#!/usr/bin/env python3

"""This script will take the output of cactus-align-batch (which was run on the output of cactus-graphmap-split).  This would be
a set of .vg files, one for each chromosome.  It will do the following in order to make a single output graph
- clip unmapped regions out of each graph
- join the ids of the clipped graph
- convert each clipped graph to gfa
- merge the gfas into a whole-genome graph
- convert the gfa into a gbwt/gbwt-graph pair
- export the gbwt to xg
- export the xg/gbwt to vcf

so the final output will be
- set of clipped vg graphs
- graph.gfa.gz
- graph.gbwt
- graph.gg
- graph.xg
- graph.snarls
- graph.vcf.gz
- graph.vcf.gz.tbi
"""
import os
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.pipeline.cactus_workflow import CactusWorkflowArguments
from cactus.pipeline.cactus_workflow import addCactusWorkflowOptions
from cactus.pipeline.cactus_workflow import CactusTrimmingBlastPhase
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import unzip_gz, write_s3
from cactus.preprocessor.fileMasking import get_mask_bed_from_fasta
from toil.job import Job
from toil.common import Toil
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions
from toil.realtimeLogger import RealtimeLogger
from toil.lib.threading import cpu_count

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory, getTempFile, catFiles

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    parser.add_argument("--vg", required=True, nargs='+',  help = "Input vg files")
    parser.add_argument("--outDir", required=True, type=str, help = "Output directory")
    parser.add_argument("--outName", required=True, type=str, help = "Basename of all output files")
    parser.add_argument("--reference", required=True, type=str, help = "Reference event name")
    parser.add_argument("--rename", nargs='+', default = [], help = "Path renaming, each of form src>dest (see clip-vg -r)")
    parser.add_argument("--clipLength", type=int, default=0, help = "clip out unaligned sequences longer than this")
    parser.add_argument("--wlineSep", type=str, help = "wline separator for vg convert")
    parser.add_argument("--indexCores", type=int, default=1, help = "cores for indexing processes")
    
    #Progressive Cactus Options
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)

    options = parser.parse_args()

    setupBinaries(options)
    setLoggingFromOptions(options)
    enableDumpStack()

    if options.outDir and not options.outDir.startswith('s3://'):
        if not os.path.isdir(options.outDir):
            os.makedirs(options.outDir)
        
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    start_time = timeit.default_timer()
    runCactusGraphMapJoin(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-graphmap-join has finished after {} seconds".format(run_time))

def runCactusGraphMapJoin(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            wf_output = toil.restart()
        else:
            options.cactusDir = getTempDirectory()

            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals()

            # load up the vgs
            vg_ids = []
            for vg_path in options.vg:
                logger.info("Importing {}".format(vg_path))
                vg_ids.append(toil.importFile(makeURL(vg_path)))

            # run the workflow
            wf_output = toil.start(Job.wrapJobFn(graphmap_join_workflow, options, config, vg_ids))
                
        #export the split data
        export_join_data(toil, options, wf_output[0], wf_output[1])

def graphmap_join_workflow(job, options, config, vg_ids):

    root_job = Job()
    job.addChild(root_job)

    # run clip-vg on each input
    clipped_vg_ids = []
    for vg_path, vg_id in zip(options.vg, vg_ids):
        clip_job = root_job.addChildJobFn(clip_vg, options, config, vg_path, vg_id,
                                          disk=vg_id.size * 2, memory=vg_id.size * 4)
        clipped_vg_ids.append(clip_job.rv())

    # join the ids
    join_job = root_job.addFollowOnJobFn(join_vg, options, config, clipped_vg_ids,
                                         disk=sum([f.size for f in vg_ids]))
    clipped_vg_ids = join_job.rv()

    # make a gfa for each
    gfa_root_job = Job()
    join_job.addFollowOn(gfa_root_job)
    clipped_gfa_ids = []
    for i in range(len(options.vg)):
        vg_path = options.vg[i]
        clipped_id = join_job.rv(i)
        vg_id = vg_ids[i]
        gfa_job = gfa_root_job.addChildJobFn(vg_to_gfa, options, config, vg_path, clipped_id,
                                             disk=vg_id.size * 5)
        clipped_gfa_ids.append(gfa_job.rv())

    # merge up the gfas and make the various vg indexes
    gfa_merge_job = gfa_root_job.addFollowOnJobFn(vg_indexes, options, config, clipped_gfa_ids,
                                                  cores=options.indexCores,
                                                  disk=sum(f.size for f in vg_ids) * 10)
    
    return clipped_vg_ids, gfa_merge_job.rv()

def clip_vg(job, options, config, vg_path, vg_id):
    """ run clip-vg 
    """
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    out_path = vg_path + '.clip'

    cmd = ['clip-vg', vg_path, '-u', str(options.clipLength), '-f']
    for rs in options.rename:
        cmd += ['-r', rs]
    if options.reference:
        cmd += ['-e', options.reference]

    # sort while we're at it
    cmd = [cmd, ['vg', 'ids', '-s', '-']]
        
    cactus_call(parameters=cmd, outfile=out_path)

    # worth it
    cactus_call(parameters=['vg', 'validate', out_path])

    return job.fileStore.writeGlobalFile(out_path)

def join_vg(job, options, config, clipped_vg_ids):
    """ run vg ids -j
    """
    work_dir = job.fileStore.getLocalTempDir()
    vg_paths = []

    for vg_path, vg_id in zip(options.vg, clipped_vg_ids):
        vg_path = os.path.join(work_dir, os.path.basename(vg_path))
        job.fileStore.readGlobalFile(vg_id, vg_path, mutable=True)
        vg_paths.append(vg_path)
        
    cactus_call(parameters=['vg', 'ids', '-j'] + vg_paths)

    return [job.fileStore.writeGlobalFile(f) for f in vg_paths]

def vg_to_gfa(job, options, config, vg_path, vg_id):
    """ run gfa conversion """
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    out_path = vg_path + '.gfa'

    cmd = ['vg', 'convert', '-f', '-Q', options.reference, vg_path, '-B']
    if options.wlineSep:
        cmd += ['-w', options.wlineSep]

    cactus_call(parameters=cmd, outfile=out_path)

    return job.fileStore.writeGlobalFile(out_path)

def vg_indexes(job, options, config, gfa_ids):
    """ merge of the gfas, then make gbwt / xg / snarls / vcf 
    """ 
    work_dir = job.fileStore.getLocalTempDir()
    vg_paths = []
    merge_gfa_path = os.path.join(work_dir, 'merged.gfa')
    with open(merge_gfa_path, 'w') as merge_gfa_file:
        merge_gfa_file.write('H\tVN:Z:1.0\n')

    # merge the gfas
    for vg_path, gfa_id in zip(options.vg, gfa_ids):
        gfa_path = os.path.join(work_dir, os.path.basename(vg_path) +  '.gfa')
        job.fileStore.readGlobalFile(gfa_id, gfa_path, mutable=True)
        cactus_call(parameters=['grep', '-v', '^H', gfa_path], outfile=merge_gfa_path, outappend=True)
        os.remove(gfa_path)

    # make the gbwt
    gbwt_path = os.path.join(work_dir, 'merged.gbwt')
    gg_path = os.path.join(work_dir, 'merged.gg')
    trans_path = os.path.join(work_dir, 'merged.trans')
    cactus_call(parameters=['vg', 'gbwt', '-G', merge_gfa_path, '-o', gbwt_path, '-g', gg_path, '--translation', trans_path])

    # zip the gfa
    cactus_call(parameters=['bgzip', merge_gfa_path, '--threads', str(job.cores)])
    gfa_path = merge_gfa_path + '.gz'

    # make the xg
    xg_path = os.path.join(work_dir, 'merged.xg')
    cactus_call(parameters=['vg', 'convert', gg_path, '-b', gbwt_path, '-x', '-t', str(job.cores)], outfile=xg_path)

    # worth it
    cactus_call(parameters=['vg', 'validate', xg_path])

    # make the snarls
    snarls_path = os.path.join(work_dir, 'merged.snarls')
    cactus_call(parameters=['vg', 'snarls', xg_path, '-t', str(job.cores)], outfile=snarls_path)

    # make the vcf
    vcf_path = os.path.join(work_dir, 'merged.vcf.gz')
    cactus_call(parameters=[['vg', 'deconstruct', xg_path, '-r', snarls_path, '-g', gbwt_path, '-T', trans_path, '-t', str(job.cores)],
                            ['bgzip', '--threads', str(job.cores)]],
                outfile=vcf_path)
    cactus_call(parameters=['tabix', '-p', 'vcf', vcf_path])
                            
    return { 'gfa.gz' : job.fileStore.writeGlobalFile(gfa_path),
             'gbwt' : job.fileStore.writeGlobalFile(gbwt_path),
             'gg' : job.fileStore.writeGlobalFile(gg_path),
             'trans' : job.fileStore.writeGlobalFile(trans_path),
             'xg' : job.fileStore.writeGlobalFile(xg_path),
             'snarls' : job.fileStore.writeGlobalFile(snarls_path),
             'vcf.gz' : job.fileStore.writeGlobalFile(vcf_path),
             'tbi' : job.fileStore.writeGlobalFile(vcf_path + '.tbi') }

def export_join_data(toil, options, clip_ids, idx_map):
    """ download all the output data
    """

    # download the clip vgs
    clip_base = os.path.join(options.outDir, 'clip')
    if not clip_base.startswith('s3://') and not os.path.isdir(clip_base):
        os.makedirs(clip_base)

    for vg_path, vg_id in zip(options.vg, clip_ids):
        toil.exportFile(vg_id, makeURL(os.path.join(clip_base, os.path.basename(vg_path))))

    # download everything else
    for ext, idx_id in idx_map.items():
        toil.exportFile(idx_id, makeURL(os.path.join(options.outDir, '{}.{}'.format(options.outName, ext))))
    
        
if __name__ == "__main__":
    main()