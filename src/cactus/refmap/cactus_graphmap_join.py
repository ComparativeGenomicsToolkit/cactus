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
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import unzip_gz, write_s3
from cactus.preprocessor.fileMasking import get_mask_bed_from_fasta
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from toil.lib.threading import cpu_count

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory, getTempFile, catFiles

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    parser.add_argument("--vg", required=True, nargs='+',  help = "Input vg files (PackedGraph or HashGraph format)")
    parser.add_argument("--outDir", required=True, type=str, help = "Output directory")
    parser.add_argument("--outName", required=True, type=str, help = "Basename of all output files")
    parser.add_argument("--reference", required=True, type=str, help = "Reference event name")
    parser.add_argument("--vcfReference", type=str, help = "Produce additional VCF for given reference event")
    parser.add_argument("--rename", nargs='+', default = [], help = "Path renaming, each of form src>dest (see clip-vg -r)")
    parser.add_argument("--clipLength", type=int, default=None, help = "clip out unaligned sequences longer than this")
    parser.add_argument("--wlineSep", type=str, help = "wline separator for vg convert")
    parser.add_argument("--indexCores", type=int, default=1, help = "cores for general indexing and VCF constructions")
    parser.add_argument("--giraffeCores", type=int, default=None, help = "cores for giraffe-specific indexing (defaults to --indexCores)")
    parser.add_argument("--decoyGraph", help= "decoy sequences vg graph to add (PackedGraph or HashGraph format)")
    parser.add_argument("--hal", nargs='+', default = [], help = "Input hal files (for merging)")
    parser.add_argument("--vcf", action="store_true", help= "make VCF")
    parser.add_argument("--giraffe", action="store_true", help= "make Giraffe-specific indexes (distance and minimizer)")
    parser.add_argument("--normalizeIterations", type=int, default=None,
                        help="Run this many iterations of vg normamlization (shared prefix zipping)")
    parser.add_argument("--gfaffix", action="store_true",
                        help="Run GFAFfix normalization")
    parser.add_argument("--vgClipOpts", nargs='+', help = "If specified, run vg clip with given options (surround in quotes; multiple allowed to chain multiple commands)")
    
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
    set_logging_from_options(options)
    enableDumpStack()

    if options.outDir and not options.outDir.startswith('s3://'):
        if not os.path.isdir(options.outDir):
            os.makedirs(options.outDir)

    if options.hal and len(options.hal) != len(options.vg):
        raise RuntimeError("If --hal and --vg should specify the same number of files")
    if not options.giraffeCores:
        options.giraffeCores = options.indexCores
        
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

            # tack on the decoys
            if options.decoyGraph:
                logger.info("Importing decoys {}".format(options.decoyGraph))
                vg_ids.append(toil.importFile(makeURL(options.decoyGraph)))
                # we'll treat it like any other graph downstream, except clipping
                # where we'll check first using the path name
                options.vg.append(options.decoyGraph)
                
            # load up the hals
            hal_ids = []
            for hal_path in options.hal:
                logger.info("Importing {}".format(hal_path))
                hal_ids.append(toil.importFile(makeURL(hal_path)))                

            # run the workflow
            wf_output = toil.start(Job.wrapJobFn(graphmap_join_workflow, options, config, vg_ids, hal_ids))
                
        #export the split data
        export_join_data(toil, options, wf_output[0], wf_output[1])

def graphmap_join_workflow(job, options, config, vg_ids, hal_ids):

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

    # optional clipping -- we do this down here after joining and normalization so
    # our graph is id-compatible with a graph that wasn't clipped (but run with same parameters otherwise)
    if options.vgClipOpts:
        clip_root_job = Job()
        join_job.addFollowOn(clip_root_job)
        clipped_vg_ids = []
        for i in range(len(vg_ids)):
            vg_clip_job = clip_root_job.addChildJobFn(vg_clip_vg, options, config, options.vg[i], join_job.rv(i),
                                                      disk=vg_ids[i].size * 2)
            join_job.addFollowOn(vg_clip_job)
            clipped_vg_ids.append(vg_clip_job.rv())
        join_job = clip_root_job

    # make a gfa for each
    gfa_root_job = Job()
    join_job.addFollowOn(gfa_root_job)
    clipped_gfa_ids = []
    for i in range(len(options.vg)):
        vg_path = options.vg[i]
        clipped_id = clipped_vg_ids[i] if options.vgClipOpts else join_job.rv(i)
        vg_id = vg_ids[i]
        gfa_job = gfa_root_job.addChildJobFn(vg_to_gfa, options, config, vg_path, clipped_id,
                                             disk=vg_id.size * 5)
        clipped_gfa_ids.append(gfa_job.rv())

    # merge up the gfas and make the various vg indexes
    gfa_merge_job = gfa_root_job.addFollowOnJobFn(vg_indexes, options, config, clipped_gfa_ids,
                                                  cores=options.indexCores,
                                                  disk=sum(f.size for f in vg_ids) * 5)
    out_dicts = [gfa_merge_job.rv()]
                                
    # optional vcf
    if options.vcf:
        deconstruct_job = gfa_merge_job.addFollowOnJobFn(make_vcf, options.outName, options.reference, out_dicts[0],
                                                         cores=options.indexCores,
                                                         disk = sum(f.size for f in vg_ids) * 2)
        out_dicts.append(deconstruct_job.rv())

    # optional vcf with different reference
    if options.vcfReference:
        ref_deconstruct_job = gfa_merge_job.addFollowOnJobFn(make_vcf, options.outName, options.vcfReference, out_dicts[0],
                                                             tag=options.vcfReference + '.',
                                                             cores=options.indexCores,
                                                             disk = sum(f.size for f in vg_ids) * 2)
        out_dicts.append(ref_deconstruct_job.rv())

    # optional giraffe
    if options.giraffe:
        giraffe_job = gfa_merge_job.addFollowOnJobFn(make_giraffe_indexes, options, out_dicts[0],
                                                     cores=options.giraffeCores,
                                                     disk = sum(f.size for f in vg_ids) * 4)
        out_dicts.append(giraffe_job.rv())

    # optional hal merge
    if hal_ids:
        hal_merge_job = root_job.addChildJobFn(merge_hal, options, hal_ids,
                                               cores = 1,
                                               disk=sum(f.size for f in hal_ids) * 2)
        out_dicts.append(hal_merge_job.rv())
    
    return clipped_vg_ids, out_dicts

def clip_vg(job, options, config, vg_path, vg_id):
    """ run clip-vg 
    """
    work_dir = job.fileStore.getLocalTempDir()
    is_decoy = vg_path == options.decoyGraph
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    clipped_path = vg_path + '.clip'
    out_path = vg_path + '.out'

    # remove masked unaligned regions with clip-vg
    cmd = ['clip-vg', vg_path, '-f']
    if options.clipLength is not None and not is_decoy:
        cmd += ['-u', str(options.clipLength)]
    for rs in options.rename:
        cmd += ['-r', rs]
    if options.reference:
        cmd += ['-e', options.reference]

    if getOptionalAttrib(findRequiredNode(config.xmlRoot, "hal2vg"), "includeMinigraph", typeFn=bool, default=False):
        # our vg file has minigraph sequences -- we'll filter them out, along with any nodes
        # that don't appear in a non-minigraph path
        graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
        cmd += ['-d', graph_event]
        
    cactus_call(parameters=cmd, outfile=clipped_path)
        
    # optional normalization.  this will (in theory) correct a lot of small underalignments due to cactus bugs
    # by zipping up redundant nodes. done before clip-vg otherwise ref paths not guaranteed to be forwardized
    # todo: would be nice if clip-vg worked on stdin
    if options.normalizeIterations:
        normalized_path = clipped_path + '.normalized'
        mod_cmd = ['vg', 'mod', '-O', '-U', str(options.normalizeIterations), clipped_path]
        if options.reference:
            mod_cmd += ['-z', options.reference]
        cactus_call(parameters=mod_cmd, outfile=normalized_path)
        # worth it
        cactus_call(parameters=['vg', 'validate', normalized_path])
        clipped_path = normalized_path

    # GFAFfix is a tool written just to normalize graphs, and alters a (faster) alternative to vg
    # (though requires GFA conversion)
    if options.gfaffix:
        normalized_path = clipped_path + '.gfaffixed'
        gfa_in_path = vg_path + '.gfa'
        gfa_out_path = normalized_path + '.gfa'
        cactus_call(parameters=['vg', 'convert', '-f', clipped_path], outfile=gfa_in_path)
        fix_cmd = ['gfaffix', gfa_in_path, '--output_refined', gfa_out_path]
        if options.reference:
            fix_cmd += ['--dont_collapse', options.reference + '*']
        cactus_call(parameters=fix_cmd)
        # GFAFfix doesn't seem to unchop that well, so we do that in vg after
        # The double convert is to work around: https://github.com/vgteam/vg/issues/3373
        cactus_call(parameters=[['vg', 'convert', '-g', '-a', gfa_out_path], ['vg', 'mod', '-u', '-'], ['vg', 'convert', '-p', '-']], outfile=normalized_path)
        clipped_path = normalized_path

    # also forwardize just in case
    if options.reference and (options.normalizeIterations or options.gfaffix):
        forward_path = clipped_path + '.forward'
        cmd = ['clip-vg', clipped_path, '-e', options.reference]
        cactus_call(parameters=cmd, outfile=forward_path)
        clipped_path = forward_path

    # sort by id
    cmd = ['vg', 'ids', '-s', clipped_path]        
    cactus_call(parameters=cmd, outfile=out_path)

    # worth it
    cactus_call(parameters=['vg', 'validate', out_path])

    return job.fileStore.writeGlobalFile(out_path)

def vg_clip_vg(job , options, config, vg_path, vg_id):
    """ run vg clip, chaining multiple invocations if desired.
    """
    work_dir = job.fileStore.getLocalTempDir()
    is_decoy = vg_path == options.decoyGraph
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    clipped_path = vg_path + '.clip'

    clip_cmd = ['vg', 'clip', vg_path, '-P', options.reference] +  options.vgClipOpts[0].split()
    if len(options.vgClipOpts) > 1:
        clip_cmd = [clip_cmd]
        for clip_opts in options.vgClipOpts[1:]:
            clip_cmd.append(['vg', 'clip', '-', '-P', options.reference] +  clip_opts.split())

    cactus_call(parameters=clip_cmd, outfile=clipped_path)

    # worth it
    cactus_call(parameters=['vg', 'validate', clipped_path])

    return job.fileStore.writeGlobalFile(clipped_path)
    
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

    cmd = ['vg', 'convert', '-f', '-Q', options.reference, os.path.basename(vg_path), '-B']
    if options.wlineSep:
        cmd += ['-w', options.wlineSep]

    # important, when options.wlineSep is ., it throws off prepareWorkDir in cactus_call
    # so important to specify the work_dir below       
    cactus_call(parameters=cmd, outfile=out_path, work_dir=work_dir)

    return job.fileStore.writeGlobalFile(out_path)

def vg_indexes(job, options, config, gfa_ids):
    """ merge of the gfas, then make gbwt / xg / snarls / trans
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

    # compress the trans
    cactus_call(parameters=['bgzip', trans_path, '--threads', str(job.cores)])
    trans_path += '.gz'                            
    
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
    cactus_call(parameters=['vg', 'snarls', xg_path, '-T', '-t', str(job.cores)], outfile=snarls_path)
                            
    return { 'gfa.gz' : job.fileStore.writeGlobalFile(gfa_path),
             'gbwt' : job.fileStore.writeGlobalFile(gbwt_path),
             'gg' : job.fileStore.writeGlobalFile(gg_path),
             'trans.gz' : job.fileStore.writeGlobalFile(trans_path),
             'xg' : job.fileStore.writeGlobalFile(xg_path),
             'snarls' : job.fileStore.writeGlobalFile(snarls_path) }

def make_vcf(job, out_name, vcf_ref, index_dict, tag=''):
    """ make the vcf
    """ 
    work_dir = job.fileStore.getLocalTempDir()
    xg_path = os.path.join(work_dir, os.path.basename(out_name) + '.xg')
    gbwt_path = os.path.join(work_dir, os.path.basename(out_name) + '.gbwt')
    snarls_path = os.path.join(work_dir, os.path.basename(out_name) + '.snarls')
    trans_path = os.path.join(work_dir, os.path.basename(out_name) + '.trans.gz')
    job.fileStore.readGlobalFile(index_dict['xg'], xg_path)
    job.fileStore.readGlobalFile(index_dict['gbwt'], gbwt_path)
    job.fileStore.readGlobalFile(index_dict['snarls'], snarls_path)
    job.fileStore.readGlobalFile(index_dict['trans.gz'], trans_path)

    # unzip the trans
    cactus_call(parameters=['gzip', '-fd', trans_path])
    trans_path = trans_path[:-3]

    # make the vcf
    vcf_path = os.path.join(work_dir, 'merged.vcf.gz')
    cactus_call(parameters=[['vg', 'deconstruct', xg_path, '-P', vcf_ref, '-a', '-r', snarls_path, '-g', gbwt_path,
                             '-T', trans_path, '-t', str(job.cores)],
                            ['bgzip', '--threads', str(job.cores)]],
                outfile=vcf_path)
    cactus_call(parameters=['tabix', '-p', 'vcf', vcf_path])

    return { '{}vcf.gz'.format(tag) : job.fileStore.writeGlobalFile(vcf_path),
             '{}vcf.gz.tbi'.format(tag) : job.fileStore.writeGlobalFile(vcf_path + '.tbi') }
    
def make_giraffe_indexes(job, options, index_dict):
    """ make giraffe-specific indexes: distance and minimaer """
    work_dir = job.fileStore.getLocalTempDir()
    xg_path = os.path.join(work_dir, os.path.basename(options.outName) + '.xg')
    gbwt_path = os.path.join(work_dir, os.path.basename(options.outName) + '.gbwt')
    snarls_path = os.path.join(work_dir, os.path.basename(options.outName) + '.snarls')
    job.fileStore.readGlobalFile(index_dict['xg'], xg_path)
    job.fileStore.readGlobalFile(index_dict['gbwt'], gbwt_path)
    job.fileStore.readGlobalFile(index_dict['snarls'], snarls_path)

    # make the distance index
    dist_path = os.path.join(work_dir, os.path.basename(options.outName) + '.dist')
    cactus_call(parameters=['vg', 'index', '-t', str(job.cores), '-j', dist_path, '-s', snarls_path, xg_path])

    # make the minimizer index
    min_path = os.path.join(work_dir, os.path.basename(options.outName) + '.min')
    cactus_call(parameters=['vg', 'minimizer', '-k', '29', '-w', '11', '-t', str(job.cores), '-g', gbwt_path,
                            '-d', dist_path, '-o', min_path, xg_path])                            
                            
    return { 'min' : job.fileStore.writeGlobalFile(min_path),
             'dist' : job.fileStore.writeGlobalFile(dist_path) }

def merge_hal(job, options, hal_ids):
    """ call halMergeChroms to make one big hal file out of the chromosome hal files """
    work_dir = job.fileStore.getLocalTempDir()
    hal_paths = []
    for in_path, hal_id in zip(options.hal, hal_ids):
        hal_path = os.path.join(work_dir, os.path.basename(in_path))
        job.fileStore.readGlobalFile(hal_id, hal_path)
        hal_paths.append(hal_path)

    merged_path = os.path.join(work_dir, '__merged__.hal')
    assert merged_path not in hal_paths

    # note: cactus_call tries to sort out relative paths by itself for docker.  but the comma-separated list
    # will likely throw it off, so we take care to specify it relative manually.
    # also note: most hal commands need --inMemory to run at scale, but the access patterns for chrom
    # merging are linear enough that it shouldn't be needed
    cmd = ['halMergeChroms',
           ','.join([os.path.basename(p) for p in hal_paths]),
           os.path.basename(merged_path),
           '--progress']
    cactus_call(parameters=cmd, work_dir = work_dir)

    return { 'hal' : job.fileStore.writeGlobalFile(merged_path) }

def export_join_data(toil, options, clip_ids, idx_maps):
    """ download all the output data
    """

    # download the clip vgs
    clip_base = os.path.join(options.outDir, 'clip-{}'.format(options.outName))
    if not clip_base.startswith('s3://') and not os.path.isdir(clip_base):
        os.makedirs(clip_base)

    for vg_path, vg_id in zip(options.vg, clip_ids):
        toil.exportFile(vg_id, makeURL(os.path.join(clip_base, os.path.basename(vg_path))))

    # download everything else
    for idx_map in idx_maps:
        for ext, idx_id in idx_map.items():
            toil.exportFile(idx_id, makeURL(os.path.join(options.outDir, '{}.{}'.format(options.outName, ext))))
        
if __name__ == "__main__":
    main()
