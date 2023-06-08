#!/usr/bin/env python3

"""
end-to-end pangenome pipeline, which chains:
cactus-minigraph
cactus-graphmap
cactus-graphmap-split
cactus-align --batch
cactus-graphmap-join

data between the stages is exported (and imported where necessary)
to make it easier to recover from errors
"""

import os, sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit, time

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.version import cactus_commit
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count

from cactus.refmap.cactus_minigraph import minigraph_construct_workflow
from cactus.refmap.cactus_minigraph import check_sample_names
from cactus.refmap.cactus_graphmap import minigraph_workflow
from cactus.refmap.cactus_graphmap_split import graphmap_split_workflow, export_split_data
from cactus.setup.cactus_align import make_batch_align_jobs, batch_align_jobs
from cactus.refmap.cactus_graphmap_join import graphmap_join_workflow, export_join_data, graphmap_join_options, graphmap_join_validate_options

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("seqFile", help = "Seq file (will be modified if necessary to include graph Fasta sequence)")
    parser.add_argument("--outDir", help = "Output directory", required=True)
    parser.add_argument("--outName", help = "Output name (without extension)", required=True)
    parser.add_argument("--reference", required=True, nargs='+', type=str, help = "Reference event name(s). The first will be the \"true\" reference and will be left unclipped and uncollapsed. It also should have been used with --reference in all upstream commands. Other names will be promoted to reference paths in vg") 

    # cactus-minigraph options
    parser.add_argument("--mgCores", type=int, help = "Number of cores for minigraph construction (defaults to the same as --maxCores).")

    # cactus-graphmap options
    parser.add_argument("--mapCores", type=int, help = "Number of cores for minigraph map.  Overrides graphmap cpu in configuration")
    
    parser.add_argument("--delFilter", type=int, help = "Filter out split-mapping-implied deletions > Nbp (default will be \"delFilter\" from the config")
    
    # cactus-graphmap-split options
    parser.add_argument("--refContigs", nargs="*", help = "Subset to these reference contigs (multiple allowed)", default=[])
    parser.add_argument("--otherContig", type=str, help = "Lump all reference contigs unselected by above options into single one with this name")
    parser.add_argument("--permissiveContigFilter", nargs='?', const='0.25', default=None, type=float, help = "If specified, override the configuration to accept contigs so long as they have at least given fraction of coverage (0.25 if no fraction specified). This can increase sensitivity of very small, fragmented and/or diverse assemblies.")

    # cactus-align options
    parser.add_argument("--maxLen", type=int, default=10000,
                        help="Only align up to this many bases (overrides <bar bandingLimit> and <caf maxRecoverableChainLength> in configuration)[default=10000]")
    parser.add_argument("--consCores", type=int, 
                        help="Number of cores for each cactus_consolidated job (defaults to all available / maxCores on single_machine)", default=None)
    parser.add_argument("--consMemory", type=int,
                        help="Memory IN GIGABYTES for each cactus_consolidated job (defaults to an estimate based on the input data size)", default=None)
    
    # cactus-graphmap-join options
    graphmap_join_options(parser)

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

    # set defaults
    # most are progressive-specific options, or pangenome options that are not longer recommended
    # they still need to be set to defaults in order to get through the different workflows
    options.pangenome = True
    options.batch = True
    options.refFromGFA = False
    options.maskFilter = None
    options.checkpointInfo = None
    options.gpu = False
    options.root = None
    options.includeRoot = False
    options.pathOverrides = None
    options.singleCopySpecies = None
    options.barMaskFilter = None
    options.pafMaskFilter = None
    options.outVG = True
    options.outGFA = False
    options.minIdentity = None
    
    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    if not options.seqFile.startswith('s3://'):
        options.seqFile = os.path.abspath(options.seqFile)
    
    if options.outDir and not options.outDir.startswith('s3://'):
        if not os.path.isdir(options.outDir):
            os.makedirs(options.outDir)
        options.outDir = os.path.abspath(options.outDir)

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    # Try to juggle --maxCores and --consCores to give some reasonable defaults where possible
    if options.batchSystem.lower() in ['single_machine', 'singlemachine']:
        if options.maxCores is not None:
            if int(options.maxCores) <= 0:
                raise RuntimeError('Cactus requires --maxCores >= 1')
        if options.consCores is None:
            if options.maxCores is not None:
                options.consCores = int(options.maxCores)
            else:
                options.consCores = cactus_cpu_count()
        elif options.maxCores is not None and options.consCores > int(options.maxCores):
            raise RuntimeError('--consCores must be <= --maxCores')
    else:
        if not options.consCores:
            raise RuntimeError('--consCores required for non single_machine batch systems')
    if options.maxCores is not None and options.consCores > int(options.maxCores):
        raise RuntimeError('--consCores must be <= --maxCores')

    # Sort out the graphmap-join options, which can be rather complex
    # pass in dummy values for now, they will get filled in later
    # (but we want to do as much error-checking upfront as possible)
    options.hal = [None]
    options.vg = [None]
    options = graphmap_join_validate_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            toil.restart()
        else:
            # load up the seqfile
            config_node = ET.parse(options.configFile).getroot()
            config_wrapper = ConfigWrapper(config_node)
            config_wrapper.substituteAllPredefinedConstantsWithLiterals()
            graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")

            # load the seqfile
            seqFile = SeqFile(options.seqFile)
            input_seq_map = seqFile.pathMap

            # validate the sample names
            check_sample_names(input_seq_map.keys(), options.reference)

            # apply cpu override                
            if options.mapCores is not None:
                findRequiredNode(config_node, "graphmap").attrib["cpu"] = str(options.mapCores)
            mg_cores = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "cpu", typeFn=int, default=1)
            if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
                mg_cores = min(mg_cores, cactus_cpu_count(), int(options.maxCores) if options.maxCores else sys.maxsize)
                findRequiredNode(config_node, "graphmap").attrib["cpu"] = str(mg_cores)

            if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
                if not options.mgCores:
                    options.mgCores = sys.maxsize
                options.mgCores = min(options.mgCores, cactus_cpu_count(), int(options.maxCores) if options.maxCores else sys.maxsize)
            else:
                if not options.mgCores:
                    raise RuntimeError("--mgCores required run *not* running on single machine batch system")
                        
            #import the sequences
            input_seq_id_map = {}
            input_path_map = {}
            input_seq_order = [options.reference[0]]
            leaves = set([seqFile.tree.getName(node) for node in seqFile.tree.getLeaves()])
            for (genome, seq) in input_seq_map.items():
                if genome != graph_event and genome in leaves:                
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    logger.info("Importing {}".format(seq))
                    input_seq_id_map[genome] = toil.importFile(seq)
                    input_path_map[genome] = seq
                    if genome != options.reference[0]:
                        input_seq_order.append(genome)
            
            toil.start(Job.wrapJobFn(pangenome_end_to_end_workflow, options, config_wrapper, input_seq_id_map, input_path_map, input_seq_order))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-pangenome has finished after {} seconds".format(run_time))

def export_minigraph_wrapper(job, options, sv_gfa_id, sv_gfa_path):
    """ export the GFA from minigraph """
    job.fileStore.exportFile(sv_gfa_id, makeURL(os.path.join(options.outDir, os.path.basename(sv_gfa_path))))

def export_graphmap_wrapper(job, options, paf_id, paf_path, gaf_id, unfiltered_paf_id, paf_filter_log):
    """ export the PAF file from minigraph """
    paf_path = os.path.join(options.outDir, os.path.basename(paf_path))
    job.fileStore.exportFile(paf_id, makeURL(paf_path))
    if gaf_id:
        gaf_name = os.path.splitext(os.path.basename(paf_path))[0]
        job.fileStore.exportFile(gaf_id, makeURL(os.path.join(options.outDir, gaf_name)))
    if unfiltered_paf_id:
        job.fileStore.exportFile(paf_id, makeURL(paf_path + '.unfiltered.gz'))
        job.fileStore.exportFile(paf_filter_log, makeURL(paf_path + '.filter.log'))

def update_seqfile(job, options, seq_id_map, seq_path_map, seq_order, gfa_fa_id, gfa_fa_path, graph_event):
    """ put the minigraph gfa.fa file into the seqfile and export both """
    work_dir = job.fileStore.getLocalTempDir()    
    seq_id_map[graph_event] = gfa_fa_id
    seq_path_map[graph_event] = gfa_fa_path
    local_seqfile = os.path.join(work_dir, 'seqfile.txt')
    with open(local_seqfile, 'w') as seqfile:
        for sample in seq_order:
            seqfile.write('{}\t{}\n'.format(sample, seq_path_map[sample]))
        seqfile.write('{}\t{}\n'.format(graph_event, gfa_fa_path))            
    job.fileStore.exportFile(job.fileStore.writeGlobalFile(local_seqfile), makeURL(os.path.join(options.outDir, os.path.basename(options.seqFile))))
    job.fileStore.exportFile(gfa_fa_id, makeURL(os.path.join(options.outDir, os.path.basename(gfa_fa_path))))

    seq_name_map = {}
    for sample,seq_path in seq_path_map.items():
        seq_name_map[sample] = os.path.basename(seq_path)
    return seq_id_map, seq_path_map, seq_name_map
    
def export_split_wrapper(job, wf_output, out_dir, config_node):
    """ toil job wrapper for cactus_graphmap_split's exporter """
    if not out_dir.startswith('s3://') and not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    export_split_data(job.fileStore, wf_output[0], wf_output[1], wf_output[2], wf_output[3], out_dir, config_node)

def make_batch_align_jobs_wrapper(job, options, chromfile_path, config_wrapper):
    """ toil job wrapper for make_batch_align_jobs from cactus_align """
    work_dir = job.fileStore.getLocalTempDir()
    local_chromfile_path = os.path.join(work_dir, 'chromfile.txt')
    chromfile_id = job.fileStore.importFile(makeURL(chromfile_path))
    job.fileStore.readGlobalFile(chromfile_id, local_chromfile_path)
    options.seqFile = local_chromfile_path    
    align_jobs = make_batch_align_jobs(options, job.fileStore, job.fileStore, config_wrapper)
    return align_jobs
    
def export_align_wrapper(job, options, results_dict, reference_list):
    """ toil job wrapper for exporting from cactus_align """
    vg_ids = []
    vg_paths = []
    hal_ids = []
    hal_paths = []
    chrom_dir = os.path.join(options.outDir, 'chrom-alignments')
    if not chrom_dir.startswith('s3://') and not os.path.isdir(chrom_dir):
        os.makedirs(chrom_dir)

    for chrom, results in results_dict.items():
        hal_path = makeURL(os.path.join(chrom_dir, '{}.hal'.format(chrom)))
        vg_path = makeURL(os.path.join(chrom_dir, '{}.vg'.format(chrom)))
        job.fileStore.exportFile(results[0], hal_path)
        job.fileStore.exportFile(results[1], vg_path)
        hal_ids.append(results[0])
        hal_paths.append(hal_path)
        vg_ids.append(results[1])
        vg_paths.append(vg_path)

    join_options = options
    join_options.hal = hal_paths
    join_options.vg = vg_paths
    join_options.reference = reference_list

    return join_options, vg_ids, hal_ids

def export_join_wrapper(job, options, wf_output):
    """ toil join wrapper for cactus_graphmap_join """
    export_join_data(job.fileStore, options, wf_output[0], wf_output[1], wf_output[2], wf_output[3], wf_output[4])

def pangenome_end_to_end_workflow(job, options, config_wrapper, seq_id_map, seq_path_map, seq_order):
    """ chain the entire workflow together, doing exports after each step to mitigate annoyance of failures """
    root_job = Job()
    job.addChild(root_job)
    config_node = config_wrapper.xmlRoot

    # it's really only graphmap_join that can accept (or use) multiple references
    # so we forget about it until then
    assert type(options.reference) == list
    reference_list = options.reference
    options.reference = options.reference[0]

    # sanitize headers (once here, skip in all workflows below)
    sanitize_job = root_job.addFollowOnJobFn(sanitize_fasta_headers, seq_id_map, pangenome=True)
    seq_id_map = sanitize_job.rv()
    
    # cactus_minigraph
    sv_gfa_path = os.path.join(options.outDir, options.outName + '.sv.gfa.gz')
    
    minigraph_job = sanitize_job.addFollowOnJobFn(minigraph_construct_workflow, options, config_node, seq_id_map, seq_order, sv_gfa_path, sanitize=False)
    sv_gfa_id = minigraph_job.rv()
    minigraph_job.addFollowOnJobFn(export_minigraph_wrapper, options, sv_gfa_id, sv_gfa_path)

    # cactus_graphmap
    paf_path = os.path.join(options.outDir, options.outName + '.paf')
    gfa_fa_path = os.path.join(options.outDir, options.outName + '.gfa.fa.gz')
    options.minigraphGFA = sv_gfa_path
    options.outputFasta = gfa_fa_path
    graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    
    graphmap_job = minigraph_job.addFollowOnJobFn(minigraph_workflow, options, config_wrapper, seq_id_map, sv_gfa_id, graph_event, sanitize=False)
    paf_id, gfa_fa_id, gaf_id, unfiltered_paf_id, paf_filter_log = graphmap_job.rv(0), graphmap_job.rv(1), graphmap_job.rv(2), graphmap_job.rv(3), graphmap_job.rv(4)
    graphmap_job.addFollowOnJobFn(export_graphmap_wrapper, options, paf_id, paf_path, gaf_id, unfiltered_paf_id, paf_filter_log)

    # we need to update the seqfile with the phonied in minigraph event
    update_seqfile_job = graphmap_job.addFollowOnJobFn(update_seqfile, options, seq_id_map, seq_path_map, seq_order, gfa_fa_id, gfa_fa_path, graph_event)
    seq_id_map, seq_path_map, seq_name_map = update_seqfile_job.rv(0), update_seqfile_job.rv(1), update_seqfile_job.rv(2)

    # cactus_graphmap_split
    split_job = update_seqfile_job.addFollowOnJobFn(graphmap_split_workflow, options, config_wrapper, seq_id_map, seq_name_map, sv_gfa_id,
                                                    sv_gfa_path, paf_id, paf_path, options.refContigs, options.otherContig, sanitize=False)

    wf_output = split_job.rv()
    split_out_path = os.path.join(options.outDir, 'chrom-subproblems')    
    split_export_job = split_job.addFollowOnJobFn(export_split_wrapper, wf_output, split_out_path, config_wrapper)
        
    # cactus_align
    chromfile_path = os.path.join(split_out_path, 'chromfile.txt')

    align_jobs_make_job = split_export_job.addFollowOnJobFn(make_batch_align_jobs_wrapper, options, chromfile_path, config_wrapper)
    align_jobs = align_jobs_make_job.rv()
    align_job = align_jobs_make_job.addFollowOnJobFn(batch_align_jobs, align_jobs)
    results_dict = align_job.rv()
    align_export_job = align_job.addFollowOnJobFn(export_align_wrapper, options, results_dict, reference_list)
    join_options, vg_ids, hal_ids = align_export_job.rv(0), align_export_job.rv(1), align_export_job.rv(2)

    # cactus_graphmap_join
    join_job = align_export_job.addFollowOnJobFn(graphmap_join_workflow, join_options, config_wrapper, vg_ids, hal_ids)
    join_wf_output = join_job.rv()
    join_export_job = join_job.addFollowOnJobFn(export_join_wrapper, join_options, join_wf_output)
        
        
