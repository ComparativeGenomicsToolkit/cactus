#!/usr/bin/env python3

"""This sciprt will take the cactus-graphmap output and input and split it up into chromosomes.  This is done after a whole-genome graphmap as we need the whole-genome minigraph alignments to do the chromosome splitting in the first place.  For each chromosome, it will make a PAF, GFA and Seqfile (pointing to chrosome fastas)
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
from cactus.progressive.schedule import Schedule
from cactus.progressive.projectWrapper import ProjectWrapper
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

    parser.add_argument("seqFile", help = "Seq file (gzipped fastas supported)")
    parser.add_argument("minigraphGFA", help = "Minigraph-compatible reference graph in GFA format (can be gzipped)")
    parser.add_argument("graphmapPAF", type=str, help = "Output pairwise alignment file in PAF format (can be gzipped)")
    parser.add_argument("--outDir", required=True, type=str, help = "Output directory")
    parser.add_argument("--fastaHeaderTable", type=str,
                        help="Fasta contig information (from cactus-preprocess) required to not use minigraph node sequences")    
    parser.add_argument("--refContigs", nargs="*", help = "Subset to these reference contigs (multiple allowed)", default=[])
    parser.add_argument("--refContigsFile", type=str, help = "Subset to (newline-separated) reference contigs in this file")
    parser.add_argument("--otherContig", type=str, help = "Lump all reference contigs unselected by above options into single one with this name")
    parser.add_argument("--reference", type=str, help = "Name of reference (in seqFile).  Ambiguity filters will not be applied to it")
    parser.add_argument("--maskFilter", type=int, help = "Ignore softmasked sequence intervals > Nbp")
    
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
        
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    start_time = timeit.default_timer()
    runCactusGraphMapSplit(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-graphmap-split has finished after {} seconds".format(run_time))

def runCactusGraphMapSplit(options):
    with Toil(options) as toil:
        importSingularityImage(options)

        #load cactus config
        configNode = ET.parse(options.configFile).getroot()
        config = ConfigWrapper(configNode)
        config.substituteAllPredefinedConstantsWithLiterals()

        #Run the workflow
        if options.restart:
            wf_output = toil.restart()
        else:
            options.cactusDir = getTempDirectory()

            # load up the contigs if any
            ref_contigs = set(options.refContigs)
            # todo: use import?
            if options.refContigsFile:
                with open(options.refContigsFile, 'r') as rc_file:
                    for line in rc_file:
                        if len(line.strip()):
                            ref_contigs.add(line.strip().split()[0])

            if options.otherContig:
                assert options.otherContig not in ref_contigs

            # get the minigraph "virutal" assembly name
            graph_event = getOptionalAttrib(findRequiredNode(configNode, "graphmap"), "assemblyName", default="_MINIGRAPH_")

            # import the fasta header table
            header_table_id = toil.importFile(makeURL(options.fastaHeaderTable)) if options.fastaHeaderTable else None

            # load the seqfile
            seqFile = SeqFile(options.seqFile)
            
            #import the graph
            gfa_id = toil.importFile(makeURL(options.minigraphGFA))

            #import the paf
            paf_id = toil.importFile(makeURL(options.graphmapPAF))

            #import the sequences (that we need to align for the given event, ie leaves and outgroups)
            seqIDMap = {}
            leaves = set([seqFile.tree.getName(node) for node in seqFile.tree.getLeaves()])

            options.post_stable = False
            if graph_event not in leaves:
                if not options.fastaHeaderTable:
                    raise RuntimeError("Minigraph name {} not found in seqfile".format(graph_event))
            elif options.fastaHeaderTable:
                options.post_stable = True
            if options.reference and options.reference not in leaves:
                raise RuntimeError("Name given with --reference {} not found in seqfile".format(options.reference))
                
            for genome, seq in seqFile.pathMap.items():
                if genome in leaves:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    logger.info("Importing {}".format(seq))
                    seqIDMap[genome] = (seq, toil.importFile(seq))

            # run the workflow
            wf_output = toil.start(Job.wrapJobFn(graphmap_split_workflow, options, config, seqIDMap,
                                                 gfa_id, options.minigraphGFA,
                                                 paf_id, options.graphmapPAF, ref_contigs, options.otherContig,
                                                 header_table_id))

        #export the split data
        export_split_data(toil, wf_output[0], wf_output[1], wf_output[2:], options.outDir, config, options.post_stable)

def graphmap_split_workflow(job, options, config, seqIDMap, gfa_id, gfa_path, paf_id, paf_path, ref_contigs, other_contig, header_table_id):

    root_job = Job()
    job.addChild(root_job)

    # get the sizes before we overwrite below
    gfa_size = gfa_id.size
    paf_size = paf_id.size
    
    # use file extension to sniff out compressed input
    if gfa_path.endswith(".gz"):
        gfa_id = root_job.addChildJobFn(unzip_gz, gfa_path, gfa_id, disk=gfa_id.size * 10).rv()
        gfa_size *= 10
    if paf_path.endswith(".gz"):
        paf_id = root_job.addChildJobFn(unzip_gz, paf_path, paf_id, disk=paf_id.size * 10).rv()
        paf_size *= 10

    mask_bed_id = None
    if options.maskFilter:
        mask_bed_id = root_job.addChildJobFn(get_mask_bed, seqIDMap, options.maskFilter).rv()
        
    # use rgfa-split to split the gfa and paf up by contig
    split_gfa_job = root_job.addFollowOnJobFn(split_gfa, config, gfa_id, [paf_id], ref_contigs,
                                              other_contig, options.reference, mask_bed_id,
                                              header_table_id if not options.post_stable else None,
                                              disk=(gfa_size + paf_size) * 5)

    # use the output of the above splitting to do the fasta splitting
    split_fas_job = split_gfa_job.addFollowOnJobFn(split_fas, seqIDMap, split_gfa_job.rv(0))

    # gather everythign up into a table
    gather_fas_job = split_fas_job.addFollowOnJobFn(gather_fas, seqIDMap, split_gfa_job.rv(0), split_fas_job.rv())

    # try splitting the ambiguous sequences using minimap2, which is more sensitive in some cases
    remap_job = gather_fas_job.addFollowOnJobFn(split_minimap_fallback, options, config, seqIDMap, gather_fas_job.rv())

    # partition these into fasta files
    split_fallback_gfa_job = remap_job.addFollowOnJobFn(split_gfa, config, None, remap_job.rv(0), ref_contigs,
                                                        other_contig, options.reference, None, None,
                                                        disk=(gfa_size + paf_size) * 5)

    # use the output of the above to split the ambiguous fastas
    split_fallback_fas_job = split_fallback_gfa_job.addFollowOnJobFn(split_fas, remap_job.rv(1), split_fallback_gfa_job.rv(0))

    # gather the fallback contigs into a table
    gather_fallback_fas_job = split_fallback_fas_job.addFollowOnJobFn(gather_fas, remap_job.rv(1), split_fallback_gfa_job.rv(0),
                                                                      split_fallback_fas_job.rv())

    # combine the split sequences with the split ambigious sequences
    combine_split_job = gather_fallback_fas_job.addFollowOnJobFn(combine_splits, options, config, seqIDMap, gather_fas_job.rv(),
                                                                 gather_fallback_fas_job.rv())
    output_map = combine_split_job.rv()

    # convert to stable coordinates.  if already stable coordinates (conversion done in graphmap), then run filter to make
    # sure no missing targets
    if options.fastaHeaderTable:
        stablefy_job = combine_split_job.addFollowOnJobFn(stablefy, gfa_id, output_map, header_table_id, not options.post_stable)
        output_map = stablefy_job.rv()
        
    # return all the files, as well as the 2 split logs
    return (seqIDMap, output_map, split_gfa_job.rv(1), split_fallback_gfa_job.rv(1))

def get_mask_bed(job, seq_id_map, min_length):
    """ make a bed file from the fastas """
    beds = []
    for event in seq_id_map.keys():
        fa_path, fa_id = seq_id_map[event]
        beds.append(job.addChildJobFn(get_mask_bed_from_fasta, event, fa_id, fa_path, min_length, disk=fa_id.size * 5).rv())
    return job.addFollowOnJobFn(cat_beds, beds).rv()

def cat_beds(job, bed_ids):
    in_beds = [job.fileStore.readGlobalFile(bed_id) for bed_id in bed_ids]
    out_bed = job.fileStore.getLocalTempFile()
    catFiles(in_beds, out_bed)
    return job.fileStore.writeGlobalFile(out_bed)
    
def split_gfa(job, config, gfa_id, paf_ids, ref_contigs, other_contig, reference_event, mask_bed_id, fasta_header_id):
    """ Use rgfa-split to divide a GFA and PAF into chromosomes.  The GFA must be in minigraph RGFA output using
    the desired reference. """

    if not paf_ids:
        # we can bypass when, ex, doing second pass on ambiguous sequences but not are present
        return [None, None]

    if not gfa_id and not getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "remap", typeFn=bool, default=False):
        # also bypass if remapping is off in the config (we know it's the second pass because gfa_id is None)
        return [None, None]

    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    paf_path = os.path.join(work_dir, "mg.paf")
    out_prefix = os.path.join(work_dir, "split_")
    bed_path = os.path.join(work_dir, "mask.bed")
    log_path = os.path.join(work_dir, "split.log")
    fasta_header_path = os.path.join(work_dir, "fasta-headers.tsv")
    
    if (mask_bed_id):
        job.fileStore.readGlobalFile(mask_bed_id, bed_path)

    if gfa_id:
        job.fileStore.readGlobalFile(gfa_id, gfa_path)
        
    paf_paths = []
    for i, paf_id in enumerate(paf_ids):
        paf_paths.append('{}.{}'.format(paf_path, i) if len(paf_ids) > 1 else paf_path)
        job.fileStore.readGlobalFile(paf_id, paf_paths[-1])
    if len(paf_paths) > 1:
        catFiles(paf_paths, paf_path)

    if fasta_header_id:
        job.fileStore.readGlobalFile(fasta_header_id, fasta_header_path)
    
    # get the minigraph "virutal" assembly name
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    # and look up its unique id prefix.  this will be needed to pick its contigs out of the list
    mg_id = graph_event

    # get the specificity filters
    query_coverage = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "minQueryCoverage", default="0")
    small_query_coverage = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "minQuerySmallCoverage", default="0")
    small_coverage_threshold = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "minQuerySmallThreshold", default="0")
    query_uniqueness = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "minQueryUniqueness", default="0")
    max_gap = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "maxGap", default="0")
    amb_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "ambiguousName", default="_AMBIGUOUS_")

    cmd = ['rgfa-split',
           '-p', paf_path,
           '-b', out_prefix,
           '-n', query_coverage,
           '-N', small_query_coverage,
           '-T', small_coverage_threshold,
           '-Q', query_uniqueness,
           '-P', max_gap,
           '-a', amb_name,
           '-l', log_path]
    if gfa_id:
        cmd += ['-g', gfa_path, '-G']
    if other_contig:
        cmd += ['-o', other_contig]
    if reference_event:
        cmd += ['-r', 'id={}|'.format(reference_event)]
    if mask_bed_id:
        cmd += ['-B', bed_path]
    min_mapq = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minMAPQ")
    if min_mapq:
        cmd += ['-A', min_mapq]
    if fasta_header_id:
        cmd += ['-L', fasta_header_path]
    # optional stuff added to second pass:
    if not gfa_id:
        remap_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "remapSplitOptions", default=None)
        if remap_opts:
            cmd += remap_opts.split(' ')        
    for contig in ref_contigs:
        cmd += ['-c', contig]

    cactus_call(parameters=cmd, work_dir=work_dir)

    output_id_map = {}
    for out_name in os.listdir(work_dir):
        file_name, ext = os.path.splitext(out_name)
        if file_name.startswith(os.path.basename(out_prefix)) and ext in [".gfa", ".paf", ".fa_contigs"] and \
           os.path.isfile(os.path.join(work_dir, file_name + ".fa_contigs")):
            name = file_name[len(os.path.basename(out_prefix)):]
            if name not in output_id_map:
                output_id_map[name] = {}
            if ext == '.paf':
                # apply the hacky naming correction so that subpaths have no special characterse in the hal (to make hubs happy)
                # this gets undone by hal2vg
                cactus_call(parameters=['sed', '-i', '-e', 's/\([^:]*\):\([0-9]*\)-\([0-9]*\)/echo "\\1_sub_$((\\2-1))_\\3"/e',
                                        '-e', 's/ /\t/g', os.path.join(work_dir, out_name)]) 
            output_id_map[name][ext[1:]] = job.fileStore.writeGlobalFile(os.path.join(work_dir, out_name))
            
    return output_id_map, job.fileStore.writeGlobalFile(log_path)

def split_fas(job, seq_id_map, split_id_map):
    """ Use samtools to split a bunch of fasta files into reference contigs, using the output of rgfa-split as a guide"""

    if (seq_id_map, split_id_map) == (None, None):
        return None

    root_job = Job()
    job.addChild(root_job)

    # map event name to dict of contgs.  ex fa_contigs["CHM13"]["chr13"] = file_id
    fa_contigs = {}
    # we do each fasta in parallel
    for event in seq_id_map.keys():
        fa_path, fa_id = seq_id_map[event]
        if fa_id.size:
            fa_contigs[event] = root_job.addChildJobFn(split_fa_into_contigs, event, fa_id, fa_path, split_id_map,
                                                       disk=fa_id.size * 3).rv()

    return fa_contigs

def split_fa_into_contigs(job, event, fa_id, fa_path, split_id_map):
    """ Use samtools turn on fasta into one for each contig. this relies on the informatino in .fa_contigs
    files made by rgfa-split """

    # download the fasta
    work_dir = job.fileStore.getLocalTempDir()
    fa_path = os.path.join(work_dir, os.path.basename(fa_path))
    is_gz = fa_path.endswith(".gz")
    job.fileStore.readGlobalFile(fa_id, fa_path, mutable=is_gz)
    if is_gz:
        # samtools can only work on bgzipped files.  so we uncompress here to be able to support gzipped too
        cactus_call(parameters=['gzip', '-fd', fa_path])
        fa_path = fa_path[:-3]

    unique_id = 'id={}|'.format(event)
                        
    contig_fa_dict = {}

    for ref_contig in split_id_map.keys():
        query_contig_list_id = split_id_map[ref_contig]['fa_contigs']
        list_path = os.path.join(work_dir, '{}.fa_contigs'.format(ref_contig))
        job.fileStore.readGlobalFile(query_contig_list_id, list_path)
        faidx_input_path = os.path.join(work_dir, '{}.fa_contigs.clean'.format(ref_contig))
        contig_count = 0
        with open(list_path, 'r') as list_file, open(faidx_input_path, 'w') as clean_file:
            for line in list_file:
                query_contig = line.strip()
                if query_contig.startswith(unique_id):
                    assert query_contig.startswith(unique_id)
                    query_contig = query_contig[len(unique_id):]
                    clean_file.write('{}\n'.format(query_contig))
                    contig_count += 1
        contig_fasta_path = os.path.join(work_dir, '{}_{}.fa'.format(event, ref_contig))
        if contig_count > 0:
            cmd = [['samtools', 'faidx', fa_path, '--region-file', faidx_input_path]]
            # transform chr1:10-15 (1-based inclusive) into chr1_sub_9_15 (0-based end open)
            # this is a format that contains no special characters in order to make assembly hubs
            # happy.  But it does require conversion going into vg which wants chr[9-15] and
            # hal2vg can get updated to do this autmatically
            cmd.append(['sed', '-e', 's/\([^:]*\):\([0-9]*\)-\([0-9]*\)/echo "\\1_sub_$((\\2-1))_\\3"/e'])
            if is_gz:
                cmd.append(['gzip'])
            cactus_call(parameters=cmd, outfile=contig_fasta_path)
        else:
            # TODO: review how cases like this are handled
            with open(contig_fasta_path, 'w') as empty_file:
                empty_file.write("")
        contig_fa_dict[ref_contig] = job.fileStore.writeGlobalFile(contig_fasta_path)

    return contig_fa_dict

def gather_fas(job, seq_id_map, output_id_map, contig_fa_map):
    """ take the split_fas output which has everything sorted by event, and move into the ref-contig-based table
    from split_gfa.  return the updated table, which can then be exported into the chromosome projects """

    if not contig_fa_map:
        return None

    for ref_contig in output_id_map.keys():
        output_id_map[ref_contig]['fa'] = {}
        for event, fa_id in contig_fa_map.items():
            if ref_contig in fa_id:
                output_id_map[ref_contig]['fa'][event] = fa_id[ref_contig]

    return output_id_map

def split_minimap_fallback(job, options, config, seqIDMap, output_id_map):
    """ take the output table from gather_fas, pull out the ambiguous sequences, remap them to the reference, and 
    add them to the events where possible"""

    # can't do anything without a reference
    if not options.reference:
        logger.info("Skipping minimap2 fallback as --reference was not specified")
        return None, None
    # todo: also skip if no ambgious sequences
    
    ref_path, ref_id = seqIDMap[options.reference]
    mm_mem = ref_id.size * 5
    if seqIDMap[options.reference][0].endswith('.gz'):
        mm_mem *= 4
    mm_index_job = job.addChildJobFn(minimap_index, ref_path, ref_id, disk=ref_id.size * 5, memory=mm_mem)
    mm_map_root_job = Job()
    mm_index_job.addFollowOn(mm_map_root_job)
    
    amb_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "ambiguousName", default="_AMBIGUOUS_")

    if amb_name not in output_id_map:
        logger.info("Skipping minmap2 fallback as no ambigious sequences found")
        return None, None

    # map every ambgiuous sequence against the reference in parallel
    paf_ids = []
    ambiguous_seq_id_map = {}
    for event, fa_id in output_id_map[amb_name]['fa'].items():
        paf_job = mm_map_root_job.addChildJobFn(minimap_map, config, mm_index_job.rv(), event, fa_id, seqIDMap[event][0],
                                                disk=ref_id.size * 3, memory=mm_mem)
        paf_ids.append(paf_job.rv())
        ambiguous_seq_id_map[event] = (seqIDMap[event][0], fa_id)

    return paf_ids, ambiguous_seq_id_map

def minimap_index(job, ref_name, ref_id):
    """ make a minimap2 index of a reference genome """

    work_dir = job.fileStore.getLocalTempDir()
    fa_path = os.path.join(work_dir, os.path.basename(ref_name))
    idx_path = fa_path + ".idx"
    job.fileStore.readGlobalFile(ref_id, fa_path)

    cactus_call(parameters=['minimap2', fa_path, '-d', idx_path, '-x', 'asm5'])

    return job.fileStore.writeGlobalFile(idx_path)

def minimap_map(job, config, minimap_index_id, event, fa_id, fa_name):
    """ run minimap2 """
    work_dir = job.fileStore.getLocalTempDir()
    idx_path = os.path.join(work_dir, "minmap2.idx")
    fa_path = os.path.join(work_dir, "ambiguous_" + os.path.basename(fa_name))
    job.fileStore.readGlobalFile(minimap_index_id, idx_path)
    job.fileStore.readGlobalFile(fa_id, fa_path)
    paf_path = fa_path + ".paf"

    # call minimap2 and stick our unique identifiers on the output right away to be consistent
    # with cactus-graphmap's paf output
    cmd = [['minimap2', idx_path, fa_path, '-c', '-x', 'asm5', '--secondary=no'],
           ['awk', 'BEGIN {{OFS="\t"}} $1="id={}|"$1'.format(event)]]
    
    min_mapq = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minMAPQ")
    if min_mapq:
        # add mapq filter used in mzgaf2paf
        cmd.append(['awk', '$12>={}'.format(min_mapq)])

    cactus_call(parameters=cmd, outfile=paf_path)

    return job.fileStore.writeGlobalFile(paf_path)    
    
def combine_splits(job, options, config, seq_id_map, original_id_map, remap_id_map):
    """ combine the output of two runs of gather_fas.  the first is the contigs determined by minigraph,
    the second from remapping the ambigious contigs with minimap2 """    
    root_job = Job()
    job.addChild(root_job)

    # no ambiguous remappings, nothing to do
    if not remap_id_map or len(remap_id_map) == 0:
        return original_id_map

    amb_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "ambiguousName", default="_AMBIGUOUS_")
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")

    # we're overwriting this, but need it below
    orig_amb_entry = original_id_map[amb_name]
    
    # note: we're not handling case where 100% of a given reference contigs are ambiguous
    for ref_contig in original_id_map:
        if ref_contig == amb_name:
            # for ambiguous sequence, we overwrite and don't combine
            if ref_contig in remap_id_map:
                original_id_map[ref_contig] = remap_id_map[ref_contig]
            else:
                original_id_map[ref_contig] = None
        elif ref_contig in remap_id_map:            
            total_size = 0
            for event in original_id_map[ref_contig]['fa']:
                total_size += original_id_map[ref_contig]['fa'][event].size
                if event in remap_id_map[ref_contig]['fa']:
                    total_size += remap_id_map[ref_contig]['fa'][event].size
            original_id_map[ref_contig] = root_job.addChildJobFn(combine_ref_contig_splits,
                                                                 original_id_map[ref_contig],
                                                                 remap_id_map[ref_contig],
                                                                 disk=total_size * 4).rv()

    # the remap/split logic without this doesn't work with stable coordintes
    # and even then results are flaky:  much better to run with post_stable == True (pass in PAFs with minigraph targets)
    # the convert at the very end 
    if options.fastaHeaderTable and not options.post_stable:
        raise RuntimeError("Splitting pafs with stable coordinates (generated with cactus-graphmap --fastaHeaderTable)"
                           " requires useMinimapPaf to be set to \"1\" in the config.  But even then, this logic has not "
                           "yet been fully tested")
        
    return root_job.addFollowOnJobFn(combine_paf_splits, options, config, seq_id_map, original_id_map, orig_amb_entry,
                                     remap_id_map, amb_name, graph_event).rv()

def combine_ref_contig_splits(job, original_ref_entry, remap_ref_entry):
    """ combine fa and paf files for splits for a ref contig """
    work_dir = job.fileStore.getLocalTempDir()

    # combine the event fas
    for event in original_ref_entry['fa']:
        if event in remap_ref_entry['fa']:
            orig_fa_path = os.path.join(work_dir, event + '.fa')
            remap_fa_path = os.path.join(work_dir, event + '.remap.fa')
            new_fa_path = os.path.join(work_dir, event + '.combine.fa')
            job.fileStore.readGlobalFile(original_ref_entry['fa'][event], orig_fa_path, mutable=True)
            job.fileStore.readGlobalFile(remap_ref_entry['fa'][event], remap_fa_path, mutable=True)
            catFiles([orig_fa_path, remap_fa_path], new_fa_path)
            original_ref_entry['fa'][event] = job.fileStore.writeGlobalFile(new_fa_path)
            os.remove(orig_fa_path)
            os.remove(remap_fa_path)
                
    return original_ref_entry

def combine_paf_splits(job, options, config, seq_id_map, original_id_map, orig_amb_entry,
                       remap_id_map, amb_name, graph_event):
    """ pull out PAF entries for contigs that were ambiguous in the first round but assigned by minimap2
    then add them to the chromosome PAFs     
    """

    if amb_name not in original_id_map:
        return original_id_map

    work_dir = job.fileStore.getLocalTempDir()
    amb_paf_path = os.path.join(work_dir, 'amb.paf')
    job.fileStore.readGlobalFile(orig_amb_entry['paf'], amb_paf_path, mutable=True)

    # use_minimap_paf = True: return the minimap2 mappings for ambiguous contigs in final output
    # use_minimap_paf = False: ambiguous contigs are assigned to chromosomes base on minimap2, but their minigraph 
    #                          alignments are returned in the final paf"""
    use_minimap_paf = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "useMinimapPAF",
                                        typeFn=bool, default=False)

    # it's simpler not to support both codepaths right now.  the main issue is that -u can cause contigs to be split
    # in which case they get renamed, so pulling them in from the existing PAF would require a pass to resolove all the
    # offsets
    if not use_minimap_paf and '-u' in getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "remapSplitOptions",
                                                         default=""):
        raise RuntimeError("useMinimapPAF must be set when -u present in remapSplitOptions")

    for ref_contig in remap_id_map.keys():
        if ref_contig != amb_name and ref_contig in original_id_map:

            # make a set of all minigraph nodes in this contig
            mg_contig_set = set()
            if not use_minimap_paf:
                mg_fa_path = os.path.join(work_dir, '{}.{}.fa'.format(graph_event, ref_contig))
                if seq_id_map[graph_event][0].endswith('.gz'):
                    mg_fa_path += '.gz'
                mg_contigs_path = os.path.join(work_dir, '{}.contigs'.format(graph_event))
                job.fileStore.readGlobalFile(original_id_map[ref_contig]['fa'][graph_event], mg_fa_path, mutable=True)
                cactus_call(parameters=[['zcat' if mg_fa_path.endswith('.gz') else 'cat', mg_fa_path],
                                        ['grep', '>'], ['cut', '-c', '2-']], outfile=mg_contigs_path)
                with open(mg_contigs_path, 'r') as mg_contigs_file:
                    for line in mg_contigs_file:
                        mg_contig_set.add('id={}|{}'.format(graph_event, line.strip()))
                os.remove(mg_fa_path)
                os.remove(mg_contigs_path)

            #make a set of all the query contigs that we want to remove from ambiguous and add to this contig
            query_contig_set = set()

            for event in remap_id_map[ref_contig]['fa']:
                if event != graph_event and remap_id_map[ref_contig]['fa'][event].size > 0:
                    # read the contigs assigned to this sample for this chromosome by scanning fasta headers
                    tmp_fa_path = os.path.join(work_dir, 'tmp.fa')
                    if seq_id_map[event][0].endswith('.gz'):
                        tmp_fa_path += '.gz'
                    if os.path.isfile(tmp_fa_path):
                        os.remove(tmp_fa_path)
                    job.fileStore.readGlobalFile(remap_id_map[ref_contig]['fa'][event], tmp_fa_path, mutable=True)
                    contigs_path = os.path.join(work_dir, '{}.contigs'.format(event))
                    cactus_call(parameters=[['zcat' if tmp_fa_path.endswith('.gz') else 'cat', tmp_fa_path],
                                            ['grep', '>'], ['cut', '-c', '2-']], outfile=contigs_path)
                    # add them to the grep
                    with open(contigs_path, 'r') as contigs_file:
                        for line in contigs_file:
                            query_contig_set.add('id={}|{}'.format(event, line.strip()))

            if query_contig_set:
                # pull out remapped contigs into this path
                new_contig_path = os.path.join(work_dir, '{}.remap.paf'.format(ref_contig))
                do_append = False
                if ref_contig in original_id_map and 'paf' in original_id_map[ref_contig]:
                    job.fileStore.readGlobalFile(original_id_map[ref_contig]['paf'], new_contig_path, mutable=True)
                    do_append = True
                    
                # make an updated ambiguous paf with the contigs removed in this path                
                temp_contig_path = os.path.join(work_dir, amb_paf_path + '.temp.remove')                    
                with open(new_contig_path, 'a' if do_append else 'w') as new_contig_file, \
                     open(amb_paf_path, 'r') as amb_paf_file, \
                     open(temp_contig_path, 'w') as temp_contig_file:
                    # scan the ambgiuous paf from minigraph
                    for line in amb_paf_file:
                        toks = line.split('\t')
                        if len(toks) > 5 and toks[0] in query_contig_set:
                            if toks[5] in mg_contig_set and not use_minimap_paf:
                                # move the contig if both the query and target belong to reference contig
                                new_contig_file.write(line)
                        else:
                            # leave the contig in ambiguous
                            temp_contig_file.write(line)
                    if use_minimap_paf:
                        # if we're taking the contigs from minigraph, append them here (as they weren't added in
                        # the loop above)
                        minimap_paf_path = os.path.join(work_dir, '{}.minimap.paf'.format(ref_contig))
                        job.fileStore.readGlobalFile(remap_id_map[ref_contig]['paf'], minimap_paf_path)
                        with open(minimap_paf_path, 'r') as minimap_paf_file:
                            for line in minimap_paf_file:
                                toks = line.split('\t')
                                if len(toks) > 5:
                                    toks[5] = 'id={}|{}'.format(options.reference, toks[5])
                                new_contig_file.write('\t'.join(toks))
                        
                # update the map
                original_id_map[ref_contig]['paf'] = job.fileStore.writeGlobalFile(new_contig_path)
                # update the ambigious paf
                cactus_call(parameters=['mv', temp_contig_path, amb_paf_path])

    # update the ambiguous paf
    if amb_name in original_id_map and original_id_map[amb_name]:
        original_id_map[amb_name]['paf'] = job.fileStore.writeGlobalFile(amb_paf_path)
    else:
        assert os.path.getsize(amb_paf_path) == 0
    
    return original_id_map

def stablefy(job, gfa_id, id_map, header_table_id, filter_only):
    """ convert all minigraph node targets into stable query coordinates using the table, 
    then run filter to make sure that all targets are present in queries"""
    for ref_contig in id_map:
        if 'paf' in id_map[ref_contig]:
            id_map[ref_contig]['paf'] = job.addChildJobFn(stablefy_and_filter_paf, gfa_id, id_map[ref_contig]['paf'],
                                                          header_table_id, filter_only,
                                                          disk = gfa_id.size + 3 * id_map[ref_contig]['paf'].size).rv()

    return id_map


def stablefy_and_filter_paf(job, gfa_id, paf_id, header_table_id, filter_only):
    """
    use the table to convert target contigs from minigraph node names to query contig (stable) space. 
    then run filter:
    with stable coordinates, we can have cases where a minigraph node translates to a reference path that itself
    is not assigned to the given contig.  this is a rare case where maybe minigraph maps the same contig to 2 chromosomes
    or maybe it only maps a little piece of the contig.  whichever the case, it will cause errors downstream, so we must remove

    todo: it's bit wasteful to reparse the gfa each time here, but won't affect bottom line much
    
    """
    work_dir = job.fileStore.getLocalTempDir()
    paf_path = os.path.join(work_dir, "contig.paf")
    job.fileStore.readGlobalFile(paf_id, paf_path)
    
    if not filter_only:
        header_table_path = os.path.join(work_dir, "header-table.tsv")
        gfa_path = os.path.join(work_dir, "graph.gfa")
        job.fileStore.readGlobalFile(header_table_id, header_table_path)
        job.fileStore.readGlobalFile(gfa_id, gfa_path)
        stable_paf_path = paf_path + ".stable"                
        cactus_call(parameters=['paf2stable', gfa_path, header_table_path, paf_path], outfile = stable_paf_path)
        return job.fileStore.writeGlobalFile(stable_paf_path)
            
    fixed_path = paf_path + '.remove_unmapped_targets'
    query_set = set()
    with open(paf_path, 'r') as paf_file:
        for line in paf_file:
            toks = line.split('\t')
            if len(toks) > 11:
                query_set.add(toks[0])
    with open(paf_path, 'r') as paf_file, open(fixed_path, 'w') as fixed_file:
        for line in paf_file:
            toks = line.split('\t')
            if len(toks) <= 11 or toks[5] in query_set:
                fixed_file.write(line)
    return job.fileStore.writeGlobalFile(fixed_path)


def export_split_data(toil, input_seq_id_map, output_id_map, split_log_ids, output_dir, config, filter_graph_event):
    """ download all the split data locally """

    amb_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "ambiguousName", default="_AMBIGUOUS_")
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    
    chrom_file_map = {}
    
    for ref_contig in output_id_map.keys():
        if output_id_map[ref_contig] is None:
            # todo: check ambigous?
            continue
        ref_contig_path = os.path.join(output_dir, ref_contig)
        if not os.path.isdir(ref_contig_path) and not ref_contig_path.startswith('s3://'):
            os.makedirs(ref_contig_path)

        # GFA: <output_dir>/<contig>/<contig>.gfa
        if 'gfa' in output_id_map[ref_contig]:
            # we do this check because no gfa made for ambiguous sequences "contig"
            toil.exportFile(output_id_map[ref_contig]['gfa'], makeURL(os.path.join(ref_contig_path, '{}.gfa'.format(ref_contig))))

        # PAF: <output_dir>/<contig>/<contig>.paf
        paf_path = os.path.join(ref_contig_path, '{}.paf'.format(ref_contig))
        if 'paf' in output_id_map[ref_contig]:
            # rgfa-split doesn't write empty pafs: todo: should it?
            toil.exportFile(output_id_map[ref_contig]['paf'], makeURL(paf_path))

        # Fasta: <output_dir>/<contig>/fasta/<event>_<contig>.fa ..
        seq_file_map = {}
        for event, ref_contig_fa_id in output_id_map[ref_contig]['fa'].items():
            fa_base = os.path.join(ref_contig_path, 'fasta')
            if not os.path.isdir(fa_base) and not fa_base.startswith('s3://'):
                os.makedirs(fa_base)
            fa_path = makeURL(os.path.join(fa_base, '{}_{}.fa'.format(event, ref_contig)))
            if input_seq_id_map[event][0].endswith('.gz'):
                fa_path += '.gz'
            seq_file_map[event] = fa_path
            toil.exportFile(ref_contig_fa_id, fa_path)

        # Seqfile: <output_dir>/seqfiles/<contig>.seqfile
        seq_file_path = os.path.join(output_dir, 'seqfiles', '{}.seqfile'.format(ref_contig))
        if seq_file_path.startswith('s3://'):
            seq_file_temp_path = getTempFile()
        else:
            seq_file_temp_path = seq_file_path
            if not os.path.isdir(os.path.dirname(seq_file_path)):
                os.makedirs(os.path.dirname(seq_file_path))
        with open(seq_file_temp_path, 'w') as seq_file:
            for event, fa_path in seq_file_map.items():
                # cactus can't handle empty fastas.  if there are no sequences for a sample for this
                # contig, just don't add it.
                if output_id_map[ref_contig]['fa'][event].size > 0:
                    seq_file.write('{}\t{}\n'.format(event, fa_path))
        if seq_file_path.startswith('s3://'):
            write_s3(seq_file_temp_path, seq_file_path)

        # Top-level seqfile
        chrom_file_map[ref_contig] = seq_file_path, paf_path

    # Chromfile : <coutput_dir>/chromfile.txt
    chrom_file_path = os.path.join(output_dir, 'chromfile.txt')
    if chrom_file_path.startswith('s3://'):
        chrom_file_temp_path = getTempFile()
    else:
        chrom_file_temp_path = chrom_file_path        
    with open(chrom_file_temp_path, 'w') as chromfile:
        for ref_contig, seqfile_paf in chrom_file_map.items():
            if ref_contig != amb_name:
                seqfile, paf = seqfile_paf[0], seqfile_paf[1]
                if seqfile.startswith('s3://'):
                    # no use to have absolute s3 reference as cactus-align requires seqfiles passed locally
                    seqfile = 'seqfiles/{}'.format(os.path.basename(seqfile))
                chromfile.write('{}\t{}\t{}\n'.format(ref_contig, seqfile, paf))
    if chrom_file_path.startswith('s3://'):
        write_s3(chrom_file_temp_path, chrom_file_path)

    toil.exportFile(split_log_ids[0], makeURL(os.path.join(output_dir, 'minigraph.split.log')))
    if split_log_ids[1]:
        toil.exportFile(split_log_ids[1], makeURL(os.path.join(output_dir, 'minimap2.ambiguous.split.log')))
        
if __name__ == "__main__":
    main()
