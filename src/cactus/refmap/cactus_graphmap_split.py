#!/usr/bin/env python3

"""This sciprt will take the cactus-graphmap output and input and split it up into chromosomes.  This is done after a whole-genome graphmap as we need the whole-genome minigraph alignments to do the chromosome splitting in the first place.  For each chromosome, it will make a PAF, GFA and Seqfile (pointing to chrosome fastas)
"""
import os, sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit
import shutil
import re
import gzip

from array import array
from collections import defaultdict
from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import setupBinaries, importSingularityImage, cactus_walltime
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles, pangenome_file_name, pangenome_seqfile_name
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options, add_cactus_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import unzip_gz, write_s3
from cactus.shared.common import get_faidx_subpath_rename_cmd
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.common import clean_jobstore_files
from cactus.shared.version import cactus_commit
from cactus.preprocessor.fileMasking import get_mask_bed_from_fasta
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers
from cactus.refmap.cactus_graphmap import filter_paf, apply_mgsplit_filter_overrides
from cactus.refmap.cactus_graphmap import FILTER_PAF_SECS_PER_GB
from cactus.refmap.cactus_minigraph import check_sample_names, minigraph_gfa_from_pansn, minigraph_gfa_to_pansn
from cactus.refmap.cactus_minigraph import bgzip_gfa_rename
from cactus.refmap.cactus_minigraph import GFA_RENAME_SECS_PER_GB, RAW_BYTES_PER_GZ_BYTE
from cactus.preprocessor.maskingCommon import longest_record_size
from toil.job import Job
from toil.lib.conversions import bytes2human
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory, getTempFile, catFiles

def main():
    parser = Job.Runner.getDefaultArgumentParser()
    add_cactus_toil_options(parser)

    parser.add_argument("seqFile", help = "Seq file (gzipped fastas supported)")
    parser.add_argument("minigraphGFA", help = "Minigraph-compatible reference graph in GFA format (can be gzipped)")
    parser.add_argument("graphmapPAF", type=str, help = "Output pairwise alignment file in PAF format (can be gzipped)")
    parser.add_argument("--outDir", required=True, type=str, help = "Output directory")
    parser.add_argument("--refContigs", nargs="*", help = "Subset to these reference contigs (multiple allowed)", default=[])
    parser.add_argument("--refContigsFile", type=str, help = "Subset to (newline-separated) reference contigs in this file")
    parser.add_argument("--otherContig", type=str, help = "Lump all reference contigs unselected by above options into single one with this name")
    parser.add_argument("--reference", required=True, nargs='+', type=str, help = "Name of reference (in seqFile).  Ambiguity filters will not be applied to it")
    parser.add_argument("--maskFilter", type=int, help = "Ignore softmasked sequence intervals > Nbp")
    parser.add_argument("--minIdentity", type=float, help = "Ignore PAF lines with identity (column 10/11) < this (overrides minIdentity in <graphmap> in config)")
    parser.add_argument("--permissiveContigFilter", nargs='?', const='0.25', default=None, type=float, help = "If specified, override the configuration to accept contigs so long as they have at least given fraction of coverage (0.25 if no fraction specified). This can increase sensitivity of very small, fragmented and/or diverse assemblies.")
    parser.add_argument("--mgSplit", action="store_true", default=False,
                        help="Use this when the input PAF was produced against a reference-only minigraph (cactus-pangenome --mgSplit / cactus-minigraph --refOnly): relaxes the block-length and overlap filters so palindromic mappings don't get axed before rgfa-split.")

    parser.add_argument("--collapse", help = "Incorporate minimap2 self-alignments.", action='store_true', default=False)
    
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

    # support but ignore multi reference
    if options.reference:
        options.reference = options.reference[0]

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()
    cactus_graphmap_split(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-graphmap-split has finished after {} seconds".format(run_time))

def cactus_graphmap_split(options):
    with Toil(options) as toil:
        importSingularityImage(options)

        #load cactus config
        config_node = ET.parse(options.configFile).getroot()
        config = ConfigWrapper(config_node)
        config.substituteAllPredefinedConstantsWithLiterals(options)

        if options.collapse:
            findRequiredNode(config_node, "graphmap").attrib["collapse"] = 'all'

        if options.mgSplit:
            apply_mgsplit_filter_overrides(config_node)

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
            options.refContigs = list(ref_contigs)

            if options.refContigs:
                max_ref_contigs = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "maxRefContigs", typeFn=int, default=128)
                if len(options.refContigs) > max_ref_contigs:
                    logger.warning('You specified {} refContigs, which is greater than the suggested limit of {}. This may cause issues downstream.'.format(len(options.refContigs), max_ref_contigs))
                

            if options.otherContig:
                assert options.otherContig not in options.refContigs

            # get the minigraph "virutal" assembly name
            graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")

            # load the seqfile
            seqFile = SeqFile(options.seqFile, defaultBranchLen=config.getDefaultBranchLen(pangenome=True), pangenome=True,
                              dot_haplotypes=config.getDotHaplotypes())

            #import the graph
            gfa_id = toil.importFile(makeURL(options.minigraphGFA))

            #import the paf
            paf_id = toil.importFile(makeURL(options.graphmapPAF))

            #import the sequences (that we need to align for the given event, ie leaves and outgroups)
            input_seq_id_map = {}
            input_name_map = {}
            leaves = set([seqFile.tree.getName(node) for node in seqFile.tree.getLeaves()])
            
            if graph_event not in leaves:
                raise RuntimeError("Minigraph name {} not found in seqfile".format(graph_event))

            # validate the sample names
            options.reference = check_sample_names(seqFile.pathMap.keys(), options.reference, config.getDotHaplotypes())
            if options.reference and options.reference not in leaves:
                raise RuntimeError("Name given with --reference {} not found in seqfile".format(options.reference))
                        
            for genome, seq in seqFile.pathMap.items():
                if genome in leaves:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in sorted(os.listdir(seq))], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    input_seq_id_map[genome] = toil.importFile(seq)
                    input_name_map[genome] = os.path.basename(seq)

            # run the workflow
            wf_output = toil.start(Job.wrapJobFn(graphmap_split_workflow, options, config, input_seq_id_map, input_name_map,
                                                 gfa_id, options.minigraphGFA,
                                                 paf_id, options.graphmapPAF, walltime=cactus_walltime()))

        #export the split data
        export_split_data(toil, wf_output[0], wf_output[1], wf_output[2], wf_output[3], options.outDir, config)

# Seconds per GB of GFA-plus-PAF for rgfa-split: 1220s to split the 8.3 GiB GFA and 34.4 GiB PAF
# of the HPRC v2.0 whole panel, ie ~27 s/GB, with the per-chromosome sed renames that follow adding
# ~130s more.  This is the heaviest job in cactus-graphmap-split, and separate_ref_contigs runs the
# same binary on the same kind of input.
RGFA_SPLIT_SECS_PER_GB = 32

def graphmap_split_workflow(job, options, config, seq_id_map, seq_name_map, gfa_id, gfa_path, paf_id, paf_path, sanitize=True,
                            pansn_gfa_input=True):

    root_job = Job(walltime=cactus_walltime())
    job.addChild(root_job)

    # can be a list coming in from cactus-pangenome, but we only need first item
    if type(options.reference) is list:
        options.reference = options.reference[0]

    #override the minIdentity
    if options.minIdentity is not None:
        findRequiredNode(config.xmlRoot, "graphmap").attrib["minIdentity"] = str(options.minIdentity)

    apply_permissive_contig_filter(config, options.permissiveContigFilter)

    # get the sizes before we overwrite below
    gfa_size = gfa_id.size
    paf_size = paf_id.size
    genome_names = set(seq_id_map.keys())    

    # fix up the headers
    if sanitize:
        sanitize_job = root_job.addChildJobFn(sanitize_fasta_headers, seq_id_map, pangenome=True, walltime=cactus_walltime())
        seq_id_map = sanitize_job.rv()
    else:
        sanitize_job = Job(walltime=cactus_walltime())
        root_job.addChild(sanitize_job)

    # auto-set --refContigs
    if not options.refContigs:
        # samtools faidx of a whole reference is 21s measured; the rest is downloading the fasta,
        # whose size is not in scope here (seq_id_map is a promise), so it goes in the constant --
        # 600s is a ~3 GB reference downloaded on a badly contended filesystem, and a bigger one
        # still fits inside the safety factor
        refcontig_job = sanitize_job.addFollowOnJobFn(detect_ref_contigs, config, options, seq_id_map,
                                                     walltime=cactus_walltime(600))
        ref_contigs = refcontig_job.rv()
        options.otherContig = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "otherContigName", typeFn=str, default="chrOther")
        sanitize_job = refcontig_job
    else:
        # move from options to avoid promise confusion
        ref_contigs = options.refContigs

    # convert the GFA from PanSN to Cactus names
    if pansn_gfa_input:
        # the renaming pass decompresses the GFA before bgzipping it back up, so it needs room for
        # the raw copy (reckoned at 10x, as below) on top of the compressed input and output
        rename_gfa_job = root_job.addChildJobFn(minigraph_gfa_from_pansn, genome_names, gfa_path, gfa_id,
                                                disk=gfa_size*12,
                                                walltime=cactus_walltime(GFA_RENAME_SECS_PER_GB * gfa_size / 1e9,
                                                                         io_bytes=2*gfa_size))
        new_root_job = Job(walltime=cactus_walltime())
        root_job.addFollowOn(new_root_job)
        root_job = new_root_job
        gfa_id = rename_gfa_job.rv(0)
            
    # use file extension to sniff out compressed input
    if gfa_path.endswith(".gz"):
        # gunzip is fast (27s for the 0.83 GiB compressed whole-panel GFA); the cost is writing the
        # ~10x bigger raw file back to the jobstore, which is what io_bytes is counting
        gfa_id = root_job.addChildJobFn(unzip_gz, gfa_path, gfa_id, delete_original=False, disk=gfa_size * 10,
                                        walltime=cactus_walltime(60, io_bytes=(1 + RAW_BYTES_PER_GZ_BYTE) * gfa_size)).rv()
        gfa_size *= 10
    if paf_path.endswith(".gz"):
        # the same, but on a PAF that reaches 34 GiB raw at panel scale
        paf_id = root_job.addChildJobFn(unzip_gz, paf_path, paf_id, delete_original=False, disk=paf_id.size * 10,
                                        walltime=cactus_walltime(120, io_bytes=(1 + RAW_BYTES_PER_GZ_BYTE) * paf_size)).rv()
        paf_size *= 10

    # do some basic paf filtering
    paf_filter_mem = max(paf_id.size * 10, 2**32)
    paf_filter_job = root_job.addFollowOnJobFn(filter_paf, paf_id, config, reference=options.reference,
                                               disk = paf_id.size * 10, memory=cactus_clamp_memory(paf_filter_mem),
                                               walltime=cactus_walltime(60 + FILTER_PAF_SECS_PER_GB * paf_size / 1e9,
                                                                        io_bytes=2*paf_size))
    paf_id = paf_filter_job.rv()
    root_job = paf_filter_job

    mask_bed_id = None
    if options.maskFilter:
        mask_bed_id = sanitize_job.addFollowOnJobFn(get_mask_bed, seq_id_map, options.maskFilter, walltime=cactus_walltime()).rv()
        
    # use rgfa-split to split the gfa and paf up by contig
    split_gfa_job = root_job.addFollowOnJobFn(split_gfa, config, gfa_id, [paf_id], ref_contigs,
                                              options.otherContig, options.reference, mask_bed_id,
                                              disk=(gfa_size + paf_size) * 5,
                                              memory=cactus_clamp_memory((gfa_size + paf_size) * 3),
                                              walltime=cactus_walltime(120 + RGFA_SPLIT_SECS_PER_GB * (gfa_size + paf_size) / 1e9,
                                                                       io_bytes=2 * (gfa_size + paf_size)))

    # use the output of the above splitting to do the fasta splitting
    split_fas_job = split_gfa_job.addFollowOnJobFn(split_fas, seq_id_map, seq_name_map, split_gfa_job.rv(0), walltime=cactus_walltime())

    # gather everythign up into a table
    gather_fas_job = split_fas_job.addFollowOnJobFn(gather_fas, split_gfa_job.rv(0), split_fas_job.rv(0), split_fas_job.rv(1), walltime=cactus_walltime())

    # lump "other" contigs together into one file (to make fewer align jobs downstream)
    bin_other_job = gather_fas_job.addFollowOnJobFn(bin_other_contigs, config, ref_contigs, options.otherContig, gather_fas_job.rv(0),
                                                    disk=(gfa_size + paf_size) * 2,
                                                    walltime=cactus_walltime(0, io_bytes=(gfa_size + paf_size) * 2))

    # return all the files, as well as the 2 split logs
    return (seq_name_map, bin_other_job.rv(), split_gfa_job.rv(1), gather_fas_job.rv(1))

def detect_ref_contigs(job, config, options, seq_id_map):
    """ automatically determine --refContigs from the data """
    work_dir = job.fileStore.getLocalTempDir()
    # it will be unzipped since we're after sanitize
    fa_path = os.path.join(work_dir, options.reference + '.fa')
    job.fileStore.readGlobalFile(seq_id_map[options.reference], fa_path)
    cactus_call(parameters=['samtools', 'faidx', fa_path])
    contigs = []
    with open(fa_path + '.fai', 'r') as fai_file:
        for line in fai_file:
            toks = line.split('\t')
            assert len(toks) > 2
            # clean out prefix that was added by sanitize
            assert toks[0].startswith('id={}|'.format(options.reference))
            clean_name = toks[0][len('id={}|'.format(options.reference)):]
            contigs.append((clean_name, float(toks[1])))

    # get the cutoff heuristics 
    max_ref_contigs = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "maxRefContigs", typeFn=int, default=128)
    ref_contig_dropoff = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "refContigDropoff", typeFn=float, default=10.0)
    ref_contig_regex = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "refContigRegex", typeFn=str, default=None)
    other_contig_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "otherContigName", typeFn=str, default="chrOther")

    # sort by decreasing size
    sorted_contigs = sorted(contigs, key=lambda x : x[1], reverse=True)
    
    ref_contigs = set()
    # pass 1: do the regex matches
    if ref_contig_regex:
        for contig, length in contigs:
            if re.fullmatch(ref_contig_regex, contig):
                ref_contigs.add(contig)
                if len(ref_contigs) >= max_ref_contigs:
                    break
        
    # pass 2: do the dropoff
    longest_length = sorted_contigs[0][1]
    for contig_len in sorted_contigs:
        if len(ref_contigs) >= max_ref_contigs:
            break
        if contig_len[0] not in ref_contigs:
            dropoff = longest_length / contig_len[1]
            if dropoff >= ref_contig_dropoff:
                break
            ref_contigs.add(contig_len[0])

    ref_contigs = sorted(list(ref_contigs))

    msg = "auto-detected --refContigs {}".format(' '.join(ref_contigs))
    if len(ref_contigs) < len(sorted_contigs):
        msg += " --otherContig {}".format(other_contig_name)

    RealtimeLogger.info(msg)

    return sorted(list(ref_contigs))
    

def get_mask_bed(job, seq_id_map, min_length):
    """ make a bed file from the fastas """
    beds = []
    for event in seq_id_map.keys():
        fa_id = seq_id_map[event]
        fa_path = '{}.fa'.format(event) 
        # cactus_softmask2hardmask scanning a vertebrate genome is p99 51s / max 89s over the 625
        # runs of the VGP 577-way; the rest is the fasta in and the BED out
        beds.append(job.addChildJobFn(get_mask_bed_from_fasta, event, fa_id, fa_path, min_length, disk=fa_id.size * 5,
                                      walltime=cactus_walltime(120, io_bytes=2*fa_id.size)).rv())
    # one BED per genome -- order tens of GiB across a whole panel -- concatenated and written back.
    # bed_ids are promises and no fasta size is in scope either, so the whole cost has to be a
    # constant: 1200s is ~30 GiB of round-trip staging at a pessimistic 30 MB/s
    return job.addFollowOnJobFn(cat_beds, beds, walltime=cactus_walltime(1200)).rv()

def cat_beds(job, bed_ids):
    in_beds = [job.fileStore.readGlobalFile(bed_id) for bed_id in bed_ids]
    out_bed = job.fileStore.getLocalTempFile()
    catFiles(in_beds, out_bed)
    return job.fileStore.writeGlobalFile(out_bed)
    
def apply_permissive_contig_filter(config, permissive_contig_filter):
    """ override the contig filter parameters with a single value taken from --permissiveContigFilter,
    also disabling the uniqueness threshold.  every rgfa-split pass needs this applied, otherwise a
    later pass undoes with the strict defaults what the user asked to keep """
    if permissive_contig_filter is not None:
        findRequiredNode(config.xmlRoot, "graphmap_split").attrib["minQueryCoverages"] = str(permissive_contig_filter)
        findRequiredNode(config.xmlRoot, "graphmap_split").attrib["minQueryCoverageThresholds"] = ""
        findRequiredNode(config.xmlRoot, "graphmap_split").attrib["minQueryUniqueness"] = "1"

def get_split_specificity_opts(config):
    """ the rgfa-split coverage / uniqueness thresholds that decide which reference contig a query
    belongs to.  returns (coverage options, uniqueness, ambiguous bin name) """
    query_coverages = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "minQueryCoverages", default=None)
    query_coverage_thresholds = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "minQueryCoverageThresholds", default=None)
    coverage_opts = []
    if query_coverages is not None:
        try:
            query_cov_vals = [float(x) for x in query_coverages.split()]
            query_thresh_vals = [int(x) for x in query_coverage_thresholds.split()]
            assert len(query_cov_vals) == len(query_thresh_vals) + 1
            assert sorted(query_thresh_vals) == query_thresh_vals
            for cv in query_cov_vals:
                coverage_opts += ['-n', str(cv)]
            for tv in query_thresh_vals:
                coverage_opts += ['-T', str(tv)]
        except:
            raise RuntimeError("minQueryCoverages and / or minQueryCoverageThresholds malspecified in config")
    query_uniqueness = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "minQueryUniqueness", default="0")
    amb_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "ambiguousName", default="_AMBIGUOUS_")
    return coverage_opts, query_uniqueness, amb_name

def split_gfa(job, config, gfa_id, paf_ids, ref_contigs, other_contig, reference_event, mask_bed_id):
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

    # get the minigraph "virutal" assembly name
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    # and look up its unique id prefix.  this will be needed to pick its contigs out of the list
    mg_id = graph_event

    # get the specificity filters
    coverage_opts, query_uniqueness, amb_name = get_split_specificity_opts(config)

    cmd = ['rgfa-split',
           '-p', paf_path,
           '-b', out_prefix,
           '-Q', query_uniqueness,
           '-a', amb_name,
           '-L', log_path]
    cmd += coverage_opts    
    if gfa_id:
        cmd += ['-g', gfa_path, '-G']
    if reference_event:
        cmd += ['-r', 'id={}|'.format(reference_event)]
    if mask_bed_id:
        cmd += ['-B', bed_path]
    min_mapq = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minMAPQ")
    if min_mapq:
        cmd += ['-A', min_mapq]
    # optional stuff added to second pass:
    if not gfa_id:
        remap_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "remapSplitOptions", default=None)
        if remap_opts:
            cmd += remap_opts.split(' ')
    if not other_contig:
        # we use this option only when subsetting to given ref contigs *without* using the "other"
        # option.  otherwise, we rely on bin_other_contigs here in python to do the job
        for contig in ref_contigs:
            cmd += ['-c', contig]
            # hack in to make sure we deal with gfas with/without prexixes
            # todo: once all merged, this needs to be simplified and cleaned
            if reference_event and not contig.startswith('id={}|'.format(reference_event)):
                cmd += ['-c', 'id={}|{}'.format(reference_event, contig)]
            if contig.startswith('id=') and contig.find('|') > 3:
                cmd += ['-c', contig[contig.find('|')+1:]]            

    cactus_call(parameters=cmd, work_dir=work_dir, job_memory=job.memory)

    output_id_map = {}
    for out_name in sorted(os.listdir(work_dir)):
        file_name, ext = os.path.splitext(out_name)
        if file_name.startswith(os.path.basename(out_prefix)) and ext in [".gfa", ".paf", ".fa_contigs"] and \
           os.path.isfile(os.path.join(work_dir, file_name + ".fa_contigs")):
            name = file_name[len(os.path.basename(out_prefix)):]

            # don't leave unique identifier in names
            if name.startswith('id=') and name.find('|') > 3:
                name=name[name.find('|') + 1:]

            if name not in output_id_map:
                output_id_map[name] = {}
            if ext == '.paf':
                # spaces become tabs.  This used to be sed -i, with the subpath rename of
                # get_faidx_subpath_rename_cmd(), which never matched; GNU sed's -i and \t do not work on Mac.
                paf_path = os.path.join(work_dir, out_name)
                with open(paf_path, 'r') as paf_file, open(paf_path + '.tabs', 'w') as tab_file:
                    for line in paf_file:
                        tab_file.write(line.replace(' ', '\t'))
                os.replace(paf_path + '.tabs', paf_path)
            output_id_map[name][ext[1:]] = job.fileStore.writeGlobalFile(os.path.join(work_dir, out_name))
            
    return output_id_map, job.fileStore.writeGlobalFile(log_path)

def ref_contig_names(gfa_path):
    """ the distinct reference contigs in a minigraph GFA, ie its rank-0 (SR:i:0) sources """
    ref_contigs = set()
    with open(gfa_path, 'r') as gfa_file:
        for line in gfa_file:
            if line.startswith('S\t'):
                toks = line.rstrip('\n').split('\t')
                sn = [t[5:] for t in toks if t.startswith('SN:Z:')]
                sr = [t[5:] for t in toks if t.startswith('SR:i:')]
                if sn and sr and sr[0] == '0':
                    ref_contigs.add(sn[0])
    return ref_contigs

def strip_event_prefix(name):
    """ turn an rgfa-split bin name like id=CHM13|chr2 into chr2, the way split_gfa does """
    if name.startswith('id=') and name.find('|') > 3:
        return name[name.find('|') + 1:]
    return name

def cat_paths(in_paths, out_path):
    """ catFiles() interpolates its arguments into a shell command unquoted, so it cannot take a path
    holding a shell metacharacter.  these paths are rgfa-split bins, named after reference contigs
    whose sanitized form is id=EVENT|CONTIG -- don't make the concatenation depend on whether the
    prefix survives into the filename """
    with open(out_path, 'wb') as out_file:
        for in_path in in_paths:
            with open(in_path, 'rb') as in_file:
                shutil.copyfileobj(in_file, out_file)

def read_fasta_names(path):
    """ the sequence names in a fasta, which may or may not be compressed """
    with open(path, 'rb') as fa_file:
        is_gzipped = fa_file.read(2) == b'\x1f\x8b'
    opener = gzip.open if is_gzipped else open
    names = set()
    with opener(path, 'rt') as fa_file:
        for line in fa_file:
            if line.startswith('>'):
                names.add(line[1:].split()[0])
    return names

def bins_for_chrom(chrom, ref_contigs, bin_ref_contigs):
    """ which of a graph's reference contigs belong to one chromosome's bin.  a normal chromosome owns
    the contig it is named after; the --otherContig lump owns whichever reference contigs its own
    fasta holds.  only needed when the graph holds more of the reference than the chromosome it was
    built for, ie under mgSplitWholeGenomeRef.

    the lump is read off its fasta rather than inferred from the surviving bin names because the two
    disagree: export_split_data drops a reference contig that fewer than two genomes align to
    ("Omitting reference contig ..."), so it names no bin, and inferring would hand it to the lump --
    whose seqfile has no fasta for it, leaving cactus-align with PAF rows for an absent query """
    stripped = set(strip_event_prefix(name) for name in ref_contigs)
    if chrom in stripped:
        return {chrom}
    return stripped & set(strip_event_prefix(name) for name in (bin_ref_contigs or set()))

# Pruning what a wholeGenomeRef construction added to the rest of the reference
# (<graphmap_split pruneOtherContigs>).
#
# With mgSplitWholeGenomeRef each chromosome's second-pass minigraph is built against the whole
# reference, so construction can attach that chromosome's haplotypes to other chromosomes'
# reference as well: non-reference segments hanging only off another chromosome, and sample links
# between another chromosome's reference segments (deletion and translocation edges, eg --ggs-cnv's).
# The separation pass cuts all of that out again before cactus-align, but the zip and graphmap's
# mapping both see it first, and it distorts them: copies of a repeat array added to another
# chromosome's array (4q D4Z4 units on chr10's, in a chr4 graph) raise its minimizers' occurrence
# towards minigraph's cap, and a chr21 haplotype's 171 kb deletion edge on chr13 gives the chr21
# p-arm blocks a path to chain onto chr13 by.  So, before anything maps to the graph, keep
# only what this chromosome's construction could have added to its own reference: every rank-0
# segment and rank-0 link (the other chromosomes stay as flat linear reference, so they still compete
# for MAPQ), the non-reference components attached to this chromosome, and of the sample (rank > 0)
# links, those that touch no other chromosome's reference.

# <graphmap_split pruneOtherContigs> for a config without it
PRUNE_OTHER_CONTIGS_DEFAULT = True

# The prune is single-threaded python bar its bgzips (in, and the two graphs out), which are all
# that more cores speed up
PRUNE_MAX_CORES = 8

# Seconds per GB of compressed GFA, for the whole job: inflating it, the two python passes, the
# python PanSN rename of the pruned graph, and bgzipping both graphs.  Measured on a 30-way CHM13
# chr4 graph (767 MB bgzipped, 3.1 GB raw, 52,819 segments) at 8 cores: 141 s, of which the two
# bgzips of the pruned graphs were 120 s, the rename 8 s and inflating, scanning and writing 13 s.
# That is 184 s/GB; this allows 1.6x it
PRUNE_SECS_PER_GB = 300

# Memory.  What the scan keeps (a few arrays and a name index per segment) is small next to the
# longest line, which the python passes, and minigraph_gfa_to_pansn most of all, copy several times
# over.  On a whole-reference graph that is the longest reference contig, as one segment that no
# haplotype of this chromosome touched: CHM13's chr1 is a single 248 MB segment in the chr4 graph
# above, where the job peaked at 1.06 GB, 4.3x it.  So the job asks for 6x the longest reference
# contig (prune_other_contigs_memory), and never less than the floor, which leaves room for ten times
# that graph's segments
PRUNE_MEMORY_PER_CONTIG_BYTE = 6
PRUNE_MIN_MEMORY = 4 * 1024**3

def prune_other_contigs_enabled(config_node):
    """ <graphmap_split pruneOtherContigs>: whether each wholeGenomeRef graph is pruned of what its
    construction added to other chromosomes' reference """
    graphmap_split_node = findRequiredNode(config_node, "graphmap_split")
    try:
        return getOptionalAttrib(graphmap_split_node, "pruneOtherContigs", typeFn=bool,
                                 default=PRUNE_OTHER_CONTIGS_DEFAULT)
    except ValueError:
        raise RuntimeError('<graphmap_split pruneOtherContigs="{}"> must be 0 or 1'.format(
            graphmap_split_node.attrib.get("pruneOtherContigs")))

def prune_other_contigs_memory(longest_ref_contig):
    """ the memory of the prune of a whole-reference graph whose longest reference contig is
    longest_ref_contig bytes (PRUNE_MEMORY_PER_CONTIG_BYTE), before cactus_clamp_memory """
    return max(PRUNE_MIN_MEMORY, PRUNE_MEMORY_PER_CONTIG_BYTE * int(longest_ref_contig or 0))

def fasta_longest_record(job, fa_id):
    """ the size in bytes of the longest record of a fasta (uncompressed, as sanitize_fasta_headers
    writes them), headers and newlines included: an upper bound on its longest sequence """
    fa_path = os.path.join(job.fileStore.getLocalTempDir(), 'reference.fa')
    job.fileStore.readGlobalFile(fa_id, fa_path)
    with open(fa_path, 'rb') as fa_file:
        if fa_file.read(2) == b'\x1f\x8b':
            raise RuntimeError('{} is compressed, where a sanitized fasta is expected'.format(fa_id))
    longest = longest_record_size(fa_path)
    RealtimeLogger.info('The longest reference contig is at most {} bytes: each whole-reference graph prune '
                        '(<graphmap_split pruneOtherContigs>) gets {}'.format(
                            longest, bytes2human(cactus_clamp_memory(prune_other_contigs_memory(longest)))))
    return longest

def gfa_ref_contig(sn):
    """ the contig an rGFA SN tag names, without its genome: chr4 from either the cactus naming,
    id=CHM13|chr4 (as strip_event_prefix reads it), or PanSN, CHM13#0#chr4 (or SAMPLE#chr4) """
    if sn.startswith('id='):
        return strip_event_prefix(sn)
    return sn.split('#', 2)[-1]

def open_gfa_bytes(gfa_path):
    """ a GFA, gzipped or not, for reading as bytes """
    with open(gfa_path, 'rb') as gfa_file:
        is_gzipped = gfa_file.read(2) == b'\x1f\x8b'
    return gzip.open(gfa_path, 'rb') if is_gzipped else open(gfa_path, 'rb')

def gfa_tag_value(tags, key):
    """ the value of the first of tags (bytes) that starts with key (eg b'SR:i:'), or None """
    for tag in tags:
        if tag.startswith(key):
            return tag[len(key):]
    return None

class RgfaScan:
    """ what the prune needs of an rGFA, read without keeping its sequences: per segment (numbered in
    the order the file first names it) its rank, reference contig and length, and per L line (in file
    order) its two segments and rank.  minigraph writes S and L lines only (and an H line), and the
    prune only knows how to drop those, so a W or P line, which could name a dropped segment, is an
    error, as is an S or L line without its rank (SR:i:) """
    def __init__(self):
        self.index = {}              # segment name (bytes) -> number
        self.rank = array('q')       # SR:i:, -1 while the segment has no S line
        self.contig = array('q')     # a rank-0 segment's number in contigs, else -1
        self.length = array('q')     # LN:i:, or the sequence length without one
        self.contigs = []            # the reference contigs, as gfa_ref_contig names them
        self.contig_index = {}
        self.link_from = array('q')
        self.link_to = array('q')
        self.link_rank = array('q')  # SR:i:

    def segment(self, name):
        """ the number of the segment called name, which is given one if it is new """
        number = self.index.get(name)
        if number is None:
            number = self.index[name] = len(self.rank)
            self.rank.append(-1)
            self.contig.append(-1)
            self.length.append(0)
        return number

def gfa_segment_name_end(line, end):
    """ where the name of an S line ends """
    name_end = line.find(b'\t', 2, end)
    if name_end < 0:
        raise RuntimeError('Malformed GFA S line: {}'.format(line[:200]))
    return name_end

def scan_rgfa(gfa_path):
    """ an RgfaScan of a minigraph rGFA (gzipped or not), streamed a line at a time """
    scan = RgfaScan()
    with open_gfa_bytes(gfa_path) as gfa_file:
        for line in gfa_file:
            if line.startswith(b'S\t'):
                end = len(line.rstrip(b'\r\n'))
                name_end = gfa_segment_name_end(line, end)
                seq_end = line.find(b'\t', name_end + 1, end)
                tags = line[seq_end + 1:end].split(b'\t') if seq_end >= 0 else []
                if seq_end < 0:
                    seq_end = end
                number = scan.segment(line[2:name_end])
                if scan.rank[number] >= 0:
                    raise RuntimeError('Segment {} has two S lines in {}'.format(line[2:name_end].decode(), gfa_path))
                rank = gfa_tag_value(tags, b'SR:i:')
                if rank is None:
                    raise RuntimeError('Segment {} in {} has no SR:i: tag: it is not an rGFA'.format(
                        line[2:name_end].decode(), gfa_path))
                scan.rank[number] = int(rank)
                length = gfa_tag_value(tags, b'LN:i:')
                scan.length[number] = int(length) if length is not None else seq_end - name_end - 1
                if scan.rank[number] == 0:
                    sn = gfa_tag_value(tags, b'SN:Z:')
                    if sn is None:
                        raise RuntimeError('Rank-0 segment {} in {} has no SN:Z: tag'.format(
                            line[2:name_end].decode(), gfa_path))
                    contig = gfa_ref_contig(sn.decode())
                    if contig not in scan.contig_index:
                        scan.contig_index[contig] = len(scan.contigs)
                        scan.contigs.append(contig)
                    scan.contig[number] = scan.contig_index[contig]
            elif line.startswith(b'L\t'):
                toks = line.rstrip(b'\r\n').split(b'\t')
                rank = gfa_tag_value(toks[6:], b'SR:i:')
                if rank is None:
                    raise RuntimeError('Link {} in {} has no SR:i: tag: it is not an rGFA'.format(
                        b' '.join(toks[1:5]).decode(), gfa_path))
                scan.link_from.append(scan.segment(toks[1]))
                scan.link_to.append(scan.segment(toks[3]))
                scan.link_rank.append(int(rank))
            elif line.startswith(b'W\t') or line.startswith(b'P\t'):
                raise RuntimeError('{} has {} lines: the prune of a whole-reference minigraph only reads an rGFA '
                                   '(S and L lines), and cannot drop the segments of a path'.format(
                                       gfa_path, line[:1].decode()))
    missing = [name.decode() for name, number in scan.index.items() if scan.rank[number] < 0]
    if missing:
        raise RuntimeError('{} segment(s) in links of {} have no S line, eg {}'.format(
            len(missing), gfa_path, ' '.join(missing[:5])))
    return scan

def prune_own_contigs(chrom, ref_contigs, bin_ref_contigs=None):
    """ the reference contigs (of ref_contigs, the graph's) that belong to chrom's bin, by
    bins_for_chrom: the contig chrom is named after, or for the --otherContig lump, those its own
    reference fasta holds, which are bin_ref_contigs.  That can be a function returning them, which is
    only called for the lump """
    if chrom not in ref_contigs and callable(bin_ref_contigs):
        bin_ref_contigs = bin_ref_contigs()
    elif callable(bin_ref_contigs):
        bin_ref_contigs = None
    own = bins_for_chrom(chrom, set(ref_contigs), set(gfa_ref_contig(name) for name in (bin_ref_contigs or ())))
    if not own:
        raise RuntimeError('None of the {} reference contigs in {}\'s graph belong to it, so there is nothing to '
                           'prune it back to'.format(len(ref_contigs), chrom))
    return own

class OtherContigPrune:
    """ what prune_other_contigs_plan decides: drop_seg (per segment number) and drop_link (per L line)
    are 0 to keep, and for a link 1 where it touches a dropped segment, 2 where it is a sample link
    onto another chromosome's reference.  summary is what prune_other_contigs_report reads """
    def __init__(self, drop_seg, drop_link, summary):
        self.drop_seg = drop_seg
        self.drop_link = drop_link
        self.summary = summary

    def dropped(self):
        return self.summary['dropped_segments'] > 0 or self.summary['dropped_links'] > 0

def prune_other_contigs_plan(scan, own_contigs):
    """ what to drop from a whole-reference rGFA (an RgfaScan) to keep only what could have been added
    to own_contigs, its own reference contigs:
      - every non-reference component (the rank > 0 segments, joined by the links between them) not
        attached to an own contig's rank-0 segment, which includes one attached to no rank-0 segment.
        a component attached to both an own contig and another is kept, and counted as bridging
      - every link touching a dropped segment
      - every sample link (SR:i: > 0) touching another contig's rank-0 segment: deletion or
        translocation edges between other chromosomes' reference, and the links onto another
        chromosome from a kept (bridging) component or from this chromosome's own reference
    Every rank-0 segment and rank-0 link is kept """
    num_segs, num_links = len(scan.rank), len(scan.link_from)
    rank, contig, length = scan.rank, scan.contig, scan.length
    own = bytearray(len(scan.contigs))
    for name in own_contigs:
        if name in scan.contig_index:
            own[scan.contig_index[name]] = 1

    # the non-reference components, by union-find
    parent = array('q', range(num_segs))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in zip(scan.link_from, scan.link_to):
        if rank[a] > 0 and rank[b] > 0:
            root_a, root_b = find(a), find(b)
            if root_a != root_b:
                parent[root_a] = root_b
    touches = defaultdict(set)
    for a, b in zip(scan.link_from, scan.link_to):
        if rank[a] > 0 and rank[b] == 0:
            touches[find(a)].add(contig[b])
        if rank[b] > 0 and rank[a] == 0:
            touches[find(b)].add(contig[a])

    def contig_key(contig_numbers):
        return ','.join(sorted(scan.contigs[c] for c in contig_numbers)) or 'none'

    # each component's fate, decided at its first segment: (keep it, the tally it counts towards)
    component = {}
    dropped_components = defaultdict(lambda: [0, 0, 0])  # contigs it touches -> components, segments, bp
    bridging_components = defaultdict(lambda: [0, 0, 0])
    drop_seg = bytearray(num_segs)
    for seg in range(num_segs):
        if rank[seg] <= 0:
            continue
        root = find(seg)
        if root not in component:
            touched = touches.get(root, ())
            keep = any(own[c] for c in touched)
            if not keep:
                tally = dropped_components[contig_key(touched)]
            elif not all(own[c] for c in touched):
                tally = bridging_components[contig_key(touched)]
            else:
                tally = None
            if tally is not None:
                tally[0] += 1
            component[root] = (keep, tally)
        keep, tally = component[root]
        if not keep:
            drop_seg[seg] = 1
        if tally is not None:
            tally[1] += 1
            tally[2] += length[seg]

    drop_link = bytearray(num_links)
    links_to_dropped = 0
    sample_links = defaultdict(int)  # the other contigs it touches -> links
    for k in range(num_links):
        a, b = scan.link_from[k], scan.link_to[k]
        if drop_seg[a] or drop_seg[b]:
            drop_link[k] = 1
            links_to_dropped += 1
        elif scan.link_rank[k] > 0:
            others = set(contig[x] for x in (a, b) if rank[x] == 0 and not own[contig[x]])
            if others:
                drop_link[k] = 2
                sample_links[contig_key(others)] += 1

    summary = {'own_contigs': sorted(own_contigs),
               'segments': num_segs,
               'rank0_segments': sum(1 for r in rank if r == 0),
               'links': num_links,
               'dropped_segments': sum(drop_seg),
               'dropped_bp': sum(length[s] for s in range(num_segs) if drop_seg[s]),
               'dropped_links': links_to_dropped + sum(sample_links.values()),
               'dropped_components': {k: tuple(v) for k, v in dropped_components.items()},
               'bridging_components': {k: tuple(v) for k, v in bridging_components.items()},
               'links_to_dropped_segments': links_to_dropped,
               'sample_links_onto_other_contigs': dict(sample_links)}
    return OtherContigPrune(drop_seg, drop_link, summary)

def prune_other_contigs_report(summary):
    """ the lines that say what prune_other_contigs_plan dropped and kept """
    lines = ['own reference contig(s) {}; graph has {} segments ({} rank 0) and {} links'.format(
        ','.join(summary['own_contigs']), summary['segments'], summary['rank0_segments'], summary['links'])]
    for key, (comps, segs, bp) in sorted(summary['dropped_components'].items()):
        lines.append('dropped {} non-reference component(s) attached {}: {} segments, {} bp'.format(
            comps, 'only to ' + key if key != 'none' else 'to no reference contig', segs, bp))
    for key, (comps, segs, bp) in sorted(summary['bridging_components'].items()):
        lines.append('kept {} bridging component(s) attached to {}: {} segments, {} bp (their sample links onto '
                     'other contigs are dropped)'.format(comps, key, segs, bp))
    lines.append('dropped {} link(s) touching dropped segments'.format(summary['links_to_dropped_segments']))
    for key, count in sorted(summary['sample_links_onto_other_contigs'].items()):
        lines.append('dropped {} sample link(s) onto the reference of {}'.format(count, key))
    lines.append('kept {} of {} segments ({} bp dropped) and {} of {} links'.format(
        summary['segments'] - summary['dropped_segments'], summary['segments'], summary['dropped_bp'],
        summary['links'] - summary['dropped_links'], summary['links']))
    return lines

def write_pruned_gfa(gfa_path, out_path, scan, prune):
    """ stream gfa_path (the file scan was read from) to out_path (uncompressed), without the
    segments and links prune drops.  Every other line is copied as it was """
    num_links = 0
    with open_gfa_bytes(gfa_path) as gfa_file, open(out_path, 'wb') as out_file:
        for line in gfa_file:
            if line.startswith(b'S\t'):
                name_end = gfa_segment_name_end(line, len(line))
                if prune.drop_seg[scan.index[line[2:name_end]]]:
                    continue
            elif line.startswith(b'L\t'):
                num_links += 1
                if prune.drop_link[num_links - 1]:
                    continue
            out_file.write(line)
    if num_links != len(prune.drop_link):
        raise RuntimeError('{} changed while it was being pruned: {} links, where it had {}'.format(
            gfa_path, num_links, len(prune.drop_link)))

def prune_other_contigs_gfa(gfa_path, out_path, chrom, bin_ref_contigs=None):
    """ prune a whole-reference rGFA (gzipped or not, cactus or PanSN named) built for chrom's bin
    (prune_other_contigs_plan) to out_path, uncompressed.  bin_ref_contigs is only read for an
    --otherContig lump (prune_own_contigs).  Returns the summary prune_other_contigs_report reads.
    This is the whole prune outside Toil """
    scan = scan_rgfa(gfa_path)
    prune = prune_other_contigs_plan(scan, prune_own_contigs(chrom, scan.contigs, bin_ref_contigs))
    write_pruned_gfa(gfa_path, out_path, scan, prune)
    return prune.summary

def prune_other_contigs_workflow(job, chrom, gfa_id, pansn_gfa_id, gfa_path, graph_names, bin_ref_fa_id=None,
                                 longest_ref_contig=None, max_cores=1):
    """ size and schedule prune_other_contigs, which takes what chrom's wholeGenomeRef construction
    added to other chromosomes' reference out of both its graphs, the cactus-named gfa_id and the PanSN
    pansn_gfa_id (as minigraph_construct returned them).  A job of its own because the graph only has
    a size once it exists.  longest_ref_contig (fasta_longest_record of the whole reference) sets the
    memory, with PRUNE_MIN_MEMORY when it is not given.  Returns the two graphs, pruned """
    cores = max(1, min(int(max_cores or 1), PRUNE_MAX_CORES))
    gfa_bytes = gfa_id.size
    raw_bytes = RAW_BYTES_PER_GZ_BYTE * gfa_bytes if gfa_path.endswith('.gz') else gfa_bytes
    prune_job = job.addChildJobFn(prune_other_contigs, chrom, gfa_id, pansn_gfa_id, gfa_path, graph_names,
                                  bin_ref_fa_id,
                                  cores=cores,
                                  memory=cactus_clamp_memory(prune_other_contigs_memory(longest_ref_contig)),
                                  # the graph compressed, inflated and pruned, and the PanSN rename's
                                  # uncompressed copy and both compressed outputs
                                  disk=4 * gfa_bytes + 2 * raw_bytes + (bin_ref_fa_id.size if bin_ref_fa_id else 0),
                                  walltime=cactus_walltime(PRUNE_SECS_PER_GB * gfa_bytes / 1e9,
                                                           io_bytes=3 * gfa_bytes))
    return prune_job.rv(0), prune_job.rv(1)

def prune_other_contigs(job, chrom, gfa_id, pansn_gfa_id, gfa_path, graph_names, bin_ref_fa_id=None):
    """ prune chrom's whole-reference cactus-named graph gfa_id (prune_other_contigs_plan) and make its
    PanSN copy from the result, the way minigraph_construct made pansn_gfa_id from the unpruned one,
    so the two stay the same graph.  Both are compressed the way gfa_path says.  Returns the pair, or
    the graphs as they came in when there is nothing to drop """
    work_dir = job.fileStore.getLocalTempDir()
    name = os.path.basename(gfa_path)
    in_path = os.path.join(work_dir, 'unpruned.' + name)
    job.fileStore.readGlobalFile(gfa_id, in_path)
    with open(in_path, 'rb') as gfa_file:
        is_gzipped = gfa_file.read(2) == b'\x1f\x8b'
    raw_path = in_path
    if is_gzipped:
        # two passes over it in python, so inflate it once, with threads, rather than twice in python
        raw_path = os.path.join(work_dir, 'unpruned.gfa')
        cactus_call(parameters=['bgzip', '-dc', '--threads', str(job.cores), os.path.basename(in_path)],
                    work_dir=work_dir, outfile=raw_path)
    scan = scan_rgfa(raw_path)

    def read_bin_ref_contigs():
        if not bin_ref_fa_id:
            return set()
        bin_ref_path = os.path.join(work_dir, 'bin_reference.fa')
        job.fileStore.readGlobalFile(bin_ref_fa_id, bin_ref_path)
        return read_fasta_names(bin_ref_path)
    prune = prune_other_contigs_plan(scan, prune_own_contigs(chrom, scan.contigs, read_bin_ref_contigs))
    report = prune_other_contigs_report(prune.summary)
    if not prune.dropped():
        RealtimeLogger.info('Nothing to prune from {}\'s whole-reference graph: {}'.format(chrom, report[0]))
        return gfa_id, pansn_gfa_id
    RealtimeLogger.info('Pruning what construction added to other reference contigs from {}\'s whole-reference '
                        'graph:\n  {}'.format(chrom, '\n  '.join(report)))

    out_path = os.path.join(work_dir, 'pruned.' + name)
    raw_out_path = out_path[:-3] if out_path.endswith('.gz') else out_path
    write_pruned_gfa(raw_path, raw_out_path, scan, prune)
    if raw_path != in_path:
        os.remove(raw_path)
    pansn_out_path = os.path.join(work_dir, 'pan-sn.pruned.' + name)
    minigraph_gfa_to_pansn(set(graph_names), raw_out_path, pansn_out_path, job.cores)
    bgzip_gfa_rename(raw_out_path, out_path, job.cores)
    # the unpruned graphs are nobody's input from here on.  deleted from a follow-on rather than
    # inline so a retry of this job doesn't trip over its own delete
    job.addFollowOnJobFn(clean_jobstore_files, file_ids=[gfa_id, pansn_gfa_id], walltime=cactus_walltime())
    return job.fileStore.writeGlobalFile(out_path), job.fileStore.writeGlobalFile(pansn_out_path)

def separate_ref_contigs_batch(job, config, graphmap_input_dict, graphmap_batch_results, reference_event,
                              permissive_contig_filter, whole_genome_ref=False):
    """ --otherContig lumps every leftover reference contig into a single bin, so with --mgSplit that
    bin gets one minigraph covering all of them.  minigraph keeps those contigs disjoint, but the
    per-chromosome graphmap that follows can align a sample contig across two, which welds them into
    one component in the final graph.  Run rgfa-split over each bin to hold the contigs apart.

    With mgSplitWholeGenomeRef every bin's graph carries the whole reference, so this pass runs on
    all of them and does the real work of the option: dropping the mappings that landed on another
    chromosome, and cutting the graph artifacts back down to this chromosome.

    Appends the rgfa-split log to each chromosome's graphmap results, or None where the pass didn't
    run, so it can be exported alongside the PAFs, and (in whole-genome-reference mode) the pruned
    PanSN GFA after it. """
    # the first pass applies this inside its own job, so its copy of the config never reaches us
    apply_permissive_contig_filter(config, permissive_contig_filter)
    # the pruned fasta's bgzip and the PanSN rename both compress a whole-genome-scale GFA here
    mg_cores = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "cpu", typeFn=int, default=1)
    output_dict = {}
    for chrom, gm_result in graphmap_batch_results.items():
        seq_id_map, gfa_id = graphmap_input_dict[chrom][0], graphmap_input_dict[chrom][1]
        paf_id = gm_result[0]
        if not paf_id or not gfa_id:
            if whole_genome_ref:
                # every artifact this chromosome has is still whole-genome, and passing them through
                # would publish a whole reference under its name and hand cactus-align the same.
                # there is no correct thing to emit here, so say so rather than emit the wrong thing
                raise RuntimeError('{} reached the separation pass with no {}, so its whole-genome '
                                   'reference cannot be pruned away'.format(chrom, 'PAF' if not paf_id else 'GFA'))
            output_dict[chrom] = tuple(gm_result) + (None,)
            continue
        # the per-chromosome GFA comes in bgzipped, and we decompress it here: budget for that the
        # same way the first-pass split does.  rgfa-split's working set is what sets the memory, so
        # the extra whole-genome-reference files below are added to the disk budget alone -- this job
        # writes the minigraph fasta and streams the reference fasta for its names, it holds neither
        tot_size = gfa_id.size * 10 + paf_id.size
        disk_size = tot_size
        ref_fa_id = None
        if whole_genome_ref:
            disk_size += (gm_result[1].size if gm_result[1] else 0) * 10
            ref_fa_id = seq_id_map.get(reference_event) if reference_event else None
            disk_size += (ref_fa_id.size if ref_fa_id else 0) * 4
        output_dict[chrom] = job.addChildJobFn(separate_ref_contigs, config, chrom, gfa_id, gm_result, reference_event,
                                               whole_genome_ref=whole_genome_ref, ref_fa_id=ref_fa_id,
                                               genome_names=sorted(seq_id_map.keys()),
                                               cores=mg_cores,
                                               disk=disk_size * 5,
                                               memory=cactus_clamp_memory(tot_size * 3),
                                               walltime=cactus_walltime(120 + RGFA_SPLIT_SECS_PER_GB * tot_size / 1e9,
                                                                        io_bytes=tot_size)).rv()
    return output_dict

def separate_ref_contigs(job, config, chrom, gfa_id, gm_result, reference_event,
                         whole_genome_ref=False, ref_fa_id=None, genome_names=None):
    """ Use rgfa-split to assign each query in one chromosome's PAF to a single reference contig, and
    drop the alignments that don't fit that assignment.  Bins with a single reference contig -- every
    normal chromosome -- are left untouched.  A query that straddles two contigs evenly is ambiguous
    and loses all its alignments; one with a clear best contig keeps those and loses the rest.  Either
    way the sequence stays in the bin's fasta, so cactus-align sees it as unaligned and the usual
    non-minigraph clipping takes it out.

    Under mgSplitWholeGenomeRef the graph holds the whole reference rather than one chromosome, so
    every bin lands here with something to do: a query that maps better to another chromosome loses
    its alignments to this one, which is the inter-chromosome competition the whole-genome pipeline
    gets for free.  In that mode this pass is required for correctness, not just quality:
    extract_paf_from_gfa synthesizes rGFA-derived rows for every reference contig genome-wide, and it
    is the binning here that keeps the ones belonging to other chromosomes out of this PAF.  Only the bins this chromosome owns are kept, and the minigraph fasta and PanSN
    GFA are rebuilt off the pruned graph so nothing off-chromosome reaches cactus-align or the
    merged .sv.gfa.gz. """
    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, 'mg.gfa')
    paf_path = os.path.join(work_dir, 'mg.paf')
    out_prefix = os.path.join(work_dir, 'separate_')
    log_path = os.path.join(work_dir, 'separate.log')
    # what a bypass returns: the results as they came in, plus the log slot, plus (in whole-genome
    # reference mode) the pruned PanSN GFA slot
    passthrough = tuple(gm_result) + (None,) + ((None,) if whole_genome_ref else ())
    job.fileStore.readGlobalFile(gfa_id, gfa_path)
    with open(gfa_path, 'rb') as gfa_file:
        is_gzipped = gfa_file.read(2) == b'\x1f\x8b'
    if is_gzipped:
        unzipped_path = os.path.join(work_dir, 'mg.unzipped.gfa')
        cactus_call(parameters=['gzip', '-dc', os.path.basename(gfa_path)], work_dir=work_dir, outfile=unzipped_path)
        gfa_path = unzipped_path

    # nothing to hold apart on a normal chromosome: get out before touching the PAF, which is big.
    # a single-contig reference has nothing off-chromosome to prune either, so this is still the
    # right bypass in whole-genome reference mode
    ref_contigs = ref_contig_names(gfa_path)
    num_ref_contigs = len(ref_contigs)
    if num_ref_contigs < 2:
        return passthrough

    job.fileStore.readGlobalFile(gm_result[0], paf_path)
    coverage_opts, query_uniqueness, amb_name = get_split_specificity_opts(config)
    cmd = ['rgfa-split',
           '-p', paf_path,
           '-g', gfa_path,
           '-b', out_prefix,
           '-Q', query_uniqueness,
           '-a', amb_name,
           '-L', log_path]
    cmd += coverage_opts
    if whole_genome_ref:
        # we need the per-bin GFAs too, to cut the whole reference back down to this chromosome
        cmd += ['-G']
    if reference_event:
        cmd += ['-r', 'id={}|'.format(reference_event)]
    min_mapq = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minMAPQ")
    if min_mapq:
        cmd += ['-A', min_mapq]
    cactus_call(parameters=cmd, work_dir=work_dir, job_memory=job.memory)

    # which bins belong to this chromosome?  normally every non-ambiguous one: the bin was a
    # chromosome going in and rgfa-split only subdivided it.  with the whole reference in the graph
    # it split by every reference contig, so keep this chromosome's own contig -- or, for the
    # --otherContig lump, every reference contig that isn't some other chromosome's
    keep_names = None
    if whole_genome_ref:
        bin_ref_contigs = None
        if chrom not in set(strip_event_prefix(name) for name in ref_contigs) and ref_fa_id:
            bin_ref_path = os.path.join(work_dir, 'bin_reference.fa')
            job.fileStore.readGlobalFile(ref_fa_id, bin_ref_path)
            bin_ref_contigs = read_fasta_names(bin_ref_path)
        keep_names = bins_for_chrom(chrom, ref_contigs, bin_ref_contigs)

    # stitch the per-contig PAFs back into one: the bin stays whole, but no query is left
    # straddling two reference contigs.  the ambiguous bin is what we're dropping
    keep_pafs, keep_gfas = [], []
    for out_name in sorted(os.listdir(work_dir)):
        file_name, ext = os.path.splitext(out_name)
        if ext in ['.paf', '.gfa'] and file_name.startswith(os.path.basename(out_prefix)) and \
           os.path.isfile(os.path.join(work_dir, file_name + '.fa_contigs')):
            name = file_name[len(os.path.basename(out_prefix)):]
            if name == amb_name:
                continue
            if keep_names is not None and strip_event_prefix(name) not in keep_names:
                continue
            (keep_pafs if ext == '.paf' else keep_gfas).append(os.path.join(work_dir, out_name))

    separated_paf_path = os.path.join(work_dir, 'separated.paf')
    cat_paths(keep_pafs, separated_paf_path)

    def line_count(path):
        with open(path, 'r') as paf_file:
            return sum(1 for line in paf_file if line.strip())
    before, after = line_count(paf_path), line_count(separated_paf_path)
    RealtimeLogger.info('Separating {} reference contigs in {}: kept {} of them and {} of {} PAF records'.format(
        num_ref_contigs, chrom, len(keep_pafs), after, before))

    out_fa_id, pansn_gfa_id = gm_result[1], None
    if whole_genome_ref:
        if not keep_gfas:
            raise RuntimeError('rgfa-split left no reference contig assigned to {}: it cannot be pruned back '
                               'to a chromosome'.format(chrom))
        # a contig we meant to keep but that rgfa-split gave no bin contributes no rank-0 backbone to
        # the pruned graph, while the bin's seqfile still carries its fasta.  every reference contig
        # should have at least its own rgfa2paf row and so a bin, but say so rather than drop it quietly
        missing = keep_names - set(strip_event_prefix(os.path.splitext(os.path.basename(g))[0][len(os.path.basename(out_prefix)):])
                                   for g in keep_gfas)
        if missing:
            RealtimeLogger.warning('No rgfa-split bin for reference contig(s) {} in {}: they are absent from its '
                                   'pruned graph and minigraph fasta'.format(sorted(missing), chrom))
        graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
        # rgfa-split's GFA splitter re-emits only S and L lines, which is all minigraph wrote: its
        # GFAs carry no header, here or in any other run
        separated_gfa_path = os.path.join(work_dir, 'separated.gfa')
        cat_paths(keep_gfas, separated_gfa_path)

        # the minigraph fasta cactus-align consumes has to describe the pruned graph, not the
        # whole-genome one it was made from.  same command as make_minigraph_fasta()
        if gm_result[1]:
            out_fa_path = os.path.join(work_dir, 'separated.sv.gfa.fa.gz')
            cactus_call(outfile=out_fa_path,
                        parameters=[['gfatools', 'gfa2fa', separated_gfa_path],
                                    ['sed', '-e', r's/^>\(.\)/>id={}|\1/g'.format(graph_event)],
                                    ['bgzip', '--threads', str(job.cores)]])
            out_fa_id = job.fileStore.writeGlobalFile(out_fa_path)
            # the whole-genome one it replaces is nobody's input from here on.  deleted from a
            # follow-on rather than inline so a retry of this job doesn't trip over its own delete
            job.addFollowOnJobFn(clean_jobstore_files, file_ids=[gm_result[1]], walltime=cactus_walltime())

        # and the published per-chromosome graph, which graphmap-join merges
        pansn_gfa_path = os.path.join(work_dir, 'separated.sv.gfa.gz')
        minigraph_gfa_to_pansn(set(genome_names or []), separated_gfa_path, pansn_gfa_path, job.cores)
        pansn_gfa_id = job.fileStore.writeGlobalFile(pansn_gfa_path)

    return (job.fileStore.writeGlobalFile(separated_paf_path), out_fa_id) + tuple(gm_result[2:]) + \
        (job.fileStore.writeGlobalFile(log_path),) + ((pansn_gfa_id,) if whole_genome_ref else ())

# Seconds per reference contig for split_fa_into_contigs.  The contig count is the driver rather
# than the fasta size because each contig costs its own faidx and bgzip, and <graphmap_split
# maxRefContigs> defaults to 128 -- five times the human case.
#
# Re-measured by grouping every timed command in the two HPRC v2.1 runs by the job work directory
# it ran in: 459 jobs per run, ~50 commands each, summing to a p99 of 646 s (sep16) and 604 s
# (sep17) against a worst of 654 s, over a 196-contig GRCh38 split.  That is 3 s/contig at the
# p99, which is the convention WALLTIME_FACTOR's comment sets for these constants.  The earlier
# 13 s/contig came from grouping by a 23-26 contig split, where the per-job fixed cost is spread
# over a fifth as many contigs and so lands in the slope.
#
# The difference is worth the re-fit because of how many of these there are: at 15 all 459 asked
# 2.1 h each, and at 3 they ask 37 min, which is the side of the one-hour partition they belong on.
SPLIT_FA_SECS_PER_REF_CONTIG = 3

def split_fas(job, seq_id_map, seq_name_map, split_id_map):
    """ Use samtools to split a bunch of fasta files into reference contigs, using the output of rgfa-split as a guide"""

    if (seq_id_map, split_id_map) == (None, None):
        return None

    root_job = Job(walltime=cactus_walltime())
    job.addChild(root_job)

    # map event name to dict of contgs.  ex fa_contigs["CHM13"]["chr13"] = file_id
    fa_contigs = {}
    # map event name to dict of contig sizes ex fa_contigs["CHM13"]["chr13"] = N
    # where N is total number of fasta bases
    fa_contig_sizes = {}
    
    # we do each fasta in parallel
    for event, fa_id in seq_id_map.items():
        fa_path = seq_name_map[event]
        if fa_id.size:
            split_job = root_job.addChildJobFn(split_fa_into_contigs, event, fa_id, fa_path, split_id_map,
                                               strip_prefix=False, 
                                               disk=fa_id.size * 3,
                                               walltime=cactus_walltime(60 + SPLIT_FA_SECS_PER_REF_CONTIG * len(split_id_map),
                                                                        io_bytes=2*fa_id.size))
            fa_contigs[event] = split_job.rv(0)
            fa_contig_sizes[event] = split_job.rv(1)

    return fa_contigs, fa_contig_sizes

def split_fa_into_contigs(job, event, fa_id, fa_name, split_id_map, strip_prefix=False):
    """ Use samtools turn on fasta into one for each contig. this relies on the informatino in .fa_contigs
    files made by rgfa-split """

    # download the fasta
    work_dir = job.fileStore.getLocalTempDir()
    fa_path = os.path.join(work_dir, '{}.fa'.format(event))
    is_gz = fa_name.endswith(".gz")
    job.fileStore.readGlobalFile(fa_id, fa_path)

    unique_id = 'id={}|'.format(event)
                        
    contig_fa_dict = {}
    contig_size_dict = {}

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
                    if strip_prefix:
                        query_contig = query_contig[len(unique_id):]                        
                    clean_file.write('{}\n'.format(query_contig))
                    contig_count += 1
        contig_fasta_path = os.path.join(work_dir, '{}_{}.fa'.format(event, ref_contig))
        if contig_count > 0:
            cmd = [['samtools', 'faidx', fa_path, '--region-file', faidx_input_path]]
            # transform chr1:10-15 (1-based inclusive) into chr1_sub_9_15 (0-based end open)
            # this is a format that contains no special characters in order to make assembly hubs
            # happy.  But it does require conversion going into vg which wants chr[9-15] and
            # hal2vg is updated to do this autmatically
            cmd.append(get_faidx_subpath_rename_cmd())
            if is_gz:
                cmd.append(['bgzip', '--threads', str(job.cores)])
            cactus_call(parameters=cmd, outfile=contig_fasta_path)

            # now get the total sequence size
            size_cmd = [['cat', contig_fasta_path], ['grep', '-v', '>'], ['wc'], ['awk', '{printf \"%f\", $3-$1}']]
            if is_gz:
                size_cmd[0] = ['bgzip', '-dc', contig_fasta_path, '--threads', str(job.cores)]
            num_bases = int(float(cactus_call(parameters=size_cmd, check_output=True).strip()))
        else:
            # TODO: review how cases like this are handled
            with open(contig_fasta_path, 'w') as empty_file:
                empty_file.write("")
            num_bases = 0
        contig_fa_dict[ref_contig] = job.fileStore.writeGlobalFile(contig_fasta_path)
        contig_size_dict[ref_contig] = num_bases
        
    return contig_fa_dict, contig_size_dict

def gather_fas(job, output_id_map, contig_fa_map, contig_size_map):
    """ take the split_fas output which has everything sorted by event, and move into the ref-contig-based table
    from split_gfa.  return the updated table, which can then be exported into the chromosome projects """

    if not contig_fa_map:
        return None

    events = set()
    for ref_contig in output_id_map.keys():
        output_id_map[ref_contig]['fa'] = {}
        for event, fa_id in contig_fa_map.items():
            if ref_contig in fa_id:
                output_id_map[ref_contig]['fa'][event] = fa_id[ref_contig]
            events.add(event)

    # go ahead and make the size table file here
    # download the fasta
    work_dir = job.fileStore.getLocalTempDir()
    size_table_path = os.path.join(work_dir, 'contig_sizes.tsv')
    events = sorted(events)
    with open(size_table_path, 'w') as size_table_file:
        #write the header
        size_table_file.write("Contig\t" + '\t'.join(events) + '\tmin\tmax\tavg\n')
        for ref_contig in output_id_map.keys():
            sizes = [contig_size_map[event][ref_contig] for event in events]
            #total size for each event
            size_table_file.write(ref_contig + '\t' + '\t'.join([str(s) for s in sizes]))
            #tack on 3 basic stats
            size_table_file.write('\t{}\t{}\t{}\n'.format(min(sizes), max(sizes), int(sum(sizes) / len(sizes)) if sizes else 0))
    contig_size_table_id = job.fileStore.writeGlobalFile(size_table_path)

    return output_id_map, contig_size_table_id

def bin_other_contigs(job, config, ref_contigs, other_contig, output_id_map):
    """ take all the other (ie non ref) contigs (in practice, unplaced bits of grch38) and merge them
    all up into a single "other" contig.  this avoids a 1000 align jobs getting created downstream. 
    Note we've moved this to the ned here (it used to be done usinc -c -o in rgfa-split) so as to
    avoid letting these contigs glom together into componenets: the final output will be identital
    to if they were kept in separate files. """
    work_dir = job.fileStore.getLocalTempDir()
    if not ref_contigs or not other_contig or not output_id_map:
        return output_id_map

    amb_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "ambiguousName", default="_AMBIGUOUS_")
    
    ref_contigs = set(ref_contigs)
    other_contigs = set()
    merged_output_path_map = {} # map[TYPE][EVENT] = local file path

    for ref_contig in output_id_map.keys():
        if ref_contig not in ref_contigs and ref_contig != amb_name:
            other_contigs.add(ref_contig)
            for type_key in output_id_map[ref_contig].keys():
                if type_key == 'fa':
                    for event in output_id_map[ref_contig][type_key].keys():
                        # we read the existing id
                        local_path = os.path.join(work_dir, '{}_{}.{}'.format(ref_contig, event, type_key))
                        job.fileStore.readGlobalFile(output_id_map[ref_contig][type_key][event], local_path)

                        # find its new home
                        other_path = os.path.join(work_dir, '{}_{}.{}'.format(other_contig, event, type_key))

                        # update the map
                        if type_key not in merged_output_path_map:
                            merged_output_path_map[type_key] = {}
                        if event not in merged_output_path_map[type_key]:
                            merged_output_path_map[type_key][event] = other_path

                        # append the file
                        with open(local_path, 'rb') as local_file, open(other_path, 'ab') as other_file:
                            shutil.copyfileobj(local_file, other_file)
                else:
                    # we read the existing id
                    local_path = os.path.join(work_dir, '{}.{}'.format(ref_contig, type_key))
                    job.fileStore.readGlobalFile(output_id_map[ref_contig][type_key], local_path)

                    # find its new home
                    other_path = os.path.join(work_dir, '{}.{}'.format(other_contig, type_key))

                    # update the map
                    if type_key not in merged_output_path_map:
                        merged_output_path_map[type_key] = other_path

                    # append the file
                    with open(local_path, 'rb') as local_file, open(other_path, 'ab') as other_file:
                        shutil.copyfileobj(local_file, other_file)
                    
    # and convert to ids
    for type_key in merged_output_path_map.keys():
        if type_key == 'fa':
            for event in merged_output_path_map[type_key].keys():
                merged_output_path_map[type_key][event] = job.fileStore.writeGlobalFile(merged_output_path_map[type_key][event])
        else:
            merged_output_path_map[type_key] = job.fileStore.writeGlobalFile(merged_output_path_map[type_key])

    # remove the other contigs from the map
    for contig in other_contigs:
        del output_id_map[contig]

    # drop in the other contig
    output_id_map[other_contig] = merged_output_path_map

    return output_id_map

def export_split_data(toil, input_name_map, output_id_map, split_log_id, contig_size_table_id, output_dir, config):
    """ download all the split data locally """

    amb_name = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"), "ambiguousName", default="_AMBIGUOUS_")
    
    chrom_file_map = {}

    # export the sizes
    contig_size_table_path = os.path.join(output_dir, 'contig_sizes.tsv')
    toil.exportFile(contig_size_table_id, makeURL(contig_size_table_path))

    # export the log
    toil.exportFile(split_log_id, makeURL(os.path.join(output_dir, 'minigraph.split.log')))

    # hack to filter out reference contigs where nothing maps to it (not even itself)
    # this is usually worked around by using --refContigs <chroms> --otherContig chrOther
    # but filtering them out here lets users who don't do that make graphs
    # 
    empty_contigs = set()
    
    for ref_contig in output_id_map.keys():        
        if output_id_map[ref_contig] is None or len(output_id_map[ref_contig]) == 0:
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
            fa_path = makeURL(os.path.join(fa_base, '{}_{}.fa'.format(pangenome_file_name(event), ref_contig)))
            if input_name_map[event].endswith('.gz'):
                fa_path += '.gz'
            seq_file_map[event] = fa_path
            toil.exportFile(ref_contig_fa_id, fa_path)

        seq_file_map_size = 0
        for event, fa_path in seq_file_map.items():
            if output_id_map[ref_contig]['fa'][event].size > 0:
                seq_file_map_size += 1
        if seq_file_map_size < 2:
            logger.warning("Omitting reference contig {} from the graph because it doesn't align well enough. If you absolutely want to include it, rerun with --refContigs and --otherContig specified".format(ref_contig))
            empty_contigs.add(ref_contig)
            continue
        
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
                    seq_file.write('{}\t{}\n'.format(pangenome_seqfile_name(event), fa_path))
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
            if ref_contig != amb_name and ref_contig not in empty_contigs:
                seqfile, paf = seqfile_paf[0], seqfile_paf[1]
                chromfile.write('{}\t{}\t{}\n'.format(ref_contig, seqfile, paf))
    if chrom_file_path.startswith('s3://'):
        write_s3(chrom_file_temp_path, chrom_file_path)
        
if __name__ == "__main__":
    main()
