#!/usr/bin/env python3

"""
This is a spiritual successor to hal2mafMP.py (from hal).  putting it here in cactus as cactus already has the setup tools and docker command stuff
needed toil autoscale, and it'll be easier to keep all dependencies managed in submodules.  
"""

import os, sys
from argparse import ArgumentParser
import copy
import timeit, time
import math
import xml.etree.ElementTree as ET

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
from cactus.shared.common import clean_jobstore_files

from cactus.shared.version import cactus_commit
from cactus.progressive.cactus_prepare import human2bytesN
from cactus.progressive.multiCactusTree import MultiCactusTree

from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count
from toil.lib.humanize import bytes2human
from sonLib.bioio import getTempDirectory
from cactus.shared.common import cactus_clamp_memory
from sonLib.nxnewick import NXNewick

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("halFile", help = "HAL file to convert to MAF")
    parser.add_argument("outputMAF", help = "Output MAF (will be gzipped if ends in .gz). Suffix with .taf or .taf.gz if you want TAF output")
    parser.add_argument("--batchSize", type=int, help = "Number of chunks for each hal2maf batch", default=None)
    parser.add_argument("--batchCount", type=int, help = "Number of hal2maf batches [default 1 unless --batchSize set]", default=None)
    parser.add_argument("--batchCores", type=int, help = "Number of cores for each hal2maf batch.")
    parser.add_argument("--batchMemory", type=human2bytesN,
                        help="Memory in bytes for each hal2maf batch (defaults to an estimate based on the input data size). "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))", default=None)
    parser.add_argument("--chunkSize", type=int, help = "Size of chunks to operate on.", required=True)
    parser.add_argument("--batchParallelHal2maf", type=int, help = "Number of hal2maf commands to be executed in parallel in batch. Use to throttle down number of concurrent jobs to save memory. [default=batchCores]", default=None)
    parser.add_argument("--batchParallelTaf", type=int, help = "Number of taf normalization command chains to be executed in parallel in batch. Use to throttle down number of concurrent jobs to save memory. [default=batchCores]", default=None)    
    parser.add_argument("--raw", action="store_true", help = "Do not run taf-based normalization on the MAF")

    # new dupe-handler option
    parser.add_argument("--dupeMode", type=str, choices=["single", "consensus", "ancestral", "all"],
                        help="Toggle how to handle duplications: single: heuristically choose single, most similar homolog; consensus: squish all duplicate rows into a single conensus row; Ancestral: keep only duplications that are also separate in ancestor; All: keep all duplications, including paralogies (self-alignments) in given genome [default=all]",
                        default="all")

    # toggle taffy indexing
    parser.add_argument("--index", action = "store_true", help = "Produce Taffy index (.tai) of the output", default=False)
    
    # pass through a subset of hal2maf options
    parser.add_argument("--refGenome", required=True,
                        help="name of reference genome (root if empty)",
                        default=None)
    parser.add_argument("--refSequence",
                        help="subset to this contig in reference genome (multiple allowed) [default=all]",
                        nargs='+',
                        default=None)
    parser.add_argument("--start",
                        help="restrict to subrange of --refSequence beginning at this 0-based position (multiple allowed)",
                        type=int,
                        nargs='+',
                        default=None),
    parser.add_argument("--length",
                        help="restrict to subrange of --refSequence using this length (multiple allowed)",
                        type=int,
                        nargs='+',
                        default=None),
    parser.add_argument("--bedRanges",
                        help="coordinates in BED format of sequence ranges in --refGenome to export")
    parser.add_argument("--rootGenome",
                        help="name of root genome (none if empty)",
                        default=None)
    parser.add_argument("--targetGenomes",
                        help="comma-separated (no spaces) list of target "
                        "genomes (others are excluded) (vist all if empty)",
                        default=None)
    parser.add_argument("--noAncestors",
                        help="don't write ancestral sequences. IMPORTANT: "
                        "Must be used in conjunction with --refGenome"
                        " to set a non-ancestral genome as the reference"
                        " because the default reference is the root.",
                        action="store_true",
                        default=False)    
    
    # pass through taffy norm options
    parser.add_argument("--maximumBlockLengthToMerge",
                        help="Only merge together any two adjacent blocks if one or both is less than this many bases long, [default: see taffy norm -h]",
                        type=int,
                        default=None)
    parser.add_argument("--maximumGapLength",
                         help="Only merge together two adjacent blocks if the total number of unaligned bases between the blocks is less than this many bases, [default: see taffy norm -h]",
                         type=int,
                         default=None)
    parser.add_argument("--fractionSharedRows",
                        help="The fraction of rows between two blocks that need to be shared for a merge, [default: see taffy norm -h]",
                        type=float,
                        default=None)
    parser.add_argument("--filterGapCausingDupes",
                        help="Turn on (experimental) taffy norm -d filter that removes duplications that would induce gaps > maximumGapLength [default=Off unless --dupeMode single]",
                        action="store_true")
    parser.add_argument("--maxRefNFrac",
                        help="Filter out MAF blocks whose reference (first) line has a greater fraction of Ns than the given amount. Should be between 0.0 (filter everything) and 1.0 (filter nothing). [default=0.95]",
                        type=float,
                        default=0.95)

    # coverage options
    parser.add_argument("--coverage",
                        help="Produce coverage TSV for output alignment",
                        action="store_true")
    parser.add_argument("--coverageSexChroms",
                        help="Add breakdown of sex chromosomes vs autosomes to output coverage, using this space-separated list of sex chromosome. Don't prefix with reference genome name. Ex --coverageSexChroms chrX chrY",
                        nargs='+')
    parser.add_argument("--coverageGapThresholds",
                        help="add breakdown of coverage by maximum gap size using given thresholds. Ex. --overageGapThresholds 10 50 1000 1000000",
                        nargs='+')
                        
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

    if options.batchSize and options.batchCount:
        raise RuntimeError('Only one of --batchSize and --batchCount can be specified')

    if not options.outputMAF.endswith('.maf') and not options.outputMAF.endswith('.maf.gz') and \
       not options.outputMAF.endswith('.taf') and not options.outputMAF.endswith('.taf.gz'):
        raise RuntimeError('Output path must end with .maf, .maf.gz, .taf or .taf.gz')

    if options.dupeMode == 'single' and options.raw:
        raise RuntimeError('--dupeMode single requires TAF normalization and therefore is incompatible with --raw')
    if options.noAncestors and not options.refGenome:
        raise RuntimeError('(non-ancestral) --refGenome required with --noAncestors')
    if options.start is not None or options.length is not None:
        if options.start is None or options.length is None or options.refGenome is None or options.refSequence is None:
            raise RuntimeError('in order to use --start or --length, all of --start, --length, --refGenome and --refSequence must be specified')
        if len(options.start) != len(options.length) or len(options.start) != len(options.refSequence):
            raise RuntimeError('--start and --length must have the same number of values as --refSequence')
        if options.bedRanges is not None:
            raise RuntimeError('--start and --length cannot be used with --bedRanges')
    if options.bedRanges and options.refSequence:
        raise RuntimeError('--refSequence cannot be used with --bedRanges')
    if options.bedRanges and not options.bedRanges.endswith('.bed'):
        raise RuntimeError('file passed to --bedRanges must end with .bed')

    if options.dupeMode == 'single':
        options.filterGapCausingDupes = True

    if options.coverageSexChroms or options.coverageGapThresholds:
        logger.warning('Toggling --coverage on because --coverageSexChroms and/or --coverageGapThresholds was set')
        options.coverage = True

    # export job needs absolute path
    if not options.outputMAF.startswith('s3://'):
        options.outputMAF = os.path.abspath(options.outputMAF)
        
    # apply cpu override                
    if options.batchCores is None:
        if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
            options.batchCores = cactus_cpu_count()
            if options.maxCores:
                options.batchCores = min(options.batchCores, int(options.maxCores))
            logger.info('Setting batchCores to {}'.format(options.batchCores))
        else:
            raise RuntimeError('--batchCores must be specified for batch systems other than singleMachine')

    if not options.batchParallelHal2maf:
        options.batchParallelHal2maf = options.batchCores
    if options.batchParallelHal2maf > options.batchCores:
        raise RuntimeError('--batchParallelHal2maf cannot exceed the number of batch cores ({})'.format(options.batchCores))
    if not options.batchParallelTaf:
        options.batchParallelTaf = options.batchCores            
    if options.batchParallelTaf > options.batchCores:
        raise RuntimeError('--batchParallelTaf cannot exceed the number of batch cores ({})'.format(options.batchCores))
        
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    if not options.batchCount and not options.batchSize:
        logger.info('Using default batch count of 1')
        options.batchCount = 1
                    
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            maf_id = toil.restart()
        else:
            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals(options)

            bed_id = toil.importFile(options.bedRanges) if options.bedRanges else None
            logger.info("Importing {}".format(options.halFile))
            hal_id = toil.importFile(options.halFile)
            maf_id = toil.start(Job.wrapJobFn(hal2maf_workflow, hal_id, bed_id, options, config))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-hal2maf has finished after {} seconds".format(run_time))

def export_file(job, file_id, out_path):
    """ run toil export in its own job in order to a) do it right away but b) in a separate job (in case it fails) """
    job.fileStore.exportFile(file_id, makeURL(out_path))
    
def hal2maf_workflow(job, hal_id, bed_id, options, config):

    hal2maf_ranges_job = job.addChildJobFn(hal2maf_ranges, hal_id, bed_id, options, cores=1, disk=hal_id.size)
    chunks, genome_list = hal2maf_ranges_job.rv(0), hal2maf_ranges_job.rv(1)
    hal2maf_all_job = hal2maf_ranges_job.addFollowOnJobFn(hal2maf_all, hal_id, chunks, genome_list, options, config)
    hal2maf_merge_job = hal2maf_all_job.addFollowOnJobFn(hal2maf_merge, hal2maf_all_job.rv(), options, disk=hal_id.size)
    hal2maf_merge_job.addFollowOnJobFn(clean_jobstore_files, file_ids=hal2maf_all_job.rv())
    maf_id = hal2maf_merge_job.rv()
    hal2maf_merge_job.addFollowOnJobFn(export_file, maf_id, options.outputMAF)

    if options.index:
        index_job = hal2maf_merge_job.addFollowOnJobFn(taffy_index, maf_id, options.outputMAF, disk=hal_id.size)
        index_job.addFollowOnJobFn(export_file, index_job.rv(), options.outputMAF + '.tai')

    if options.coverage:
        coverage_job = hal2maf_merge_job.addFollowOnJobFn(taffy_coverage, maf_id, genome_list, options, disk=hal_id.size,
                                                          memory=cactus_clamp_memory(hal_id.size / 8))
        coverage_job.addFollowOnJobFn(export_file, coverage_job.rv(), options.outputMAF + '.cov.tsv')
        

def hal2maf_ranges(job, hal_id, bed_id, options):
    """ get the ranges (in reference *sequence* coordinates) for each hal2maf job """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    bed_path = os.path.join(work_dir, 'ranges.bed')
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))    
    job.fileStore.readGlobalFile(hal_id, hal_path)
    if bed_id:
        job.fileStore.readGlobalFile(bed_id, bed_path)
    RealtimeLogger.info("Computing range information")

    res = cactus_call(parameters=['halStats', hal_path, '--sequenceStats', options.refGenome], check_output=True)
    ref_sequence_lengths = {}
    for line in res.strip().split('\n')[1:]:
        tokens = line.strip().split(",")
        if len(tokens) == 4:
            ref_sequence_lengths[tokens[0]] = int(tokens[1])

    # upstream error handling should make sure that start/length, if specified, are 1:1 for refsequences
    if options.start:
        assert options.refSequence and options.length and len(options.start) == len(options.refSequence)
    if options.length:
        assert options.start and len(options.length) == len(options.start)

    # collect the ranges to convert, which can be specified using various options
    raw_intervals = []
    if options.start:
        # from --refSequence --start --length
        for seq, start, length in zip(options.refSequence, options.start, options.length):
            raw_intervals.append((seq, start, start + length))
    elif options.bedRanges:
        # from BED file with --bedRanges
        with open(bed_path, 'r') as bed_file:
            for line in bed_file:
                toks = line.split()
                if len(toks) >= 3:
                    raw_intervals.append((toks[0], int(toks[1]), int(toks[2])))
    else:
        if options.refSequence:
            # from --refSequence (but not subrange)
            for seq in options.refSequence:
                raw_intervals.append((seq, 0, ref_sequence_lengths[seq] if seq in ref_sequence_lengths else -1))
        else:
            # default to everything
            for seq, length in ref_sequence_lengths.items():
                raw_intervals.append((seq, 0, length))

    # sanity check the invervals
    for seq, start, end in raw_intervals:
        if seq not in ref_sequence_lengths:
            raise RuntimeError('target reference sequence {} not found in genome {}'.format(seq, options.refGenome))
        if start < 0:
            raise RuntimeError('negative start value {} not allowed'.format(start))
        if end > ref_sequence_lengths[seq]:
            raise RuntimeError('invalid range start={} end={} exceeds {} sequence length of {}'.format(start, end, seq, ref_sequence_lengths[seq]))
                                                                                                          

    # chunk the intervals
    # note this breaks up, but doesn't fuse input chunks
    # todo: if we want to support fusing, then probably need to make sure normalization respects boundaries
    chunks = []
    for seq, start, end in raw_intervals:
        num_chunks = max(1, int((end - start) / options.chunkSize))
        last_chunk_size = int(end - start) - (num_chunks - 1) * options.chunkSize
        # avoid huge chunks
        if last_chunk_size > options.chunkSize * 1.5 and num_chunks:
            num_chunks += 1
            last_chunk_size = int(end - start) - (num_chunks - 1) * options.chunkSize
        for i in range(num_chunks):
            chunk_size = options.chunkSize if i < num_chunks - 1 else last_chunk_size
            chunk_start = start + i * options.chunkSize
            chunk_end = chunk_start + chunk_size
            chunks.append((seq, chunk_start, chunk_end))
    assert chunks

    # get list of genomes sorted by their distance to reference in pre-order traversal
    tree_str = cactus_call(parameters=['halStats', hal_path, '--tree'], check_output=True).strip()
    mc_tree = MultiCactusTree(NXNewick().parseString(tree_str, addImpliedRoots=False))
    genome_ranks = {}
    for i, node in enumerate(mc_tree.preOrderTraversal()):
        genome_ranks[mc_tree.getName(node)] = i
    genome_list = sorted(list(genome_ranks.keys()),
                         key = lambda genome : abs(genome_ranks[genome] - genome_ranks[options.refGenome]))
    assert genome_list[0] == options.refGenome
    return chunks, genome_list

def hal2maf_all(job, hal_id, chunks, genome_list, options, config):
    """ make a job for each batch of chunks """
    num_batches = options.batchCount
    if not num_batches:
        # default to 1 unless batch_size is set
        if options.batchSize:
            num_batches = int(len(chunks) / options.batchSize)
            remainder = len(chunks) % options.batchSize
            if remainder:
                num_batches += 1
        else:
            num_batches = 1
        RealtimeLogger.info('Setting batchCount to {}'.format(num_batches))
            
    batch_size = options.batchSize
    if not batch_size:
        batch_size = math.ceil(len(chunks) / num_batches)
        RealtimeLogger.info('Setting batchSize to {}'.format(batch_size))

    # try to figure out memory
    if options.batchMemory:
        batch_memory = options.batchMemory
    else:
        maf_mem = math.ceil(hal_id.size / 200)
        taf_mem = math.ceil(hal_id.size / 33)
        batch_memory = cactus_clamp_memory(max(maf_mem * options.batchParallelHal2maf, taf_mem * options.batchParallelTaf))
        
    chunks_left = len(chunks)
    batch_results = []        
    for i in range(num_batches):
        cur_chunk = i * batch_size
        cur_batch_size = min(chunks_left, batch_size)
        if cur_batch_size:
            batch_results.append(job.addChildJobFn(hal2maf_batch, hal_id, chunks[cur_chunk:cur_chunk+cur_batch_size],
                                                   genome_list, options, config,
                                                   disk=math.ceil((1 + 1.5 / num_batches)*hal_id.size), cores=options.batchCores,
                                                   memory=batch_memory).rv())
        chunks_left -= cur_batch_size
    assert chunks_left == 0
    
    return batch_results

def chunk_name(chunk_num, options, tag=''):
    """ get a readable maf chunk name """
    bname, ext = os.path.splitext(os.path.basename(options.outputMAF))
    assert ext in ['.maf', '.taf', '.gz']
    if ext == '.gz':
        bname, ext2 = os.path.splitext(bname)
        assert ext2 in ['.maf', '.taf']
        ext = ext2 + ext
    # note, we're stuck with keeping chunks in MAF because TAF merging somehow isn't working:
    return '{}_chunk_{}{}{}'.format(bname, chunk_num, tag, ext.replace('taf', 'maf'))
    
def hal2maf_cmd(hal_path, chunk, chunk_num, options, config):
    """ make a hal2maf command for a chunk """

    # we are going to run this relative to work_dir
    hal_path = os.path.basename(hal_path)

    time_cmd = '/usr/bin/time -vp' if os.environ.get("CACTUS_LOG_MEMORY") else '(time -p '
    time_end = '' if os.environ.get("CACTUS_LOG_MEMORY") else ')'
    
    cmd = 'set -eo pipefail && {} hal2maf {} stdout --refGenome {} --refSequence {} --start {} --length {} --maxBlockLen 10000'.format(
        time_cmd, hal_path, options.refGenome, chunk[0], chunk[1], chunk[2]-chunk[1])
    if options.rootGenome:
        cmd += ' --rootGenome {}'.format(options.rootGenome)
    if options.targetGenomes:
        cmd += ' --targetGenomes {}'.format(options.targetGenomes)
    if options.dupeMode == 'ancestral':
        cmd += ' --onlyOrthologs'
    if options.noAncestors:
        cmd += ' --noAncestors'
    cmd += '{} 2> {}.h2m.time'.format(time_end, chunk_num)
    # todo: it would be better to filter this out upstream, but shouldn't make any difference here since we're taf normalizing anyway
    cmd += ' | grep -v "^s	{}"'.format(getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_"))
    if chunk_num != 0 and options.raw:
        cmd += ' | grep -v ^#'        
    if options.outputMAF.endswith('.gz'):
        cmd += ' | bgzip'        
    cmd += ' > {}'.format(chunk_name(chunk_num, options))
    return cmd

def taf_cmd(hal_path, chunk, chunk_num, genome_list_path, sed_script_paths, options):
    """ make a taf normalization command for the chunk """
    # we are going to run this relative to work_dir
    hal_path = os.path.basename(hal_path)

    time_cmd = '/usr/bin/time -vp' if os.environ.get("CACTUS_LOG_MEMORY") else '(time -p '
    time_end = '' if os.environ.get("CACTUS_LOG_MEMORY") else ')'
    read_cmd = 'gzip -dc' if options.outputMAF.endswith ('.gz') else 'cat'

    # we don't pipe directly from hal2maf because add_gap_bases (now norm) uses even more memory in hal
    cmd = 'set -eo pipefail && {} {}'.format(read_cmd, chunk_name(chunk_num, options))
    if sed_script_paths:
        # rename to alphanumeric names
        cmd += ' | sed -f {}'.format(os.path.basename(sed_script_paths[0]))
    if genome_list_path:
        # note, using mafRowOrderer instead of taffy sort because latter does not keep reference row
        # (also, using taffy in this spot requires roundtrip maf-taf-maf tho that can be refactored away)
        cmd += ' | {} mafRowOrderer -m - --order-file {}{} 2> {}.sort.time'.format(time_cmd, os.path.basename(genome_list_path), time_end, chunk_num)
    if sed_script_paths:
        # rename to original, since it needs to be compatible with hal in next step
        cmd += ' | sed -f {}'.format(os.path.basename(sed_script_paths[1]))
    cmd += ' | {} taffy view {} 2> {}.m2t.time'.format(time_cmd, time_end, chunk_num)        
    norm_opts = ''
    if options.maximumBlockLengthToMerge is not None:
        norm_opts += '-m {}'.format(options.maximumBlockLengthToMerge)
    if options.maximumGapLength is not None:
        norm_opts += ' -n {}'.format(options.maximumGapLength)
    if options.filterGapCausingDupes:
        norm_opts += ' -d '
    if options.fractionSharedRows is not None:
        norm_opts += '-q {}'.format(options.fractionSharedRows)
    cmd += ' | {} taffy norm -a {} -k {}{} 2> {}.tn.time'.format(time_cmd, hal_path, norm_opts, time_end, chunk_num)
    if sed_script_paths:
        # rename to alphanumeric names
        cmd += ' | sed -f {}'.format(os.path.basename(sed_script_paths[0]))    
    if options.maxRefNFrac:
        cmd += ' | mafFilter -m - -N {}'.format(options.maxRefNFrac)
    # get rid of single-row (ie ref-only) blcks while we're filtering
    cmd += ' | mafFilter -m - -l 2'
    if options.dupeMode == 'single':
        cmd += ' | mafDuplicateFilter -m - -k'                                               
    elif options.dupeMode == 'consensus':
        cmd += ' | maf_stream merge_dups consensus'
        # need to resort after merge_dups
        if genome_list_path:
            cmd += ' | mafRowOrderer -m - --order-file {}'.format(os.path.basename(genome_list_path))
    if sed_script_paths:
        #rename back to original names
        cmd += ' | sed -f {}'.format(os.path.basename(sed_script_paths[1]))
    if chunk_num != 0:
        cmd += ' | grep -v ^#'
    if options.outputMAF.endswith('.gz'):
        cmd += ' | bgzip'        
    cmd += ' > {}'.format(chunk_name(chunk_num, options, tag='.norm'))
    cmd += ' && mv {} {}'.format(chunk_name(chunk_num, options, tag='.norm'), chunk_name(chunk_num, options))
    return cmd
    
def read_time_mem(timefile_path):
    ''' return a string with some resource usage from the given logfile '''
    mem_usage = '?'
    time_usage = '?'
    do_mem = os.environ.get("CACTUS_LOG_MEMORY")

    with open(timefile_path, 'r') as timefile:
        if do_mem:
            # parse /usr/bin/time -v default output kind of like common.py        
            for line in timefile:
                if 'Maximum resident set size (kbytes):' in line:
                    mem_usage = bytes2human(int(line.split()[-1]) * 1024)
                elif 'Elapsed (wall clock) time (h:mm:ss or m:ss):' in line:
                    time_toks = line.split()[-1].split(':')
                    if len(time_toks) == 2:
                        time_usage = 60 * int(time_toks[0]) + float(time_toks[1])
                    else:
                        time_usage = 3600 * int(time_toks[0]) + 60 * int(time_toks[1]) + float(time_toks[2])
        else:
            time_usage = timefile.readline().strip().split()[1]

    msg = 'in time: {}'.format(time_usage)
    if do_mem:
        msg += ' and memory: {}'.format(mem_usage)
    return msg

def get_rename_map(genome_list, prefix='g'):
    """ map genome names to new alphanumeric placeholder names """

def get_sed_rename_scripts(work_dir, genome_list, out_bed = False, prefix='g'):
    """ get a pair of sciprts to run before/after mafDuplicateFilter or hgMafSummary to
    make sure all genomes alphanumeric """

    # make a map of genome names to unique alphanumeric genome names
    genome_set = set(genome_list)
    name_map = {}
    counter = 0
    for genome in genome_list:
        if not genome.isalnum():
            assert ' ' not in genome and '\t' not in genome
            new_name = '{}{}'.format(prefix, counter)
            while new_name in genome_set or new_name in name_map:
                counter += 1
                new_name = '{}{}'.format(prefix, counter)
            name_map[genome] = new_name
            counter += 1

    # make a sed script to map from genome to alphanumeric
    if len(name_map) > 0:
        before_script_path = os.path.join(work_dir, 'genome_to_alnum.sed')
        after_script_path = os.path.join(work_dir, 'alnum_to_genome.sed')    
        with open(before_script_path, 'w') as before_script_file, open(after_script_path, 'w') as after_script_file:
            # maftools uses spaces and taffy tabs, we make sure consistent here
            before_script_file.write('s/[ ]\\+/\\t/g\n')
            if not out_bed:
                after_script_file.write('s/[ ]\\+/\\t/g\n')
            for genome in reversed(sorted(name_map.keys())):
                before_script_file.write('s/^s\\t{}\\./s\\t{}\\./\n'.format(genome, name_map[genome]))
                if out_bed:
                    after_script_file.write('s/\\t{}\\t/\\t{}\\t/\n'.format(name_map[genome], genome))
                else:
                    after_script_file.write('s/^s\\t{}\\./s\\t{}\\./\n'.format(name_map[genome], genome))

        return (before_script_path, after_script_path, name_map)

    return None
    

def hal2maf_batch(job, hal_id, batch_chunks, genome_list, options, config):
    """ run hal2maf on a batch of chunks in parallel """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))
    job.fileStore.readGlobalFile(hal_id, hal_path)
    # apply any applicable renaming to the genome names
    sed_script_paths = get_sed_rename_scripts(work_dir, genome_list)
    if sed_script_paths:
        renamed_genome_list = []
        for genome in genome_list:
            renamed_genome = sed_script_paths[2][genome] if genome in sed_script_paths[2] else genome
            renamed_genome_list.append(renamed_genome)
    else:
        renamed_genome_list = genome_list
            
    # serialize the list
    genome_list_path = os.path.join(work_dir, 'genome.list')
    with open(genome_list_path, 'w') as genome_list_file:
        for genome in renamed_genome_list:
            genome_list_file.write('{}\n'.format(genome))
        
    h2m_cmds = [hal2maf_cmd(hal_path, chunk, i, options, config) for i, chunk in enumerate(batch_chunks)]

    # do it with parallel
    h2m_cmd_path = os.path.join(work_dir, 'hal2maf_cmds.txt')
    with open(h2m_cmd_path, 'w') as h2m_cmd_file:
        for cmd in h2m_cmds:
            h2m_cmd_file.write(cmd + '\n')
    parallel_cmd = [['cat', h2m_cmd_path],
                    ['parallel', '-j', str(options.batchParallelHal2maf), '{}']]
    RealtimeLogger.info('First of {} commands in parallel batch: {}'.format(len(h2m_cmds), h2m_cmds[0]))
    try:
        cactus_call(parameters=parallel_cmd, work_dir=work_dir)
    except Exception as e:
        logger.error("Parallel hal2maf command failed, dumping all stderr")
        for chunk_num in range(len(batch_chunks)):
            for tag, cmd in [('h2m', 'hal2maf')]:
                stderr_file_path = os.path.join(work_dir, '{}.{}.time'.format(chunk_num, tag))
                if os.path.isfile(stderr_file_path):
                    with open(stderr_file_path, 'r') as stderr_file:
                        for line in stderr_file:
                            logger.error(line.rstrip())
        raise e

    # realtime log the running times in format similar to cactus_call
    # (but breakign up the big piped mess into something more readable)
    for chunk_num in range(len(batch_chunks)):
        chunk_msg = "Successfully ran "
        cmd_toks = h2m_cmds[chunk_num].split(' ')
        for i in range(len(cmd_toks)):
            cmd_toks[i] = cmd_toks[i].rstrip(')')
        for tag, cmd in [('h2m', 'hal2maf')]:
            tag_start = cmd_toks.index(cmd)
            tag_end = cmd_toks.index('2>')        
            tag_cmd = ' '.join(cmd_toks[tag_start:tag_end])
            chunk_msg += "{}{} {} ".format('' if cmd == 'hal2maf' else 'and ', tag_cmd, read_time_mem(os.path.join(work_dir, '{}.{}.time'.format(chunk_num, tag))))
            cmd_toks = cmd_toks[tag_end + 1:]
        RealtimeLogger.info(chunk_msg)

    if not options.raw:
        # and a separate batch for the taf commands. potentially much less efficient as all will block on the
        # slowest h2m command, but we need to be able to specify fewer cores due to higher memory usage
        taf_cmds = [taf_cmd(hal_path, chunk, i, genome_list_path, sed_script_paths, options) for i, chunk in enumerate(batch_chunks)]

        # do it with parallel
        taf_cmd_path = os.path.join(work_dir, 'taf_cmds.txt')
        with open(taf_cmd_path, 'w') as taf_cmd_file:
            for cmd in taf_cmds:
                taf_cmd_file.write(cmd + '\n')
        parallel_cmd = [['cat', taf_cmd_path],
                        ['parallel', '-j', str(options.batchParallelTaf), '{}']]
        RealtimeLogger.info('First of {} commands in parallel batch: {}'.format(len(taf_cmds), taf_cmds[0]))
        try:
            cactus_call(parameters=parallel_cmd, work_dir=work_dir)
        except Exception as e:
            logger.error("Parallel taffy command failed, dumping all stderr")
            for chunk_num in range(len(batch_chunks)):
                for tag, cmd in [('m2t', 'view'), ('sort', 'sort'), ('tn', 'norm')]:                
                    stderr_file_path = os.path.join(work_dir, '{}.{}.time'.format(chunk_num, tag))
                    if os.path.isfile(stderr_file_path):
                        with open(stderr_file_path, 'r') as stderr_file:
                            for line in stderr_file:
                                logger.error(line.rstrip())
            raise e

        # realtime log the running times in format similar to cactus_call
        # (but breakign up the big piped mess into something more readable)
        for chunk_num in range(len(batch_chunks)):
            chunk_msg = "Successfully ran "
            cmd_toks = taf_cmds[chunk_num].split(' ')
            for i in range(len(cmd_toks)):
                cmd_toks[i] = cmd_toks[i].rstrip(')')
            for tag, cmd in [('m2t', 'view'), ('sort', 'sort'), ('tn', 'norm')]:
                tag_start = None
                if tag == 'sort' and 'mafRowOrderer' in cmd_toks:
                    tag_start = cmd_toks.index('mafRowOrderer')
                    tag_end = tag_start + cmd_toks[tag_start:].index('2>')
                    tag_cmd = 'mafRowOrderer'
                elif cmd in cmd_toks:
                    tag_start = cmd_toks.index(cmd) - 1
                    tag_end = tag_start + cmd_toks[tag_start:].index('2>')        
                    tag_cmd = ' '.join(cmd_toks[tag_start:tag_end])
                if tag_start is not None:
                    chunk_msg += "{}{} {} ".format('' if cmd == 'hal2maf' else 'and ', tag_cmd, read_time_mem(os.path.join(work_dir, '{}.{}.time'.format(chunk_num, tag))))
            RealtimeLogger.info(chunk_msg)

    # merge up the results
    maf_path = os.path.join(work_dir, 'out.' + os.path.basename(options.outputMAF)).replace('.taf', '.maf')
    catFiles([os.path.join(work_dir, '{}'.format(chunk_name(i, options))) for i in range(len(batch_chunks))], maf_path)

    return job.fileStore.writeGlobalFile(maf_path)
            

def hal2maf_merge(job, maf_ids, options):
    """ just cat the results """
    work_dir = job.fileStore.getLocalTempDir()
    merged_path = os.path.join(work_dir, 'merged.' + os.path.basename(options.outputMAF)).replace('.taf', '.maf')
    in_paths = []
    for i, maf_id in enumerate(maf_ids):
        in_paths.append(job.fileStore.readGlobalFile(maf_id))
    catFiles(in_paths, merged_path)
    # todo: figure out why converting to taf at the chunk level doesn't work (result can't be parsed)
    if options.outputMAF.endswith('.taf') or options.outputMAF.endswith('.taf.gz'):
        taf_path = os.path.join(work_dir, os.path.basename(options.outputMAF))
        view_cmd = ['taffy', 'view', '-u', '-i', merged_path]
        if taf_path.endswith('.gz'):
            view_cmd += ['-c']
        cactus_call(parameters=view_cmd, outfile=taf_path)
        merged_path = taf_path
        
    return job.fileStore.writeGlobalFile(merged_path)

def taffy_index(job, maf_id, output_path):
    """ run taffy index """
    work_dir = job.fileStore.getLocalTempDir()
    maf_path = os.path.join(work_dir, os.path.basename(output_path))
    job.fileStore.readGlobalFile(maf_id, maf_path)
    
    cactus_call(parameters=['taffy', 'index', '-i', maf_path])

    return job.fileStore.writeGlobalFile(maf_path + '.tai')
    
def taffy_coverage(job, maf_id, genome_list, options):
    """ compute coverage with taffy coverage """
    work_dir = job.fileStore.getLocalTempDir()
    maf_path = os.path.join(work_dir, os.path.basename(options.outputMAF))
    job.fileStore.readGlobalFile(maf_id, maf_path)
    
    cmd = ['taffy', 'coverage', '-r', options.refGenome, '-g', ' '.join(genome_list)]
    if options.coverageSexChroms:
        for contig in options.coverageSexChroms:
            cmd += ['-s', '{}.{}'.format(options.refGenome, contig)]
    if options.coverageGapThresholds:
        for gap in options.coverageGapThresholds:
            cmd += ['-a', str(gap)]

    if maf_path.endswith('.maf') or maf_path.endswith('.maf.gz'):
        cmd = [['taffy', 'view', '-i', maf_path], cmd]
    else:
        cmd += ['-i', maf_path]

    tsv_path = maf_path + '.cov.tsv'
    cactus_call(parameters=cmd, outfile=tsv_path)

    return job.fileStore.writeGlobalFile(tsv_path)
