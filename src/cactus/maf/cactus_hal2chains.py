#!/usr/bin/env python3

"""
This is a script to convert a HAL file into a set of pairwise alignment in UCSC chains format.
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
from cactus.shared.version import cactus_commit
from cactus.progressive.cactus_prepare import human2bytesN
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.paf.paf import get_distances

from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count
from cactus.shared.common import cactus_clamp_memory
from toil.lib.humanize import bytes2human
from sonLib.bioio import getTempDirectory
from sonLib.nxnewick import NXNewick
from sonLib.bioio import newickTreeParser

def main():
    parser = Job.Runner.getDefaultArgumentParser()

    parser.add_argument("halFile", help = "HAL file to convert to MAF")
    parser.add_argument("outDir", help = "Output directory")
    parser.add_argument("--targetGenomes", nargs='*',
                        help="name(s) of target genomes (all leaves if empty).",
                        default=[])
    parser.add_argument("--queryGenomes", nargs='*',
                        help="name(s) of query genomes (all leaves if empty).",
                        default=[])
    parser.add_argument("--includeSelfAlignments", action="store_true",
                        help="include genome-vs-self chains (which are trivial one-block alignments)",
                        default=False)
    parser.add_argument("--bigChain", action="store_true",
                        help="output bigChain.bb and bigChain.link.bb files as well",
                        default=False)
    parser.add_argument("--useHalSynteny",
                        help="use halSynteny instead of halLiftover. halLiftover is default because that's what CAT uses",
                        action="store_true")
    parser.add_argument("--maxAnchorDistance", type=int,
                        help="set --maxAnchorDistance in halSynteny")
    parser.add_argument("--minBlockSize", type=int,
                        help="set --minBlockSize in halSynteny")
    parser.add_argument("--linearGapThreshold", type=float, default=1.0,
                        help="Use axtChain -linearGap=medium if the tree distance between query and target is less than this value, and -linearGap=loose otherwise. UCSC's guidance (from doBlastzChainNet.pl) is to use medium for pairs within the same taxonomic class (e.g. mammal-mammal) and loose across class boundaries (e.g. mammal-bird) [default=1.0]")

    parser.add_argument("--inMemory",
                       help="use --inMemory for halLiftover and/or halSynteny",
                       action="store_true")

    # batching options (mirror cactus-hal2maf): each batch is a single toil job that reads the hal
    # once and runs multiple hal2chains invocations in parallel via GNU parallel. this dramatically
    # reduces the number of hal copies when running on a cluster.
    parser.add_argument("--batchSize", type=int, help="Number of (query,target) pairs in each hal2chains batch", default=None)
    parser.add_argument("--batchCount", type=int, help="Number of hal2chains batches [default 1 unless --batchSize set]", default=None)
    parser.add_argument("--batchCores", type=int, help="Number of cores for each hal2chains batch.")
    parser.add_argument("--batchMemory", type=human2bytesN,
                        help="Memory in bytes for each hal2chains batch (defaults to an estimate based on the input hal size). "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes)", default=None)
    parser.add_argument("--batchParallelHal2chains", type=int,
                        help="Number of hal2chains commands to execute in parallel within a batch. Use to throttle down concurrency "
                        "to save memory. [default=batchCores]", default=None)

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

    if options.maxAnchorDistance and not options.useHalSynteny:
        raise RuntimeError('--maxAnchorDistance can only be used with --useHalSynteny')
    if options.minBlockSize and not options.useHalSynteny:
        raise RuntimeError('--minBlockSize can only be used with --useHalSynteny')
    if options.includeSelfAlignments and options.useHalSynteny:
        raise RuntimeError('--includeSelfAlignments cannot be used with --useHalSynteny')

    if options.batchSize and options.batchCount:
        raise RuntimeError('Only one of --batchSize and --batchCount can be specified')

    # batchCores default: full cpu count on single_machine, required otherwise
    if options.batchCores is None:
        if options.batchSystem.lower() in ['single_machine', 'singlemachine']:
            options.batchCores = cactus_cpu_count()
            if options.maxCores:
                options.batchCores = min(options.batchCores, int(options.maxCores))
            logger.info('Setting batchCores to {}'.format(options.batchCores))
        else:
            raise RuntimeError('--batchCores must be specified for batch systems other than singleMachine')

    if not options.batchParallelHal2chains:
        options.batchParallelHal2chains = options.batchCores
    if options.batchParallelHal2chains > options.batchCores:
        raise RuntimeError('--batchParallelHal2chains cannot exceed the number of batch cores ({})'.format(options.batchCores))

    if not options.batchCount and not options.batchSize:
        logger.info('Using default batch count of 1')
        options.batchCount = 1

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    # make sure we have an outdir
    if not options.outDir.startswith('s3://') and not os.path.isdir(options.outDir):
        os.makedirs(options.outDir)

    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            chains_id_dict = toil.restart()
        else:
            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals(options)
            
            hal_id = toil.importFile(options.halFile)            
            chains_id_dict = toil.start(Job.wrapJobFn(hal2chains_workflow, config, options, hal_id))

        #export the chains
        for query_genome in chains_id_dict.keys():
            clean_query_name = query_genome.replace('#', '.').replace(' ', '.')
            for target_genome in chains_id_dict[query_genome].keys():                
                clean_target_name = target_genome.replace('#', '.').replace(' ', '.')
                chain_id = chains_id_dict[query_genome][target_genome]['chains']
                out_base = makeURL(os.path.join(options.outDir, clean_target_name + '_vs_' + clean_query_name))
                toil.exportFile(chain_id, out_base + ".chain.gz")
                if options.bigChain:
                    bigchain_id = chains_id_dict[query_genome][target_genome]['bigChain']
                    toil.exportFile(bigchain_id, out_base + ".bigChain.bb")
                    biglink_id = chains_id_dict[query_genome][target_genome]['bigLink']
                    toil.exportFile(biglink_id, out_base + ".bigChain.link.bb")

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-hal2chains has finished after {} seconds".format(run_time))


def hal2chains_workflow(job, config, options, hal_id):
    root_job = Job()
    job.addChild(root_job)
    check_tools_job = root_job.addChildJobFn(hal2chains_check_tools, options)
    get_genomes_job = check_tools_job.addFollowOnJobFn(hal2chains_get_genomes, config, options, hal_id,
                                                       disk=int(hal_id.size * 1.2))
    leaf_genomes = get_genomes_job.rv(0)
    distance_matrix = get_genomes_job.rv(1)
    chrom_info_job = get_genomes_job.addFollowOnJobFn(hal2chains_chrom_info_all, config, options, hal_id, leaf_genomes)
    hal2chains_all_job = chrom_info_job.addFollowOnJobFn(hal2chains_all, config, options, hal_id, chrom_info_job.rv(), distance_matrix)
    return hal2chains_all_job.rv()

def hal2chains_check_tools(job, options):
    """ make sure we have the required ucsc commands available on the PATH"""
    tools = ['axtChain', 'faToTwoBit', 'pslPosTarget']
    if options.bigChain:
        tools += ['hgLoadChain', 'bedToBigBed']
    for tool in tools:
        try:
            cactus_call(parameters=[tool])
        except Exception as e:
            # hacky way to see if the usage info came out of running the command with no arguments
            if "usage:" in str(e):
                continue
        raise RuntimeError('Required tool, {}, not found in PATH. Please download it with wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/{}'.format(tool, tool))

def hal2chains_get_genomes(job, config, options, hal_id):    
    """ get the genomes in the hal file, return it along with a distance matrix between them"""
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))    
    job.fileStore.readGlobalFile(hal_id, hal_path)
    tree_str = cactus_call(parameters=['halStats', hal_path, '--tree'], check_output=True).strip()
    mc_tree = MultiCactusTree(NXNewick().parseString(tree_str, addImpliedRoots=False))
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    leaf_genomes = []
    anc_genomes = []
    for node in mc_tree.preOrderTraversal():
        if mc_tree.getName(node) == graph_event:
            continue
        if mc_tree.isLeaf(node):
            leaf_genomes.append(mc_tree.getName(node))
        else:
            anc_genomes.append(mc_tree.getName(node))

    if options.targetGenomes:
        for input_genome in options.targetGenomes:
            if input_genome not in leaf_genomes and input_genome not in anc_genomes:
                raise RuntimeError('Input --targetGenomes {} not found in HAL file'.format(input_genome))

    if options.queryGenomes:
        for input_genome in options.queryGenomes:
            if input_genome not in leaf_genomes and input_genome not in anc_genomes:
                raise RuntimeError('Input --queryGenomes {} not found in HAL file'.format(input_genome))

    # Distances between all pairs of nodes
    event_tree = newickTreeParser(tree_str)
    distances = get_distances(event_tree)
    distance_matrix = {}
    for k,v in distances.items():
        g1,g2 = k[0].iD, k[1].iD
        if g1 not in distance_matrix:
            distance_matrix[g1] = {}
        distance_matrix[g1][g2] = v
    
    return leaf_genomes, distance_matrix
    
def compute_batches(num_items, options):
    """ figure out the number of batches and items per batch given --batchSize/--batchCount """
    if num_items <= 0:
        return 0, 0
    if options.batchSize:
        num_batches = math.ceil(num_items / options.batchSize)
        batch_size = options.batchSize
    else:
        num_batches = max(1, min(options.batchCount or 1, num_items))
        batch_size = math.ceil(num_items / num_batches)
    return num_batches, batch_size

def lpt_assign(items, num_batches, cost_fn):
    """ Longest-Processing-Time-first greedy assignment: sort items by cost descending, then
    repeatedly place the next item into the batch with the smallest running total. Produces
    a list of num_batches lists of items with roughly equal total cost — good for balancing
    toil-job completion times when per-item cost varies widely. """
    import heapq
    if num_batches <= 0 or not items:
        return []
    scored = sorted(items, key=lambda x: -cost_fn(x))
    heap = [[0.0, i, []] for i in range(num_batches)]
    heapq.heapify(heap)
    for item in scored:
        total, idx, bucket = heapq.heappop(heap)
        bucket.append(item)
        heapq.heappush(heap, [total + cost_fn(item), idx, bucket])
    heap.sort(key=lambda x: x[1])
    return [bucket for _, _, bucket in heap]

def chain_pair_cost(q, t, chrom_info_dict, distance_matrix, epsilon=0.05):
    """ Cost proxy for a hal2chains pair: alignment volume is bounded by the smaller 2bit
    (the column of the alignment) and decreases with tree distance. The small epsilon keeps
    very-close-sister pairs from blowing up (matches the 577-way data, where d<0.01 pairs are
    only somewhat slower than d=0.1-0.3 pairs despite being 30x closer). """
    q_size = chrom_info_dict[q]['2bit'].size
    t_size = chrom_info_dict[t]['2bit'].size
    d = distance_matrix[q][t]
    return min(q_size, t_size) / (d + epsilon)

def estimate_batch_memory(options, hal_id, pair_2bit_max=0):
    """ pick a memory request for a hal2chains batch job.

    Tuning data:
      Per-pair runs on a 577-way HAL (~900 GiB), one pair per toil job:
        target=hg38:     worst pair used 42 GiB (ratio 26x q+t); old formula requested 106 GiB.
        target=chicken:  worst pair used 20 GiB (ratio 10x q+t).
        target=catshark: worst pair used  7.6 GiB (ratio 3.8x q+t).
      Batched run on same HAL with k=3 parallel mammal-to-hg38 pipelines:
        batch peak was 41 GiB — virtually identical to the single-pair peak (~40 GiB).
        Parallel tasks share the HAL page cache; only the per-task axtChain working set adds.

    Model (non-inMemory):
      total ≈ first_task_peak + (k-1) × per_task_axtchain_working_set
      first_task_peak absorbs the shared HAL page cache + one task's full working set.
      We estimate first_task_peak as 30 × pair_2bit_max (covers the 26x observed worst case).
      Additional tasks add ~15 × pair_2bit_max each (the axtChain DP state without the cache).

    inMemory mode: each task loads its own full HAL copy, so all terms scale with k.

    Users can override via --batchMemory. """
    if options.batchMemory:
        return options.batchMemory
    k = options.batchParallelHal2chains
    first_task_peak = 30 * pair_2bit_max  # covers shared hal cache + one task's axtChain
    additional_per_task = 15 * pair_2bit_max  # axtChain working set only (cache already paid for)
    if options.inMemory:
        # each task has its own hal copy — nothing is shared
        total = k * (hal_id.size + 30 * pair_2bit_max)
    else:
        total = first_task_peak + max(0, k - 1) * additional_per_task
    return cactus_clamp_memory(total)

def hal2chains_chrom_info_all(job, config, options, hal_id, genomes):
    """ get the 2bit and chrom sizes of every relevant genome, in batches that each copy the hal once """

    # default query / target genome sets to all leaves if they weren't input
    if not options.queryGenomes:
        options.queryGenomes = genomes
    if not options.targetGenomes:
        options.targetGenomes = genomes

    all_genomes = sorted(set(options.queryGenomes + options.targetGenomes))
    num_batches, batch_size = compute_batches(len(all_genomes), options)

    # no 2bit sizes available yet — the batch's own job is what generates them.
    # chrom_info is lighter than chain-building (halStats + hal2fasta|faToTwoBit), so assume
    # per-task memory is bounded by hal page cache + a little overhead.
    batch_memory = estimate_batch_memory(options, hal_id)

    chrom_info_dict = {}
    for i in range(num_batches):
        batch_genomes = all_genomes[i * batch_size : (i + 1) * batch_size]
        if not batch_genomes:
            continue
        # disk: hal copy + headroom for all the 2bits/beds we'll generate (2bits ~ fasta size ~ hal/ngenomes)
        batch_disk = int(hal_id.size * 1.2) + int(hal_id.size * 1.5 * len(batch_genomes) / max(1, len(all_genomes)))
        batch_job = job.addChildJobFn(hal2chains_chrom_info_batch, config, options, hal_id, batch_genomes,
                                      disk=batch_disk,
                                      cores=options.batchCores,
                                      memory=batch_memory)
        for g in batch_genomes:
            chrom_info_dict[g] = batch_job.rv(g)
    return chrom_info_dict

def hal2chains_chrom_info_batch(job, config, options, hal_id, batch_genomes):
    """ extract chrom info (bed, 2bit, optionally chrom.sizes) for a batch of genomes,
    running each extraction in parallel via GNU parallel after a single hal copy """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))
    job.fileStore.readGlobalFile(hal_id, hal_path)
    hal_base = os.path.basename(hal_path)

    cmds = []
    for g in batch_genomes:
        sizes = '{}.chrom.sizes'.format(g)
        bed = '{}.bed'.format(g)
        tbit = '{}.2bit'.format(g)
        err = '{}.chrominfo.stderr'.format(g)
        cmd = ('set -eo pipefail && '
               "halStats {hal} --chromSizes {g} > {sizes} && "
               "awk '{{print $1\"\\t0\\t\"$2}}' {sizes} > {bed} && "
               "hal2fasta {hal} {g} | faToTwoBit stdin {tbit}"
               ).format(hal=hal_base, g=g, sizes=sizes, bed=bed, tbit=tbit)
        cmd += ' 2> {}'.format(err)
        cmds.append(cmd)

    cmd_path = os.path.join(work_dir, 'chrominfo_cmds.txt')
    with open(cmd_path, 'w') as f:
        for c in cmds:
            f.write(c + '\n')

    parallel_cmd = [['cat', cmd_path],
                    ['parallel', '-j', str(options.batchParallelHal2chains), '{}']]
    RealtimeLogger.info('First of {} commands in chrom_info batch: {}'.format(len(cmds), cmds[0]))
    try:
        cactus_call(parameters=parallel_cmd, work_dir=work_dir)
    except Exception as e:
        logger.error("Parallel chrom_info command failed, dumping all stderr")
        for g in batch_genomes:
            err_path = os.path.join(work_dir, '{}.chrominfo.stderr'.format(g))
            if os.path.isfile(err_path) and os.path.getsize(err_path) > 0:
                logger.error('--- stderr for chrom_info of {} ---'.format(g))
                with open(err_path, 'r') as ef:
                    for line in ef:
                        logger.error(line.rstrip())
        raise

    out = {}
    for g in batch_genomes:
        bed_path = os.path.join(work_dir, '{}.bed'.format(g))
        tbit_path = os.path.join(work_dir, '{}.2bit'.format(g))
        info = {'bed': job.fileStore.writeGlobalFile(bed_path),
                '2bit': job.fileStore.writeGlobalFile(tbit_path)}
        if options.bigChain:
            sizes_path = os.path.join(work_dir, '{}.chrom.sizes'.format(g))
            info['sizes'] = job.fileStore.writeGlobalFile(sizes_path)
        out[g] = info
    return out

def hal2chains_all(job, config, options, hal_id, chrom_info_dict, distance_matrix):
    """ build chains for all (query, target) pairs, in batches that each copy the hal once """

    # default query / target genome sets to all leaves if they weren't input
    if not options.queryGenomes:
        options.queryGenomes = list(chrom_info_dict.keys())
    if not options.targetGenomes:
        options.targetGenomes = list(chrom_info_dict.keys())

    pairs = []
    for q in options.queryGenomes:
        for t in options.targetGenomes:
            if options.includeSelfAlignments or t != q:
                pairs.append((q, t))

    num_batches, batch_size = compute_batches(len(pairs), options)

    # Balance batch completion times: LPT-greedy assignment using a cost proxy that tracks
    # observed runtime on a 577-way vertebrate alignment (r=-0.63 between tree distance and
    # chain pipeline runtime for target=hg38). Also sort the pairs within each batch by cost
    # descending so GNU parallel fills its slots with the longest jobs first — this keeps the
    # tail short (the stragglers finish while later, small jobs are still arriving).
    cost_fn = lambda p: chain_pair_cost(p[0], p[1], chrom_info_dict, distance_matrix)
    batches = lpt_assign(pairs, num_batches, cost_fn)
    for b in batches:
        b.sort(key=lambda p: -cost_fn(p))
    RealtimeLogger.info('Assigned {} pairs to {} batches via LPT (cost = min(q_2bit,t_2bit)/(distance+0.05))'.format(
        len(pairs), len([b for b in batches if b])))

    output_dict = {}
    for batch_pairs in batches:
        if not batch_pairs:
            continue
        batch_genomes = set()
        for q, t in batch_pairs:
            batch_genomes.add(q)
            batch_genomes.add(t)
        batch_chrom_info = {g: chrom_info_dict[g] for g in batch_genomes}
        batch_distances = {(q, t): distance_matrix[q][t] for q, t in batch_pairs}

        # size memory/disk using actual 2bit sizes. worst-case pair: max(q_2bit + t_2bit) over batch_pairs
        pair_2bit_max = max(chrom_info_dict[q]['2bit'].size + chrom_info_dict[t]['2bit'].size
                             for q, t in batch_pairs)
        total_2bit = sum(chrom_info_dict[g]['2bit'].size for g in batch_genomes)
        batch_memory = estimate_batch_memory(options, hal_id, pair_2bit_max=pair_2bit_max)
        # disk: hal copy + 2bits + bed files + per-pair scratch/outputs (~10x 2bit per pair, as in old code)
        batch_disk = int(hal_id.size * 1.2) + total_2bit + len(batch_pairs) * 10 * pair_2bit_max

        batch_job = job.addChildJobFn(hal2chains_batch, config, options, hal_id,
                                      batch_pairs, batch_chrom_info, batch_distances,
                                      disk=batch_disk,
                                      cores=options.batchCores,
                                      memory=batch_memory)

        for q, t in batch_pairs:
            if q not in output_dict:
                output_dict[q] = {}
            output_dict[q][t] = {'chains': batch_job.rv(q, t)}
            if options.bigChain:
                # bigChain doesn't use the HAL — just the chain.gz + chrom.sizes → bigChain.bb.
                # disk: chain file + intermediates (~10x target 2bit is generous).
                # memory: hgLoadChain peak is ~2x uncompressed chain size ≈ a few GiB for big mammals.
                t_2bit_size = chrom_info_dict[t]['2bit'].size
                bigchains_job = batch_job.addFollowOnJobFn(chain2bigchain, options, q, t,
                                                           chrom_info_dict[t], batch_job.rv(q, t),
                                                           disk=max(20 * t_2bit_size, 1024**3),
                                                           memory=cactus_clamp_memory(max(5 * t_2bit_size, 2 * 1024**3)))
                output_dict[q][t]['bigChain'] = bigchains_job.rv(0)
                output_dict[q][t]['bigLink'] = bigchains_job.rv(1)

    return output_dict


def hal2chains_batch(job, config, options, hal_id, batch_pairs, batch_chrom_info, batch_distances):
    """ build chains for a batch of (query, target) pairs: copy hal + 2bit/bed files once,
    then run all halLiftover|...|axtChain pipelines in parallel via GNU parallel """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))
    job.fileStore.readGlobalFile(hal_id, hal_path)
    hal_base = os.path.basename(hal_path)

    # pull down the per-genome chrom-info files (beds and 2bits) — one copy per unique genome
    for g, info in batch_chrom_info.items():
        job.fileStore.readGlobalFile(info['bed'], os.path.join(work_dir, '{}.bed'.format(g)))
        job.fileStore.readGlobalFile(info['2bit'], os.path.join(work_dir, '{}.2bit'.format(g)))

    cmds = []
    out_names = {}
    for q, t in batch_pairs:
        out_name = '{}_vs_{}.chain.gz'.format(t, q)
        err_name = '{}_vs_{}.chain.stderr'.format(t, q)
        out_names[(q, t)] = out_name

        distance = batch_distances[(q, t)]
        gap = 'medium' if distance < options.linearGapThreshold else 'loose'

        if options.useHalSynteny:
            first = 'halSynteny {hal} /dev/stdout --queryGenome {q} --targetGenome {t}'.format(hal=hal_base, q=q, t=t)
            if options.inMemory:
                first += ' --inMemory'
            if options.maxAnchorDistance:
                first += ' --maxAnchorDistance {}'.format(options.maxAnchorDistance)
            if options.minBlockSize:
                first += ' --minBlockSize {}'.format(options.minBlockSize)
        else:
            first = 'halLiftover {hal} {q} {q}.bed {t} /dev/stdout --outPSL'.format(hal=hal_base, q=q, t=t)
            if options.inMemory:
                first += ' --inMemory'

        cmd = ('set -eo pipefail && '
               '{first} | pslPosTarget /dev/stdin /dev/stdout | '
               'axtChain -psl -verbose=0 -linearGap={gap} /dev/stdin {t}.2bit {q}.2bit /dev/stdout | '
               'gzip > {out}').format(first=first, gap=gap, t=t, q=q, out=out_name)
        cmd += ' 2> {}'.format(err_name)
        cmds.append(cmd)

    cmd_path = os.path.join(work_dir, 'chain_cmds.txt')
    with open(cmd_path, 'w') as f:
        for c in cmds:
            f.write(c + '\n')

    parallel_cmd = [['cat', cmd_path],
                    ['parallel', '-j', str(options.batchParallelHal2chains), '{}']]
    RealtimeLogger.info('First of {} commands in chain batch: {}'.format(len(cmds), cmds[0]))
    try:
        cactus_call(parameters=parallel_cmd, work_dir=work_dir)
    except Exception as e:
        logger.error("Parallel hal2chains command failed, dumping all stderr")
        for q, t in batch_pairs:
            err_path = os.path.join(work_dir, '{}_vs_{}.chain.stderr'.format(t, q))
            if os.path.isfile(err_path) and os.path.getsize(err_path) > 0:
                logger.error('--- stderr for {} vs {} ---'.format(q, t))
                with open(err_path, 'r') as ef:
                    for line in ef:
                        logger.error(line.rstrip())
        raise

    out = {}
    for q, t in batch_pairs:
        chain_path = os.path.join(work_dir, out_names[(q, t)])
        if q not in out:
            out[q] = {}
        out[q][t] = job.fileStore.writeGlobalFile(chain_path)
    return out

# https://raw.githubusercontent.com/ucscGenomeBrowser/kent/cb3cd0e9c5b9006b7d316df7905faa516f880c6b/src/hg/lib/bigChain.as
bigChain_as = \
'''table bigChain
"bigChain pairwise alignment"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name or ID of item, ideally both human readable and unique"
    uint score;         "Score (0-1000)"
    char[1] strand;     "+ or - for strand"
    uint tSize;         "size of target sequence"
    string qName;       "name of query sequence"
    uint qSize;         "size of query sequence"
    uint qStart;        "start of alignment on query sequence"
    uint qEnd;          "end of alignment on query sequence"
    double chainScore;    "score from chain"
    )
'''

# https://genome.ucsc.edu/goldenPath/help/examples/bigLink.as
bigLink_as = \
'''
table bigLink
"bigLink pairwise alignment"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name or ID of item, ideally both human readable and unique"
    uint qStart;        "start of alignment on query sequence"
    )
'''

def chain2bigchain(job, options, query_genome, target_genome, target_info, chain_id):
    """ convert the chain to big chain. from https://genome.ucsc.edu/goldenPath/help/bigChain.html """
    work_dir = job.fileStore.getLocalTempDir()
    chain_path = os.path.join(work_dir, target_genome + '_vs_' + query_genome + '.chain.gz')
    job.fileStore.readGlobalFile(chain_id, chain_path)
    target_sizes_path = os.path.join(work_dir, target_genome + '.chrom.sizes')
    job.fileStore.readGlobalFile(target_info['sizes'], target_sizes_path)

    # unzip the chains
    uz_chain_path = chain_path[:-3]
    cactus_call(parameters=['gzip', '-dc', chain_path], outfile=uz_chain_path)

    # make the bigChain.as file
    bigchain_as_path = os.path.join(work_dir, 'bigChain.as')
    with open(bigchain_as_path, 'w') as as_file:
        as_file.write(bigChain_as)

    # Use the hgLoadChain utility to generate the chain.tab and link.tab files needed to create the bigChain file:
    cactus_call(parameters=['hgLoadChain', '-noBin', '-test', target_genome, 'bigChain', uz_chain_path])

    # Create the bigChain file from your input chain file using a combination of sed, awk and the bedToBigBed utility:
    bigchain_path = os.path.join(work_dir, target_genome + '_vs_' + query_genome + '.bigChain')
    cactus_call(parameters=[['sed', 's/\\.000000//', 'chain.tab'],
                            ['awk', 'BEGIN {OFS=\"\\t\"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}']],
                outfile=bigchain_path)

    bigchain_bb_path = bigchain_path + '.bb'
    cactus_call(parameters=['bedToBigBed', '-type=bed6+6', '-as={}'.format(bigchain_as_path),
                            '-tab', bigchain_path, target_sizes_path, bigchain_bb_path])

    # make the bigLink.as file
    biglink_as_path = os.path.join(work_dir, 'bigLink.as')
    with open(biglink_as_path, 'w') as as_file:
        as_file.write(bigLink_as)

    # To display your date in the Genome Browser, you must also create a binary indexed link file to accompany your bigChain file:
    biglink_path = bigchain_path + '.link'
    cactus_call(parameters=[['awk', 'BEGIN {OFS=\"\\t\"} {print $1, $2, $3, $5, $4}', 'link.tab'],
                            ['sort', '-k1,1', '-k2,2n']],
                outfile=biglink_path)

    biglink_bb_path = biglink_path + '.bb'
    cactus_call(parameters=['bedToBigBed', '-type=bed4+1', 'as={}'.format(biglink_as_path),
                            '-tab', biglink_path, target_sizes_path, biglink_bb_path])

    return (job.fileStore.writeGlobalFile(bigchain_bb_path), job.fileStore.writeGlobalFile(biglink_bb_path))
    
         


