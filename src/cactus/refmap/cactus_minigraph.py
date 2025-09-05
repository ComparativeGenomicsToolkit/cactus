#!/usr/bin/env python3

"""
build a minigraph in Toil, using a cactus seqfile as input
"""

import os, sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit, time
import math
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from operator import itemgetter
import gzip

from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles, write_s3
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import clean_jobstore_files
from cactus.shared.version import cactus_commit
from cactus.progressive.cactus_prepare import human2bytesN
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers
from cactus.paf.last_scoring import last_train
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from toil.lib.conversions import bytes2human
from cactus.shared.common import cactus_cpu_count
from cactus.shared.common import cactus_clamp_memory
from cactus.progressive.multiCactusTree import MultiCactusTree
from sonLib.bioio import getTempDirectory, getTempFile

def main():
    parser = Job.Runner.getDefaultArgumentParser()

    parser.add_argument("seqFile", help = "Seq file (or chromfile with --batch)")
    parser.add_argument("outputGFA", help = "Output Minigraph GFA (or directory in --batch mode)")
    parser.add_argument("--reference", required=True, nargs='+', type=str,
                        help = "Reference genome name(s) (added to minigraph first). Mash distance to 1st reference to determine order of other genomes (use minigraphSortInput in the config xml to toggle this behavior).")
    parser.add_argument("--mgCores", type=int, help = "Number of cores for minigraph construction (defaults to the same as --maxCores).")
    parser.add_argument("--mgMemory", type=human2bytesN,
                        help="Memory in bytes for the minigraph construction job (defaults to an estimate based on the input data size). "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))", default=None)
    parser.add_argument("--lastTrain", action="store_true",
                        help="Use last-train to estimate scoring matrix from input data", default=False)
    parser.add_argument("--refOnly", action="store_true",
                        help="Only build the graph out of reference genome(s). Can be used when it will only be used for chromosome-splitting, for example")
    parser.add_argument("--batch", action="store_true",
                        help="Run independently on set of chromosomea inputs (chromfile as from cactus-graphmap-split). Note that the output will be a directory and not a GFA")
        
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

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    if options.batch:
        # the output gfa is a directory, make sure it's there
        if not os.path.isdir(options.outputGFA):
            os.makedirs(options.outputGFA)

    # map chrom name to seqFile
    input_seqfiles = {}
    if options.batch:
        input_seqfiles = read_chromfile(options.seqFile)
    else:
        input_seqfiles['all'] = options.seqFile
            
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            output_dict = toil.restart()
        else:
            # load up the config
            config_node = ET.parse(options.configFile).getroot()
            config_wrapper = ConfigWrapper(config_node)
            config_wrapper.substituteAllPredefinedConstantsWithLiterals(options)

            # apply cpu override
            if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
                if not options.mgCores:
                    options.mgCores = sys.maxsize
                options.mgCores = min(options.mgCores, cactus_cpu_count(), int(options.maxCores) if options.maxCores else sys.maxsize)
            else:
                if not options.mgCores:
                    raise RuntimeError("--mgCores required run *not* running on single machine batch system")

            if '://' not in options.outputGFA:
                options.outputGFA = os.path.abspath(options.outputGFA)

            # maps name -> input_seq_id_map, input_seq_order
            input_dict = minigraph_construct_import_sequences(options, config_wrapper, input_seqfiles, toil)
                
            # output_dict:  chrom-> (gfa_id, pansn_gfa_id, train_id)
            output_dict = toil.start(Job.wrapJobFn(minigraph_construct_batch_workflow, options, config_node, input_dict, options.outputGFA))

        export_minigraph_construct_output(options, input_seqfiles, output_dict, toil)
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-minigraph has finished after {} seconds".format(run_time))

def read_chromfile(chromfile_path):
    """ read tsv into map """
    # map chrom name to seqFile
    chromfile = {}
    with open(chromfile_path, 'r') as chrom_file:
        for line in chrom_file:
            toks = line.strip().split()
            if len(toks):
                assert len(toks) >= 2
                assert toks[0] not in chromfile
                chromfile[toks[0]] = toks[1:]
    return chromfile

def minigraph_construct_import_sequences(options, config_wrapper, input_seqfiles, file_store):
    """ import the files and return a map """
    # maps name -> input_seq_id_map, input_seq_order
    input_dict = {}
    config_node = config_wrapper.xmlRoot
    graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")

    # load the seqfiles
    for chrom, seqfile_path in input_seqfiles.items():
        if type(seqfile_path) is list:
            seqfile_path = seqfile_path[0]
        seqFile = SeqFile(seqfile_path, defaultBranchLen=config_wrapper.getDefaultBranchLen(pangenome=True))
        input_seq_map = seqFile.pathMap
        raw_input_seq_order = seqFile.seqOrder

        # make sure the reference is first
        input_seq_order = [options.reference[0]]
        for seq in raw_input_seq_order:
            if seq != options.reference[0]:
                input_seq_order.append(seq)

        # hack out everything but reference
        if options.refOnly:
            input_seq_order = options.reference
            ref_seq_map = {}
            for sample in options.reference:
                ref_seq_map[sample] = input_seq_map[sample]
            input_seq_map = ref_seq_map

        # validate the sample names
        check_sample_names(input_seq_map.keys(), options.reference)

        #import the sequences
        input_seq_id_map = {}
        leaves = set([seqFile.tree.getName(node) for node in seqFile.tree.getLeaves()])
        for (genome, seq) in input_seq_map.items():
            if genome != graph_event and genome in leaves:                
                if os.path.isdir(seq):
                    tmpSeq = getTempFile()
                    catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                    seq = tmpSeq
                seq = makeURL(seq)
                input_seq_id_map[genome] = file_store.importFile(seq)
            elif genome in input_seq_order:
                input_seq_order.remove(genome)

        input_dict[chrom] = (input_seq_id_map, input_seq_order)
        
    return input_dict

def export_minigraph_construct_output(options, input_seqfiles, output_dict, toil):
    if options.batch:
        chrom_file_path = os.path.join(options.outputGFA, 'chromfile.mg.txt')
        if chrom_file_path.startswith('s3://'):
            chrom_file_temp_path = getTempFile()
        else:
            chrom_file_temp_path = chrom_file_path                    
        chromfile = open(chrom_file_temp_path, 'w')
    for chrom, output_ids in output_dict.items():
        gfa_id, pansn_gfa_id, train_id = output_ids
        if options.batch:
            gfa_path = os.path.join(options.outputGFA, chrom + '.sv.gfa.gz')
        else:
            gfa_path = options.outputGFA
        if train_id:
            train_path = gfa_path.replace('.gfa.gz', '.gfa').replace('.gfa', '.train')
        else:
            train_path = None
        #export the gfa            
        toil.exportFile(pansn_gfa_id, makeURL(gfa_path))
        if train_path:
            # export the scoring model (.train)
            toil.exportFile(train_id, makeURL(train_path))
        if options.batch:
            chromfile.write('{}\t{}\t{}\t{}\n'.format(chrom, input_seqfiles[chrom][0], gfa_path,
                                                      train_path if train_path else '*'))
    if options.batch:
        chromfile.close()
        if chrom_file_path.startswith('s3://'):
            write_s3(chrom_file_temp_path, chrom_file_path)

def check_sample_names(sample_names, references):
    """ make sure we have a workable set of sample names """

    # make sure we have the reference
    if references:
        assert type(references) in [list, str]
        if type(references) is str:
            references = [references]
        if references[0] not in sample_names:
            raise RuntimeError("Specified reference, \"{}\" not in seqfile".format(references[0]))
        for reference in references:
            # graphmap-join uses reference names as prefixes, so make sure we don't get into trouble with that
            reference_base = os.path.splitext(reference)[0]
            for sample in sample_names:
                sample_base = os.path.splitext(sample)[0]
                if sample != reference and sample_base.startswith(reference_base):
                    raise RuntimeError("Input sample {} is prefixed by given reference {}. ".format(sample_base, reference_base) +    
                                       "This is not supported by this version of Cactus, " +
                                       "so one of these samples needs to be renamed to continue")

    # the "." character is overloaded to specify haplotype, make sure that it makes sense
    for sample in sample_names:
        sample_base, sample_ext = os.path.splitext(sample)
        if not sample_base or (not sample_ext and sample_base.startswith(".")):
            raise RuntimeError("Sample name {} invalid because it begins with \".\"".format(sample))
        if sample_ext and (len(sample_ext) == 1 or not sample_ext[1:].isnumeric()):
            raise RuntimeError("Sample name {} with \"{}\" suffix is not supported. You must either remove this suffix or use .N where N is an integer to specify haplotype".format(sample, sample_ext))

def minigraph_construct_batch_workflow(job, options, config_node, input_dict, gfa_path, sanitize=True):
    """ run the construction workflow on individual chromosomes """
    output_dict = {}
    for chrom, input_info in input_dict.items():
        seq_id_map, seq_order = input_info
        if options.batch:
            gfa_path = os.path.join(options.outputGFA, '{}.gfa.gz'.format(chrom))
        else:
            gfa_path = options.outputGFA
        mgwf_job = job.addChildJobFn(minigraph_construct_workflow, options, config_node, seq_id_map, seq_order, gfa_path, sanitize)
        output_dict[chrom] = mgwf_job.rv()
    return output_dict
                                    
def minigraph_construct_workflow(job, options, config_node, seq_id_map, seq_order, gfa_path, sanitize=True):
    """ minigraph can handle bgzipped files but not gzipped; so unzip everything in case before running"""
    assert type(options.reference) is list
    assert options.reference[0] == seq_order[0]
    if options.refOnly:
        refonly_seq_id_map = {}
        refonly_seq_order = []
        for seq in seq_order:
            if seq in options.reference:
                refonly_seq_order.append(seq)
                refonly_seq_id_map[seq] = seq_id_map[seq]
        seq_id_map, seq_order = refonly_seq_id_map, refonly_seq_order
    ref_size = seq_id_map[options.reference[0]].size
    if sanitize:
        sanitize_job = job.addChildJobFn(sanitize_fasta_headers, seq_id_map, pangenome=True)
        sanitized_seq_id_map = sanitize_job.rv()
    else:
        sanitized_seq_id_map = seq_id_map
        sanitize_job = Job()
        job.addChild(sanitize_job)
    xml_node = findRequiredNode(config_node, "graphmap")
    sort_type = getOptionalAttrib(xml_node, "minigraphSortInput", str, default=None)
    if sort_type == "mash" and len(seq_id_map) > 2:
        sort_job = sanitize_job.addFollowOnJobFn(sort_minigraph_input_with_mash, options, config_node, sanitized_seq_id_map, seq_order)
        seq_order = sort_job.rv()
        prev_job = sort_job
    else:
        prev_job = sanitize_job
    minigraph_job = prev_job.addFollowOnJobFn(minigraph_construct_in_batches, options, config_node, sanitized_seq_id_map, seq_order, gfa_path)
    if options.lastTrain:
        last_train_job = prev_job.addFollowOnJobFn(last_train, config_node, seq_order, sanitized_seq_id_map, 
                                                   cores=options.mgCores,
                                                   disk=8*ref_size,
                                                   memory=cactus_clamp_memory(8*ref_size))
        
    return minigraph_job.rv(0), minigraph_job.rv(1), last_train_job.rv() if options.lastTrain else None

def sort_minigraph_input_with_mash(job, options, config_node, seq_id_map, seq_order):
    """ Sort the input """
    # (dist, length) pairs which will be sorted decreasing on dist, breaking ties with increasing on length
    # assumption : reference is first
    mash_dists = [(0, sys.maxsize)]
    # start by sketching the reference to avoid a bunch of recomputation
    sketch_job = job.addChildJobFn(mash_sketch, seq_order[0], seq_id_map,
                                   disk = seq_id_map[seq_order[0]].size * 2)
    ref_sketch_id = sketch_job.rv()

    dist_root_job = Job()
    sketch_job.addFollowOn(dist_root_job)

    xml_node = findRequiredNode(config_node, "graphmap")
    sort_by_sample = getOptionalAttrib(xml_node, "minigraphSortBySample", str, default='0')
    sort_by_sample = True if sort_by_sample == '1' or (sort_by_sample == 'nonbatch' and not options.batch) else False
    # in the case of diploid inputs, the distance will be driven by the sex of the haplotype
    # ie in human, the haplotype with chrX is always going to be much closer to the reference
    # than the one with chrY, which makes the resulting order pretty meaningless.
    #
    # so to get around this, we group the haplotypes by sample, and compute their distances together
    # so, for example HG002.1 and HG002.2 would be always be consecutive in the order and their distance
    # will be determined jointly (by just concatenating them).
    seq_by_sample = defaultdict(list)
    # note: seq_order[0] is (first) reference, so we don't include it
    for seq_name in seq_order[1:]:
        sample_name = seq_name[:seq_name.rfind('.')] if '.' in seq_name and sort_by_sample else seq_name
        seq_by_sample[sample_name].append(seq_name)

    # list of dictionary (promises) that map genome name to mash distance output
    dist_maps = []
    for sample, names in seq_by_sample.items():
        dist_map = dist_root_job.addChildJobFn(mash_dist, names, seq_order[0], seq_id_map, ref_sketch_id,
                                               disk = 2 * sum(seq_id_map[x].size for x in names) + seq_id_map[seq_order[0]].size).rv()
        dist_maps.append(dist_map)
            
    return dist_root_job.addFollowOnJobFn(mash_distance_order, options, config_node, seq_order, dist_maps).rv()

def mash_sketch(job, ref_seq, seq_id_map):
    """ get the sketch """
    work_dir = job.fileStore.getLocalTempDir()
    ref_path = os.path.join(ref_seq + '.fa')
    job.fileStore.readGlobalFile(seq_id_map[ref_seq], ref_path)

    cactus_call(parameters=['mash', 'sketch', ref_path])

    return job.fileStore.writeGlobalFile(ref_path + '.msh')
    
def mash_dist(job, query_seqs, ref_seq, seq_id_map, ref_sketch_id):
    """ get the mash distance
    returns a map from genome name to -> (sample distance, distance, size)
    where sample_distance is the concatentation of all sequences from the same sample (ie HG002.1 and HG002.2)
    """
    work_dir = job.fileStore.getLocalTempDir()
    ref_sketch_path = os.path.join(ref_seq + '.fa.msh')
    query_paths = [os.path.join(query_seq + '.fa') for query_seq in query_seqs]
    job.fileStore.readGlobalFile(ref_sketch_id, ref_sketch_path)
    for query_seq, query_path in zip(query_seqs, query_paths):
        job.fileStore.readGlobalFile(seq_id_map[query_seq], query_path)

    def parse_mash_output(mash_output):
        return float(mash_output.strip().split()[2])

    output_dist_map = {}
    
    # make the concatenated distance
    cat_mash_dist = None
    if len(query_seqs) > 1:        
        cat_path = query_paths[0] + '.cat'
        catFiles(query_paths, cat_path)
        cat_mash_output = cactus_call(parameters=['mash', 'dist', cat_path, ref_sketch_path], check_output=True)
        cat_mash_dist = parse_mash_output(cat_mash_output)
        # we want samples to stay together in the event of ties (which happen).  So add a a little bit to make unique
        cat_mash_dist += float(sorted(seq_id_map.keys()).index(query_seqs[0])) * sys.float_info.epsilon

    # make the individual distance
    for query_seq, query_path in zip(query_seqs, query_paths):
        mash_output = cactus_call(parameters=['mash', 'dist', query_path, ref_sketch_path], check_output=True)
        dist = parse_mash_output(mash_output)
        sample_dist = cat_mash_dist if cat_mash_dist is not None else dist
        size = sum([len(r.seq) for r in SeqIO.parse(query_path, 'fasta')])
        output_dist_map[query_seq] = (sample_dist, dist, size)
        log_msg = 'mash distance of {} (size = {}) to reference {} = {}'.format(query_seq, size, ref_seq, dist)
        if cat_mash_dist is not None:
            log_msg += ' sample distance = {}'.format(cat_mash_dist)
        RealtimeLogger.info(log_msg)
        
    return output_dist_map
    
def mash_distance_order(job, options, config_node, seq_order, mash_output_maps):
    """ get the sequence order from the mash distance"""

    # we first orient the list of dicts along seq_order
    mash_output_map = mash_output_maps[0] if mash_output_maps else {}
    for output_map in mash_output_maps[1:]:
        mash_output_map.update(output_map)    
    mash_dists = []
    for seq in seq_order:
        if seq in mash_output_map:
            mash_dists.append(mash_output_map[seq])
        else:
            assert seq == options.reference[0]

    # we want to sort reverse on size, so make them negative
    mash_dists = [(x, y, -z) for x,y,z in mash_dists]
    seq_to_dist = {}
    for seq, md in zip(seq_order[1:], mash_dists):
        seq_to_dist[seq] = md

    # sanity check
    max_seq_dist = max(seq_to_dist.items(), key = lambda x : x[1][1])
    if max_seq_dist[1][1] > 0.02:
        job.fileStore.logToMaster('\n\nWARNING: Sample {} has mash distance {} from the reference. A value this high likely means your data is too diverse to construct a useful pangenome graph from.\n'.format(max_seq_dist[0], max_seq_dist[1][1]))

    # sort by mash distance
    mash_order = [seq_order[0]] + sorted(seq_order[1:], key = lambda x : seq_to_dist[x])
        
    # optionally fix secondary references back to their original positions in the seqfile
    sort_refs  = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "minigraphSortReference", typeFn=bool, default=True)
    if not sort_refs and len(options.reference) > 1:
        fixed_order = copy.deepcopy(seq_order)
        empty_slots = []
        for i, seq in enumerate(seq_order):
            if seq not in options.reference:
                empty_slots.append(i)
        j = 0
        for seq in mash_order:
            if seq not in options.reference:
                fixed_order[empty_slots[j]] = seq
                j += 1
        for ref in options.reference[1:]:
            mash_pos = mash_order.index(ref)
            fix_pos = fixed_order.index(ref)
            assert fix_pos == seq_order.index(ref)
            if mash_pos != fix_pos:
                RealtimeLogger.info('Secondary reference {}, which would have mash rank {}, fixed at input rank {} because minigraphSortReference is disabled'.format(ref, mash_pos, fix_pos))
        mash_order = fixed_order

    return mash_order
            
def minigraph_construct_in_batches(job, options, config_node, seq_id_map, seq_order, gfa_path):
    """ Make minigraph in sequential batches"""

    max_size = max([x.size for x in seq_id_map.values()])
    total_size = sum([x.size for x in seq_id_map.values()])
    disk = total_size * 2
    mem = cactus_clamp_memory(60 * max_size + int(total_size / 4))
    if options.mgMemory is not None:
        RealtimeLogger.info('Overriding minigraph_construct memory estimate of {} with {} value {} from --mgMemory'.format(bytes2human(mem), 'greater' if options.mgMemory > mem else 'lesser', bytes2human(options.mgMemory)))     
        mem = options.mgMemory

    # parse options from the config
    xml_node = findRequiredNode(config_node, "graphmap")
    max_batch_size = getOptionalAttrib(xml_node, "minigraphConstructBatchSize", int, default=-1)
    assert max_batch_size > 0
    num_batches = int(math.ceil(len(seq_order) / max_batch_size))
    if num_batches >= 990:
        # Toil will fail with a python recursion error if we try to chain too many jobs
        new_max_batch_size = int(len(seq_order) / 990 + 1)
        job.fileStore.logToMaster('WARNING: Increasing minigraphConstructBatchSize from {} to {} to avoid Toil error from chaining too many jobs'.format(max_batch_size, new_max_batch_size))
        max_batch_size = new_max_batch_size
        num_batches = int(math.ceil(len(seq_order) / max_batch_size))
        assert num_batches > 0 and num_batches <= 991
    prev_job = None
    prev_gfa_path = None
    for i in range(num_batches):        
        batch_size = len(seq_order) - i * max_batch_size if i == num_batches - 1 else max_batch_size
        input_seq_order = seq_order[i * max_batch_size : (i * max_batch_size) + batch_size]
        assert input_seq_order and len(input_seq_order) <= max_batch_size
        out_gfa_path = gfa_path
        pan_sn_output = True
        if i < num_batches - 1:
            if out_gfa_path.endswith('.gz'):
                out_gfa_path = '{}.{}.gz'.format(gfa_path[:-3], i)
            else:
                out_gfa_path = '{}.{}'.format(gfa_path, i)
            pan_sn_output = False
        minigraph_job = Job.wrapJobFn(minigraph_construct, options, config_node, seq_id_map, input_seq_order, out_gfa_path,
                                      prev_job.rv() if prev_job else None, prev_gfa_path,
                                      pan_sn_output,
                                      disk=disk, memory=mem, cores=options.mgCores)
        if prev_job:
            prev_job.addFollowOn(minigraph_job)
            # delete the output of the previous batch from the job store            
            minigraph_job.addFollowOnJobFn(clean_jobstore_files, file_ids=[prev_job.rv()])
        else:
            job.addChild(minigraph_job)
        prev_job = minigraph_job
        prev_gfa_path = out_gfa_path

    return prev_job.rv()

def minigraph_construct(job, options, config_node, seq_id_map, seq_order, gfa_path, prev_gfa_id, prev_gfa_path, pan_sn_output):
    """ Make minigraph """

    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, os.path.basename(gfa_path))
    if prev_gfa_id:
        prev_gfa_path = os.path.join(work_dir, os.path.basename(prev_gfa_path))
        job.fileStore.readGlobalFile(prev_gfa_id, prev_gfa_path)

    # parse options from the config
    xml_node = findRequiredNode(config_node, "graphmap")
    minigraph_opts = getOptionalAttrib(xml_node, "minigraphConstructOptions", str, default="")     
    opts_list = minigraph_opts.split()
    if '-t' not in opts_list:
        opts_list += ['-t', str(job.cores)]
    
    # download the sequences
    local_fa_paths = {}
    for event in seq_order:
        fa_id = seq_id_map[event]
        fa_path = os.path.join(work_dir, '{}.fa'.format(event))
        job.fileStore.readGlobalFile(fa_id, fa_path)
        local_fa_paths[event] = fa_path
        assert os.path.getsize(local_fa_paths[event]) > 0

    mg_cmd = ['minigraph'] + opts_list
    if prev_gfa_id:
        mg_cmd += [os.path.basename(prev_gfa_path)]
    for event in seq_order:
        mg_cmd += [os.path.basename(local_fa_paths[event])]

    if gfa_path.endswith('.gz'):
        mg_cmd = [mg_cmd, ['bgzip', '--threads', str(job.cores)]]

    if options.batch:
        prefix = '[minigraph-{}]'.format(os.path.basename(gfa_path).replace('.gz', '').replace('.gfa', ''))
    else:
        prefix = '[minigraph]'
    cactus_call(parameters=mg_cmd, outfile=gfa_path, work_dir=work_dir, realtimeStderrPrefix=prefix, job_memory=job.memory)

    gfa_out_id = job.fileStore.writeGlobalFile(gfa_path)
    if pan_sn_output:
        # rename to pan-sn before serializing, so it's more useful (ie for anything except cactus)
        pansn_gfa_path = os.path.join(work_dir, 'pan-sn.' + os.path.basename(gfa_path))
        minigraph_gfa_to_pansn(set(seq_id_map.keys()), gfa_path, pansn_gfa_path)
        pansn_gfa_out_id = job.fileStore.writeGlobalFile(pansn_gfa_path)
        return gfa_out_id, pansn_gfa_out_id
    else:
        return gfa_out_id

def minigraph_gfa_to_pansn(names, gfa_path, out_gfa_path):
    """ hack to convert cactus names like id=simChimp.0|simChimp.chr6 to PanSN simChimp#0#simpChimp.chr6
    so that minigraph GFA file can be used outside of catus

    todo: Cactus should probably be changed to just use PanSN internally as well, but that's a much
    bigger lift
    """
    if gfa_path.endswith('.gz'):
        in_file = gzip.open(gfa_path, 'rb')
        out_file = gzip.open(out_gfa_path, 'wb')
    else:
        in_file = open(gfa_path, 'rb')
        out_file = open(out_gfa_path, 'wb')

    for line in in_file:
        line = line.decode()
        if line.startswith('S'):
            toks = line.strip().split('\t')
            for i, tok in enumerate(toks[4:]):
                if tok.startswith('SN:Z:id='):
                    barpos = tok.find('|')
                    assert barpos > 8
                    name = tok[8:barpos]
                    assert name in names
                    dotpos = name.rfind('.')
                    if dotpos > 0:
                        hap = name[dotpos+1:]
                        name = name[:dotpos]
                    else:
                        # add #0 to names without haplotype to make valid PAN-SN
                        hap = '0'
                    toks[4+i] = 'SN:Z:{}#{}#{}'.format(name, hap, tok[barpos+1:])
                    break
            out_file.write(('\t'.join(toks) + '\n').encode())
        else:
            out_file.write(line.encode())

    in_file.close()
    out_file.close()

def minigraph_gfa_from_pansn(job, names, gfa_path, gfa_id):
    """ hack to convert PanSN names like simChimp#0#simpChimp.chr6 to Cactus names like id=simChimp.0|simChimp.chr6
    so that a minigrpah GFA (as converted panSN by minigraph_gfa_to_pansn() above) can be read back into Cactus

    todo: Cactus should probably be changed to just use PanSN internally as well, but that's a much
    bigger lift
    """
    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, os.path.basename(gfa_path))
    job.fileStore.readGlobalFile(gfa_id, gfa_path)
    out_gfa_path = os.path.join(work_dir, 'cactus.' + os.path.basename(gfa_path))
    
    if gfa_path.endswith('.gz'):
        in_file = gzip.open(gfa_path, 'rb')
        out_file = gzip.open(out_gfa_path, 'wb')
    else:
        in_file = open(gfa_path, 'rb')
        out_file = open(out_gfa_path, 'wb')

    for line in in_file:
        line = line.decode()
        if line.startswith('S'):
            toks = line.strip().split('\t')
            for i, tok in enumerate(toks[4:]):
                if tok.startswith('SN:Z:'):
                    hashpos = tok.find('#')
                    if hashpos < 0:
                        # no prefix found: do nothing and hope for the best
                        continue
                    hashpos2 = hashpos + 1 + tok[hashpos+1:].find('#')
                    if hashpos2 < 0:
                        name = tok[5:]
                        hap = "0"
                    else:
                        name = tok[5:hashpos]
                        hap = tok[hashpos+1:hashpos2]
                    # minigraph_to_pansn() will add #0 to names without any dots
                    # we untangle that here using the names list                    
                    if name not in names:
                        name = '{}.{}'.format(name, hap)
                        assert name in names                    
                    toks[4+i] = 'SN:Z:id={}|{}'.format(name, tok[hashpos2+1:])
                    break
            out_file.write(('\t'.join(toks) + '\n').encode())
        else:
            out_file.write(line.encode())

    in_file.close()
    out_file.close()

    return job.fileStore.writeGlobalFile(out_gfa_path)

    
