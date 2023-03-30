#!/usr/bin/env python3

"""
build a minigraph in Toil, using a cactus seqfile as input
"""

import os, sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit, time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
from cactus.refmap.cactus_graphmap_join import unzip_seqfile
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count
from cactus.progressive.multiCactusTree import MultiCactusTree
from sonLib.bioio import getTempDirectory

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("seqFile", help = "Seq file (will be modified if necessary to include graph Fasta sequence)")
    parser.add_argument("outputGFA", help = "Output Minigraph GFA")
    parser.add_argument("--reference", type=str, required=True,
                        help = "Reference genome name (added to minigraph first). Order in seqfile used otherwise")
    parser.add_argument("--mapCores", type=int, help = "Number of cores for minigraph.  Overrides graphmap cpu in configuration")
    
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

    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            gfa_id = toil.restart()
        else:
            # load up the seqfile
            config_node = ET.parse(options.configFile).getroot()
            config_wrapper = ConfigWrapper(config_node)
            config_wrapper.substituteAllPredefinedConstantsWithLiterals()
            graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")

            # load the seqfile
            seqFile = SeqFile(options.seqFile)
            input_seq_map = seqFile.pathMap
            
            if options.reference not in input_seq_map:
                raise RuntimeError("Specified reference not in seqfile")
            for sample in input_seq_map.keys():
                if sample != options.reference and sample.startswith(options.reference):
                    raise RuntimeError("Input sample {} is prefixed by given reference {}. ".format(sample, options.reference) +                        
                                       "This is not supported by this version of Cactus, " +
                                       "so one of these samples needs to be renamed to continue")

            # apply cpu override                
            if options.mapCores is not None:
                findRequiredNode(config_node, "graphmap").attrib["cpu"] = str(options.mapCores)
            mg_cores = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "cpu", typeFn=int, default=1)
            if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
                mg_cores = min(mg_cores, cactus_cpu_count(), int(options.maxCores) if options.maxCores else sys.maxsize)
                findRequiredNode(config_node, "graphmap").attrib["cpu"] = str(mg_cores)
            
            #import the sequences
            input_seq_id_map = {}
            input_seq_order = [options.reference]
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
                    if genome != options.reference:
                        input_seq_order.append(genome)
            
            gfa_id = toil.start(Job.wrapJobFn(minigraph_construct_workflow, config_node, input_seq_id_map, input_seq_order, options.outputGFA))

        #export the gfa
        toil.exportFile(gfa_id, makeURL(options.outputGFA))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-minigraph has finished after {} seconds".format(run_time))
        
def minigraph_construct_workflow(job, config_node, seq_id_map, seq_order, gfa_path, sanitize=True):
    """ minigraph can handle bgzipped files but not gzipped; so unzip everything in case before running"""
    if sanitize:
        sanitize_job = job.addChildJobFn(sanitize_fasta_headers, seq_id_map)
        sanitized_seq_id_map = sanitize_job.rv()
    else:
        sanitized_seq_id_map = seq_id_map
        sanitize_job = Job()
        job.addChild(sanitize_job)
    xml_node = findRequiredNode(config_node, "graphmap")
    mg_cores = getOptionalAttrib(xml_node, "cpu", typeFn=int, default=1)
    sort_type = getOptionalAttrib(xml_node, "minigraphSortInput", str, default=None)
    if sort_type == "mash":
        sort_job = sanitize_job.addFollowOnJobFn(sort_minigraph_input_with_mash, sanitized_seq_id_map, seq_order)
        seq_order = sort_job.rv()
        prev_job = sort_job
    else:
        prev_job = sanitize_job
    minigraph_job = prev_job.addFollowOnJobFn(minigraph_construct, config_node, sanitized_seq_id_map, seq_order, gfa_path,
                                              cores = mg_cores,
                                              disk = 5 * sum([seq_id.size for seq_id in seq_id_map.values()]))
    return minigraph_job.rv()

def sort_minigraph_input_with_mash(job, seq_id_map, seq_order):
    """ Sort the input """
    # (dist, length) pairs which will be sorted decreasing on dist, breaking ties with increasing on length
    root_job = Job()
    job.addChild(root_job)
    # assumption : reference is first
    mash_dists = [(0, sys.maxsize)]
    for seq in seq_order[1:]:
        dist = root_job.addChildJobFn(mash_dist, seq, seq_order[0], seq_id_map,
                                      disk = seq_id_map[seq].size + seq_id_map[seq_order[0]].size).rv()
        mash_dists.append(dist)                    

    return root_job.addFollowOnJobFn(mash_distance_order, seq_order, mash_dists).rv()

def mash_dist(job, query_seq, ref_seq, seq_id_map):
    """ get the mash distance """
    work_dir = job.fileStore.getLocalTempDir()
    ref_path = os.path.join(ref_seq + '.fa')
    query_path = os.path.join(query_seq + '.fa')
    job.fileStore.readGlobalFile(seq_id_map[ref_seq], ref_path)
    job.fileStore.readGlobalFile(seq_id_map[query_seq], query_path)

    mash_output = cactus_call(parameters=['mash', 'dist', query_path, ref_path], check_output=True)
    dist = float(mash_output.strip().split()[2])
    size = sum([len(r.seq) for r in SeqIO.parse(query_path, 'fasta')])
    RealtimeLogger.info('mash distance of {} (size = {}) to reference {} = {}'.format(query_seq, size, ref_seq, dist))
    return (dist, size)
    
def mash_distance_order(job, seq_order, mash_dists):
    """ get the sequence order from the mash distance"""
    seq_to_dist = { }
    # we want to sort reverse on size, so make them negative
    # and we round the dists to ignore tiny changes
    mash_dists = [(round(x, 4),-y) for x,y in mash_dists]
    for seq, md in zip(seq_order, mash_dists):
        seq_to_dist[seq] = md
    return sorted(seq_order, key = lambda x : seq_to_dist[x])
    
def minigraph_construct(job, config_node, seq_id_map, seq_order, gfa_path):
    """ Make minigraph """
    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, os.path.basename(gfa_path))

    # parse options from the config
    xml_node = findRequiredNode(config_node, "graphmap")
    minigraph_opts = getOptionalAttrib(xml_node, "minigraphConstructOptions", str, default="")     
    opts_list = minigraph_opts.split()
    if '-t' not in opts_list:
        opts_list += ['-t', str(job.cores)]
    
    # download the sequences
    local_fa_paths = {}
    for event, fa_id in seq_id_map.items():
        fa_id = seq_id_map[event]
        fa_path = os.path.join(work_dir, '{}.fa'.format(event))
        job.fileStore.readGlobalFile(fa_id, fa_path)
        local_fa_paths[event] = fa_path
        assert os.path.getsize(local_fa_paths[event]) > 0

    sort_type = getOptionalAttrib(xml_node, "minigraphSortInput", str, default=None)
    if sort_type == "size":
        # don't touch the reference
        sorted_order = [seq_order[0]]
        # sort the rest by size, biggest first, which is usually how we help poa
        seq_to_size = {}
        for seq in seq_order[1:]:
            seq_to_size[seq] = sum([len(r.seq) for r in SeqIO.parse(local_fa_paths[seq], 'fasta')])
        sorted_order += sorted(seq_order[1:], key = lambda e : seq_to_size[e], reverse=True)
        seq_order = sorted_order

    mg_cmd = ['minigraph'] + opts_list
    for event in seq_order:
        mg_cmd += [os.path.basename(local_fa_paths[event])]

    if gfa_path.endswith('.gz'):
        mg_cmd = [mg_cmd, ['bgzip', '--threads', str(job.cores)]]

    cactus_call(parameters=mg_cmd, outfile=gfa_path, work_dir=work_dir, realtimeStderrPrefix='[minigraph]')

    return job.fileStore.writeGlobalFile(gfa_path)
        
        
        

    
