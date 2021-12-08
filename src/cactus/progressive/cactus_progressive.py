#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

"""Wrapper to run the cactus_workflow progressively, using the input species tree as a guide

tree.
"""

import logging
import os
import xml.etree.ElementTree as ET
import timeit
from argparse import ArgumentParser
from base64 import b64encode

from toil.lib.bioio import getTempFile
from toil.lib.threading import cpu_count
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from toil.job import Job
from toil.common import Toil

from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import findRequiredNode
from cactus.shared.common import makeURL
from cactus.shared.common import catFiles
from cactus.shared.common import cactus_call
from cactus.shared.version import cactus_commit
from cactus.shared.common import cactusRootPath
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import write_s3

from cactus.pipeline.cactus_workflow import cactus_cons_with_resources
from cactus.progressive.progressive_decomposition import compute_outgroups, compute_schedule, parse_seqfile, get_subtree
from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor
from cactus.preprocessor.dnabrnnMasking import loadDnaBrnnModel
from cactus.paf.local_alignment import make_paf_alignments
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.schedule import Schedule
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import setupBinaries, importSingularityImage

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory

def logAssemblyStats(job, message, name, sequenceID, preemptable=True):
    sequenceFile = job.fileStore.readGlobalFile(sequenceID)
    analysisString = cactus_call(parameters=["cactus_analyseAssembly", sequenceFile], check_output=True)
    job.fileStore.logToMaster("%s, got assembly stats for genome %s: %s" % (message, name, analysisString))

def preprocess_all(job, options, config_node, input_seq_id_map):
    ''' run prepreprocessor on every input sequence '''
    root_job = Job()
    job.addChild(root_job)
    events = list(input_seq_id_map.keys())
    seq_ids = list(input_seq_id_map.values())
    preprocessor_job = root_job.addChild(CactusPreprocessor(seq_ids, config_node, eventNames=events))
    pp_seq_ids = {}
    for i, event in enumerate(events):
        pp_seq_ids[event] = preprocessor_job.rv(i)

    # do the logging and checkpointing
    root_job.addFollowOnJobFn(save_preprocessed_files, options, config_node, input_seq_id_map)
    
    return pp_seq_ids

def save_preprocessed_files(job, options, config_node, seq_id_map):
    ''' checkpoint the the intermediate url; also log the stats '''
    # Save preprocessed sequences
    if options.intermediateResultsUrl is not None:
        for genome, seqID in list(seq_id_map.items()):
            job.fileStore.exportFile(seqID, options.intermediateResultsUrl + '-preprocessed-' + genome)

    # Log the stats for the preprocessed assemblies
    for name, sequence in list(seq_id_map.items()):
        job.addChildJobFn(logAssemblyStats, "After preprocessing", name, sequence)

def progressive_step(job, options, config_node, seq_id_map, tree, schedule, og_map, event):
    ''' use the schedule to run the events in order '''

    RealtimeLogger.info("STEP {}".format(event))
    
    root_job = Job()
    job.addChild(root_job)

    # each results are dicst that map event name -> consolidated results
    results_list = []
    if event in schedule.depTree:
        for dep in schedule.deps(event):
            child_job = root_job.addChildJobFn(progressive_step, options, config_node, seq_id_map, tree, schedule, og_map, dep)
            results_list.append(child_job.rv())

    next_job = root_job.addFollowOnJobFn(progressive_next, options, config_node, seq_id_map, tree, schedule, og_map, event, results_list)

    return next_job.rv()

def progressive_next(job, options, config_node, seq_id_map, tree, schedule, og_map, event, results_lists):
    ''' compute the alignment and make the hal '''

    RealtimeLogger.info("NEXT {}".format(event))
    results_lists = flatten_lists(results_lists)
    results_dict = results_list_to_dict(results_lists)
    for ev, res in results_dict.items():
        seq_id_map[ev] = res['fa']

    # get our subtree
    subtree = get_subtree(tree, event, ConfigWrapper(config_node), og_map)

    # get our events
    subtree_eventmap = {}
    for leaf in subtree.getLeaves():
        subtree_eventmap[subtree.getName(leaf)] = seq_id_map[subtree.getName(leaf)]
    
    # do the blast
    paf_job = job.addChildJobFn(make_paf_alignments, NXNewick().writeString(subtree), subtree_eventmap, event, config_node)

    # do the consolidated
    consolidated_job = paf_job.addFollowOnJobFn(cactus_cons_with_resources, subtree, event, config_node, subtree_eventmap, og_map, paf_job.rv(),
                                                cons_cores = options.consCores, intermediate_results_url = options.intermediateResultsUrl)

    return results_lists + [consolidated_job.rv()]

def flatten_lists(results_list):
    ''' make sure our list is just 1d '''
    def flatten(l):
        flat_list = []
        did_something = False
        for x in l:
            if isinstance(x, list):
                did_something = True
                for y in x:
                    flat_list.append(y)
            else:
                flat_list.append(x)
        return (flat_list, did_something)

    while True:
        (results_list, did_something) = flatten(results_list)
        if not did_something:
            break
    return results_list

def results_list_to_dict(results_list):
    ''' convert results list of tuples into a dcit'''
    results_dict = {}
    for res in results_list:
        assert isinstance(res, tuple)
        # note, this has to be consistent with interface in cactus_cons() in cactus_workflow.py
        results_dict[res[0]] = { 'c2h' : res[1], 'c2h.fa' : res[2], 'fa' : res[3] }
        
    return results_dict

def export_hal(job, tree, config_node, seq_id_map, og_map, results, event=None, cacheBytes=None,
               cacheMDC=None, cacheRDC=None, cacheW0=None, chunk=None, deflate=None, inMemory=True,
               checkpointInfo=None, acyclicEvent=None):

    # todo: going through list nonsense because (i think) it helps with promises, should at least clean up
    results = results_list_to_dict(flatten_lists(results))
    print(results)
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, '{}.hal'.format(event if event else tree.getRootName()))

    # make the multicactus tree
    # todo: use throughout?
    tree = MultiCactusTree(tree)
    tree.nameUnlabeledInternalNodes(ConfigWrapper(config_node).getDefaultInternalNodePrefix())
    tree.computeSubtreeRoots()
    subtree_roots = set(tree.getSubtreeRootNames())

    # find subtree if event specified
    root_node = None
    if event is not None:
        assert event in tree.nameToId and not tree.isLeaf(tree.nameToId[event])
        root_node = tree.nameToId[event]

    hal_path = os.path.join(work_dir, '{}.hal'.format(event if event else tree.getRootName()))

    for node in tree.breadthFirstTraversal(root_node):
        genome_name = tree.getName(node)
        if genome_name in subtree_roots:
            outgroups = og_map[genome_name] if genome_name in og_map else []
            subtree = get_subtree(tree, genome_name, ConfigWrapper(config_node), og_map, include_outgroups=False)
            tree_string = NXNewick().writeString(subtree)
            sub_hal_path = os.path.join(work_dir, '{}.hal.c2h'.format(genome_name))
            hal_fasta_path = os.path.join(work_dir, '{}.hal.fa'.format(genome_name))
            assert genome_name in results
            job.fileStore.readGlobalFile(results[genome_name]['c2h'], sub_hal_path)
            job.fileStore.readGlobalFile(results[genome_name]['c2h.fa'], hal_fasta_path)

            args = ['halAppendCactusSubtree', sub_hal_path, hal_fasta_path, tree_string, hal_path]
            
            if len(outgroups) > 0:
                args += ["--outgroups", ",".join(outgroups)]
            if cacheBytes is not None:
                args += ["--cacheBytes", cacheBytes]
            if cacheMDC is not None:
                args += ["--cacheMDC", cacheMDC]
            if cacheRDC is not None:
                args += ["--cacheRDC", cacheRDC]
            if cacheW0 is not None:
                args += ["--cacheW0", cacheW0]
            if chunk is not None:
                args += ["--chunk", chunk]
            if deflate is not None:
                args += ["--deflate", deflate]
            if inMemory is True:
                args += ["--inMemory"]

            cactus_call(parameters=args)

    cactus_call(parameters=["halSetMetadata", hal_path, "CACTUS_COMMIT", cactus_commit])
    config_path = os.path.join(work_dir, 'config.xml')
    ConfigWrapper(config_node).writeXML(config_path)
    with open(config_path, 'rb') as configFile:
        cactus_call(parameters=["halSetMetadata", hal_path, "CACTUS_CONFIG", b64encode(configFile.read()).decode()])

    if acyclicEvent:
        cactus_call(parameters=["halRemoveDupes", hal_path, acyclicEvent])

    if checkpointInfo:
        write_s3(hal_path, checkpointInfo[1], region=checkpointInfo[0])

    return job.fileStore.writeGlobalFile(hal_path)
    
def progressive_workflow(job, options, config_node, tree, schedule, og_map, input_seq_id_map):
    ''' run the entire progressive workflow '''

    # start with the preprocessor
    if not options.skipPreprocessor:
        pp_job = job.addChildJobFn(preprocess_all, options, config_node, input_seq_id_map)
        seq_id_map = pp_job.rv()
    else:
        pp_job = Job()
        job.addChild(pp_job)
        seq_id_map = input_seq_id_map

    # then do the progressive workflow
    root_event = options.root if options.root else tree.getRootName()
    progressive_job = pp_job.addFollowOnJobFn(progressive_step, options, config_node, seq_id_map, tree, schedule, og_map, root_event)

    # then do the hal export
    hal_export_job = progressive_job.addFollowOnJobFn(export_hal, tree, config_node, seq_id_map, og_map, progressive_job.rv(), event=root_event,
                                                      disk=ConfigWrapper(config_node).getExportHalDisk(), preemptable=False)

    return hal_export_job.rv()

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("outputHal", type=str, help = "Output HAL file")

    #Progressive Cactus Options
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--root", dest="root", help="Name of ancestral node (which"
                      " must appear in NEWICK tree in <seqfile>) to use as a "
                      "root for the alignment.  Any genomes not below this node "
                      "in the tree may be used as outgroups but will never appear"
                      " in the output.  If no root is specifed then the root"
                      " of the tree is used. ", default=None)
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)
    parser.add_argument("--gpu", action="store_true",
                        help="Enable GPU acceleration by using Segaling instead of lastz")
    parser.add_argument("--consCores", type=int, 
                        help="Number of cores for each cactus_consolidated job (defaults to all available / maxCores on single_machine)", default=None)
    parser.add_argument("--intermediateResultsUrl",
                        help="URL prefix to save intermediate results like DB dumps to (e.g. "
                        "prefix-dump-caf, prefix-dump-avg, etc.)", default=None)
    parser.add_argument("--skipPreprocessor",
                        help="Do not run any preprocessor jobs",
                        action="store_true")

    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    # Try to juggle --maxCores and --consCores to give some reasonable defaults where possible
    if options.batchSystem.lower() in ['single_machine', 'singlemachine']:
        if options.maxCores is not None:
            if int(options.maxCores) <= 0:
                raise RuntimeError('Cactus requires --maxCores >= 1')
        if options.consCores is None:
            if options.maxCores is not None:
                options.consCores = int(options.maxCores)
            else:
                options.consCores = cpu_count()
        elif options.maxCores is not None and options.consCores > int(options.maxCores):
            raise RuntimeError('--consCores must be <= --maxCores')
    else:
        if not options.consCores:
            raise RuntimeError('--consCores required for non single_machine batch systems')
    if options.maxCores is not None and options.consCores > int(options.maxCores):
        raise RuntimeError('--consCores must be <= --maxCores')

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    start_time = timeit.default_timer()

    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            halID = toil.restart()
        else:

            # load up the seqfile and figure out the outgroups and schedule
            config_node = ET.parse(options.configFile).getroot()
            ConfigWrapper(config_node).substituteAllPredefinedConstantsWithLiterals()
            tree, input_seq_map, og_candidates = parse_seqfile(options.seqFile)
            og_map = compute_outgroups(tree, ConfigWrapper(config_node), set(og_candidates), options.root)
            schedule = compute_schedule(tree, ConfigWrapper(config_node), og_map)
            event_set = set([tree.getName(node) for node in tree.postOrderTraversal(options.root)]).union(set(og_map.keys()))

            #hack
            for i in ["Anc0", "Anc1", "Anc2", "mr"]:
                sb = get_subtree(tree, i, ConfigWrapper(config_node), og_map, False)
                print("SUB", i, NXNewick().writeString(sb))
                        
            #import the sequences
            input_seq_id_map = {}
            for (genome, seq) in input_seq_map.items():
                if genome in event_set:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    logger.info("Importing {}".format(seq))
                    input_seq_id_map[genome] = toil.importFile(seq)
                
            # Make sure we have the dna-brnn model in the filestore if we need it
            loadDnaBrnnModel(toil, config_node)

            # run the whole workflow
            hal_id = toil.start(Job.wrapJobFn(progressive_workflow, options, config_node, tree, schedule, og_map, input_seq_id_map))

        toil.exportFile(hal_id, makeURL(options.outputHal))
    
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("Cactus has finished after {} seconds".format(run_time))

if __name__ == '__main__':
    main()
