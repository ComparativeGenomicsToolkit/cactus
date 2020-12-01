#!/usr/bin/env python3

"""jobs and command for running the minigraph pan-genome building pipeline, which takes
   input: usual list of fasta files + a minigraph-compatible GFA graph
   output: PAF file containing alignment of each input Fasta to the contig sequences of the graph
           (these contig sequences are also output in their own Fasta, which can be treated as an "assembly"
            by the rest of cactus)

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
from cactus.pipeline.cactus_workflow import CactusTrimmingBlastPhase
from cactus.shared.common import makeURL
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from toil.job import Job
from toil.common import Toil
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions
from toil.realtimeLogger import RealtimeLogger

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    parser.add_argument("seqFile", help = "Seq file (will be modified if necessary to include graph Fasta sequence)")
    parser.add_argument("minigraphGFA", help = "Minigraph-compatible reference graph in GFA format (can be gzipped)")
    parser.add_argument("outputPAF", type=str, help = "Output pairwise alignment file in PAF format")
    parser.add_argument("--outputFasta", type=str, help = "Output graph sequence file in FASTA format (required if not present in seqFile)")
    parser.add_argument("--ignoreSoftmasked", action="store_true", help = "Ignore softmasked sequence (by hardmasking before sending to minigraph)")

    #WDL hacks
    parser.add_argument("--pathOverrides", nargs="*", help="paths (multiple allowed) to override from seqFile")
    parser.add_argument("--pathOverrideNames", nargs="*", help="names (must be same number as --pathOverrides) of path overrides")

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

    if (options.pathOverrides or options.pathOverrideNames):
        if not options.pathOverrides or not options.pathOverrideNames or \
           len(options.pathOverrideNames) != len(options.pathOverrides):
            raise RuntimeError('same number of values must be passed to --pathOverrides and --pathOverrideNames')

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    start_time = timeit.default_timer()
    runCactusGraphMap(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-graphmap has finished after {} seconds".format(run_time))

def runCactusGraphMap(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            alignmentID = toil.restart()
        else:
            options.cactusDir = getTempDirectory()

            # apply path overrides.  this was necessary for wdl which doesn't take kindly to
            # text files of local paths (ie seqfile).  one way to fix would be to add support
            # for s3 paths and force wdl to use it.  a better way would be a more fundamental
            # interface shift away from files of paths throughout all of cactus
            if options.pathOverrides:
                seqFile = SeqFile(options.seqFile)
                configNode = ET.parse(options.configFile).getroot()
                config = ConfigWrapper(configNode)
                tree = MultiCactusTree(seqFile.tree)
                tree.nameUnlabeledInternalNodes(prefix = config.getDefaultInternalNodePrefix())
                for name, override in zip(options.pathOverrideNames, options.pathOverrides):
                    seqFile.pathMap[name] = override
                override_seq = os.path.join(options.cactusDir, 'seqFile.override')
                with open(override_seq, 'w') as out_sf:
                    out_sf.write(str(seqFile))
                options.seqFile = override_seq

            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals()

            # get the minigraph "virutal" assembly name
            graph_event = getOptionalAttrib(findRequiredNode(configNode, "graphmap"), "assemblyName", default="__MINIGRAPH_SEQUENCES__")

            # load the seqfile
            seqFile = SeqFile(options.seqFile)
            
            logger.info("Genomes for graphmap, {}".format(seqFile.pathMap))

            if not options.outputFasta and graph_event not in seqFile.pathMap:
                raise RuntimeError("{} assembly not found in seqfile so it must be specified with --outputFasta".format(graph_event))

            #import the graph
            gfa_id = toil.importFile(makeURL(options.minigraphGFA))

            #import the sequences (that we need to align for the given event, ie leaves and outgroups)
            seqIDMap = {}
            leaves = set([seqFile.tree.getName(node) for node in seqFile.tree.getLeaves()])
            for genome, seq in seqFile.pathMap.items():
                if genome != graph_event and genome in leaves:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    seqIDMap[genome] = toil.importFile(seq)

            # run the workflow
            paf_id, gfa_fa_id = toil.start(Job.wrapJobFn(minigraph_workflow, options, config, seqIDMap, gfa_id, graph_event))

        #export the paf
        toil.exportFile(paf_id, makeURL(options.outputPAF))
        if gfa_fa_id:
            toil.exportFile(gfa_fa_id, makeURL(options.outputFasta))

        # update the input seqfile (in place!)
        if options.outputFasta:
            add_genome_to_seqfile(options.seqFile, makeURL(options.outputFasta), graph_event)

def minigraph_workflow(job, options, config, seq_id_map, gfa_id, graph_event):
    """ Overall workflow takes command line options and returns (paf-id, (optional) fa-id) """
    fa_id = None
    if options.outputFasta:
        fa_job = job.addChildJobFn(make_minigraph_fasta, gfa_id, graph_event)
        fa_id = fa_job.rv()
    
    paf_job = job.addChildJobFn(minigraph_map_all, config, gfa_id, seq_id_map, options.ignoreSoftmasked)

    return paf_job.rv(), fa_id                                
        
def make_minigraph_fasta(job, gfa_file_id, name):
    """ Use gfatools to make the minigraph "assembly" """

    # note: using the toil-vg convention of naming working files manually so that logging is more readable
    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    fa_path = os.path.join(work_dir, "minigraph_sequences.fa")
    
    job.fileStore.readGlobalFile(gfa_file_id, gfa_path)

    cactus_call(work_dir=work_dir, outfile=fa_path,
                parameters=["gfatools", "gfa2fa", os.path.basename(gfa_path)])

    return job.fileStore.writeGlobalFile(fa_path)

def minigraph_map_all(job, config, gfa_id, fa_id_map, ignore_softmasked):
    """ top-level job to run the minigraph mapping in parallel, returns paf """
    
    # hang everything on this job, to self-contain workflow
    top_job = Job()
    job.addChild(top_job)
    
    # do the mapping
    gaf_ids = []
    for event, fa_id in fa_id_map.items():
        RealtimeLogger.info("adding child event={} faid={} gfaid={}".format(event, fa_id, gfa_id))
        minigraph_map_job = top_job.addChildJobFn(minigraph_map_one, config, event, fa_id, gfa_id,
                                                  ignore_softmasked,
                                                  cores=1, disk=5*(fa_id.size + gfa_id.size))
        gaf_ids.append(minigraph_map_job.rv())

    # convert to paf
    paf_job = top_job.addFollowOnJobFn(merge_gafs_into_paf, config, gaf_ids)

    return paf_job.rv()

def minigraph_map_one(job, config, event_name, fa_file_id, gfa_file_id, ignore_softmasked):
    """ Run minigraph to map a Fasta file to a GFA graph, producing a GAF output """

    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    fa_path = os.path.join(work_dir, "{}.fa".format(event_name))
    gaf_path = os.path.join(work_dir, "{}.gaf".format(event_name))
    
    job.fileStore.readGlobalFile(gfa_file_id, gfa_path)
    job.fileStore.readGlobalFile(fa_file_id, fa_path)

    # parse options from the config
    xml_node = findRequiredNode(config.xmlRoot, "graphmap")
    minigraph_opts = getOptionalAttrib(xml_node, "minigraphMapOptions", str, default="")     
    opts_list = minigraph_opts.split()
    # add required options if not present
    if "-S" not in opts_list:
        opts_list += ["-S"]
    if "--write-mz" not in opts_list:
        opts_list += ["--write-mz"]
    if "-t" not in opts_list:
        opts_list += ["-t", str(int(job.cores))]

    cmd = ["minigraph",
           os.path.basename(gfa_path),
           os.path.basename(fa_path),
           "-o", os.path.basename(gaf_path)] + opts_list

    if ignore_softmasked:
        cmd[2] = '-'
        cmd = [['cactus_softmask2hardmask', os.path.basename(fa_path)], cmd]
    
    # todo: pipe into gzip directly as these files can be huge!!! (requires gzip support be added to mzgaf2paf)
    cactus_call(work_dir=work_dir, parameters=cmd)

    return job.fileStore.writeGlobalFile(gaf_path)

def merge_gafs_into_paf(job, config, gaf_file_ids):
    """ Merge GAF alignments into a single PAF, applying some filters """

    work_dir = job.fileStore.getLocalTempDir()
    paf_path = os.path.join(work_dir, "mz_alignments.paf")
    gaf_paths = []
    for i, gaf_id in enumerate(gaf_file_ids):
        gaf_paths.append("mz_alignment_{}.gaf".format(i))
        job.fileStore.readGlobalFile(gaf_id, os.path.join(work_dir, gaf_paths[-1]))

    xml_node = findRequiredNode(config.xmlRoot, "graphmap")
    mzgaf2paf_opts = []
    mz_filter = getOptionalAttrib(xml_node, "universalMZFilter", float)
    if mz_filter:
        mzgaf2paf_opts += ['-u', str(mz_filter)]
    if getOptionalAttrib(xml_node, "nodeBasedUniversal", typeFn=bool, default=False):
        mzgaf2paf_opts += ['-n']
    if getOptionalAttrib(xml_node, "strictUniversal", typeFn=bool, default=False):
        mzgaf2paf_opts += ['-i']
    min_mz = getOptionalAttrib(xml_node, "minMZBlockLength", int)
    if min_mz:
        mzgaf2paf_opts += ['-m', str(min_mz)]
    mapq = getOptionalAttrib(xml_node, "minMAPQ", int)
    if mapq:
        mzgaf2paf_opts += ['-q', str(mapq)]
    gaf_block = getOptionalAttrib(xml_node, "minGAFBlockLength", int)
    if gaf_block:
        mzgaf2paf_opts += ['-b', str(gaf_block)]
    gaf_node = getOptionalAttrib(xml_node, "minGAFNodeLength", int)
    if gaf_node:
        mzgaf2paf_opts += ['-s', str(gaf_node)]        

    cactus_call(work_dir=work_dir, outfile=paf_path, parameters=["mzgaf2paf"] + gaf_paths + mzgaf2paf_opts)

    # these are big, get rid of them as soon as we can (which is now)
    for gaf_id in gaf_file_ids:
        job.fileStore.deleteGlobalFile(gaf_id)

    return job.fileStore.writeGlobalFile(paf_path)
                                               

def add_genome_to_seqfile(seqfile_path, fasta_path, name):
    """ hack the auto-generated minigraph assembly back into the seqfile for future use """
    seq_file = SeqFile(seqfile_path)

    # add the genome to the tree (branching off root)
    in_tree = False
    max_id = 0
    for node in seq_file.tree.preOrderTraversal():
        max_id = max(max_id, node)
        if seq_file.tree.getName(node) == name:
            in_tree = True
            break
    if not in_tree:
        label = max_id + 1
        seq_file.tree.nxDg.add_edge(0, label)
        seq_file.tree.setName(label, name)
        seq_file.tree.setWeight(0, label, SeqFile.branchLen)

    # add the sequence to the map
    seq_file.pathMap[name] = fasta_path

    # write the seq file back to disk
    with open(seqfile_path, 'w') as seqfile_handle:
        seqfile_handle.write(str(seq_file))



if __name__ == "__main__":
    main()
