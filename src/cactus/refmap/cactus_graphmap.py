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
from cactus.pipeline.cactus_workflow import prependUniqueIDs
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from toil.job import Job
from toil.common import Toil
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions
from toil.realtimeLogger import RealtimeLogger
from toil.lib.threading import cpu_count

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
    parser.add_argument("--maskFilter", type=int, help = "Ignore softmasked sequence intervals > Nbp (overrides config option of same name)")
    parser.add_argument("--outputGAFDir", type=str, help = "Output GAF alignments (raw minigraph output before PAF conversion) to this directory")
    parser.add_argument("--refFromGFA", type=str, help = "Do not align given genome from seqfile, and instead extract its alignment from the rGFA tags (must have been used as reference for minigraph GFA construction)")

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

    if options.outputGAFDir:
        if not os.path.isdir(options.outputGAFDir):
            os.makedirs(options.outputGAFDir)

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
            paf_id, gfa_fa_id, gaf_id_map = toil.restart()
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

            #apply the maskfilter override
            if options.maskFilter is not None:
                findRequiredNode(configNode, "graphmap").attrib["maskFilter"] = str(options.maskFilter)

            # get the minigraph "virutal" assembly name
            graph_event = getOptionalAttrib(findRequiredNode(configNode, "graphmap"), "assemblyName", default="_MINIGRAPH_")

            # load the seqfile
            seqFile = SeqFile(options.seqFile)

            if options.refFromGFA:
                if options.refFromGFA not in seqFile.pathMap:
                    raise RuntimeError("{}, specified with --refFromGFA, was not found in the seqfile".format(options.refFromGFA))
                # we're not going to need the fasta for anything, so forget about it now
                del seqFile.pathMap[options.refFromGFA]
                
            if not options.outputFasta and graph_event not in seqFile.pathMap:
                raise RuntimeError("{} assembly not found in seqfile so it must be specified with --outputFasta".format(graph_event))

            #import the graph
            gfa_id = toil.importFile(makeURL(options.minigraphGFA))

            #import the sequences (that we need to align for the given event, ie leaves and outgroups)
            seqIDMap = {}
            leaves = set([seqFile.tree.getName(node) for node in seqFile.tree.getLeaves()])
            for genome, seq in seqFile.pathMap.items():
                if genome != graph_event and genome in leaves and genome != options.refFromGFA:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    logger.info("Importing {}".format(seq))
                    seqIDMap[genome] = (seq, toil.importFile(seq))

            # run the workflow
            paf_id, gfa_fa_id, gaf_id_map = toil.start(Job.wrapJobFn(minigraph_workflow, options, config, seqIDMap, gfa_id, graph_event))

        #export the paf
        toil.exportFile(paf_id, makeURL(options.outputPAF))
        if gfa_fa_id:
            toil.exportFile(gfa_fa_id, makeURL(options.outputFasta))

        #export the gafs
        if options.outputGAFDir:
            for event, gaf_id in gaf_id_map.items():
                gaf_path = os.path.join(options.outputGAFDir, '{}.gaf.gz'.format(event))
                toil.exportFile(gaf_id, makeURL(gaf_path))

        # update the input seqfile (in place!)
        if options.outputFasta:
            add_genome_to_seqfile(options.seqFile, makeURL(options.outputFasta), graph_event)

def minigraph_workflow(job, options, config, seq_id_map, gfa_id, graph_event):
    """ Overall workflow takes command line options and returns (paf-id, (optional) fa-id) """
    fa_id = None
        
    if options.outputFasta:
        # convert GFA to fasta
        fa_job = job.addChildJobFn(make_minigraph_fasta, gfa_id, graph_event)
        fa_id = fa_job.rv()

    paf_job = job.addChildJobFn(minigraph_map_all, config, gfa_id, seq_id_map, graph_event, options.outputGAFDir is not None)
    
    if options.refFromGFA:
        # extract a PAF directly from the rGFAs tag for the given reference
        gfa2paf_job = job.addChildJobFn(extract_paf_from_gfa, gfa_id, options.minigraphGFA, options.refFromGFA, graph_event,
                                        disk=gfa_id.size * 12)

        merge_paf_job = Job.wrapJobFn(merge_pafs, {"1" : paf_job.rv(0), "2" : gfa2paf_job.rv()}, disk=gfa_id.size * 12)
        paf_job.addFollowOn(merge_paf_job)
        gfa2paf_job.addFollowOn(merge_paf_job)
        out_paf_id = merge_paf_job.rv()
    else:
        out_paf_id = paf_job.rv(0)

    return out_paf_id, fa_id, paf_job.rv(1)
        
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

def minigraph_map_all(job, config, gfa_id, fa_id_map, graph_event, keep_gaf):
    """ top-level job to run the minigraph mapping in parallel, returns paf """
    
    # hang everything on this job, to self-contain workflow
    top_job = Job()
    job.addChild(top_job)

    mg_cores = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "cpu", typeFn=int, default=1)
    mg_cores = min(mg_cores, cpu_count())

    # doing the paf conversion is more efficient when done separately for each genome.  we can get away
    # with doing this if the universal filter (which needs to process everything at once) is disabled
    xml_node = findRequiredNode(config.xmlRoot, "graphmap")
    paf_per_genome = not getOptionalAttrib(xml_node, "universalMZFilter", float)
    
    # do the mapping
    gaf_id_map = {}
    paf_id_map = {}
                
    for event, fa_path_fa_id in fa_id_map.items():
        fa_path = fa_path_fa_id[0]
        fa_id = fa_path_fa_id[1]
        minigraph_map_job = top_job.addChildJobFn(minigraph_map_one, config, event, fa_path, fa_id, gfa_id,
                                                  keep_gaf or not paf_per_genome, paf_per_genome,
                                                  # todo: estimate RAM
                                                  cores=mg_cores, disk=5*(fa_id.size + gfa_id.size))
        gaf_id_map[event] = minigraph_map_job.rv(0)
        paf_id_map[event] = minigraph_map_job.rv(1)

    # convert to paf
    if paf_per_genome:
        paf_job = top_job.addFollowOnJobFn(merge_pafs, paf_id_map)
    else:
        paf_job = top_job.addFollowOnJobFn(merge_gafs_into_paf, config, gaf_id_map)

    if not keep_gaf:
        gaf_id_map = None
    else:
        gaf_id_map = paf_job.addFollowOnJobFn(compress_gafs, gaf_id_map).rv()
        
    return paf_job.rv(), gaf_id_map

def minigraph_map_one(job, config, event_name, fa_path, fa_file_id, gfa_file_id, gaf_output, paf_output):
    """ Run minigraph to map a Fasta file to a GFA graph, producing a GAF output """

    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    fa_dir = job.fileStore.getLocalTempDir()
    fa_path = os.path.join(fa_dir, os.path.basename(fa_path))
    gaf_path = os.path.join(work_dir, "{}.gaf".format(event_name))
    
    job.fileStore.readGlobalFile(gfa_file_id, gfa_path)
    job.fileStore.readGlobalFile(fa_file_id, fa_path)

    if fa_path.endswith('.gz'):
        fa_path = fa_path[:-3]
        cactus_call(parameters = ['gzip', '-d', '-c', fa_path + '.gz'], outfile=fa_path)

    # prepend the unique id before mapping so the GAF has cactus-compatible event names
    fa_path = prependUniqueIDs({event_name : fa_path}, work_dir, eventNameAsID=True)[event_name]

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

    mask_filter = getOptionalAttrib(xml_node, "maskFilter", int, default=-1)
    if mask_filter >= 0:
        cmd[2] = '-'
        cmd = [['cactus_softmask2hardmask', os.path.basename(fa_path), '-m', str(mask_filter)], cmd]

    cactus_call(work_dir=work_dir, parameters=cmd)

    paf_id, gaf_id = None, None
    if paf_output:
        # optional gaf->paf step.  we are not piping directly out of minigraph because mzgaf2paf's overlap filter
        # (which is usually on) requires 2 passes so it won't read stdin when it's enabled
        paf_id =  merge_gafs_into_paf(job, config, None, [gaf_path])
    if gaf_output:
        gaf_id = job.fileStore.writeGlobalFile(gaf_path)

    return gaf_id, paf_id

def merge_gafs_into_paf(job, config, gaf_file_id_map, gaf_paths = []):
    """ Merge GAF alignments into a single PAF, applying some filters """

    work_dir = job.fileStore.getLocalTempDir()
    paf_path = os.path.join(work_dir, "mz_alignments.paf")
    if not gaf_paths:
        for event, gaf_id in gaf_file_id_map.items():
            gaf_paths.append("{}.gaf".format(event))
            job.fileStore.readGlobalFile(gaf_id, os.path.join(work_dir, gaf_paths[-1]))

    xml_node = findRequiredNode(config.xmlRoot, "graphmap")
    mzgaf2paf_opts = []
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    # this must be consistent with prependUniqueIDs() in cactus_workflow.py
    mzgaf2paf_opts += ['-p', 'id={}|'.format(graph_event)]
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
    overlap_filter_len = getOptionalAttrib(xml_node, "minGAFQueryOverlapFilter", int)
    if overlap_filter_len:
        mzgaf2paf_opts += ['-o', str(overlap_filter_len)]

    cactus_call(outfile=paf_path, parameters=["mzgaf2paf"] + gaf_paths + mzgaf2paf_opts)

    return job.fileStore.writeGlobalFile(paf_path)
                                               
def merge_pafs(job, paf_file_id_map):
    """ merge up some pafs """
    paf_paths = [job.fileStore.readGlobalFile(paf_id) for paf_id in paf_file_id_map.values()]
    merged_path = job.fileStore.getLocalTempFile()
    catFiles(paf_paths, merged_path)
    return job.fileStore.writeGlobalFile(merged_path)

def compress_gafs(job, gaf_file_id_map):
    for event, file_id in gaf_file_id_map.items():
        gaf_file_id_map[event] = job.addChildJobFn(compress_gaf, file_id, disk=int(1.5 * file_id.size)).rv()
    return gaf_file_id_map

def compress_gaf(job, gaf_file_id):
    gaf_path = job.fileStore.readGlobalFile(gaf_file_id)
    zip_path = job.fileStore.getLocalTempFile()
    cactus_call(parameters=['gzip', gaf_path, '-c', ], outfile=zip_path)
    job.fileStore.deleteGlobalFile(gaf_file_id)
    return job.fileStore.writeGlobalFile(zip_path)

def extract_paf_from_gfa(job, gfa_id, gfa_path, ref_event, graph_event):
    """ make a paf directly from the rGFA tags.  rgfa2paf supports other ranks, but we're only
    using rank=0 here to produce an alignment for the reference genome """
    work_dir = job.fileStore.getLocalTempDir()

    # download the gfa
    gfa_path = os.path.join(work_dir, os.path.basename(gfa_path))
    job.fileStore.readGlobalFile(gfa_id, gfa_path, mutable=True)
    # unzip if needed
    if gfa_path.endswith(".gz"):
        cactus_call(parameters=['gzip', '-fd', gfa_path])
        gfa_path = gfa_path[:-3]
    # make the paf
    paf_path = job.fileStore.getLocalTempFile()
    cactus_call(parameters=['rgfa2paf', gfa_path, '-T', 'id={}|'.format(graph_event), '-P', 'id={}|'.format(ref_event)], outfile=paf_path)
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
