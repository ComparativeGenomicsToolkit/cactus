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
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from toil.lib.threading import cpu_count
from cactus.progressive.progressive_decomposition import compute_outgroups, parse_seqfile, get_subtree, get_spanning_subtree, get_event_set
from cactus.progressive.multiCactusTree import MultiCactusTree
from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory
from cactus.setup.cactus_align import cactus_align

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("seqFile", help = "Seq file (will be modified if necessary to include graph Fasta sequence)")
    parser.add_argument("minigraphGFA", help = "Minigraph-compatible reference graph in GFA format (can be gzipped)")
    parser.add_argument("outputPAF", type=str, help = "Output pairwise alignment file in PAF format")
    parser.add_argument("--outputFasta", type=str, help = "Output graph sequence file in FASTA format (required if not present in seqFile)")
    parser.add_argument("--maskFilter", type=int, help = "Ignore softmasked sequence intervals > Nbp (overrides config option of same name)")
    parser.add_argument("--delFilter", type=int, help = "Filter out split-mapping-implied deletions > Nbp")
    parser.add_argument("--outputGAFDir", type=str, help = "Output GAF alignments (raw minigraph output before PAF conversion) to this directory")
    parser.add_argument("--reference", type=str, help = "Reference genome name.  MAPQ filter will not be applied to it")
    parser.add_argument("--refFromGFA", action="store_true", help = "Do not align reference (--reference) from seqfile, and instead extract its alignment from the rGFA tags (must have been used as reference for minigraph GFA construction)")
    parser.add_argument("--base", action="store_true", help = "Use cactus to fill in (pairwise) base alignments between minigraph minimizers")
    parser.add_argument("--mapCores", type=int, help = "Number of cores for each minigraph (and base-alignment job).  Overrides graphmap cpu in configuration")

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
    set_logging_from_options(options)
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
    graph_map(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-graphmap has finished after {} seconds".format(run_time))

def graph_map(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        config_node = ET.parse(options.configFile).getroot()
        # get the minigraph "virutal" assembly name
        graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")
        #Run the workflow
        if options.restart:
            paf_id, gfa_fa_id, gaf_id_map = toil.restart()
        else:
            # load up the seqfile and figure out the outgroups and schedule
            config_wrapper = ConfigWrapper(config_node)
            config_wrapper.substituteAllPredefinedConstantsWithLiterals()
            mc_tree, input_seq_map, og_candidates = parse_seqfile(options.seqFile, config_wrapper)
            og_map = compute_outgroups(mc_tree, config_wrapper, set(og_candidates))
            event_set = get_event_set(mc_tree, config_wrapper, og_map, mc_tree.getRootName())

            # apply path overrides.  this was necessary for wdl which doesn't take kindly to
            # text files of local paths (ie seqfile).  one way to fix would be to add support
            # for s3 paths and force wdl to use it.  a better way would be a more fundamental
            # interface shift away from files of paths throughout all of cactus
            if options.pathOverrides:
                for name, override in zip(options.pathOverrideNames, options.pathOverrides):
                    input_seq_map[name] = override

            #apply the maskfilter override
            if options.maskFilter is not None:
                findRequiredNode(config_node, "graphmap").attrib["maskFilter"] = str(options.maskFilter)
            if options.delFilter is not None:
                findRequiredNode(config_node, "graphmap").attrib["delFilter"] = str(options.delFilter)

            # apply cpu override
            if options.mapCores is not None:
                findRequiredNode(config_node, "graphmap").attrib["cpu"] = str(options.mapCores)

            if graph_event in event_set:
                # dont need to import this
                event_set.remove(graph_event)

            # check --reference input
            if options.reference:
                leaves = [mc_tree.getName(leaf) for leaf in mc_tree.getLeaves()]
                if options.reference not in leaves:
                    raise RuntimeError("Genome specified with --reference, {}, not found in tree leaves".format(options.reference))
            
            if options.refFromGFA:
                if not options.reference:
                    raise RuntimeError("--reference must be used with --refFromGFA")
                # ugly, but this option used to be a string
                # todo: probably best to eventually get rid of this option entirely.
                options.refFromGFA = options.reference
                # we're not going to need the fasta for anything, so forget about it now
                del input_seq_map[options.refFromGFA]
                
            if not options.outputFasta and graph_event not in input_seq_map:
                raise RuntimeError("{} assembly not found in seqfile so it must be specified with --outputFasta".format(graph_event))

            #import the graph
            gfa_id = toil.importFile(makeURL(options.minigraphGFA))

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
            
            # run the workflow
            paf_id, gfa_fa_id, gaf_id_map, unfiltered_paf_id, paf_filter_log = toil.start(Job.wrapJobFn(
                minigraph_workflow, options, config_wrapper, input_seq_map, input_seq_id_map, gfa_id, graph_event))

        #export the paf
        toil.exportFile(paf_id, makeURL(options.outputPAF))
        if unfiltered_paf_id:
            toil.exportFile(unfiltered_paf_id, makeURL(options.outputPAF + ".unfiltered"))
            toil.exportFile(paf_filter_log, makeURL(options.outputPAF + ".filter.log"))
                        
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

def minigraph_workflow(job, options, config, seq_path_map, seq_id_map, gfa_id, graph_event):
    """ Overall workflow takes command line options and returns (paf-id, (optional) fa-id) """
    fa_id = None

    root_job = Job()
    job.addChild(root_job)
    
    if options.outputFasta or options.base:
        # convert GFA to fasta
        fa_job = root_job.addChildJobFn(make_minigraph_fasta, gfa_id, options.outputFasta, graph_event)
        fa_id = fa_job.rv()

    paf_job = Job.wrapJobFn(minigraph_map_all, config, gfa_id, seq_path_map, seq_id_map, graph_event, options.outputGAFDir is not None,
                            options.base, fa_id if options.base else None)
    if options.base:
        root_job.addFollowOn(paf_job)
    else:
        root_job.addChild(paf_job)

    if options.reference:
        # extract a PAF directly from the rGFAs tag for the given reference
        # if --refFromGFA is specified, we get the entire alignment from that, otherwise we just take contigs
        # that didn't get mapped by anything else
        gfa2paf_job = Job.wrapJobFn(extract_paf_from_gfa, gfa_id, options.minigraphGFA, options.reference, graph_event, paf_job.rv(0) if not options.refFromGFA else None,
                                    disk=gfa_id.size * 12)
        if options.refFromGFA:
            root_job.addChild(gfa2paf_job)
        else:
            paf_job.addFollowOn(gfa2paf_job)
        merge_paf_job = Job.wrapJobFn(merge_pafs, {"1" : paf_job.rv(0), "2" : gfa2paf_job.rv()}, disk=gfa_id.size * 12)
        paf_job.addFollowOn(merge_paf_job)
        gfa2paf_job.addFollowOn(merge_paf_job)
        out_paf_id = merge_paf_job.rv()
        prev_job = merge_paf_job
    else:
        out_paf_id = paf_job.rv(0)
        prev_job = paf_job

    # apply the optional deletion filter
    unfiltered_paf_id = None
    filtered_paf_log = None
    del_filter = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "delFilter", int, default=-1)
    if del_filter > 0:
        mg_cores = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "cpu", typeFn=int, default=1)
        mg_cores = min(mg_cores, cpu_count())
        del_filter_job = prev_job.addFollowOnJobFn(filter_paf_deletions, out_paf_id, gfa_id, options.minigraphGFA, del_filter,
                                                   disk=gfa_id.size * 12, cores=mg_cores)
        unfiltered_paf_id = out_paf_id
        out_paf_id = del_filter_job.rv(0)
        filtered_paf_log = del_filter_job.rv(1)

    return out_paf_id, fa_id if options.outputFasta else None, paf_job.rv(1), unfiltered_paf_id, filtered_paf_log
        
def make_minigraph_fasta(job, gfa_file_id, gfa_fa_file_path, name):
    """ Use gfatools to make the minigraph "assembly. if name specified, use it as unique prefix """

    # note: using the toil-vg convention of naming working files manually so that logging is more readable
    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    fa_path = os.path.join(work_dir, "minigraph_sequences.fa")
    
    job.fileStore.readGlobalFile(gfa_file_id, gfa_path)
    cmd = [["gfatools", "gfa2fa", gfa_path]]
    if name:
        cmd.append(["sed", "-e", "s/^>\(.\)/>id={}|\\1/g".format(name)])
    if gfa_fa_file_path and gfa_fa_file_path.endswith('.gz'):
        cmd.append(['gzip'])
        fa_path += '.gz'
    if len(cmd) == 1:
        cmd = cmd[0]
    cactus_call(outfile=fa_path, parameters=cmd)

    return job.fileStore.writeGlobalFile(fa_path)
    
def minigraph_map_all(job, config, gfa_id, fa_path_map, fa_id_map, graph_event, keep_gaf, base_alignment, gfa_fa_id):
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
                
    for event, fa_id in fa_id_map.items():
        fa_path = fa_path_map[event]
        minigraph_map_job = top_job.addChildJobFn(minigraph_map_one, config, event, fa_path, fa_id, gfa_id,
                                                  keep_gaf or not paf_per_genome, paf_per_genome,
                                                  # todo: estimate RAM
                                                  cores=mg_cores, disk=5*(fa_id.size + gfa_id.size))
        gaf_id_map[event] = minigraph_map_job.rv(0)
        if base_alignment:
            base_alignment_job = minigraph_map_job.addFollowOnJobFn(minigraph_base_align, config, event,
                                                                    fa_path, fa_id, gfa_id, gfa_fa_id, graph_event, minigraph_map_job.rv(1))
            paf_id_map[event] = base_alignment_job.rv()
        else:
            paf_id_map[event] = minigraph_map_job.rv(1)

    # convert to paf
    if paf_per_genome:
        paf_job = top_job.addFollowOnJobFn(merge_pafs, paf_id_map)
    else:
        paf_job = top_job.addFollowOnJobFn(merge_gafs_into_paf, config, gaf_id_map, [])

    if not keep_gaf:
        gaf_id_map = None
    else:
        gaf_id_map = paf_job.addFollowOnJobFn(compress_gafs, gaf_id_map).rv()
        
    return paf_job.rv(), gaf_id_map

def minigraph_map_one(job, config, event_name, fa_path, fa_file_id, gfa_file_id, gaf_output, paf_output):
    """ Run minigraph to map a Fasta file to a GFA graph, producing a GAF output """

    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    fa_path = os.path.join(work_dir, os.path.basename(fa_path))
    if fa_path == gfa_path or fa_path == gfa_path + ".gz":
        gfa_path += ".1"
    gaf_path = os.path.join(work_dir, "{}.gaf".format(event_name))
    
    job.fileStore.readGlobalFile(gfa_file_id, gfa_path)
    job.fileStore.readGlobalFile(fa_file_id, fa_path)

    if fa_path.endswith('.gz'):
        fa_path = fa_path[:-3]
        cactus_call(parameters = ['gzip', '-d', '-c', fa_path + '.gz'], outfile=fa_path)

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

    cmd = ["minigraph", gfa_path, fa_path, "-o", gaf_path] + opts_list

    mask_filter = getOptionalAttrib(xml_node, "maskFilter", int, default=-1)
    if mask_filter >= 0:
        cmd[2] = '-'
        cmd = [['cactus_softmask2hardmask', fa_path, '-m', str(mask_filter)], cmd]
    
    cactus_call(parameters=cmd)

    paf_id, gaf_id = None, None
    if paf_output:
        # optional gaf->paf step.  we are not piping directly out of minigraph because mzgaf2paf's overlap filter
        # (which is usually on) requires 2 passes so it won't read stdin when it's enabled
        paf_id =  merge_gafs_into_paf(job, config, None, [gaf_path])
    if gaf_output:
        gaf_id = job.fileStore.writeGlobalFile(gaf_path)

    return gaf_id, paf_id

def merge_gafs_into_paf(job, config, gaf_file_id_map, gaf_paths):
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
    # this must be consistent with checkUniqueHeaders.py
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

    cmd = [["mzgaf2paf"] + gaf_paths + mzgaf2paf_opts, ['sort', '-k' , '1,1', '-k', '3,3n']]
    cactus_call(outfile=paf_path, parameters=cmd)

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

def extract_paf_from_gfa(job, gfa_id, gfa_path, ref_event, graph_event, ignore_paf_id):
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
    # optional paf whose queries we ignore
    ignore_paf_path = os.path.join(work_dir, os.path.basename(gfa_path) + ".tofilter.paf")
    if ignore_paf_id:
        job.fileStore.readGlobalFile(ignore_paf_id, ignore_paf_path)
    # make the paf
    paf_path = job.fileStore.getLocalTempFile()
    cmd = ['rgfa2paf', gfa_path, '-T', 'id={}|'.format(graph_event), '-P', 'id={}|'.format(ref_event)]
    if ignore_paf_id:
        cmd += ['-i', ignore_paf_path]
    cactus_call(parameters=cmd, outfile=paf_path)
    return job.fileStore.writeGlobalFile(paf_path)

def minigraph_base_align(job, config, event, fa_path, fa_id, gfa_id, gfa_fa_id, graph_event, paf_id):
    """ send the paf from map_one through cactus-align then hal2paf to get a paf with full base alignment
    """
    work_dir = job.fileStore.getLocalTempDir()
    # download the gfa
    gfa_path = os.path.join(work_dir, 'mg.gfa')
    job.fileStore.readGlobalFile(gfa_id, gfa_path)

    # download the fa
    fa_path = os.path.join(work_dir, os.path.basename(fa_path))
    job.fileStore.readGlobalFile(fa_id, fa_path)

    # download the minigraph gfa fa (_MINIGRAPH_ sequences)
    gfa_fa_path = os.path.join(work_dir, 'mg.gfa.fa')
    job.fileStore.readGlobalFile(gfa_fa_id, gfa_fa_path)

    # put together input for cactus (pairwise)
    tree = NXNewick().parseString('({}:0.01){}:0.01;'.format(event, graph_event));
    mc_tree = MultiCactusTree(tree)
    mc_tree.nameUnlabeledInternalNodes(config.getDefaultInternalNodePrefix())
    mc_tree.computeSubtreeRoots()

    input_seq_map = {event : fa_path, graph_event : 'mg.gfa.fa'}
    input_seq_id_map = {event : fa_id, graph_event : gfa_fa_id}
    og_map = {}

    # use same cores as minigraph
    mg_cores = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "cpu", typeFn=int, default=1)
    mg_cores = min(mg_cores, cpu_count())

    # run cactus
    align_job = job.addChildJobFn(cactus_align, config, mc_tree, input_seq_map, input_seq_id_map, paf_id, graph_event, og_map, None, False, False, cons_cores=mg_cores)

    # run hal2paf
    hal2paf_job = align_job.addFollowOnJobFn(hal2paf, align_job.rv(0), graph_event, event, disk=10*gfa_fa_id.size)

    return hal2paf_job.rv()

def hal2paf(job, hal_id, graph_event, event, add_unique_prefix = True):
    """ run hal2paf """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, 'h2p.hal')
    job.fileStore.readGlobalFile(hal_id, hal_path)

    paf_path = os.path.join(work_dir, 'h2p.paf')

    cmd = ['hal2paf', hal_path, '--rootGenome', graph_event, '--onlySequenceNames']
    if add_unique_prefix:
        # unique prefixes will have gotten stripped by the cactus workflow.  best to add them back
        cmd = [cmd, ['awk', '{{$1="id={}|"$1; $6="id={}|"$6}}1'.format(event, graph_event),  'OFS=\t']]

    cactus_call(parameters=cmd, outfile=paf_path)

    return job.fileStore.writeGlobalFile(paf_path)

def filter_paf_deletions(job, paf_id, gfa_id, gfa_path, max_deletion):
    """ run filter-paf-deletions on a paf to break out giant-snarl-making edges """
    work_dir = job.fileStore.getLocalTempDir()
    paf_path = os.path.join(work_dir, 'mg.paf')
    gfa_path = os.path.join(work_dir, os.path.basename(gfa_path))
    if paf_path == gfa_path:
        paf_path += ".1"                            

    job.fileStore.readGlobalFile(paf_id, paf_path)
    job.fileStore.readGlobalFile(gfa_id, gfa_path)

    # make the vg graph
    if gfa_path.endswith('.gz'):
        cactus_call(parameters = ['gzip', '-f', '-d', gfa_path])
        gfa_path = gfa_path[:-3]

    vg_path = gfa_path + '.vg'
    trans_path = gfa_path + '.trans'

    cactus_call(parameters = ['vg', 'convert', '-r', '0', '-g', gfa_path, '-p', '-T', trans_path],
                outfile=vg_path)

    # call filter-paf-deletionts
    filter_paf_path = paf_path + ".filter"
    filter_log_path = paf_path + ".filter.log"
    cactus_call(parameters=['filter-paf-deletions', vg_path, trans_path, paf_path, str(max_deletion), '-v', '-t', str(job.cores)],
                outfile=filter_paf_path, errfile=filter_log_path)

    # return the results
    return (job.fileStore.writeGlobalFile(filter_paf_path), job.fileStore.writeGlobalFile(filter_log_path))

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
        seq_file.tree.setWeight(0, label, seq_file.branchLen)

    # add the sequence to the map
    seq_file.pathMap[name] = fasta_path

    # write the seq file back to disk
    with open(seqfile_path, 'w') as seqfile_handle:
        seqfile_handle.write(str(seq_file))

    

if __name__ == "__main__":
    main()
