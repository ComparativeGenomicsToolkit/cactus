#!/usr/bin/env python3

"""jobs and command for running the minigraph pan-genome building pipeline, which takes
   input: usual list of fasta files + a minigraph-compatible GFA graph
   output: PAF file containing alignment of each input Fasta to the contig sequences of the graph
           (these contig sequences are also output in their own Fasta, which can be treated as an "assembly"
            by the rest of cactus)

"""
import os, sys
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
from cactus.shared.common import unzip_gz, zip_gz
from cactus.shared.version import cactus_commit
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count
from cactus.shared.common import cactus_clamp_memory
from cactus.progressive.progressive_decomposition import compute_outgroups, parse_seqfile, get_subtree, get_spanning_subtree, get_event_set
from cactus.refmap.cactus_minigraph import check_sample_names
from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("seqFile", help = "Seq file (will be modified if necessary to include graph Fasta sequence)")
    parser.add_argument("minigraphGFA", help = "Minigraph-compatible reference graph in GFA format (can be gzipped)")
    parser.add_argument("outputPAF", type=str, help = "Output pairwise alignment file in PAF format")
    parser.add_argument("--outputFasta", type=str, help = "Output graph sequence file in FASTA format (required if not present in seqFile)")
    parser.add_argument("--maskFilter", type=int, help = "Ignore softmasked sequence intervals > Nbp (overrides config option of same name)")
    parser.add_argument("--delFilter", type=int, help = "Filter out split-mapping-implied deletions > Nbp (default will be \"delFilter\" from the config")
    parser.add_argument("--outputGAFDir", type=str, help = "Output GAF alignments (raw minigraph output before PAF conversion) to this directory")
    parser.add_argument("--reference", nargs='+', type=str, help = "Reference genome name.  MAPQ filter will not be applied to it")
    parser.add_argument("--refFromGFA", action="store_true", help = "Do not align reference (--reference) from seqfile, and instead extract its alignment from the rGFA tags (must have been used as reference for minigraph GFA construction)")
    parser.add_argument("--mapCores", type=int, help = "Number of cores for minigraph.  Overrides graphmap cpu in configuration")    

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

    # support but ignore multi reference
    if options.reference:
        options.reference = options.reference[0]

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()
    graph_map(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-graphmap has finished after {} seconds".format(run_time))

def graph_map(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        config_node = ET.parse(options.configFile).getroot()
        config_wrapper = ConfigWrapper(config_node)        
        graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")        
        if options.restart:
            paf_id, gfa_fa_id, gaf_id, unfiltered_paf_id, paf_filter_log = toil.restart()            
        else:
            # load up the seqfile and figure out the outgroups and schedule
            config_wrapper.substituteAllPredefinedConstantsWithLiterals(options)
            mc_tree, input_seq_map, og_candidates = parse_seqfile(options.seqFile, config_wrapper, pangenome=True)
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
            mg_cores = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "cpu", typeFn=int, default=1)
            if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
                mg_cores = min(mg_cores, cactus_cpu_count(), int(options.maxCores) if options.maxCores else sys.maxsize)
                findRequiredNode(config_node, "graphmap").attrib["cpu"] = str(mg_cores)
                
            # get the minigraph "virutal" assembly name
            graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")
            if graph_event in event_set:
                # dont need to import this
                event_set.remove(graph_event)

            # validate the sample names
            check_sample_names(input_seq_map.keys(), options.reference)
                
            # check --reference input (a bit redundant to above, but does additional leaf check)
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
            logger.info("Importing {}".format(options.minigraphGFA))
            gfa_id = toil.importFile(makeURL(options.minigraphGFA))

            #import the sequences (that we need to align for the given event, ie leaves and outgroups)
            seq_id_map = {}
            fa_id_map = {}
            for (genome, seq) in input_seq_map.items():
                if genome in event_set:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    logger.info("Importing {}".format(seq))
                    seq_id_map[genome] = toil.importFile(seq)
                    fa_id_map[genome] = seq

            # run the workflow
            paf_id, gfa_fa_id, gaf_id, unfiltered_paf_id, paf_filter_log, paf_was_filtered = toil.start(Job.wrapJobFn(
                minigraph_workflow, options, config_wrapper, seq_id_map, gfa_id, graph_event))

        #export the paf
        toil.exportFile(paf_id, makeURL(options.outputPAF))
        output_gaf = options.outputPAF[:-4] if options.outputPAF.endswith('.paf') else options.outputPAF
        output_gaf += '.gaf.gz'
        toil.exportFile(gaf_id, makeURL(output_gaf))
        if paf_was_filtered:
            toil.exportFile(unfiltered_paf_id, makeURL(options.outputPAF + ".unfiltered.gz"))
            toil.exportFile(paf_filter_log, makeURL(options.outputPAF + ".filter.log"))
                        
        if gfa_fa_id:
            toil.exportFile(gfa_fa_id, makeURL(options.outputFasta))

        # update the input seqfile (in place!)
        if options.outputFasta:
            add_genome_to_seqfile(options.seqFile, makeURL(options.outputFasta), graph_event)

def minigraph_workflow(job, options, config, seq_id_map, gfa_id, graph_event, sanitize=True):
    """ Overall workflow takes command line options and returns (paf-id, (optional) fa-id) """
    fa_id = None
    gfa_id_size = gfa_id.size

    # can be a list coming in from cactus-pangenome, but we only need first item
    if type(options.reference) is list:
        options.reference = options.reference[0]

    root_job = Job()
    job.addChild(root_job)

    mg_cores = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "cpu", typeFn=int, default=1)

    # enforce unique prefixes and unzip fastas
    if sanitize:
        sanitize_job = root_job.addChildJobFn(sanitize_fasta_headers, seq_id_map, pangenome=True)
        seq_id_map = sanitize_job.rv()

    zipped_gfa = options.minigraphGFA.endswith('.gz')
    if options.outputFasta:
        # convert GFA to fasta
        scale = 5 if zipped_gfa else 1
        fa_job = root_job.addChildJobFn(make_minigraph_fasta, gfa_id, options.outputFasta, graph_event,
                                        disk=scale*2*gfa_id.size, memory=cactus_clamp_memory(2*scale*gfa_id.size))
        fa_id = fa_job.rv()

    if zipped_gfa:
        # gaf2paf needs unzipped gfa, so we take care of that upfront
        gfa_unzip_job = root_job.addChildJobFn(unzip_gz, options.minigraphGFA, gfa_id, disk=5*gfa_id.size)
        gfa_id = gfa_unzip_job.rv()
        gfa_id_size *= 10
        options.minigraphGFA = options.minigraphGFA[:-3]
    paf_job = Job.wrapJobFn(minigraph_map_all, config, gfa_id, seq_id_map, graph_event)
    root_job.addFollowOn(paf_job)

    if options.reference:
        # extract a PAF directly from the rGFAs tag for the given reference
        # if --refFromGFA is specified, we get the entire alignment from that, otherwise we just take contigs
        # that didn't get mapped by anything else
        gfa2paf_job = Job.wrapJobFn(extract_paf_from_gfa, gfa_id, options.minigraphGFA, options.reference, graph_event, paf_job.rv(0) if not options.refFromGFA else None,
                                    disk=gfa_id_size, memory=cactus_clamp_memory(gfa_id_size))
        if options.refFromGFA:
            root_job.addChild(gfa2paf_job)
        else:
            paf_job.addFollowOn(gfa2paf_job)
        merge_paf_job = Job.wrapJobFn(merge_pafs, {"1" : paf_job.rv(0), "2" : gfa2paf_job.rv()}, disk=gfa_id_size)
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
        del_filter_threshold = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "delFilterThreshold", float, default=None)
        del_size_threshold = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "delFilterQuerySizeThreshold", float, default=None)
        del_filter_job = prev_job.addFollowOnJobFn(filter_paf_deletions, out_paf_id, gfa_id, del_filter, del_filter_threshold,
                                                   del_size_threshold,
                                                   disk=8*gfa_id_size, cores=mg_cores,
                                                   memory=cactus_clamp_memory(16*gfa_id_size))
        unfiltered_paf_id = prev_job.addFollowOnJobFn(zip_gz, 'mg.paf.unfiltered', out_paf_id, disk=gfa_id_size).rv()
        out_paf_id = del_filter_job.rv(0)
        filtered_paf_log = del_filter_job.rv(1)
        paf_was_filtered = del_filter_job.rv(2)

    return out_paf_id, fa_id if options.outputFasta else None, paf_job.rv(1), unfiltered_paf_id, filtered_paf_log, paf_was_filtered
                    
def make_minigraph_fasta(job, gfa_file_id, gfa_file_path, name):
    """ Use gfatools to make the minigraph "assembly" """

    # note: using the toil-vg convention of naming working files manually so that logging is more readable
    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    fa_path = os.path.join(work_dir, "minigraph_sequences.fa")
    
    job.fileStore.readGlobalFile(gfa_file_id, gfa_path)
    cmd = [["gfatools", "gfa2fa", gfa_path]]
    if name:
        cmd.append(["sed", "-e", "s/^>\(.\)/>id={}|\\1/g".format(name)])
    if gfa_file_path.endswith('.gz'):
        cmd.append(['bgzip', '--threads', str(job.cores)])
        fa_path += '.gz'
    if len(cmd) == 1:
        cmd = cmd[0]
    cactus_call(outfile=fa_path, parameters=cmd)

    return job.fileStore.writeGlobalFile(fa_path)

def minigraph_map_all(job, config, gfa_id, fa_id_map, graph_event):
    """ top-level job to run the minigraph mapping in parallel, returns paf """
    
    # hang everything on this job, to self-contain workflow
    top_job = Job()
    job.addChild(top_job)

    mg_cores = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "cpu", typeFn=int, default=1)

    # doing the paf conversion is more efficient when done separately for each genome.  we can get away
    # with doing this if the universal filter (which needs to process everything at once) is disabled
    xml_node = findRequiredNode(config.xmlRoot, "graphmap")
    paf_per_genome = not getOptionalAttrib(xml_node, "universalMZFilter", float)
    
    # do the mapping
    gaf_id_map = {}
    paf_id_map = {}
                
    for event, fa_id in fa_id_map.items():
        minigraph_map_job = top_job.addChildJobFn(minigraph_map_one, config, event, fa_id, gfa_id,
                                                  # todo: estimate RAM
                                                  cores=mg_cores, disk=5*fa_id.size + gfa_id.size,
                                                  memory=cactus_clamp_memory(72*fa_id.size + 2*gfa_id.size))
        gaf_id_map[event] = minigraph_map_job.rv(0)
        paf_id_map[event] = minigraph_map_job.rv(1)

    # merge up
    paf_merge_job = top_job.addFollowOnJobFn(merge_pafs, paf_id_map)
    gaf_merge_job = top_job.addFollowOnJobFn(merge_pafs, gaf_id_map, gzip=True)
    
    return paf_merge_job.rv(), gaf_merge_job.rv()

def minigraph_map_one(job, config, event_name, fa_file_id, gfa_file_id):
    """ Run minigraph to map a Fasta file to a GFA graph, producing a GAF output """

    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    fa_path = os.path.join(work_dir, "{}.fa".format(event_name))
    if fa_path == gfa_path or fa_path == gfa_path + ".gz":
        gfa_path += ".1"
    gaf_path = os.path.join(work_dir, "{}.gaf".format(event_name))
    
    job.fileStore.readGlobalFile(gfa_file_id, gfa_path)
    job.fileStore.readGlobalFile(fa_file_id, fa_path)

    # parse options from the config
    xml_node = findRequiredNode(config.xmlRoot, "graphmap")
    minigraph_opts = getOptionalAttrib(xml_node, "minigraphMapOptions", str, default="")     
    opts_list = minigraph_opts.split()
    # add required options if not present
    if "-c" not in opts_list:
        opts_list += ["-c"]
    if "-t" not in opts_list:
        opts_list += ["-t", str(int(job.cores))]

    cmd = []

    # optional hardmasking of softmasked fasta input (to ignore masked sequence)
    mask_filter = getOptionalAttrib(xml_node, "maskFilter", int, default=-1)
    if mask_filter >= 0:
        cmd += [['cactus_softmask2hardmask', fa_path, '-m', str(mask_filter)]]
        fa_path = '-'

    # run minigraph mapping
    cmd += [["minigraph", gfa_path, fa_path, "-o", gaf_path] + opts_list]

    cactus_call(parameters=cmd, job_memory=job.memory)

    # convert the gaf into unstable gaf (targets are node sequences)
    # note: the gfa needs to be uncompressed for this tool to work
    mg_lengths_path = gfa_path + '.node_lengths.tsv'
    unstable_gaf_path = gaf_path + '.unstable'
    cmd = ['gaf2unstable', gaf_path, '-g', gfa_path, '-o', mg_lengths_path]

    # optional gaf overlap filter
    overlap_ratio = getOptionalAttrib(xml_node, "GAFOverlapFilterRatio", typeFn=float, default=0)
    length_ratio = getOptionalAttrib(xml_node, "GAFOverlapFilterMinLengthRatio", typeFn=float, default=0)
    overlap_filter_len = getOptionalAttrib(xml_node, "minGAFQueryOverlapFilter", int, default=0)    
    min_block = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minGAFBlockLength", typeFn=int, default=0)
    min_mapq = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minMAPQ", typeFn=int, default=0)
    min_ident = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minIdentity", typeFn=float, default=0)    
    if overlap_ratio or overlap_filter_len:
        cmd = [cmd, ['gaffilter', '-', '-r', str(overlap_ratio), '-m', str(length_ratio), '-q', str(min_mapq),
                     '-b', str(min_block), '-o', str(overlap_filter_len), '-i', str(min_ident)]]
    cactus_call(parameters=cmd, outfile=unstable_gaf_path, job_memory=job.memory)

    # convert the unstable gaf into unstable paf, which is what cactus expects
    # also tack on the unique id to the target column
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    unstable_paf_path = unstable_gaf_path + '.paf'
    unstable_paf_cmd = [['gaf2paf', unstable_gaf_path, '-l', mg_lengths_path],
                        ['awk', 'BEGIN{{OFS=\"	\"}} {{$6="id={}|"$6; print}}'.format(graph_event)]]
    cactus_call(parameters=unstable_paf_cmd, outfile=unstable_paf_path, job_memory=job.memory)

    # return the stable gaf (minigraph output) and the unstable paf
    return job.fileStore.writeGlobalFile(gaf_path), job.fileStore.writeGlobalFile(unstable_paf_path)

def merge_pafs(job, paf_file_id_map, gzip=False):
    """ merge up some pafs """
    paf_paths = [job.fileStore.readGlobalFile(paf_id) for paf_id in paf_file_id_map.values()]
    merged_path = job.fileStore.getLocalTempFile()
    catFiles(paf_paths, merged_path)
    if gzip:
        cactus_call(parameters=['bgzip', merged_path, '--threads', str(job.cores)])
        merged_path += '.gz'                    
    return job.fileStore.writeGlobalFile(merged_path)

def extract_paf_from_gfa(job, gfa_id, gfa_path, ref_event, graph_event, ignore_paf_id):
    """ make a paf directly from the rGFA tags.  rgfa2paf supports other ranks, but we're only
    using rank=0 here to produce an alignment for the reference genome """
    work_dir = job.fileStore.getLocalTempDir()
    # download the gfa
    gfa_path = os.path.join(work_dir, os.path.basename(gfa_path))
    job.fileStore.readGlobalFile(gfa_id, gfa_path, mutable=True)
    # unzip if needed
    if gfa_path.endswith(".gz"):
        cactus_call(parameters=['bgzip', '-fd', gfa_path, '--threads', str(job.cores)])
        gfa_path = gfa_path[:-3]
    # optional paf whose queries we ignore
    ignore_paf_path = os.path.join(work_dir, os.path.basename(gfa_path) + ".tofilter.paf")
    if ignore_paf_id:
        job.fileStore.readGlobalFile(ignore_paf_id, ignore_paf_path)
    # make the paf
    paf_path = job.fileStore.getLocalTempFile()
    cmd = ['rgfa2paf', gfa_path, '-T', 'id={}|'.format(graph_event)]
    if ref_event:
        cmd += ['-P', 'id={}|'.format(ref_event)]
    if ignore_paf_id:
        cmd += ['-i', ignore_paf_path]
    cactus_call(parameters=cmd, outfile=paf_path)
    return job.fileStore.writeGlobalFile(paf_path)

def filter_paf(job, paf_id, config, reference=None):
    """ run basic paf-filtering.  these are quick filters that are best to do on-the-fly when reading the paf and 
        as such, they are called by cactus-graphmap-split and cactus-align, not here """
    work_dir = job.fileStore.getLocalTempDir()
    paf_path = os.path.join(work_dir, 'mg.paf')
    filter_paf_path = os.path.join(work_dir, 'mg.paf.filter')
    job.fileStore.readGlobalFile(paf_id, paf_path)

    min_block = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minGAFBlockLength", typeFn=int, default=0)
    min_mapq = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minMAPQ", typeFn=int, default=0)
    min_ident = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minIdentity", typeFn=float, default=0)
    RealtimeLogger.info("Running PAF filter with minBlock={} minMAPQ={} minIdentity={}".format(min_block, min_mapq, min_ident))
    with open(paf_path, 'r') as paf_file, open(filter_paf_path, 'w') as filter_paf_file:
        for line in paf_file:
            toks = line.split('\t')
            is_ref = reference and toks[0].startswith('id={}|'.format(reference))
            mapq = int(toks[11])
            query_len = int(toks[1])
            ident = float(toks[9]) / (float(toks[10]) + 0.00000001)
            bl = None
            for tok in toks[12:]:
                # this is a special tag that was written by gaf2paf in order to preserve the original gaf block length
                # we use it to be able to filter by the gaf block even after it's been broken in the paf
                if tok.startswith('gl:i:'):
                    bl = int(tok[5:])
                # we can also get the identity of the parent gaf block 
                if tok.startswith('gi:i:'):
                    ident = min(ident, float(toks[5:]))
            if is_ref or (mapq >= min_mapq and (bl is None or query_len <= min_block or bl >= min_block) and ident >= min_ident):
                filter_paf_file.write(line)

    overlap_ratio = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "PAFOverlapFilterRatio", typeFn=float, default=0)
    length_ratio = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "PAFOverlapFilterMinLengthRatio", typeFn=float, default=0)    
    if overlap_ratio:
        overlap_filter_paf_path = filter_paf_path + ".overlap"
        cactus_call(parameters=['gaffilter', filter_paf_path, '-p', '-r', str(overlap_ratio), '-m', str(length_ratio),
                                '-b', str(min_block), '-q', str(min_mapq), '-i', str(min_ident)],
                    outfile=overlap_filter_paf_path, job_memory=job.memory)
        filter_paf_path = overlap_filter_paf_path

    return job.fileStore.writeGlobalFile(filter_paf_path)    

def filter_paf_deletions(job, paf_id, gfa_id, max_deletion, filter_threshold, filter_query_size_threshold):
    """ run filter-paf-deletions on a paf to break out giant-snarl-making edges """
    work_dir = job.fileStore.getLocalTempDir()
    paf_path = os.path.join(work_dir, 'mg.paf')
    gfa_path = os.path.join(work_dir, 'mg.gfa')
    job.fileStore.readGlobalFile(paf_id, paf_path)
    job.fileStore.readGlobalFile(gfa_id, gfa_path)

    # make the vg graph
    vg_path = gfa_path + '.vg'
    trans_path = gfa_path + '.trans'

    cactus_call(parameters = ['vg', 'convert', '-r', '0', '-g', gfa_path, '-p', '-T', trans_path],
                outfile=vg_path, job_memory=job.memory)

    # call filter-paf-deletionts
    filter_paf_path = paf_path + ".filter"
    filter_log_path = paf_path + ".filter.log"
    filter_paf_cmd = ['filter-paf-deletions', vg_path, trans_path, paf_path, '-d', str(max_deletion), '-v', '-p', '-t', str(job.cores)]
    if filter_threshold:
        filter_paf_cmd += ['-m', str(filter_threshold)]
    if filter_query_size_threshold:
        filter_paf_cmd += ['-s', str(filter_query_size_threshold)]
    filter_stdout, filter_stderr = cactus_call(parameters=filter_paf_cmd, check_output=True, returnStdErr=True, job_memory=job.memory)
    with open(filter_log_path, 'w') as filter_log_file:
        for line in filter_stderr:
            filter_log_file.write(line)            
    with open(filter_paf_path, 'w') as filter_paf_file:
        for line in filter_stdout:
            filter_paf_file.write(line)
    unfiltered_paf_lines = int(cactus_call(parameters=['wc', '-l', paf_path], check_output=True).strip().split()[0])            
    filtered_paf_lines = int(cactus_call(parameters=['wc', '-l', filter_paf_path], check_output=True).strip().split()[0])
    assert filtered_paf_lines <= unfiltered_paf_lines
    was_filtered = filtered_paf_lines < unfiltered_paf_lines

    # return the results
    return (job.fileStore.writeGlobalFile(filter_paf_path), job.fileStore.writeGlobalFile(filter_log_path), was_filtered)  

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
