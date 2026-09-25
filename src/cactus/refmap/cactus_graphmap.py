#!/usr/bin/env python3

"""jobs and command for running the minigraph pan-genome building pipeline, which takes
   input: usual list of fasta files + a minigraph-compatible GFA graph
   output: PAF file containing alignment of each input Fasta to the contig sequences of the graph
           (these contig sequences are also output in their own Fasta, which can be treated as an "assembly"
            by the rest of cactus)

"""
import os, sys, re
import gzip
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit
import shutil

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import setupBinaries, importSingularityImage, cactus_walltime
from cactus.shared.common import GZIP_COMPRESS_BYTES_PER_SEC
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options, add_cactus_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import unzip_gz, zip_gz
from cactus.shared.version import cactus_commit
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers
from cactus.refmap.pangenome_exclusions import event_to_pansn_prefix
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count
from cactus.shared.common import cactus_clamp_memory
from cactus.progressive.progressive_decomposition import compute_outgroups, parse_seqfile, get_subtree, get_spanning_subtree, get_event_set
from cactus.refmap.cactus_minigraph import check_sample_names, minigraph_gfa_from_pansn, read_chromfile
from cactus.refmap.cactus_minigraph import GFA_RENAME_SECS_PER_GB, RAW_BYTES_PER_GZ_BYTE
from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory, getTempFile

def main():
    parser = Job.Runner.getDefaultArgumentParser()
    add_cactus_toil_options(parser)

    parser.add_argument("seqFile", help = "Seq file (will be modified if necessary to include graph Fasta sequence) (or chromfile with --batch)")
    parser.add_argument("minigraphGFA", nargs='?', default='', type=str,
                        help = "Minigraph-compatible reference graph in GFA format (can be gzipped) (don't specify when using --batch)")
    parser.add_argument("outputPAF", type=str, help = "Output pairwise alignment file in PAF format (or directory in --batch mode)")
    parser.add_argument("--outputFasta", type=str, help = "Output graph sequence file in FASTA format (required if not present in seqFile)")
    parser.add_argument("--maskFilter", type=int, help = "Ignore softmasked sequence intervals > Nbp (overrides config option of same name)")
    parser.add_argument("--delFilter", type=int, help = "Filter out split-mapping-implied deletions > Nbp (default will be \"delFilter\" from the config")
    parser.add_argument("--minIdentity", type=float, help = "Ignore PAF lines with identity (column 10/11) < this (overrides minIdentity in <graphmap> in config)")
    parser.add_argument("--reference", nargs='+', type=str, help = "Reference genome name.  MAPQ filter will not be applied to it")
    parser.add_argument("--refFromGFA", action="store_true", help = "Do not align reference (--reference) from seqfile, and instead extract its alignment from the rGFA tags (must have been used as reference for minigraph GFA construction)")
    parser.add_argument("--mapCores", type=int, help = "Number of cores for minigraph.  Overrides graphmap cpu in configuration")
    parser.add_argument("--collapse", help = "Incorporate minimap2 self-alignments.", action='store_true', default=False)
    parser.add_argument("--collapseRefPAF", help ="Incorporate given (reference-only) self-alignments in PAF format [Experimental]")
    parser.add_argument("--inGAF", type=str, default=None,
                        help = "Reuse the mappings in this GAF (as published by a previous cactus-graphmap or cactus-pangenome run) "
                        "instead of re-running minigraph for the genomes it covers. Minigraph GAF is in stable coordinates, which "
                        "node splitting does not change, so these mappings are re-derived against the given graph instead. Only "
                        "genomes the GAF does not cover are mapped. Intended for use with cactus-minigraph --inGFA")
    parser.add_argument("--remap", action="store_true", default=False,
                        help = "Map every genome with minigraph even if --inGAF already covers it. Slower, but the existing "
                        "genomes then see the nodes contributed by the newly added ones, as they would in a from-scratch run")

    parser.add_argument("--batch", action="store_true",
                        help="Run independently on set of chromosomea inputs (chromfile as from cactus-minigraph --batch). Note that the output will be a directory and not a PAF")
    parser.add_argument("--mgSplit", action="store_true", default=False,
                        help="Relax block-length, overlap and deletion filters because this PAF will only feed cactus-graphmap-split for chromosome binning (matches the cactus-pangenome --mgSplit / cactus-minigraph --refOnly pipeline). Do not use on per-chromosome (--batch) runs.")

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

    # support but ignore multi reference
    if options.reference:
        options.reference = options.reference[0]    

    if options.collapseRefPAF:
        if not options.collapseRefPAF.endswith('.paf'):
            raise RuntimeError('file passed to --collapseRefPAF must end with .paf')
        if not options.reference:
            raise RuntimeError('--reference must be used with --collapseRefPAF')
        if options.collapse:
            raise RuntimeError('--collapseRefPAF cannot be used with --collapse')
    
    if options.batch and options.minigraphGFA:
        raise RuntimeError("minigraphGFA argument must *not* be specified when using --batch")
    if not options.batch and not options.minigraphGFA:
        raise RuntimeError("minigraphGFA argument must be specified when *not* using --batch")

    if options.batch and options.outputFasta:
        raise RuntimeError("--outputFasta cannot be used with --batch")

    if options.mgSplit and options.batch:
        raise RuntimeError("--mgSplit is for the whole-genome splitting pass and cannot be used with --batch")

    if options.inGAF:
        if options.batch:
            raise RuntimeError("--inGAF cannot be used with --batch")
        if options.collapse or options.collapseRefPAF:
            raise RuntimeError("--inGAF cannot be used with --collapse or --collapseRefPAF: collapse PAFs are minimap2 "
                               "self-alignments and are not derived from the GAF")
    elif options.remap:
        raise RuntimeError("--remap only means something with --inGAF, which is what it overrides")
    
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    if options.batch:
        # the output paf is a directory, make sure it's there
        if not os.path.isdir(options.outputPAF):
            os.makedirs(options.outputPAF)

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

        # map chrom name to seqFile, gfa
        input_map = {}
        if options.batch:
            input_map = read_chromfile(options.seqFile)
        else:
            input_map['all'] = options.seqFile, options.minigraphGFA
        
        if options.restart:
            # output_dict maps chrom -> (paf_id, gfa_fa_id, gaf_id, unfiltered_paf_id, paf_filter_log, paf_was_filtered)
            output_dict = toil.restart()
        else:
            # load the config
            config_wrapper.substituteAllPredefinedConstantsWithLiterals(options)

            # relax splitting-stage filters first; explicit --delFilter / --maskFilter below override
            if options.mgSplit:
                apply_mgsplit_filter_overrides(config_node)

            #apply the maskfilter override
            if options.maskFilter is not None:
                findRequiredNode(config_node, "graphmap").attrib["maskFilter"] = str(options.maskFilter)
            if options.delFilter is not None:
                findRequiredNode(config_node, "graphmap").attrib["delFilter"] = str(options.delFilter)
            if options.minIdentity is not None:
                findRequiredNode(config_node, "graphmap").attrib["minIdentity"] = str(options.minIdentity)

            # apply cpu override                
            if options.mapCores is not None:
                findRequiredNode(config_node, "graphmap").attrib["cpu"] = str(options.mapCores)
            mg_cores = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "cpu", typeFn=int, default=1)
            if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
                mg_cores = min(mg_cores, cactus_cpu_count(), int(options.maxCores) if options.maxCores else sys.maxsize)
                findRequiredNode(config_node, "graphmap").attrib["cpu"] = str(mg_cores)

            # apply the collapse overrides
            if options.collapse:
                findRequiredNode(config_node, "graphmap").attrib["collapse"] = 'all'
            if options.collapseRefPAF:
                assert options.reference
                if options.batch:
                    raise RuntimeError('--collapseRefPAF not supported with --batch')
                findRequiredNode(config_node, "graphmap").attrib["collapse"] = 'reference'

            if '://' not in options.outputPAF:
                options.outputPAF = os.path.abspath(options.outputPAF)
                
            # get the minigraph "virutal" assembly name
            graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")
                                
            # load up the chromfile / seqfiles
            input_dict = {} # chrom -> seq_id_map, gfa_id, ref_collapse_paf_id, seqfile_path, gfa_path
            for chrom, inputs in input_map.items():
                input_seqfile, input_gfa = inputs[0], inputs[1]
            
                seqFile = SeqFile(input_seqfile, defaultBranchLen=config_wrapper.getDefaultBranchLen(pangenome=True))
                input_seq_map = seqFile.pathMap

                # validate the sample names
                check_sample_names(input_seq_map.keys(), options.reference)
                
                # check --reference input (a bit redundant to above, but does additional leaf check)
                if options.reference:
                    leaves = [seqFile.tree.getName(leaf) for leaf in seqFile.tree.getLeaves()]
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
                
                if not options.outputFasta and not options.batch and graph_event not in input_seq_map:
                    raise RuntimeError("{} assembly not found in seqfile so it must be specified with --outputFasta".format(graph_event))

                #import the graph
                gfa_id = toil.importFile(makeURL(input_gfa))

                #import the reference collapse paf
                ref_collapse_paf_id = None
                if options.collapseRefPAF:
                    ref_collapse_paf_id = toil.importFile(options.collapseRefPAF)

                #import the sequences (that we need to align for the given event, ie leaves and outgroups)
                seq_id_map = {}
                fa_id_map = {}
                for (genome, seq) in input_seq_map.items():
                    if genome != graph_event:
                        if os.path.isdir(seq):
                            tmpSeq = getTempFile()
                            catFiles([os.path.join(seq, subSeq) for subSeq in sorted(os.listdir(seq))], tmpSeq)
                            seq = tmpSeq
                        seq = makeURL(seq)
                        seq_id_map[genome] = toil.importFile(seq)
                        fa_id_map[genome] = seq

                input_dict[chrom] = seq_id_map, gfa_id, ref_collapse_paf_id, input_map[chrom][0], input_map[chrom][1]
                
            #import the mappings to reuse
            in_gaf_id = toil.importFile(makeURL(options.inGAF)) if options.inGAF and not options.remap else None

            # run the workflow
            # output_dict is chrom -> paf_id, gfa_fa_id, gaf_id, unfiltered_paf_id, paf_filter_log, paf_was_filtered
            output_dict = toil.start(Job.wrapJobFn(minigraph_batch_separate_workflow, options, config_wrapper, input_dict, graph_event, True,
                                                   in_gaf_id=in_gaf_id, walltime=cactus_walltime()))

        export_graphmap_output(options, config_node, input_map, output_dict, toil)

def export_graphmap_output(options, config_node, input_map, output_dict, toil):
    graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    if options.batch:
        chrom_file_path = os.path.join(options.outputPAF, 'chromfile.gm.txt')
        if chrom_file_path.startswith('s3://'):
            chrom_file_temp_path = getTempFile()
        else:
            chrom_file_temp_path = chrom_file_path                    
        chromfile = open(chrom_file_temp_path, 'w')
        construct_chromfile = read_chromfile(options.seqFile)
    for chrom, output_ids in output_dict.items():
        paf_id, gfa_fa_id, gaf_id, unfiltered_paf_id, paf_filter_log, paf_was_filtered = output_ids[:6]
        # the batch path appends the log of the pass that holds multi-reference-contig bins apart
        separate_log_id = output_ids[6] if len(output_ids) > 6 else None
        if options.batch:
            paf_path = os.path.join(options.outputPAF, chrom + '.paf')
        else:
            paf_path = options.outputPAF

        #export the paf / gaf
        toil.exportFile(paf_id, makeURL(paf_path))
        output_gaf = paf_path[:-4] if paf_path.endswith('.paf') else paf_path
        output_gaf += '.gaf.gz'
        toil.exportFile(gaf_id, makeURL(output_gaf))
        if paf_was_filtered:
            toil.exportFile(unfiltered_paf_id, makeURL(paf_path + ".unfiltered.gz"))
            toil.exportFile(paf_filter_log, makeURL(paf_path + ".filter.log"))
        if separate_log_id:
            toil.exportFile(separate_log_id, makeURL(os.path.join(os.path.dirname(paf_path),
                                                                 'mgSplit.{}.log'.format(chrom))))

        #export the fa and add it to the seqfile
        out_seqfile_path = getTempFile() if options.batch else input_map[chrom][0]
        if gfa_fa_id:
            if options.batch:
                fa_path = paf_path[:-4] + '.sv.gfa.fa.gz'
            else:
                fa_path = options.outputFasta                
            toil.exportFile(gfa_fa_id, makeURL(fa_path))

            # update the input seqfile (in place!)
            if options.batch:
                shutil.copyfile(input_map[chrom][0], out_seqfile_path)
            else:
                assert input_map[chrom][0] == options.seqFile
            add_genome_to_seqfile(out_seqfile_path, makeURL(fa_path), graph_event)

        #update the chromfile and copy the seqfile
        if options.batch:
            toil.exportFile(out_seqfile_path, paf_path[:-4] + '.gm.seqfile')
            chromfile.write('{}\t{}\t{}\t{}\n'.format(chrom, paf_path[:-4] + '.gm.seqfile', paf_path, construct_chromfile[chrom][2]))

    if options.batch:
        chromfile.close()
        if chrom_file_path.startswith('s3://'):
            write_s3(chrom_file_temp_path, chrom_file_path)

def minigraph_batch_workflow(job, options, config, input_dict, graph_event, sanitize, pansn_gfa_input=True, in_gaf_id=None):
    """ Batch wrapper to run grpahmap independently at the chromosome level."""
    output_dict = {}
    options.mg_chrom_name = None
    for chrom, input_info in input_dict.items():
        seq_id_map, gfa_id, ref_collapse_paf_id, seqfile, gfa = input_info
        if options.batch:
            chrom_options = copy.deepcopy(options)
            chrom_options.minigraphGFA = gfa
            chrom_options.seqFile = seqfile
            chrom_options.outputFasta = 'sv.gfa.fa.gz'
            chrom_options.mg_chrom_name = chrom
        else:
            chrom_options = options
        mgwf_job = job.addChildJobFn(minigraph_workflow, chrom_options, config, seq_id_map, gfa_id, graph_event,
                                     sanitize, ref_collapse_paf_id, pansn_gfa_input, in_gaf_id=in_gaf_id,
                                     walltime=cactus_walltime())
        output_dict[chrom] = mgwf_job.rv()
    return output_dict

def add_separate_ref_contigs_job(batch_job, options, config, input_dict):
    """ chain the pass that holds multi-reference-contig bins apart onto a minigraph_batch_workflow
    job, and return it.  it has to hang off the batch job rather than being added inside it: whoever
    reads these results must run after the separation pass's own children, and sibling follow-ons of the
    same job run concurrently, so a follow-on added inside would race with the caller's """
    # imported here because cactus_graphmap_split imports this module
    from cactus.refmap.cactus_graphmap_split import separate_ref_contigs_batch
    reference = options.reference[0] if type(options.reference) is list else options.reference
    return batch_job.addFollowOnJobFn(separate_ref_contigs_batch, config, input_dict, batch_job.rv(), reference,
                                      getattr(options, 'permissiveContigFilter', None),
                                      whole_genome_ref=getattr(options, 'mgSplitWholeGenomeRef', False),
                                      walltime=cactus_walltime())

def minigraph_batch_separate_workflow(job, options, config, input_dict, graph_event, sanitize, pansn_gfa_input=True, in_gaf_id=None):
    """ minigraph_batch_workflow followed by the separation pass, for callers that just want the final
    result and add nothing after it """
    batch_job = job.addChildJobFn(minigraph_batch_workflow, options, config, input_dict, graph_event, sanitize,
                                  pansn_gfa_input, in_gaf_id=in_gaf_id, walltime=cactus_walltime())
    return add_separate_ref_contigs_job(batch_job, options, config, input_dict).rv()

# Walltime estimates for the graphmap jobs.  Everything below was measured on the two HPRC
# pangenome runs, the biggest graphmaps we have logs for: the whole-panel pass of v2.0 (one 8.3
# GiB raw GFA, a 34.4 GiB merged PAF) and the per-chromosome pass of v2.1 (~0.6 GiB of GFA and
# ~2.2 GiB of PAF each).

# The merged minigraph PAF comes out at about 4x the raw GFA it was mapped against: 34.4 GiB
# against 8.3 GiB whole-panel, 2.2 GiB against 0.6 GiB per chromosome.  The PAF is only a promise
# while the workflow is being built, so this is how its size reaches the walltimes of the jobs
# that make it, read it and copy it.
PAF_BYTES_PER_GFA_BYTE = 4

# Seconds per GB of raw GFA for rgfa2paf: 402s on the 8.3 GiB whole-panel GFA against a 112s worst
# case on the 0.6 GiB per-chromosome ones, ie ~40 s/GB on top of a ~90s fixed cost.
RGFA2PAF_SECS_PER_GB = 40

# Seconds per GB of raw GFA for filter_paf_deletions: 4337s for filter-paf-deletions plus ~500s
# for the vg convert that precedes it on the whole-panel GFA, against 802s worst case on the
# per-chromosome ones.
#
# This looks far too high under --mgSplitWholeGenomeRef, where the whole job came to 563 s over
# 50 jobs against a graph term of 4386 s.  Do not cut the rate to close that gap: the rate is
# right and the size it is applied to is wrong.  gfa_id_size is the compressed GFA expanded by
# the hardcoded 10 above, which the whole-panel measurement supports (0.83 GiB gz to 8.3 GiB raw)
# but these runs do not -- their 27 unzip_gz jobs report 3.12 GB of raw GFA from a 757 MB input,
# a ratio of 4.12.  Cutting the rate to 120 fits the split runs and leaves the whole-panel case
# asking 1669 s for work that was measured at 4837 s, which is the wrong trade.  The expansion
# ratio is the thing to fix, once it is known why the two graphs differ by 2.4x.
FILTER_PAF_DELETIONS_SECS_PER_GB = 500


def minigraph_workflow(job, options, config, seq_id_map, gfa_id, graph_event, sanitize, ref_collapse_paf_id, pansn_gfa_input=True,
                       in_gaf_id=None):
    """ Overall workflow takes command line options and returns (paf-id, (optional) fa-id) """
    fa_id = None
    gfa_id_size = gfa_id.size
    genome_names = set(seq_id_map.keys())

    # can be a list coming in from cactus-pangenome, but we only need first item
    if type(options.reference) is list:
        options.reference = options.reference[0]

    root_job = Job(walltime=cactus_walltime())
    job.addChild(root_job)

    mg_cores = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "cpu", typeFn=int, default=1)

    # enforce unique prefixes and unzip fastas
    if sanitize:
        sanitize_job = root_job.addChildJobFn(sanitize_fasta_headers, seq_id_map, pangenome=True, walltime=cactus_walltime())
        seq_id_map = sanitize_job.rv()

    # add unique prefixes to the input PAF
    if ref_collapse_paf_id:
        # one awk pass over the PAF, then staging it in and out.  --collapseRefPAF was not used in
        # any of the runs we have logs for, so the ~50 MB/s awk rate is inferred, not measured
        ref_collapse_paf_id = root_job.addChildJobFn(add_paf_prefixes, ref_collapse_paf_id, options.reference,
                                                     disk=2*ref_collapse_paf_id.size,
                                                     walltime=cactus_walltime(20 * ref_collapse_paf_id.size / 1e9,
                                                                              io_bytes=2*ref_collapse_paf_id.size)).rv()

    # convert the GFA from PanSN to Cactus names
    if pansn_gfa_input:
        # the renaming pass decompresses the GFA before bgzipping it back up, so it needs room for
        # the raw copy (reckoned at 10x, as elsewhere) on top of the compressed input and output
        rename_gfa_job = root_job.addChildJobFn(minigraph_gfa_from_pansn, genome_names, options.minigraphGFA, gfa_id,
                                                disk=gfa_id.size*12,
                                                walltime=cactus_walltime(GFA_RENAME_SECS_PER_GB * gfa_id.size / 1e9,
                                                                         io_bytes=2*gfa_id.size))
        new_root_job = Job(walltime=cactus_walltime())
        root_job.addFollowOn(new_root_job)
        root_job = new_root_job
        gfa_id = rename_gfa_job.rv(0)

    # split up any mappings we've been given to reuse, so each genome's re-derivation is its own
    # job just as its mapping would have been
    in_gaf_map = None
    if in_gaf_id:
        # one pass over the reused GAF, splitting it per genome: I/O rather than compute
        split_gaf_job = root_job.addChildJobFn(split_gaf_by_event, in_gaf_id, genome_names, options.inGAF,
                                               disk=12*in_gaf_id.size,
                                               walltime=cactus_walltime(0, io_bytes=RAW_BYTES_PER_GZ_BYTE * in_gaf_id.size))
        in_gaf_map = split_gaf_job.rv()

    zipped_gfa = options.minigraphGFA.endswith('.gz')
    if options.outputFasta:
        # convert GFA to fasta
        scale = 5 if zipped_gfa else 1
        # gfatools barely scales with GFA size -- 199s on the whole-panel GFA against a 239s worst
        # case per chromosome -- so a flat compute term plus the staging is the honest shape
        fa_job = root_job.addChildJobFn(make_minigraph_fasta, gfa_id, options.outputFasta, graph_event,
                                        disk=scale*2*gfa_id_size, memory=cactus_clamp_memory(2*scale*gfa_id_size),
                                        walltime=cactus_walltime(300, io_bytes=2*gfa_id_size))
        fa_id = fa_job.rv()

    if zipped_gfa:
        # gaf2paf needs unzipped gfa, so we take care of that upfront
        # gunzip itself is fast (27s for the 0.83 GiB compressed whole-panel GFA); what costs is
        # writing the ~10x bigger raw GFA back to the jobstore
        gfa_unzip_job = root_job.addChildJobFn(unzip_gz, options.minigraphGFA, gfa_id, delete_original=False, disk=5*gfa_id_size,
                                               walltime=cactus_walltime(60, io_bytes=(1 + RAW_BYTES_PER_GZ_BYTE) * gfa_id_size))
        gfa_id = gfa_unzip_job.rv()
        gfa_id_size *= 10
        options.minigraphGFA = options.minigraphGFA[:-3]

    # size of the merged PAF every job below either makes, reads or copies
    paf_bytes = PAF_BYTES_PER_GFA_BYTE * gfa_id_size

    if in_gaf_id:
        # resolving a reused GAF is the same work for every genome, so a GAF that does not belong
        # to this graph fails identically in all of them -- once per genome, after the fan-out, and
        # again on every Toil retry.  Checking a sample up front turns that into one quick failure
        # with something actionable in it.  chained onto the unzip (when there is one) because it
        # needs the same uncompressed graph the per-genome jobs use
        check_parent = gfa_unzip_job if zipped_gfa else root_job
        check_parent.addFollowOnJobFn(check_reusable_gaf, config, in_gaf_id, gfa_id, genome_names,
                                      options.inGAF, options.minigraphGFA,
                                      disk=4*gfa_id_size, memory=cactus_clamp_memory(2*gfa_id_size),
                                      walltime=cactus_walltime(GAF_CHECK_SECS, io_bytes=gfa_id_size + in_gaf_id.size))

    paf_job = Job.wrapJobFn(minigraph_map_all, options, config, gfa_id, seq_id_map, graph_event, in_gaf_map,
                            walltime=cactus_walltime())
    root_job.addFollowOn(paf_job)

    collapse_paf_id = ref_collapse_paf_id
    if options.reference:
        # extract a PAF directly from the rGFAs tag for the given reference
        # if --refFromGFA is specified, we get the entire alignment from that, otherwise we just take contigs
        # that didn't get mapped by anything else
        gfa2paf_job = Job.wrapJobFn(extract_paf_from_gfa, gfa_id, options.minigraphGFA, options.reference, graph_event, paf_job.rv(0) if not options.refFromGFA else None,
                                    disk=gfa_id_size, memory=cactus_clamp_memory(gfa_id_size),
                                    walltime=cactus_walltime(120 + RGFA2PAF_SECS_PER_GB * gfa_id_size / 1e9,
                                                             io_bytes=gfa_id_size + (0 if options.refFromGFA else paf_bytes)))
        if options.refFromGFA:
            root_job.addChild(gfa2paf_job)
        else:
            paf_job.addFollowOn(gfa2paf_job)
        collapse_mode = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "collapse", typeFn=str, default="none")
        if collapse_mode in ['reference', 'all', 'nonref']:
            collapse_job = paf_job.addChildJobFn(self_align_all, config, seq_id_map, options.reference, collapse_mode, walltime=cactus_walltime())
            
            if ref_collapse_paf_id:
                collapse_paf_id = collapse_job.addFollowOnJobFn(merge_pafs, {"1":collapse_job.rv(), "2":ref_collapse_paf_id},
                                                                disk=gfa_id_size,
                                                                walltime=merge_pafs_walltime(paf_bytes)).rv()
            else:
                collapse_paf_id = collapse_job.rv()
        merge_paf_job = Job.wrapJobFn(merge_pafs,  {"1" : paf_job.rv(0), "2" : gfa2paf_job.rv()}, disk=gfa_id_size,
                                      walltime=merge_pafs_walltime(paf_bytes))
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
    paf_was_filtered = False
    del_filter = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "delFilter", int, default=-1)
    if del_filter > 0:
        del_filter_threshold = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "delFilterThreshold", float, default=None)
        del_size_threshold = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "delFilterQuerySizeThreshold", float, default=None)
        del_max_remove = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "delFilterMaxRemove", int, default=None)
        del_min_support = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "delFilterMinSupport", int, default=None)
        del_filter_job = prev_job.addFollowOnJobFn(filter_paf_deletions, out_paf_id, gfa_id, del_filter, del_filter_threshold,
                                                   del_size_threshold, del_max_remove, del_min_support,
                                                   disk=8*gfa_id_size, cores=mg_cores,
                                                   memory=cactus_clamp_memory(30*gfa_id_size),
                                                   walltime=cactus_walltime(600 + FILTER_PAF_DELETIONS_SECS_PER_GB * gfa_id_size / 1e9,
                                                                            io_bytes=gfa_id_size + 2*paf_bytes))
        # the 600s floor is for the tail: per chromosome this gzip has a 79s median but a 1390s
        # worst case, which is contention on the shared filesystem, not PAF size
        unfiltered_paf_id = prev_job.addFollowOnJobFn(zip_gz, 'mg.paf.unfiltered', out_paf_id, delete_original=False,
                                                      disk=gfa_id_size,
                                                      walltime=cactus_walltime(600 + paf_bytes / GZIP_COMPRESS_BYTES_PER_SEC,
                                                                               io_bytes=2*paf_bytes)).rv()
        out_paf_id = del_filter_job.rv(0)
        filtered_paf_log = del_filter_job.rv(1)
        paf_was_filtered = del_filter_job.rv(2)
        prev_job = del_filter_job

    if collapse_paf_id:
        # note: the collapse paf doesn't get merged into unfiltered_paf
        merge_collapse_job = prev_job.addFollowOnJobFn(merge_pafs, {"1" : out_paf_id, "2" : collapse_paf_id}, disk=gfa_id_size,
                                                       walltime=merge_pafs_walltime(paf_bytes))
        out_paf_id = merge_collapse_job.rv()

    return out_paf_id, fa_id if options.outputFasta else None, paf_job.rv(1), unfiltered_paf_id, filtered_paf_log, paf_was_filtered

def add_paf_prefixes(job, paf_id, name):
    """ Add prefixes to paf """
    work_dir = job.fileStore.getLocalTempDir()
    paf_path = os.path.join(work_dir, name + ".paf")
    job.fileStore.readGlobalFile(paf_id, paf_path)
    renamed_paf_path = os.path.join(work_dir, name + '.prefixed.paf')
    cmd = ['awk', 'BEGIN{{OFS=\"	\"}} {{$1="id={}|"$1; $6="id={}|"$6; print}}'.format(name, name), paf_path]
    cactus_call(parameters=cmd, outfile=renamed_paf_path)
    return job.fileStore.writeGlobalFile(renamed_paf_path)    
    
def make_minigraph_fasta(job, gfa_file_id, gfa_file_path, name):
    """ Use gfatools to make the minigraph "assembly" """
    # note: using the toil-vg convention of naming working files manually so that logging is more readable
    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    fa_path = os.path.join(work_dir, "minigraph_sequences.fa")
    
    job.fileStore.readGlobalFile(gfa_file_id, gfa_path)
    cmd = [["gfatools", "gfa2fa", gfa_path]]
    if name:
        cmd.append(["sed", "-e", r"s/^>\(.\)/>id={}|\1/g".format(name)])
    if gfa_file_path.endswith('.gz'):
        cmd.append(['bgzip', '--threads', str(job.cores)])
        fa_path += '.gz'
    if len(cmd) == 1:
        cmd = cmd[0]
    cactus_call(outfile=fa_path, parameters=cmd)

    return job.fileStore.writeGlobalFile(fa_path)

# minigraph mapping, per job and per GB of sanitized fasta.  Both re-fitted on the 27,097
# mappings of an HPRC v2.1 run, the first at scale with the faster minigraph, which separate
# cleanly into the two populations the line was drawn through: 456 whole-genome haplotypes
# (~3.1 GB of fasta, p50 712 s, p90 1053 s, max 1499 s) and 26,641 per-chromosome ones
# (0.05-0.26 GB, p50 100 s, p99 440 s, max 713 s).  Through those the slope is 250 s/GB at the
# p90 and 330 s/GB at the p99, so 400 sits above the measurement everywhere; 1200 came from the
# slower fork, where the same whole-genome mapping ran to a p90 of 4062 s.
#
# The intercept is the tail allowance, and it used to be sized against a per-chromosome max of
# 3521 s over 11,390 invocations -- 3x its own p99, taken to be cluster contention rather than
# anything the fasta size can see.  That tail is now 713 s against a p99 of 440 s, 1.6x rather
# than 3x, so 2000 was buying a cushion that costs more than it is worth: it put all 18,927
# mapping jobs of that run above an hour, which on a cluster whose shortest partition is an hour
# is the entire scheduling decision.  600 holds every per-chromosome mapping under the hour while
# still covering the observed worst by 3x, and leaves the whole-genome pass -- where the per-GB
# term dominates anyway -- around 2 h against a 1499 s worst.
MINIGRAPH_MAP_SECS = 600
MINIGRAPH_MAP_SECS_PER_GB = 400

# Re-deriving one genome's PAF from a GAF it already has (--inGAF): the same gaf2unstable/gaffilter/
# gaf2paf chain minigraph_map_one runs, without the minigraph.  Those three came to p99 18s and max
# 29s per genome across the HPRC runs (gaf2unstable|gaffilter n=11380), so this is the GAF read and
# the graph load rather than the chain itself, keyed off the shard the way the memory request is.
TRANSLATE_GAF_SECS_PER_GB = 400

# check_reusable_gaf loads the graph and resolves a sample of records against it.
GAF_CHECK_SECS = 600

def minigraph_map_all(job, options, config, gfa_id, fa_id_map, graph_event, in_gaf_map=None):
    """ top-level job to run the minigraph mapping in parallel, returns paf.

    a genome that in_gaf_map already has mappings for has its PAF re-derived from them rather
    than being mapped again -- see translate_gaf_one() """
    # hang everything on this job, to self-contain workflow
    top_job = Job(walltime=cactus_walltime())
    job.addChild(top_job)

    mg_cores = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "cpu", typeFn=int, default=1)

    # do the mapping
    gaf_id_map = {}
    paf_id_map = {}
                
    # the estimate below is anchored on the query, which holds while the graph is no bigger than the
    # chromosome the query came from.  the --mgSplitWholeGenomeRef second pass breaks that -- the graph
    # is whole-genome while the query stays one chromosome, so the index dominates and the graph term
    # has to carry it: measured at ~5.5x the (already decompressed) GFA on HPRC, against the 2x below.
    # it must be gated on batch as well as the option: the option's own first pass maps whole-genome
    # queries against a whole-genome graph, where the query anchor still holds and 2x is right
    gfa_coefficient = 6 if options.batch and getattr(options, 'mgSplitWholeGenomeRef', False) else 2

    # every genome whose contigs can name a step in the GAF's paths, which is more than the
    # genomes being mapped: --refFromGFA takes the reference out of the sequence map
    genome_names = set(fa_id_map.keys())
    if options.reference:
        genome_names.add(options.reference if type(options.reference) is str else options.reference[0])
    for event, fa_id in fa_id_map.items():
        mem = 72*fa_id.size + gfa_coefficient*gfa_id.size
        event_name = event
        if options.batch:
            # the memory heuristc seems to drastically underestimate some chromosomes in batch mode...
            mem *= 2
            event_name = '{}.{}'.format(event, options.mg_chrom_name)
        if in_gaf_map and event in in_gaf_map:
            # no minigraph, and no input fasta: gaf2unstable/gaffilter/gaf2paf against the new graph
            # is the whole job.  gaffilter reads its input into memory, as it does when mapping
            gaf_shard_id = in_gaf_map[event]
            map_job = top_job.addChildJobFn(translate_gaf_one, config, event_name, gaf_shard_id, gfa_id, genome_names,
                                            disk=12*gaf_shard_id.size + 2*gfa_id.size,
                                            memory=cactus_clamp_memory(24*gaf_shard_id.size + 4*gfa_id.size),
                                            walltime=cactus_walltime(TRANSLATE_GAF_SECS_PER_GB * gaf_shard_id.size / 1e9,
                                                                     io_bytes=2*gaf_shard_id.size + gfa_id.size))
        else:
            map_job = top_job.addChildJobFn(minigraph_map_one, config, event_name, fa_id, gfa_id,
                                            cores=mg_cores, disk=5*fa_id.size + gfa_id.size,
                                            memory=cactus_clamp_memory(mem),
                                            walltime=cactus_walltime(MINIGRAPH_MAP_SECS + MINIGRAPH_MAP_SECS_PER_GB * fa_id.size / 1e9,
                                                                     io_bytes=2*fa_id.size + gfa_id.size))
        gaf_id_map[event] = map_job.rv(0)
        paf_id_map[event] = map_job.rv(1)

    # merge up.  these two are the merges whose inputs scale with the number of genomes, so they get
    # sized off them rather than taking the default; the GAF one also bgzips, so give it the mapping
    # cores instead of leaving bgzip single-threaded.  merge_pafs_sized resolves the promises and
    # sets the real disk and walltime from them, so these two are just the coordination job
    merge_name = getattr(options, 'mg_chrom_name', None) if options.batch else None
    merge_name = merge_name if merge_name else 'merged'
    paf_merge_job = top_job.addFollowOnJobFn(merge_pafs_sized, paf_id_map,
                                             merged_name='{}.paf'.format(merge_name), walltime=cactus_walltime())
    gaf_merge_job = top_job.addFollowOnJobFn(merge_pafs_sized, gaf_id_map, gzip=True,
                                             merged_name='{}.gaf'.format(merge_name),
                                             gzip_cores=mg_cores, walltime=cactus_walltime())

    return paf_merge_job.rv(), gaf_merge_job.rv()

# id=EVENT|CONTIG, as it appears in a stable GAF's query column and in each of its path segments
# anchored to a field start (line start or tab) or a path-segment orientation mark, because
# "id=" is only a name prefix in those positions.  unanchored it also fires inside a contig
# name that happens to contain id=...|, rewriting a name that exists in no graph and no input
gaf_pansn_re = re.compile(r'(^|[\t><])id=([^|\t\n<>]+)\|')

def pansn_to_event_map(names):
    """ SAMPLE#HAP -> seqfile event, for the genomes in names.  event_to_pansn_prefix is lossy
    (both S288C and S288C.0 give S288C#0), so going back needs the event names to hand """
    prefix_map = {}
    for event in names:
        prefix_map.setdefault(event_to_pansn_prefix(event), event)
    return prefix_map

# SAMPLE#HAP, as it appears in a published (PanSN) GAF's query column and in each of its path
# segments.  anchored like gaf_pansn_re above, and stopping at the second '#' so that a PanSN
# phase block (SAMPLE#HAP#CONTIG#PHASEBLOCK) is left in the contig part where it belongs
pansn_gaf_re = re.compile(r'(^|[\t><])([^|\t\n<>#]+#[^|\t\n<>#]+)#')

def gaf_from_pansn(names, gaf_path, out_path):
    """ the inverse of gaf_to_pansn(): rewrite a published GAF's PanSN SAMPLE#HAP#CONTIG names back
    to cactus's id=EVENT|CONTIG, so it can be resolved against a cactus-named GFA again.

    a prefix that is not one of the seqfile's genomes is left alone rather than mangled, which
    covers both an already-cactus-named GAF (from a cactus old enough to have published one) and
    any contig name that happens to look like a PanSN prefix """
    prefix_map = pansn_to_event_map(names)

    def replace(m):
        event = prefix_map.get(m.group(2))
        if event is None:
            return m.group(0)
        return '{}id={}|'.format(m.group(1), event)

    with open(gaf_path, 'r') as in_file, open(out_path, 'w') as out_file:
        for line in in_file:
            out_file.write(pansn_gaf_re.sub(replace, line))

def split_gaf_by_event(job, gaf_id, names, gaf_path):
    """ split a published (merged) GAF into one file per genome, returning {event: file id} """
    work_dir = job.fileStore.getLocalTempDir()
    local_gaf_path = os.path.join(work_dir, 'extend.gaf.gz' if gaf_path.endswith('.gz') else 'extend.gaf')
    job.fileStore.readGlobalFile(gaf_id, local_gaf_path)
    shard_dir = os.path.join(work_dir, 'shards')
    os.makedirs(shard_dir)

    shard_paths, dropped = split_gaf_file_by_event(local_gaf_path, names, shard_dir)

    if dropped:
        RealtimeLogger.info('Ignoring mappings in {} for {}name(s) not in the seqfile: {}'.format(
            gaf_path, 'at least ' if len(dropped) >= MAX_DROPPED_NAMES_REPORTED else '{} '.format(len(dropped)),
            ' '.join(sorted(dropped))))
    RealtimeLogger.info('Reusing mappings for {} genome(s) from {}'.format(len(shard_paths), gaf_path))

    return {event: job.fileStore.writeGlobalFile(shard_path) for event, shard_path in shard_paths.items()}

# how many unrecognised names a split reports before it stops collecting them.  the genome part of
# a name is normally one of a handful, but a name in neither naming falls back to its contig
MAX_DROPPED_NAMES_REPORTED = 20

def split_gaf_file_by_event(gaf_path, names, shard_dir):
    """ split a published (merged) GAF into one file per genome, returning ({event: path}, dropped names).

    the merged GAF is a concatenation of the per-genome files, so this puts each genome's mappings
    back exactly as minigraph_map_one() left them -- which is what lets the reused mappings be
    re-derived by the very same code that produced them in the first place.

    genomes in the GAF that are not in names are dropped, not an error: --refFromGFA legitimately
    takes the reference out of the sequence map.  cactus-minigraph --inGFA is where a genome
    missing from the seqfile is caught, because there it is unrecoverable """
    prefix_map = pansn_to_event_map(names)

    def event_of(query_name):
        """ the seqfile event a GAF query column belongs to, or None """
        if query_name.startswith('id='):
            barpos = query_name.find('|')
            event = query_name[3:barpos] if barpos > 3 else None
            return event if event in names else None
        hashpos = query_name.find('#')
        if hashpos < 0:
            return None
        hashpos2 = query_name.find('#', hashpos + 1)
        if hashpos2 < 0:
            return None
        return prefix_map.get(query_name[:hashpos2])

    shard_paths = {}
    dropped = set()
    # the merged GAF groups each genome's records together, so one handle at a time is enough.
    # append mode keeps it correct even if some other producer interleaved them
    cur_event, cur_file = None, None
    opener = gzip.open if gaf_path.endswith('.gz') else open
    with opener(gaf_path, 'rt') as gaf_file:
        for line in gaf_file:
            tab = line.find('\t')
            if tab < 0:
                continue
            query_name = line[:tab]
            event = event_of(query_name)
            if event is None:
                # report the genome part of the name, in whichever naming it is in.  a name that
                # carries no genome at all falls back to the whole contig, so this is bounded by
                # contig count rather than genome count and needs a cap of its own
                if len(dropped) < MAX_DROPPED_NAMES_REPORTED:
                    dropped.add(query_name[3:query_name.find('|')] if query_name.startswith('id=') and '|' in query_name
                                else query_name.split('#')[0])
                continue
            if event != cur_event:
                if cur_file:
                    cur_file.close()
                if event not in shard_paths:
                    shard_paths[event] = os.path.join(shard_dir, '{}.gaf'.format(event))
                cur_file = open(shard_paths[event], 'a')
                cur_event = event
            cur_file.write(line)
    if cur_file:
        cur_file.close()

    return shard_paths, dropped

def gaf_to_pansn(gaf_path, out_path):
    """ rewrite cactus's internal id=EVENT|CONTIG names as PanSN SAMPLE#HAP#CONTIG

    minigraph_gfa_to_pansn() converts the GFA on the way out, but the GAF was left in cactus
    naming, so the two files cactus publishes sat in different namespaces.  Nothing inside
    cactus was affected -- gaf2unstable runs above on the cactus-named GAF against the
    cactus-named GFA -- but a reader of the published pair could not resolve a GAF path
    segment against the graph it came from, nor recover sample and haplotype without knowing
    that cactus encodes the haplotype as a .N suffix.

    Both the query column and the path column are rewritten: id=EVENT| appears in each and
    nowhere else in the record, so a single substitution over the line covers it. """
    with open(gaf_path, 'r') as in_file, open(out_path, 'w') as out_file:
        for line in in_file:
            out_file.write(gaf_pansn_re.sub(lambda m: m.group(1) + event_to_pansn_prefix(m.group(2)) + '#', line))

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

    return stable_gaf_to_paf(job, config, gaf_path, gfa_path)

# how many reused GAF records the up-front check resolves before trusting the rest
GAF_REUSE_CHECK_RECORDS = 1000

def check_reusable_gaf(job, config, gaf_file_id, gfa_file_id, genome_names, gaf_path, gfa_path):
    """ Resolve the first few reused mappings against the graph, and fail with a diagnosis if they
    do not fit.

    The pair has to come from the same graphmap run.  The trap is that a --mgSplit run publishes
    <outName>.sv.gfa.gz and <outName>.gaf.gz that look like a pair and are not: the GAF is the
    whole-genome pass against the reference-only first-pass graph, while the GFA is the merged
    per-chromosome graphs.  Same stable coordinates, different node boundaries, so gaf2unstable
    fails on the tiling rather than on a name and the mismatch is not obvious from the error. """
    work_dir = job.fileStore.getLocalTempDir()
    gfa_local = os.path.join(work_dir, 'mg.gfa')
    job.fileStore.readGlobalFile(gfa_file_id, gfa_local)
    pansn_head = os.path.join(work_dir, 'check.pansn.gaf')
    head = os.path.join(work_dir, 'check.gaf')
    job.fileStore.readGlobalFile(gaf_file_id, pansn_head + '.full')
    with open(pansn_head, 'w') as out_file:
        opener = gzip.open if gaf_path.endswith('.gz') else open
        with opener(pansn_head + '.full', 'rt') as in_file:
            for i, line in enumerate(in_file):
                if i >= GAF_REUSE_CHECK_RECORDS: break
                out_file.write(line)
    os.remove(pansn_head + '.full')
    gaf_from_pansn(genome_names, pansn_head, head)

    try:
        cactus_call(parameters=['gaf2unstable', head, '-g', gfa_local, '-o', os.path.join(work_dir, 'lens.tsv')],
                    outfile=os.path.join(work_dir, 'check.unstable.gaf'), job_memory=job.memory)
    except RuntimeError as e:
        # name the genomes each side actually talks about: the mismatch above shows up as a GAF
        # whose paths only ever name the references while the graph is full of samples
        step_genomes = set()
        with open(head) as in_file:
            for line in in_file:
                toks = line.split('\t')
                if len(toks) > 5:
                    for step in gaf_step_re.findall(toks[5]):
                        name = step[1:].rsplit(':', 1)[0] if ':' in step else step[1:]
                        step_genomes.add(name[3:name.find('|')] if name.startswith('id=') and '|' in name else name)
        raise RuntimeError(
            'The mappings in {} do not resolve against {}: the two must come from the same graphmap run.\n'
            'The first {} records only ever walk through: {}.\n'
            'A --mgSplit run is the usual way to get a mismatched pair that looks like a matching one: its '
            '<outName>.gaf.gz is the whole-genome pass against the reference-only first-pass graph, while its '
            '<outName>.sv.gfa.gz is the merged per-chromosome graphs.  Drop --inGAF to map every genome '
            'against the extended graph instead -- --inGFA still saves the construction, which is the '
            'expensive half.\nUnderlying error: {}'.format(
                gaf_path, gfa_path, GAF_REUSE_CHECK_RECORDS,
                ' '.join(sorted(step_genomes)[:12]) or '(nothing)', e))

    RealtimeLogger.info('Reused mappings from {} resolve against the graph'.format(gaf_path))

def translate_gaf_one(job, config, event_name, gaf_file_id, gfa_file_id, genome_names):
    """ Re-derive one genome's PAF from mappings it already has, against a (possibly extended) graph.

    minigraph GAF is in stable coordinates -- rGFA SN/SO names and offsets -- which adding genomes
    to a graph does not change: new nodes are appended and existing ones are only ever split, so
    the stable sequence a node covers stays exactly where it was.  That makes re-deriving the PAF a
    matter of running the same gaf2unstable/gaf2paf chain minigraph_map_one() runs, against the new
    graph, which gaf2unstable resolves into the new (finer) node ids for free.

    the graph the genome was originally mapped to is not needed, and neither is minigraph """

    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    gaf_path = os.path.join(work_dir, "{}.gaf".format(event_name))
    job.fileStore.readGlobalFile(gfa_file_id, gfa_path)

    # the published GAF is PanSN, but gaf2unstable resolves it against the cactus-named GFA.
    # gaf_from_pansn passes through anything already in cactus naming.  note the input cannot be
    # named <gaf>.pansn: that is where stable_gaf_to_paf() writes the copy it publishes
    in_gaf_path = os.path.join(work_dir, "{}.in.gaf".format(event_name))
    job.fileStore.readGlobalFile(gaf_file_id, in_gaf_path)
    gaf_from_pansn(genome_names, in_gaf_path, gaf_path)

    # the reused GAF is published back (PanSN in, PanSN out, unchanged), so a run that extends a
    # pangenome can itself be extended
    return stable_gaf_to_paf(job, config, gaf_path, gfa_path, regranulated=True)

# a GAF path step: an orientation mark followed by a name that runs to the next mark
gaf_step_re = re.compile(r'[<>][^<>]+')

def trim_unstable_gaf(gaf_path, out_path, node_lengths_path):
    """ drop the path steps of each record that carry none of its alignment, moving the path
    offsets along with them.

    gaf2paf reads a record's path start as an offset into its *first* step.  That holds for a GAF
    gaf2unstable resolved against the graph it was mapped to, where each of minigraph's stable
    steps is one node.  Against a graph that has since been extended, the same stable step resolves
    into the several finer nodes it was split into, and the offset can now reach past the first of
    them -- which gaf2paf asserts on rather than handles.

    The steps it reaches past hold no aligned bases, and gaf2paf emits nothing for them even when it
    does cope, so taking them off the front and back restores the shape gaf2paf expects without
    changing the alignment at all.  When nothing needs trimming -- every mapping that was made
    against the graph it is being resolved against -- every line is passed through untouched.

    Nothing in cactus-gfa-tools does this today.  gaffilter's rebase_path() is the same algorithm in
    C++, but it runs only on the records -t actually trims, and only when trimming is on at all, so
    it never sees the reuse case.  gaf2unstable, which is the thing changing the granularity, is
    where this belongs if it ever moves. """
    node_len = {}
    with open(node_lengths_path) as lengths_file:
        for line in lengths_file:
            toks = line.split()
            if len(toks) >= 2:
                node_len[toks[0]] = int(toks[1])

    def first_step_len(path):
        end = path.find('>', 1)
        alt = path.find('<', 1)
        if alt != -1 and (end == -1 or alt < end):
            end = alt
        return node_len[path[1:] if end == -1 else path[1:end]]

    def last_step_len(path):
        start = max(path.rfind('>'), path.rfind('<'))
        return node_len[path[start + 1:]]

    trimmed_records = 0
    with open(gaf_path) as in_file, open(out_path, 'w') as out_file:
        for line in in_file:
            toks = line.rstrip('\n').split('\t')
            if len(toks) < 12 or not toks[5] or toks[5][0] not in '<>':
                out_file.write(line)
                continue
            path, path_len, path_start, path_end = toks[5], int(toks[6]), int(toks[7]), int(toks[8])
            # the overwhelming majority of records need nothing done, and deciding that needs only
            # the two end steps: with the offsets inside them, nothing in between can be outside
            if path_start < first_step_len(path) and path_end > path_len - last_step_len(path):
                out_file.write(line)
                continue
            steps = gaf_step_re.findall(path)
            lo, hi = 0, len(steps)
            while lo < hi - 1 and node_len[steps[lo][1:]] <= path_start:
                dropped = node_len[steps[lo][1:]]
                path_start -= dropped
                path_end -= dropped
                path_len -= dropped
                lo += 1
            while hi - 1 > lo and path_len - node_len[steps[hi - 1][1:]] >= path_end:
                path_len -= node_len[steps[hi - 1][1:]]
                hi -= 1
            if lo == 0 and hi == len(steps):
                # nothing was outside the alignment after all, so the record stands as it is
                out_file.write(line)
                continue
            toks[5], toks[6], toks[7], toks[8] = ''.join(steps[lo:hi]), str(path_len), str(path_start), str(path_end)
            out_file.write('\t'.join(toks) + '\n')
            trimmed_records += 1

    return trimmed_records

def stable_gaf_to_paf(job, config, gaf_path, gfa_path, regranulated=False):
    """ Turn a stable-coordinate (ie minigraph output) GAF into the node-coordinate PAF cactus
    consumes, returning (published PanSN gaf id, paf id).  Shared by mapping and by reuse of an
    existing mapping, so that the two produce identical output for identical input.

    regranulated says the GAF was made against a coarser version of this graph, so its path offsets
    have to be brought back inside their first and last steps -- see trim_unstable_gaf().  It is a
    no-op when the graph has not changed, but it is only asked for on the reuse path so that
    mapping keeps running exactly the commands it always has """

    xml_node = findRequiredNode(config.xmlRoot, "graphmap")

    # convert the gaf into unstable gaf (targets are node sequences)
    # note: the gfa needs to be uncompressed for this tool to work
    mg_lengths_path = gfa_path + '.node_lengths.tsv'
    unstable_gaf_path = gaf_path + '.unstable'
    cmd = ['gaf2unstable', gaf_path, '-g', gfa_path, '-o', mg_lengths_path]

    # optional gaf overlap filter
    overlap_ratio = getOptionalAttrib(xml_node, "GAFOverlapFilterRatio", typeFn=float, default=0)
    length_ratio = getOptionalAttrib(xml_node, "GAFOverlapFilterMinLengthRatio", typeFn=float, default=0)
    min_block = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minGAFBlockLength", typeFn=int, default=0)
    min_mapq = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minMAPQ", typeFn=int, default=0)
    min_ident = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minIdentity", typeFn=float, default=0)    
    overlap_trim = getOptionalAttrib(xml_node, "GAFOverlapFilterTrim", typeFn=bool, default=False)
    trim_edge = getOptionalAttrib(xml_node, "GAFOverlapFilterTrimEdge", typeFn=int, default=5000)
    trim_min_mapq = getOptionalAttrib(xml_node, "GAFOverlapFilterTrimMinMAPQ", typeFn=int, default=20)
    if overlap_ratio:
        if overlap_trim:
            # GAFOverlapFilterMinLengthRatio exists because deleting a record is expensive: it stops
            # a small overlap from destroying a whole one.  A trim costs only the contested span, so
            # the guard has nothing left to protect and only leaves small conflicts between long
            # contigs unadjudicated -- measured on a 12-sample chr15 graph, trimming with it at 0.25
            # doubly places 409,584 bp on one segmental-duplication pair that 0 does not.  There is
            # no setting of it that helps while trimming, so it is not a knob here.
            length_ratio = 0
        overlap_cmd = ['gaffilter', '-', '-r', str(overlap_ratio), '-m', str(length_ratio), '-q', str(min_mapq),
                       '-b', str(min_block), '-i', str(min_ident)]
        if overlap_trim:
            # cut the contested span out of a losing record instead of deleting the record whole.
            # -l: the unstable GAF names bare nodes, so gaffilter needs their lengths to shorten a
            # path.  gaf2unstable writes that file in full before it emits its first GAF line, so
            # it is complete by the time gaffilter, which reads all of its input first, opens it.
            # gaffilter's -g/--close-holes are deliberately not exposed: --close-holes is off and
            # cannot be justified (no claimant to a hole can meet the bar -r sets), and with it off
            # -g provably does not change the output.
            overlap_cmd += ['-t', '-e', str(trim_edge), '-Q', str(trim_min_mapq),
                            '-l', mg_lengths_path]
        cmd = [cmd, overlap_cmd]
    try:
        cactus_call(parameters=cmd, outfile=unstable_gaf_path, job_memory=job.memory)
    except RuntimeError as e:
        if not regranulated:
            raise
        # gaf2unstable asserts, rather than reporting, when a record names stable sequence the
        # graph does not have.  Mapping cannot reach that -- the GAF came from this graph -- but
        # reuse can, and the assertion on its own says nothing about which input is wrong
        raise RuntimeError('Failed to resolve reused mappings against this graph. If the GAF names sequence the graph '
                           'does not have, it was made against a different pangenome: the GAF must come from the run '
                           'that produced the graph being reused. Underlying error: {}'.format(e))

    if regranulated:
        trimmed_path = unstable_gaf_path + '.trimmed'
        trimmed = trim_unstable_gaf(unstable_gaf_path, trimmed_path, mg_lengths_path)
        RealtimeLogger.info('Moved the path offsets of {} reused GAF record(s) back inside their end steps'.format(trimmed))
        os.replace(trimmed_path, unstable_gaf_path)

    # convert the unstable gaf into unstable paf, which is what cactus expects
    # also tack on the unique id to the target column
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    unstable_paf_path = unstable_gaf_path + '.paf'
    unstable_paf_cmd = [['gaf2paf', unstable_gaf_path, '-l', mg_lengths_path],
                        ['awk', 'BEGIN{{OFS=\"	\"}} {{$6="id={}|"$6; print}}'.format(graph_event)]]
    cactus_call(parameters=unstable_paf_cmd, outfile=unstable_paf_path, job_memory=job.memory)

    # the gaf is published as-is, so put it in PanSN to match the exported minigraph GFA
    # (the filtering chain above has already consumed gaf_path in its cactus-named form)
    pansn_gaf_path = gaf_path + '.pansn'
    gaf_to_pansn(gaf_path, pansn_gaf_path)
    # nothing reads the cactus-named copy after this, and the job's disk request did not grow
    # to hold both
    os.remove(gaf_path)

    # return the stable gaf (minigraph output) and the unstable paf
    return job.fileStore.writeGlobalFile(pansn_gaf_path), job.fileStore.writeGlobalFile(unstable_paf_path)

# What is left of a merge_pafs job once its staging is accounted for: catFiles runs no command, so
# this is worker startup and the python copy loop.
MERGE_PAF_SECS = 60

def merge_pafs_walltime(merged_bytes, gzip=False, cores=1):
    """ walltime for a merge_pafs job whose output comes to roughly merged_bytes.  The job is all
    I/O -- every input is read out of the jobstore and the concatenation written back -- except
    with gzip=True, which bgzips the result on the way out, threaded when it is given cores.

    io_bytes is 4x rather than 2x because the bytes move twice: once staging in and out of the
    jobstore, and again locally, where catFiles reads every input back and writes the joined
    file.  At 2x the whole-panel merge of the HPRC v2.0 PAF came to under half an hour. """
    secs = MERGE_PAF_SECS
    if gzip:
        secs += merged_bytes / (GZIP_COMPRESS_BYTES_PER_SEC * max(1, cores))
    return cactus_walltime(secs, io_bytes=4*merged_bytes)

def merge_pafs(job, paf_file_id_map, gzip=False, merged_name=None):
    """ merge up some pafs.  merged_name is what the merged file is called on disk: getLocalTempFile()
    would give it an anonymous .tmp, which is all anyone reading the log of the bgzip below would see.

    it is NOT called "name": toil's FunctionWrappingJob pops memory/cores/disk/accelerators/
    preemptible/checkpoint/name out of the kwargs of addChildJobFn before the function is called, and
    uses "name" as the job's unitName.  a parameter with any of those names is silently swallowed """
    paf_paths = [job.fileStore.readGlobalFile(paf_id) for paf_id in paf_file_id_map.values()]
    merged_path = os.path.join(job.fileStore.getLocalTempDir(), merged_name if merged_name else 'merged.paf')
    catFiles(paf_paths, merged_path)
    if gzip:
        cactus_call(parameters=['bgzip', merged_path, '--threads', str(job.cores)])
        merged_path += '.gz'                    
    return job.fileStore.writeGlobalFile(merged_path)

def merge_pafs_sized(job, paf_file_id_map, gzip=False, merged_name=None, gzip_cores=1):
    """ merge_pafs, sized off its inputs.  callers upstream of the mapping jobs hold promises and so
    cannot measure them; by the time this job runs they are resolved.  the merge holds every input
    plus the merged copy, and bgzip then writes a compressed copy alongside.

    gzip_cores is not called "cores" for the reason in merge_pafs above: toil would eat it as this
    job's own resource request and the function would keep its default, leaving the bgzip on one
    thread.  it is passed on as the child's cores, where toil eating it is exactly what we want """
    total_size = sum(paf_id.size for paf_id in paf_file_id_map.values() if paf_id)
    return job.addChildJobFn(merge_pafs, paf_file_id_map, gzip=gzip, merged_name=merged_name,
                             cores=gzip_cores, disk=max(total_size * 3, 2**31),
                             walltime=merge_pafs_walltime(total_size, gzip=gzip, cores=gzip_cores)).rv()

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

# Seconds per GB of fasta for a minimap2 -xasm5 self-alignment.  There is no measurement behind
# this one: <graphmap collapse> defaults to "none", so self_align ran in none of the runs we have
# logs for.  It is minigraph's own mapping rate standing in, and wants replacing with a
# measurement the first time --collapse is used at panel scale.
SELF_ALIGN_SECS_PER_GB = 1200

def self_align_all(job, config, seq_id_map, reference, collapse_mode):
    """ run self-alignment. if reference event given, just run on that, otherwise do all genomes """
    assert collapse_mode in ['reference', 'all', 'nonref']
    assert reference or collapse_mode == 'all'
    root_job = Job(walltime=cactus_walltime())
    job.addChild(root_job)
    events = []
    for event in seq_id_map.keys():
        if collapse_mode == 'all' or (collapse_mode == 'reference' and event == reference) or \
           (collapse_mode == 'nonref' and event != reference):
            events.append(event)
    mg_cores = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "cpu", typeFn=int, default=1)
    paf_dict = {}
    for event in events:
        collapse_job = root_job.addChildJobFn(self_align, config, event, seq_id_map[event],
                                              disk=4*seq_id_map[event].size,
                                              memory=4*seq_id_map[event].size,
                                              cores=mg_cores,
                                              walltime=cactus_walltime(600 + SELF_ALIGN_SECS_PER_GB * seq_id_map[event].size / 1e9,
                                                                       io_bytes=2*seq_id_map[event].size))
        paf_dict[event] = collapse_job.rv()

    # every input here is a promise, so the merge is sized the same way the mapping merges are: by
    # a coordination job that runs once they have resolved and can read their real sizes.  what a
    # self-alignment PAF comes to as a fraction of the fasta it came from has never been measured
    # -- <graphmap collapse> defaults to "none" -- and this way it does not have to be guessed
    merge_paf_job = root_job.addFollowOnJobFn(merge_pafs_sized, paf_dict, walltime=cactus_walltime())
    return merge_paf_job.rv()

def self_align(job, config, seq_name, seq_id):
    """ run minimap2 self alignment on a single genome """
    work_dir = job.fileStore.getLocalTempDir()
    seq_name = seq_name.split('|')[-1]
    fa_path = os.path.join(work_dir, seq_name + '.fa')
    job.fileStore.readGlobalFile(seq_id, fa_path)
    cactus_call(parameters=['samtools', 'faidx', fa_path])
    contigs = []
    with open(fa_path + '.fai', 'r') as fai_file:
        for line in fai_file:
            if line.strip():
                contigs.append(line.split()[0])
    paf_path = os.path.join(work_dir, seq_name + '.self.paf')
    xml_node = findRequiredNode(config.xmlRoot, "graphmap")
    minimap_opts = getOptionalAttrib(xml_node, "minimapCollapseOptions", str, default="")     
    for contig in contigs:
        contig_path = os.path.join(work_dir, '{}.{}.fa'.format(seq_name, contig))
        cactus_call(parameters=['samtools', 'faidx', fa_path, contig, '-o', contig_path])
        cmd = ['minimap2', '-t', str(job.cores)] + minimap_opts.split() + [contig_path, contig_path]
        cactus_call(parameters=cmd, outfile=paf_path, outappend=True)
    return job.fileStore.writeGlobalFile(paf_path)        

def apply_mgsplit_filter_overrides(config_node):
    """ Disable the overlap/block-length/deletion filters in the <graphmap> config.
        These filters tidy up alignment blocks but mangle the signal rgfa-split needs
        when the splitting graph contains only the reference (eg --mgSplit / --refOnly):
        palindromic mappings against a single-haplotype chrY all look similarly sized
        and gaffilter axes them.  rgfa-split has its own interval merging.
        Caller deep-copies the config first if it must be preserved for other branches. """
    graphmap_node = findRequiredNode(config_node, "graphmap")
    graphmap_node.attrib["GAFOverlapFilterRatio"] = "0"
    graphmap_node.attrib["PAFOverlapFilterRatio"] = "0"
    graphmap_node.attrib["minGAFBlockLength"] = "0"
    graphmap_node.attrib["delFilter"] = "-1"

# Seconds per GB of PAF for filter_paf: gaffilter is measured at 1475s on the 34.4 GiB whole-panel
# PAF and 270s worst case on the 2.2 GiB per-chromosome ones, and the python line-by-line pass that
# always runs adds about as much again at ~60 MB/s.  Lives here, next to the job, because
# cactus-graphmap-split and cactus-align both schedule it.
FILTER_PAF_SECS_PER_GB = 60

def filter_paf(job, paf_id, config, reference=None):
    """ run basic paf-filtering.  these are quick filters that are best to do on-the-fly when reading the paf and
        as such, they are called by cactus-graphmap-split and cactus-align, not here

        both callers must pass the same arguments: depending on --noSplit and on how a step-by-step
        run is driven, this may run before the split, before cactus-align, or both, and those three
        need to agree """
    work_dir = job.fileStore.getLocalTempDir()
    paf_path = os.path.join(work_dir, 'mg.paf')
    filter_paf_path = os.path.join(work_dir, 'mg.paf.filter')
    job.fileStore.readGlobalFile(paf_id, paf_path)

    min_block = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minGAFBlockLength", typeFn=int, default=0)
    min_mapq = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minMAPQ", typeFn=int, default=0)
    min_ident = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minIdentity", typeFn=float, default=0)
    min_score = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "minScore", typeFn=float, default=0)
    max_collapse_ratio = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "maxCollapseDistanceRatio", typeFn=float, default=-1)
    RealtimeLogger.info("Running PAF filter with minBlock={} minMAPQ={} minIdentity={} minScore={} maxCollapseDistanceRatio={}".format(min_block, min_mapq, min_ident, min_score, max_collapse_ratio))
    # a contig whose alignments are all filtered out here vanishes from everything downstream: it
    # reaches neither rgfa-split's coverage map nor any chromosome, so it is never assigned, never
    # binned ambiguous, and never logged.  count what goes in and what survives so it can be said
    # out loud rather than inferred from a missing output fasta at the end of the run.
    query_line_counts = {}
    query_kept_counts = {}
    with open(paf_path, 'r') as paf_file, open(filter_paf_path, 'w') as filter_paf_file:
        for line in paf_file:
            toks = line.split('\t')
            query_name = toks[0]
            query_line_counts[query_name] = query_line_counts.get(query_name, 0) + 1
            target_name = toks[5]
            is_ref = reference and query_name.startswith('id={}|'.format(reference)) and query_name != target_name
            mapq = int(toks[11])
            query_len = int(toks[1])
            ident = float(toks[9]) / (float(toks[10]) + 0.00000001)
            bl = None
            score = None
            collapse_ratio = None
            for tok in toks[12:]:
                # this is a special tag that was written by gaf2paf in order to preserve the original gaf block length
                # we use it to be able to filter by the gaf block even after it's been broken in the paf
                if tok.startswith('gl:i:'):
                    bl = int(tok[5:])
                # we can also get the identity of the parent gaf block
                # (gaf2paf writes this as a float, gi:f:, not gi:i:)
                if tok.startswith('gi:f:'):
                    ident = min(ident, float(tok[5:]))
                if tok.startswith('AS:i:'):
                    score = int(tok[5:])
            if query_name == target_name and max_collapse_ratio >= 0:
                # compute the distance between the minimap2 self alignment intervals
                # in order to see if they are too far apart via the max ratio
                query_start, query_end = int(toks[2]), int(toks[3])
                target_start, target_end = int(toks[7]), int(toks[8])
                block_length = int(toks[10])
                dist = target_start - query_end if query_start < target_start else query_start - target_end
                if block_length > 0 and dist > 0:
                    collapse_ratio = dist / block_length
            if is_ref or (mapq >= min_mapq and (bl is None or query_len <= min_block or bl >= min_block) and ident >= min_ident and \
                          (score is None or score >= min_score) and (collapse_ratio is None or collapse_ratio <= max_collapse_ratio)):
                filter_paf_file.write(line)
                query_kept_counts[query_name] = query_kept_counts.get(query_name, 0) + 1

    eliminated = sorted(q for q, n in query_line_counts.items() if not query_kept_counts.get(q))
    if eliminated:
        RealtimeLogger.warning(
            'PAF filter removed every alignment of {} query contig(s), which will therefore be '
            'absent from all output graphs without appearing in any chromosome or in the ambiguous '
            'bin. The usual cause is minGAFBlockLength={} exceeding the contig\'s longest alignment '
            'block (contigs shorter than that are exempt). Contigs: {}'.format(
                len(eliminated), min_block, ', '.join(eliminated[:20])))

    overlap_ratio = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "PAFOverlapFilterRatio", typeFn=float, default=0)
    length_ratio = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "PAFOverlapFilterMinLengthRatio", typeFn=float, default=0)
    allow_collapse = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "collapse", typeFn=str, default="none") != "none"

    if overlap_ratio and not allow_collapse:
        overlap_filter_paf_path = filter_paf_path + ".overlap"
        cactus_call(parameters=['gaffilter', filter_paf_path, '-p', '-r', str(overlap_ratio), '-m', str(length_ratio),
                                '-b', str(min_block), '-q', str(min_mapq), '-i', str(min_ident)],
                    outfile=overlap_filter_paf_path, job_memory=job.memory)
        filter_paf_path = overlap_filter_paf_path

    return job.fileStore.writeGlobalFile(filter_paf_path)    

def filter_paf_deletions(job, paf_id, gfa_id, max_deletion, filter_threshold, filter_query_size_threshold, max_remove=None, min_support=None):
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
    if max_remove is not None and max_remove >= 0:
        # the absolute budget for deletions fewer than min_support contigs assert (see the config)
        filter_paf_cmd += ['-M', str(max_remove)]
        if min_support is not None:
            filter_paf_cmd += ['-S', str(min_support)]
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
