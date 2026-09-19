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
from cactus.refmap.pangenome_exclusions import event_to_pansn_prefix
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
from toil.job import PromisedRequirement
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
    parser.add_argument("--inGFA", type=str, default=None,
                        help="Start from this existing minigraph GFA (as made by a previous cactus-minigraph or cactus-pangenome run) "
                        "instead of building from scratch. Only the seqFile genomes that are not already in the graph get added, in "
                        "mash-distance order among themselves; if there are none, the graph is passed through untouched. The seqFile "
                        "must still contain every genome in the graph: genomes cannot be removed from a minigraph")
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

            in_gfa_id = None
            if options.inGFA:
                if options.batch:
                    raise RuntimeError('--inGFA cannot be used with --batch')
                if options.refOnly:
                    raise RuntimeError('--inGFA cannot be used with --refOnly')
                if '://' not in options.inGFA:
                    options.inGFA = os.path.abspath(options.inGFA)
                in_gfa_id = toil.importFile(makeURL(options.inGFA))

            # maps name -> input_seq_id_map, input_seq_order
            input_dict = minigraph_construct_import_sequences(options, config_wrapper, input_seqfiles, toil)
                
            # output_dict:  chrom-> (gfa_id, pansn_gfa_id, uncollapsed_pansn_gfa_id, collapse_report_id, train_id)
            output_dict = toil.start(Job.wrapJobFn(minigraph_construct_batch_workflow, options, config_node, input_dict, options.outputGFA,
                                                   in_gfa_id=in_gfa_id))

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
                    catFiles([os.path.join(seq, subSeq) for subSeq in sorted(os.listdir(seq))], tmpSeq)
                    seq = tmpSeq
                seq = makeURL(seq)
                input_seq_id_map[genome] = file_store.importFile(seq)
            elif genome in input_seq_order:
                input_seq_order.remove(genome)

        input_dict[chrom] = (input_seq_id_map, input_seq_order)
        
    return input_dict

def collapse_artifact_paths(gfa_path):
    """ the two side artifacts a collapse run leaves beside <x>.sv.gfa.gz:
    <x>.sv.uncollapsed.gfa.gz (the graph as minigraph built it) and <x>.sv.collapse.tsv """
    uncollapsed_path = gfa_path.replace('.gfa', '.uncollapsed.gfa') if '.gfa' in gfa_path \
                       else gfa_path + '.uncollapsed'
    report_path = gfa_path.replace('.gfa.gz', '').replace('.gfa', '') + '.collapse.tsv'
    return uncollapsed_path, report_path

def export_collapse_artifacts(exporter, gfa_path, uncollapsed_pansn_gfa_id, collapse_report_id):
    """ write the pre-collapse graph and the per-call report beside the collapsed graph.  `exporter`
    is anything with exportFile: a Toil object at the top level, or job.fileStore inside a job """
    uncollapsed_path, report_path = collapse_artifact_paths(gfa_path)
    if uncollapsed_pansn_gfa_id:
        exporter.exportFile(uncollapsed_pansn_gfa_id, makeURL(uncollapsed_path))
    if collapse_report_id:
        exporter.exportFile(collapse_report_id, makeURL(report_path))

def export_minigraph_construct_output(options, input_seqfiles, output_dict, toil):
    if options.batch:
        chrom_file_path = os.path.join(options.outputGFA, 'chromfile.mg.txt')
        if chrom_file_path.startswith('s3://'):
            chrom_file_temp_path = getTempFile()
        else:
            chrom_file_temp_path = chrom_file_path                    
        chromfile = open(chrom_file_temp_path, 'w')
    for chrom, output_ids in output_dict.items():
        gfa_id, pansn_gfa_id, uncollapsed_pansn_gfa_id, collapse_report_id, train_id = output_ids
        if options.batch:
            gfa_path = os.path.join(options.outputGFA, chrom + '.sv.gfa.gz')
        else:
            gfa_path = options.outputGFA
        if train_id:
            train_path = gfa_path.replace('.gfa.gz', '.gfa').replace('.gfa', '.train')
        else:
            train_path = None
        #export the gfa
        # with --mgSplitWholeGenomeRef this graph is still whole-genome: it only becomes
        # chromosome-only after the post-mapping prune, which writes it to this same path
        # (export_pruned_minigraph_gfa_wrapper in cactus_pangenome.py).  the chromfile below is
        # still written now -- downstream only ever reads its .train column
        if not getattr(options, 'mgSplitWholeGenomeRef', False):
            # pansn_gfa_id is already the collapsed graph when collapseInversions is on, so the
            # main export is unconditional; the pre-collapse graph and the per-call report are
            # extra artifacts kept beside it for comparison
            toil.exportFile(pansn_gfa_id, makeURL(gfa_path))
            export_collapse_artifacts(toil, gfa_path, uncollapsed_pansn_gfa_id, collapse_report_id)
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

def minigraph_construct_batch_workflow(job, options, config_node, input_dict, gfa_path, sanitize=True,
                                       construct_ref_id_map=None, in_gfa_id=None):
    """ run the construction workflow on individual chromosomes.  construct_ref_id_map, if given,
    swaps the whole-genome reference fastas in for the chromosome's own slice of them (--mgSplit
    --mgSplitWholeGenomeRef).  the merge happens here rather than inside minigraph_construct_workflow
    because both dicts are resolved job arguments at this point, whereas the sanitized map down there
    can still be an unresolved promise """
    output_dict = {}
    for chrom, input_info in input_dict.items():
        seq_id_map, seq_order = input_info
        construct_seq_id_map = None
        if construct_ref_id_map:
            construct_seq_id_map = dict(seq_id_map)
            construct_seq_id_map.update(construct_ref_id_map)
        if options.batch:
            gfa_path = os.path.join(options.outputGFA, '{}.gfa.gz'.format(chrom))
        else:
            gfa_path = options.outputGFA
        mgwf_job = job.addChildJobFn(minigraph_construct_workflow, options, config_node, seq_id_map, seq_order, gfa_path, sanitize,
                                     construct_seq_id_map=construct_seq_id_map, in_gfa_id=in_gfa_id)
        output_dict[chrom] = mgwf_job.rv()
    return output_dict
                                    
def minigraph_construct_workflow(job, options, config_node, seq_id_map, seq_order, gfa_path, sanitize=True,
                                 construct_seq_id_map=None, in_gfa_id=None):
    """ minigraph can handle bgzipped files but not gzipped; so unzip everything in case before running

    construct_seq_id_map, when given, replaces seq_id_map for the graph construction alone.  it is how
    --mgSplitWholeGenomeRef hands each chromosome the whole reference: mash sorting and last-training
    (and the size estimate for last-training's own job) keep using the chromosome's own much smaller
    reference slice, since training against a whole genome would be both expensive and
    cross-chromosome contaminated -- and would silently produce no model at all, as last_train()
    requires its partner sequence to be at least half the size of the database it trains against.
    The construction job itself is sized off the substituted map, since that is what it runs on.

    with a graph to extend, which genomes still need constructing is not known until that graph's
    SN tags have been read, so the rest of the workflow is deferred behind the job that reads them """
    if not in_gfa_id:
        return minigraph_construct_run(job, options, config_node, seq_id_map, seq_order, gfa_path, sanitize,
                                       construct_seq_id_map=construct_seq_id_map)

    # the renaming pass decompresses the GFA before bgzipping it back up, so it needs room for
    # the raw copy (reckoned at 10x, as elsewhere) on top of the compressed input and output
    rename_job = job.addChildJobFn(minigraph_gfa_from_pansn, set(seq_id_map.keys()), options.inGFA, in_gfa_id,
                                   disk=in_gfa_id.size*12)
    run_job = rename_job.addFollowOnJobFn(minigraph_construct_run, options, config_node, seq_id_map, seq_order, gfa_path,
                                          sanitize, construct_seq_id_map,
                                          rename_job.rv(0), rename_job.rv(1), in_gfa_id)
    # all five slots: minigraph_construct_run returns
    # (gfa, pansn_gfa, uncollapsed_pansn_gfa, collapse_report, train).  Truncating here to three
    # left cactus_pangenome's rv(3)/rv(4) reading off the end of the tuple, which Toil reports as
    # "IndexError: tuple index out of range" from _fulfillPromises, nowhere near the cause.
    return run_job.rv(0), run_job.rv(1), run_job.rv(2), run_job.rv(3), run_job.rv(4)

def minigraph_construct_run(job, options, config_node, seq_id_map, seq_order, gfa_path, sanitize=True,
                            construct_seq_id_map=None, seed_gfa_id=None, seed_events=None, seed_pansn_gfa_id=None):
    assert type(options.reference) is list
    # the substituted map is a plain dict here, but sanitized_seq_id_map below is a promise when
    # sanitize is on, so the two can't be reconciled in this job
    assert not (construct_seq_id_map and sanitize)
    ref_size = seq_id_map[options.reference[0]].size
    # the PanSN rename at the end of construction has to resolve every SN tag in the finished
    # graph, which on the extend path is more genomes than minigraph is being given
    graph_names = set(seq_id_map.keys())
    # last-training is over fastas, so it wants every genome, not just the ones still to construct
    train_seq_id_map, train_seq_order = seq_id_map, seq_order
    if seed_events is not None:
        if options.reference[0] not in seed_events:
            # it would otherwise be constructed in last, at the highest rGFA rank rather than rank 0,
            # and every rank-0 assumption downstream would be reading the wrong genome
            raise RuntimeError('Reference {} is not in {}, whose genomes are: {}. A graph can only be reused with the '
                               'reference it was built on'.format(options.reference[0], options.inGFA,
                                                                  ' '.join(sorted(seed_events))))
        # everything already in the seed graph is left alone: minigraph only gets the genomes that
        # are new to it, appended after the ones the graph was built from.  the reference is kept in
        # the sequence map (but not the order) because the mash sort below still sketches against it
        seq_order = [seq for seq in seq_order if seq not in seed_events]
        seq_id_map = {name: fa_id for name, fa_id in seq_id_map.items()
                      if name not in seed_events or name == options.reference[0]}
        if seq_order:
            RealtimeLogger.info('Extending the {} genomes in {} with {}: {}'.format(
                len(seed_events), options.inGFA, len(seq_order), ' '.join(seq_order)))
        else:
            # the seqfile asks for exactly the genomes the graph already holds, so there is nothing
            # to construct and the run resumes from it.  the graph handed back is the one --inGFA
            # was given, but its compression follows the *input* name while everything downstream
            # reads the output name to decide whether to unzip, so it is re-emitted to match
            RealtimeLogger.info('Resuming from the {} genomes in {}: nothing left to construct'.format(
                len(seed_events), options.inGFA))
            match_job = job.addChildJobFn(match_gfa_compression, seed_gfa_id, seed_pansn_gfa_id,
                                          options.inGFA, gfa_path,
                                          disk=12 * (seed_gfa_id.size if hasattr(seed_gfa_id, 'size') else 0))
            # last_train reads fastas, not the graph, so resuming is no reason to skip it: without
            # this --lastTrain would quietly fall back to the default scoring matrix
            train_id = None
            if options.lastTrain and len(train_seq_id_map) > 1:
                train_job = job.addChildJobFn(last_train, config_node, train_seq_order, train_seq_id_map,
                                              ref_name=options.reference[0],
                                              cores=options.mgCores, disk=8*ref_size,
                                              memory=cactus_clamp_memory(max(8*ref_size, 12*10**9)))
                train_id = train_job.rv()
            # same five slots as the constructing path below.  nothing was constructed, so there
            # is no pre-collapse graph and no report: the graph handed back is the seed as it was.
            return match_job.rv(0), match_job.rv(1), None, None, train_id
    else:
        assert options.reference[0] == seq_order[0]
    if options.refOnly:
        refonly_seq_id_map = {}
        refonly_seq_order = []
        for seq in seq_order:
            if seq in options.reference:
                refonly_seq_order.append(seq)
                refonly_seq_id_map[seq] = seq_id_map[seq]
        seq_id_map, seq_order = refonly_seq_id_map, refonly_seq_order
        if construct_seq_id_map:
            construct_seq_id_map = {seq: construct_seq_id_map[seq] for seq in refonly_seq_order}
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
        sort_job = sanitize_job.addFollowOnJobFn(sort_minigraph_input_with_mash, options, config_node, sanitized_seq_id_map, seq_order,
                                                 ref_name=options.reference[0] if seed_events is not None else None)
        seq_order = sort_job.rv()
        prev_job = sort_job
    else:
        prev_job = sanitize_job
    minigraph_job = prev_job.addFollowOnJobFn(minigraph_construct_in_batches, options, config_node,
                                              construct_seq_id_map if construct_seq_id_map else sanitized_seq_id_map,
                                              seq_order, gfa_path,
                                              whole_genome_ref=bool(construct_seq_id_map),
                                              seed_gfa_id=seed_gfa_id, graph_names=graph_names)

    # optionally rewrite inverted alleles stored as novel sequence into inversion edges.  when this
    # runs, the collapsed graph REPLACES the constructed one in both namings: it is what graphmap,
    # rgfa-split and the join all have to see, since a graph they map against that still carries the
    # uncollapsed alt would defeat the point.  the originals are returned alongside so the export can
    # keep them for comparison.
    #
    # rgfa-collapse needs the PanSN graph, not the cactus-named one: it locates sites with
    # `vg snarls -n -P <reference sample>`, and cactus names (id=EVENT|CONTIG) carry no PanSN sample
    # for -P to match -- with the wrong -P, vg skips every snarl as a non-reference boundary and the
    # tool silently finds nothing.  So collapse the PanSN graph and rename the result back.
    #
    # skipped on the --mgSplit reference-only first pass: that graph has no non-reference nodes, so
    # there is nothing to collapse and the snarl decomposition would be pure cost.  it still runs on
    # the per-chromosome all-sample graphs, which is where the alleles are.
    uncollapsed_pansn_gfa_id, collapse_report_id = None, None
    if getOptionalAttrib(xml_node, "collapseInversions", typeFn=bool, default=False) and \
       not getattr(options, 'refOnly', False):
        collapse_job = minigraph_job.addFollowOnJobFn(collapse_inversions, options, config_node,
                                                      minigraph_job.rv(1), gfa_path,
                                                      cores=options.mgCores,
                                                      # sized from the graph it collapses, not the
                                                      # reference: rgfa-collapse's cost is almost all
                                                      # fixed.  Over 24 HPRC chromosomes it used
                                                      # 5.9-6.0 GiB of memory and up to 6.57 GiB of
                                                      # disk while the gfa spanned 6000x (0.03 MB to
                                                      # 190 MB).  Both are flat in graph size, so the
                                                      # floors are what is measured -- roughly 2x the
                                                      # observed peak each -- and the per-byte terms
                                                      # only bite above them, as insurance for graphs
                                                      # larger than any seen.
                                                      disk=PromisedRequirement(
                                                          lambda gfa: max(24 * gfa.size, 16 * 2**30),
                                                          minigraph_job.rv(1)),
                                                      memory=PromisedRequirement(
                                                          lambda gfa: cactus_clamp_memory(max(64 * gfa.size, 12 * 2**30)),
                                                          minigraph_job.rv(1)))
        # graph_names, the same set the forward rename uses, and for the same reason its comment
        # gives: it has to resolve every SN tag in the finished graph, which on the --inGFA extend
        # path is more genomes than minigraph is given.  It is captured above before seq_id_map is
        # trimmed to the new genomes, so using seq_id_map here would drop the seed graph's events
        # and minigraph_gfa_from_pansn would assert on the first one it met.
        rename_job = collapse_job.addFollowOnJobFn(minigraph_gfa_from_pansn, graph_names,
                                                   gfa_path, collapse_job.rv(0),
                                                   disk=12*ref_size,
                                                   memory=cactus_clamp_memory(4*ref_size))
        uncollapsed_pansn_gfa_id = minigraph_job.rv(1)
        collapse_report_id = collapse_job.rv(1)
        # rv(0), not rv(): master's minigraph_gfa_from_pansn returns
        # (converted gfa id, set of genomes the graph was built from), and handing the whole tuple
        # on made graphmap fail with "'tuple' object has no attribute 'size'"
        gfa_ids = (rename_job.rv(0), collapse_job.rv(0))
    else:
        gfa_ids = (minigraph_job.rv(0), minigraph_job.rv(1))

    train_id = None
    if options.lastTrain and len(seq_id_map) > 1:
        # note: somehow last training memory overruns don't seem to be detected by slurm so we
        # give 12G at least whenever possible, as --doubleMem won't help...
        last_train_job = prev_job.addFollowOnJobFn(last_train, config_node, seq_order, sanitized_seq_id_map,
                                                   ref_name=options.reference[0] if seed_events is not None else None,
                                                   cores=options.mgCores,
                                                   disk=8*ref_size,
                                                   memory=cactus_clamp_memory(max(8*ref_size, 12*10**9)))
        train_id = last_train_job.rv()
        
    # (cactus-named graph, PanSN graph, the PanSN graph before collapsing or None, per-call report
    #  or None, LAST scoring model or None).  slots 0 and 1 are the graph the pipeline uses.
    return gfa_ids[0], gfa_ids[1], uncollapsed_pansn_gfa_id, collapse_report_id, train_id

def collapse_inversions(job, options, config_node, pansn_gfa_id, gfa_path):
    """ rewrite inverted alleles that minigraph stored as novel sequence into proper inversion
    edges through the reference.

    minigraph sometimes cannot align a haplotype's inverted copy of a region back to the
    reference copy and keeps it as an alt node.  Downstream that is an insertion, not an
    inversion: cactus-align cannot merge two separate backbone nodes and vcfbub flattens the
    giant allele, so the inversion never reaches the VCF.  rgfa-collapse finds these and rewires
    them.  Sites come from vg's snarl decomposition, hence the vg calls here.

    Returns (collapsed_gfa_id, report_id).  The uncollapsed graph is exported alongside so the
    two can be compared.
    """
    work_dir = job.fileStore.getLocalTempDir()
    gzipped = gfa_path.endswith('.gz')
    in_gfa = os.path.join(work_dir, 'in.gfa')
    job.fileStore.readGlobalFile(pansn_gfa_id, in_gfa + ('.gz' if gzipped else ''))
    if gzipped:
        cactus_call(parameters=['bgzip', '-d', '--threads', str(job.cores), in_gfa + '.gz'],
                    work_dir=work_dir)

    xml_node = findRequiredNode(config_node, "graphmap")
    opts = getOptionalAttrib(xml_node, "collapseOptions", str, default="")

    # snarls: -P orients the snarl tree along the reference.  A minigraph rGFA imports with only
    # the reference as a path so this is a no-op there, but it is correct for graphs that carry
    # more, and vg is already a pipeline dependency.
    snarls = os.path.join(work_dir, 'snarls.json')
    cactus_call(parameters=[['vg', 'snarls', '-n', '-P', options.reference[0], '-t', str(job.cores),
                             os.path.basename(in_gfa)],
                            ['vg', 'view', '-Rj', '-']],
                outfile=snarls, work_dir=work_dir)

    out_gfa = os.path.join(work_dir, 'collapsed.gfa')
    report = os.path.join(work_dir, 'collapse.tsv')
    # Bias toward threads-per-minimap2 rather than concurrent invocations.  Runtime is dominated
    # by a few enormous single alignments -- on CHM13 chr9 individual minimap2 calls ran over
    # three hours at -t 2 -- and those cannot be spread across concurrent chunks, so a wide -j
    # leaves cores idle waiting on the tail (measured average utilisation 3.25 of 8 requested at
    # -j 4 -t 2).  A chunk holds ~200 query sequences, so minimap2 keeps a high -t busy on its
    # own.  Note minimap2 runs about two threads more than -t asks for, so this deliberately
    # leaves headroom rather than saturating.
    jobs = min(8, max(1, int(job.cores) // 8))
    threads = max(1, int(job.cores) // jobs)
    cmd = ['rgfa-collapse'] + opts.split() + \
          ['-j', str(jobs), '-t', str(threads), '-r', os.path.basename(report),
           os.path.basename(in_gfa), os.path.basename(snarls)]
    prefix = '[rgfa-collapse-{}]'.format(os.path.basename(gfa_path).replace('.gz', '').replace('.gfa', ''))
    cactus_call(parameters=cmd, outfile=out_gfa, work_dir=work_dir,
                realtimeStderrPrefix=prefix, job_memory=job.memory)

    if gzipped:
        cactus_call(parameters=['bgzip', '--threads', str(job.cores)], infile=out_gfa,
                    outfile=out_gfa + '.gz')
        out_gfa += '.gz'
    return job.fileStore.writeGlobalFile(out_gfa), job.fileStore.writeGlobalFile(report)

def match_gfa_compression(job, gfa_id, pansn_gfa_id, in_path, out_path):
    """ re-emit an unchanged seed graph at the compression its output path asks for.

    every path out of construction bgzips iff the output path ends in .gz, and the stages after it
    read that same suffix to decide whether the file needs unzipping.  Extending a graph by nothing
    is the one case that returns a graph nobody wrote, so it is the one case that can arrive with
    the wrong compression for where it is going. """
    want_gz = out_path.endswith('.gz')
    if want_gz == in_path.endswith('.gz'):
        return gfa_id, pansn_gfa_id

    work_dir = job.fileStore.getLocalTempDir()
    out_ids = []
    for i, file_id in enumerate([gfa_id, pansn_gfa_id]):
        in_gfa_path = os.path.join(work_dir, 'seed.{}.gfa{}'.format(i, '.gz' if not want_gz else ''))
        job.fileStore.readGlobalFile(file_id, in_gfa_path)
        out_gfa_path = os.path.join(work_dir, 'out.{}.gfa{}'.format(i, '.gz' if want_gz else ''))
        if want_gz:
            cactus_call(parameters=['bgzip', '--threads', str(job.cores), '-c', in_gfa_path], outfile=out_gfa_path)
        else:
            cactus_call(parameters=['gzip', '-dc', in_gfa_path], outfile=out_gfa_path)
        out_ids.append(job.fileStore.writeGlobalFile(out_gfa_path))
    return out_ids[0], out_ids[1]

def sort_minigraph_input_with_mash(job, options, config_node, seq_id_map, seq_order, ref_name=None):
    """ Sort the input.

    ref_name names the genome to measure distance against when it is not seq_order[0], as when
    extending an existing minigraph and the order holds only the genomes being added.  It is put
    back at the front for the sort and taken off again on the way out, so the genomes being added
    are still ordered by their distance to the reference """
    trim_ref = ref_name is not None and ref_name not in seq_order
    if trim_ref:
        seq_order = [ref_name] + list(seq_order)
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
            
    return dist_root_job.addFollowOnJobFn(mash_distance_order, options, config_node, seq_order, dist_maps, trim_ref).rv()

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
    
def mash_distance_order(job, options, config_node, seq_order, mash_output_maps, trim_ref=False):
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
            if ref not in mash_order:
                continue
            mash_pos = mash_order.index(ref)
            fix_pos = fixed_order.index(ref)
            assert fix_pos == seq_order.index(ref)
            if mash_pos != fix_pos:
                RealtimeLogger.info('Secondary reference {}, which would have mash rank {}, fixed at input rank {} because minigraphSortReference is disabled'.format(ref, mash_pos, fix_pos))
        mash_order = fixed_order

    return mash_order[1:] if trim_ref else mash_order
            
def minigraph_construct_in_batches(job, options, config_node, seq_id_map, seq_order, gfa_path, whole_genome_ref=False,
                                   seed_gfa_id=None, graph_names=None):
    """ Make minigraph in sequential batches.

    seed_gfa_id is an existing graph to extend: it is fed to the first batch exactly as each batch
    already feeds its output to the next, so extending a graph and constructing one in batches are
    the same operation """

    max_size = max([x.size for x in seq_id_map.values()])
    total_size = sum([x.size for x in seq_id_map.values()])
    # a seed graph can dwarf the genomes being added to it (one sample onto an HPRC-scale graph),
    # so it needs to be in the estimates rather than lost in the headroom the way a batch-to-batch
    # intermediate is
    seed_size = seed_gfa_id.size if seed_gfa_id else 0
    disk = total_size * 2 + seed_size * 12
    mem = cactus_clamp_memory(60 * max_size + int(total_size / 4) + seed_size * 12)
    if whole_genome_ref:
        # with --mgSplitWholeGenomeRef the largest input is the whole reference, so the estimate above
        # is already the whole-genome one, calibrated against a graph carrying far more sample
        # material than one chromosome's.  applying the batch multiple on top would ask for more
        # memory than most machines have.  --mgMemory overrides either way
        RealtimeLogger.info('Sizing whole-genome-reference minigraph_construct for {} at {}'.format(
            os.path.basename(gfa_path), bytes2human(mem)))
    elif options.batch:
        # the memory heuristc seems to drastically underestimate some chromosomes in batch mode...
        mem *= 3
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
    seed_gfa_path = None
    if seed_gfa_id:
        # minigraph_construct() only uses this to name its local copy, but keep the compression
        # suffix honest since that is what says whether the file it reads is bgzipped
        seed_gfa_path = 'extend.gfa.gz' if options.inGFA.endswith('.gz') else 'extend.gfa'
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
                                      prev_job.rv() if prev_job else seed_gfa_id,
                                      prev_gfa_path if prev_job else seed_gfa_path,
                                      pan_sn_output, graph_names,
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

def minigraph_construct(job, options, config_node, seq_id_map, seq_order, gfa_path, prev_gfa_id, prev_gfa_path, pan_sn_output,
                        graph_names=None):
    """ Make minigraph.

    graph_names is every genome the finished graph can contain, which is only seq_id_map's keys
    when the graph is built from scratch: extending one leaves the genomes already in it out of
    seq_id_map, but their SN tags still have to be renamed on the way out """

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
        minigraph_gfa_to_pansn(graph_names if graph_names else set(seq_id_map.keys()), gfa_path, pansn_gfa_path, job.cores)
        pansn_gfa_out_id = job.fileStore.writeGlobalFile(pansn_gfa_path)
        return gfa_out_id, pansn_gfa_out_id
    else:
        return gfa_out_id

def open_gfa_for_rename(gfa_path, out_gfa_path):
    """ open a GFA and its renamed copy.  the copy is left uncompressed here and bgzipped by
    bgzip_gfa_rename() below: python's gzip module is single threaded and defaults to the slowest
    compression level, which otherwise takes well over 90% of the runtime of these renaming passes """
    if gfa_path.endswith('.gz'):
        in_file = gzip.open(gfa_path, 'rb')
    else:
        in_file = open(gfa_path, 'rb')
    raw_out_path = out_gfa_path[:-3] if out_gfa_path.endswith('.gz') else out_gfa_path
    return in_file, open(raw_out_path, 'wb'), raw_out_path

def bgzip_gfa_rename(raw_out_path, out_gfa_path, threads):
    """ compress what open_gfa_for_rename() left uncompressed, if the output wanted compressing """
    if raw_out_path != out_gfa_path:
        cactus_call(parameters=['bgzip', '--threads', str(threads)], infile=raw_out_path, outfile=out_gfa_path)
        os.remove(raw_out_path)

def minigraph_gfa_to_pansn(names, gfa_path, out_gfa_path, threads=1):
    """ hack to convert cactus names like id=simChimp.0|simChimp.chr6 to PanSN simChimp#0#simpChimp.chr6
    so that minigraph GFA file can be used outside of catus

    todo: Cactus should probably be changed to just use PanSN internally as well, but that's a much
    bigger lift
    """
    in_file, out_file, raw_out_path = open_gfa_for_rename(gfa_path, out_gfa_path)

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
                    # one splitter for every artifact: this GFA, the GAF that indexes into it,
                    # and hal2vg's final graph.  splitting here on the last '.' verbatim instead
                    # gave HG002.01 -> HG002#01# against the GAF's HG002#1#, and HG002.pat ->
                    # HG002#pat#, which is not a PanSN haplotype at all
                    toks[4+i] = 'SN:Z:{}#{}'.format(event_to_pansn_prefix(name), tok[barpos+1:])
                    break
            out_file.write(('\t'.join(toks) + '\n').encode())
        else:
            out_file.write(line.encode())

    in_file.close()
    out_file.close()
    bgzip_gfa_rename(raw_out_path, out_gfa_path, threads)

def minigraph_gfa_from_pansn(job, names, gfa_path, gfa_id):
    """ hack to convert PanSN names like simChimp#0#simpChimp.chr6 to Cactus names like id=simChimp.0|simChimp.chr6
    so that a minigrpah GFA (as converted panSN by minigraph_gfa_to_pansn() above) can be read back into Cactus

    returns (converted gfa id, set of genomes the graph was built from).  the genome set is what
    --inGFA needs to work out which of the seqfile's genomes are new to the graph, and it comes
    free with the pass that has to read every SN tag anyway.

    a GFA that is already in cactus naming -- from a cactus old enough to have published one, or
    handed straight from one stage to the next -- is read for its genome set and returned as-is.

    todo: Cactus should probably be changed to just use PanSN internally as well, but that's a much
    bigger lift
    """
    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, os.path.basename(gfa_path))
    job.fileStore.readGlobalFile(gfa_id, gfa_path)
    out_gfa_path = os.path.join(work_dir, 'cactus.' + os.path.basename(gfa_path))

    in_file, out_file, raw_out_path = open_gfa_for_rename(gfa_path, out_gfa_path)

    events = set()
    unresolved = set()
    already_cactus = False
    for line in in_file:
        line = line.decode()
        if line.startswith('S'):
            toks = line.strip().split('\t')
            for i, tok in enumerate(toks[4:]):
                if tok.startswith('SN:Z:'):
                    if tok.startswith('SN:Z:id='):
                        # already cactus-named: nothing to rewrite, just harvest the genome
                        already_cactus = True
                        barpos = tok.find('|')
                        if barpos > 8:
                            events.add(tok[8:barpos])
                        break
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
                        if '{}.{}'.format(name, hap) in names:
                            name = '{}.{}'.format(name, hap)
                        else:
                            # collected rather than asserted on: with --inGFA this is usually
                            # the user leaving a genome out of the seqfile, which deserves to be
                            # named.  the whole tag goes in the message because the other way to
                            # land here is an SN tag that is not SAMPLE#HAP#CONTIG at all, and
                            # then the prefix alone says nothing about what is wrong
                            unresolved.add(tok)
                            break
                    events.add(name)
                    toks[4+i] = 'SN:Z:id={}|{}'.format(name, tok[hashpos2+1:])
                    break
            out_file.write(('\t'.join(toks) + '\n').encode())
        else:
            out_file.write(line.encode())

    in_file.close()
    out_file.close()

    if unresolved:
        raise RuntimeError('{} sequence name(s) in {} could not be matched to a seqfile genome: {}. Genomes cannot '
                           'be removed from a minigraph, so every genome in the graph must be in the seqfile -- '
                           'unless these names are not SAMPLE#HAP#CONTIG, in which case the graph was not written '
                           'by cactus and cannot be read back into it'.format(
                               len(unresolved), gfa_path, ' '.join(sorted(unresolved)[:10])))

    if already_cactus:
        return gfa_id, events

    bgzip_gfa_rename(raw_out_path, out_gfa_path, job.cores)

    return job.fileStore.writeGlobalFile(out_gfa_path), events

    
