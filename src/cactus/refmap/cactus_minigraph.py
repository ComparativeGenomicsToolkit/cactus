#!/usr/bin/env python3

"""
build a minigraph in Toil, using a cactus seqfile as input
"""

import os, sys, re
import signal
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
from cactus.shared.common import setupBinaries, importSingularityImage, cactus_walltime
from cactus.refmap.pangenome_exclusions import event_to_pansn_prefix
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles, write_s3
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options, add_cactus_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import clean_jobstore_files, unzip_gz
from cactus.shared.version import cactus_commit
from cactus.progressive.cactus_prepare import human2bytesN
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers
from cactus.paf.last_scoring import last_train, last_train_enabled
from toil.job import Job
from toil.job import PromisedRequirement
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from toil.lib.conversions import bytes2human
from cactus.shared.common import cactus_cpu_count
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.common import GZIP_COMPRESS_BYTES_PER_SEC
from cactus.progressive.multiCactusTree import MultiCactusTree
from sonLib.bioio import getTempDirectory, getTempFile

def main():
    parser = Job.Runner.getDefaultArgumentParser()
    add_cactus_toil_options(parser)

    parser.add_argument("seqFile", help = "Seq file (or chromfile with --batch)")
    parser.add_argument("outputGFA", help = "Output Minigraph GFA (or directory in --batch mode)")
    parser.add_argument("--reference", required=True, nargs='+', type=str,
                        help = "Reference genome name(s) (added to minigraph first). Mash distance to 1st reference to determine order of other genomes (use minigraphSortInput in the config xml to toggle this behavior).")
    parser.add_argument("--mgCores", type=int, help = "Number of cores for minigraph construction (defaults to the same as --maxCores).")
    parser.add_argument("--mgMemory", type=human2bytesN,
                        help="Memory in bytes for the minigraph construction job (defaults to an estimate based on the input data size). "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))", default=None)
    parser.add_argument("--lastTrain", action="store_true",
                        help="Deprecated: last-train is now on by default. Use the lastTrain attribute of <graphmap> "
                        "in the config to turn it off", default=False)
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
            # a zip setting that cannot run should fail here, not after construction
            check_graph_rewrite_config(config_node)
            # from here on options.lastTrain is the config's say, not the deprecated flag's.  a
            # --refOnly graph is only used for splitting, so no alignment will read its model
            options.lastTrain = last_train_enabled(options, config_node) and not options.refOnly

            # apply cpu override
            if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
                if not options.mgCores:
                    options.mgCores = sys.maxsize
                options.mgCores = min(options.mgCores, cactus_cpu_count(), int(options.maxCores) if options.maxCores else sys.maxsize)
                # zipWalks="gaf" maps with graphmap's per-job cores, so they get the clamp cactus-graphmap gives them
                graphmap_node = findRequiredNode(config_node, "graphmap")
                graphmap_node.attrib["cpu"] = str(min(getOptionalAttrib(graphmap_node, "cpu", typeFn=int, default=1),
                                                      options.mgCores))
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
                
            # output_dict:  chrom-> (gfa_id, pansn_gfa_id, input_pansn_gfa_id, rewrite_artifacts, train_id)
            output_dict = toil.start(Job.wrapJobFn(minigraph_construct_batch_workflow, options, config_node, input_dict, options.outputGFA,
                                                   in_gfa_id=in_gfa_id, walltime=cactus_walltime()))

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

def rewrite_artifact_paths(gfa_path):
    """ the side artifacts the zip (rgfa-zip, <graphmap zipAlleles>) of <x>.sv.gfa.gz leaves beside it:
        input_gfa  <x>.sv.unzipped.gfa.gz   the graph as minigraph built it
        report     <x>.sv.zip.tsv           rgfa-zip's report, one row per candidate and skipped site
        walks_gaf  <x>.sv.unzipped.gaf.gz   zipWalks="gaf" only: every genome mapped to the unzipped
                                            graph, which is where rgfa-zip read its walks

    only the file name is rewritten, so a '.gfa' in a directory name is left alone """
    out_dir, name = gfa_path[:gfa_path.rfind('/') + 1], gfa_path[gfa_path.rfind('/') + 1:]
    stem, ext = name, ''
    for suffix in ('.gfa.gz', '.gfa'):
        if name.endswith(suffix):
            stem, ext = name[:-len(suffix)], suffix
            break
    return {'input_gfa': out_dir + stem + '.unzipped' + ext,
            'report': out_dir + stem + '.zip.tsv',
            'walks_gaf': out_dir + stem + '.unzipped.gaf.gz'}

def export_rewrite_artifacts(exporter, gfa_path, input_pansn_gfa_id, rewrite_artifacts):
    """ write what the zip leaves beside the graph it rewrote: the graph as minigraph built it, the
    report, and with zipWalks="gaf" the mappings the zip read its walks from.  rewrite_artifacts is
    the fourth slot minigraph_construct_run returns: None when nothing was zipped, else a dict of the
    zip's files (its 'report', and its 'walks_gaf' or None).  `exporter` is anything with
    exportFile: a Toil object at the top level, or job.fileStore inside a job """
    if not rewrite_artifacts:
        return
    paths = rewrite_artifact_paths(gfa_path)
    if input_pansn_gfa_id:
        exporter.exportFile(input_pansn_gfa_id, makeURL(paths['input_gfa']))
    for key in ('report', 'walks_gaf'):
        if rewrite_artifacts.get(key):
            exporter.exportFile(rewrite_artifacts[key], makeURL(paths[key]))

def export_minigraph_construct_output(options, input_seqfiles, output_dict, toil):
    if options.batch:
        chrom_file_path = os.path.join(options.outputGFA, 'chromfile.mg.txt')
        if chrom_file_path.startswith('s3://'):
            chrom_file_temp_path = getTempFile()
        else:
            chrom_file_temp_path = chrom_file_path                    
        chromfile = open(chrom_file_temp_path, 'w')
    for chrom, output_ids in output_dict.items():
        gfa_id, pansn_gfa_id, input_pansn_gfa_id, rewrite_artifacts, train_id = output_ids
        if options.batch:
            gfa_path = os.path.join(options.outputGFA, chrom + '.sv.gfa.gz')
        else:
            gfa_path = options.outputGFA
        if train_id:
            train_path = gfa_path.replace('.gfa.gz', '.gfa').replace('.gfa', '.train')
        else:
            train_path = None
        #export the gfa
        # with mgSplitWholeGenomeRef this graph is still whole-genome: it only becomes
        # chromosome-only after the post-mapping prune, which writes it to this same path
        # (export_pruned_minigraph_gfa_wrapper in cactus_pangenome.py).  the chromfile below is
        # still written now -- downstream only ever reads its .train column
        if not getattr(options, 'mgSplitWholeGenomeRef', False):
            # pansn_gfa_id is already the zipped graph when zipAlleles is on, so the main export is
            # unconditional
            toil.exportFile(pansn_gfa_id, makeURL(gfa_path))
        # the pre-zip graph, the report and (zipWalks="gaf") the walks GAF are extra artifacts kept
        # beside the graph for comparison.  They go out even when the main export above is deferred
        # to the prune: the prune never touches them, and keeping such artifacts inside that branch
        # once left a 460-haplotype mgSplitWholeGenomeRef run with 25 rewritten chromosomes and not
        # one report on disk.  In that mode the unzipped graph is the unpruned, whole-reference one,
        # and the walks GAF is against it.
        export_rewrite_artifacts(toil, gfa_path, input_pansn_gfa_id, rewrite_artifacts)
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
                                       construct_ref_id_map=None, in_gfa_id=None, scores_id=None):
    """ run the construction workflow on individual chromosomes.  construct_ref_id_map, if given,
    swaps the whole-genome reference fastas in for the chromosome's own slice of them (--mgSplit
    mgSplitWholeGenomeRef).  the merge happens here rather than inside minigraph_construct_workflow
    because both dicts are resolved job arguments at this point, whereas the sanitized map down there
    can still be an unresolved promise.  scores_id is a --scoresFile model, which the mapping that
    zipWalks="gaf" runs before the zip uses as graphmap will """
    # everything hangs off a child, so that borrow_last_train_models, a follow-on, is done before
    # any follow-on of this job starts.  As a follow-on of this job it ran alongside cactus-pangenome's
    # export_minigraph_batch_wrapper, the follow-on that reads its result, which then failed on the
    # unresolved promise whenever it started first ("This job was passed promise ... that wasn't yet
    # resolved"), and only got through by being retried
    root_job = Job(walltime=cactus_walltime())
    job.addChild(root_job)
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
        mgwf_job = root_job.addChildJobFn(minigraph_construct_workflow, options, config_node, seq_id_map, seq_order, gfa_path,
                                          sanitize, construct_seq_id_map=construct_seq_id_map, in_gfa_id=in_gfa_id,
                                          scores_id=scores_id, walltime=cactus_walltime())
        output_dict[chrom] = mgwf_job.rv()
    if options.lastTrain and len(input_dict) > 1:
        # sized off the chromosome's own reference slice, which is what it trains on
        ref_sizes = {chrom: input_info[0][options.reference[0]].size for chrom, input_info in input_dict.items()
                     if options.reference[0] in input_info[0]}
        return root_job.addFollowOnJobFn(borrow_last_train_models, output_dict, ref_sizes,
                                         walltime=cactus_walltime()).rv()
    return output_dict

def borrow_last_train_models(job, output_dict, ref_sizes):
    """ a chromosome that couldn't train its own scoring model (too small, no partner, or last-train
    failed) gets another chromosome's instead of falling back to the default scores: a model trained
    on the same genomes is much closer to right than HOXD70 is.  The donor is the chromosome with a
    model whose reference size is the median of those that have one, as neither the smallest (least
    data) nor the largest is typical.  The borrowed model is the same file id, so it is exported
    under each chromosome's own name """
    trained = sorted([chrom for chrom, val in output_dict.items() if val[4]], key=lambda c: (ref_sizes.get(c, 0), c))
    if not trained or len(trained) == len(output_dict):
        return output_dict
    donor = trained[len(trained) // 2]
    borrowers = sorted(chrom for chrom, val in output_dict.items() if not val[4])
    RealtimeLogger.info('Borrowing the scoring model of {} for {} chromosome(s) that could not train their own: {}'.format(
        donor, len(borrowers), ' '.join(borrowers)))
    output_dict = dict(output_dict)
    for chrom in borrowers:
        output_dict[chrom] = tuple(output_dict[chrom][:4]) + (output_dict[donor][4],)
    return output_dict
                                    
def minigraph_construct_workflow(job, options, config_node, seq_id_map, seq_order, gfa_path, sanitize=True,
                                 construct_seq_id_map=None, in_gfa_id=None, scores_id=None):
    """ minigraph can handle bgzipped files but not gzipped; so unzip everything in case before running

    construct_seq_id_map, when given, replaces seq_id_map for the graph construction alone.  it is how
    mgSplitWholeGenomeRef hands each chromosome the whole reference: mash sorting and last-training
    (and the size estimate for last-training's own job) keep using the chromosome's own much smaller
    reference slice, since training against a whole genome would be both expensive and
    cross-chromosome contaminated -- and would silently produce no model at all, as last_train()
    requires its partner sequence to be at least half the size of the database it trains against.
    The construction job itself is sized off the substituted map, since that is what it runs on.

    with a graph to extend, which genomes still need constructing is not known until that graph's
    SN tags have been read, so the rest of the workflow is deferred behind the job that reads them

    scores_id is a --scoresFile model, only read by zipAlleles' zipWalks="gaf" mapping """
    if not in_gfa_id:
        return minigraph_construct_run(job, options, config_node, seq_id_map, seq_order, gfa_path, sanitize,
                                       construct_seq_id_map=construct_seq_id_map, scores_id=scores_id)

    # the renaming pass decompresses the GFA before bgzipping it back up, so it needs room for
    # the raw copy (reckoned at 10x, as elsewhere) on top of the compressed input and output
    rename_job = job.addChildJobFn(minigraph_gfa_from_pansn, set(seq_id_map.keys()), options.inGFA, in_gfa_id,
                                   disk=in_gfa_id.size*12,
                                   walltime=cactus_walltime(GFA_RENAME_SECS_PER_GB * in_gfa_id.size / 1e9,
                                                            io_bytes=RAW_BYTES_PER_GZ_BYTE * in_gfa_id.size))
    run_job = rename_job.addFollowOnJobFn(minigraph_construct_run, options, config_node, seq_id_map, seq_order, gfa_path,
                                          sanitize, construct_seq_id_map,
                                          rename_job.rv(0), rename_job.rv(1), in_gfa_id,
                                          scores_id=scores_id, walltime=cactus_walltime())
    # all five slots: minigraph_construct_run returns
    # (gfa, pansn_gfa, input_pansn_gfa, rewrite_artifacts, train).  Truncating here to three
    # left cactus_pangenome's rv(3)/rv(4) reading off the end of the tuple, which Toil reports as
    # "IndexError: tuple index out of range" from _fulfillPromises, nowhere near the cause.
    return run_job.rv(0), run_job.rv(1), run_job.rv(2), run_job.rv(3), run_job.rv(4)

def minigraph_construct_run(job, options, config_node, seq_id_map, seq_order, gfa_path, sanitize=True,
                            construct_seq_id_map=None, seed_gfa_id=None, seed_events=None, seed_pansn_gfa_id=None,
                            scores_id=None):
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
            # zipping the extended graph splits nodes at alignment-block boundaries, deletes the alt
            # it replaces and rewires the edges around them, so the mappings cactus-pangenome --inGAF
            # reuses no longer tile it.  Left alone this surfaces only after construction, as a
            # gaf2unstable tiling assertion in check_reusable_gaf.  Resuming (no seq_order) zips
            # nothing, so the reused mappings still fit.  getattr: cactus-minigraph has no --inGAF
            if getattr(options, 'inGAF', None) and not getattr(options, 'remap', False) and \
               zip_alleles_enabled(config_node):
                raise RuntimeError(
                    'zipAlleles cannot be used with --inGAF when adding genomes to a graph: the zip rewrites the '
                    'node boundaries of the extended graph, so the mappings in {} no longer resolve against it.  '
                    'Either drop --inGAF (--inGFA still saves the construction, which is the expensive half), add '
                    '--remap to map every genome against the zipped graph, or set zipAlleles="0" in <graphmap> in '
                    'the config'.format(options.inGAF))
            RealtimeLogger.info('Extending the {} genomes in {} with {}: {}'.format(
                len(seed_events), options.inGFA, len(seq_order), ' '.join(seq_order)))
        else:
            # the seqfile asks for exactly the genomes the graph already holds, so there is nothing
            # to construct and the run resumes from it.  the graph handed back is the one --inGFA
            # was given, but its compression follows the *input* name while everything downstream
            # reads the output name to decide whether to unzip, so it is re-emitted to match
            RealtimeLogger.info('Resuming from the {} genomes in {}: nothing left to construct'.format(
                len(seed_events), options.inGFA))
            # one bgzip or gunzip of each of the two seed graphs, so it is the compression rate
            # over a raw GFA ~10x the compressed input
            seed_gfa_size = seed_gfa_id.size if hasattr(seed_gfa_id, 'size') else 0
            match_job = job.addChildJobFn(match_gfa_compression, seed_gfa_id, seed_pansn_gfa_id,
                                          options.inGFA, gfa_path,
                                          disk=12 * seed_gfa_size,
                                          walltime=cactus_walltime(
                                              2 * RAW_BYTES_PER_GZ_BYTE * seed_gfa_size / GZIP_COMPRESS_BYTES_PER_SEC,
                                              io_bytes=4 * seed_gfa_size))
            # last_train reads fastas, not the graph, so resuming is no reason to skip it: without
            # this --lastTrain would quietly fall back to the default scoring matrix
            train_id = None
            if options.lastTrain and len(train_seq_id_map) > 1:
                train_job = job.addChildJobFn(last_train, config_node, train_seq_order, train_seq_id_map,
                                              ref_name=options.reference[0],
                                              cores=options.mgCores, disk=8*ref_size,
                                              memory=cactus_clamp_memory(max(8*ref_size, 12*10**9)),
                                              walltime=cactus_walltime(LAST_TRAIN_SECS + ref_size / 1e6,
                                                                       io_bytes=3 * ref_size))
                train_id = train_job.rv()
            # same five slots as the constructing path below.  nothing was constructed, so nothing
            # is zipped and there is no report: the graph handed back is the seed as it was.  a
            # seed published with zipAlleles on is already zipped
            if zip_alleles_enabled(config_node):
                RealtimeLogger.info('Not running rgfa-zip on {}: resuming constructs nothing, so the graph is passed '
                                    'on as it was given'.format(options.inGFA))
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
        sanitize_job = job.addChildJobFn(sanitize_fasta_headers, seq_id_map, pangenome=True, walltime=cactus_walltime())
        sanitized_seq_id_map = sanitize_job.rv()
    else:
        sanitized_seq_id_map = seq_id_map
        sanitize_job = Job(walltime=cactus_walltime())
        job.addChild(sanitize_job)
    xml_node = findRequiredNode(config_node, "graphmap")
    sort_type = getOptionalAttrib(xml_node, "minigraphSortInput", str, default=None)
    if sort_type == "mash" and len(seq_id_map) > 2:
        sort_job = sanitize_job.addFollowOnJobFn(sort_minigraph_input_with_mash, options, config_node, sanitized_seq_id_map, seq_order,
                                                 ref_name=options.reference[0] if seed_events is not None else None,
                                                 walltime=cactus_walltime())
        seq_order = sort_job.rv()
        prev_job = sort_job
    else:
        prev_job = sanitize_job
    minigraph_job = prev_job.addFollowOnJobFn(minigraph_construct_in_batches, options, config_node,
                                              construct_seq_id_map if construct_seq_id_map else sanitized_seq_id_map,
                                              seq_order, gfa_path,
                                              whole_genome_ref=bool(construct_seq_id_map),
                                              seed_gfa_id=seed_gfa_id, graph_names=graph_names,
                                              walltime=cactus_walltime())

    # last-train is scheduled ahead of the zip below only because zipWalks="gaf" maps with its
    # model, and so has to wait for it.  It still runs alongside construction either way
    train_id, last_train_job = None, None
    if options.lastTrain and len(seq_id_map) > 1:
        # note: somehow last training memory overruns don't seem to be detected by slurm so we
        # give 12G at least whenever possible, as --doubleMem won't help...
        # lastdb dominates the runtime and is erratic: over the 24 per-chromosome runs of the
        # HPRC v2.1 pangenome (2.0e8-byte reference, 8 cores) it took 55 s at the p50 but
        # 3934 s at the worst, while the 16x bigger whole-genome reference of HPRC v2.0
        # (32 cores) took 2821 s.  last-train itself adds 409-680 s on top.  So the estimate is
        # mostly a flat allowance for that tail, with a small linear term so that small inputs
        # aren't over-provisioned.  The I/O is the reference plus the training partner, which
        # last_train() only picks inside the job but constrains to at least half the reference.
        last_train_job = prev_job.addFollowOnJobFn(last_train, config_node, seq_order, sanitized_seq_id_map,
                                                   ref_name=options.reference[0] if seed_events is not None else None,
                                                   cores=options.mgCores,
                                                   disk=8*ref_size,
                                                   memory=cactus_clamp_memory(max(8*ref_size, 12*10**9)),
                                                   walltime=cactus_walltime(LAST_TRAIN_SECS + ref_size / 1e6,
                                                                            io_bytes=3 * ref_size))
        train_id = last_train_job.rv()

    # optionally zip the alleles minigraph stored as novel sequence where they duplicate the
    # reference (or each other) with rgfa-zip (zipAlleles).  when it runs, the zipped graph
    # REPLACES the constructed one in both namings: it is what graphmap, rgfa-split and the join
    # all have to see, since a graph they map against that still carries the duplicated alt would
    # defeat the point.  the originals are returned alongside so the export can keep them for
    # comparison.
    #
    # rgfa-zip needs the PanSN graph, not the cactus-named one: it locates sites with
    # `vg snarls -n -P <reference sample>`, and cactus names (id=EVENT|CONTIG) carry no PanSN sample
    # for -P to match -- with the wrong -P, vg skips every snarl as a non-reference boundary and the
    # tool silently finds nothing.  So zip the PanSN graph and rename the result back.
    #
    # skipped on the --mgSplit reference-only first pass: that graph has no non-reference nodes, so
    # there is nothing to zip and the snarl decomposition would be pure cost.  it still runs on
    # the per-chromosome all-sample graphs, which is where the alleles are.  With
    # mgSplitWholeGenomeRef those graphs hold the whole reference until the prune after graphmap,
    # so that is the graph zipped: graphmap maps to it, and the prune then cuts the zipped graph
    # back to the chromosome.
    input_pansn_gfa_id, rewrite_artifacts = None, None
    if not getattr(options, 'refOnly', False) and zip_alleles_enabled(config_node):
        # Size on the reference the graph was actually built from.  ref_size above is the
        # chromosome's slice, but with mgSplitWholeGenomeRef construct_seq_id_map swaps in the
        # whole-genome reference, so every per-chromosome graph carries all of it.  On a 460-
        # haplotype human run 8x the slice came to 2 GiB, and vg snarls alone peaks at 2.0-2.4
        # GiB: every chromosome OOM'd (exit 137 in vg snarls) and passed only after Toil had
        # doubled it to 8 GiB, or 16 GiB for chr6, chr17 and chrX.
        zip_ref_size = (construct_seq_id_map or seq_id_map)[options.reference[0]].size
        gaf_seq_id_maps, zip_scores_id, wait_for_training = None, None, False
        if zip_walks_mode(config_node) == 'gaf':
            # every genome in the graph is mapped, the way graphmap will map them: rgfa-zip trusts
            # these walks to be all the haplotypes there are.  Extending a graph adds genomes to it
            # without the seed's own genomes passing through construction, so those are added back
            # (sanitizing them first if this run sanitizes, which construction never needed to)
            gaf_seq_id_maps = [sanitized_seq_id_map]
            seed_only = {name: fa_id for name, fa_id in train_seq_id_map.items()
                         if name not in seq_id_map} if seed_events is not None else {}
            if seed_only and sanitize:
                # a child of construction, so it is done before construction's follow-on, the zip
                seed_sanitize_job = minigraph_job.addChildJobFn(sanitize_fasta_headers, seed_only, pangenome=True,
                                                                walltime=cactus_walltime())
                gaf_seq_id_maps.append(seed_sanitize_job.rv())
            elif seed_only:
                gaf_seq_id_maps.append(seed_only)
            # the model graphmap will map with: --scoresFile's, or the one trained on this graph's
            # genomes, which means waiting for last-train.  A --mgSplit chromosome that trains no
            # model of its own borrows one only once every chromosome is built, which is after its
            # zip, so its walks are mapped with minigraph's default penalties.  They only change the
            # base-level alignment, not the chaining that picks the walk
            if getOptionalAttrib(xml_node, "lastTrainMap", typeFn=bool, default=False):
                if scores_id:
                    zip_scores_id = scores_id
                elif last_train_job:
                    zip_scores_id = train_id
                    wait_for_training = True
        zip_job = minigraph_job.addFollowOnJobFn(zip_alleles_workflow, options, config_node,
                                                 minigraph_job.rv(0), minigraph_job.rv(1), gfa_path, graph_names,
                                                 zip_ref_size, gaf_seq_id_maps=gaf_seq_id_maps,
                                                 scores_id=zip_scores_id, walltime=cactus_walltime())
        if wait_for_training:
            # a second predecessor: the zip starts once both construction and training are done
            last_train_job.addFollowOn(zip_job)
        input_pansn_gfa_id = minigraph_job.rv(1)
        rewrite_artifacts = zip_job.rv(2)
        gfa_ids = (zip_job.rv(0), zip_job.rv(1))
    else:
        gfa_ids = (minigraph_job.rv(0), minigraph_job.rv(1))

    # (cactus-named graph, PanSN graph, the PanSN graph before zipping or None, the zip's artifacts
    #  or None, LAST scoring model or None).  slots 0 and 1 are the graph the pipeline uses; slot 3
    #  is the dict export_rewrite_artifacts() reads
    return gfa_ids[0], gfa_ids[1], input_pansn_gfa_id, rewrite_artifacts, train_id

# Fixed seconds for a last_train job.  last-train's cost is set by how divergent the pair it picks
# is, not by how big the reference is, and the two are close to anti-correlated here: over the 27
# per-chromosome buckets of two HPRC v2.1 runs the median job took 130 s and the largest reference
# (chr1) finished in 188 s, while the *second smallest* -- the unplaced/unlocalized bucket -- ran
# 9807 s in one run and 6931 s in the other.  last_train() picks the furthest genome in the mash
# order that clears its size floor, and that bucket's mash distances span the whole range, so its
# pick sits at 0.106 against 0.0007-0.0037 for a real chromosome.
#
# ref_size is therefore the wrong regressor and this flat term is what carries the pathology.  The
# old 3000 asked 2:05:30 for that job; it ran 2:43:26 and survived only because the partition it
# was in did not enforce the limit (its two retries doubled memory, correctly leaving the walltime
# alone -- both failures were MEMLIMIT).  6000 asks 4:10:32, which clears the worse of the two runs
# by 1.5x, enough for the 41% they differed by.  It costs the other 26 jobs nothing that matters:
# they already ask over two hours, so none of them changes partition.
LAST_TRAIN_SECS = 6000

# Zipping alleles (<graphmap zipAlleles>): rgfa-zip aligns every allele minigraph stored as novel
# sequence to the reference window it bypasses, on both strands, and (in its alt-vs-alt pass) to a
# parallel allele that shares its bounding handles, and where the homology is confident (5 kb at
# 95% gap-compressed identity, a unique placement, an unfragmented chain) replaces the allele by
# the sequence it duplicates.  The three versions compared are options of it:
#   v1  the default: zipOptions with --no-alt, the reference pass only
#   v2  zipOptions without --no-alt: the reference pass, then alt-vs-alt
#   v3  zipWalks="gaf" with v2's options: walks observed by mapping every genome to the unzipped graph
# On a CHM13 30-way pangenome the alt-vs-alt pass was a small but consistent loss (raw GT F1 lower
# than v1's for 13 of its 14 samples), hence the default

# v1's zipOptions, as in the shipped config, for a config that turns zipAlleles on without giving
# them: rgfa-zip's own defaults would be v2.  An explicit zipOptions="" still means those
ZIP_DEFAULT_OPTIONS = '-b 5000 -i 0.95 -G 50 -x asm20 --no-alt'

# rgfa-zip's walk sources.  creator (the default) and witness are rebuilt from the graph alone; gaf
# reads them from a GAF that zip_alleles_workflow makes by mapping every genome to the graph
ZIP_WALKS = ('creator', 'witness', 'gaf')

# zipOptions that cactus owns: the output files, and the walk source, which is zipWalks
ZIP_RESERVED_OPTIONS = ('-o', '-r', '--walks')

# rgfa-zip's exit codes other than 0, as its CLI contract defines them.  Nothing is written on any
# of them, so each fails the job
ZIP_EXIT_CODES = {2: 'invalid input: an option it does not take (see zipOptions), or the graph, its snarls or the walks '
                     'failed its input checks',
                  3: 'invariant failure: an edit failed validation and could not be undone',
                  4: 'systemic aligner failure: minimap2 failed on more than 5% of windows, or was killed twice '
                     'on one (an out-of-memory kill is a signal too)',
                  5: 'I/O error writing the output'}

# How zip_alleles runs rgfa-zip: bash -c ZIP_LOG_TEE rgfa-zip <the rgfa-zip command>.  Its log still
# reaches the job log line by line as it is written (it is the job's stderr), and tee keeps a copy
# in zip.log, which the summary reads rgfa-zip's totals from.  Its stdout is passed through as it
# was, and pipefail hands back its exit status, or 128 + the signal if one killed it
ZIP_LOG_TEE = 'exec 3>&1; set -o pipefail; "$@" 2>&1 >&3 3>&- | tee zip.log >&2'

# Memory.  rgfa-zip's own peak (its RUSAGE_SELF, which leaves out the aligners it runs) against the
# uncompressed GFA it reads.  With the options cactus runs it with: 1.08x on a 30-way CHM13
# mgSplitWholeGenomeRef graph, which carries the whole reference (chr1's: 3.15 GiB on 3.13 GB),
# and 3.4x and 4.3x on the GRCh38-464 chr1 and chr20 graphs.  The sweep runs, with --dump keeping
# every alignment record, give 1.1-1.9x on the 30-way chromosome graphs, 1.2-3.4x on most of
# GRCh38-464's, and up to 5.4-6.7x on chr1, chr18, chr20 and chrY.  What pushes it up is the big
# sites rather than the graph (chr1:2.65's alt-vs-alt pass alone holds about 1 GB), so the ratio is
# highest on a small graph with a tangle in it.  6x is the worst of the runs without --dump with 40%
# to spare
ZIP_GFA_MEMORY_FACTOR = 6
# each minimap2 it runs peaked at 1.34 GiB at most over the same runs (30-way chr17, with rgfa-zip
# itself at 0.16 GiB: a child's peak RSS counts its parent's at the spawn, so where rgfa-zip is big
# the minimap2 peaks it logs are really rgfa-zip's), and rgfa-zip lowers -j until j x 1.5e9 bytes
# fits --mem, so that is what each aligner process gets
ZIP_ALIGNER_MEMORY = int(1.5e9)
# vg snarls runs first, on its own, and gets a floor of 8x the reference.  It peaks at 2.3-2.5 GiB
# on a 30-way CHM13 mgSplitWholeGenomeRef chromosome graph, where this floor is 25 GB
ZIP_SNARLS_MEMORY_PER_REF_BYTE = 8
# The job takes --mgCores, up to this many.  Every core adds a site thread and an aligner process,
# and so 1.5 GB to the request, while the work is small: production GRCh38-464 chr1 is 2,042 CPU-s,
# and past a few cores its wall time is set by what does not divide: chr1:2.65's single-threaded
# alt-vs-alt pass (270 s), vg snarls (single-threaded in practice: 457 s on a whole-reference graph
# at -t 4) and reading and writing the graph.  The --mgCores 64 of the GRCh38-464 runs would
# otherwise ask 96 GB for aligners that never all run at once
ZIP_MAX_CORES = 16

# Walltime, from the spec's cost tests: a fixed 900 s, 600 s per GB of reference (vg snarls, reading
# and writing the graph), and the CPU the array tangles cost spread over the cores.  The tangles are
# what dominates: 2-3.5 CPU-hours over all of GRCh38-464, nearly all of it at 23 sites (chr1:2.65 alone
# 1,628 s, chr20:28.95 1,022 s), so the top of that range is allowed for every graph, which is
# exactly right for a whole-genome one and leaves a chromosome room to hold most of the tangles.
# There is no single-alignment tail to allow for: rgfa-zip caps every query and window at 5 Mb and
# screens out the satellite windows minimap2 can spend hours on.  Measured since, the worst
# whole chromosome is GRCh38-464 chr1: 3,056 CPU-s and 17:35 wall at -t 3 -j 3 with --check, 3,296
# CPU-s with zipWalks="gaf" (which aligns more units, up to 8x as many on chr21, but adds little CPU
# where the tangles are), and vg snarls takes 460-515 s on a whole-reference chromosome graph, so
# even at 16 cores this asks several times what the slowest chromosome needs
ZIP_ALLELES_BASE_SECS = 900
ZIP_ALLELES_SECS_PER_GB = 600
ZIP_TANGLE_CPU_SECS = 12600

# rgfa-zip's own totals, from the summary lines that end its log (summarize_zip_log):
#   accepted N chain(s), M bp zipped; R reverted[; K of the chain(s) zip bp of their own, ...]
#                                                            ("would zip ..." under --detect-only)
#   alt-vs-alt: ...; X zipped (Y bp), Z refused, W reverted; V near-parallel    (not with --no-alt)
#   done in T s, peak RSS P GB                               (P is GiB: its kB / 1048576)
ZIP_LOG_TOTALS_RE = re.compile(r'\b(accepted|would zip) (\d+) chain\(s\), (-?\d+) bp(?: zipped)?; (\d+) reverted'
                               r'(?:; (\d+) of the chain\(s\) zip bp of their own)?')
ZIP_LOG_ALT_RE = re.compile(r'\balt-vs-alt: .*; (\d+) zipped \((-?\d+) bp\), (\d+) refused, (\d+) reverted; '
                            r'(\d+) near-parallel')
ZIP_LOG_DONE_RE = re.compile(r'\bdone in ([0-9.]+) s, peak RSS ([0-9.]+) GB')

def zip_alleles_enabled(config_node):
    """ whether <graphmap zipAlleles> asks for rgfa-zip after construction """
    return getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "zipAlleles", typeFn=bool, default=False)

def zip_walks_mode(config_node):
    """ <graphmap zipWalks>: where rgfa-zip gets its haplotype walks.  'creator' (the default) """
    walks = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "zipWalks", str, default="creator")
    walks = walks.strip() if walks else "creator"
    if walks not in ZIP_WALKS:
        raise RuntimeError('<graphmap zipWalks="{}"> is not one of {}'.format(walks, ', '.join(ZIP_WALKS)))
    return walks

def zip_options_list(config_node):
    """ <graphmap zipOptions>, split (v1's when the attribute is missing).  Passed through to rgfa-zip
    as they are, except that the output files and the walk source belong to cactus """
    opts = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "zipOptions", str,
                             default=ZIP_DEFAULT_OPTIONS).split()
    for opt in opts:
        if opt in ZIP_RESERVED_OPTIONS or opt.startswith('--walks='):
            raise RuntimeError('<graphmap zipOptions> cannot contain {}: cactus sets {}, and the walk source with '
                               'zipWalks'.format(opt, ' '.join(ZIP_RESERVED_OPTIONS[:2])))
    if opts and opts[-1] == '-m':
        raise RuntimeError('<graphmap zipOptions> ends in -m, which needs the path of a minimap2 after it')
    return opts

def check_graph_rewrite_config(config_node):
    """ check the zip settings up front, so that a config that cannot run fails before the hours of
    construction that come before the zip.  Returns whether zipAlleles is on """
    zip_on = zip_alleles_enabled(config_node)
    if zip_on:
        zip_walks_mode(config_node)
        zip_options_list(config_node)
    return zip_on

def graph_rewrite_name(gfa_path):
    """ what a graph is called in log prefixes and saved inputs: the chromosome with --mgSplit
    (whose construction paths are <chrom>.gfa.gz), the output name otherwise """
    return os.path.basename(gfa_path).replace('.gz', '').replace('.gfa', '')

def zip_alleles_workflow(job, options, config_node, gfa_id, pansn_gfa_id, gfa_path, graph_names, ref_size,
                         gaf_seq_id_maps=None, scores_id=None):
    """ zip the alleles of a freshly constructed graph, returning (cactus-named zipped graph, PanSN
    zipped graph, artifacts), where artifacts is the dict export_rewrite_artifacts() reads.

    gaf_seq_id_maps, given only with zipWalks="gaf", hold the (sanitized) fasta of every genome in
    the graph.  Each is then mapped to the unzipped graph exactly the way cactus-graphmap maps it to
    the zipped one afterwards -- minigraph_map_all, the same options, the same last-train model
    (scores_id) -- and rgfa-zip reads its walks from the merged GAF.  The graph mapped to is always
    the one this zip is about to rewrite, so the GAF and the graph match in every mode:
      whole-genome           the whole-genome graph and every genome
      --mgSplit              this chromosome's graph as built, before graphmap's prune: the whole
                             reference plus the genomes binned to the chromosome, of which the
                             reference is mapped as its chromosome slice, as graphmap maps it
      --mgSplit, with <graphmap_split wholeGenomeRef="0">
                             this chromosome's graph and the genomes binned to it """
    gaf_id = None
    if gaf_seq_id_maps is not None:
        # cactus_graphmap imports this module
        from cactus.refmap.cactus_graphmap import minigraph_map_all
        seq_id_map = {}
        for id_map in gaf_seq_id_maps:
            seq_id_map.update(id_map)
        graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"), "assemblyName", default="_MINIGRAPH_")
        # with --batch graphmap names its per-genome GAFs and their merge after the chromosome
        map_options = copy.deepcopy(options)
        map_options.mg_chrom_name = graph_rewrite_name(gfa_path)
        RealtimeLogger.info('Mapping {} genomes to the unzipped {} for rgfa-zip\'s walks{}'.format(
            len(seq_id_map), graph_rewrite_name(gfa_path), ' (with the last-train model)' if scores_id else ''))
        if gfa_path.endswith('.gz'):
            # graphmap maps to the decompressed graph and sizes its mapping jobs from it, so do the same
            unzip_job = job.addChildJobFn(unzip_gz, gfa_path, gfa_id, delete_original=False, disk=5*gfa_id.size,
                                          walltime=cactus_walltime(60, io_bytes=(1 + RAW_BYTES_PER_GZ_BYTE) * gfa_id.size))
            map_job = unzip_job.addFollowOnJobFn(minigraph_map_all, map_options, ConfigWrapper(config_node),
                                                 unzip_job.rv(), seq_id_map, graph_event, scores_id=scores_id,
                                                 gaf_only=True, walltime=cactus_walltime())
            map_job.addFollowOnJobFn(clean_jobstore_files, file_ids=[unzip_job.rv()], walltime=cactus_walltime())
        else:
            map_job = job.addChildJobFn(minigraph_map_all, map_options, ConfigWrapper(config_node),
                                        gfa_id, seq_id_map, graph_event, scores_id=scores_id,
                                        gaf_only=True, walltime=cactus_walltime())
        gaf_id = map_job.rv(1)
        schedule_job = map_job.addFollowOnJobFn(zip_alleles_schedule, options, config_node, pansn_gfa_id, gfa_path,
                                                graph_names, ref_size, gaf_id, walltime=cactus_walltime())
    else:
        schedule_job = job.addChildJobFn(zip_alleles_schedule, options, config_node, pansn_gfa_id, gfa_path,
                                         graph_names, ref_size, None, walltime=cactus_walltime())
    return schedule_job.rv(0), schedule_job.rv(1), schedule_job.rv(2)

def zip_alleles_walltime(ref_size, cores):
    """ estimated seconds for one zip_alleles job, from the reference its graph carries """
    return (ZIP_ALLELES_BASE_SECS + ZIP_ALLELES_SECS_PER_GB * ref_size / 1e9 +
            ZIP_TANGLE_CPU_SECS / max(1, int(cores or 1)))

def zip_alleles_cores(mg_cores):
    """ cores for one zip_alleles job: --mgCores, up to ZIP_MAX_CORES """
    return max(1, min(int(mg_cores or 1), ZIP_MAX_CORES))

def zip_alleles_memory(raw_gfa_bytes, ref_size, cores):
    """ memory for one zip_alleles job, before clamping: vg snarls, which runs first and on its own,
    or rgfa-zip holding the uncompressed graph plus an aligner process per core, whichever is more """
    return max(ZIP_SNARLS_MEMORY_PER_REF_BYTE * ref_size,
               ZIP_GFA_MEMORY_FACTOR * raw_gfa_bytes + ZIP_ALIGNER_MEMORY * cores)

def zip_aligner_memory(job_memory, raw_gfa_bytes):
    """ rgfa-zip's --mem: what the job has left for aligner processes once rgfa-zip holds the graph,
    from the size of the graph it was actually given.  The job is sized from the same size, so this
    leaves the aligner per core it asked for unless the request was clamped; rgfa-zip lowers -j to
    fit whatever it gets, never below one aligner """
    return max(int(job_memory) - ZIP_GFA_MEMORY_FACTOR * int(raw_gfa_bytes), ZIP_ALIGNER_MEMORY)

def gfa_raw_bytes(job, gfa_id, gzipped):
    """ how big a GFA in the jobstore is once decompressed, counted by streaming it through zlib.
    Memory is sized off 6x this, and the RAW_BYTES_PER_GZ_BYTE guess the disk requests use is 10x the
    compressed size where minigraph graphs inflate only 4.1-4.6x (GRCh38-464 and 30-way CHM13 graphs;
    yeast 3.8x), which would double the request.  Takes about 1 s per 400 MB it inflates to """
    if not gzipped:
        return gfa_id.size
    raw_bytes = 0
    with job.fileStore.readGlobalFileStream(gfa_id) as gfa_stream:
        # GzipFile reads every gzip member in turn, and a bgzipped file is many of them
        with gzip.GzipFile(fileobj=gfa_stream) as gfa_file:
            while True:
                chunk = gfa_file.read(1 << 24)
                if not chunk:
                    return raw_bytes
                raw_bytes += len(chunk)

def zip_alleles_schedule(job, options, config_node, pansn_gfa_id, gfa_path, graph_names, ref_size, gaf_id=None):
    """ size and schedule rgfa-zip and the rename back to cactus naming.  A job of its own because the
    graph and the walks GAF only have sizes once they exist.  Returns what zip_alleles_workflow does """
    cores = zip_alleles_cores(options.mgCores)
    gfa_bytes = pansn_gfa_id.size
    raw_gfa_bytes = gfa_raw_bytes(job, pansn_gfa_id, gfa_path.endswith('.gz'))
    gaf_bytes = gaf_id.size if gaf_id else 0
    raw_gaf_bytes = RAW_BYTES_PER_GZ_BYTE * gaf_bytes
    memory = cactus_clamp_memory(zip_alleles_memory(raw_gfa_bytes, ref_size, cores))
    # the graph compressed and not, the snarls, the zipped graph compressed and not, and the GAF
    # compressed and not, with a floor of 8x the reference
    disk = max(8 * ref_size, 2 * gfa_bytes + 3 * raw_gfa_bytes + gaf_bytes + raw_gaf_bytes)
    RealtimeLogger.info('Scheduling rgfa-zip on {}: {} cores, {} memory, {} disk, for a {} GFA ({} compressed)'.format(
        graph_rewrite_name(gfa_path), cores, bytes2human(memory), bytes2human(disk), bytes2human(raw_gfa_bytes),
        bytes2human(gfa_bytes)))
    zip_job = job.addChildJobFn(zip_alleles, options, config_node, pansn_gfa_id, gfa_path, gaf_id,
                                cores=cores, memory=memory, disk=disk,
                                walltime=cactus_walltime(zip_alleles_walltime(ref_size, cores),
                                                         io_bytes=2 * gfa_bytes + gaf_bytes))
    # graph_names, the same set the forward rename uses: it has to resolve every SN tag in the
    # finished graph, which on the --inGFA extend path is more genomes than minigraph is given
    rename_job = zip_job.addFollowOnJobFn(minigraph_gfa_from_pansn, graph_names, gfa_path, zip_job.rv(0),
                                          disk=12 * gfa_bytes,
                                          walltime=cactus_walltime(GFA_RENAME_SECS_PER_GB * gfa_bytes / 1e9,
                                                                   io_bytes=RAW_BYTES_PER_GZ_BYTE * gfa_bytes))
    # rv(0): minigraph_gfa_from_pansn returns (gfa id, genome set)
    return rename_job.rv(0), zip_job.rv(0), {'report': zip_job.rv(1), 'walks_gaf': gaf_id}

def zip_minimap2(work_dir):
    """ the minimap2 rgfa-zip aligns with, as an absolute path where the binaries run (here, or in
    the container).  rgfa-zip has no PATH default because placements differ between minimap2
    versions -- on chr19 2.17 made 1,280 calls where 2.30 made 753 -- so this is cactus's own copy,
    the one build-tools/downloadPangenomeTools installs beside cactus_consolidated.  A minimap2 that
    is only on the PATH is used with a warning, and the version is logged either way (rgfa-zip also
    writes it into its report) """
    script = ('b=$(command -v cactus_consolidated || true); '
              'if [ -n "$b" ] && [ -x "$(dirname "$b")/minimap2" ]; then '
              'echo "bundled $(cd "$(dirname "$b")" && pwd)/minimap2"; '
              'else echo "path $(command -v minimap2 || true)"; fi')
    found = cactus_call(parameters=['bash', '-c', script], check_output=True, work_dir=work_dir).strip().split(None, 1)
    if len(found) < 2 or not found[1]:
        raise RuntimeError('rgfa-zip needs minimap2, and there is none beside cactus_consolidated or on the PATH')
    source, path = found
    if source != 'bundled':
        RealtimeLogger.warning('No minimap2 beside cactus_consolidated, so rgfa-zip is aligning with {} from the '
                               'PATH: zip placements depend on the minimap2 version'.format(path))
    version = cactus_call(parameters=[path, '--version'], check_output=True, work_dir=work_dir).strip()
    RealtimeLogger.info('rgfa-zip aligner: {} (minimap2 {})'.format(path, version))
    return path

def summarize_zip_report(report_path):
    """ (rows zipped, rows reverted) in an rgfa-zip report.  Read off the outcome values themselves
    rather than a column number, so a column added to the report cannot break it.  These count
    chains, one row each, not what the zip did to the graph: chains that agree on an edit each get a
    row for it, which is why the bp are taken from rgfa-zip's log instead (summarize_zip_log) """
    zipped, reverted = 0, 0
    with open(report_path) as report_file:
        for line in report_file:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if 'zipped' in fields:
                zipped += 1
            if any(field.startswith('reverted:') for field in fields):
                reverted += 1
    return zipped, reverted

def summarize_zip_log(log_path):
    """ rgfa-zip's own totals, from the summary lines that end its log (ZIP_LOG_*_RE), as a dict
    holding whichever were found: 'chains', 'bp' and 'reverted' over both passes ('detect_only' if
    it only said what it would zip, and 'own_chains', the chains that zip bp no earlier chain did,
    where rgfa-zip says), 'alt_chains', 'alt_bp' and 'alt_reverted' for the alt-vs-alt pass, and
    'secs' and 'peak_rss' (bytes) for the whole run.  bp is what the zip took out of the graph.
    The last of each line wins """
    totals = {}
    with open(log_path, errors='replace') as log_file:
        for line in log_file:
            match = ZIP_LOG_TOTALS_RE.search(line)
            if match:
                totals.update(detect_only=match.group(1) == 'would zip', chains=int(match.group(2)),
                              bp=int(match.group(3)), reverted=int(match.group(4)))
                if match.group(5) is not None:
                    totals['own_chains'] = int(match.group(5))
                continue
            match = ZIP_LOG_ALT_RE.search(line)
            if match:
                totals.update(alt_chains=int(match.group(1)), alt_bp=int(match.group(2)),
                              alt_reverted=int(match.group(4)))
                continue
            match = ZIP_LOG_DONE_RE.search(line)
            if match:
                totals.update(secs=float(match.group(1)), peak_rss=int(float(match.group(2)) * 2**30))
    return totals

def zip_failure_meaning(error_text):
    """ what a failed rgfa-zip call means, from cactus_call's error: rgfa-zip's own exit status
    (ZIP_EXIT_CODES), or the signal that killed it, which bash (ZIP_LOG_TEE) reports as 128 + the
    signal number """
    exited = re.search(r'\bexited (\d+)', error_text)
    if not exited:
        return 'killed by a signal' if re.search(r'\bsignaled \w+', error_text) else 'failed'
    status = int(exited.group(1))
    if status in ZIP_EXIT_CODES:
        return ZIP_EXIT_CODES[status]
    if status in (126, 127):
        return 'rgfa-zip could not be run (exit status {}: not found, or not executable)'.format(status)
    if status > 128:
        try:
            sig_name = signal.Signals(status - 128).name
        except ValueError:
            sig_name = 'signal {}'.format(status - 128)
        return 'killed by {}{}'.format(sig_name, ' (out of memory, most likely)' if status - 128 == signal.SIGKILL else '')
    return 'unexpected exit status {}'.format(status)

def zip_summary_line(name, totals, report_zipped, report_reverted, raw_gfa_bytes):
    """ the line the job log gets for one zip: chains and bp from rgfa-zip's own totals, split into
    the reference and alt-vs-alt passes, and its peak memory as a multiple of the GFA it read.  Falls
    back to the report's row counts when the log has no totals line """
    if 'chains' not in totals:
        return ('rgfa-zip on {}: {} zipped and {} reverted row(s) in its report (no totals line in its log, so '
                'no bp)'.format(name, report_zipped, report_reverted))
    line = 'rgfa-zip on {}: {} {} chain(s), {} bp'.format(
        name, 'would zip' if totals['detect_only'] else 'zipped', totals['chains'], totals['bp'])
    if 'alt_chains' in totals:
        line += ' ({} chain(s), {} bp onto the reference, {} chain(s), {} bp alt-vs-alt)'.format(
            totals['chains'] - totals['alt_chains'], totals['bp'] - totals['alt_bp'],
            totals['alt_chains'], totals['alt_bp'])
    else:
        line += ' (all onto the reference: no alt-vs-alt pass)'
    line += '; {} reverted'.format(totals['reverted'])
    if 'own_chains' in totals and totals['own_chains'] != totals['chains']:
        line += '; {} chain(s) only agree with an earlier one\'s edit'.format(totals['chains'] - totals['own_chains'])
    if 'peak_rss' in totals:
        line += '; peak RSS {} in {:.0f} s'.format(bytes2human(totals['peak_rss']), totals['secs'])
        if raw_gfa_bytes:
            # the job allows it ZIP_GFA_MEMORY_FACTOR x; zip_alleles warns when that mattered
            line += ', {:.2f}x its {} GFA'.format(totals['peak_rss'] / raw_gfa_bytes, bytes2human(raw_gfa_bytes))
    return line

def zip_alleles(job, options, config_node, pansn_gfa_id, gfa_path, gaf_id=None):
    """ run rgfa-zip on a PanSN minigraph GFA (see the notes above ZIP_WALKS).  Its sites are vg's
    top-level snarls with reference boundaries, hence the vg calls.  gaf_id is the walks GAF
    (zipWalks="gaf").

    Returns (zipped gfa id, report id).  Any non-zero exit fails the job, after saving the inputs
    and rgfa-zip's log under <outDir>/zip-failed/ so the failure can be reproduced; chains that
    rgfa-zip reverts are counted and logged as warnings.  The job log gets rgfa-zip's log as it runs
    and then one summary line, with the chains and bp zipped (zip_summary_line) """
    work_dir = job.fileStore.getLocalTempDir()
    gzipped = gfa_path.endswith('.gz')
    name = graph_rewrite_name(gfa_path)
    in_gfa = os.path.join(work_dir, 'in.gfa')
    job.fileStore.readGlobalFile(pansn_gfa_id, in_gfa + ('.gz' if gzipped else ''))
    if gzipped:
        cactus_call(parameters=['bgzip', '-d', '--threads', str(job.cores), in_gfa + '.gz'], work_dir=work_dir)

    # -P orients the snarl tree along the reference.  A minigraph rGFA imports with only the
    # reference as a path, so this is a no-op there, but it is correct for graphs that carry more
    snarls = os.path.join(work_dir, 'snarls.json')
    cactus_call(parameters=[['vg', 'snarls', '-n', '-P', options.reference[0], '-t', str(job.cores), 'in.gfa'],
                            ['vg', 'view', '-Rj', '-']],
                outfile=snarls, work_dir=work_dir)

    walks = zip_walks_mode(config_node)
    if walks == 'gaf':
        assert gaf_id
        # merged and bgzipped by minigraph_map_all; rgfa-zip reads it plain
        job.fileStore.readGlobalFile(gaf_id, os.path.join(work_dir, 'walks.gaf.gz'))
        cactus_call(parameters=['bgzip', '-d', '--threads', str(job.cores), 'walks.gaf.gz'], work_dir=work_dir)
        walks = 'gaf:walks.gaf'

    # zipOptions pass through.  cactus adds the aligner and the job's resources unless they are
    # set there: -t sites in parallel and -j aligner processes are the job's cores, and --mem is
    # what is left of its memory once rgfa-zip holds the graph, which rgfa-zip lowers -j to fit
    opts = zip_options_list(config_node)
    cores = max(1, int(job.cores))
    raw_gfa_bytes = os.path.getsize(in_gfa)
    if '-m' in opts:
        RealtimeLogger.info('rgfa-zip aligner for {}: {}, from zipOptions'.format(name, opts[opts.index('-m') + 1]))
    else:
        opts = ['-m', zip_minimap2(work_dir)] + opts
    if '-t' not in opts:
        opts += ['-t', str(cores)]
    if '-j' not in opts:
        opts += ['-j', str(cores)]
    if '--mem' not in opts:
        aligner_memory = zip_aligner_memory(job.memory, raw_gfa_bytes)
        opts += ['--mem', str(aligner_memory)]
        RealtimeLogger.info('rgfa-zip on {}: {} of job memory, less {}x the {} GFA, leaves --mem {}, room for {} aligner '
                            'process(es)'.format(name, bytes2human(job.memory), ZIP_GFA_MEMORY_FACTOR,
                                                 bytes2human(raw_gfa_bytes), aligner_memory,
                                                 aligner_memory // ZIP_ALIGNER_MEMORY))
    if '--tmpdir' not in opts:
        os.makedirs(os.path.join(work_dir, 'zip-tmp'))
        opts += ['--tmpdir', 'zip-tmp']
    out_gfa = os.path.join(work_dir, 'zipped.gfa')
    report = os.path.join(work_dir, 'zip.tsv')
    log_path = os.path.join(work_dir, 'zip.log')
    cmd = ['rgfa-zip'] + opts + ['--walks', walks, 'in.gfa', 'snarls.json', '-o', 'zipped.gfa', '-r', 'zip.tsv']
    try:
        cactus_call(parameters=['bash', '-c', ZIP_LOG_TEE, 'rgfa-zip'] + cmd, work_dir=work_dir,
                    realtimeStderrPrefix='[rgfa-zip-{}]'.format(name), job_memory=job.memory)
    except RuntimeError as e:
        # Keep the exact inputs so the failure can be reproduced, and rgfa-zip's log, which the job
        # log has interleaved with every other job's: by the time anyone looks, the job's temp dir
        # is gone.  Under --mgSplit options.outputGFA is '' and gfa_path is a bare filename, so
        # anchor on --outDir when the pipeline has one
        meaning = zip_failure_meaning(str(e))
        out_root = getattr(options, 'outDir', None) or os.path.dirname(gfa_path) or '.'
        debug_dir = os.path.join(out_root, 'zip-failed')
        if '://' not in debug_dir:
            os.makedirs(debug_dir, exist_ok=True)
        saved = []
        cmd_path = os.path.join(work_dir, 'zip-input.cmd')
        with open(cmd_path, 'w') as cmd_file:
            cmd_file.write('# in.gfa is {0}.zip-input.gfa{1} and snarls.json {0}.zip-input.snarls.json{2}, all '
                           'decompressed\n'.format(name, '.gz' if gzipped else '',
                                                   ', walks.gaf {}.zip-input.gaf.gz'.format(name) if gaf_id else ''))
            cmd_file.write(' '.join(cmd) + '\n')
        for file_id, local, dest_name in ((pansn_gfa_id, None, name + '.zip-input.gfa' + ('.gz' if gzipped else '')),
                                          (None, snarls, name + '.zip-input.snarls.json'),
                                          (gaf_id, None, name + '.zip-input.gaf.gz'),
                                          (None, cmd_path, name + '.zip-input.cmd'),
                                          (None, log_path if os.path.isfile(log_path) else None, name + '.zip.log')):
            if file_id is None and local is None:
                continue
            dest = os.path.join(debug_dir, dest_name)
            job.fileStore.exportFile(file_id if file_id else job.fileStore.writeGlobalFile(local), makeURL(dest))
            saved.append(dest)
        raise RuntimeError('rgfa-zip failed on {} ({}).\nInputs and log saved for reproduction: {}\n{}'.format(
            name, meaning, ' '.join(saved), e))
    missing = [os.path.basename(path) for path in (out_gfa, report) if not os.path.isfile(path)]
    if missing:
        raise RuntimeError('rgfa-zip exited 0 on {} without writing {}'.format(name, ' and '.join(missing)))

    # chains and bp from rgfa-zip's own totals: the report has a row per chain, but no bp total
    totals = summarize_zip_log(log_path) if os.path.isfile(log_path) else {}
    report_zipped, report_reverted = summarize_zip_report(report)
    RealtimeLogger.info(zip_summary_line(name, totals, report_zipped, report_reverted, raw_gfa_bytes))
    reverted = totals.get('reverted', report_reverted)
    if reverted:
        # each one cost only its own chain, which stayed alt: a warning, not a failure
        job.fileStore.logToMaster('WARNING: rgfa-zip reverted {} chain(s) on {} after a validator refused them; '
                                  'they stay unzipped.  See the reverted:<check> rows of its report'.format(reverted, name))
    if totals.get('peak_rss', 0) > max(ZIP_GFA_MEMORY_FACTOR * raw_gfa_bytes, ZIP_ALIGNER_MEMORY):
        # nothing failed, but the next graph like this one could: the job is sized on that factor.
        # Below an aligner's 1.5 GB it is only rgfa-zip's fixed cost showing on a small graph
        job.fileStore.logToMaster('WARNING: rgfa-zip on {} peaked at {}, {:.1f}x its {} GFA, above the {}x cactus '
                                  'sizes its job for'.format(name, bytes2human(totals['peak_rss']),
                                                             totals['peak_rss'] / raw_gfa_bytes,
                                                             bytes2human(raw_gfa_bytes), ZIP_GFA_MEMORY_FACTOR))

    if gzipped:
        cactus_call(parameters=['bgzip', '--threads', str(job.cores)], infile=out_gfa, outfile=out_gfa + '.gz')
        out_gfa += '.gz'
    return job.fileStore.writeGlobalFile(out_gfa), job.fileStore.writeGlobalFile(report)

# Bytes of reference fasta `mash sketch` gets through per second.  It took 225 s on the whole
# 3.15e9-byte CHM13 of HPRC v2.0 (1.4e7 B/s) and a p50 of 3.4 s on the ~1.3e8-byte chromosome
# references of HPRC v2.1 (3.9e7 B/s); this is the slow end of that.
MASH_SKETCH_BYTES_PER_SEC = 1e7

# Bytes of query fasta one mash_dist job gets through per second.  `mash dist` alone runs at
# 7e6-1.9e7 B/s (694 calls over the 232 mash_dist jobs of HPRC v2.0: p50 168 s, p99 429 s for a
# 3.1e9-byte haplotype), but the job also concatenates the sample's haplotypes and counts every
# base of each with Bio.SeqIO.parse, neither of which shows up in the logs as its own command.
# Budgeting those two at ~2e7 and ~2e8 B/s respectively lands the whole job here.
MASH_DIST_BYTES_PER_SEC = 3e6

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
    ref_bytes = seq_id_map[seq_order[0]].size
    sketch_job = job.addChildJobFn(mash_sketch, seq_order[0], seq_id_map,
                                   disk = ref_bytes * 2,
                                   walltime=cactus_walltime(ref_bytes / MASH_SKETCH_BYTES_PER_SEC,
                                                            io_bytes=ref_bytes))
    ref_sketch_id = sketch_job.rv()

    dist_root_job = Job(walltime=cactus_walltime())
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
        sample_bytes = sum(seq_id_map[x].size for x in names)
        dist_map = dist_root_job.addChildJobFn(mash_dist, names, seq_order[0], seq_id_map, ref_sketch_id,
                                               disk = 2 * sample_bytes + ref_bytes,
                                               walltime=cactus_walltime(sample_bytes / MASH_DIST_BYTES_PER_SEC,
                                                                        io_bytes=sample_bytes)).rv()
        dist_maps.append(dist_map)
            
    return dist_root_job.addFollowOnJobFn(mash_distance_order, options, config_node, seq_order, dist_maps, trim_ref,
                                          walltime=cactus_walltime()).rv()

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
            
# Bytes of input fasta `minigraph -xggs` gets through per second per core.  Re-fitted on the
# 207 batch constructs of an HPRC v2.1 run (64 cores, 24 chromosomes, 2.7-13.0 GB per batch),
# which is the first at-scale run of the faster minigraph: dividing the batch bytes by the
# observed seconds gives p50 7.5e5, p10 3.3e5 and a single worst point of 7.9e4 B/s/core.
#
# The old value of 4e4 came from the slower fork (8- and 32-core runs of HPRC v2.0/v2.1), and
# carrying it over asked 29 h for a chr1 batch that now takes 2.3 h and 72 h for the last one --
# every construct into the longest partition, which is the opposite of the point.
#
# The value is set by chrY and nothing else.  Solving both runs for the rate at which each batch's
# ask would exactly equal its observed time, the four tightest points of 307 are all chrY (1.8e5
# to 2.4e5) and the next is chr16 at 6.3e5, against a p50 of 4.1e6 -- so any value that covers
# chrY leaves every other chromosome several times over-provisioned, and there is no way to tell
# them apart from bytes.  chrY is slow because it parallelises badly: it holds a CPU factor of 3.9
# against a p50 of 7.6, and chr2, the next worst, holds 3.2.
#
# 1.5e5 clears chrY's worst observed run by 1.3x, which is the margin that matters because chrY
# swung 14% between the two runs.  The cost of covering it is what pushes the median construct to
# 8 h and 46 of 207 batches past 12 h; 2e5 would leave only 14 past 12 h but clears chrY by 1.0x,
# which is to say it times out.  Nothing reaches even the shorter of the two ceilings these runs
# saw (84 h), and --doubleTime covers a third run worse than either of these.
MINIGRAPH_CONSTRUCT_BYTES_PER_SEC_PER_CORE = 1.5e5

# ...but not linearly.  minigraph parallelises over query contigs, and the faster fork adds more
# parallel sections on top of that, but parts of construction remain single-threaded -- so it is
# Amdahl's law rather than a hard ceiling, and how far it scales depends on the data.  Measured
# from the CPU factor minigraph reports itself (cputime/elapsed, the `*N` in its
# [M::ggen_map::T*N] lines), time-weighted over every construct process of an HPRC run: at
# --mgCores 64 on the current fork it averages 7.96 cores busy, which is a serial fraction of
# 0.112 and an asymptote near nine cores.
#
# Fitted at one core count, so the asymptote is measured and the shape is Amdahl's assumption.
# It matters most in the middle: at --mgCores 8 this gives 4.5 effective cores where crediting
# the request in full would give 8, and simply capping at the asymptote would too.
MINIGRAPH_SERIAL_FRACTION = 0.112

def minigraph_effective_cores(cores):
    """ cores minigraph construction can actually keep busy at this core count, by Amdahl's law """
    cores = max(1, int(cores or 1))
    s = MINIGRAPH_SERIAL_FRACTION
    return 1.0 / (s + (1.0 - s) / cores)

# Fixed cost of a construct batch on top of the alignment itself: staging in the previous
# batch's GFA (or the seed graph when extending), which is a promise here and so cannot be
# sized, and -- on the final batch only -- the in-python PanSN rename plus its bgzip, which is
# under 310 s even for the 830 MB whole-genome GFA of HPRC v2.0.
MINIGRAPH_CONSTRUCT_OVERHEAD_SECS = 600

# Each batch aligns its genomes against everything the batches before it already put in the
# graph, so the same amount of new sequence costs more the later it arrives.  Measured over 23
# chromosomes of an HPRC run, as batch i's wall time against batch 0's: the p90 runs 1.23, 1.58,
# 1.81, 1.92, 2.08 for i = 1..5 and then flattens, so the growth saturates rather than
# compounding.  Keyed off the bytes already in the graph rather than the batch index, so uneven
# batches and the --inGFA seed graph are handled the same way.
#
# The graph itself is the better predictor, but it reaches this point as prev_job.rv(), a promise
# with no size, so the sequence that went into it is the closest thing in scope.
MINIGRAPH_GRAPH_GROWTH = 0.25
MINIGRAPH_GRAPH_GROWTH_MAX = 2.5

def minigraph_graph_growth(prior_bytes, batch_bytes):
    """ how much slower this batch is than the first, for the graph already built ahead of it """
    if batch_bytes <= 0:
        return 1.0
    return min(1.0 + MINIGRAPH_GRAPH_GROWTH * (prior_bytes / batch_bytes), MINIGRAPH_GRAPH_GROWTH_MAX)

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
        # with mgSplitWholeGenomeRef the largest input is the whole reference, so the estimate above
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
    def batch_construct_work(i):
        """ (bytes batch i adds to the graph, seconds minigraph spends adding them) """
        start = i * max_batch_size
        batch_size = len(seq_order) - start if i == num_batches - 1 else max_batch_size
        batch_bytes = sum(seq_id_map[e].size for e in seq_order[start:start + batch_size])
        # everything already in the graph this batch has to align against
        prior_bytes = sum(seq_id_map[e].size for e in seq_order[:start])
        return batch_bytes, (batch_bytes * minigraph_graph_growth(prior_bytes, batch_bytes) /
                             (MINIGRAPH_CONSTRUCT_BYTES_PER_SEC_PER_CORE *
                              minigraph_effective_cores(options.mgCores)))

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
        # batch 0 pays for batch 1 as well.  Toil chains a job into its predecessor's allocation
        # whenever the successor's memory, cores and disk all fit -- walltime is not among the
        # things nextChainable() looks at -- so a chained pair runs under the *first* job's Slurm
        # time limit.  Every batch of a chromosome is issued with identical requirements, and
        # batch 0 is the only one whose sole successor is the next batch: from batch 1 on, each
        # also carries the clean_jobstore_files follow-on below, and two successors end the chain.
        # So batches 0 and 1 always share one allocation, and asking only for batch 0 is what
        # killed chrY's first batch in both of the runs these estimates are drawn from.
        chained = [i] + ([1] if i == 0 and num_batches > 1 else [])
        work = [batch_construct_work(j) for j in chained]
        minigraph_job = Job.wrapJobFn(minigraph_construct, options, config_node, seq_id_map, input_seq_order, out_gfa_path,
                                      prev_job.rv() if prev_job else seed_gfa_id,
                                      prev_gfa_path if prev_job else seed_gfa_path,
                                      pan_sn_output, graph_names,
                                      disk=disk, memory=mem, cores=options.mgCores,
                                      walltime=cactus_walltime(
                                          MINIGRAPH_CONSTRUCT_OVERHEAD_SECS * len(chained) +
                                          sum(secs for _, secs in work),
                                          io_bytes=sum(nbytes for nbytes, _ in work)))
        if prev_job:
            prev_job.addFollowOn(minigraph_job)
            # delete the output of the previous batch from the job store            
            minigraph_job.addFollowOnJobFn(clean_jobstore_files, file_ids=[prev_job.rv()], walltime=cactus_walltime())
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

# A bgzipped GFA or PAF decompresses to about 10x its size, the figure the disk requests here are
# already reckoned at.
RAW_BYTES_PER_GZ_BYTE = 10

# Seconds per GB of *compressed* GFA to rename it between PanSN and Cactus.  Nothing measures this
# one directly -- cactus-pangenome passes pansn_gfa_input=False, so it only runs from the
# standalone cactus-graphmap and cactus-graphmap-split entry points and it fired in none of the
# runs we have logs for.  It rewrites every S-line of the decompressed GFA in python (~40 MB/s) and
# bgzips the result back up (~25 MB/s, the rate of the whole-panel bgzips that were measured), both
# single-threaded and both over a raw GFA ~10x the compressed input it is handed.
GFA_RENAME_SECS_PER_GB = 700

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

    
