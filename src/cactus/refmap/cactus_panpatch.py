#!/usr/bin/env python3

"""
end-to-end assembly patching pipeline, which chains:
cactus-pangenome
panpatch

panpatch (https://github.com/glennhickey/panpatch) uses a pangenome graph to patch a (slightly)
fragmented assembly into T2T chromosomes: it fills N-gaps, scaffolds disconnected contigs and, with
--requireTelomeres, completes missing terminal telomeres, taking the patch sequence from one or more
donor assemblies.

the pangenome is built with --chrom-vg full and nothing else, since the per-chromosome "full" vg
graphs are the only thing panpatch needs.  cactus-pangenome is the overwhelming majority of the
compute here: panpatch itself takes a couple of minutes on a whole human pangenome.  with
--keepGraphs its input graphs are left in <outDir>/<name>.chroms/, so panpatch's parameters can be
re-tuned without rebuilding the pangenome
"""

import os, sys
import copy
import glob
import gzip
import shlex
import shutil
import timeit
import xml.etree.ElementTree as ET

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.cactus_prepare import human2bytesN
from cactus.shared.common import importSingularityImage
from cactus.shared.common import makeURL
from cactus.shared.common import cactus_call
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.common import delete_directory
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import setupBinaries
from cactus.shared.common import enableDumpStack
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.version import cactus_commit
from cactus.refmap.cactus_minigraph import check_sample_names
from cactus.refmap.cactus_graphmap_join import graphmap_join_defaults
from cactus.refmap.cactus_pangenome import pangenome_options, pangenome_validate_options
from cactus.refmap.cactus_pangenome import pangenome_config_overrides, pangenome_end_to_end_workflow
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger, set_logging_from_options
from toil.realtimeLogger import RealtimeLogger

def main():
    parser = Job.Runner.getDefaultArgumentParser()

    parser.add_argument("seqFile", nargs='?', default=None,
                        help = "Seq file (as with cactus-pangenome), or, with --batch, a chromfile with one "
                        "\"<name> <seqfile>\" line per sample")
    parser.add_argument("--outDir", help = "Output directory", required=True)
    parser.add_argument("--outName", help = "Output name (without extension) [default: the name of the sample being patched]")
    parser.add_argument("--reference", nargs='+', type=str,
                        help = "Reference event name(s), as in cactus-pangenome. Must be haploid (no .1/.2 haplotype suffix). "
                        "If omitted, patching is done reference-free: each haplotype of the sample is its own reference, so a "
                        "separate pangenome (and panpatch) is built for each haplotype -- N times the cost of a normal run -- "
                        "and the seqfile must contain only the sample to patch (first) and its donors")
    parser.add_argument("--sample", nargs='+',
                        help = "Sample to patch, followed by the donor samples in priority order. "
                        "[default: every non-reference sample in the seqfile, in seqfile order]")
    # like cactus-align --batch, the (positional) seqFile is read as a chromfile when this is set.
    # the dest is not "batch" because the pangenome pipeline uses options.batch for its own purposes
    parser.add_argument("--batch", dest="panpatchBatch", action="store_true", default=False,
                        help = "Patch many samples at once. The seqFile argument is instead a chromfile with one "
                        "\"<name> <seqfile>\" line per sample. Every seqfile is assumed to have the same layout, differing "
                        "only in its sample names. Cannot be used with --sample or --outName")

    # panpatch options
    parser.add_argument("--requireTelomeres", action="store_true", default=False,
                        help = "Require a telomere at both ends of each patched haplotype (and none inside). A missing terminal "
                        "telomere is patched in from a donor when possible, otherwise the contig reverts to its input (panpatch -T)")
    parser.add_argument("--maxTelomerePatch", type=int, default=None,
                        help = "Maximum bp of target sequence that a --requireTelomeres telomere patch may replace at a contig end "
                        "(panpatch -M) [default=500000]")
    parser.add_argument("--excludeBed", default=None,
                        help = "BED file of target regions to leave untouched (panpatch -b)")
    parser.add_argument("--defaultSample", default=None,
                        help = "Use this sample's contig when a patch is rejected (panpatch -e)")
    parser.add_argument("--panpatchOptions", type=str, default=None,
                        help = "Options to pass through to panpatch, ex --panpatchOptions \"--min-cover 0.9 --graft-recovery 60\" "
                        "(don't forget to wrap in quotes)")
    parser.add_argument("--keepPangenome", action="store_true", default=False,
                        help = "Keep the intermediate cactus-pangenome output in <outDir>/<name>.cactus-scratch/ instead of "
                        "deleting it once panpatch has succeeded")
    parser.add_argument("--keepGraphs", action="store_true", default=False,
                        help = "Keep the per-chromosome vg graphs panpatch ran on in <outDir>/<name>.chroms/, so panpatch's "
                        "parameters can be re-tuned (via --panpatchOptions) in minutes without rebuilding the pangenome")

    # cactus-pangenome options.  note that we deliberately do *not* pull in graphmap_join_options():
    # cactus-panpatch wants nothing but the chromosome vg graphs out of the pangenome, so there is no
    # output for the user to select.  graphmap_join_defaults() fills those options in below
    pangenome_options(parser)

    # cactus-graphmap-join uses this to cap the memory of the vg jobs that we do run
    parser.add_argument("--indexMemory", type=human2bytesN,
                        help="Memory in bytes to upper-bound the per-chromosome vg jobs. "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))", default=None)

    # hal2vg (which converts each chromosome alignment to vg) is sized from the alignment's memory,
    # which badly under-estimates it for the large unclipped reference-free graphs here (measured up
    # to ~27Gi), so floor it.  standard suffixes supported (default=bytes)
    parser.add_argument("--halExportMemory", type=human2bytesN, default=48 * 2**30,
                        help="Minimum memory in bytes for each per-chromosome hal2vg job [default=48Gi]")

    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    panpatch_validate_options(options)

    # we need the minigraph event name to filter it out of the seqfiles below
    graph_event = getOptionalAttrib(findRequiredNode(ET.parse(options.configFile).getroot(), "graphmap"),
                                    "assemblyName", default="_MINIGRAPH_")

    # work out the runs (one cactus-pangenome + panpatch each) before validating, since with
    # --referenceFree the reference is derived from the sample being patched
    runs = make_runs(options, graph_event)
    if not options.reference:
        # --referenceFree: each run has its own reference (its own haplotype). the validation below
        # just needs a representative one -- every run gets its own via pg_options further down
        options.reference = list(runs[0]['reference'])

    # fill in the graphmap-join options that we don't expose, then ask for the one output we want.
    # --chrom-vg full is not just a convenience: it's what gets the _MINIGRAPH_ paths taken out of
    # the full graphs (see drop_graph_event() in cactus_graphmap_join)
    options = graphmap_join_defaults(options)
    options.chrom_vg = ['full']
    # nothing in a chrom-vg-only workflow runs with --indexCores, but graphmap-join insists on
    # having it when we're not on a single machine
    options.indexCores = 1

    options = pangenome_validate_options(options)

    # graphmap-join's validation switches a default GFA back on whenever no whole-genome output was
    # requested, which in turn keeps the clip phase alive.  we want neither: the "full" graphs come
    # out of the join workflow before any clipping or filtering
    disable_pangenome_outputs(options)

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    for run in runs:
        logger.info('cactus-panpatch will patch {} using donor(s) {} against reference {} (output: {})'.format(
            run['samples'][0], ', '.join(run['samples'][1:]) if len(run['samples']) > 1 else '<none>',
            run['panpatchReference'], run['name']))
    start_time = timeit.default_timer()

    with Toil(options) as toil:
        importSingularityImage(options)
        if options.restart:
            toil.restart()
        else:
            config_node = ET.parse(options.configFile).getroot()
            config_wrapper = ConfigWrapper(config_node)
            config_wrapper.substituteAllPredefinedConstantsWithLiterals(options)
            pangenome_config_overrides(options, config_node)

            # floor the per-chromosome hal2vg memory (read by export_vg in cactus_align): its default
            # estimate under-provisions the large reference-free graphs we make
            findRequiredNode(config_node, "hal2vg").attrib["minMemory"] = str(int(options.halExportMemory))

            if options.collapseRefPAF:
                findRequiredNode(config_node, "graphmap").attrib["collapse"] = 'reference'

            exclude_bed_id = None
            if options.excludeBed:
                exclude_bed_id = toil.importFile(makeURL(options.excludeBed))

            # import the input for each run, and make each run its own cactus-pangenome options,
            # pointed at its own working directory.
            #
            # every run imports its own copy of everything it needs, even when several of them read
            # the same file (the reference is shared by every sample in a --batch, and the donors by
            # every haplotype of a --referenceFree sample).  it is tempting to import each path once
            # and hand the same file ID to each run, but the pangenome workflow *consumes* the
            # fastas it is given: sanitize_fasta_header deletes each input from the jobstore once it
            # has sanitized it (checkUniqueHeaders.py:75).  a shared file ID would therefore be
            # deleted out from under every run but the first
            run_inputs = []
            for run in runs:
                seq_id_map = {}
                seq_path_map = {}
                seq_order = []
                for sample, seq_path in run['seq_rows']:
                    seq_url = makeURL(seq_path)
                    seq_id_map[sample] = toil.importFile(seq_url)
                    seq_path_map[sample] = seq_url
                    seq_order.append(sample)

                ref_collapse_paf_id = toil.importFile(makeURL(options.collapseRefPAF)) if options.collapseRefPAF else None
                last_scores_id = toil.importFile(makeURL(options.scoresFile)) if options.scoresFile else None

                # a *second* import of the target fastas, kept aside for the dropped-contig rescue in
                # run_panpatch.  it can't share the pangenome's copy above: that one is consumed
                # (deleted) by sanitize_fasta_header during construction
                target_fasta_ids = {hap: toil.importFile(makeURL(path)) for hap, path in run['target_haps']}

                pg_options = copy.deepcopy(options)
                pg_options.outDir = pangenome_dir(options, run)
                pg_options.outName = run['name']
                pg_options.reference = list(run['reference'])
                pg_options.seqFile = run['seqFile']
                # we only want the chromosome vgs: don't make a hal we'd throw away, and don't
                # export the vgs into the pangenome directory when we're about to export them
                # ourselves (they're the biggest output there is)
                pg_options.noHal = True
                pg_options.noJoinExport = not options.keepPangenome
                if not pg_options.outDir.startswith('s3://'):
                    os.makedirs(pg_options.outDir, exist_ok=True)

                run_inputs.append((run, pg_options, seq_id_map, seq_path_map, seq_order,
                                   ref_collapse_paf_id, last_scores_id, target_fasta_ids))

            toil.start(Job.wrapJobFn(panpatch_batch_workflow, options, config_wrapper, run_inputs,
                                     exclude_bed_id))

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-panpatch has finished after {} seconds".format(run_time))

def sample_base(sample):
    """ strip the .N haplotype suffix off a sample name to get the name it has in the graph, which is
    the name panpatch wants.  only a *numeric* extension is a haplotype: a base name is allowed to
    contain dots of its own (ex HG002.verkko.1 is haplotype 1 of HG002.verkko), and the names given
    to --sample / --defaultSample are base names already """
    base, ext = os.path.splitext(sample)
    return base if ext and ext[1:].isdigit() else sample

def sample_hap(sample):
    """ the haplotype number of a seqfile sample name, or None if it's haploid (no .N suffix) """
    ext = os.path.splitext(sample)[1]
    return int(ext[1:]) if ext and ext[1:].isdigit() else None

def sanitize_contig_name(header):
    """ mirror cactus_sanitizeFastaHeaders -p on a single fasta header, so an input contig name
    matches the CONTIG field of the graph path names (SAMPLE#HAP#CONTIG).  cut at whitespace, strip
    through the last '#', convert a trailing ':start-end' range to '_sub_(start-1)_end', and replace
    any other ':' with '_' (kept consistent with cactus_sanitizeFastaHeaders.c) """
    name = header.split()[0]
    if '#' in name:
        name = name.rsplit('#', 1)[-1]
    colon = name.rfind(':')
    if colon != -1:
        dash = name.rfind('-')
        if dash > colon + 1 and name[colon+1:dash].isdigit() and name[dash+1:].isdigit() \
           and int(name[dash+1:]) >= int(name[colon+1:dash]):
            start, end = int(name[colon+1:dash]), int(name[dash+1:])
            name = name[:colon] + '_sub_{}_{}'.format(start - 1 if start > 0 else start, end)
        name = name.replace(':', '_')
    return name

def panpatch_validate_options(options):
    """ check the cactus-panpatch options before doing anything expensive """

    # like cactus-align --batch, --batch reads the (positional) seqFile as a chromfile
    if not options.seqFile:
        what = 'chromfile' if options.panpatchBatch else 'seqFile'
        raise RuntimeError('{} argument is required'.format(what))
    if options.panpatchBatch:
        for opt, val in [('--sample', options.sample), ('--outName', options.outName)]:
            if val:
                raise RuntimeError('{} cannot be used with --batch (the chromfile provides it for each sample)'.format(opt))

    # reference-free (each haplotype patched against itself) is the default; giving --reference
    # switches to external-reference mode
    options.referenceFree = not options.reference
    if options.reference:
        # a Minigraph-Cactus reference must be haploid: a .1/.2/... suffix means it is one haplotype
        # of a diploid sample, which cannot be the reference (use reference-free for that)
        for ref in options.reference:
            if sample_hap(ref) not in (None, 0):
                raise RuntimeError('--reference {} must be haploid: it cannot have a .1/.2/... haplotype suffix '
                                   '(use no suffix, or .0)'.format(ref))

    if options.noSplit:
        # a single unsplit graph has every reference chromosome in it, and panpatch skips any graph
        # that doesn't have exactly one reference path -- it would run to completion and hand back
        # the unpatched input
        raise RuntimeError('--noSplit cannot be used with cactus-panpatch: panpatch needs one graph per reference chromosome')

def disable_pangenome_outputs(options):
    """ switch off every graphmap-join output except the chromosome vgs, so that we don't spend
    hours building indexes, VCFs and GFAs that get thrown away.  must be called *after*
    graphmap_join_validate_options(), which turns a default GFA back on if it finds no whole-genome
    output selected """
    for opt in ['gfa', 'unchopped_gfa', 'gbz', 'xg', 'odgi', 'viz', 'draw', 'chrom_og',
                'vcf', 'giraffe', 'lrGiraffe', 'haplo', 'snarlStats', 'panacus']:
        setattr(options, opt, [])
    options.gref = None
    options.grefL = None
    # with no whole-genome output left, nothing downstream of the "full" phase is needed
    options.clip = None
    options.filter = None
    # these two are the only reason graphmap-join would still build a whole-genome gbz (it needs a
    # reference fasta to normalize a VCF with, and checks the config rather than --vcf)
    options.vcfbub = 0
    options.vcfwave = False
    assert options.chrom_vg == ['full']

def read_seqfile(seqfile_path, graph_event):
    """ the (sample, fasta) rows of a seqfile, in order.  the minigraph event is skipped: a seqfile
    that cactus-pangenome exported from a previous run will have it """
    seq_file = SeqFile(seqfile_path)
    leaves = set(seq_file.tree.getName(node) for node in seq_file.tree.getLeaves())
    return [(sample, seq_file.pathMap[sample]) for sample in seq_file.seqOrder
            if sample in leaves and sample != graph_event]

def read_batch_file(batch_path):
    """ read the --batch file: one "<name> <seqfile>" line per sample """
    rows = []
    names = set()
    with open(batch_path, 'r') as batch_file:
        for line in batch_file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            toks = line.split()
            if len(toks) != 2:
                raise RuntimeError('--batch file lines must be "<name> <seqfile>", got: {}'.format(line))
            if toks[0] in names:
                raise RuntimeError('--batch file has more than one row named {}'.format(toks[0]))
            names.add(toks[0])
            rows.append((toks[0], toks[1]))
    if not rows:
        raise RuntimeError('No samples found in --batch file {}'.format(batch_path))
    return rows

def make_runs(options, graph_event):
    """ work out the runs implied by the options: one cactus-pangenome + panpatch each.  a normal
    invocation is a single run, --referenceFree is one run per haplotype of the sample being patched,
    and --batch multiplies all that by the number of samples """
    if options.panpatchBatch:
        # with --batch the seqFile argument is the chromfile (as in cactus-align --batch)
        seqfiles = read_batch_file(options.seqFile)
    else:
        seqfiles = [(options.outName, options.seqFile)]

    runs = []
    for name, seqfile_path in seqfiles:
        runs += make_seqfile_runs(options, name, seqfile_path, graph_event)

    names = [run['name'] for run in runs]
    if len(set(names)) != len(names):
        raise RuntimeError('Output names are not unique: {}'.format(names))
    for run_name in names:
        # the names get pasted into output paths (and into a directory we later delete), so don't
        # let them wander out of --outDir
        if os.path.basename(run_name) != run_name or run_name.startswith('.'):
            raise RuntimeError('Invalid output name "{}": it cannot contain "/" or begin with "."'.format(run_name))
    return runs

def make_seqfile_runs(options, name, seqfile_path, graph_event):
    """ the runs implied by a single seqfile """
    rows = read_seqfile(seqfile_path, graph_event)
    if not rows:
        raise RuntimeError('No sequences found in {}'.format(seqfile_path))
    check_sample_names([row[0] for row in rows], options.reference)

    # the sample names as they appear in the graph, in seqfile order
    samples = []
    for sample, _ in rows:
        if sample_base(sample) not in samples:
            samples.append(sample_base(sample))

    # the reference is never inferred from the seqfile: it's either given with --reference, or (with
    # --referenceFree) it *is* the sample being patched, which is named just below.  note that cactus
    # takes the reference as it's named in the seqfile (which can carry a .N haplotype suffix),
    # whereas panpatch wants the name it has in the graph, which never does
    reference = None
    references = set()
    if not options.referenceFree:
        reference = sample_base(options.reference[0])
        # --reference takes several names (the extra ones get promoted to reference paths in vg).
        # none of them is a candidate for patching
        references = set(sample_base(ref) for ref in options.reference)
        for ref in references:
            if ref not in samples:
                raise RuntimeError('Reference {} not found in {}'.format(ref, seqfile_path))

    # the sample to patch, followed by its donors in priority order
    if options.sample:
        patch_samples = [sample_base(sample) for sample in options.sample]
        for sample in patch_samples:
            if sample not in samples:
                raise RuntimeError('--sample {} not found in {}'.format(sample, seqfile_path))
    else:
        patch_samples = [sample for sample in samples if sample not in references]
    if not patch_samples:
        raise RuntimeError('No sample to patch found in {}'.format(seqfile_path))

    target, donors = patch_samples[0], patch_samples[1:]
    target_rows = [row for row in rows if sample_base(row[0]) == target]
    donor_rows = [row for row in rows if sample_base(row[0]) in donors]

    if not options.referenceFree:
        # the graph is the whole seqfile, with the reference first (as cactus-pangenome wants it)
        seq_rows = [row for row in rows if sample_base(row[0]) == reference]
        seq_rows += [row for row in rows if sample_base(row[0]) != reference]
        # target_haps maps each target haplotype's graph number to its input fasta -- used to rescue
        # dropped contigs (see run_panpatch).  panpatch names its output for haplotype N after the
        # PanSN haplotype of the target path: the seqfile's .N suffix, or 0 for a suffix-less
        # (haploid) target -- which is the number panpatch actually writes (.hap0)
        target_haps = [(sample_hap(row[0]) if sample_hap(row[0]) is not None else 0, row[1])
                       for row in target_rows]
        return [{'name' : name if name else target,
                 'seqFile' : seqfile_path,
                 'seq_rows' : seq_rows,
                 'reference' : options.reference,
                 'panpatchReference' : reference,
                 'samples' : [target] + donors,
                 'defaultSample' : resolve_default_sample(options, target, donors, target),
                 'target_haps' : target_haps,
                 'ploidy' : len(target_rows),
                 'referenceFree' : False,
                 'hap' : None}]

    # --referenceFree: a Minigraph-Cactus reference must be haploid, so we can't just point
    # --reference at TARGET.1 (cactus would read that as haplotype 1 of a diploid TARGET). we rename
    # it to the suffix-free TARGET_1 instead, which makes it a haploid sample in its own right, and
    # build a separate graph -- and run panpatch -- for each haplotype of the target. the other
    # haplotypes of the target (and anything else in the seqfile) are left out of the graph
    runs = []
    for target_row in target_rows:
        hap = sample_hap(target_row[0])
        # the "." is what cactus uses to mark a haplotype, so it can't survive into the new name:
        # SAMPLE.1 -> SAMPLE_1.  a base name is allowed to contain dots of its own (only the last
        # one is the haplotype), and those have to go too, or check_sample_names would read the
        # result as a bogus haplotype suffix (ex HG002.verkko.1 -> HG002_verkko_1)
        graph_name = target.replace('.', '_') if hap is None else '{}_{}'.format(target.replace('.', '_'), hap)
        # the renamed target is a suffix-free (haploid) sample, so it is PanSN haplotype 0 in the
        # graph, and panpatch names its single output file hap0 (mapped to the user haplotype in
        # run_panpatch).  target_haps carries the input fasta for the contig rescue there
        runs.append({'name' : (name if name else target) + ('' if hap is None else '.hap{}'.format(hap)),
                     'seqFile' : seqfile_path,
                     'seq_rows' : [(graph_name, target_row[1])] + donor_rows,
                     'reference' : [graph_name],
                     'panpatchReference' : graph_name,
                     'samples' : [graph_name] + donors,
                     'defaultSample' : resolve_default_sample(options, target, donors, graph_name),
                     'target_haps' : [(0, target_row[1])],
                     'ploidy' : 1,
                     'referenceFree' : True,
                     'hap' : hap})
    return runs

def resolve_default_sample(options, target, donors, target_graph_name):
    """ --defaultSample as it's named in the graph.  it has to be one of the samples that's actually
    in there, and it has to survive the --referenceFree rename if it's the sample being patched --
    panpatch doesn't check -e (unlike -r and -s), it just quietly matches no paths """
    if not options.defaultSample:
        return None
    default_sample = sample_base(options.defaultSample)
    if default_sample == target:
        return target_graph_name
    if default_sample in donors:
        return default_sample
    raise RuntimeError('--defaultSample {} must be the sample being patched ({}) or one of its donors ({})'.format(
        options.defaultSample, target, ', '.join(donors) if donors else '<none>'))

def pangenome_dir(options, run):
    """ where a run's intermediate cactus-pangenome output goes (deleted on success) """
    return os.path.join(options.outDir, run['name'] + '.cactus-scratch')

def panpatch_batch_workflow(job, options, config_wrapper, run_inputs, exclude_bed_id):
    """ patch each sample.  the runs are completely independent of each other: they share no input
    files in the jobstore (see the import loop in main()) """
    for run, pg_options, seq_id_map, seq_path_map, seq_order, ref_collapse_paf_id, last_scores_id, target_fasta_ids in run_inputs:
        job.addChildJobFn(panpatch_run_workflow, options, pg_options, config_wrapper, run, seq_id_map,
                          seq_path_map, seq_order, ref_collapse_paf_id, last_scores_id, exclude_bed_id,
                          target_fasta_ids)

def panpatch_run_workflow(job, options, pg_options, config_wrapper, run, seq_id_map, seq_path_map, seq_order,
                          ref_collapse_paf_id, last_scores_id, exclude_bed_id, target_fasta_ids):
    """ cactus-pangenome, then panpatch on the chromosome graphs it made """
    pangenome_job = job.addChildJobFn(pangenome_end_to_end_workflow, pg_options, config_wrapper, seq_id_map,
                                      seq_path_map, seq_order, ref_collapse_paf_id, last_scores_id)

    # a follow-on of the job hosting the pangenome workflow only runs once that job's entire child
    # subtree is done, which is what makes this safe
    pangenome_job.addFollowOnJobFn(panpatch_workflow, options, run, pangenome_job.rv(0), pangenome_job.rv(1),
                                   pangenome_job.rv(2), exclude_bed_id, target_fasta_ids)

def panpatch_workflow(job, options, run, join_options, join_wf_output, seq_id_map, exclude_bed_id, target_fasta_ids):
    """ the pangenome's promises have resolved by now.  run one single-threaded panpatch job per
    chromosome graph -- panpatch handles each graph independently, so this matches running it over
    all of them at once, but lets toil spread the chromosomes over the cluster -- then gather the
    per-chromosome outputs into the final assembly """

    # the "full" chromosome graphs (cactus-pangenome was run with --chrom-vg full), named the same
    # way cactus-graphmap-join names them on disk
    full_vg_ids = join_wf_output[0]
    vg_names = [os.path.splitext(os.path.basename(vg_path))[0] + '.full.vg' for vg_path in join_options.vg]

    if not full_vg_ids:
        # nothing to patch: the pangenome produced no chromosome graphs (eg the reference had no
        # sequence binned to any chromosome)
        raise RuntimeError('cactus-pangenome produced no chromosome graphs for {}: nothing to patch'.format(run['name']))

    # one single-threaded panpatch job per chromosome graph.  each needs memory for just its own
    # graph, and --indexMemory (if given) caps the estimate
    chrom_jobs = []
    for vg_id, vg_name in zip(full_vg_ids, vg_names):
        mem = cactus_clamp_memory(max(2**32, vg_id.size * 8))
        if options.indexMemory:
            mem = min(mem, options.indexMemory)
        cj = job.addChildJobFn(run_panpatch_chrom, options, run, vg_id, vg_name, exclude_bed_id,
                               cores=1, memory=mem, disk=vg_id.size * 4 + 2**30)
        chrom_jobs.append(cj)

    # gather the per-chromosome outputs, rescue dropped contigs, and bgzip.  the concatenated fastas
    # are ~ploidy * genome; seq_id_map holds the sanitized (uncompressed) inputs, so the reference's
    # size is the genome size
    ref_size = seq_id_map[run['reference'][0]].size
    gather_disk = int(run['ploidy'] * ref_size * 3) + sum(fid.size for fid in target_fasta_ids.values()) + 2**30
    gather_job = job.addFollowOnJobFn(gather_panpatch, options, run, [cj.rv() for cj in chrom_jobs],
                                      target_fasta_ids, memory=cactus_clamp_memory(max(2**32, ref_size * 2)),
                                      disk=gather_disk)
    gather_job.addFollowOnJobFn(export_panpatch_wrapper, options, run, gather_job.rv(), full_vg_ids, vg_names,
                                disk=sum(vg_id.size for vg_id in full_vg_ids) * 2 + 2**30)

def run_panpatch_chrom(job, options, run, vg_id, vg_name, exclude_bed_id):
    """ run single-threaded panpatch on one chromosome graph.  returns the per-haplotype fastas, the
    bed and report, and the set of target contigs this graph held (for the dropped-contig rescue) """
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, vg_name)
    job.fileStore.readGlobalFile(vg_id, vg_path)

    # which samples this graph contains, and which of the target's contigs.  a donor need not have
    # contigs on every chromosome, and panpatch aborts if a -s sample is absent from the graph, so we
    # pass only the samples actually here (a graph with no donors is fine -- panpatch just passes the
    # target through).  this matches what panpatch does over all graphs at once: it only ever uses a
    # sample's paths in the graphs that have them.  the graph path names are PanSN
    # (SAMPLE#HAP#CONTIG[#OFFSET]); sanitize strips '#' out of contig names, so the CONTIG field
    # (index 2) is exactly the first token of the input fasta header
    target_sample = run['samples'][0]
    present_samples = set()
    # the target contigs this graph holds, keyed "HAP#CONTIG": a diploid target's two haplotypes
    # share the graphs and can reuse a contig name, so the rescue below must not treat hap2's dropped
    # copy as "seen" just because hap1's copy is here
    target_contigs = set()
    for line in cactus_call(parameters=['vg', 'paths', '-L', '-x', vg_path], check_output=True).split('\n'):
        fields = line.split('#')
        if len(fields) >= 3:
            present_samples.add(fields[0])
            if fields[0] == target_sample:
                target_contigs.add('{}#{}'.format(fields[1], fields[2]))

    if target_sample not in present_samples:
        # no target sequence in this graph (eg a bin of only donor or unplaced contigs): nothing to
        # patch or emit here (any target contig that belongs here but was dropped is rescued later)
        return {'fastas' : {}, 'bed' : None, 'report' : None, 'seen' : sorted(target_contigs)}

    samples = [s for s in run['samples'] if s in present_samples]   # target first, order preserved

    # panpatch's -f is a template, not a prefix: it puts .hap<N> in front of the extension, so the
    # .fa is what makes the output panpatch.hap1.fa rather than an extension-less panpatch.hap1
    fasta_path = os.path.join(work_dir, 'panpatch.fa')
    bed_path = os.path.join(work_dir, 'panpatch.bed')
    report_path = os.path.join(work_dir, 'panpatch.tsv')

    cmd = ['panpatch', vg_path, '-r', run['panpatchReference'], '-t', '1']
    for sample in samples:
        cmd += ['-s', sample]
    cmd += ['-f', fasta_path, '--bed', bed_path]
    if options.requireTelomeres:
        cmd += ['-T']
    if options.maxTelomerePatch is not None:
        cmd += ['-M', str(options.maxTelomerePatch)]
    if run['defaultSample'] and run['defaultSample'] in present_samples:
        cmd += ['-e', run['defaultSample']]
    if exclude_bed_id:
        exclude_bed_path = os.path.join(work_dir, 'exclude.bed')
        job.fileStore.readGlobalFile(exclude_bed_id, exclude_bed_path)
        cmd += ['-b', exclude_bed_path]
    if options.panpatchOptions:
        cmd += shlex.split(options.panpatchOptions)

    cactus_call(parameters=cmd, outfile=report_path, job_memory=job.memory)

    # one fasta per haplotype panpatch found target paths for, keyed by the PanSN haplotype number
    fastas = {}
    for hap_fasta in glob.glob(os.path.join(work_dir, 'panpatch.hap*.fa')):
        hap = int(os.path.basename(hap_fasta)[len('panpatch.hap'):-len('.fa')])
        fastas[hap] = job.fileStore.writeGlobalFile(hap_fasta)

    return {'fastas' : fastas,
            'bed' : job.fileStore.writeGlobalFile(bed_path),
            'report' : job.fileStore.writeGlobalFile(report_path),
            'seen' : sorted(target_contigs)}

def gather_panpatch(job, options, run, chrom_results, target_fasta_ids):
    """ concatenate the per-chromosome panpatch outputs (in chromosome order, as panpatch itself
    would), rescue any target contig the pangenome dropped, bgzip, and return the output file map """
    work_dir = job.fileStore.getLocalTempDir()

    # every target contig panpatch saw across all the chromosome graphs.  a contig the pangenome
    # dropped (to the excluded _AMBIGUOUS_ bin, or unmapped) is in none of them, so it is rescued below
    seen_contigs = set()
    for res in chrom_results:
        seen_contigs.update(res['seen'])

    # the haplotypes to emit: every one panpatch produced a fasta for, plus every target haplotype we
    # have an input for -- so a haplotype whose contigs were *all* dropped still gets an output (an
    # empty combined fasta that the rescue fills), rather than silently vanishing
    emitted = set(hap for res in chrom_results for hap in res['fastas'])
    haps = sorted(set(target_fasta_ids) | emitted)
    if not emitted:
        RealtimeLogger.warning('panpatch patched nothing for {}: output will be the input contigs unchanged'.format(
            run['name']))
    if run['referenceFree'] and haps != [0]:
        # the reference-free target is a haploid sample (PanSN hap 0) in its own graph
        raise RuntimeError('reference-free panpatch produced haplotype(s) {} for {}, expected only hap0'.format(
            haps, run['name']))

    output_id_map = {}
    for hap in haps:
        combined = os.path.join(work_dir, 'hap{}.fa'.format(hap))
        with open(combined, 'wb') as out:
            for res in chrom_results:   # chromosome order
                if hap in res['fastas']:
                    chrom_fa = os.path.join(work_dir, 'chrom.fa')
                    job.fileStore.readGlobalFile(res['fastas'][hap], chrom_fa)
                    with open(chrom_fa, 'rb') as f:
                        shutil.copyfileobj(f, out)
                    os.remove(chrom_fa)

        # rescue any target contig the pangenome dropped, appending it verbatim from the input assembly
        if hap in target_fasta_ids:
            rescue_dropped_contigs(job, combined, target_fasta_ids[hap], seen_contigs, run['samples'][0], hap)

        # in reference-free mode the single (hap0) output is named for the user's haplotype (already
        # in the run name); otherwise keep panpatch's haplotype numbering
        output_name = (run['name'] + '.fa.gz') if run['referenceFree'] else (run['name'] + '.hap{}.fa.gz'.format(hap))
        cactus_call(parameters=['bgzip', combined, '--threads', str(int(job.cores))])
        output_id_map[output_name] = job.fileStore.writeGlobalFile(combined + '.gz')

    # concatenate the beds (in chromosome order; a skipped graph has no bed)
    bed_out = os.path.join(work_dir, 'patched.bed')
    with open(bed_out, 'wb') as out:
        for res in chrom_results:
            if res['bed'] is None:
                continue
            p = os.path.join(work_dir, 'c.bed')
            job.fileStore.readGlobalFile(res['bed'], p)
            with open(p, 'rb') as f:
                shutil.copyfileobj(f, out)
            os.remove(p)
    output_id_map[run['name'] + '.bed'] = job.fileStore.writeGlobalFile(bed_out)

    # concatenate the reports: the TSV header from the first chromosome that has one, then every
    # chromosome's rows (a skipped graph has no report)
    report_out = os.path.join(work_dir, 'patched.tsv')
    with open(report_out, 'w') as out:
        wrote_header = False
        for res in chrom_results:
            if res['report'] is None:
                continue
            p = os.path.join(work_dir, 'c.tsv')
            job.fileStore.readGlobalFile(res['report'], p)
            with open(p) as f:
                lines = f.readlines()
            out.writelines(lines if not wrote_header else lines[1:])
            wrote_header = True
            os.remove(p)
    output_id_map[run['name'] + '.tsv'] = job.fileStore.writeGlobalFile(report_out)

    return output_id_map

def rescue_dropped_contigs(job, fasta_path, target_fasta_id, seen_contigs, target_sample, hap):
    """ append to fasta_path (a panpatch output haplotype) any contig from the input target assembly
    that panpatch never saw -- ie that the pangenome dropped (binned to _AMBIGUOUS_ or unmapped).
    those contigs were never patched or scaffolded, so they are safe to carry through verbatim """
    work_dir = job.fileStore.getLocalTempDir()
    in_fa = os.path.join(work_dir, 'target_input.fa')
    job.fileStore.readGlobalFile(target_fasta_id, in_fa)
    with open(in_fa, 'rb') as fh:
        is_gzipped = fh.read(2) == b'\x1f\x8b'
    opener = gzip.open if is_gzipped else open

    rescued = 0
    with open(fasta_path, 'a') as out, opener(in_fa, 'rt') as inp:
        keep = False
        for line in inp:
            if line.startswith('>'):
                # sanitize the input header the same way cactus did before building the graph, so the
                # name matches what "seen" holds (and the passthrough contigs already in the output)
                contig = sanitize_contig_name(line[1:])
                keep = '{}#{}'.format(hap, contig) not in seen_contigs
                if keep:
                    # match panpatch's passthrough header style (SAMPLE#HAP#CONTIG)
                    out.write('>{}#{}#{}\n'.format(target_sample, hap, contig))
                    rescued += 1
            elif keep:
                out.write(line)
    RealtimeLogger.info('Rescued {} dropped contig(s) into {}'.format(rescued, os.path.basename(fasta_path)))

def export_panpatch_wrapper(job, options, run, output_id_map, full_vg_ids, vg_names):
    """ write the panpatch output to the output directory, and (with --keepGraphs) the graphs it was
    made from """
    for output_name, output_id in output_id_map.items():
        job.fileStore.exportFile(output_id, makeURL(os.path.join(options.outDir, output_name)))

    # optionally keep panpatch's input around: re-running it with different parameters takes minutes,
    # whereas rebuilding the pangenome takes hours
    if options.keepGraphs:
        chrom_dir = os.path.join(options.outDir, run['name'] + '.chroms')
        if not chrom_dir.startswith('s3://') and not os.path.isdir(chrom_dir):
            os.makedirs(chrom_dir)
        for vg_id, vg_name in zip(full_vg_ids, vg_names):
            job.fileStore.exportFile(vg_id, makeURL(os.path.join(chrom_dir, vg_name)))

    if not options.keepPangenome:
        job.addFollowOnJobFn(cleanup_pangenome_wrapper, options, run)

def cleanup_pangenome_wrapper(job, options, run):
    """ the cactus-pangenome output is only ever scratch for us.  everything the user asked for has
    been exported by the time we get here, so a failure anywhere upstream leaves it alone (and
    leaves it ready for --restart) """
    pg_dir = pangenome_dir(options, run)
    RealtimeLogger.info('Removing intermediate cactus-pangenome output: {}'.format(pg_dir))
    delete_directory(pg_dir)

if __name__ == "__main__":
    main()
