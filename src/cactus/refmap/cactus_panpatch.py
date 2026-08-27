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
    parser.add_argument("--assemblyErrorBeds", default=None,
                        help = "Manifest of per-assembly error BEDs: one \"<seqfile-event-name>\\t<bed-path>\" line "
                        "per assembly (a subset is fine). Intervals are treated almost like runs of Ns: in a DONOR "
                        "assembly they are never used to patch; in the TARGET they are patched like gaps when possible "
                        "and otherwise left as the original sequence (never N). BED contig names are the assembly's own "
                        "fasta-header first token, coordinates 0-based half-open in that assembly's frame. Distinct from "
                        "--excludeBed (which protects target regions and uses full graph path names)")
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
    # to ~27Gi with hal2vg 1.1.8, less since), so floor it.  standard suffixes supported (default=bytes)
    parser.add_argument("--halExportMemory", type=human2bytesN, default=32 * 2**30,
                        help="Minimum memory in bytes for each per-chromosome hal2vg job [default=32Gi]")

    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    panpatch_validate_options(options)

    # per-assembly error BEDs (treated almost like runs of Ns): parse the manifest up front so run
    # construction (make_seqfile_runs) can attach each BED to its assembly, and so a typo'd event name
    # is caught before anything expensive runs
    options.errorBedMap = read_error_bed_manifest(options.assemblyErrorBeds) if options.assemblyErrorBeds else {}

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

                # per-assembly error BEDs: import each path once, then key it by seq_rows sample name
                # (for masking the pangenome inputs) and by graph haplotype for the target (for report
                # tagging in run_panpatch_chrom)
                error_bed_id_cache = {}
                for bed_path in set(run['mask_beds'].values()) | set(run['target_error_beds'].values()):
                    error_bed_id_cache[bed_path] = toil.importFile(makeURL(bed_path))
                mask_bed_ids = {sample: error_bed_id_cache[p] for sample, p in run['mask_beds'].items()}
                target_error_bed_ids = {hap: error_bed_id_cache[p] for hap, p in run['target_error_beds'].items()}

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
                                   ref_collapse_paf_id, last_scores_id, target_fasta_ids,
                                   mask_bed_ids, target_error_bed_ids))

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

def read_error_bed_manifest(path):
    """ read the --assemblyErrorBeds manifest: one "<seqfile-event-name>\t<bed-path>" line per assembly
    (a subset of the assemblies is fine).  returns {event_name: absolute_bed_path} """
    error_beds = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            toks = line.split()
            if len(toks) != 2:
                raise RuntimeError('--assemblyErrorBeds lines must be "<event> <bedpath>", got: {}'.format(line))
            event, bed_path = toks[0], toks[1]
            if event in error_beds:
                raise RuntimeError('--assemblyErrorBeds has more than one entry for {}'.format(event))
            if not os.path.isfile(bed_path):
                raise RuntimeError('--assemblyErrorBeds file for {} not found: {}'.format(event, bed_path))
            error_beds[event] = os.path.abspath(bed_path)
    return error_beds

def parse_error_bed(bed_path, sanitize=False):
    """ read a 3-column error BED into {contig: [(start, end), ...]} (0-based half-open).  with
    sanitize=True the contig names are put through sanitize_contig_name so they line up with the graph
    CONTIG field (used when comparing against panpatch's report/output coordinates) """
    intervals = {}
    with open(bed_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line[0] == '#' or line.startswith('track') or line.startswith('browser'):
                continue
            toks = line.split()
            if len(toks) < 3:
                continue
            contig, start, end = toks[0], int(toks[1]), int(toks[2])
            if start < 0 or end < 0 or start >= end:
                continue
            if sanitize:
                contig = sanitize_contig_name(contig)
            intervals.setdefault(contig, []).append((start, end))
    return intervals

_REVCOMP = bytes.maketrans(b'ACGTNacgtnMRWSYKVHDBmrwsykvhdb', b'TGCANtgcanKYWSRMBDHVkywsrmbdhv')
def revcomp(seq):
    """ reverse-complement a bytes sequence (IUPAC-aware) """
    return seq.translate(_REVCOMP)[::-1]

def intervals_overlap(a_start, a_end, intervals):
    """ does [a_start, a_end) overlap any [s, e) in the list of half-open intervals? """
    return any(s < a_end and a_start < e for s, e in intervals)

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

    # per-assembly error BEDs default to none, so make_runs / make_seqfile_runs can read errorBedMap
    # even when main has not parsed a manifest (e.g. the offline unit tests); main overwrites this with
    # the parsed manifest right after calling us
    if not hasattr(options, 'errorBedMap'):
        options.errorBedMap = {}

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
    options.vcfL = None
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

    # every event named in --assemblyErrorBeds must have matched an assembly in some run, else it's a
    # typo (or a reference, which is never masked) -- catch it before anything expensive runs.  compare
    # by the consumed event names (run['error_events']), not by bed path: two events may share one bed
    # file, and the reference-free target's mask_beds key is its graph_name, not its manifest event name
    if options.errorBedMap:
        matched = set()
        for run in runs:
            matched.update(run['error_events'])
        unmatched = sorted(set(options.errorBedMap) - matched)
        if unmatched:
            raise RuntimeError('--assemblyErrorBeds event(s) did not match any patchable assembly (the target or a '
                               'donor) in any seqfile -- check the names (the reference is never masked): {}'.format(
                                   ', '.join(unmatched)))
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
        # attach per-assembly error BEDs (see make_seqfile_runs's --referenceFree branch for the keying
        # rationale).  in reference mode the seq_rows keys are the original event names, so mask_beds is
        # keyed by those; target_error_beds is keyed by graph haplotype for report tagging
        errmap = options.errorBedMap
        mask_beds = {row[0]: errmap[row[0]] for row in target_rows + donor_rows if row[0] in errmap}
        error_events = {row[0] for row in target_rows + donor_rows if row[0] in errmap}
        target_error_beds = {(sample_hap(row[0]) if sample_hap(row[0]) is not None else 0): errmap[row[0]]
                             for row in target_rows if row[0] in errmap}
        return [{'name' : name if name else target,
                 'seqFile' : seqfile_path,
                 'seq_rows' : seq_rows,
                 'reference' : options.reference,
                 'panpatchReference' : reference,
                 'samples' : [target] + donors,
                 'defaultSample' : resolve_default_sample(options, target, donors, target),
                 'target_haps' : target_haps,
                 'mask_beds' : mask_beds,
                 'error_events' : error_events,
                 'target_error_beds' : target_error_beds,
                 'ploidy' : len(target_rows),
                 'referenceFree' : False,
                 'hap' : None}]

    # --referenceFree: a Minigraph-Cactus reference must be haploid, so we can't just point
    # --reference at TARGET.1 (cactus would read that as haplotype 1 of a diploid TARGET). we rename
    # it to the suffix-free TARGET_1 instead, which makes it a haploid sample in its own right, and
    # build a separate graph -- and run panpatch -- for each haplotype of the target. the other
    # haplotypes of the target (and anything else in the seqfile) are left out of the graph
    runs = []
    errmap = options.errorBedMap
    for target_row in target_rows:
        hap = sample_hap(target_row[0])
        # the "." is what cactus uses to mark a haplotype, so it can't survive into the new name:
        # SAMPLE.1 -> SAMPLE_1.  a base name is allowed to contain dots of its own (only the last
        # one is the haplotype), and those have to go too, or check_sample_names would read the
        # result as a bogus haplotype suffix (ex HG002.verkko.1 -> HG002_verkko_1)
        graph_name = target.replace('.', '_') if hap is None else '{}_{}'.format(target.replace('.', '_'), hap)
        # per-assembly error BEDs.  the target is renamed to graph_name in seq_rows so its mask BED is
        # keyed by graph_name; donors keep their original names.  the reference-free target is graph
        # haplotype 0, so its report-tagging BED is keyed by 0
        mask_beds = {dr[0]: errmap[dr[0]] for dr in donor_rows if dr[0] in errmap}
        error_events = {dr[0] for dr in donor_rows if dr[0] in errmap}
        target_error_beds = {}
        if target_row[0] in errmap:
            mask_beds[graph_name] = errmap[target_row[0]]
            target_error_beds[0] = errmap[target_row[0]]
            error_events.add(target_row[0])
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
                     'mask_beds' : mask_beds,
                     'error_events' : error_events,
                     'target_error_beds' : target_error_beds,
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
    summaries = []
    for run, pg_options, seq_id_map, seq_path_map, seq_order, ref_collapse_paf_id, last_scores_id, target_fasta_ids, mask_bed_ids, target_error_bed_ids in run_inputs:
        run_job = job.addChildJobFn(panpatch_run_workflow, options, pg_options, config_wrapper, run, seq_id_map,
                                    seq_path_map, seq_order, ref_collapse_paf_id, last_scores_id, exclude_bed_id,
                                    target_fasta_ids, mask_bed_ids, target_error_bed_ids)
        summaries.append((run['name'], run_job.rv()))
    # once every sample is done, roll the per-sample reports up into one cross-sample summary
    job.addFollowOnJobFn(write_batch_summary, options, summaries)

def mask_assembly_errors(job, fasta_id, bed_id):
    """ return a copy of an input assembly fasta with each error-BED interval replaced by an equal-length
    run of N (length-preserving, so downstream coordinates are unchanged).  operates on the raw fasta --
    headers untouched, sequence only -- before cactus sanitizes it, matching each BED contig to the fasta
    header's first whitespace token (the user's own contig names).  masking a target region makes it a gap
    panpatch can fill (and revert if it can't); masking a donor region breaks that donor's path so its
    error sequence can never be selected """
    work_dir = job.fileStore.getLocalTempDir()
    in_fa = os.path.join(work_dir, 'in.fa')
    job.fileStore.readGlobalFile(fasta_id, in_fa)
    bed_path = os.path.join(work_dir, 'errors.bed')
    job.fileStore.readGlobalFile(bed_id, bed_path)
    intervals = parse_error_bed(bed_path)   # {raw_contig_name: [(start, end)]}

    with open(in_fa, 'rb') as fh:
        is_gzipped = fh.read(2) == b'\x1f\x8b'
    in_opener = gzip.open if is_gzipped else open

    # the masked file's only consumer is the pangenome's sanitize step, which detects compression by
    # magic bytes (checkUniqueHeaders.py sniffs the first two bytes), so we always emit plaintext: it is
    # simpler and avoids recompressing a multi-GB assembly that sanitize immediately decompresses again
    out_fa = os.path.join(work_dir, 'masked.fa')
    width = 80
    masked_bp = [0]

    def emit(out, header, seq):
        if header is None:
            return
        hparts = header[1:].split()
        contig = hparts[0] if hparts else ''
        for s, e in intervals.get(contig, []):
            e = min(e, len(seq))
            if 0 <= s < e:
                seq[s:e] = b'N' * (e - s)
                masked_bp[0] += e - s
        out.write(header.encode('utf-8') + b'\n')
        for i in range(0, len(seq), width):
            out.write(bytes(seq[i:i + width]) + b'\n')

    with in_opener(in_fa, 'rt') as inp, open(out_fa, 'wb') as out:
        header, seq = None, bytearray()
        for line in inp:
            if line.startswith('>'):
                emit(out, header, seq)
                header, seq = line.rstrip('\n'), bytearray()
            else:
                seq += line.strip().encode('ascii')
        emit(out, header, seq)

    RealtimeLogger.info('mask_assembly_errors: masked {} bp across {} contig(s)'.format(masked_bp[0], len(intervals)))
    return job.fileStore.writeGlobalFile(out_fa)

def panpatch_run_workflow(job, options, pg_options, config_wrapper, run, seq_id_map, seq_path_map, seq_order,
                          ref_collapse_paf_id, last_scores_id, exclude_bed_id, target_fasta_ids,
                          mask_bed_ids, target_error_bed_ids):
    """ (optionally mask assembly-error intervals to N, then) cactus-pangenome, then panpatch on the
    chromosome graphs it made """
    # mask any assembly-error intervals to N in the pangenome's inputs, before it is built.  masking a
    # target error region turns it into a genuine gap (filled if a donor spans it, otherwise reverted to
    # the original sequence in gather_panpatch); masking a donor error region breaks that donor's path
    # there so its error sequence can never be selected.  the pristine target (target_fasta_ids) is not
    # masked -- it is what the revert reaches back to
    masked = dict(seq_id_map)
    for sample, bed_id in mask_bed_ids.items():
        mask_job = job.addChildJobFn(mask_assembly_errors, seq_id_map[sample], bed_id,
                                     disk=seq_id_map[sample].size * 4 + 2**30,
                                     memory=cactus_clamp_memory(max(2**32, seq_id_map[sample].size * 2)))
        masked[sample] = mask_job.rv()

    # the pangenome runs as a follow-on so any masking children have completed first (their promises
    # resolved); a follow-on of the job hosting the pangenome then runs once its entire child subtree is
    # done, which is what makes panpatch_workflow safe
    pangenome_job = job.addFollowOnJobFn(pangenome_end_to_end_workflow, pg_options, config_wrapper, masked,
                                         seq_path_map, seq_order, ref_collapse_paf_id, last_scores_id)
    # return the run's output-file map (its report id feeds the cross-sample batch summary)
    return pangenome_job.addFollowOnJobFn(panpatch_workflow, options, run, pangenome_job.rv(0), pangenome_job.rv(1),
                                          pangenome_job.rv(2), exclude_bed_id, target_fasta_ids,
                                          target_error_bed_ids).rv()

def panpatch_workflow(job, options, run, join_options, join_wf_output, seq_id_map, exclude_bed_id, target_fasta_ids,
                      target_error_bed_ids):
    """ the pangenome's promises have resolved by now.  run one single-threaded panpatch job per
    chromosome graph -- panpatch handles each graph independently, so this matches running it over
    all of them at once, but lets toil spread the chromosomes over the cluster -- then gather the
    per-chromosome outputs into the final assembly """

    # the "full" chromosome graphs (cactus-pangenome was run with --chrom-vg full), named the same
    # way cactus-graphmap-join names them on disk
    full_vg_ids = join_wf_output[0]
    full_vg_empty = join_wf_output[7]
    vg_names = [os.path.splitext(os.path.basename(vg_path))[0] + '.full.vg' for vg_path in join_options.vg]

    # drop empty chromosome graphs: a "clean" reference-free haplotype can leave an unplaced
    # (chrOther) bin that is only minigraph, so once drop_graph_event strips it the graph is
    # empty/unloadable.  filter it here so it is neither handed to panpatch (which would crash
    # loading it) nor written out by the --keepGraphs export
    kept = [(vid, name) for vid, name, empty in zip(full_vg_ids, vg_names, full_vg_empty) if not empty]
    if len(kept) < len(full_vg_ids):
        RealtimeLogger.info('cactus-panpatch: skipping {} empty chromosome graph(s) for {}'.format(
            len(full_vg_ids) - len(kept), run['name']))
    full_vg_ids = [k[0] for k in kept]
    vg_names = [k[1] for k in kept]

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
                               target_error_bed_ids,
                               cores=1, memory=mem, disk=vg_id.size * 4 + 2**30)
        chrom_jobs.append(cj)

    # gather the per-chromosome outputs, rescue dropped contigs, and bgzip.  the concatenated fastas
    # are ~ploidy * genome; seq_id_map holds the sanitized (uncompressed) inputs, so the reference's
    # size is the genome size
    ref_size = seq_id_map[run['reference'][0]].size
    gather_disk = int(run['ploidy'] * ref_size * 3) + sum(fid.size for fid in target_fasta_ids.values()) + 2**30
    # reverting target bytes to the original loads one pristine target haplotype into memory, so give
    # gather more headroom when a revert will actually run
    gather_mem = ref_size * (3 if target_error_bed_ids else 2)
    gather_job = job.addFollowOnJobFn(gather_panpatch, options, run, [cj.rv() for cj in chrom_jobs],
                                      target_fasta_ids, target_error_bed_ids,
                                      memory=cactus_clamp_memory(max(2**32, gather_mem)),
                                      disk=gather_disk)
    gather_job.addFollowOnJobFn(export_panpatch_wrapper, options, run, gather_job.rv(), full_vg_ids, vg_names,
                                disk=sum(vg_id.size for vg_id in full_vg_ids) * 2 + 2**30)
    # surface this run's output-file map (incl. its report id) up to the cross-sample batch summary
    return gather_job.rv()

def tag_report_gap_origin(job, report_path, target_error_bed_ids):
    """ add a 'gap_origin' column to panpatch's per-chromosome report, distinguishing a genuine N-gap
    fill ('N-gap') from an error-BED-induced fill ('bed-gap'); '.' for non-gap-fill rows.  panpatch
    cannot tell them apart (both are just Ns in the graph) -- we can, because we masked the bed-gaps: a
    gap-fill row is a bed-gap iff its replaced target span overlaps a target error interval.  the column
    is appended so existing column positions are unchanged """
    # target error intervals per graph haplotype, in sanitized (graph) contig coordinates
    err_by_hap = {}
    if target_error_bed_ids:
        work_dir = job.fileStore.getLocalTempDir()
        for hap, bed_id in target_error_bed_ids.items():
            p = os.path.join(work_dir, 'target_err_{}.bed'.format(hap))
            job.fileStore.readGlobalFile(bed_id, p)
            err_by_hap[hap] = parse_error_bed(p, sanitize=True)

    with open(report_path) as f:
        lines = f.readlines()
    out_lines = []
    cols = None
    idx = {}
    for line in lines:
        stripped = line.rstrip('\n')
        if stripped.startswith('#') or not stripped:
            out_lines.append(line)
            continue
        if cols is None:
            # the header (first non-comment line)
            cols = stripped.split('\t')
            idx = {c: i for i, c in enumerate(cols)}
            out_lines.append(stripped + '\tgap_origin\n')
            continue
        fields = stripped.split('\t')
        origin = '.'
        if len(fields) == len(cols) and fields[idx['type']] == 'gap-fill':
            origin = 'N-gap'
            try:
                hap = int(fields[idx['hap']])
                tgt = fields[idx['target']]
                contig = tgt.split('#')[2] if tgt.count('#') >= 2 else tgt
                ts, te = fields[idx['target_start']], fields[idx['target_end']]
                if ts not in ('.', '') and te not in ('.', '') \
                   and intervals_overlap(int(ts), int(te), err_by_hap.get(hap, {}).get(contig, [])):
                    origin = 'bed-gap'
            except (ValueError, KeyError, IndexError):
                pass
        out_lines.append(stripped + '\t' + origin + '\n')
    with open(report_path, 'w') as f:
        f.writelines(out_lines)

def run_panpatch_chrom(job, options, run, vg_id, vg_name, exclude_bed_id, target_error_bed_ids):
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

    # tag each gap-fill row as a genuine N-gap or an error-BED-induced ('bed-gap') fill
    tag_report_gap_origin(job, report_path, target_error_bed_ids)

    # one fasta per haplotype panpatch found target paths for, keyed by the PanSN haplotype number
    fastas = {}
    for hap_fasta in glob.glob(os.path.join(work_dir, 'panpatch.hap*.fa')):
        hap = int(os.path.basename(hap_fasta)[len('panpatch.hap'):-len('.fa')])
        fastas[hap] = job.fileStore.writeGlobalFile(hap_fasta)

    return {'fastas' : fastas,
            'bed' : job.fileStore.writeGlobalFile(bed_path),
            'report' : job.fileStore.writeGlobalFile(report_path),
            'seen' : sorted(target_contigs)}

def parse_bed_blocks(bed_path):
    """ parse a panpatch --bed into {(locus, hap): [(path_name, start, end, strand), ...]}, one entry per
    '#Patched assembly on <locus> for <sample>#<hap>:' block.  the target revert walks a patched record's
    intervals with this; passthrough / reverted records are handled whole-contig by name and their blocks
    are simply never looked up """
    blocks = {}
    cur = None
    with open(bed_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#Patched assembly on '):
                rest = line[len('#Patched assembly on '):]      # "<locus> for <sample>#<hap>:"
                if ' for ' in rest:
                    locus, after = rest.rsplit(' for ', 1)
                    cur = (locus, after.rstrip(':').split('#')[-1])
                    blocks.setdefault(cur, [])
                else:
                    cur = None
            elif line.startswith('#') or not line:
                cur = None
            elif cur is not None:
                toks = line.split('\t')
                if len(toks) >= 4:
                    blocks[cur].append((toks[0], int(toks[1]), int(toks[2]), toks[3]))
    return blocks

def revert_target_to_original(job, combined_fasta, patched_blocks, pristine_target_id, target_sample, hap):
    """ rewrite a panpatch output haplotype so every target-sample byte comes from the original (unmasked)
    assembly, while donor bytes are kept verbatim.  this restores a masked-but-unpatched error region (and
    any reverted / passthrough target contig, whose graph copy may carry masked Ns) to its original
    sequence -- never N.  driven by the panpatch --bed, which is a byte-exact description of the output:
    output = original_contig[start:end] (revcomp iff strand '-'), concatenated in order """
    work_dir = job.fileStore.getLocalTempDir()
    pristine_fa = os.path.join(work_dir, 'pristine_target.fa')
    job.fileStore.readGlobalFile(pristine_target_id, pristine_fa)
    with open(pristine_fa, 'rb') as fh:
        is_gzipped = fh.read(2) == b'\x1f\x8b'
    opener = gzip.open if is_gzipped else open

    # load the pristine target haplotype, keyed by the graph CONTIG name (sanitized the way cactus did),
    # so it lines up with the PanSN path names in the panpatch bed / fasta (exactly as rescue does)
    pristine = {}
    with opener(pristine_fa, 'rt') as inp:
        cname, parts = None, []
        for line in inp:
            if line.startswith('>'):
                if cname is not None:
                    pristine[cname] = b''.join(parts)
                cname, parts = sanitize_contig_name(line[1:]), []
            else:
                parts.append(line.strip().encode('ascii'))
        if cname is not None:
            pristine[cname] = b''.join(parts)

    def contig_of(path_name):
        # SAMPLE#HAP#CONTIG -> CONTIG (sanitize stripped through the last '#', so CONTIG has none)
        return path_name.split('#', 2)[2] if path_name.count('#') >= 2 else path_name

    def is_target(path_name):
        # sample must match AND the contig must be one we hold, so a donor contig that happens to share a
        # target contig's bare name is never mistaken for the target
        return path_name.split('#')[0] == target_sample and contig_of(path_name) in pristine

    out_path = combined_fasta + '.reverted'
    width, reverted_bp = 80, [0]
    with open(combined_fasta) as inp, open(out_path, 'wb') as out:
        name, seq_parts = None, []

        def flush():
            if name is None:
                return
            if '#' in name:
                # a reverted / passthrough record: one whole target contig, restored from the original
                if is_target(name):
                    new = pristine[contig_of(name)]
                    reverted_bp[0] += len(new)
                else:
                    new = ''.join(seq_parts).encode('ascii')
            else:
                # a patched record "<locus>_hap_<hap>": rebuild from its bed block -- target spans from the
                # pristine assembly, donor spans kept verbatim from panpatch's own emitted bytes
                raw = ''.join(seq_parts).encode('ascii')
                lines = patched_blocks.get((name[:name.rfind('_hap_')], str(hap)))
                if lines is None:
                    new = raw
                else:
                    cursor, chunks = 0, []
                    for (path_name, s, e, strand) in lines:
                        length = e - s
                        if is_target(path_name):
                            seg = pristine[contig_of(path_name)][s:e]
                            chunks.append(revcomp(seg) if strand == '-' else seg)
                            reverted_bp[0] += length
                        else:
                            chunks.append(raw[cursor:cursor + length])
                        cursor += length
                    if cursor != len(raw):
                        raise RuntimeError('panpatch revert: bed length {} != fasta record length {} for {} (hap {})'.format(
                            cursor, len(raw), name, hap))
                    new = b''.join(chunks)
            out.write(('>' + name + '\n').encode('utf-8'))
            for i in range(0, len(new), width):
                out.write(new[i:i + width] + b'\n')

        for line in inp:
            if line.startswith('>'):
                flush()
                name, seq_parts = line[1:].rstrip('\n'), []
            else:
                seq_parts.append(line.strip())
        flush()
    os.replace(out_path, combined_fasta)
    RealtimeLogger.info('panpatch revert: restored {} target bp from the original assembly in hap {}'.format(
        reverted_bp[0], hap))

def gather_panpatch(job, options, run, chrom_results, target_fasta_ids, target_error_bed_ids):
    """ concatenate the per-chromosome panpatch outputs (in chromosome order, as panpatch itself would),
    revert target bytes to the original assembly where the target was masked, rescue any target contig the
    pangenome dropped, bgzip, and return the output file map """
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

    # concatenate the beds first (chromosome order; a skipped graph has no bed).  the bed is a byte-exact
    # description of the output and drives the per-haplotype target revert below
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
    patched_blocks = parse_bed_blocks(bed_out) if target_error_bed_ids else {}

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

        # revert target-sample bytes to the original (unmasked) assembly, so a masked-but-unpatched error
        # region comes back as its original sequence -- never N.  only for a haplotype whose target was
        # masked; donor bytes stay exactly as panpatch emitted them
        if hap in target_error_bed_ids and hap in target_fasta_ids:
            revert_target_to_original(job, combined, patched_blocks, target_fasta_ids[hap],
                                      run['samples'][0], hap)

        # rescue any target contig the pangenome dropped, appending it verbatim from the input assembly
        if hap in target_fasta_ids:
            rescue_dropped_contigs(job, combined, target_fasta_ids[hap], seen_contigs, run['samples'][0], hap)

        # in reference-free mode the single (hap0) output is named for the user's haplotype (already
        # in the run name); otherwise keep panpatch's haplotype numbering
        output_name = (run['name'] + '.fa.gz') if run['referenceFree'] else (run['name'] + '.hap{}.fa.gz'.format(hap))
        cactus_call(parameters=['bgzip', combined, '--threads', str(int(job.cores))])
        output_id_map[output_name] = job.fileStore.writeGlobalFile(combined + '.gz')

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

def summarize_report(report_path):
    """ tally one sample's panpatch report: accepted patches by category (and how many gap-fills were
    error-BED 'bed-gap' fills), plus telomere-to-telomere output contigs from the '#Contig' cap lines
    (both tips capped).  rejected candidates and passthrough (unpatched) rows are not patches """
    counts = {'gap-fill': 0, 'bed-gap': 0, 'scaffold': 0, 'telomere': 0}
    t2t, contigs = 0, 0
    idx = None
    with open(report_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#Contig'):                 # per-output-contig telomere cap status
                contigs += 1
                if 'left=YES' in line and 'right=YES' in line:
                    t2t += 1
                continue
            if not line or line.startswith('#'):
                continue
            if idx is None:                                # the header (first non-comment line)
                idx = {c: i for i, c in enumerate(line.split('\t'))}
                continue
            f2 = line.split('\t')
            d, t = idx.get('decision'), idx.get('type')
            if d is None or d >= len(f2) or f2[d] != 'accepted':
                continue
            typ = f2[t] if t is not None and t < len(f2) else ''
            if typ in ('gap-fill', 'scaffold', 'telomere'):
                counts[typ] += 1
                g = idx.get('gap_origin')
                if typ == 'gap-fill' and g is not None and g < len(f2) and f2[g] == 'bed-gap':
                    counts['bed-gap'] += 1
    return {'counts': counts, 't2t': t2t, 'contigs': contigs}

def write_batch_summary(job, options, summaries):
    """ roll the per-sample panpatch reports up into one cross-sample summary
    (outDir/panpatch-summary.tsv): a row of accepted patches by category + T2T output contigs per sample,
    then a TOTAL row.  the per-sample <name>.tsv reports are left unchanged.  summaries is a list of
    (run name, that run's output-file-id map) """
    work_dir = job.fileStore.getLocalTempDir()
    rows = []
    for name, output_id_map in summaries:
        report_id = output_id_map.get(name + '.tsv') if output_id_map else None
        if report_id is None:                              # a run that produced no report (skipped)
            continue
        p = os.path.join(work_dir, 'r.tsv')
        job.fileStore.readGlobalFile(report_id, p)
        rows.append((name, summarize_report(p)))
        os.remove(p)

    total = {'gap-fill': 0, 'bed-gap': 0, 'scaffold': 0, 'telomere': 0, 't2t': 0, 'contigs': 0}
    summary_path = os.path.join(work_dir, 'panpatch-summary.tsv')
    with open(summary_path, 'w') as out:
        out.write('# cactus-panpatch batch summary -- {} sample(s); "patches" = accepted gap_fill + '
                  'scaffold + telomere (bed_gap is the error-BED subset of gap_fill)\n'.format(len(rows)))
        out.write('sample\tgap_fill\tbed_gap\tscaffold\ttelomere\tpatches\tt2t_contigs\tcontigs_analyzed\n')
        for name, s in rows:
            c = s['counts']
            patches = c['gap-fill'] + c['scaffold'] + c['telomere']
            out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                name, c['gap-fill'], c['bed-gap'], c['scaffold'], c['telomere'], patches, s['t2t'], s['contigs']))
            for k in ('gap-fill', 'bed-gap', 'scaffold', 'telomere'):
                total[k] += c[k]
            total['t2t'] += s['t2t']
            total['contigs'] += s['contigs']
        out.write('TOTAL\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            total['gap-fill'], total['bed-gap'], total['scaffold'], total['telomere'],
            total['gap-fill'] + total['scaffold'] + total['telomere'], total['t2t'], total['contigs']))

    RealtimeLogger.info('cactus-panpatch summary: {} sample(s), {} accepted patches (gap-fill {}, '
                        'scaffold {}, telomere {}; {} bed-gap), {} T2T contigs'.format(
                            len(rows), total['gap-fill'] + total['scaffold'] + total['telomere'],
                            total['gap-fill'], total['scaffold'], total['telomere'], total['bed-gap'], total['t2t']))
    job.fileStore.exportFile(job.fileStore.writeGlobalFile(summary_path),
                             makeURL(os.path.join(options.outDir, 'panpatch-summary.tsv')))

def cleanup_pangenome_wrapper(job, options, run):
    """ the cactus-pangenome output is only ever scratch for us.  everything the user asked for has
    been exported by the time we get here, so a failure anywhere upstream leaves it alone (and
    leaves it ready for --restart) """
    pg_dir = pangenome_dir(options, run)
    RealtimeLogger.info('Removing intermediate cactus-pangenome output: {}'.format(pg_dir))
    delete_directory(pg_dir)

if __name__ == "__main__":
    main()
