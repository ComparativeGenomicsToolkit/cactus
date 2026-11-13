#!/usr/bin/env python3

"""
end-to-end assembly patching pipeline, which chains:
cactus-pangenome
panpatch

panpatch (https://github.com/glennhickey/panpatch) uses a pangenome graph to patch a (slightly)
fragmented assembly into T2T chromosomes: it fills N-gaps, scaffolds disconnected contigs,
completes missing terminal telomeres (--requireTelomeres makes them mandatory) and replaces regions
flagged as assembly errors (--assemblyErrorBeds), taking the patch sequence from one or more donor
assemblies.

the pangenome is built with --chrom-vg full and nothing else, since the per-chromosome "full" vg
graphs are the only thing panpatch needs.  cactus-pangenome is the overwhelming majority of the
compute here: panpatch itself takes a couple of minutes on a whole human pangenome.  with
--keepGraphs its input graphs are left in <outDir>/<name>.chroms/, and --fromGraphs patches those
again without rebuilding the pangenome: the assemblies are aligned unmasked and their error BEDs
(--assemblyErrorBeds) only go to panpatch, so other BEDs or panpatch parameters can be tried in minutes
"""

import os, sys
import re
import copy
import glob
import gzip
import shlex
import shutil
import timeit
import collections
import xml.etree.ElementTree as ET

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.cactus_prepare import human2bytesN
from cactus.shared.common import importSingularityImage, cactus_walltime, add_cactus_toil_options
from cactus.shared.common import makeURL
from cactus.shared.common import normalize_pangenome_name, pangenome_name_sample
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
from cactus.refmap.cactus_graphmap_join import graphmap_join_defaults, graphmap_join_resolve_references
from cactus.refmap.cactus_pangenome import pangenome_options, pangenome_validate_options
from cactus.refmap.cactus_pangenome import pangenome_config_overrides, pangenome_end_to_end_workflow
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger, set_logging_from_options
from toil.realtimeLogger import RealtimeLogger

def main():
    parser = Job.Runner.getDefaultArgumentParser()
    add_cactus_toil_options(parser)

    parser.add_argument("seqFile", nargs='?', default=None,
                        help = "Seq file (as with cactus-pangenome), or, with --batch, a chromfile with one "
                        "\"<name> <seqfile>\" line per sample")
    parser.add_argument("--outDir", help = "Output directory", required=True)
    parser.add_argument("--outName", help = "Output name (without extension) [default: the name of the sample being patched]")
    parser.add_argument("--reference", nargs='+', type=str,
                        help = "Reference event name(s), as in cactus-pangenome. Must be haploid (no #1/#2 haplotype). "
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
                        help = "Maximum bp of target sequence that a telomere patch, or a tip graft (a donor carrying the "
                        "patch past the target's end), may replace at a contig end (panpatch -M) [default=20000]")
    parser.add_argument("--excludeBed", default=None,
                        help = "BED file of target regions to leave untouched, named by full graph path (panpatch -b)")
    parser.add_argument("--noPatchFasta", default=None,
                        help = "FASTA of more sequences no graft may replace or put in (panpatch --no-patch-fasta): high-copy "
                        "sequence whose copies share their k-mers, so a swap of it cannot be verified.  panpatch already "
                        "guards the human rDNA unit (--panpatchOptions --patch-rdna lets grafts take rDNA)")
    parser.add_argument("--assemblyErrorBeds", default=None,
                        help = "Manifest of per-assembly BEDs of suspected errors: one \"<seqfile-event-name> <bed-path>\" line "
                        "per assembly (a subset is fine), each BED in that assembly's own contig names and coordinates. "
                        "The assemblies are aligned unmasked and the intervals passed to panpatch (--error-bed): a donor's "
                        "are never used to patch; the target's are replaced by donor sequence where a donor differs and "
                        "can patch them, and otherwise left as they are")
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
                        "parameters can be re-tuned (via --fromGraphs) in minutes without rebuilding the pangenome")
    parser.add_argument("--fromGraphs", default=None,
                        help = "Don't build the pangenome: patch the chromosome graphs an earlier run kept with --keepGraphs, "
                        "from <DIR>/<name>.chroms/ (DIR is that run's --outDir).  Give the same seqfile or chromfile it "
                        "was built from.  The graphs don't depend on --assemblyErrorBeds, so other BEDs (or "
                        "--panpatchOptions) can be tried without rebuilding")

    # cactus-pangenome options.  note that we deliberately do *not* pull in graphmap_join_options():
    # cactus-panpatch wants nothing but the chromosome vg graphs out of the pangenome, so there is no
    # output for the user to select.  graphmap_join_defaults() fills those options in below
    pangenome_options(parser)

    # cactus-graphmap-join uses this to cap the memory of the vg jobs that we do run
    parser.add_argument("--indexMemory", type=human2bytesN,
                        help="Memory in bytes to upper-bound the per-chromosome vg and panpatch jobs. "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes))", default=None)

    # hal2vg (which converts each chromosome alignment to vg) is sized from the alignment's memory,
    # which badly under-estimates it for the large unclipped reference-free graphs here (measured up
    # to ~27Gi with hal2vg 1.1.8, less since), so floor it.  standard suffixes supported (default=bytes)
    parser.add_argument("--halExportMemory", type=human2bytesN, default=24 * 2**30,
                        help="Minimum memory in bytes for each per-chromosome hal2vg job [default=24Gi]")

    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    panpatch_validate_options(options)

    # per-assembly error BEDs (treated almost like runs of Ns): parse the manifest up front so run
    # construction (make_seqfile_runs) can attach each BED to its assembly, and so a typo'd event name
    # is caught before anything expensive runs
    options.errorBedMap = read_error_bed_manifest(options.assemblyErrorBeds) if options.assemblyErrorBeds else {}

    # we need the minigraph event name to filter it out of the seqfiles below, and how to read their names
    config_node = ET.parse(options.configFile).getroot()
    graph_event = getOptionalAttrib(findRequiredNode(config_node, "graphmap"),
                                    "assemblyName", default="_MINIGRAPH_")
    dot_haplotypes = ConfigWrapper(config_node).getDotHaplotypes()

    # work out the runs (one cactus-pangenome + panpatch each) before validating, since with
    # --referenceFree the reference is derived from the sample being patched
    runs = make_runs(options, graph_event, dot_haplotypes)
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

            # don't cut the reference's ends into tips in the full graphs.  that is only for vg haplotypes
            # (--haplo full), which we never run, and it breaks patching: telomere anchors align every
            # haplotype's array to the reference's end, which they overhang by a few bp, so the cut splits
            # the target's own paths there.  an unpatched T2T contig then comes out as two records (a few-bp
            # stub and the rest), and a tip cut off a patched contig is dropped without a report row
            findRequiredNode(config_node, "graphmap_join").attrib["fullStubbify"] = "0"

            if options.collapseRefPAF:
                findRequiredNode(config_node, "graphmap").attrib["collapse"] = 'reference'

            exclude_bed_id = None
            if options.excludeBed:
                exclude_bed_id = toil.importFile(makeURL(options.excludeBed))
            # the --noPatchFasta sequences, which every chromosome job hands to panpatch
            options.noPatchFastaID = toil.importFile(makeURL(options.noPatchFasta)) if options.noPatchFasta else None

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
                    if options.fromGraphs and sample != run['reference'][0]:
                        # no pangenome to build: only the reference is read (its size scales gather)
                        continue
                    seq_url = makeURL(seq_path)
                    seq_id_map[sample] = toil.importFile(seq_url)
                    seq_path_map[sample] = seq_url
                    seq_order.append(sample)

                # --fromGraphs: the chromosome graphs an earlier --keepGraphs run left for this one
                graph_ids = None
                if options.fromGraphs:
                    chrom_dir = os.path.join(options.fromGraphs, run['name'] + '.chroms')
                    vg_paths = sorted(glob.glob(os.path.join(chrom_dir, '*.vg')))
                    if not vg_paths:
                        raise RuntimeError('--fromGraphs: no chromosome graphs (*.vg) for {} in {}'.format(run['name'], chrom_dir))
                    graph_ids = [(toil.importFile(makeURL(p)), os.path.basename(p)) for p in vg_paths]

                ref_collapse_paf_id = toil.importFile(makeURL(options.collapseRefPAF)) if options.collapseRefPAF else None
                last_scores_id = toil.importFile(makeURL(options.scoresFile)) if options.scoresFile else None

                # a *second* import of the target fastas, kept aside for the dropped-contig rescue in
                # run_panpatch.  it can't share the pangenome's copy above: that one is consumed
                # (deleted) by sanitize_fasta_header during construction
                target_fasta_ids = {hap: toil.importFile(makeURL(path)) for hap, path in run['target_haps']}

                # per-assembly error BEDs: import each path once, then key it by seq_rows sample name
                # (for panpatch's --error-bed, named as the graph names its paths) and by graph haplotype
                # for the target (for report tagging in run_panpatch_chrom)
                error_bed_id_cache = {}
                for bed_path in set(run['error_beds'].values()) | set(run['target_error_beds'].values()):
                    error_bed_id_cache[bed_path] = toil.importFile(makeURL(bed_path))
                error_bed_ids = {sample: error_bed_id_cache[p] for sample, p in run['error_beds'].items()}
                target_error_bed_ids = {hap: error_bed_id_cache[p] for hap, p in run['target_error_beds'].items()}

                pg_options = copy.deepcopy(options)
                pg_options.outDir = pangenome_dir(options, run)
                pg_options.outName = run['name']
                pg_options.reference = list(run['reference'])
                graphmap_join_resolve_references(pg_options, dot_haplotypes, references_resolved=True)
                pg_options.seqFile = run['seqFile']
                # we only want the chromosome vgs: don't make a hal we'd throw away, and don't
                # export the vgs into the pangenome directory when we're about to export them
                # ourselves (they're the biggest output there is)
                pg_options.noHal = True
                pg_options.noJoinExport = not options.keepPangenome
                if not pg_options.outDir.startswith('s3://') and not options.fromGraphs:
                    os.makedirs(pg_options.outDir, exist_ok=True)

                run_inputs.append((run, pg_options, seq_id_map, seq_path_map, seq_order,
                                   ref_collapse_paf_id, last_scores_id, target_fasta_ids,
                                   error_bed_ids, target_error_bed_ids, graph_ids))

            toil.start(Job.wrapJobFn(panpatch_batch_workflow, options, config_wrapper, run_inputs,
                                     exclude_bed_id, walltime=cactus_walltime()))

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-panpatch has finished after {} seconds".format(run_time))

def sample_base(sample):
    """ strip the #N haplotype off a genome name to get the (PanSN) sample name it has in the graph,
    which is the name panpatch wants.  a sample name is allowed to contain dots of its own (ex
    HG002.verkko#1 is haplotype 1 of HG002.verkko), and the names given to --sample / --defaultSample
    are sample names already.  the old SAMPLE.N spelling of a haplotype is only read in the seqfile
    (see normalize_pangenome_name), so a '.' is never a haplotype here """
    return pangenome_name_sample(sample)

def sample_hap(sample):
    """ the haplotype number of a genome name (as normalized from the seqfile), or None if it
    doesn't give one (haploid) """
    hap = sample.partition('#')[2]
    return int(hap) if hap else None

def sanitize_contig_name(header):
    """ mirror cactus_sanitizeFastaHeaders -p on a single fasta header, so an input contig name
    matches the CONTIG field of the graph path names (SAMPLE#HAP#CONTIG).  cut at whitespace, drop
    an id=EVENT| prefix (which cactus_consolidated strips, and whose EVENT can contain a '#'), strip
    through the last '#', convert a trailing ':start-end' range to '_sub_(start-1)_end', and replace
    any other ':' with '_' (kept consistent with cactus_sanitizeFastaHeaders.c) """
    name = header.split()[0]
    if name.startswith('id=') and '|' in name:
        name = name[name.find('|') + 1:]
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
    # switches to external-reference mode (whose reference make_seqfile_runs checks is haploid, once
    # it knows which genome the name is)
    options.referenceFree = not options.reference

    if options.noSplit:
        # a single unsplit graph has every reference chromosome in it, and panpatch skips any graph
        # that doesn't have exactly one reference path -- it would run to completion and hand back
        # the unpatched input
        raise RuntimeError('--noSplit cannot be used with cactus-panpatch: panpatch needs one graph per reference chromosome')

    # getattr because the unit tests build a minimal namespace rather than going through the parser
    if getattr(options, 'fromGraphs', None):
        if getattr(options, 'keepGraphs', False) or getattr(options, 'keepPangenome', False):
            raise RuntimeError('--keepGraphs / --keepPangenome cannot be used with --fromGraphs: no pangenome is built, '
                               'and its graphs are already kept')
        if options.fromGraphs.startswith('s3://') or not os.path.isdir(options.fromGraphs):
            raise RuntimeError('--fromGraphs {} is not a local directory'.format(options.fromGraphs))
        out_dir = getattr(options, 'outDir', None)
        if out_dir and not out_dir.startswith('s3://') and os.path.realpath(out_dir) == os.path.realpath(options.fromGraphs):
            raise RuntimeError('--outDir cannot be the --fromGraphs directory {}: its run\'s outputs would be '
                               'overwritten'.format(options.fromGraphs))
    if getattr(options, 'inGFA', None) or getattr(options, 'inGAF', None) or getattr(options, 'remap', False):
        # these come in via pangenome_options(), which panpatch shares.  the graph panpatch builds
        # is a throwaway, built per sample out of that sample and its donors, so there is nothing
        # an earlier run could usefully be extended from
        raise RuntimeError('--inGFA / --inGAF / --remap cannot be used with cactus-panpatch: it builds a fresh '
                           'graph per sample being patched, which is not an extension of anything')

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

def read_seqfile(seqfile_path, graph_event, dot_haplotypes):
    """ the (sample, fasta) rows of a seqfile, in order.  the minigraph event is skipped: a seqfile
    that cactus-pangenome exported from a previous run will have it """
    seq_file = SeqFile(seqfile_path, pangenome=True, dot_haplotypes=dot_haplotypes)
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

def make_runs(options, graph_event, dot_haplotypes):
    """ work out the runs implied by the options: one cactus-pangenome + panpatch each.  a normal
    invocation is a single run, --referenceFree is one run per haplotype of the sample being patched,
    and --batch multiplies all that by the number of samples.  dot_haplotypes is the config's
    <graphmap dotHaplotypes>, for reading the seqfiles """
    if options.panpatchBatch:
        # with --batch the seqFile argument is the chromfile (as in cactus-align --batch)
        seqfiles = read_batch_file(options.seqFile)
    else:
        seqfiles = [(options.outName, options.seqFile)]

    # the --assemblyErrorBeds events are named as in the seqfile, so they're normalized the same way
    # (TGT.1 is TGT#1 with dotHaplotypes).  given_events maps them back to how they were given, to
    # report them
    given_events = {}
    for event in options.errorBedMap:
        norm_event = normalize_pangenome_name(event, dot_haplotypes)
        if norm_event in given_events and given_events[norm_event] != event:
            raise RuntimeError('--assemblyErrorBeds has more than one entry for {}: {} and {}'.format(
                norm_event, given_events[norm_event], event))
        given_events[norm_event] = event
    options.errorBedMap = {normalize_pangenome_name(event, dot_haplotypes) : bed for event, bed in options.errorBedMap.items()}

    runs = []
    for name, seqfile_path in seqfiles:
        runs += make_seqfile_runs(options, name, seqfile_path, graph_event, dot_haplotypes)

    names = [run['name'] for run in runs]
    if len(set(names)) != len(names):
        raise RuntimeError('Output names are not unique: {}'.format(names))
    for run_name in names:
        # the names get pasted into output paths (and into a directory we later delete), so don't
        # let them wander out of --outDir.  nor can they have a '#', which Toil can't export to S3
        if os.path.basename(run_name) != run_name or run_name.startswith('.') or '#' in run_name:
            raise RuntimeError('Invalid output name "{}": it cannot contain "/" or "#" or begin with "."'.format(run_name))

    # every event named in --assemblyErrorBeds must have matched an assembly in some run, else it's a
    # typo (or a reference, which is never patched) -- catch it before anything expensive runs.  compare
    # by the consumed event names (run['error_events']), not by bed path: two events may share one bed
    # file, and the reference-free target's error_beds key is its graph_name, not its manifest event name
    if options.errorBedMap:
        matched = set()
        for run in runs:
            matched.update(run['error_events'])
        unmatched = sorted(given_events[event] for event in set(options.errorBedMap) - matched)
        if unmatched:
            raise RuntimeError('--assemblyErrorBeds event(s) did not match any patchable assembly (the target or a '
                               'donor) in any seqfile -- check the names (the reference is never patched): {}'.format(
                                   ', '.join(unmatched)))
    return runs

def make_seqfile_runs(options, name, seqfile_path, graph_event, dot_haplotypes):
    """ the runs implied by a single seqfile """
    rows = read_seqfile(seqfile_path, graph_event, dot_haplotypes)
    if not rows:
        raise RuntimeError('No sequences found in {}'.format(seqfile_path))
    options.reference = check_sample_names([row[0] for row in rows], options.reference, dot_haplotypes)
    if options.reference:
        # a Minigraph-Cactus reference must be haploid: a #1/#2/... haplotype means it is one haplotype
        # of a diploid sample, which cannot be the reference (use reference-free for that)
        for ref in options.reference:
            if sample_hap(ref):
                raise RuntimeError('--reference {} must be haploid: it cannot be haplotype {} of {} '
                                   '(use no haplotype, or #0)'.format(ref, sample_hap(ref), sample_base(ref)))

    # the sample names as they appear in the graph, in seqfile order
    samples = []
    for sample, _ in rows:
        if sample_base(sample) not in samples:
            samples.append(sample_base(sample))

    # the reference is never inferred from the seqfile: it's either given with --reference, or (with
    # --referenceFree) it *is* the sample being patched, which is named just below.  note that cactus
    # takes the reference as it's named in the seqfile (which can carry a #N haplotype),
    # whereas panpatch wants the sample name it has in the graph, which never does
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
        # PanSN haplotype of the target path: the seqfile's #N, or 0 for a haploid target
        # that gives none -- which is the number panpatch actually writes (.hap0)
        target_haps = [(sample_hap(row[0]) if sample_hap(row[0]) is not None else 0, row[1])
                       for row in target_rows]
        # attach per-assembly error BEDs (see make_seqfile_runs's --referenceFree branch for the keying
        # rationale).  in reference mode the seq_rows keys are the original event names, so error_beds is
        # keyed by those; target_error_beds is keyed by graph haplotype for report tagging
        errmap = options.errorBedMap
        error_beds = {row[0]: errmap[row[0]] for row in target_rows + donor_rows if row[0] in errmap}
        error_events = {row[0] for row in target_rows + donor_rows if row[0] in errmap}
        target_error_beds = {(sample_hap(row[0]) if sample_hap(row[0]) is not None else 0): errmap[row[0]]
                             for row in target_rows if row[0] in errmap}
        return [{'name' : name if name else target,
                 'seqFile' : seqfile_path,
                 'seq_rows' : seq_rows,
                 'reference' : options.reference,
                 'panpatchReference' : reference,
                 'samples' : [target] + donors,
                 'target_sample' : target,
                 'defaultSample' : resolve_default_sample(options, target, donors, target),
                 'target_haps' : target_haps,
                 'error_beds' : error_beds,
                 'error_events' : error_events,
                 'target_error_beds' : target_error_beds,
                 'ploidy' : len(target_rows),
                 'referenceFree' : False,
                 'hap' : None}]

    # --referenceFree: a Minigraph-Cactus reference must be haploid, so we can't just point
    # --reference at TARGET#1 (haplotype 1 of a diploid TARGET). we rename it to the haplotype-free
    # TARGET_1 instead, which makes it a haploid sample in its own right, and
    # build a separate graph -- and run panpatch -- for each haplotype of the target. the other
    # haplotypes of the target (and anything else in the seqfile) are left out of the graph
    runs = []
    errmap = options.errorBedMap
    for target_row in target_rows:
        hap = sample_hap(target_row[0])
        # SAMPLE#1 -> SAMPLE_1.  the dots of a sample name are flattened as well (HG002.verkko#1 ->
        # HG002_verkko_1): that was needed when a '.' marked the haplotype, and it keeps the names
        # the same as they were
        graph_name = target.replace('.', '_') if hap is None else '{}_{}'.format(target.replace('.', '_'), hap)
        # per-assembly error BEDs.  the target is renamed to graph_name in seq_rows so its error BED is
        # keyed by graph_name; donors keep their original names.  the reference-free target is graph
        # haplotype 0, so its report-tagging BED is keyed by 0
        error_beds = {dr[0]: errmap[dr[0]] for dr in donor_rows if dr[0] in errmap}
        error_events = {dr[0] for dr in donor_rows if dr[0] in errmap}
        target_error_beds = {}
        if target_row[0] in errmap:
            error_beds[graph_name] = errmap[target_row[0]]
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
                     'target_sample' : target,
                     'defaultSample' : resolve_default_sample(options, target, donors, graph_name),
                     'target_haps' : [(0, target_row[1])],
                     'error_beds' : error_beds,
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

# Walltime estimates for the jobs below.  cactus-panpatch is newer than any run we have logs for,
# so neither of these two rates was measured on it: each is the slow end of the nearest thing that
# was measured, and each covers one job per sample or one per chromosome graph.

# Bytes of chromosome graph per second for one single-threaded panpatch job.  panpatch itself takes
# a couple of minutes on a whole human pangenome (see the module docstring); what dominates is
# loading the graph, and it gets loaded twice -- once by the `vg paths -L -x` that lists the samples
# in it, once by panpatch.  20 MB/s is the slow end of the vg graph loads measured on the HPRC v2.1
# chromosome graphs (`vg paths` p90 165s, `vg paths|awk` p90 471s, on graphs of ~24 GiB).
PANPATCH_BYTES_PER_SEC = 20e6

# Bytes of output assembly per second for the gather pass: bgzip of each output haplotype (single-
# threaded it runs 21-27 MB/s, measured on the HPRC PAFs -- see GZIP_COMPRESS_BYTES_PER_SEC in
# cactus.shared.common), plus whole-genome python passes over it (renaming and the dropped-contig rescue)
# plus panpatch's telomere report over each finished fasta.  the rate is still
# the single-threaded one, so the walltime it gives is an upper bound
GATHER_BYTES_PER_SEC = 8e6

# Cores for the gather pass.  only its bgzip uses more than one (bgzip --threads job.cores): one
# thread took ~7 minutes per haplotype
GATHER_CORES = 4

def panpatch_batch_workflow(job, options, config_wrapper, run_inputs, exclude_bed_id):
    """ patch each sample.  the runs are completely independent of each other: they share no input
    files in the jobstore (see the import loop in main()) """
    summaries = []
    for run, pg_options, seq_id_map, seq_path_map, seq_order, ref_collapse_paf_id, last_scores_id, target_fasta_ids, error_bed_ids, target_error_bed_ids, graph_ids in run_inputs:
        run_job = job.addChildJobFn(panpatch_run_workflow, options, pg_options, config_wrapper, run, seq_id_map,
                                    seq_path_map, seq_order, ref_collapse_paf_id, last_scores_id, exclude_bed_id,
                                    target_fasta_ids, error_bed_ids, target_error_bed_ids, graph_ids,
                                    walltime=cactus_walltime())
        summaries.append((run['name'], run_job.rv()))
    # once every sample is done, roll the per-sample reports up into one cross-sample summary
    job.addFollowOnJobFn(write_batch_summary, options, summaries, walltime=cactus_walltime())

def make_panpatch_error_bed(job, error_bed_ids, label=None):
    """ write the error intervals of every assembly in a run into one BED for panpatch --error-bed, each
    contig named as the graph names its path: SAMPLE#HAP#CONTIG, with the sample and haplotype of the
    assembly's seqfile name (a haploid one, such as the reference-free target, is haplotype 0) and the
    contig name as cactus sanitizes the fasta header.  error_bed_ids maps seq_rows sample name -> BED.

    the assemblies go into the pangenome unmasked, so where a donor agrees with the target it shares the
    target's nodes.  panpatch then blocks each path only at its own intervals: the target's flagged
    sequence is replaced where a donor differs and can patch it, and a donor's never goes into a patch """
    work_dir = job.fileStore.getLocalTempDir()
    out_path = os.path.join(work_dir, 'errors.bed')
    n_ivs = 0
    with open(out_path, 'w') as out:
        for i, sample in enumerate(sorted(error_bed_ids)):
            bed_path = os.path.join(work_dir, 'in.{}.bed'.format(i))
            job.fileStore.readGlobalFile(error_bed_ids[sample], bed_path)
            hap = sample_hap(sample)
            prefix = '{}#{}#'.format(sample_base(sample), hap if hap is not None else 0)
            for contig, ivs in parse_error_bed(bed_path, sanitize=True).items():
                for s, e in ivs:
                    out.write('{}{}\t{}\t{}\n'.format(prefix, contig, s, e))
                    n_ivs += 1
    RealtimeLogger.info('{}: {:,} error-BED intervals on {} assemblies for panpatch --error-bed'.format(
        label or 'make_panpatch_error_bed', n_ivs, len(error_bed_ids)))
    return job.fileStore.writeGlobalFile(out_path)

def panpatch_run_workflow(job, options, pg_options, config_wrapper, run, seq_id_map, seq_path_map, seq_order,
                          ref_collapse_paf_id, last_scores_id, exclude_bed_id, target_fasta_ids,
                          error_bed_ids, target_error_bed_ids, graph_ids=None):
    """ cactus-pangenome, then panpatch on the chromosome graphs it made (or, with --fromGraphs, on
    graph_ids: the [(id, name)] graphs an earlier run kept) """
    # the assemblies are aligned as they are: their error intervals go to panpatch (--error-bed), not into
    # the pangenome's inputs as N's.  masking cut a flagged region out of the alignment, so a donor that
    # agreed with the target there no longer shared its nodes, every flagged region became a gap to fill,
    # and the anchors either side of a donor's own masked copy could land on the wrong ends of it
    error_bed_id = None
    if error_bed_ids:
        error_bed_id = job.addChildJobFn(make_panpatch_error_bed, error_bed_ids, run['name'],
                                         disk=sum(b.size for b in error_bed_ids.values()) * 3 + 2**28,
                                         walltime=cactus_walltime()).rv()

    if graph_ids is not None:
        # --fromGraphs: nothing to build.  gather is sized by the reference (the sanitized input's size on
        # the pangenome path); a gzipped fasta is ~4x smaller than the sequence it holds
        ref = run['reference'][0]
        ref_size = seq_id_map[ref].size * (4 if seq_path_map[ref].endswith('.gz') else 1)
        return job.addFollowOnJobFn(panpatch_graphs_workflow, options, run, [g[0] for g in graph_ids],
                                    [g[1] for g in graph_ids], ref_size, exclude_bed_id, target_fasta_ids,
                                    target_error_bed_ids, error_bed_id, walltime=cactus_walltime()).rv()

    # a follow-on of the job hosting the pangenome runs once its entire child subtree is done, which is
    # what makes panpatch_workflow safe (and resolves the error BED's promise)
    pangenome_job = job.addFollowOnJobFn(pangenome_end_to_end_workflow, pg_options, config_wrapper, seq_id_map,
                                         seq_path_map, seq_order, ref_collapse_paf_id, last_scores_id,
                                         walltime=cactus_walltime())
    # return the run's output-file map (its report id feeds the cross-sample batch summary)
    return pangenome_job.addFollowOnJobFn(panpatch_workflow, options, run, pangenome_job.rv(0), pangenome_job.rv(1),
                                          pangenome_job.rv(2), exclude_bed_id, target_fasta_ids,
                                          target_error_bed_ids, error_bed_id, walltime=cactus_walltime()).rv()

def panpatch_workflow(job, options, run, join_options, join_wf_output, seq_id_map, exclude_bed_id, target_fasta_ids,
                      target_error_bed_ids, error_bed_id=None):
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

    # seq_id_map holds the sanitized (uncompressed) inputs, so the reference's size is the genome size
    return add_panpatch_jobs(job, options, run, full_vg_ids, vg_names, seq_id_map[run['reference'][0]].size,
                             exclude_bed_id, target_fasta_ids, target_error_bed_ids, error_bed_id)

def panpatch_graphs_workflow(job, options, run, full_vg_ids, vg_names, ref_size, exclude_bed_id, target_fasta_ids,
                             target_error_bed_ids, error_bed_id=None):
    """ --fromGraphs: patch the chromosome graphs an earlier --keepGraphs run left, as panpatch_workflow
    patches the ones the pangenome just made """
    return add_panpatch_jobs(job, options, run, full_vg_ids, vg_names, ref_size, exclude_bed_id, target_fasta_ids,
                             target_error_bed_ids, error_bed_id)

def add_panpatch_jobs(job, options, run, full_vg_ids, vg_names, ref_size, exclude_bed_id, target_fasta_ids,
                      target_error_bed_ids, error_bed_id):
    """ add the per-chromosome panpatch jobs as children of job, and the gather and export as its
    follow-ons.  ref_size is the genome size gather is scaled by.  returns gather's output-file map """
    # one single-threaded panpatch job per chromosome graph.  each needs memory for just its own
    # graph, and --indexMemory (if given) caps the estimate
    chrom_jobs = []
    for vg_id, vg_name in zip(full_vg_ids, vg_names):
        mem = cactus_clamp_memory(max(2**32, vg_id.size * 8))
        if options.indexMemory:
            mem = min(mem, options.indexMemory)
        cj = job.addChildJobFn(run_panpatch_chrom, options, run, vg_id, vg_name, exclude_bed_id,
                               target_error_bed_ids, error_bed_id,
                               cores=1, memory=mem, disk=vg_id.size * 4 + 2**30,
                               walltime=cactus_walltime(vg_id.size / PANPATCH_BYTES_PER_SEC,
                                                        io_bytes=vg_id.size))
        chrom_jobs.append(cj)

    # gather the per-chromosome outputs, rescue dropped contigs, and bgzip.  the concatenated fastas
    # are ~ploidy * genome
    gather_disk = int(run['ploidy'] * ref_size * 3) + sum(fid.size for fid in target_fasta_ids.values()) + 2**30
    gather_mem = ref_size * 2
    gather_job = job.addFollowOnJobFn(gather_panpatch, options, run, [cj.rv() for cj in chrom_jobs],
                                      target_fasta_ids, target_error_bed_ids,
                                      cores=min(GATHER_CORES, int(options.maxCores)) if options.maxCores else GATHER_CORES,
                                      memory=cactus_clamp_memory(max(2**32, gather_mem)),
                                      disk=gather_disk,
                                      walltime=cactus_walltime(run['ploidy'] * ref_size / GATHER_BYTES_PER_SEC,
                                                               io_bytes=gather_disk))
    # this copies the patched assemblies out, and with --keepGraphs every chromosome graph too,
    # so it is bounded by the same bytes its disk request is
    export_bytes = sum(vg_id.size for vg_id in full_vg_ids) * 2 + 2**30
    gather_job.addFollowOnJobFn(export_panpatch_wrapper, options, run, gather_job.rv(), full_vg_ids, vg_names,
                                disk=export_bytes,
                                walltime=cactus_walltime(0, io_bytes=export_bytes))
    # surface this run's output-file map (incl. its report id) up to the cross-sample batch summary
    return gather_job.rv()

def tag_report_gap_origin(job, report_path, target_error_bed_ids):
    """ add a 'gap_origin' column to panpatch's per-chromosome report, distinguishing a genuine N-gap
    fill ('N-gap') from a fill of a target error interval ('bed-gap'); '.' for rows with no gap origin.
    panpatch reports both as gap-fill rows: a gap-fill OR telomere row is 'bed-gap' iff its replaced
    target span overlaps a target error interval, which is how a telomere completed over a flagged tip
    becomes self-evident in the report.  a gap-fill with no overlap is a genuine 'N-gap'; a telomere with
    no overlap capped an unflagged end and stays '.'.  the column is appended so existing column positions
    are unchanged """
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
        rtype = fields[idx['type']] if len(fields) == len(cols) else None
        if rtype in ('gap-fill', 'telomere'):
            # both a gap-fill and a telomere completion replace a target span; the span is 'bed-gap' iff it
            # overlaps a target error interval.  a gap-fill not over one sits on an N gap, so it is a
            # genuine 'N-gap'; a telomere completion may instead cap an unflagged (truncated/absent-
            # telomere) end, which is no gap at all, so with no overlap it stays '.'
            origin = 'N-gap' if rtype == 'gap-fill' else '.'
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

def run_panpatch_chrom(job, options, run, vg_id, vg_name, exclude_bed_id, target_error_bed_ids, error_bed_id=None):
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
    if error_bed_id:
        error_bed_path = os.path.join(work_dir, 'errors.bed')
        job.fileStore.readGlobalFile(error_bed_id, error_bed_path)
        cmd += ['--error-bed', error_bed_path]
    if getattr(options, 'noPatchFastaID', None):
        no_patch_path = os.path.join(work_dir, 'no_patch.fa')
        job.fileStore.readGlobalFile(options.noPatchFastaID, no_patch_path)
        cmd += ['--no-patch-fasta', no_patch_path]
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

PANPATCH_REPORT_COLUMNS = ('chrom\thap\ttype\ttarget\ttarget_bp\tdonor\tdonor_bp\treplaced_bp\tkmer%\tflankL%\t'
                           'flankR%\tdecision\treason\ttarget_start\ttarget_end')

def run_renames(run):
    """ {graph sample: (sample, haplotype)} for the names panpatch's graph gives a run's target that the
    user never wrote: a reference-free target is the haploid TARGET_1 in its graph, which is TARGET#1 """
    if not run.get('referenceFree'):
        return {}
    return {run['samples'][0]: (run['target_sample'], run['hap'] if run['hap'] is not None else '0')}

def clean_path_name(path, renames):
    """ a graph path name (SAMPLE#HAP#CONTIG, with the trailing #0 phase block of a full graph and maybe a
    [START-END] subrange) as the user's SAMPLE#HAP#CONTIG: the phase block dropped, a renamed reference-free
    target mapped back to its own name.  anything else is returned as it is """
    m = re.match(r'^(.*?)(\[\d+(?:-\d+)?\])?$', path)
    base, sub = m.group(1), m.group(2) or ''
    fields = base.split('#')
    if len(fields) < 3 or not fields[1].isdigit():
        return path
    sample, hap, contig = fields[0], fields[1], fields[2]
    if sample in renames:
        sample, hap = renames[sample]
    return '{}#{}#{}{}'.format(sample, hap, contig, sub)

def clean_names_in_text(text, renames):
    """ clean_path_name applied to every path name in a free-text field (a report row's reason) """
    return ' '.join(clean_path_name(t, renames) if t.count('#') >= 2 else t for t in text.split(' '))

def clean_report_names(report_path, renames, hap=None):
    """ rewrite a run's report in the user's names: the target and donor columns and any path in a
    reason, and (reference-free) the hap column, which panpatch gives as the haploid graph sample's 0 """
    with open(report_path) as f:
        lines = f.readlines()
    idx = None
    out = []
    for line in lines:
        s = line.rstrip('\n')
        if not s or s.startswith('#'):
            out.append(line)
            continue
        if idx is None:
            idx = {c: i for i, c in enumerate(s.split('\t'))}
            out.append(line)
            continue
        f = s.split('\t')
        for col in ('target', 'donor'):
            if col in idx and idx[col] < len(f):
                f[idx[col]] = clean_path_name(f[idx[col]], renames)
        if 'reason' in idx and idx['reason'] < len(f):
            f[idx['reason']] = clean_names_in_text(f[idx['reason']], renames)
        if hap is not None and 'hap' in idx and idx['hap'] < len(f):
            f[idx['hap']] = str(hap)
        out.append('\t'.join(f) + '\n')
    with open(report_path, 'w') as fo:
        fo.writelines(out)

def clean_bed_names(bed_path, renames):
    """ rewrite a run's bed in the user's names: each block header's SAMPLE#HAP and each interval's path """
    with open(bed_path) as f:
        lines = f.readlines()
    out = []
    for line in lines:
        s = line.rstrip('\n')
        m = re.match(r'^(#.* for )(.*)(:)$', s)
        if m:
            sample, _, hap = m.group(2).rpartition('#')
            if sample in renames:
                sample, hap = renames[sample]
            out.append('{}{}#{}{}\n'.format(m.group(1), sample, hap, m.group(3)))
        elif s and not s.startswith('#'):
            f = s.split('\t')
            f[0] = clean_path_name(f[0], renames)
            out.append('\t'.join(f) + '\n')
        else:
            out.append(line)
    with open(bed_path, 'w') as fo:
        fo.writelines(out)

def output_record_name(name, reference_free):
    """ the name a record of the patched assembly is delivered under.  a record that is one input contig
    (patched reference-free, or written as it is) takes that contig's name, as cactus read it from the
    fasta header; a reference-based patched record (one chromosome, maybe scaffolded from several contigs)
    keeps panpatch's <chrom>_hap_<N> """
    if '#' in name:
        return name.split('[')[0].split('#')[2] if name.count('#') >= 2 else name
    if reference_free and '_hap_' in name:
        return name[:name.rfind('_hap_')]
    return name

def rename_fasta_records(fasta_path, reference_free):
    """ rename every record of a patched haplotype (see output_record_name), in place.  returns the new
    names; two records ending up with one name is an error """
    out_path = fasta_path + '.renamed'
    names = []
    with open(fasta_path) as inp, open(out_path, 'w') as out:
        for line in inp:
            if line.startswith('>'):
                names.append(output_record_name(line[1:].split()[0], reference_free))
                out.write('>' + names[-1] + '\n')
            else:
                out.write(line)
    dups = sorted(n for n, c in collections.Counter(names).items() if c > 1)
    if dups:
        raise RuntimeError('cactus-panpatch: output record name(s) {} would be used more than once in {}'.format(
            ', '.join(dups[:10]), os.path.basename(fasta_path)))
    os.replace(out_path, fasta_path)
    return names

def fasta_contig_names(fasta_path):
    """ the contig names of a (maybe gzipped) fasta, sanitized as cactus sanitizes them for the graph """
    with open(fasta_path, 'rb') as fh:
        opener = gzip.open if fh.read(2) == b'\x1f\x8b' else open
    names = []
    with opener(fasta_path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                names.append(sanitize_contig_name(line[1:]))
    return names

def missing_target_contigs(input_contigs, record_names, bed_path, report_path, target_sample, hap):
    """ the input target contigs (graph haplotype hap) that a patched haplotype does not account for.  a
    contig is accounted for when it is written as a record of its own (SAMPLE#HAP#CONTIG), is a target
    interval of a patched record in the bed, or is named by an accepted patch (a fragment a scaffold or a
    telomere graft took over).  anything else would be silently lost """
    accounted = set()
    def take(path):
        fields = path.split('[')[0].split('#')
        if len(fields) >= 3 and fields[0] == target_sample and fields[1] == str(hap):
            accounted.add(fields[2])
    for name in record_names:
        take(name)
    if bed_path and os.path.exists(bed_path):
        for (locus, block_hap), segs in parse_bed_blocks(bed_path).items():
            if block_hap == str(hap):
                for seg in segs:
                    take(seg[0])
    if report_path and os.path.exists(report_path):
        idx = None
        with open(report_path) as f:
            for line in f:
                s = line.rstrip('\n')
                if not s or s.startswith('#'):
                    continue
                if idx is None:
                    idx = {c: i for i, c in enumerate(s.split('\t'))}
                    continue
                f2 = s.split('\t')
                if f2[idx['decision']] != 'accepted' or f2[idx['hap']] != str(hap):
                    continue
                take(f2[idx['target']])
                take(f2[idx['donor']])
    return [c for c in input_contigs if c not in accounted]

def gather_panpatch(job, options, run, chrom_results, target_fasta_ids, target_error_bed_ids):
    """ concatenate the per-chromosome panpatch outputs (in chromosome order, as panpatch itself would),
    rescue any target contig the pangenome dropped, bgzip, and return the output file map """
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

    # panpatch's telomere threshold, if the user overrode it, so our post-revert measurement below calls
    # a telomere exactly the way the run that produced these patches did
    telo_opts = []
    if options.panpatchOptions:
        popts = shlex.split(options.panpatchOptions)
        for i, o in enumerate(popts):
            if o == '--telomere-threshold' and i + 1 < len(popts):
                telo_opts = ['--telomere-threshold', popts[i + 1]]
    cap_lines = []   # "#Contig" cap status of the assembly we actually deliver (filled per haplotype)

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

    # and the reports: the TSV header from the first chromosome that has one, then every chromosome's
    # rows (a skipped graph has no report).  the completeness check below reads this one; the delivered
    # report is rewritten at the end with the measured cap lines and the user's names
    report_out = os.path.join(work_dir, 'patched.tsv')
    chrom_report_paths = []
    for i, res in enumerate(chrom_results):
        if res['report'] is None:
            continue
        p = os.path.join(work_dir, 'c{}.tsv'.format(i))
        job.fileStore.readGlobalFile(res['report'], p)
        chrom_report_paths.append(p)
    concat_reports(chrom_report_paths, [], report_out)

    def t2t_count(lines):
        return sum(1 for l in lines if 'left=YES' in l and 'right=YES' in l)

    summary_lines = []   # "#Summary" lines, one per haplotype, for the batch summary
    rescue_rows = []     # passthrough rows for the contigs the rescue writes
    for hap in haps:
        # the user's haplotype: reference-free runs are the haploid graph sample's 0
        user_hap = (run['hap'] if run['hap'] is not None else '0') if run['referenceFree'] else str(hap)
        label = '{}{}'.format(run['name'], '' if run['referenceFree'] else ' hap{}'.format(hap))
        combined = os.path.join(work_dir, 'hap{}.fa'.format(hap))
        with open(combined, 'wb') as out:
            for res in chrom_results:   # chromosome order
                if hap in res['fastas']:
                    chrom_fa = os.path.join(work_dir, 'chrom.fa')
                    job.fileStore.readGlobalFile(res['fastas'][hap], chrom_fa)
                    with open(chrom_fa, 'rb') as f:
                        shutil.copyfileobj(f, out)
                    os.remove(chrom_fa)

        # rescue any target contig the pangenome dropped, appending it verbatim from the input assembly, and
        # describe it in the report and the bed like any other record written as it is
        if hap in target_fasta_ids:
            rescued = rescue_dropped_contigs(job, combined, target_fasta_ids[hap], seen_contigs,
                                             run['samples'][0], hap, label=label)
            with open(bed_out, 'a') as bed:
                for contig, length in rescued:
                    path = '{}#{}#{}'.format(run['samples'][0], hap, contig)
                    rescue_rows.append('\t'.join(['.', str(hap), 'passthrough', path, str(length), '.', '0', '.', '.',
                                                  '.', '.', 'accepted', 'not in any pangenome graph; written as it is',
                                                  '.', '.', '.']))
                    bed.write('#Passthrough (not in any graph) for {}#{}:\n{}\t0\t{}\t+\n\n'.format(
                        run['samples'][0], hap, path, length))

        # every input contig must be accounted for: written as a record of its own, a target interval of a
        # patched record, or taken over by an accepted patch.  anything else is a contig the output would
        # silently lack, so stop here rather than deliver it
        input_contigs, input_caps = [], []
        if hap in target_fasta_ids:
            target_in = os.path.join(work_dir, 'target_input.fa')
            job.fileStore.readGlobalFile(target_fasta_ids[hap], target_in)
            input_contigs = fasta_contig_names(target_in)
            input_caps = [l for l in cactus_call(parameters=['panpatch', '--telomere-report', target_in] + telo_opts,
                                                 check_output=True).split('\n') if l.startswith('#Contig')]
            os.remove(target_in)
            with open(combined) as f:
                record_names = [l[1:].split()[0] for l in f if l.startswith('>')]
            missing = missing_target_contigs(input_contigs, record_names, bed_out, report_out,
                                             run['samples'][0], hap)
            if missing:
                msg = '{} target contig(s) of {} are missing from the patched assembly: {}{}'.format(
                    len(missing), label, ', '.join(missing[:20]), ' ...' if len(missing) > 20 else '')
                if run['defaultSample']:
                    # --defaultSample replaces a reverted target contig with the default sample's on purpose
                    RealtimeLogger.warning('{} (replaced by --defaultSample {})'.format(msg, run['defaultSample']))
                else:
                    raise RuntimeError('cactus-panpatch: ' + msg)

        # the delivered names: each record that is one input contig takes that contig's name
        record_names = rename_fasta_records(combined, run['referenceFree'])

        # Measure telomere completeness on the assembly we are actually delivering: panpatch's own
        # "#Contig" lines are measured before this function renames and rescues, so re-measure the
        # finished FASTA with panpatch's own detector (--telomere-report, which shares seq_has_telomere
        # with the in-graph check) and use those lines.
        cap_out = cactus_call(parameters=['panpatch', '--telomere-report', combined] + telo_opts,
                              check_output=True)
        hap_caps = [l for l in cap_out.split('\n') if l.startswith('#Contig')]
        cap_lines += hap_caps
        summary_lines.append('#Summary hap={} input_contigs={} input_t2t={} output_records={} t2t={}'.format(
            user_hap, len(input_contigs), t2t_count(input_caps), len(record_names), t2t_count(hap_caps)))

        # in reference-free mode the single (hap0) output is named for the user's haplotype (already
        # in the run name); otherwise keep panpatch's haplotype numbering
        output_name = (run['name'] + '.fa.gz') if run['referenceFree'] else (run['name'] + '.hap{}.fa.gz'.format(hap))
        cactus_call(parameters=['bgzip', combined, '--threads', str(int(job.cores))])
        output_id_map[output_name] = job.fileStore.writeGlobalFile(combined + '.gz')

    # the delivered report and bed, in the user's names (the reference-free target's own sample and
    # haplotype, no #0 phase blocks): rows, then a "#Summary" line per haplotype, then the cap lines
    renames = run_renames(run)
    if not chrom_report_paths:
        # panpatch wrote no report at all (every target contig was left out of every graph): give the
        # rescued contigs' rows the header panpatch's report would have had
        p = os.path.join(work_dir, 'header.tsv')
        with open(p, 'w') as f:
            f.write(PANPATCH_REPORT_COLUMNS + '\tgap_origin\n')
        chrom_report_paths.append(p)
    concat_reports(chrom_report_paths, rescue_rows + summary_lines + cap_lines, report_out)
    for p in chrom_report_paths:
        os.remove(p)
    clean_report_names(report_out, renames,
                       hap=((run['hap'] if run['hap'] is not None else '0') if run['referenceFree'] else None))
    output_id_map[run['name'] + '.tsv'] = job.fileStore.writeGlobalFile(report_out)
    clean_bed_names(bed_out, renames)
    output_id_map[run['name'] + '.bed'] = job.fileStore.writeGlobalFile(bed_out)

    s = summarize_report(report_out)
    RealtimeLogger.info('{}: {} output records from {} input contigs, T2T {} -> {}; {} patches accepted, {} '
                        'rejected'.format(run['name'], s['output_records'], s['input_contigs'], s['input_t2t'],
                                          s['t2t'], sum(v for k, v in s['counts'].items() if k != 'bed-gap'),
                                          s['rejected']))
    return output_id_map

def rescue_dropped_contigs(job, fasta_path, target_fasta_id, seen_contigs, target_sample, hap, label=None):
    """ append to fasta_path (a panpatch output haplotype) any contig from the input target assembly
    that panpatch never saw -- ie that the pangenome dropped (binned to _AMBIGUOUS_ or unmapped).
    those contigs were never patched or scaffolded, so they are safe to carry through verbatim """
    work_dir = job.fileStore.getLocalTempDir()
    in_fa = os.path.join(work_dir, 'target_input.fa')
    job.fileStore.readGlobalFile(target_fasta_id, in_fa)
    with open(in_fa, 'rb') as fh:
        is_gzipped = fh.read(2) == b'\x1f\x8b'
    opener = gzip.open if is_gzipped else open

    rescued = []    # (contig, length)
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
                    rescued.append([contig, 0])
            elif keep:
                out.write(line)
                rescued[-1][1] += len(line.strip())
    if rescued:
        RealtimeLogger.info('{}: {} contig(s) the pangenome left out of every graph written as they are: {}{}'.format(
            label or os.path.basename(fasta_path), len(rescued), ', '.join(r[0] for r in rescued[:10]),
            ' ...' if len(rescued) > 10 else ''))
    return [tuple(r) for r in rescued]

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

    if not options.keepPangenome and not getattr(options, 'fromGraphs', None):
        # deleting the scratch tree costs one operation per file (~n_chroms * n_genomes of them),
        # not one per byte: seconds of rmtree locally, but a serial paginated delete loop on an
        # s3:// outDir, which is why it is off the coordination tier
        job.addFollowOnJobFn(cleanup_pangenome_wrapper, options, run, walltime=cactus_walltime(600))

def concat_reports(chrom_report_paths, cap_lines, out_path):
    """ concatenate the per-chromosome panpatch reports (in the order given) into one sample report:
    the TSV header once, then every chromosome's patch rows, then cap_lines: the caller's own lines (rows
    for the contigs it rescued, '#Summary' lines and the '#Contig' telomere cap lines).

    panpatch's own cap lines are dropped and replaced by cap_lines, which the caller measured on the
    finished assembly.  panpatch emits one line per graph path that contributed intervals -- donor paths
    included, with len= the sum of that path's intervals -- and measures them before the wrapper reverts
    masked target regions, so they describe neither the records nor the sequence the user receives. """
    with open(out_path, 'w') as out:
        wrote_header = False
        for p in chrom_report_paths:
            with open(p) as f:
                lines = f.readlines()
            if wrote_header:
                lines = lines[1:]                      # header only once
            out.writelines([l for l in lines if not l.startswith('#Contig')])
            wrote_header = True
        for l in cap_lines:
            out.write(l.rstrip('\n') + '\n')

def summarize_report(report_path):
    """ tally one run's panpatch report: accepted patches by category (and how many gap-fills were
    error-BED 'bed-gap' fills), rejected candidates, and from the '#Summary' lines gather_panpatch writes
    (one per haplotype) the input contigs, output records and telomere-to-telomere contigs before and after.
    passthrough (written-as-is) rows are neither.  a report with no '#Summary' line (an older one) counts
    output records and T2T from its '#Contig' cap lines """
    counts = {'gap-fill': 0, 'bed-gap': 0, 'scaffold': 0, 'telomere': 0}
    rejected = 0
    caps, cap_t2t = 0, 0
    totals = {'input_contigs': 0, 'input_t2t': 0, 'output_records': 0, 't2t': 0}
    have_summary = False
    idx = None
    with open(report_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#Contig'):                 # per-output-record telomere cap status
                caps += 1
                if 'left=YES' in line and 'right=YES' in line:
                    cap_t2t += 1
                continue
            if line.startswith('#Summary'):
                have_summary = True
                for kv in line.split()[1:]:
                    k, _, v = kv.partition('=')
                    if k in totals:
                        totals[k] += int(v)
                continue
            if not line or line.startswith('#'):
                continue
            if idx is None:                                # the header (first non-comment line)
                idx = {c: i for i, c in enumerate(line.split('\t'))}
                continue
            f2 = line.split('\t')
            d, t = idx.get('decision'), idx.get('type')
            if d is None or d >= len(f2):
                continue
            typ = f2[t] if t is not None and t < len(f2) else ''
            if f2[d] == 'rejected' and typ != 'passthrough':
                rejected += 1
            if f2[d] != 'accepted':
                continue
            if typ in ('gap-fill', 'scaffold', 'telomere'):
                counts[typ] += 1
                g = idx.get('gap_origin')
                if typ in ('gap-fill', 'telomere') and g is not None and g < len(f2) and f2[g] == 'bed-gap':
                    counts['bed-gap'] += 1
    if not have_summary:
        totals['output_records'], totals['t2t'] = caps, cap_t2t
    return dict(totals, counts=counts, rejected=rejected, contigs=caps)

SUMMARY_COLUMNS = ['run', 'input_contigs', 'output_records', 'input_t2t', 't2t', 'patches', 'gap_fill', 'bed_gap',
                   'scaffold', 'telomere', 'rejected']

def write_batch_summary(job, options, summaries):
    """ roll the per-run panpatch reports up into one summary (outDir/panpatch-summary.tsv): a row per run
    (a sample, or one haplotype of it when reference-free), then a TOTAL row.  summaries is a list of (run
    name, that run's output-file-id map) """
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

    def values(s):
        c = s['counts']
        return [s['input_contigs'], s['output_records'], s['input_t2t'], s['t2t'],
                c['gap-fill'] + c['scaffold'] + c['telomere'], c['gap-fill'], c['bed-gap'], c['scaffold'],
                c['telomere'], s['rejected']]
    total = [0] * (len(SUMMARY_COLUMNS) - 1)
    summary_path = os.path.join(work_dir, 'panpatch-summary.tsv')
    with open(summary_path, 'w') as out:
        out.write('# cactus-panpatch summary of {} run(s).  input_t2t / t2t: contigs with a telomere at both ends '
                  'before / after patching.  patches: accepted gap_fill + scaffold + telomere (bed_gap: the gap '
                  'fills and telomere patches over error-BED intervals).  rejected: candidate patches turned down (see each run\'s '
                  '.tsv for why)\n'.format(len(rows)))
        out.write('\t'.join(SUMMARY_COLUMNS) + '\n')
        for name, s in rows:
            v = values(s)
            out.write('\t'.join([name] + [str(x) for x in v]) + '\n')
            total = [a + b for a, b in zip(total, v)]
        out.write('\t'.join(['TOTAL'] + [str(x) for x in total]) + '\n')

    t = dict(zip(SUMMARY_COLUMNS[1:], total))
    RealtimeLogger.info('cactus-panpatch summary: {} run(s), {} output records from {} input contigs, T2T {} -> {}, '
                        '{} patches accepted (gap-fill {}, scaffold {}, telomere {}; {} bed-gap), {} rejected'.format(
                            len(rows), t['output_records'], t['input_contigs'], t['input_t2t'], t['t2t'],
                            t['patches'], t['gap_fill'], t['scaffold'], t['telomere'], t['bed_gap'], t['rejected']))
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
