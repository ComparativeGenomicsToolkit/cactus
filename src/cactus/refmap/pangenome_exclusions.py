#!/usr/bin/env python3

"""
Accounting for input sequence that does not make it into the output pangenome graphs.

Sequence goes missing between the input FASTAs and the output graphs in several ways: whole contigs
binned to _AMBIGUOUS_ (or dropped without any log record) during chromosome splitting, fragments
removed by clip-vg and the vg clip stages that follow it, and fragments removed by the allele
frequency (depth) filter.  None of these were reported completely, and clip-vg's own --out-bed is
wrong whenever paths carry subranges, which in cactus they always do after the "full" phase.

The method here does not trust any of those tools' self-reports.  For every phase graph the workflow
actually produced, path coverage is read back with `vg paths -E`, projected into input-contig
coordinates, and subtracted from a baseline built from the sanitized input FASTAs.  Because the
phases are strictly nested (filter subset of clip subset of full subset of input), each lost base is
attributed to the phase boundary at which it disappeared, by set subtraction alone -- no heuristics.

This module holds the pure, Toil-free half: name parsing, interval algebra and file parsing.  The
job functions that drive it live alongside in cactus_graphmap_join.py.
"""

import gzip
import logging
import os
import re
import shutil

from cactus.shared.common import cactus_call, getOptionalAttrib, findRequiredNode
from toil.realtimeLogger import RealtimeLogger

logger = logging.getLogger(__name__)

# every baseline table starts with exactly this line.  it is checked on read because
# <outDir>/chrom-subproblems/contig_sizes.tsv is a completely different (chromosome x sample matrix)
# file living in the same output tree, and feeding it to this parser would otherwise silently
# produce a nonsense baseline and a 100%-loss report
BASELINE_HEADER = '#event\tpansn_prefix\tbase_contig\toffset\tlength'

# written above the header so the table explains itself
BASELINE_PREAMBLE = [
    '# input contig lengths, as cactus saw them after sanitizing the fasta headers.',
    '# pansn_prefix    SAMPLE#HAP, the prefix this genome carries in the graph',
    '# base_contig     contig name in the graph, and in the exclusion BED files',
    '# offset, length  where this piece sits in base_contig. a contig that arrived whole is one',
    '#                 row at offset 0; one that arrived in fragments is several rows',
]


def open_maybe_gzip(path, mode='rt'):
    """ the baseline table is written gzipped, but accept it either way: users hand this file back
    with --inputContigSizes and may well have decompressed it """
    with open(path, 'rb') as probe:
        gzipped = probe.read(2) == b'\x1f\x8b'
    return gzip.open(path, mode) if gzipped else open(path, mode)

# reason tokens used in column 4 of the emitted BED files.  bare tokens, no whitespace, so that
# column 4 remains a legal BED name field
REASON_AMBIGUOUS = 'ambiguous'
REASON_UNASSIGNED = 'unassigned'
REASON_NO_CHROM_GRAPH = 'no_chromosome_graph'
REASON_UNALIGNED = 'unaligned'
REASON_CLIP = 'clip'
REASON_FILTER = 'filter'
REASON_REFGAP = 'refgap'

# stand-in chromosome for contigs that reached no chromosome graph at all, so that per-chromosome
# tables still sum to the per-genome totals
NO_CHROMOSOME = '_NONE_'

# a phase's name doubles as the reason token for sequence it removed.  spelled out rather than left
# implicit so that renaming a phase cannot silently rename a column-4 value users filter on
PHASE_REASON = {'clip': REASON_CLIP, 'filter': REASON_FILTER}

# sequence lost before any graph was built.  absent from every graph, whichever one you look at
DROPPED_REASONS = (REASON_AMBIGUOUS, REASON_UNASSIGNED, REASON_NO_CHROM_GRAPH, REASON_UNALIGNED)

# each graph is built from the previous one, so exclusions accumulate: what is missing from the
# clipped graph is everything missing from the full graph plus what clipping took, and so on
PHASE_REASONS = {'full': frozenset(DROPPED_REASONS),
                 'clip': frozenset(DROPPED_REASONS) | {REASON_CLIP},
                 'filter': frozenset(DROPPED_REASONS) | {REASON_CLIP, REASON_FILTER}}


def phase_tag(phase, filter_threshold=None):
    """ the infix cactus uses for a phase's output files: clip is the unsuffixed default, the full
    graph gets .full, and the frequency-filtered graph gets .dN """
    if phase == 'clip':
        return ''
    if phase == 'full':
        return '.full'
    return '.d{}'.format(filter_threshold)

# vg depth prints one line per reference base as `path <TAB> 1-based-pos <TAB> depth`, which is far
# too much to hand to python for a whole genome.  this reduces it to the zero-depth positions in the
# pipe, as refgaps.sh does, and converts 1-based positions to 0-based half-open BED in the same
# step: base p becomes [p-1, p), so a run p..q merges to [p-1, q) -- width preserved, and it can
# never extend past the end of the contig
DEPTH_TO_BED_AWK = '$3==0 {print $1 "\t" ($2-1) "\t" $2}'

_SUBRANGE_RE = re.compile(r'^(.*)\[(\d+)-(\d+)\]$')
_SPLIT_ASSIGNED_RE = re.compile(r'^Assigned (?:ref-)?contig to (\S+): id=([^|]+)\|(\S+)\s')
_SPLIT_AMBIGUOUS_RE = re.compile(r'^Query contig is ambiguous: id=([^|]+)\|(\S+)\s')


def resolve_subpath_naming(name):
    """ python port of hal2vg's resolve_subpath_naming (hal2vg.cpp).  cactus encodes a contig
    fragment as <CONTIG>_sub_<START>_<END>, and hal2vg turns that back into a path on <CONTIG> with
    subrange [START-END].  the encoding nests, in which case the offsets sum and the length is the
    innermost (first-parsed) one.

    returns (base_name, offset, length) with length None when there was no _sub_ suffix at all """
    first_length = None
    start_offset = 0
    while True:
        sp = name.rfind('_sub_')
        if sp == -1:
            break
        up = name.rfind('_')
        if up <= sp + 1:
            break
        start_tok, end_tok = name[sp + 5:up], name[up + 1:]
        if not start_tok.isdigit() or not end_tok.isdigit():
            # not a coordinate suffix -- a contig genuinely called e.g. "scaf_sub_unit".  hal2vg
            # would abort on such a name, so it cannot reach a graph; treat it as an opaque name
            break
        start, end = int(start_tok), int(end_tok)
        start_offset += start
        if first_length is None:
            first_length = end - start
        name = name[:sp]
    return name, (start_offset if first_length is not None else 0), first_length


def event_to_pansn_prefix(event):
    """ the SAMPLE#HAP prefix that graph paths for a given seqfile event carry.  mirrors hal2vg's
    resolve_haplotype_naming: a trailing .N with N all digits is the haplotype, otherwise 0 """
    dot = event.rfind('.')
    if dot != -1 and event[dot + 1:].isdigit():
        return '{}#{}'.format(event[:dot], int(event[dot + 1:]))
    return '{}#0'.format(event)


def pansn_prefix_to_event(pansn_prefix):
    """ inverse of event_to_pansn_prefix, for the degraded mode where the seqfile event names are
    not available.  SAMPLE#0 -> SAMPLE, SAMPLE#1 -> SAMPLE.1, so the two haplotypes of a diploid
    still land in separate files instead of being merged """
    toks = pansn_prefix.split('#')
    if len(toks) != 2:
        return pansn_prefix.replace('#', '.')
    return toks[0] if toks[1] == '0' else '{}.{}'.format(toks[0], toks[1])


def parse_path_name(name, path_length=None):
    """ split a vg path name into (pansn_prefix, base_contig, start, end).

    the grammar, verified against real cactus output:
        SAMPLE#HAP#CONTIG                 reference-sense    S288C#0#chrI
        SAMPLE#HAP#CONTIG#PHASEBLOCK      haplotype-sense    Y12#0#chrI#0
        ...either, with an optional [START-END] subrange   SK1#0#chrI#0[14059-228861]

    contig names contain no '#' (the sanitizer strips through the last one) and hal2vg has already
    resolved any _sub_ encoding into the subrange, so the third field is always a base contig.  the
    phase block, when present, is discarded.

    returns None for anything that is not PanSN, which the caller counts as an orphan """
    m = _SUBRANGE_RE.match(name)
    if m:
        base, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    else:
        base, start, end = name, 0, None
    toks = base.split('#')
    if len(toks) < 3:
        return None
    if end is None:
        end = start + path_length if path_length is not None else None
    return '#'.join(toks[:2]), toks[2], start, end


def merge_intervals(intervals):
    """ sort and merge, returning non-overlapping half-open [start, end) tuples """
    out = []
    for start, end in sorted(intervals):
        if start >= end:
            continue
        if out and start <= out[-1][1]:
            if end > out[-1][1]:
                out[-1][1] = end
        else:
            out.append([start, end])
    return [(s, e) for s, e in out]


def subtract_intervals(a, b):
    """ a \\ b, both merged on the way in """
    a, b = merge_intervals(a), merge_intervals(b)
    out = []
    j = 0
    for start, end in a:
        cur = start
        while j < len(b) and b[j][1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] < end:
            if b[k][0] > cur:
                out.append((cur, min(b[k][0], end)))
            cur = max(cur, b[k][1])
            if cur >= end:
                break
            k += 1
        if cur < end:
            out.append((cur, end))
    return out


def intersect_intervals(a, b):
    a, b = merge_intervals(a), merge_intervals(b)
    out = []
    i = j = 0
    while i < len(a) and j < len(b):
        start, end = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if start < end:
            out.append((start, end))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def total_bp(intervals):
    return sum(end - start for start, end in intervals)


def merge_sorted_bed_stream(line_iter, min_length, path_offsets=None):
    """ merge an already-sorted stream of `name <TAB> start <TAB> end` BED lines into maximal runs,
    keeping only those at least min_length long.  streaming, so a reference with no alignment at all
    costs O(1) memory instead of one tuple per base.

    yields (name, start, end), with any subrange offset of the named path added to both ends """
    path_offsets = path_offsets or {}
    cur_name, cur_start, cur_end = None, None, None

    def emit():
        if cur_name is None or cur_end - cur_start < min_length:
            return None
        offset = path_offsets.get(cur_name, 0)
        return (cur_name, cur_start + offset, cur_end + offset)

    for line in line_iter:
        toks = line.split('\t')
        if len(toks) < 3:
            continue
        name, start, end = toks[0], int(toks[1]), int(toks[2])
        if name == cur_name and start <= cur_end:
            cur_end = max(cur_end, end)
            continue
        got = emit()
        if got:
            yield got
        cur_name, cur_start, cur_end = name, start, end
    got = emit()
    if got:
        yield got


def contig_sizes_from_fai(fai_path, event):
    """ read a samtools faidx index of a *sanitized* cactus fasta into baseline rows.

    returns a list of (event, pansn_prefix, base_contig, offset, length) """
    pansn_prefix = event_to_pansn_prefix(event)
    rows = []
    with open(fai_path, 'r') as fai_file:
        for line in fai_file:
            toks = line.split('\t')
            if len(toks) < 2:
                continue
            name, length = toks[0], int(toks[1])
            # the sanitizer passes through a pre-existing id=X| prefix even when X is not the event,
            # so cut at the first '|' generically rather than matching id=<event>|
            if name.startswith('id='):
                bar = name.find('|')
                if bar > 0:
                    name = name[bar + 1:]
            base_contig, offset, sub_length = resolve_subpath_naming(name)
            if sub_length is not None and sub_length != length:
                logger.warning('subpath name %s in %s implies length %d but the fasta has %d; '
                               'using the fasta length', name, fai_path, sub_length, length)
            rows.append((event, pansn_prefix, base_contig, offset, length))
    return rows


def write_baseline_tsv(rows, out_path):
    opener = gzip.open if out_path.endswith('.gz') else open
    with opener(out_path, 'wt') as out_file:
        for line in BASELINE_PREAMBLE:
            out_file.write(line + '\n')
        out_file.write(BASELINE_HEADER + '\n')
        for row in rows:
            out_file.write('{}\t{}\t{}\t{}\t{}\n'.format(*row))


def check_baseline_header(path):
    """ raise unless path looks like a baseline contig-size table.  cheap enough to call during
    option validation, so a wrong file fails immediately rather than after the whole join has run """
    try:
        with open_maybe_gzip(path) as in_file:
            header = in_file.readline().rstrip('\n')
            while header.startswith('#') and header != BASELINE_HEADER:
                header = in_file.readline().rstrip('\n')
    except OSError as e:
        raise RuntimeError('could not read --inputContigSizes {}: {}'.format(path, e))
    if header != BASELINE_HEADER:
        raise RuntimeError(
            '{} is not a cactus baseline contig-size table: expected the header line\n  {}\n'
            'but found\n  {}\nThis file is written by cactus-pangenome as '
            '<outDir>/<outName>.input-contig-sizes.tsv.gz. Note it is NOT the same file as '
            '<outDir>/chrom-subproblems/contig_sizes.tsv, which is a chromosome-by-sample '
            'matrix and cannot be used here.'.format(path, BASELINE_HEADER, header))


def read_baseline_tsv(path):
    """ inverse of write_baseline_tsv.  hard-fails on a wrong header rather than silently parsing
    some other tsv (chrom-subproblems/contig_sizes.tsv is the one that would otherwise get passed) """
    rows = []
    check_baseline_header(path)
    with open_maybe_gzip(path) as in_file:
        for line in in_file:
            if not line.strip() or line.startswith('#'):
                continue
            toks = line.rstrip('\n').split('\t')
            if len(toks) != 5:
                raise RuntimeError('malformed line in {}: {}'.format(path, line.rstrip('\n')))
            rows.append((toks[0], toks[1], toks[2], int(toks[3]), int(toks[4])))
    return rows


def parse_split_log(path):
    """ read minigraph.split.log, which rgfa-split writes one line per contig it made a decision
    about.  returns (assigned, ambiguous) where assigned maps (event, contig, offset) -> chromosome
    and ambiguous is the set of (event, contig, offset) binned to _AMBIGUOUS_.

    the log names contigs the way cactus named them internally, so a contig that arrived as a
    fragment appears as e.g. chr1_sub_1000_2000.  those are resolved here into the same
    (contig, offset) space the graph and the input contig sizes use, so the three can be joined
    without anything having to carry the internal name around.

    contigs rgfa-split made no decision about appear in neither -- that silence is exactly the case
    this feature exists to surface, so it must not be conflated with either bucket """
    assigned, ambiguous = {}, set()
    with open(path, 'r') as log_file:
        for line in log_file:
            m = _SPLIT_ASSIGNED_RE.match(line)
            if m:
                contig, offset, _length = resolve_subpath_naming(m.group(3))
                assigned[(m.group(2), contig, offset)] = m.group(1)
                continue
            m = _SPLIT_AMBIGUOUS_RE.match(line)
            if m:
                contig, offset, _length = resolve_subpath_naming(m.group(2))
                ambiguous.add((m.group(1), contig, offset))
    return assigned, ambiguous


def baseline_by_key(rows):
    """ group baseline rows into {(pansn_prefix, base_contig): merged intervals}, which is the key
    graph coverage joins on """
    by_key = {}
    for _event, pansn_prefix, base_contig, offset, length in rows:
        by_key.setdefault((pansn_prefix, base_contig), []).append((offset, offset + length))
    return {key: merge_intervals(ivs) for key, ivs in by_key.items()}


def event_by_pansn_prefix(rows):
    """ {pansn_prefix: event}, so a graph path can be routed to the right output file """
    return {pansn_prefix: event for event, pansn_prefix, _c, _o, _l in rows}


def safe_event_filename(event):
    """ event names become file basenames; check_sample_names only constrains a leading '.' and a
    numeric suffix, so guard against a name that would escape the output directory """
    if os.sep in event or '/' in event or event in ('.', '..'):
        raise RuntimeError('genome name {} cannot be used as a file name'.format(event))
    return event


def synthetic_samples(config):
    """ the sample names that exist only inside the graph and have no input fasta: the minigraph
    "assembly" and, if ancestors are ever included, the cactus ancestral genomes.  read from the
    config rather than hardcoded -- assemblyName is a config attribute, and includeAncestor is
    shipped as 0 today but flipping it would otherwise silently reintroduce the same bug """
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"),
                                    "assemblyName", default="_MINIGRAPH_")
    return {graph_event, config.getDefaultInternalNodePrefix() + '0'}


###############################################################################
# toil job functions
###############################################################################

def path_coverage_job(job, config, vg_path, vg_id, chrom_name, phase):
    """ read back what a phase graph actually contains, as intervals in input-contig coordinates.

    returns (file_id, key_list, fragment_count) where file_id is a sorted, merged TSV of
    pansn_prefix / base_contig / start / end.  merging here keeps the bytes crossing the jobstore
    proportional to the number of surviving fragments rather than to the number of paths """
    work_dir = job.fileStore.getLocalTempDir()
    local_vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, local_vg_path)

    synthetic = synthetic_samples(config)
    coverage = {}
    orphan_names = []
    lines = cactus_call(parameters=['vg', 'paths', '-E', '-v', local_vg_path],
                        check_output=True, job_memory=job.memory)
    for line in lines.split('\n'):
        toks = line.split('\t')
        if len(toks) < 2:
            continue
        path_name, path_length = toks[0], int(toks[1])
        parsed = parse_path_name(path_name, path_length)
        if parsed is None:
            orphan_names.append(path_name)
            continue
        pansn_prefix, base_contig, start, end = parsed
        # drop synthetic events at parse time, not on the output: the degraded (no baseline table)
        # mode derives its baseline *from* this coverage, so a _MINIGRAPH_ path left in here would
        # become its own "input genome"
        if pansn_prefix.split('#')[0] in synthetic:
            continue
        if end - start != path_length:
            raise RuntimeError(
                'path {} in {} has subrange width {} but length {}; the assumption that a vg '
                'subrange spans exactly the path it names no longer holds'.format(
                    path_name, vg_path, end - start, path_length))
        coverage.setdefault((pansn_prefix, base_contig), []).append((start, end))

    if orphan_names:
        RealtimeLogger.warning('{} {}: {} path(s) are not PanSN and were skipped, e.g. {}'.format(
            chrom_name, phase, len(orphan_names), orphan_names[:3]))

    out_path = os.path.join(work_dir, '{}.{}.coverage.tsv'.format(chrom_name, phase))
    fragment_count = 0
    with open(out_path, 'w') as out_file:
        for (pansn_prefix, base_contig) in sorted(coverage.keys()):
            for start, end in merge_intervals(coverage[(pansn_prefix, base_contig)]):
                out_file.write('{}\t{}\t{}\t{}\n'.format(pansn_prefix, base_contig, start, end))
                fragment_count += 1

    return (job.fileStore.writeGlobalFile(out_path), sorted(coverage.keys()),
            fragment_count, len(orphan_names))


def read_coverage_tsv(path):
    """ inverse of path_coverage_job's output """
    coverage = {}
    with open(path, 'r') as in_file:
        for line in in_file:
            toks = line.rstrip('\n').split('\t')
            if len(toks) != 4:
                continue
            coverage.setdefault((toks[0], toks[1]), []).append((int(toks[2]), int(toks[3])))
    return coverage


def ref_gaps_job(job, config, options, vg_path, vg_id, chrom_name, ref_event, drop_graph_event):
    """ the reference is never clipped, so its exclusions cannot be measured by absence.  measure
    them the way refgaps.sh does instead: runs along the reference with no other assembly aligned.

    `vg depth -m0 -P <ref>#` prints one line per reference base as `path <TAB> 1-based-pos <TAB>
    depth`, where depth excludes the reference path itself.  It is streamed, never materialised """
    work_dir = job.fileStore.getLocalTempDir()
    local_vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, local_vg_path)

    min_length = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"),
                                   "refGapMinLength", typeFn=int, default=10000)
    if min_length <= 0:
        return None

    depth_input_path = local_vg_path
    if drop_graph_event:
        # falling back to the full graph, which still carries minigraph fragments.  they would
        # count towards depth and mask real gaps, so take them out first
        graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"),
                                        "assemblyName", default="_MINIGRAPH_")
        depth_input_path = local_vg_path + '.nomg'
        cactus_call(parameters=['vg', 'paths', '-d', '-S', graph_event, '-x', local_vg_path],
                    outfile=depth_input_path, job_memory=job.memory)

    # note: vg depth already reports positions in original contig coordinates -- it strips the
    # subrange from the path name and adds subrange.first to the offset itself (vg's
    # algorithms/coverage_depth.cpp path_depths).  So nothing here may add the offset a second time.
    # vg depth emits one line per reference base, so it is filtered down to the zero-depth
    # positions in the pipe (exactly as refgaps.sh does) before anything touches disk.  awk also
    # does the 1-based -> 0-based half-open conversion, emitting one [pos-1, pos) singleton per
    # zero base, which merge_sorted_bed_stream then joins back into runs
    zeros_path = os.path.join(work_dir, '{}.zeros.bed'.format(chrom_name))
    cactus_call(parameters=[['vg', 'depth', depth_input_path, '-m0', '-P', ref_event + '#'],
                            ['awk', DEPTH_TO_BED_AWK]],
                outfile=zeros_path, job_memory=job.memory)

    out_path = os.path.join(work_dir, '{}.refgaps.bed'.format(chrom_name))
    rows = 0
    with open(zeros_path, 'r') as zeros_file, open(out_path, 'w') as out_file:
        for path_name, start, end in merge_sorted_bed_stream(zeros_file, min_length):
            parsed = parse_path_name(path_name)
            contig = parsed[1] if parsed else path_name
            out_file.write('{}\t{}\t{}\t{}\n'.format(contig, start, end, REASON_REFGAP))
            rows += 1
    os.remove(zeros_path)

    if not rows:
        return None
    return job.fileStore.writeGlobalFile(out_path)


def _strip_subrange(name):
    m = _SUBRANGE_RE.match(name)
    return m.group(1) if m else name


PHASE_ORDER = ('full', 'clip', 'filter')


def resolve_assigned_graph(where, chrom_names, other_contig):
    """ the graph a split-log assignment actually landed in.  rgfa-split logs the fine-grained
    reference contig (e.g. chr14_GL000009v2_random), but cactus then folds every reference contig
    not given its own graph into the otherContig graph (chrOther by default).  So an assignment
    naming a contig with no graph of its own resolves to otherContig when that was built, and to
    None only when the sequence really reached no graph at all.  Returns the graph name or None """
    if where in chrom_names:
        return where
    if other_contig and other_contig in chrom_names:
        return other_contig
    return None


def assigned_chromosome(key_event, base_contig, offsets, split_log, chrom_names, other_contig):
    """ the graph the splitting stage put this contig in, if any """
    if split_log is None:
        return None
    assigned, _ambiguous = split_log
    for offset in offsets:
        key = (key_event, base_contig, offset)
        if key in assigned:
            return resolve_assigned_graph(assigned[key], chrom_names, other_contig)
    return None


def sublabel_dropped(key_event, base_contig, offsets, full_cov_bp, split_log, chrom_names,
                     other_contig):
    """ which of the four split-stage fates a wholly-or-partly absent contig met.

    the three split-level labels are gated on *zero* full-phase coverage: a contig with any
    coverage demonstrably reached a chromosome graph, so whatever it lost there is alignment-level,
    not split-level, no matter what the log says about it """
    if full_cov_bp > 0 or split_log is None:
        return REASON_UNALIGNED
    assigned, ambiguous = split_log
    for offset in offsets:
        if (key_event, base_contig, offset) in ambiguous:
            return REASON_AMBIGUOUS
    for offset in offsets:
        key = (key_event, base_contig, offset)
        if key in assigned:
            return (REASON_UNALIGNED
                    if resolve_assigned_graph(assigned[key], chrom_names, other_contig)
                    else REASON_NO_CHROM_GRAPH)
    return REASON_UNASSIGNED


def compute_exclusions(baseline_rows, coverage_paths_by_phase, phases_present, split_log,
                       chrom_names, out_dir, filter_threshold=None, other_contig=None):
    # chrom_names is order-significant: entry i names the chromosome whose coverage file is at
    # index i in each phase's list
    """ the pure core of the exclusion report: turn per-phase coverage into per-genome BED files
    plus a summary, and return the counters.

    baseline_rows is None in degraded mode (no --inputContigSizes), in which case the baseline is
    taken from the widest present phase and the 'dropped' class cannot be reported at all.

    processing is chromosome by chromosome.  each input contig lives in exactly one chromosome
    subproblem, so peak memory scales with the largest chromosome rather than the whole panel """
    order = [p for p in PHASE_ORDER if p in phases_present]
    if not order:
        raise RuntimeError('exclusion report needs at least one phase graph')
    widest = order[0]
    explicit = baseline_rows is not None


    if explicit:
        baseline = baseline_by_key(baseline_rows)
        event_of = event_by_pansn_prefix(baseline_rows)
        offsets_of = {}
        for event, pansn_prefix, base_contig, offset, _l in baseline_rows:
            offsets_of.setdefault((pansn_prefix, base_contig), []).append(offset)
        input_bp = {}
        for event, _p, _c, _o, length in baseline_rows:
            input_bp[event] = input_bp.get(event, 0) + length
    else:
        baseline, event_of, offsets_of, input_bp = {}, {}, {}, {}
        for path in coverage_paths_by_phase[widest]:
            for key, ivs in read_coverage_tsv(path).items():
                baseline.setdefault(key, []).extend(ivs)
        baseline = {k: merge_intervals(v) for k, v in baseline.items()}
        for pansn_prefix, _contig in baseline:
            event_of.setdefault(pansn_prefix, pansn_prefix_to_event(pansn_prefix))
        for key, ivs in baseline.items():
            event = event_of[key[0]]
            input_bp[event] = input_bp.get(event, 0) + total_bp(ivs)

    # sequence missing from every graph can only be seen against a real input baseline: without one
    # the baseline IS the widest graph, so that class would measure zero by construction
    report_dropped = explicit and 'full' in phases_present

    # one BED per (graph, genome).  writing them in one pass would need phases x genomes open
    # handles, which overruns the usual 1024 limit on a big panel, so a single reason-tagged file
    # per genome is written here and split per graph afterwards
    work_dir = os.path.join(out_dir, '.by-reason')
    os.makedirs(work_dir, exist_ok=True)
    writers = {}

    def writer_for(event):
        if event not in writers:
            writers[event] = open(
                os.path.join(work_dir, '{}.bed'.format(safe_event_filename(event))), 'w')
        return writers[event]

    counts = {}            # (event, reason) -> [intervals, bp]
    whole_dropped = {}     # event -> [(contig, bp, reason)] absent from every graph
    # one row per (chromosome, genome, contig): how much of that contig survived into each graph
    contig_rows = []
    outside_baseline_bp = 0
    orphan_keys = []
    orphan_bp = 0
    seen_keys = set()

    def tally(event, reason, intervals):
        if not intervals:
            return
        slot = counts.setdefault((event, reason), [0, 0])
        slot[0] += len(intervals)
        slot[1] += total_bp(intervals)

    def emit(event, contig, intervals, reason):
        if not intervals:
            return
        out_file = writer_for(event)
        for start, end in intervals:
            out_file.write('{}\t{}\t{}\t{}\n'.format(contig, start, end, reason))
        tally(event, reason, intervals)

    n_chroms = len(coverage_paths_by_phase[widest])
    for i in range(n_chroms):
        cov = {}
        for phase in order:
            cov[phase] = read_coverage_tsv(coverage_paths_by_phase[phase][i])
        keys = set()
        for phase in order:
            keys.update(cov[phase].keys())

        for key in sorted(keys):
            pansn_prefix, base_contig = key
            if key not in baseline:
                bp = total_bp(merge_intervals(cov[widest].get(key, [])))
                orphan_keys.append(key)
                orphan_bp += bp
                continue
            seen_keys.add(key)
            base_ivs = baseline[key]
            event = event_of[pansn_prefix]

            clamped = {}
            for phase in order:
                raw = merge_intervals(cov[phase].get(key, []))
                clamped[phase] = intersect_intervals(raw, base_ivs)
                outside_baseline_bp += total_bp(raw) - total_bp(clamped[phase])

            missing = {p: subtract_intervals(base_ivs, clamped[p]) for p in order}

            if report_dropped:
                reason = sublabel_dropped(event, base_contig, offsets_of.get(key, [0]),
                                          total_bp(clamped['full']), split_log, chrom_names,
                                          other_contig)
                emit(event, base_contig, missing['full'], reason)
                if reason != REASON_NO_CHROM_GRAPH and not total_bp(clamped['full']):
                    whole_dropped.setdefault(event, []).append(
                        (base_contig, total_bp(base_ivs), reason))
            for prev, cur in zip(order, order[1:]):
                gained = subtract_intervals(missing[cur], missing[prev])
                emit(event, base_contig, gained, PHASE_REASON[cur])

            contig_rows.append((chrom_names[i], event, base_contig, total_bp(base_ivs),
                                dict((p, total_bp(clamped[p])) for p in order),
                                dict((p, len(clamped[p])) for p in order)))

    # contigs that reached no chromosome graph at all never appear in any coverage file
    if report_dropped:
        for key, base_ivs in sorted(baseline.items()):
            if key in seen_keys:
                continue
            pansn_prefix, base_contig = key
            event = event_of[pansn_prefix]
            reason = sublabel_dropped(event, base_contig, offsets_of.get(key, [0]), 0,
                                      split_log, chrom_names, other_contig)
            emit(event, base_contig, base_ivs, reason)
            if reason != REASON_NO_CHROM_GRAPH:
                whole_dropped.setdefault(event, []).append(
                    (base_contig, total_bp(base_ivs), reason))
            # in no chromosome graph, so it gets the synthetic chromosome _NONE_ rather than being
            # left out -- otherwise the per-chromosome rows would not sum to the genome total
            # a contig assigned to a chromosome nobody built still belongs to that chromosome:
            # attributing it keeps the per-chromosome view complete under --refContigs
            where = assigned_chromosome(event, base_contig, offsets_of.get(key, [0]),
                                        split_log, chrom_names, other_contig) or NO_CHROMOSOME
            contig_rows.append((where, event, base_contig, total_bp(base_ivs),
                                dict((p, 0) for p in order), dict((p, 0) for p in order)))

    for out_file in writers.values():
        out_file.close()

    # split into one directory per graph.  every input genome gets a file in every directory,
    # including the genomes that lost nothing
    phase_dirs = {}
    for phase in order:
        keep = PHASE_REASONS[phase]
        phase_dir = os.path.join(out_dir, 'clipped' + phase_tag(phase, filter_threshold) + '.beds')
        os.makedirs(phase_dir, exist_ok=True)
        phase_dirs[phase] = phase_dir
        for event in sorted(input_bp):
            src = os.path.join(work_dir, '{}.bed'.format(safe_event_filename(event)))
            dest = os.path.join(phase_dir, '{}.bed'.format(safe_event_filename(event)))
            with open(dest, 'w') as out_file:
                if os.path.exists(src):
                    with open(src, 'r') as in_file:
                        for line in in_file:
                            if line.rstrip('\n').split('\t')[3] in keep:
                                out_file.write(line)
    shutil.rmtree(work_dir)

    return {'counts': counts, 'input_bp': input_bp, 'phase_dirs': phase_dirs,
            'whole_dropped': whole_dropped, 'contig_rows': contig_rows,
            'phases': order, 'explicit_baseline': explicit, 'report_dropped': report_dropped,
            'outside_baseline_bp': outside_baseline_bp, 'orphan_paths': len(orphan_keys),
            'orphan_bp': orphan_bp, 'orphan_keys': orphan_keys[:10]}


# a sample losing this much of itself to sequence that reached no graph at all is not normal
# divergence, it is a sign the assembly did not fit the reference's chromosome structure
WARN_DROPPED_FRACTION = 0.05
# ...and a sample losing this much more than the rest of the panel stands out even if the cause
# is spread across classes.  both must hold, so a uniformly lossy run does not flag every sample
WARN_OUTLIER_RATIO = 3.0
WARN_OUTLIER_FRACTION = 0.10
# a reference with this much of itself unaligned to anything means the panel barely agrees with it
WARN_REFGAP_FRACTION = 0.20


def detect_problems(result, refgap_bp=None, ref_input_bp=None):
    """ look over the finished accounting for the shapes that mean something went wrong, as opposed
    to sequence being legitimately clipped or filtered.

    the two things that separate a real problem from a configuration choice are *which* class the
    loss falls into and whether it is *uniform*.  sequence dropped before any graph was built is
    always worth knowing about; sequence lost to a chromosome nobody asked to build
    (no_chromosome_graph, what --refContigs does on purpose) is not, and neither is a panel where
    everything loses about the same amount.

    returns a list of human-readable warnings, empty when nothing looks wrong """
    warnings = []
    counts, input_bp = result['counts'], result['input_bp']

    # 1. sequence that reached no graph at all.  no_chromosome_graph is deliberately excluded: it
    #    is what --refContigs does, and it applies uniformly to every sample
    never_placed = (REASON_AMBIGUOUS, REASON_UNASSIGNED, REASON_UNALIGNED)
    for event in sorted(input_bp):
        total_in = input_bp[event]
        if not total_in:
            continue
        bp = sum(counts.get((event, r), [0, 0])[1] for r in never_placed)
        if bp > total_in * WARN_DROPPED_FRACTION:
            by_reason = ', '.join('{}={} bp'.format(r, counts[(event, r)][1])
                                  for r in never_placed if (event, r) in counts)
            warnings.append(
                '{}: {:.1f}% of its input ({} bp) is in no graph at all ({}). Sequence in these '
                'classes never reached a chromosome graph, so it was lost during chromosome '
                'splitting rather than by clipping. A large fraction usually means the assembly '
                'disagrees with the reference chromosome structure -- inter-chromosomal '
                'rearrangements, a mis-joined scaffold, or the wrong reference.'.format(
                    event, 100.0 * bp / total_in, bp, by_reason))

    # 2. whole input contigs absent from every graph
    for event in sorted(result.get('whole_dropped', {})):
        contigs = result['whole_dropped'][event]
        bp = sum(c[1] for c in contigs)
        shown = ', '.join('{} ({} bp, {})'.format(*c) for c in contigs[:5])
        warnings.append(
            '{}: {} whole contig(s) totalling {} bp are absent from every graph: {}{}'.format(
                event, len(contigs), bp, shown, ', ...' if len(contigs) > 5 else ''))

    # 3. one sample much worse than the rest.  needs both a ratio and a floor so that a panel which
    #    is uniformly lossy (or uniformly fine) does not flag anything
    fractions = {}
    for event, total_in in input_bp.items():
        if total_in:
            lost = sum(bp for (ev, _r), (_n, bp) in counts.items() if ev == event)
            fractions[event] = lost / total_in
    if len(fractions) > 2:
        ordered = sorted(fractions.values())
        median = ordered[len(ordered) // 2]
        for event in sorted(fractions):
            frac = fractions[event]
            if frac > WARN_OUTLIER_FRACTION and median > 0 and frac > median * WARN_OUTLIER_RATIO:
                warnings.append(
                    '{}: loses {:.1f}% of its input, against {:.1f}% for the typical genome in this '
                    'panel. One genome behaving very differently from the rest usually means a '
                    'problem with that assembly rather than with the alignment.'.format(
                        event, 100.0 * frac, 100.0 * median))

    # 4. the reference barely agrees with the panel
    for ref_event, (_n, bp) in (refgap_bp or {}).items():
        if ref_input_bp and bp > ref_input_bp * WARN_REFGAP_FRACTION:
            warnings.append(
                '{}: {:.1f}% of the reference has no other assembly aligned to it. Either the '
                'panel is very distant from this reference, or the alignment largely '
                'failed.'.format(ref_event, 100.0 * bp / ref_input_bp))

    # 5. the accounting did not add up, which means a bug rather than a data problem
    if result['outside_baseline_bp']:
        warnings.append(
            'INTERNAL: {} bp of graph path coverage falls outside the input contig bounds. The '
            'exclusion report may be wrong; please report this.'.format(result['outside_baseline_bp']))
    if result['orphan_paths']:
        warnings.append(
            'INTERNAL: {} graph path group(s) have no matching input contig (e.g. {}). The '
            'exclusion report may be wrong; please report this.'.format(
                result['orphan_paths'], result['orphan_keys']))

    return warnings


def write_summary(result, out_path, refgap_bp=None, notes=None):
    """ the concise view: one block per input genome, always, including zero-loss ones """
    refgap_bp = refgap_bp or {}
    counts, input_bp = result['counts'], result['input_bp']
    with open(out_path, 'w') as out_file:
        out_file.write('# input sequence that is not in the output graph, by genome.\n')
        out_file.write('# NOTE "clipped" in these filenames means any input base absent from the '
                       'graph, whatever removed it.\n')
        out_file.write('#      the reason column below says which stage did, and one of its values '
                       "is the narrower 'clip'\n")
        out_file.write('# phases: {}\n'.format(','.join(result['phases'])))
        out_file.write('# baseline: {}\n'.format(
            '--inputContigSizes' if result['explicit_baseline'] else
            'derived from {}-phase coverage; the dropped class is unavailable'.format(
                result['phases'][0])))
        out_file.write('# TOTAL = input bp absent from the deepest graph that was built. refgap '
                       'is not counted in it\n')
        out_file.write('# every reason measured for this run gets a row for every genome, so a 0 '
                       'means measured-and-nothing-lost\n')
        out_file.write('# lost at chromosome assignment, so absent from EVERY graph:\n')
        out_file.write('#   ambiguous           = contig binned to _AMBIGUOUS_ during splitting\n')
        out_file.write('#   unassigned          = contig that splitting made no decision about\n')
        out_file.write('#   no_chromosome_graph = contig assigned to a chromosome no graph was '
                       'built for (e.g. left out of --refContigs)\n')
        out_file.write('#   unaligned           = contig reached its chromosome graph, but some of '
                       'its sequence did not make the full graph\n')
        out_file.write('# lost later, so present in the earlier graphs:\n')
        out_file.write('#   clip   = removed by the clip phase (clip-vg -u/-a, vg clip -d1/-sS'
                       ', vg clip -D with --delEdgeFilter)\n')
        out_file.write('#   filter = removed by the filter phase (vg clip -d <filter> -m, vg clip -sS)\n')
        out_file.write('# refgap = reference bp with zero non-reference depth. PRESENT in the graph.'
                       ' NOT lost. Excluded from TOTAL lost.\n')
        out_file.write('# outside_baseline_bp: {}\n'.format(result['outside_baseline_bp']))
        out_file.write('# orphan_paths: {}\n'.format(result['orphan_paths']))
        for note in (notes or []):
            out_file.write('# NOTE: {}\n'.format(note))
        out_file.write('#genome\treason\tintervals\tbp\tpct_of_input\n')

        def pct(bp, total):
            return 'NA' if not total else '{:.3f}'.format(100.0 * bp / total)

        # the reasons this run could measure at all: the dropped classes only when there is a
        # baseline and a full graph, and each later phase only when it ran
        measured = []
        if result['report_dropped']:
            measured += list(DROPPED_REASONS)
        for phase in result['phases'][1:]:
            measured.append(PHASE_REASON[phase])

        grand_lost = grand_refgap = 0
        for event in sorted(input_bp):
            total_in = input_bp[event]
            lost_n = lost_bp = 0
            for reason in measured:
                n, bp = counts.get((event, reason), (0, 0))
                out_file.write('{}\t{}\t{}\t{}\t{}\n'.format(event, reason, n, bp,
                                                                pct(bp, total_in)))
                lost_n += n
                lost_bp += bp
            if not result['report_dropped']:
                out_file.write('{}\tdropped\tNA\tNA\tNA\n'.format(event))
            out_file.write('{}\tTOTAL\t{}\t{}\t{}\n'.format(event, lost_n, lost_bp,
                                                             pct(lost_bp, total_in)))
            if event in refgap_bp:
                n, bp = refgap_bp[event]
                out_file.write('{}\trefgap\t{}\t{}\t{}\n'.format(event, n, bp,
                                                                   pct(bp, total_in)))
                grand_refgap += bp
            grand_lost += lost_bp
        total_in = sum(input_bp.values())
        out_file.write('TOTAL\tTOTAL\t\t{}\t{}\n'.format(grand_lost, pct(grand_lost, total_in)))
        if refgap_bp:
            out_file.write('TOTAL\trefgap\t\t{}\t{}\n'.format(grand_refgap,
                                                                pct(grand_refgap, total_in)))


def write_per_contig_stats(result, contig_path, chrom_path):
    """ the same accounting as the exclusion report, rolled up two ways: one row per input contig,
    and one row per genome per reference contig.  both answer "how much of this survived", at the
    two zoom levels people actually ask about """
    phases = result['phases']
    pairs = []
    for phase in phases:
        pairs += ['{}_bp'.format(phase), '{}_frags'.format(phase)]
    preamble = [
        '# how much of the input survived into each graph that was built.',
        '# NOTE "clipped" in these filenames means any input base absent from the graph, whatever',
        '#      removed it -- chromosome splitting, the clip phase, or the frequency filter.',
        '# ref_chrom   the reference contig whose graph it landed in, or {} if it reached none'.format(
            NO_CHROMOSOME),
        '# genome      the input genome, as named in the seqfile',
        '# input_bp    its length in the input fasta',
        '# <graph>_bp    bases of it present in that graph',
        '# <graph>_frags contiguous pieces it is broken into in that graph',
        '# bases lost at a stage are the drop between two adjacent _bp columns; everything lost',
        '# overall is input_bp - {}_bp'.format(phases[-1]),
    ]

    def row(prefix, in_bp, in_frags):
        cells = []
        for phase in phases:
            cells += [str(in_bp[phase]), str(in_frags[phase])]
        return prefix + '\t' + '\t'.join(cells) + '\n'

    with gzip.open(contig_path, 'wt') as out_file:
        for line in preamble:
            out_file.write(line + '\n')
        out_file.write('#genome\tcontig\tref_chrom\tinput_bp\t{}\n'.format('\t'.join(pairs)))
        for chrom, event, contig, input_bp, in_bp, in_frags in sorted(
                result['contig_rows'], key=lambda r: (r[1], r[2], r[0])):
            out_file.write(row('{}\t{}\t{}\t{}'.format(event, contig, chrom, input_bp),
                               in_bp, in_frags))

    rolled = {}
    for chrom, event, _contig, input_bp, in_bp, in_frags in result['contig_rows']:
        slot = rolled.setdefault((chrom, event),
                                 [0, 0, dict((p, 0) for p in phases), dict((p, 0) for p in phases)])
        slot[0] += 1
        slot[1] += input_bp
        for phase in phases:
            slot[2][phase] += in_bp[phase]
            slot[3][phase] += in_frags[phase]

    chrom_opener = gzip.open if chrom_path.endswith('.gz') else open
    with chrom_opener(chrom_path, 'wt') as out_file:
        for line in preamble[:1] + preamble[2:]:
            out_file.write(line + '\n')
        out_file.write('#ref_chrom\tgenome\tcontigs\tinput_bp\t{}\n'.format('\t'.join(pairs)))
        for (chrom, event) in sorted(rolled):
            contigs, input_bp, in_bp, in_frags = rolled[(chrom, event)]
            out_file.write(row('{}\t{}\t{}\t{}'.format(chrom, event, contigs, input_bp),
                               in_bp, in_frags))


def compute_exclusions_job(job, config, options, phase_coverage, baseline_id, split_log_id,
                           refgap_ids, chrom_names):
    """ join every phase's coverage against the baseline and write the report.

    phase_coverage is {phase: [(file_id, key_list, fragment_count, orphan_count), ...]} covering
    only the phases that actually ran -- a phase that did not run contributes no rows at all, rather
    than being charged its predecessor's survivors """
    work_dir = job.fileStore.getLocalTempDir()
    out_dir = os.path.join(work_dir, 'exclusions')
    os.makedirs(out_dir, exist_ok=True)

    coverage_paths = {}
    for phase, results in phase_coverage.items():
        coverage_paths[phase] = []
        for i, (file_id, _keys, _frags, _orphans) in enumerate(results):
            local = os.path.join(work_dir, '{}.{}.cov'.format(i, phase))
            job.fileStore.readGlobalFile(file_id, local)
            coverage_paths[phase].append(local)

    notes = []
    # without the input contig lengths there is nothing to measure the graphs against.  the graphs
    # could be measured against each other, but the result would carry the same filenames and column
    # names while meaning something else -- missing.full would be empty by construction rather than
    # by measurement, and pct_of_input would be a fraction of the graph rather than of the input.
    # Reference gaps do not depend on the baseline, so those are still produced.
    baseline_rows = None
    if baseline_id:
        baseline_path = os.path.join(work_dir, 'input-contig-sizes.tsv')
        job.fileStore.readGlobalFile(baseline_id, baseline_path)
        baseline_rows = read_baseline_tsv(baseline_path)
    else:
        RealtimeLogger.warning(
            'No exclusion report: it needs the input contig lengths, which are only available when '
            'the whole pipeline runs. Pass --inputContigSizes <outName>.input-contig-sizes.tsv.gz from '
            'the cactus-pangenome run that produced these graphs to get one.')

    split_log = None
    if split_log_id:
        split_log_path = os.path.join(work_dir, 'minigraph.split.log')
        job.fileStore.readGlobalFile(split_log_id, split_log_path)
        split_log = parse_split_log(split_log_path)
    elif baseline_id:
        notes.append('no split log was available, so sequence missing from every graph is reported '
                     'as {} without distinguishing the split-stage causes.'.format(REASON_UNALIGNED))

    other_contig = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_split"),
                                     "otherContigName", typeFn=str, default="chrOther")
    result = None
    if baseline_rows is not None:
        result = compute_exclusions(baseline_rows, coverage_paths, list(phase_coverage.keys()),
                                    split_log, chrom_names, out_dir,
                                    filter_threshold=getattr(options, 'filter', None),
                                    other_contig=other_contig)

    # reference gaps: concatenated per reference genome, kept well away from the loss totals
    refgap_bp = {}
    refgap_out = {}
    for ref_event, ids in (refgap_ids or {}).items():
        ids = [i for i in ids if i]
        if not ids:
            continue
        name = 'refgaps.bed.gz' if len(refgap_ids) == 1 else 'refgaps.{}.bed.gz'.format(
            safe_event_filename(ref_event))
        path = os.path.join(out_dir, name)
        rows = bp = 0
        with gzip.open(path, 'wt') as out_file:
            for file_id in ids:
                local = job.fileStore.readGlobalFile(file_id)
                with open(local, 'r') as in_file:
                    for line in in_file:
                        toks = line.split('\t')
                        if len(toks) < 3:
                            continue
                        out_file.write(line)
                        rows += 1
                        bp += int(toks[2]) - int(toks[1])
        refgap_bp[ref_event] = (rows, bp)
        refgap_out[name] = job.fileStore.writeGlobalFile(path)

    if result is None:
        return {'summary.tsv': None, 'missing': {}, 'refgaps': refgap_out, 'WARNING': None}

    summary_path = os.path.join(out_dir, 'summary.tsv')
    write_summary(result, summary_path, refgap_bp=refgap_bp, notes=notes)

    contig_stats_path = os.path.join(out_dir, 'clipped-by-input-contig.tsv.gz')
    chrom_stats_path = os.path.join(out_dir, 'clipped-by-reference-contig.tsv')
    write_per_contig_stats(result, contig_stats_path, chrom_stats_path)

    # a run can succeed and still have gone badly wrong.  surface the shapes that mean that, in a
    # file that only exists when there is something to say
    ref_input_bp = None
    if refgap_bp:
        ref_input_bp = result['input_bp'].get(sorted(refgap_bp)[0])
    problems = detect_problems(result, refgap_bp=refgap_bp, ref_input_bp=ref_input_bp)
    warning_id = None
    if problems:
        warning_path = os.path.join(out_dir, 'WARNING')
        with open(warning_path, 'w') as warning_file:
            warning_file.write('cactus-pangenome finished, but the exclusion accounting looks '
                               'wrong. Each note below is a heuristic, not a failure -- check it '
                               'against your data before acting.\n\n')
            for problem in problems:
                warning_file.write('* ' + problem + '\n\n')
            warning_file.write('Full detail: the clipped-by-*.tsv tables and clipped*.tar.gz '
                               'archives in <outName>.stats/.\n')
        warning_id = job.fileStore.writeGlobalFile(warning_path)
        # logged from inside the job: RealtimeLogger is the only channel that reaches the console
        # both here and under cactus-pangenome, where export runs in a worker
        for problem in problems:
            RealtimeLogger.warning('exclusion report: ' + problem)
        RealtimeLogger.warning(
            'The above means {}.WARNING was written. This run finished, but its exclusion '
            'accounting looks wrong -- please read that file.'.format(
                getattr(options, 'outName', '<outName>')))

    if result['orphan_paths']:
        RealtimeLogger.warning(
            'exclusion report: {} graph path group(s) have no input contig, e.g. {}'.format(
                result['orphan_paths'], result['orphan_keys']))
    if result['outside_baseline_bp']:
        RealtimeLogger.warning('exclusion report: {} bp of graph coverage fall outside the input '
                               'contig bounds'.format(result['outside_baseline_bp']))
    with open(summary_path, 'r') as summary_file:
        RealtimeLogger.info('Exclusion report:\n' + summary_file.read())

    # one tarball per graph, so the common case -- "give me what is missing from the graph I am
    # using" -- is a single download with no post-processing
    missing_out = {}
    for phase, phase_dir in result['phase_dirs'].items():
        base = os.path.basename(phase_dir)
        tar_path = os.path.join(work_dir, base + '.tar.gz')
        cactus_call(parameters=['tar', 'czf', tar_path, '-C', out_dir, base])
        missing_out[base + '.tar.gz'] = job.fileStore.writeGlobalFile(tar_path)

    return {'summary.tsv': job.fileStore.writeGlobalFile(summary_path),
            'clipped-by-input-contig.tsv.gz': job.fileStore.writeGlobalFile(contig_stats_path),
            'clipped-by-reference-contig.tsv': job.fileStore.writeGlobalFile(chrom_stats_path),
            'missing': missing_out,
            'refgaps': refgap_out,
            'WARNING': warning_id}


def contig_sizes_job(job, seq_id_map, graph_event):
    """ build the exclusion report's baseline from the *sanitized* input fastas, whose contig names
    are already in the namespace the graph paths and the split log use.

    this has to run early, right after sanitization, because the fastas do not survive to the join:
    sanitize_fasta_header deletes its input from the jobstore, and the sanitized copies are cleaned
    up once splitting is done.  only this small table is carried forward.

    one child per genome, then a follow-on to concatenate """
    per_event = {}
    for event in sorted(seq_id_map.keys()):
        if event == graph_event:
            continue
        safe_event_filename(event)
        fa_id = seq_id_map[event]
        per_event[event] = job.addChildJobFn(contig_sizes_for_event, fa_id, event,
                                             disk=fa_id.size * 3).rv()
    return job.addFollowOnJobFn(merge_contig_sizes, per_event).rv()


def contig_sizes_for_event(job, fa_id, event):
    """ faidx one sanitized genome and return its baseline rows """
    work_dir = job.fileStore.getLocalTempDir()
    fa_path = os.path.join(work_dir, '{}.fa'.format(event))
    job.fileStore.readGlobalFile(fa_id, fa_path)
    cactus_call(parameters=['samtools', 'faidx', fa_path])
    return contig_sizes_from_fai(fa_path + '.fai', event)


def merge_contig_sizes(job, per_event_rows):
    work_dir = job.fileStore.getLocalTempDir()
    rows = []
    for event in sorted(per_event_rows.keys()):
        rows += per_event_rows[event]
    out_path = os.path.join(work_dir, 'input-contig-sizes.tsv.gz')
    write_baseline_tsv(rows, out_path)
    RealtimeLogger.info('exclusion report baseline: {} contigs over {} genomes, {} bp'.format(
        len(rows), len(set(r[0] for r in rows)), sum(r[4] for r in rows)))
    return job.fileStore.writeGlobalFile(out_path)
