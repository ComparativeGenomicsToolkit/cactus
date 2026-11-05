#!/usr/bin/env python3

"""
Telomere anchors for a pangenome alignment: every contig's telomere array is aligned to one synthetic telomere per
chromosome arm, which joins the minigraph genome, so that cactus aligns the telomeres of all haplotypes to each other
whatever subtelomere they hang off.

Minigraph's mapping never covers telomere arrays, not even the reference's own, and BAR only reaches a telomere from
the inner end of its block.  So without anchors the arrays of different haplotypes stay unaligned, or are smeared over
a short reference array, and clipping removes them, along with any subtelomere block (often one swapped in from
another chromosome's end) that only the telomere tied to the graph.  Per alignment job (a chromosome, or several
reference contigs), just before cactus_consolidated:

0. The telomere repeat is the job's motif, a unit of 5-12bp written on its G strand, which runs 5' to 3' toward every
   chromosome end: TTAGGG in vertebrates and most animals, TTTAGGG in most plants, TTAGG in most insects.  Read from a
   chromosome end inward, the telomere is its reverse complement, the unit (CCCTAA, p = 6bp, in the examples below).
   Of several candidate motifs, the job takes the one with a perfect run (in any phase) of at least max(min_run, 12
   units) starting within max_outer of the most contig ends, read inward, if those are at least 2, at least 4 times
   as many as have one read the other way (a satellite reads the same both ways) and at least 3 times as many as any
   other candidate has; otherwise it anchors nothing.  A run only counts at the end whose half of the contig it starts
   in, and one read the other way not at all if it is in the telomere the other end reads inward (its run, and its
   array if that starts in its half, even where the run doesn't): a short contig's telomere, or one made of telomere,
   is read the other way from its other end.  An end counts for one candidate, the one whose run starts nearest it:
   in a perfect array of one candidate the longest perfect run of another is 12bp, but real arrays hold variant
   repeats, and HPRC telomeres have perfect runs of TTGGGG (a candidate, and a TTAGGG variant) of up to 102bp further
   in.  Each on its own, CHM13 and the chromosome-level assemblies of three Arabidopsis thaliana, maize Mo17 and
   tomato VF36 chose TTAGGG or TTTAGGG at 8 to 48 ends, with no other candidate at more than 1.
1. Each contig end is read from the end inward, on the strand that reads CCCTAA, for a telomere array.  Its core is
   the first 100bp window (starting at most max_outer + 100bp from the end) that is at least half telomeric units
   (CCCTAA, and for TTAGGG its common variants) and whose walk inward, in 100bp windows with one failing window
   allowed, reaches min_array.  Its outer edge is the outermost telomeric k-mer (p - 1bp, but 5 to 6bp) reachable
   from the core in steps of at most 40bp, so degenerate telomere stays in the array and other junk (such as
   (TTAAAG)n) stays out.  An edge under two repeat units from the contig end is moved to the end; an array whose edge
   is further in than max_outer is not anchored.  Its inner edge is the end of the last telomeric unit in its last
   window that starts within two units (12bp) of the one before it: the window can reach up to 50bp into the
   subtelomere, which would be forced onto other haplotypes' telomeres, and a unit the subtelomere has by chance is
   on its own.
   Nor is an array without min_run bp of perfect repeat ((CCCTAA)n): telomerase makes long perfect runs, while an
   interstitial telomeric repeat (ITS), which a contig that breaks next to it would leave at its end, is degenerate,
   and anchoring it would join the telomere to the middle of the chromosome.  Over CHM13 and the chromosome ends of 21
   HPRC haplotypes, the longest perfect run of an ITS of at least 500bp was 54bp, and the shortest of a telomere 84bp
   (91% are over 300bp); 5 of the 998 arrays at the contig ends of 22 other, mostly not T2T, HPRC haplotypes fell
   under 72bp, all of them telomeres.  These are human numbers: of the 44 chromosome ends of the maize Mo17 and tomato
   VF36 T2T assemblies, 11 have arrays too short or too degenerate to anchor, and tomato has ITS with perfect runs of
   2kb.  Another motif's array has no run at all (a plant's (CCCTAAA)n has 6bp of (CCCTAA)n).  A core whose array has
   no such run is passed over for a later one, if any has it: a degenerate stretch at the very end, whose walk can
   reach min_array too, would otherwise hide the telomere behind it.
2. Each contig is anchored to the targets of one reference contig, its key: the rc:Z: tag (the reference contig of
   the minigraph segments) of a strict majority of the bases its tagged PAF lines align, or, for a contig of the
   reference event, its own name, whatever its lines' tags (cactus-graphmap's rgfa2paf lines have none).  A contig
   without tagged lines, or without a majority, is not anchored, nor is the reference if the PAF has no tags at all
   (one not made by cactus-graphmap).  So a job can hold any number of reference contigs (cactus-pangenome's
   whole-genome --noSplit job, the chrOther bin of the small ones, a plain cactus-align), each with its own targets.
   The key is the whole contig's, not each end's, as in a chromosome's job, whose lines all carry the chromosome's
   tag: so in a job of several, a telomere on a piece of another chromosome (a translocation or a misjoin) goes to the
   targets of the chromosome that most of its contig lies on.  The contig's orientation is the PAF strand that holds
   most of the bases its lines on its key align (all its lines, if none is tagged with it), and puts each array on the
   p or q arm.  A contig without PAF lines can't be oriented, and is not anchored.
3. Target id=<minigraph event>|tel_<key>p is the unit's repeat, (CCCTAA)n, with the chromosome end at position 0, and
   tel_<key>q is its reverse complement, (TTAGGG)n with the end at its last base; a _sub_ in tel_<key> is written
   _sub-, as hal2vg would read tel_X_sub_0_100p and tel_X_sub_0_100q as the same subpath of X.  Each is the longest
   array anchored to it + 1000bp long, rounded up to a whole number of units, and is only made for an arm with an
   anchor.
4. Each array is aligned to its target with the outer ends pinned together, at the repeat offset (0 to p - 1) that
   matches best over the outermost 600bp, and the inner end free (edlib, mode SHW): a slip in the repeat costs an
   indel rather than putting the rest of the array out of phase.  These anchor lines are tagged ct:i:1, which makes
   them constraints for cactus_consolidated, kept whatever chain they end up in.
   The reference's array is pinned one repeat unit further in (p to 2p - 1bp into the target, against 0 to p - 1bp
   for the others), at the offset whose alignment of its outermost 600bp has the fewest edits and starts with an
   aligned pair.
   vg haplotypes only samples a chromosome end from the haplotypes that pass through the reference's outermost node,
   and the arrays' outer ends differ by a few bp (the phase they end in, and slips near the end): pinned with the
   others, the reference's end was reached by as few as 3 of 40 haplotypes, whose ends were then the only ones
   sampled.  One unit further in, the other arrays overhang it by up to 14bp, which the join clips as stubs, so that
   their paths start on the reference's end, unless a slip puts a deletion or a mismatch on its outermost base.
5. The existing PAF lines of an anchored contig are trimmed out of its anchored ends (the array and anything outside
   it), cutting their cigars, so that nothing else pairs the telomere up some other way.  A line left with no aligned
   pair is dropped.  Every other line is written unchanged.

On 40 HPRC haplotypes x 8 chromosomes (with minigraph end events), anchoring kept 99.9% of the telomere sequence in
the clipped graph (all but the few bp of each haplotype outside the reference's end), against 28.5% without it, and 98%
of the swapped subtelomere blocks.  618 of the 640 haplotype ends reached the reference's end, and vg haplotypes could
sample 37-40 of the 40 haplotypes at every end (3-39 with the reference pinned like the others).

Copyright (C) 2026 by Glenn Hickey

Released under the MIT license, see LICENSE.txt
"""

import collections
import functools
import gzip
import logging
import multiprocessing
import os
import re
import shutil
import time
import zlib
from concurrent.futures import Future, ProcessPoolExecutor

import edlib

logger = logging.getLogger(__name__)

WINDOW = 100     # arrays are measured in windows of this many bp
SCAN = 100000    # bp read beyond max_outer at each contig end, to measure an array
OUTER_GAP = 40   # the outer edge steps outward over telomeric k-mers at most this far apart
PHASE_BP = 600   # the anchor's repeat offset is the one that matches best over this many outer bp
TARGET_PAD = 1000  # a target is this much longer than the longest array anchored to it
WORKER_BYTES = 1 << 30  # a process per this many bytes of fasta to scan
MAX_WORKERS = 16   # the scan is bound by reading the fastas well before this many processes

# the telomere motifs (G strand units) a job chooses from: vertebrates and most animals, plants, most insects,
# nematodes, beetles (TCAGG, TTGGG), green algae, Cestrum, ciliates (TTGGGG, TTTTGGGG) and Allium
DEFAULT_MOTIFS = ['TTAGGG', 'TTTAGGG', 'TTAGG', 'TTAGGC', 'TCAGG', 'TTGGG', 'TTTTAGGG', 'TTTTTTAGGG', 'TTGGGG',
                  'TTTTGGGG', 'CTCGGTTATGGG']
# variant repeat units (G strand) that are telomeric too in a motif's arrays: TTAGGG's common ones, in human arrays
VARIANTS = {'TTAGGG': ['TCAGGG', 'TGAGGG', 'TTGGGG', 'TTCGGG']}
RUN_UNITS = 12   # a run that chooses a motif is at least this many units long, as well as min_run bp
MIN_ENDS = 2     # the chosen motif has a run at this many contig ends at least,
STRAND_RATIO = 4     # this many times as many as have one read the other way,
RUNNER_UP_RATIO = 3  # and this many times as many as any other candidate

# a contig end's arm, from the end (5' or 3') and the contig's orientation relative to the reference
ARM = {('5p', '+'): 'p', ('5p', '-'): 'q', ('3p', '+'): 'q', ('3p', '-'): 'p'}
ANCHOR_MAPQ = 60
CONSTRAINT_TAG = 'ct:i:1'

COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
CIGAR_OPS = re.compile(r'(\d+)([=XIDM])')
WHITESPACE = b'\r\n\t '

# why a contig end was not anchored, in the order the log lists them
REASONS = ['no_motif', 'no_array', 'short_array', 'array_beyond_max_outer', 'degenerate_array', 'duplicate_name',
           'no_paf_lines', 'bad_paf_line', 'no_reference_contig', 'ambiguous_reference_contig', 'length_mismatch',
           'arrays_overlap', 'already_anchored']


def reverse_complement(s):
    return s.translate(COMPLEMENT)[::-1]


# a telomere motif, as make_motif() derives it from its G strand unit
Motif = collections.namedtuple('Motif', ['g', 'unit', 'p', 'repeat_re', 'repeat_start_re', 'longest', 'edge_kmer_re',
                                         'edge_reach', 'perfect_re', 'tel_unit', 'inner_gap', 'partial',
                                         'reference_inset'])


@functools.lru_cache(maxsize=None)
def make_motif(g):
    """ the Motif of telomere repeat unit g (on the G strand, as conventionally written: TTAGGG), everything the scan,
    the anchors and the targets take from it.  unit: g's reverse complement (CCCTAA), read from a chromosome end
    inward, p bp long.  repeat_re: the telomeric repeats, the unit and its VARIANTS (reverse complemented, but for any
    that is a piece of the unit's cycle), as an alternation of the longest first (it takes the first that matches),
    longest the bp of that one.  The edge k-mers: those of the unit's cycle (p - 1bp, as TTAGGG's 5, but at
    least 5bp, as random sequence has 4-mers of a cycle often enough to carry the edge over every other 40bp step, and
    at most 6bp, as degenerate telomere has fewer longer ones), and the repeats that don't start with one; edge_reach,
    how far past a position one that starts before it reaches.  inner_gap, partial (an outer edge this close to the
    end is moved to it) and reference_inset: two units, under two units and one unit """
    unit, p = reverse_complement(g), len(g)
    repeats = [unit] + [v for v in map(reverse_complement, VARIANTS.get(g, [])) if v not in unit * (len(v) // p + 2)]
    alternation = '|'.join(sorted(repeats, key=len, reverse=True))
    k = max(5, min(p - 1, 6))
    kmers = sorted(set((unit * (k // p + 2))[i:i + k] for i in range(p)))
    edge = kmers + [r for r in repeats if r[:k] not in kmers]
    return Motif(g=g, unit=unit, p=p, repeat_re=re.compile(alternation),
                 repeat_start_re=re.compile('(?=({}))'.format(alternation)), longest=max(map(len, repeats)),
                 edge_kmer_re=re.compile('(?=(?:{}))'.format('|'.join(edge))), edge_reach=max(map(len, edge)) - 1,
                 perfect_re=re.compile('(?:{})+'.format(unit)), tel_unit={'p': unit, 'q': g}, inner_gap=2 * p,
                 partial=2 * p - 1, reference_inset=p)


TTAGGG = make_motif('TTAGGG')   # the motif of the functions that take one, unless given


def check_motifs(motifs):
    """ the telomere motifs (G strand units) upper case, or a ValueError if there are none, or one is not 5-12bp of A,
    C, G and T, is a repeat of a shorter unit, or is a rotation of another or of one's reverse complement (they would
    find each other's telomeres, or one read the other way) """
    motifs = [g.strip().upper() for g in motifs]
    if not motifs:
        raise ValueError('no telomere motifs')
    for i, g in enumerate(motifs):
        if not re.fullmatch('[ACGT]{5,12}', g):
            raise ValueError('telomere motif {!r} is not 5-12bp of A, C, G and T'.format(g))
        if any(g == g[:d] * (len(g) // d) for d in range(1, len(g)) if len(g) % d == 0):
            raise ValueError('telomere motif {} is a repeat of a shorter unit'.format(g))
        for j, h in enumerate(motifs[:i + 1]):
            if len(h) != len(g):
                continue
            if j < i and h == g:
                raise ValueError('telomere motif {} is listed twice'.format(g))
            if j < i and h in g + g:
                raise ValueError('telomere motifs {} and {} are rotations of one unit'.format(h, g))
            if reverse_complement(h) in g + g:
                raise ValueError('telomere motif {} is a rotation of the reverse complement of {}'.format(
                    g, h if j < i else 'itself'))
    return motifs


def config_motifs(graphmap_node):
    """ the telomere motifs of <graphmap telomereMotifs> (separated by commas or spaces; DEFAULT_MOTIFS if it is
    missing), checked by check_motifs(), so that a pipeline can fail on a bad one before it aligns anything """
    motifs = graphmap_node.attrib.get('telomereMotifs') if graphmap_node is not None else None
    try:
        return check_motifs(motifs.replace(',', ' ').split() if motifs is not None else DEFAULT_MOTIFS)
    except ValueError as e:
        raise RuntimeError('<graphmap telomereMotifs="{}">: {}'.format(motifs, e)) from e


def telomere_fraction(w, motif):
    """ fraction of w covered by telomeric repeats, taken left to right without overlaps """
    return sum(map(len, motif.repeat_re.findall(w))) / max(1, len(w))


def array_walk(s, x0, motif):
    """ bp of telomere array in s from x0 inward: whole windows that are at least half telomeric, with one failing
    window allowed between two that pass """
    x = run = x0
    fails = 0
    while x + WINDOW <= len(s):
        if telomere_fraction(s[x:x + WINDOW], motif) >= 0.5:
            x += WINDOW
            run = x
            fails = 0
        else:
            fails += 1
            if fails > 1:
                break
            x += WINDOW
    return run - x0


def outer_edge(s, c, motif):
    """ the array's outer edge in s, from c (a telomeric repeat of its core) outward: the outermost telomeric k-mer
    reachable in steps of at most OUTER_GAP; 0 if that is within motif.partial (under two units) of the end """
    e = c
    for o in reversed([m.start() for m in motif.edge_kmer_re.finditer(s, 0, c + motif.edge_reach) if m.start() < c]):
        if e - o > OUTER_GAP:
            break
        e = o
    return 0 if e <= motif.partial else e


def inner_edge(s, x, motif):
    """ the array's inner edge in s, given x, where its last passing window ends: the end of the last telomeric
    repeat that starts in that window, within the motif's inner gap of the end of the one before it, so that a repeat
    the subtelomere has by chance doesn't count """
    hits = [(m.start(), m.end(1)) for m in motif.repeat_start_re.finditer(s, max(0, x - 2 * WINDOW),
                                                                           x + motif.longest - 1) if m.start() < x]
    for (a, a_end), (b, b_end) in zip(reversed(hits[:-1]), reversed(hits)):
        if b - a_end <= motif.inner_gap:
            return b_end
    return x


def array_at_end(s, min_array, max_outer, min_run=0, motif=TTAGGG):
    """ s: a contig end read from the end inward on the motif's C strand (as CCCTAA), upper case.  (outer_bp,
    array_bp, status, open_end), the array being s[outer_bp : outer_bp + array_bp].  status is ok; beyond (its edge is
    further in than max_outer: not to be anchored); short (no core reaches min_array: this is the first telomeric
    run); or none (no telomeric window).  The array is the first core's whose walk reaches min_array and that has
    min_run bp of perfect repeat, or if none has, the first whose walk reaches min_array (one without the run, not to
    be anchored).  open_end: its walk ran out of s, so the array may go on past it """
    short = degenerate = None
    walked = 0
    for x0 in range(0, max_outer + WINDOW + 1, WINDOW):
        # a core within an earlier core's walk would end where it does: its array would be that one or a part of it
        if x0 < walked or telomere_fraction(s[x0:x0 + WINDOW], motif) < 0.5:
            continue
        L = array_walk(s, x0, motif)
        walked = x0 + L
        e = outer_edge(s, motif.repeat_start_re.search(s, x0).start(), motif)
        # a walk only stops short of the end of s at two failing windows
        open_end = x0 + L + 2 * WINDOW > len(s)
        if L < min_array:
            if short is None:
                short = (e, inner_edge(s, x0 + L, motif) - e, 'short', open_end)
            continue
        x = inner_edge(s, x0 + L, motif)
        array = (e, x - e, ('ok' if e <= max_outer else 'beyond'), open_end)
        if perfect_run(s[e:x], motif) >= min_run:
            return array
        degenerate = degenerate or array
    return degenerate or short or (0, 0, 'none', False)


def perfect_run(array, motif=TTAGGG):
    """ bp of the longest run of the motif's unit (CCCTAA) perfectly repeated in array """
    return max((len(m) for m in motif.perfect_re.findall(array)), default=0)


def run_finder(unit, need):
    """ what telomere_run() looks for: a perfect run of unit (a C or G strand unit), in any phase, of at least need bp.
    (need, p, the whole units every such run holds, [(the first need bp of a run in a phase, the regex of whole units
    in that phase)] per phase) """
    p = len(unit)
    phases = [unit[i:] + unit[:i] for i in range(p)]
    return (need, p, unit * ((need - p + 1) // p),
            [((r * (need // p + 2))[:need], re.compile('(?:{})+'.format(r))) for r in phases])


def telomere_run(s, finder, limit):
    """ (start, end) in s of the outermost perfect run of finder's unit (run_finder()), in any phase (a stretch every
    position of which starts a rotation of the unit), that is at least finder's need bp long and starts at most limit
    bp from the start of s; None if there is none """
    need, p, whole, phases = finder
    if s.find(whole, 0, limit + need) < 0:   # most ends, quickly
        return None
    found = [(s.find(first, 0, limit + need), run_re) for first, run_re in phases]
    found = [f for f in found if f[0] >= 0]
    if not found:
        return None
    start, run_re = min(found, key=lambda f: f[0])
    end = run_re.match(s, start).end()
    while end < len(s) and s[end] == s[end - p]:   # a partial unit after the whole ones
        end += 1
    return start, end


def tally_contig(tally, length, runs5, runs3):
    """ add a contig of length bp to tally (choose_motif()'s), given per end (5' and 3'), per motif: its run read
    inward (telomere_run()'s, or None), the extent of the telomere it reads inward (that run, and its array if that
    starts in its half of the contig; None if neither), and its run read the other way.  An end with runs read inward
    of several motifs counts for the one whose run starts nearest the end (then the longest).  A run read the other
    way that lies in the telomere the contig's other end reads inward is not counted: a contig that is mostly
    telomere reads it the other way from that end, which says nothing of its strand (so one that is mostly satellite
    counts as a telomere at one end: the two can't be told apart) """
    for this, other in ((runs5, runs3), (runs3, runs5)):
        # an end counts for one motif, the one whose run starts nearest it (then the longest): a telomere's variant
        # repeats can make a perfect run of another candidate further in (TTAGGG's TTGGGG, up to 102bp in HPRC)
        best = min((i for i, (inward, _, _) in enumerate(this) if inward),
                   key=lambda i: (this[i][0][0], this[i][0][0] - this[i][0][1], i), default=None)
        for i, (t, (inward, _, wrong), (_, extent, _)) in enumerate(zip(tally, this, other)):
            if i == best:
                t[0] += 1
                t[2] += inward[1] - inward[0]
            # [a, b) from this end is [length - b, length - a) from the other
            if wrong and not (extent and length - wrong[1] < extent[1] and extent[0] < length - wrong[0]):
                t[1] += 1


def choose_motif(motifs, tally, n_ends):
    """ the motif of a job of several candidate motifs, given tally, per motif [contig ends with a perfect run of it
    read inward (telomere_run(); an end counts for one motif), ends with one read the other way, bp of the first ones'
    runs], and how many ends were read: the one at the most ends (then with the most bp), if they are at least
    MIN_ENDS, STRAND_RATIO times as many as read it the other way and RUNNER_UP_RATIO times as many as the runner-up
    has.  Returns (its index, or None if it isn't the job's motif; a line for the log; whether the line is a
    warning) """
    order = sorted(range(len(motifs)), key=lambda i: (-tally[i][0], -tally[i][2], i))
    (S, Wr, _), (S2, _, _) = tally[order[0]], tally[order[1]]
    counts = '{} at {} of {} contig ends ({} read the other way), then {} at {}'.format(motifs[order[0]], S, n_ends,
                                                                                         Wr, motifs[order[1]], S2)
    if S < MIN_ENDS or S < STRAND_RATIO * Wr:
        return None, 'no telomere motif ({}): nothing anchored'.format(counts), False
    if S < RUNNER_UP_RATIO * S2:
        return None, ('no telomere motif ({}): nothing anchored; to anchor one of them, give it alone in '
                      'telomereMotifs'.format(counts)), True
    return order[0], 'telomere motif {}'.format(counts), False


def open_maybe_gzip(path):
    """ a plain or gzipped file opened for reading bytes, by its first bytes, as cactus_cons names every fasta .fa """
    with open(path, 'rb') as f:
        gzipped = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rb') if gzipped else open(path, 'rb')


def strip_whitespace(b):
    """ b without line breaks (fast) or any other whitespace """
    b = b.replace(b'\n', b'')
    if b'\r' in b or b'\t' in b or b' ' in b:
        b = b.translate(None, WHITESPACE)
    return b


def add_sequence(record, piece, n):
    """ add the raw sequence piece to record, [name, length, first n bases, last n bases] """
    piece = strip_whitespace(piece)
    record[1] += len(piece)
    if len(record[2]) < n:
        record[2] += piece[:n - len(record[2])]
    record[3] = piece[-n:] if len(piece) >= n else (record[3] + piece)[-n:]


def fasta_ends(path, n, block_size=1 << 23):
    """ (name, length, head, tail) per record of a plain or gzipped fasta: head and tail are its first and last n
    bases (bytes, whitespace removed; each the whole sequence if it is shorter).  The file is read once, a block at a
    time, and only those ends of a record are kept, so a whole genome takes little memory """
    assert n > 0
    record = None
    carry = b''
    line_start = True   # whether the next block starts a line
    with open_maybe_gzip(path) as f:
        while True:
            block = f.read(block_size)
            if not block:
                if not carry:
                    break
                block = b'\n'   # a header line without a line break at the end of the file
            if carry:
                block, carry = carry + block, b''
            i = 0
            while True:
                # a header is a '>' at the start of a line; sequence lines never hold one
                j = block.find(b'>', i)
                while j > 0 and block[j - 1] != 10 or j == 0 and not line_start:
                    j = block.find(b'>', j + 1)
                if j < 0:
                    if record is not None and i < len(block):
                        add_sequence(record, block[i:], n)
                    break
                if record is not None and j > i:
                    add_sequence(record, block[i:j], n)
                e = block.find(b'\n', j)
                if e < 0:   # the header goes on in the next block
                    carry = block[j:]
                    break
                if record is not None:
                    yield tuple(record)
                fields = block[j + 1:e].split()
                record = [fields[0].decode('utf-8', 'surrogateescape') if fields else '', 0, b'', b'']
                i = e + 1
            line_start = bool(carry) or block.endswith(b'\n')
    if record is not None:
        yield tuple(record)


def scan_genome(event, path, min_array, max_outer, min_run, motifs):
    """ (event, contigs, tally, error) for one genome's fasta and each of motifs (G strand units).  contigs: (name,
    length, ends), ends being its 5' then 3' end as (end, arrays), arrays holding per motif None (no array) or
    (outer_bp, array_bp, status, array, reaches_scan_limit, perfect_run): array the sequence of an ok array with min_run
    bp of perfect repeat (any other is not anchored), read from the end inward as CCCTAA, reaches_scan_limit whether it
    may go on past the max_outer + SCAN bp read, and perfect_run that of the ok array.  tally, of several motifs (None
    for one), is the genome's part of choose_motif()'s: per motif [contig ends with a perfect run of at least
    max(min_run, RUN_UNITS units) starting within max_outer read inward, ends with one read the other way, bp of the
    first ones' runs], a run counting only at the end whose half of the contig it starts in (and see tally_contig()).
    A fasta that can't be read gives (event, None, None, the error) """
    motifs = [make_motif(g) for g in motifs]
    finders = []   # with one motif there is nothing to choose
    if len(motifs) > 1:
        finders = [[run_finder(u, max(min_run, RUN_UNITS * m.p)) for u in (m.unit, m.g)] for m in motifs]
    tally = [[0, 0, 0] for _ in finders] or None
    contigs = []
    try:
        for name, length, head, tail in fasta_ends(path, max_outer + SCAN):
            ends, runs = [], []
            # a run that starts in the other half of a short contig is that end's
            limit = min(max_outer, (length - 1) // 2)
            for end, s in (('5p', head.decode('latin-1').upper()),
                           ('3p', reverse_complement(tail.decode('latin-1').upper()))):
                arrays, end_runs = [], []
                for i, motif in enumerate(motifs):
                    outer, L, status, open_end = array_at_end(s, min_array, max_outer, min_run, motif)
                    array = s[outer:outer + L] if status == 'ok' else None
                    run = perfect_run(array, motif) if array else 0
                    arrays.append((outer, L, status, array if run >= min_run else None,
                                   open_end and len(s) == max_outer + SCAN, run) if L > 0 else None)
                    if finders:
                        inward, other_way = (telomere_run(s, finder, limit) for finder in finders[i])
                        # the telomere this end reads inward: its run, and its array if that starts in this end's half
                        # (to the end of s if its walk ran out of s), as a short contig's can when its run doesn't
                        extent = inward
                        if L > 0 and outer <= limit:
                            hi = len(s) if open_end else outer + L
                            extent = (min(inward[0], outer), max(inward[1], hi)) if inward else (outer, hi)
                        end_runs.append((inward, extent, other_way))
                ends.append((end, arrays))
                runs.append(end_runs)
            contigs.append((name, length, ends))
            if tally:
                tally_contig(tally, length, *runs)
    except (OSError, EOFError, zlib.error) as e:
        return event, None, None, '{}: {}'.format(path, e)
    return event, contigs, tally, None


def fasta_names(path):
    """ the record names of a plain or gzipped fasta """
    return [name for name, _, _, _ in fasta_ends(path, 1)]


def cigar_ops(cigar):
    return [(int(n), op) for n, op in CIGAR_OPS.findall(cigar)]


def cigar_string(ops):
    """ ops as a cigar, adjacent ops of one kind merged and empty ones left out """
    out = []
    for n, op in ops:
        if n == 0:
            continue
        if out and out[-1][1] == op:
            out[-1] = (out[-1][0] + n, op)
        else:
            out.append((n, op))
    return ''.join('{}{}'.format(n, op) for n, op in out)


def trim_paf(t, lo, hi):
    """ the PAF line t (a list of fields) trimmed to the query interval [lo, hi), cutting its cigar at op level, so
    that it starts and ends on aligned pairs.  Returns t itself if nothing changes, None if no aligned pair is left.
    Raises ValueError if it has no cigar or its cigar doesn't fit its coordinates """
    qs, qe, strand = int(t[2]), int(t[3]), t[4]
    lo, hi = max(qs, lo), min(qe, hi)
    if hi <= lo:
        return None
    if (lo, hi) == (qs, qe):
        return t
    cg_i = next((i for i in range(12, len(t)) if t[i].startswith('cg:Z:')), None)
    if cg_i is None:
        raise ValueError('no cigar')
    fwd = strand == '+'
    q, tp = (qs if fwd else qe), int(t[7])
    pieces = []  # (op, length, target start, query start, query end) in target order
    aligned = 0  # bp of aligned pairs (=, X and M)
    for n, op in cigar_ops(t[cg_i][5:]):
        if op in '=XM':
            aligned += n
        if op == 'D':
            pieces.append(('D', n, tp, None, None))
            tp += n
            continue
        # an op that takes query covers [q, q+n) on '+', and [q-n, q) walking down on '-'
        k0, k1 = (max(0, lo - q), min(n, hi - q)) if fwd else (max(0, q - hi), min(n, q - lo))
        if k1 > k0:
            q0, q1 = (q + k0, q + k1) if fwd else (q - k1, q - k0)
            pieces.append((op, k1 - k0, tp + k0 if op != 'I' else None, q0, q1))
        if op != 'I':
            tp += n
        q = q + n if fwd else q - n
    if q != (qe if fwd else qs) or tp != int(t[8]):
        raise ValueError('cigar does not fit the coordinates')
    while pieces and pieces[0][0] in 'ID':
        pieces.pop(0)
    while pieces and pieces[-1][0] in 'ID':
        pieces.pop()
    if not pieces:
        return None
    first, last = pieces[0], pieces[-1]
    t = list(t)
    t[2], t[3] = (str(first[3]), str(last[4])) if fwd else (str(last[3]), str(first[4]))
    t[7], t[8] = str(first[2]), str(last[2] + last[1])
    t[cg_i] = 'cg:Z:' + cigar_string([(p[1], p[0]) for p in pieces])
    matches = sum(p[1] for p in pieces if p[0] == '=')
    m_bp = sum(p[1] for p in pieces if p[0] == 'M')
    if m_bp:
        # an M is a match or a mismatch: count its bp at the line's identity, so that an rgfa2paf line (all M, all
        # matches) keeps every bp as a match
        matches += round(m_bp * int(t[9]) / aligned)
    t[9] = str(matches)
    t[10] = str(sum(p[1] for p in pieces))
    return t


def leading_gap(ops):
    """ bp of indels before the first aligned pair of ops """
    n = 0
    for k, op in ops:
        if op in '=X':
            break
        n += k
    return n


def pin_outer_base(array, target, start, ops):
    """ (start, ops) of an alignment of array to target from start, with any indels before its first aligned pair
    replaced by aligned pairs of the array's bases before that pair and the target bases just outside it, so that the
    array's outermost base is aligned """
    k = next((i for i, (n, op) in enumerate(ops) if op in '=X'), len(ops))
    i_bp = sum(n for n, op in ops[:k] if op == 'I')
    s = start + sum(n for n, op in ops[:k] if op == 'D') - i_bp
    if k == 0 or k == len(ops) or s < 0:
        return start, ops
    return s, [(1, '=' if a == b else 'X') for a, b in zip(array[:i_bp], target[s:s + i_bp])] + ops[k:]


def anchor_alignment(array, target_length, reference=False, motif=TTAGGG):
    """ align array (read from its outer edge inward as the motif's unit, CCCTAA) to the unit's repeat, (CCCTAA)n, of
    target_length (a whole number of units) with the outer ends pinned: (repeat offset r, target start, ops outer end
    first, target bp used from the start).  The array starts at r, the offset that matches best over its outermost
    PHASE_BP.  The reference's starts a repeat unit further in, at r + the motif's reference inset, r being the offset
    whose alignment of those bp has the fewest edits and, of those, the fewest before its first aligned pair, and
    starts with an aligned pair """
    target = motif.unit * (target_length // motif.p)
    if reference:
        def cost(r):
            a = edlib.align(array[:PHASE_BP], target[r:r + 2 * PHASE_BP], mode='SHW', task='path')
            return a['editDistance'], leading_gap(cigar_ops(a['cigar'])), r
        r = min(cost(r) for r in range(motif.p))[2]
        start = r + motif.reference_inset
    else:
        best = None
        for r in range(motif.p):
            m = sum(1 for a, b in zip(array[:PHASE_BP], target[r:r + PHASE_BP]) if a == b)
            if best is None or m > best[0]:
                best = (m, r)
        r = start = best[1]
    ops = cigar_ops(edlib.align(array, target[start:], mode='SHW', task='path')['cigar'])
    if reference:
        start, ops = pin_outer_base(array, target, start, ops)
    return r, start, ops, sum(n for n, op in ops if op in '=XD')


def target_name(minigraph_event, ref_contig, arm):
    """ the name of the target of reference contig ref_contig's arm (p or q), with any _sub_ in tel_<ref_contig> written
    _sub-, as hal2vg would read tel_X_sub_0_100p and tel_X_sub_0_100q as the same subpath of X (and tel_sub_0_100p as
    one of tel) """
    return 'id={}|{}{}'.format(minigraph_event, ('tel_' + ref_contig).replace('_sub_', '_sub-'), arm)


def read_paf_strands(paf_path, target_prefix):
    """ one pass over the PAF.  Returns (queries, n_lines, n_bad, n_to_targets): queries maps each query to [bp
    aligned on +, bp aligned on -, its length, whether all its lines parse and agree on its length, {reference
    contig: [bp aligned on +, bp aligned on -]}] (by the lines' rc:Z: tags, None for lines without one); n_bad counts
    the lines that don't parse, n_to_targets those already aligned to a target (a name starting with target_prefix,
    bytes) """
    queries = {}
    n_lines = n_bad = n_to_targets = 0
    with open(paf_path, 'rb') as f:
        for line in f:
            n_lines += 1
            t = line.split(b'\t', 6)
            if len(t) < 2:
                n_bad += 1
                continue
            q = queries.get(t[0])
            if q is None:
                q = queries[t[0]] = [0, 0, None, True, {}]
            try:
                ql, qs, qe = int(t[1]), int(t[2]), int(t[3])
                good = len(t) == 7 and t[4] in (b'+', b'-') and 0 <= qs < qe <= ql
            except (ValueError, IndexError):
                good = False
            if not good:
                n_bad += 1
                q[3] = False
                continue
            if t[5].startswith(target_prefix):
                n_to_targets += 1
            strand = 0 if t[4] == b'+' else 1
            q[strand] += qe - qs
            rc = t[6].find(b'\trc:Z:')
            ref_contig = None
            if rc >= 0:
                rc_end = t[6].find(b'\t', rc + 6)
                ref_contig = t[6][rc + 6:rc_end if rc_end >= 0 else None].rstrip() or None
            by_ref = q[4].get(ref_contig)
            if by_ref is None:
                by_ref = q[4][ref_contig] = [0, 0]
            by_ref[strand] += qe - qs
            if q[2] is None:
                q[2] = ql
            elif q[2] != ql:
                q[3] = False
    for q in queries.values():
        q[4] = {r if r is None else r.decode('utf-8', 'surrogateescape'): bp for r, bp in q[4].items()}
    return ({name.decode('utf-8', 'surrogateescape'): q for name, q in queries.items()}, n_lines, n_bad,
            n_to_targets)


def contig_key(q, own_name=None):
    """ the reference contig whose targets a contig's arrays are anchored to, given its read_paf_strands() entry q:
    own_name if given (a contig of the reference), otherwise the rc:Z: tag of a strict majority of the bp its tagged
    lines align.  Returns (it, None), or (None, why not) """
    if own_name is not None:
        return own_name, None
    tagged = {r: sum(bp) for r, bp in q[4].items() if r is not None}
    if not tagged:
        return None, 'no_reference_contig'
    key = max(tagged, key=tagged.get)
    if 2 * tagged[key] <= sum(tagged.values()):
        return None, 'ambiguous_reference_contig'
    return key, None


def write_trimmed_paf(paf_path, out, zones):
    """ write the PAF to the binary file out, each line of a query in zones (bytes) trimmed to its zone, [lo, hi), and
    every other line unchanged.  The queries in zones must be ones whose lines all parse (read_paf_strands).  Returns
    (n_trimmed, n_dropped, n_untrimmable, whether the last line written ended with a line break); an untrimmable line
    (no cigar, or one that doesn't fit its coordinates) is written unchanged """
    n_trimmed = n_dropped = n_untrimmable = 0
    ended = True
    with open(paf_path, 'rb') as f:
        for line in f:
            tab = line.find(b'\t')
            zone = zones.get(line[:tab]) if tab > 0 else None
            if zone is not None:
                t = line.split(b'\t', 4)
                zone = None if zone[0] <= int(t[2]) and int(t[3]) <= zone[1] else zone   # within it: unchanged
            if zone is not None:
                t = line.rstrip(b'\n').decode('utf-8', 'surrogateescape').split('\t')
                try:
                    trimmed = trim_paf(t, *zone)
                except (ValueError, IndexError):
                    trimmed = t
                    n_untrimmable += 1
                if trimmed is None:
                    n_dropped += 1
                    continue
                if trimmed is not t:
                    n_trimmed += 1
                    line = ('\t'.join(trimmed) + '\n').encode('utf-8', 'surrogateescape')
            out.write(line)
            ended = line.endswith(b'\n')
    return n_trimmed, n_dropped, n_untrimmable, ended


def write_minigraph_fasta(in_path, out_path, records):
    """ the (plain or gzipped) minigraph fasta, decompressed, then the (name, sequence) records """
    with open_maybe_gzip(in_path) as f, open(out_path, 'wb') as out:
        last = b'\n'
        while True:
            block = f.read(1 << 22)
            if not block:
                break
            out.write(block)
            last = block[-1:]
        if records and last != b'\n':
            out.write(b'\n')
        for name, seq in records:
            out.write('>{}\n{}\n'.format(name, seq).encode())


class InlineExecutor:
    """ runs each task when it is submitted, here: the pool for one thread """
    def submit(self, fn, *args):
        future = Future()
        try:
            future.set_result(fn(*args))
        except Exception as e:
            future.set_exception(e)
        return future

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def process_pool(threads, tasks, input_bytes):
    """ up to threads processes for tasks tasks and input_bytes of fasta (one per WORKER_BYTES, at most MAX_WORKERS),
    or the InlineExecutor if that comes to one """
    workers = min(int(threads), tasks, MAX_WORKERS, max(1, input_bytes // WORKER_BYTES))
    if workers <= 1:
        return InlineExecutor()
    # spawned, not forked: the caller is a toil worker, which can have threads of its own.  Each process imports the
    # caller's __main__ (toil's worker), hence the cap
    return ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn'))


# what anchor_telomeres() returns per anchored contig end
ANCHOR_KEYS = ['event', 'contig', 'end', 'orientation', 'arm', 'outer_bp', 'array_bp', 'offset', 'identity', 'strand',
               'q_start', 'q_end', 'target', 'target_length', 't_start', 't_end']


def anchor_line(d, alignment, target, target_length):
    """ the PAF line anchoring contig end d to its target, given its anchor_alignment(); fills in d's offset,
    identity, strand, target and target coordinates.  It is tagged with d's reference contig """
    r, start, ops, used = alignment
    if d['arm'] == 'p':
        ts, te = start, start + used
    else:   # the target is the reverse complement, so the alignment, outer end first, is read backwards
        ts, te, ops = target_length - start - used, target_length - start, ops[::-1]
    strand = '+' if (d['end'] == '5p') == (d['arm'] == 'p') else '-'
    matches, block = sum(n for n, op in ops if op == '='), sum(n for n, op in ops)
    d.update(offset=r, identity=matches / block, strand=strand, target=target, target_length=target_length,
             t_start=ts, t_end=te)
    fields = [d['contig'], d['length'], d['q_start'], d['q_end'], strand, target, target_length, ts, te, matches, block,
              ANCHOR_MAPQ, 'cg:Z:' + cigar_string(ops), 'gi:f:{:.4f}'.format(matches / block),
              'rc:Z:' + d['reference_contig'], 'tp:A:P', CONSTRAINT_TAG]
    return '\t'.join(map(str, fields)) + '\n'


def anchor_telomeres(label, fasta_paths, paf_path, out_paf_path, out_minigraph_path,
                     minigraph_event='_MINIGRAPH_', min_array=500, max_outer=3000, min_run=72, threads=1,
                     reference=None, motifs=DEFAULT_MOTIFS, log_choice=None):
    """ anchor the telomere arrays of the contigs to targets tel_<reference contig>p and tel_<reference contig>q of the
    minigraph genome, a pair for each reference contig with anchors (see the top of this file).  label names the job
    in the log, if not None.  fasta_paths maps each event, minigraph_event included, to its fasta, plain or gzipped,
    whose headers are the PAF's names (id=<event>|<contig>); paf_path holds the contigs' alignments to the minigraph
    segments, tagged rc:Z:<reference contig> as cactus-graphmap writes them (a contig is anchored to the targets of
    the reference contig that most of it lies on).  The contigs of event reference, if given, are anchored to the
    targets of their own names, a repeat unit further in than the others.  motifs: the telomere motifs (G strand units,
    see check_motifs()) the job chooses its motif from, or the one it takes; log_choice, if given, is called with the
    log's line about that choice and whether it is a warning (for a log that is kept: toil drops the log of a worker
    whose job succeeds, which is where the lines here go).  Writes out_paf_path, the PAF with the anchored ends trimmed
    out of its lines and one anchor line (tagged ct:i:1) per anchored end appended, every other line unchanged; and
    out_minigraph_path, the minigraph fasta (plain) with the targets of the arms that have anchors appended.  If nothing
    is anchored, the two are copies of the inputs.  A contig end that can't be anchored is skipped, and counted in the
    log.  The fastas are scanned, and the anchors aligned, by up to threads processes (fewer for small inputs; spawned,
    so a script that calls this with threads > 1 needs the usual if __name__ == '__main__' guard).  Returns a dict per
    anchored contig end, in the order of the anchor lines """
    start_time = time.time()
    motifs = check_motifs(motifs)
    log_name = 'telomere anchors' if label is None else 'telomere anchors for {}'.format(label)
    mg_path = fasta_paths[minigraph_event]
    target_prefix = target_name(minigraph_event, '', '')   # that every target's name starts with
    genomes = sorted(event for event in fasta_paths if event != minigraph_event)
    reasons = collections.Counter()   # of the contig ends not anchored
    input_bytes = sum(os.path.getsize(fasta_paths[event]) for event in genomes if os.path.isfile(fasta_paths[event]))

    with process_pool(threads, len(genomes) + 1, input_bytes) as pool:
        # the pool scans the fastas, for every motif, while the PAF is read here
        scans = [pool.submit(scan_genome, event, fasta_paths[event], min_array, max_outer, min_run, motifs)
                 for event in genomes]
        mg_names = pool.submit(fasta_names, mg_path)
        queries, n_lines, n_bad, n_to_targets = read_paf_strands(paf_path,
                                                                 target_prefix.encode('utf-8', 'surrogateescape'))

        # the job's motif: the only one, or the one its contig ends choose
        scanned = []
        tally = [[0, 0, 0] for _ in motifs]
        for scan in scans:
            event, contigs, genome_tally, error = scan.result()
            if error is not None:
                logger.warning('{}: {}; its contigs are not anchored'.format(log_name, error))
                continue
            scanned.append((event, contigs))
            for t, genome_t in zip(tally, genome_tally or []):
                t[:] = [a + b for a, b in zip(t, genome_t)]
        chosen = 0
        if len(motifs) > 1:
            chosen, text, warn = choose_motif(motifs, tally, 2 * sum(len(contigs) for _, contigs in scanned))
            (logger.warning if warn else logger.info)('{}: {}'.format(log_name, text))
            if log_choice is not None:
                log_choice(text, warn)
        motif = make_motif(motifs[chosen]) if chosen is not None else None

        # contig ends with an array of it are candidates; the others are only counted
        candidates = collections.defaultdict(list)
        contig_length = {}
        duplicates = set()
        for event, contigs in scanned:
            for name, length, ends in contigs:
                if name in contig_length:
                    duplicates.add(name)
                contig_length[name] = length
                for end, arrays in ends:
                    if motif is None:
                        reasons['no_motif'] += 1
                    elif arrays[chosen] is None:
                        reasons['no_array'] += 1
                    else:
                        outer, L, status, array, at_limit, run = arrays[chosen]
                        candidates[name].append(dict(event=event, contig=name, length=length, end=end,
                                                     outer_bp=outer, array_bp=L, status=status, array=array,
                                                     at_scan_limit=at_limit, perfect_run=run))
        existing_targets = [name for name in mg_names.result() if name.startswith(target_prefix)]

        # which candidates to anchor, and to which reference contig's targets: the reference's contigs to their own
        # names', unless no PAF line has a rc:Z: tag (a PAF not made by cactus-graphmap, which anchors nothing)
        tagged = any(r is not None for q in queries.values() for r in q[4])
        reference_prefix = 'id={}|'.format(reference)
        anchored, degenerate = [], []
        for name, ds in candidates.items():
            q = queries.get(name)
            own_name = None
            if tagged and ds[0]['event'] == reference:
                own_name = name[len(reference_prefix):] if name.startswith(reference_prefix) else name
            key, unkeyed = contig_key(q, own_name) if q is not None else (None, None)
            # two arrays that meet (a contig made of telomere): anchor neither
            overlap = len(ds) == 2 and sum(d['outer_bp'] + d['array_bp'] for d in ds) > contig_length[name]
            for d in ds:
                if d['status'] != 'ok':
                    reason = 'short_array' if d['status'] == 'short' else 'array_beyond_max_outer'
                elif d['perfect_run'] < min_run:
                    reason = 'degenerate_array'
                    degenerate.append(d)
                elif name in duplicates:
                    reason = 'duplicate_name'
                elif q is None:
                    reason = 'no_paf_lines'
                elif not q[3]:
                    reason = 'bad_paf_line'
                elif key is None:
                    reason = unkeyed
                elif q[2] != d['length']:
                    reason = 'length_mismatch'
                elif overlap:
                    reason = 'arrays_overlap'
                else:
                    # oriented by its lines on its reference contig, or by all of them if none is tagged with it
                    plus, minus = q[4].get(key) or q[:2]
                    d['reference_contig'] = key
                    d['orientation'] = '+' if plus >= minus else '-'
                    d['arm'] = ARM[(d['end'], d['orientation'])]
                    anchored.append(d)
                    continue
                reasons[reason] += 1
        if n_to_targets or existing_targets:
            logger.warning('{}: the inputs already have targets ({} in the minigraph fasta{}, {} PAF lines to them); '
                           'nothing anchored'.format(log_name, len(existing_targets),
                                                     ', such as ' + min(existing_targets) if existing_targets else '',
                                                     n_to_targets))
            reasons['already_anchored'] += len(anchored)
            anchored = []

        # targets, by (reference contig, arm), and the zone of each anchored contig its other lines may keep
        target_length = {}
        for d in anchored:
            k = (d['reference_contig'], d['arm'])
            target_length[k] = max(target_length.get(k, 0), d['array_bp'] + TARGET_PAD)
        target_length = {k: T + (-T) % motif.p for k, T in target_length.items()}
        anchored.sort(key=lambda d: (d['event'], d['contig'], d['end']))
        zones = {}
        for d in anchored:
            ql, L = d['length'], d['array_bp']
            d['q_start'], d['q_end'] = ((d['outer_bp'], d['outer_bp'] + L) if d['end'] == '5p' else
                                        (ql - d['outer_bp'] - L, ql - d['outer_bp']))
            zone = zones.setdefault(d['contig'].encode('utf-8', 'surrogateescape'), [0, ql])
            if d['end'] == '5p':
                zone[0] = d['q_end']
            else:
                zone[1] = d['q_start']

        # the pool aligns the anchors while the PAF is trimmed here
        alignments = [pool.submit(anchor_alignment, d['array'], target_length[d['reference_contig'], d['arm']],
                                  d['event'] == reference, motif) for d in anchored]
        n_trimmed = n_dropped = n_untrimmable = 0
        if not anchored:
            shutil.copyfile(paf_path, out_paf_path)
        else:
            with open(out_paf_path, 'wb') as out:
                n_trimmed, n_dropped, n_untrimmable, ended = write_trimmed_paf(paf_path, out, zones)
                if not ended:
                    out.write(b'\n')
                for d, alignment in zip(anchored, alignments):
                    k = (d['reference_contig'], d['arm'])
                    out.write(anchor_line(d, alignment.result(), target_name(minigraph_event, *k),
                                          target_length[k]).encode('utf-8', 'surrogateescape'))

    write_minigraph_fasta(mg_path, out_minigraph_path, [(target_name(minigraph_event, *k),
                                                         motif.tel_unit[k[1]] * (T // motif.p))
                                                        for k, T in sorted(target_length.items())])

    missing = set(queries) - set(contig_length)
    if missing:
        logger.warning('{}: {} PAF queries are not in the fastas, such as {}'.format(
            log_name, len(missing), sorted(missing)[0]))
    if n_bad or n_untrimmable:
        logger.warning('{}: {} PAF lines could not be parsed, and {} could not be trimmed'.format(
            log_name, n_bad, n_untrimmable))
    if duplicates:
        logger.warning('{}: {} contig names are in the fastas more than once, such as {}'.format(
            log_name, len(duplicates), sorted(duplicates)[0]))
    if degenerate:
        d = max(degenerate, key=lambda d: d['array_bp'])
        logger.info('{}: {} arrays without {}bp of perfect repeat are not anchored (interstitial telomeric repeats '
                    'where a contig breaks, or not {} telomeres), the longest of them {}bp at {} {} (its longest '
                    'perfect run {}bp)'.format(log_name, len(degenerate), min_run, motif.g, d['array_bp'], d['contig'],
                                               d['end'], d['perfect_run']))
    at_limit = [d for d in anchored if d['at_scan_limit']]
    if at_limit:
        logger.warning('{}: {} arrays reach the {}bp scan limit, such as {}'.format(
            log_name, len(at_limit), max_outer + SCAN, at_limit[0]['contig']))
    n_arm = collections.Counter(d['arm'] for d in anchored)
    logger.info('{}: {} contigs in {} genomes; anchored {} p ends and {} q ends ({}), {} past outer junk; not '
                'anchored: {}; {} of {} PAF lines trimmed, {} dropped; {:.1f}s'.format(
                    log_name, len(contig_length), len(genomes), n_arm['p'], n_arm['q'],
                    targets_text(target_length.values()), sum(1 for d in anchored if d['outer_bp'] > 0),
                    ', '.join('{} {}'.format(r, reasons[r]) for r in REASONS if reasons[r]) or 'none', n_trimmed,
                    n_lines, n_dropped, time.time() - start_time))

    return [{k: d[k] for k in ANCHOR_KEYS} for d in anchored]


def targets_text(lengths):
    """ how many targets of these lengths, and their range, for the log """
    lengths = sorted(lengths)
    if not lengths:
        return 'no targets'
    return '{} target{} of {}bp'.format(len(lengths), '' if len(lengths) == 1 else 's',
                                         lengths[0] if lengths[0] == lengths[-1] else
                                         '{}-{}'.format(lengths[0], lengths[-1]))


def anchor_summary(anchors):
    """ a line for the log about anchor_telomeres()'s anchors """
    if not anchors:
        return 'no telomere anchors'
    n_arm = collections.Counter(a['arm'] for a in anchors)
    return '{} telomere anchors ({} p, {} q) to {}, {} past outer junk, identity {:.3f}-{:.3f}'.format(
        len(anchors), n_arm['p'], n_arm['q'], targets_text({a['target']: a['target_length'] for a in anchors}.values()),
        sum(1 for a in anchors if a['outer_bp'] > 0), min(a['identity'] for a in anchors),
        max(a['identity'] for a in anchors))
