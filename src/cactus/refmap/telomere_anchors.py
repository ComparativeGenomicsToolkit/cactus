#!/usr/bin/env python3

"""
Telomere anchors for a pangenome chromosome: every contig's telomere array is aligned to one synthetic telomere per
chromosome arm, which joins the minigraph genome, so that cactus aligns the telomeres of all haplotypes to each other
whatever subtelomere they hang off.

Minigraph's mapping never covers telomere arrays, not even the reference's own, and BAR only reaches a telomere from
the inner end of its block.  So without anchors the arrays of different haplotypes stay unaligned, or are smeared over
a short reference array, and clipping removes them, along with any subtelomere block (often one swapped in from
another chromosome's end) that only the telomere tied to the graph.  Per chromosome, just before cactus_consolidated:

1. Each contig end is read from the end inward, on the strand that reads CCCTAA, for a telomere array.  Its core is
   the first 100bp window (starting at most max_outer + 100bp from the end) that is at least half telomeric hexamers
   (CCCTAA and its common variants) and whose walk inward, in 100bp windows with one failing window allowed, reaches
   min_array.  Its outer edge is the outermost telomeric k-mer reachable from the core in steps of at most 40bp, so
   degenerate telomere stays in the array and other junk (such as (TTAAAG)n) stays out.  An edge under two repeat
   units from the contig end is moved to the end; an array whose edge is further in than max_outer is not anchored.
   Nor is an array without min_run bp of perfect (CCCTAA)n: telomerase makes long perfect runs, while an interstitial
   telomeric repeat (ITS), which a contig that breaks next to it would leave at its end, is degenerate, and anchoring
   it would join the telomere to the middle of the chromosome.  Over CHM13 and the chromosome ends of 21 HPRC
   haplotypes, the longest perfect run of an ITS of at least 500bp was 54bp, and the shortest of a telomere 84bp (91%
   are over 300bp); 5 of the 998 arrays at the contig ends of 22 other, mostly not T2T, HPRC haplotypes fell under
   72bp, all of them telomeres.  A plant's (CCCTAAA)n has no run at all.
2. The contig's orientation is the PAF strand that holds most of its aligned bases, and puts each array on the p or q
   arm.  A contig without PAF lines can't be oriented, and is not anchored.  Nor is a contig whose PAF lines (by their
   rc:Z: tags, the reference contig of the minigraph segments) mostly lie on a reference contig other than chrom: in a
   job that holds several reference contigs (cactus-pangenome's whole-genome --noSplit job, named "all", or the
   chrOther bin of the small ones) every telomere would otherwise go to the same target per arm.
3. Target id=<minigraph event>|tel_<chrom>p is (CCCTAA)n with the chromosome end at position 0, and tel_<chrom>q is
   its reverse complement, (TTAGGG)n with the end at its last base.  Each is the arm's longest anchored array + 1000bp
   long, rounded up to a multiple of 6, and is only made for an arm with an anchor.
4. Each array is aligned to its target with the outer ends pinned together, at the repeat offset (0-5) that matches
   best over the outermost 600bp, and the inner end free (edlib, mode SHW): a slip in the repeat costs an indel
   rather than putting the rest of the array out of phase.  These anchor lines are tagged ct:i:1, which makes them
   constraints for cactus_consolidated, kept whatever chain they end up in.
5. The existing PAF lines of an anchored contig are trimmed out of its anchored ends (the array and anything outside
   it), cutting their cigars, so that nothing else pairs the telomere up some other way.  A line left with no aligned
   pair is dropped.  Every other line is written unchanged.

On 40 HPRC haplotypes x 8 chromosomes (with minigraph end events), anchoring kept all telomere sequence in the clipped
graph, against 28.5% without it, and 98% of the swapped subtelomere blocks.

Copyright (C) 2026 by Glenn Hickey

Released under the MIT license, see LICENSE.txt
"""

import collections
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
PARTIAL = 11     # an outer edge this close to the contig end is moved to the end (under two repeat units)
PHASE_BP = 600   # the anchor's repeat offset is the one that matches best over this many outer bp
TARGET_PAD = 1000  # a target is this much longer than the longest array anchored to it
WORKER_BYTES = 1 << 30  # a process per this many bytes of fasta to scan
MAX_WORKERS = 16   # the scan is bound by reading the fastas well before this many processes

# telomeric hexamers read from a chromosome end inward: CCCTAA, and the reverse complements of the TTAGGG variants
# TCAGGG, TGAGGG, TTGGGG and TTCGGG
MOTIFS = ['CCCTAA', 'CCCTGA', 'CCCTCA', 'CCCCAA', 'CCCGAA']
MOTIF_RE = re.compile('|'.join(MOTIFS))   # matched left to right without overlaps
MOTIF_START_RE = re.compile('(?=(?:{}))'.format('|'.join(MOTIFS)))
# the k-mers that keep the outer edge going: the six 5-mers of (CCCTAA)n, and the variant hexamers not starting with one
TEL5 = sorted(set(('CCCTAA' * 2)[i:i + 5] for i in range(6)))
EDGE_KMER_RE = re.compile('(?=(?:{}))'.format('|'.join(TEL5 + [m for m in MOTIFS if m[:5] not in TEL5])))
PERFECT_RE = re.compile('(?:CCCTAA)+')

# a contig end's arm, from the end (5' or 3') and the contig's orientation relative to the reference
ARM = {('5p', '+'): 'p', ('5p', '-'): 'q', ('3p', '+'): 'q', ('3p', '-'): 'p'}
TEL_UNIT = {'p': 'CCCTAA', 'q': 'TTAGGG'}
ANCHOR_MAPQ = 60
CONSTRAINT_TAG = 'ct:i:1'

COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
CIGAR_OPS = re.compile(r'(\d+)([=XIDM])')
WHITESPACE = b'\r\n\t '

# why a contig end was not anchored, in the order the log lists them
REASONS = ['no_array', 'short_array', 'array_beyond_max_outer', 'degenerate_array', 'duplicate_name', 'no_paf_lines',
           'bad_paf_line', 'other_chromosome', 'length_mismatch', 'arrays_overlap', 'already_anchored']


def reverse_complement(s):
    return s.translate(COMPLEMENT)[::-1]


def telomere_fraction(w):
    """ fraction of w covered by telomeric hexamers, taken left to right without overlaps """
    return 6 * len(MOTIF_RE.findall(w)) / max(1, len(w))


def array_walk(s, x0):
    """ bp of telomere array in s from x0 inward: whole windows that are at least half telomeric, with one failing
    window allowed between two that pass """
    x = run = x0
    fails = 0
    while x + WINDOW <= len(s):
        if telomere_fraction(s[x:x + WINDOW]) >= 0.5:
            x += WINDOW
            run = x
            fails = 0
        else:
            fails += 1
            if fails > 1:
                break
            x += WINDOW
    return run - x0


def outer_edge(s, c):
    """ the array's outer edge in s, from c (a telomeric hexamer of its core) outward: the outermost telomeric k-mer
    reachable in steps of at most OUTER_GAP; 0 if that is within PARTIAL of the end """
    e = c
    for o in reversed([m.start() for m in EDGE_KMER_RE.finditer(s, 0, c + 5) if m.start() < c]):
        if e - o > OUTER_GAP:
            break
        e = o
    return 0 if e <= PARTIAL else e


def array_at_end(s, min_array, max_outer):
    """ s: a contig end read from the end inward as CCCTAA, upper case.  (outer_bp, array_bp, status), the array being
    s[outer_bp : outer_bp + array_bp].  status is ok; beyond (its edge is further in than max_outer: not to be
    anchored); short (no core reaches min_array: this is the first telomeric run); or none (no telomeric window) """
    short = None
    for x0 in range(0, max_outer + WINDOW + 1, WINDOW):
        if telomere_fraction(s[x0:x0 + WINDOW]) < 0.5:
            continue
        L = array_walk(s, x0)
        e = outer_edge(s, MOTIF_START_RE.search(s, x0).start())
        if L < min_array:
            if short is None:
                short = (e, x0 + L - e, 'short')
            continue
        return e, x0 + L - e, ('ok' if e <= max_outer else 'beyond')
    return short or (0, 0, 'none')


def perfect_run(array):
    """ bp of the longest run of perfect CCCTAA repeat in array """
    return max((len(m) for m in PERFECT_RE.findall(array)), default=0)


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


def scan_genome(event, path, min_array, max_outer):
    """ (event, contigs, error) for one genome's fasta.  contigs: (name, length, ends), ends being its 5' then 3' end as
    (end, outer_bp, array_bp, status, array, reaches_scan_limit, perfect_run), array the ok array's sequence read
    from the end inward as CCCTAA, and perfect_run its perfect_run().  A fasta that can't be read gives (event, None,
    the error) """
    contigs = []
    try:
        for name, length, head, tail in fasta_ends(path, max_outer + SCAN):
            ends = []
            for end, s in (('5p', head.decode('latin-1').upper()),
                           ('3p', reverse_complement(tail.decode('latin-1').upper()))):
                outer, L, status = array_at_end(s, min_array, max_outer)
                array = s[outer:outer + L] if status == 'ok' else None
                ends.append((end, outer, L, status, array, outer + L >= len(s), perfect_run(array) if array else 0))
            contigs.append((name, length, ends))
    except (OSError, EOFError, zlib.error) as e:
        return event, None, '{}: {}'.format(path, e)
    return event, contigs, None


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
    for n, op in cigar_ops(t[cg_i][5:]):
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
    t[9] = str(sum(p[1] for p in pieces if p[0] == '='))
    t[10] = str(sum(p[1] for p in pieces))
    return t


def anchor_alignment(array, target_length):
    """ align array (read from its outer edge inward as CCCTAA) to (CCCTAA)n of target_length with the outer ends
    pinned: (repeat offset r, ops outer end first, target bp used from r) """
    target = TEL_UNIT['p'] * (target_length // 6)
    best = None
    for r in range(6):
        m = sum(1 for a, b in zip(array[:PHASE_BP], target[r:r + PHASE_BP]) if a == b)
        if best is None or m > best[0]:
            best = (m, r)
    r = best[1]
    ops = cigar_ops(edlib.align(array, target[r:], mode='SHW', task='path')['cigar'])
    return r, ops, sum(n for n, op in ops if op in '=XD')


def read_paf_strands(paf_path, target_names, chrom):
    """ one pass over the PAF.  Returns (queries, n_lines, n_bad, n_to_targets): queries maps each query to [bp
    aligned on +, bp aligned on -, its length, whether all its lines parse and agree on its length, bp aligned to
    reference contig chrom, bp aligned to other reference contigs] (by the lines' rc:Z: tags); n_bad counts the lines
    that don't parse, n_to_targets those already aligned to one of target_names (bytes) """
    queries = {}
    n_lines = n_bad = n_to_targets = 0
    chrom = chrom.encode('utf-8', 'surrogateescape')
    with open(paf_path, 'rb') as f:
        for line in f:
            n_lines += 1
            t = line.split(b'\t', 6)
            if len(t) < 2:
                n_bad += 1
                continue
            q = queries.get(t[0])
            if q is None:
                q = queries[t[0]] = [0, 0, None, True, 0, 0]
            try:
                ql, qs, qe = int(t[1]), int(t[2]), int(t[3])
                good = len(t) == 7 and t[4] in (b'+', b'-') and 0 <= qs < qe <= ql
            except (ValueError, IndexError):
                good = False
            if not good:
                n_bad += 1
                q[3] = False
                continue
            if t[5] in target_names:
                n_to_targets += 1
            q[0 if t[4] == b'+' else 1] += qe - qs
            rc = t[6].find(b'\trc:Z:')
            if rc >= 0:
                rc_end = t[6].find(b'\t', rc + 6)
                q[4 if t[6][rc + 6:rc_end if rc_end >= 0 else None].rstrip() == chrom else 5] += qe - qs
            if q[2] is None:
                q[2] = ql
            elif q[2] != ql:
                q[3] = False
    return ({name.decode('utf-8', 'surrogateescape'): q for name, q in queries.items()}, n_lines, n_bad,
            n_to_targets)


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


def anchor_line(d, alignment, target, target_length, chrom):
    """ the PAF line anchoring contig end d to its arm's target, given its anchor_alignment(); fills in d's offset,
    identity, strand, target and target coordinates """
    r, ops, used = alignment
    if d['arm'] == 'p':
        ts, te = r, r + used
    else:   # the target is the reverse complement, so the alignment, outer end first, is read backwards
        ts, te, ops = target_length - r - used, target_length - r, ops[::-1]
    strand = '+' if (d['end'] == '5p') == (d['arm'] == 'p') else '-'
    matches, block = sum(n for n, op in ops if op == '='), sum(n for n, op in ops)
    d.update(offset=r, identity=matches / block, strand=strand, target=target, target_length=target_length,
             t_start=ts, t_end=te)
    fields = [d['contig'], d['length'], d['q_start'], d['q_end'], strand, target, target_length, ts, te, matches, block,
              ANCHOR_MAPQ, 'cg:Z:' + cigar_string(ops), 'gi:f:{:.4f}'.format(matches / block), 'rc:Z:' + chrom,
              'tp:A:P', CONSTRAINT_TAG]
    return '\t'.join(map(str, fields)) + '\n'


def anchor_telomeres(chrom, fasta_paths, paf_path, out_paf_path, out_minigraph_path,
                     minigraph_event='_MINIGRAPH_', min_array=500, max_outer=3000, min_run=72, threads=1):
    """ anchor the telomere arrays of chromosome chrom's contigs to targets tel_<chrom>p and tel_<chrom>q of the
    minigraph genome (see the top of this file).  fasta_paths maps each event, minigraph_event included, to its
    fasta, plain or gzipped, whose headers are the PAF's names (id=<event>|<contig>); paf_path holds the contigs'
    alignments to the minigraph segments, tagged rc:Z:<reference contig> as cactus-graphmap writes them (only contigs
    mostly on reference contig chrom are anchored).  Writes out_paf_path, the PAF with the anchored ends trimmed out
    of its lines and one anchor line (tagged ct:i:1) per anchored end appended, every other line unchanged; and
    out_minigraph_path, the minigraph fasta (plain) with the targets of the arms that have anchors appended.  If
    nothing is anchored, the two are copies of the inputs.  A contig end that can't be anchored is skipped, and
    counted in the log.  The fastas are scanned, and the anchors aligned, by up to threads processes (fewer for
    small inputs; spawned, so a script that calls this with threads > 1 needs the usual if __name__ == '__main__'
    guard).  Returns a dict per anchored contig end, in the order of the anchor lines """
    start_time = time.time()
    mg_path = fasta_paths[minigraph_event]
    target_names = {arm: 'id={}|tel_{}{}'.format(minigraph_event, chrom, arm) for arm in 'pq'}
    genomes = sorted(event for event in fasta_paths if event != minigraph_event)
    reasons = collections.Counter()   # of the contig ends not anchored
    input_bytes = sum(os.path.getsize(fasta_paths[event]) for event in genomes if os.path.isfile(fasta_paths[event]))

    with process_pool(threads, len(genomes) + 1, input_bytes) as pool:
        # the pool scans the fastas while the PAF is read here
        scans = [pool.submit(scan_genome, event, fasta_paths[event], min_array, max_outer) for event in genomes]
        mg_names = pool.submit(fasta_names, mg_path)
        queries, n_lines, n_bad, n_to_targets = read_paf_strands(
            paf_path, set(name.encode('utf-8', 'surrogateescape') for name in target_names.values()), chrom)

        # contig ends with an array are candidates; the others are only counted
        candidates = collections.defaultdict(list)
        contig_length = {}
        duplicates = set()
        for scan in scans:
            event, contigs, error = scan.result()
            if error is not None:
                logger.warning('telomere anchors for {}: {}; its contigs are not anchored'.format(chrom, error))
                continue
            for name, length, ends in contigs:
                if name in contig_length:
                    duplicates.add(name)
                contig_length[name] = length
                for end, outer, L, status, array, at_limit, run in ends:
                    if L > 0:
                        candidates[name].append(dict(event=event, contig=name, length=length, end=end,
                                                     outer_bp=outer, array_bp=L, status=status, array=array,
                                                     at_scan_limit=at_limit, perfect_run=run))
                    else:
                        reasons['no_array'] += 1
        existing_targets = set(target_names.values()) & set(mg_names.result())

        # which candidates to anchor
        anchored, degenerate = [], []
        for name, ds in candidates.items():
            q = queries.get(name)
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
                elif q[4] <= q[5]:
                    reason = 'other_chromosome'
                elif q[2] != d['length']:
                    reason = 'length_mismatch'
                elif overlap:
                    reason = 'arrays_overlap'
                else:
                    d['orientation'] = '+' if q[0] >= q[1] else '-'
                    d['arm'] = ARM[(d['end'], d['orientation'])]
                    anchored.append(d)
                    continue
                reasons[reason] += 1
        if n_to_targets or existing_targets:
            logger.warning('telomere anchors for {}: the inputs already have {} ({} PAF lines); nothing '
                           'anchored'.format(chrom, ' and '.join(sorted(existing_targets)) or 'targets', n_to_targets))
            reasons['already_anchored'] += len(anchored)
            anchored = []

        # targets, and the zone of each anchored contig its other lines may keep
        target_length = {}
        for arm in 'pq':
            Ls = [d['array_bp'] for d in anchored if d['arm'] == arm]
            if Ls:
                target_length[arm] = max(Ls) + TARGET_PAD
                target_length[arm] += (6 - target_length[arm] % 6) % 6
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
        alignments = [pool.submit(anchor_alignment, d['array'], target_length[d['arm']]) for d in anchored]
        n_trimmed = n_dropped = n_untrimmable = 0
        if not anchored:
            shutil.copyfile(paf_path, out_paf_path)
        else:
            with open(out_paf_path, 'wb') as out:
                n_trimmed, n_dropped, n_untrimmable, ended = write_trimmed_paf(paf_path, out, zones)
                if not ended:
                    out.write(b'\n')
                for d, alignment in zip(anchored, alignments):
                    out.write(anchor_line(d, alignment.result(), target_names[d['arm']], target_length[d['arm']],
                                          chrom).encode('utf-8', 'surrogateescape'))

    write_minigraph_fasta(mg_path, out_minigraph_path, [(target_names[arm], TEL_UNIT[arm] * (target_length[arm] // 6))
                                                        for arm in sorted(target_length)])

    missing = set(queries) - set(contig_length)
    if missing:
        logger.warning('telomere anchors for {}: {} PAF queries are not in the fastas, such as {}'.format(
            chrom, len(missing), sorted(missing)[0]))
    if n_bad or n_untrimmable:
        logger.warning('telomere anchors for {}: {} PAF lines could not be parsed, and {} could not be trimmed'.format(
            chrom, n_bad, n_untrimmable))
    if duplicates:
        logger.warning('telomere anchors for {}: {} contig names are in the fastas more than once, such as {}'.format(
            chrom, len(duplicates), sorted(duplicates)[0]))
    if degenerate:
        d = max(degenerate, key=lambda d: d['array_bp'])
        logger.info('telomere anchors for {}: {} arrays without {}bp of perfect repeat are not anchored (interstitial '
                    'telomeric repeats where a contig breaks, or not TTAGGG telomeres), the longest of them {}bp at {} '
                    '{} (its longest perfect run {}bp)'.format(chrom, len(degenerate), min_run, d['array_bp'],
                                                               d['contig'], d['end'], d['perfect_run']))
    if reasons['other_chromosome'] and not anchored:
        logger.info('telomere anchors for {}: no contig with an array is mostly on reference contig {} by the rc:Z: '
                    'tags of its PAF lines, so nothing is anchored (as for a job of several reference contigs)'.format(
                        chrom, chrom))
    at_limit = [d for d in anchored if d['at_scan_limit']]
    if at_limit:
        logger.warning('telomere anchors for {}: {} arrays reach the {}bp scan limit, such as {}'.format(
            chrom, len(at_limit), max_outer + SCAN, at_limit[0]['contig']))
    n_arm = collections.Counter(d['arm'] for d in anchored)
    logger.info('telomere anchors for {}: {} contigs in {} genomes; anchored {} p ends and {} q ends (targets: {}), {} '
                'past outer junk; not anchored: {}; {} of {} PAF lines trimmed, {} dropped; {:.1f}s'.format(
                    chrom, len(contig_length), len(genomes), n_arm['p'], n_arm['q'],
                    ', '.join('{} {}bp'.format(arm, target_length[arm]) for arm in sorted(target_length)) or 'none',
                    sum(1 for d in anchored if d['outer_bp'] > 0),
                    ', '.join('{} {}'.format(r, reasons[r]) for r in REASONS if reasons[r]) or 'none', n_trimmed,
                    n_lines, n_dropped, time.time() - start_time))

    return [{k: d[k] for k in ANCHOR_KEYS} for d in anchored]


def anchor_summary(anchors):
    """ a line for the log about anchor_telomeres()'s anchors """
    if not anchors:
        return 'no telomere anchors'
    by_arm = {arm: [a for a in anchors if a['arm'] == arm] for arm in 'pq'}
    return '{} telomere anchors ({}), {} past outer junk, identity {:.3f}-{:.3f}'.format(
        len(anchors), ', '.join('{} {} to {}bp'.format(len(a), arm, a[0]['target_length'])
                                for arm, a in by_arm.items() if a),
        sum(1 for a in anchors if a['outer_bp'] > 0), min(a['identity'] for a in anchors),
        max(a['identity'] for a in anchors))
