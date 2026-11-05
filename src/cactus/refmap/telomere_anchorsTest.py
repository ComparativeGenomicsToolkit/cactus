#!/usr/bin/env python3

"""
Unit tests for the telomere anchors of pangenome chromosomes (telomere_anchors.py): array detection and its outer
edge, orientation and arm, the target, the anchor alignment and the trimming of existing PAF lines, on small
synthetic chromosomes, in both contig orientations.  Every run's outputs go through verify_outputs(), which reads the
sequences back and checks each anchor and each changed line against them, independently of the module.

These are fast and offline.
"""

import gzip
import os
import random
import re
import shutil
import tempfile
import unittest
import zlib
from unittest import mock

import edlib

from cactus.refmap import telomere_anchors
from cactus.refmap.telomere_anchors import anchor_alignment, anchor_telomeres, anchor_summary, array_at_end, \
    fasta_ends, perfect_run, pin_outer_base, trim_paf

MG = '_MINIGRAPH_'
CHROM = 'chrT'

# the verifier's own view of the telomere, rather than the module's: the six 5-mers of (CCCTAA)n and the TTAGGG
# variants (TCAGGG, TGAGGG, TTGGGG, TTCGGG) read from the chromosome end inward
TEL5 = set(('CCCTAA' * 2)[i:i + 5] for i in range(6))
VARIANTS = {'CCCTGA', 'CCCTCA', 'CCCCAA', 'CCCGAA'}
OUTER_GAP = 40
PARTIAL = 11


def rc(s):
    return s.translate(str.maketrans('ACGTacgtNn', 'TGCAtgcaNn'))[::-1]


def random_dna(rng, n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def ops_of(cigar):
    return [(int(n), op) for n, op in re.findall(r'(\d+)([=XIDM])', cigar)]


def merge_ops(ops):
    out = []
    for n, op in ops:
        if out and out[-1][1] == op:
            out[-1] = (out[-1][0] + n, op)
        else:
            out.append((n, op))
    return ''.join('{}{}'.format(n, op) for n, op in out)


def mirror_line(t):
    """ a PAF line (fields) for the reverse complement of its query """
    t = list(t)
    ql, qs, qe = int(t[1]), int(t[2]), int(t[3])
    t[2], t[3] = str(ql - qe), str(ql - qs)
    t[4] = '-' if t[4] == '+' else '+'
    return t


def reference_trim(t, lo, hi):
    """ the prototype's per-base trim of a PAF line to the query interval [lo, hi), which trim_paf() must equal """
    qs, qe, strand = int(t[2]), int(t[3]), t[4]
    lo, hi = max(qs, lo), min(qe, hi)
    if hi <= lo:
        return None
    if (lo, hi) == (qs, qe):
        return t
    cg_i = next(i for i, x in enumerate(t) if x.startswith('cg:Z:'))
    step = 1 if strand == '+' else -1
    q, tpos = (qs if strand == '+' else qe), int(t[7])
    keep = []   # (op, target base, query base) per step
    for n, op in ops_of(t[cg_i][5:]):
        for _ in range(n):
            qb = q if strand == '+' else q - 1
            if op in '=XM':
                if lo <= qb < hi:
                    keep.append((op, tpos, qb))
                tpos += 1
                q += step
            elif op == 'I':
                if lo <= qb < hi:
                    keep.append((op, None, qb))
                q += step
            else:
                keep.append((op, tpos, None))
                tpos += 1
    while keep and keep[0][0] in 'ID':
        keep.pop(0)
    while keep and keep[-1][0] in 'ID':
        keep.pop()
    aligned = [k for k in keep if k[0] in '=XM']
    if not aligned:
        return None
    t = list(t)
    t[2], t[3] = str(min(k[2] for k in aligned)), str(max(k[2] for k in aligned) + 1)
    t[7], t[8] = str(min(k[1] for k in aligned)), str(max(k[1] for k in aligned) + 1)
    t[cg_i] = 'cg:Z:' + merge_ops([(1, k[0]) for k in keep])
    t[9] = str(sum(1 for k in keep if k[0] == '='))
    t[10] = str(len(keep))
    return t


def paf_line(name, contig, q0, q1, seg_name, segment, t0, t1, strand='+', ref_contig=CHROM):
    """ a PAF line of contig[q0:q1] (reverse complemented on '-') against segment[t0:t1], with edlib's global
    alignment, on reference contig ref_contig (rc:Z:, left out if None) """
    query = contig[q0:q1].upper() if strand == '+' else rc(contig[q0:q1].upper())
    cigar = edlib.align(query, segment[t0:t1].upper(), mode='NW', task='path')['cigar']
    ops = ops_of(cigar)
    m, b = sum(n for n, op in ops if op == '='), sum(n for n, op in ops)
    return [name, str(len(contig)), str(q0), str(q1), strand, seg_name, str(len(segment)), str(t0), str(t1), str(m),
            str(b), '60', 'cg:Z:' + cigar, 'gi:f:{:.4f}'.format(m / b)] + \
        (['rc:Z:' + ref_contig] if ref_contig else []) + ['tp:A:P']


def write_fasta(path, records, gz=False, width=60):
    data = ''.join('>{}\n{}'.format(name, ''.join(seq[i:i + width] + '\n' for i in range(0, len(seq), width)))
                   for name, seq in records).encode()
    with (gzip.open if gz else open)(path, 'wb') as f:
        f.write(data)


def read_fasta(path):
    """ [(name, sequence)] of a plain or gzipped fasta """
    with open(path, 'rb') as f:
        data = f.read()
    if data[:2] == b'\x1f\x8b':
        data = gzip.decompress(data)
    out = []
    for rec in data.decode().split('>')[1:]:
        lines = rec.split('\n')
        out.append((lines[0].split()[0] if lines[0].split() else '', ''.join(l.strip() for l in lines[1:])))
    return out


class Chromosome:
    """ a synthetic chromosome: a random body tiled by three minigraph segments, the first starting and the last ending
    in 600bp of reference telomere, and contigs that copy the body with a few SNPs between telomere arrays """

    def __init__(self, seed=1, body_bp=12000):
        rng = random.Random(seed)
        self.body = random_dna(rng, body_bp)
        self.bounds = [0, body_bp // 3, 2 * body_bp // 3, body_bp]
        pieces = [self.body[a:b] for a, b in zip(self.bounds, self.bounds[1:])]
        pieces[0] = 'CCCTAA' * 100 + pieces[0]
        pieces[-1] = pieces[-1] + 'TTAGGG' * 100
        self.segments = [('id={}|s{}'.format(MG, i + 1), s) for i, s in enumerate(pieces)]

    def contig(self, name, p_array='', q_array='', outer5='', outer3='', reverse=False, array_line=False):
        """ (sequence, PAF lines as field lists) of outer5 + p_array + body + q_array + outer3, where the arrays and
        outer junk are given as read from their contig end inward (as CCCTAA...).  The lines align the body to the
        segments, the s1 and s3 lines reaching up to 600bp into the arrays, which the reference telomeres match;
        array_line adds a line within the p array.  reverse: the contig is reverse complemented, its lines mirrored """
        rng = random.Random(zlib.crc32(name.encode()))
        body = list(self.body)
        for i in range(250, len(body), 500):   # SNPs, so lines have mismatches
            body[i] = rng.choice([c for c in 'ACGT' if c != body[i]])
        body = ''.join(body)
        seq = outer5 + p_array + body + rc(q_array) + rc(outer3)
        P = len(outer5) + len(p_array)
        a5, a3 = min(len(p_array), 600), min(len(q_array), 600)
        (n1, s1), (n2, s2), (n3, s3) = self.segments
        b = self.bounds
        lines = [paf_line(name, seq, P - a5, P + b[1], n1, s1, 600 - a5, len(s1)),
                 paf_line(name, seq, P + b[1], P + b[2], n2, s2, 0, len(s2)),
                 paf_line(name, seq, P + b[2], P + b[3] + a3, n3, s3, 0, b[3] - b[2] + a3)]
        if array_line:
            x = len(outer5) + 600
            lines.append(paf_line(name, seq, x, x + 300, n1, s1, 0, 300))
        if reverse:
            seq, lines = rc(seq), [mirror_line(t) for t in lines]
        return seq, lines


def write_subproblem(work_dir, contigs, mg_records, gz=False, tag=''):
    """ contigs: [(event, name, sequence, lines)].  Writes one fasta per event (gzipped if gz, but always named .fa,
    as cactus_cons names them), the minigraph fasta and the PAF; returns (fasta_paths, paf_path) """
    fasta_paths, records = {}, {}
    for event, name, seq, _ in contigs:
        records.setdefault(event, []).append((name, seq))
    for event, recs in list(records.items()) + [(MG, mg_records)]:
        fasta_paths[event] = os.path.join(work_dir, '{}{}.fa'.format(event.replace('#', '.'), tag))
        write_fasta(fasta_paths[event], recs, gz=gz)
    paf_path = os.path.join(work_dir, 'in{}.paf'.format(tag))
    with open(paf_path, 'w') as f:
        for _, _, _, lines in contigs:
            for t in lines:
                f.write('\t'.join(t) + '\n')
    return fasta_paths, paf_path


def is_edge(s, e):
    """ s, a contig end read inward as CCCTAA: whether e is a telomere/junk boundary as the outer edge rule draws
    it: a telomeric k-mer starts at e and none within OUTER_GAP outside it, or e is the end and one starts within
    PARTIAL """
    def kmer_at(i):
        return s[i:i + 5] in TEL5 or s[i:i + 6] in VARIANTS
    if e == 0:
        return any(kmer_at(i) for i in range(PARTIAL + 1))
    return e > PARTIAL and kmer_at(e) and not any(kmer_at(i) for i in range(max(0, e - OUTER_GAP), e))


def verify_line(t, qseq, tseq, tag):
    """ coordinates, cigar lengths, = and X against the sequences, columns 10 and 11 """
    ql, qs, qe, strand, tl, ts, te = int(t[1]), int(t[2]), int(t[3]), t[4], int(t[6]), int(t[7]), int(t[8])
    assert len(qseq) == ql, '{}: query length {} != {}'.format(tag, ql, len(qseq))
    assert len(tseq) == tl, '{}: target length {} != {}'.format(tag, tl, len(tseq))
    assert 0 <= qs < qe <= ql and 0 <= ts < te <= tl and strand in '+-', '{}: coordinates out of range'.format(tag)
    cg = [x for x in t[12:] if x.startswith('cg:Z:')]
    assert len(cg) == 1 and re.fullmatch(r'(\d+[=XID])+', cg[0][5:]), '{}: bad cigar'.format(tag)
    ops = ops_of(cg[0][5:])
    assert sum(n for n, op in ops if op in '=XI') == qe - qs, '{}: cigar does not fit the query interval'.format(tag)
    assert sum(n for n, op in ops if op in '=XD') == te - ts, '{}: cigar does not fit the target interval'.format(tag)
    assert int(t[9]) == sum(n for n, op in ops if op == '=') and int(t[10]) == sum(n for n, op in ops), \
        '{}: columns 10/11 inconsistent with the cigar'.format(tag)
    # the query in target order: forward on '+', reverse complemented on '-'
    qa = qseq[qs:qe].upper() if strand == '+' else rc(qseq[qs:qe].upper())
    ta = tseq[ts:te].upper()
    i = j = 0
    for n, op in ops:
        if op == '=':
            assert qa[i:i + n] == ta[j:j + n], '{}: "=" run of {} at target {} is not a match'.format(tag, n, ts + j)
        elif op == 'X':
            assert all(a != b for a, b in zip(qa[i:i + n], ta[j:j + n])), \
                '{}: "X" run of {} at target {} holds a match'.format(tag, n, ts + j)
        i += n if op in '=XI' else 0
        j += n if op in '=XD' else 0


def verify_outputs(fasta_paths, paf_path, out_paf_path, out_mg_path, max_outer=3000, chrom=CHROM, min_run=72,
                   reference=None):
    """ the checks of the prototype's --check, raising AssertionError on the first failure.  Anchors: tagged ct:i:1,
    one per end at most, to tel_<chrom>p/q, of contigs mostly on chrom (by the input lines' rc:Z: tags), on the
    strand that holds most of the contig's aligned bases, pinned within 5bp of the target's outer end (the reference
    event's 6-11bp in, starting with an aligned pair), their query edge a telomere/junk boundary at most max_outer from
    the end, ending with a telomeric hexamer, holding min_run bp of perfect repeat, and the anchors of a contig not
    overlapping.  Other lines: those of unanchored contigs byte-identical, those of anchored ones the
    reference trim to the zone between the anchors (or dropped where it leaves nothing).  Anchors and changed lines:
    '=' real matches and 'X' real mismatches.  Minigraph fasta: the input's records unchanged, then exact repeat
    targets.  Returns (anchor lines, changed lines, dropped lines, unchanged lines) """
    with open(paf_path, 'rb') as f:
        in_lines = f.read().split(b'\n')
    with open(out_paf_path, 'rb') as f:
        out_lines = f.read().split(b'\n')
    assert in_lines[-1] == b'' or out_lines[-1] == b'', 'output PAF does not end with a line break'
    in_lines, out_lines = [l for l in in_lines if l], [l for l in out_lines if l]

    # each contig's aligned bp per strand, and on chrom or not, by the input lines
    strand_bp, chrom_bp = {}, {}
    for line in in_lines:
        t = line.decode().split('\t')
        bp = int(t[3]) - int(t[2])
        strand_bp.setdefault(t[0], [0, 0])[t[4] == '-'] += bp
        ref_contig = [x[5:] for x in t[12:] if x.startswith('rc:Z:')]
        if ref_contig:
            chrom_bp.setdefault(t[0], [0, 0])[ref_contig[0] != chrom] += bp

    def is_anchor(line):
        t = line.split(b'\t')
        return len(t) > 5 and b'|tel_' in t[5]
    anchors = [l.decode().split('\t') for l in out_lines if is_anchor(l)]
    zone, targets_of = {}, {}
    for t in anchors:
        name, ql, qs, qe, arm = t[0], int(t[1]), int(t[2]), int(t[3]), t[5][-1]
        # the end is the one the strand and arm imply; the edge checks below catch an anchor at the other end
        end = '5p' if (t[4] == '+') == (arm == 'p') else '3p'
        tag = '{} anchor {}'.format(name, end)
        assert 'ct:i:1' in t[12:], '{}: no ct:i:1'.format(tag)
        assert t[5] == 'id={}|tel_{}{}'.format(MG, chrom, arm) and arm in 'pq', '{}: target {}'.format(tag, t[5])
        assert t[4] == ('+' if strand_bp[name][0] >= strand_bp[name][1] else '-'), \
            '{}: strand {} against the contig\'s orientation'.format(tag, t[4])
        assert chrom_bp.get(name, [0, 0])[0] > chrom_bp.get(name, [0, 0])[1], '{}: not mostly on {}'.format(tag, chrom)
        T, ts, te = int(t[6]), int(t[7]), int(t[8])
        inset = 6 if reference is not None and name.startswith('id={}|'.format(reference)) else 0
        assert inset <= (ts if arm == 'p' else T - te) <= inset + 5, \
            '{}: not pinned {}-{}bp from the target outer end ({}-{} of {})'.format(tag, inset, inset + 5, ts, te, T)
        if inset:
            ops = ops_of([x for x in t[12:] if x.startswith('cg:Z:')][0][5:])
            assert (ops[0] if arm == 'p' else ops[-1])[1] in '=X', '{}: the reference starts with an indel'.format(tag)
        assert (qs if end == '5p' else ql - qe) <= max_outer, '{}: further inward than max_outer'.format(tag)
        z = zone.setdefault(name, [0, ql, set()])
        assert end not in z[2], '{}: two anchors at one end'.format(tag)
        z[2].add(end)
        z[0 if end == '5p' else 1] = qe if end == '5p' else qs
        assert t[5] not in targets_of.setdefault(name, set()), '{}: both ends anchored to {}'.format(name, t[5])
        targets_of[name].add(t[5])
    for name, z in zone.items():
        assert z[0] <= z[1], '{}: anchors overlap'.format(name)

    # input vs output, line by line
    to_verify = [(t, 'anchor') for t in anchors]
    others = [l for l in out_lines if not is_anchor(l)]
    assert out_lines[len(others):] == [l for l in out_lines if is_anchor(l)], 'anchors are not after the other lines'
    n_changed = n_dropped = n_same = 0
    k = 0
    for line in in_lines:
        name = line.split(b'\t')[0].decode()
        if name in zone:
            t = line.decode().split('\t')
            expected = reference_trim(t, zone[name][0], zone[name][1])
            if expected is None:
                n_dropped += 1
                continue
            assert int(expected[2]) >= zone[name][0] and int(expected[3]) <= zone[name][1], \
                '{}: line intrudes on an anchored end'.format(name)
            if expected is not t:
                n_changed += 1
                to_verify.append((expected, 'changed'))
            else:
                n_same += 1
            expected = '\t'.join(expected).encode()
        else:
            n_same += 1
            expected = line
        assert k < len(others) and others[k] == expected, 'output line {} is not the expected {}'.format(
            others[k][:100] if k < len(others) else None, expected[:100])
        k += 1
    assert k == len(others), 'extra output lines from {}'.format(others[k][:100] if k < len(others) else None)

    # sequences
    mg_in, mg_out = read_fasta(fasta_paths[MG]), read_fasta(out_mg_path)
    assert mg_out[:len(mg_in)] == mg_in, 'minigraph records changed'
    tseqs = dict(mg_out)
    targets = mg_out[len(mg_in):]
    assert sorted(name for name, _ in targets) == sorted(set(t[5] for t in anchors)), 'targets do not match anchors'
    for name, s in targets:
        arm = name[-1]
        assert len(s) % 6 == 0 and s == telomere_anchors.TEL_UNIT[arm] * (len(s) // 6), \
            '{} is not an exact repeat'.format(name)
    qseqs = {}
    for event, path in fasta_paths.items():
        if event != MG:
            qseqs.update(read_fasta(path))
    for t, kind in to_verify:
        tag = '{} {}:{}-{}({}) {}:{}-{}'.format(kind, t[0], t[2], t[3], t[4], t[5], t[7], t[8])
        verify_line(t, qseqs[t[0]], tseqs[t[5]], tag)
        if kind == 'anchor':
            qseq, qs, qe = qseqs[t[0]], int(t[2]), int(t[3])
            five = (t[4] == '+') == (t[5][-1] == 'p')
            e = qs if five else len(qseq) - qe
            s = (qseq[:e + 200] if five else rc(qseq[max(0, qe - 200):])).upper()
            assert is_edge(s, e), '{}: outer edge {}bp from the end is not a telomere/junk boundary'.format(tag, e)
            array = qseq[qs:qe].upper() if five else rc(qseq[qs:qe].upper())
            assert array[-6:] in VARIANTS | {'CCCTAA'}, \
                '{}: inner edge {} is not the end of a telomeric hexamer'.format(tag, array[-6:])
            run = max([len(r) for r in re.findall('(?:CCCTAA)+', qseq[qs:qe].upper())] +
                      [len(r) for r in re.findall('(?:TTAGGG)+', qseq[qs:qe].upper())])
            assert run >= min_run, '{}: longest perfect run {}bp'.format(tag, run)
    return len(anchors), n_changed, n_dropped, n_same


# outer junk read from the contig end inward, as (kind, bp) pieces, and the outer edge each one gives (from the
# prototype's tests): junk units read on the CCCTAA strand, degenerate telomere (chr16p HG00099#2's end), telomere
JUNK_UNITS = ['CTTTAA', 'TTAAAG', 'TTTAAC']
DEGENERATE = 'TCTACTAACCCTACCTACTAACCTAAACCTAACCTAACTAACCTACCCTAACCCTACCTACCCTAAACTAACTAACC'
JUNK_CASES = [[], [('junk', 5)], [('junk', 11)], [('junk', 12)], [('junk', 17)], [('junk', 30)], [('junk', 49)],
              [('junk', 60)], [('junk', 99)], [('junk', 120)], [('junk', 149)], [('junk', 150)], [('junk', 300)],
              [('junk', 1600)], [('junk', 2990)], [('junk', 2998)], [('junk', 3000)], [('junk', 3003)],
              [('junk', 3040)], [('junk', 3060)], [('junk', 5000)], [('degen', 80)], [('junk', 200), ('degen', 80)],
              [('tel', 300), ('junk', 500)], [('tel', 300), ('degen', 80), ('tel', 100)],
              [('tel', 400), ('degen', 300), ('tel', 100)]]


def degenerate_repeat(rng, n):
    """ n bp of an interstitial telomeric repeat like CHM13's: CCCTAA, perfect for at most 8 units at a time, broken
    by variants and indels """
    out = ''
    while len(out) < n:
        out += 'CCCTAA' * rng.randint(1, 8) + rng.choice(['CCCCAA', 'CCCTA', 'CCCAA', 'CCCTAAA', 'CCCTGA', 'CCCAAC'])
    return out[:n]


def noisy_telomere(rng, n, error):
    """ n bp of (CCCTAA)n with substitutions and deletions at rate error """
    out = []
    while len(out) < n:
        for c in 'CCCTAA':
            x = rng.random()
            if x < error:
                out.append(rng.choice('ACGT'))
            elif x >= error * 1.3:
                out.append(c)
    return ''.join(out[:n])


def junk_sequence(pieces, k):
    out = ''
    for kind, n in pieces:
        if kind == 'junk':
            out += (JUNK_UNITS[k % len(JUNK_UNITS)] * (n // 6 + 1))[:n]
        elif kind == 'degen':
            out += DEGENERATE[:n]
        else:
            out += ('CCCTAA' * (n // 6 + 1))[:n]
    return out


def expected_edge(pieces, s, max_outer=3000):
    """ s: the contig end read inward as CCCTAA.  The first telomeric k-mer starting at or after 5bp before the end of
    the innermost junk piece (junk ending in e.g. TTTAA reads on as TAACC into the array), 0 if none or within 11bp;
    None (not anchored) if beyond max_outer """
    lo = x = 0
    for kind, n in pieces:
        x += n
        if kind == 'junk':
            lo = max(0, x - 5)
    p = next(i for i in range(lo, len(s)) if s[i:i + 5] in TEL5 or s[i:i + 6] in VARIANTS)
    p = 0 if p <= PARTIAL else p
    return None if p > max_outer else p


class AnchorTestCase(unittest.TestCase):

    def setUp(self):
        self.work_dir = tempfile.mkdtemp(prefix='telomere_anchorsTest_')
        self.chrom = Chromosome()

    def tearDown(self):
        shutil.rmtree(self.work_dir)

    def run_anchors(self, contigs, gz=False, tag='', verify=True, mg_records=None, **kwargs):
        """ write the subproblem, anchor it and verify the outputs: (anchors, fasta_paths, paf, out paf, out mg) """
        fasta_paths, paf_path = write_subproblem(self.work_dir, contigs, mg_records or self.chrom.segments, gz=gz,
                                                 tag=tag)
        out_paf = os.path.join(self.work_dir, 'out{}.paf'.format(tag))
        out_mg = os.path.join(self.work_dir, 'out{}.mg.fa'.format(tag))
        anchors = anchor_telomeres(CHROM, fasta_paths, paf_path, out_paf, out_mg, minigraph_event=MG, **kwargs)
        if verify:
            verify_outputs(fasta_paths, paf_path, out_paf, out_mg, max_outer=kwargs.get('max_outer', 3000),
                           min_run=kwargs.get('min_run', 72), reference=kwargs.get('reference'))
        return anchors, fasta_paths, paf_path, out_paf, out_mg

    def assertSameFile(self, a, b, gunzip_a=False):
        with open(a, 'rb') as f:
            da = f.read()
        if gunzip_a and da[:2] == b'\x1f\x8b':
            da = gzip.decompress(da)
        with open(b, 'rb') as f:
            self.assertEqual(da, f.read())

    def paf_fields(self, path):
        with open(path) as f:
            return [l.rstrip('\n').split('\t') for l in f if l.strip()]

    def plain_contigs(self, reverse=lambda i: False):
        """ six haplotypes with arrays at both ends (some long, some short, one variant-rich), and a reference
        without arrays """
        contigs = []
        for i, (p, q) in enumerate([(3000, 2400), (1800, 4200), (6000, 600), (900, 3600), (2400, 2400), (4800, 1200)]):
            event = 'H{}#1'.format(i)
            name = 'id={}|{}'.format(event, CHROM)
            p_array = 'CCCTAA' * (p // 6)
            if i == 4:   # variant repeats and a slip
                p_array = ('CCCTAA' * 20 + 'CCCTGA' + 'CCCTAAA' + 'CCCCAA') * (p // 138)
            seq, lines = self.chrom.contig(name, p_array, 'CCCTAA' * (q // 6), reverse=reverse(i), array_line=i == 2)
            contigs.append((event, name, seq, lines))
        seq, lines = self.chrom.contig('id=REF|' + CHROM)
        contigs.append(('REF', 'id=REF|' + CHROM, seq, lines))
        return contigs


class TestArrayAtEnd(unittest.TestCase):

    def setUp(self):
        self.body = random_dna(random.Random(3), 2000)

    def test_plain(self):
        self.assertEqual(array_at_end('CCCTAA' * 500 + self.body, 500, 3000)[:3], (0, 3000, 'ok'))

    def test_partial_unit_at_the_end(self):
        # an edge under two repeat units from the end is moved to the end (the array's inner end is that of the last
        # hexamer starting in its last passing 100bp window: the window ends at 3000)
        self.assertEqual(array_at_end('TAA' + 'CCCTAA' * 500 + self.body, 500, 3000)[:3], (0, 3003, 'ok'))
        self.assertEqual(array_at_end('GGGGGGG' + 'CCCTAA' * 500 + self.body, 500, 3000)[:3], (0, 3001, 'ok'))
        self.assertEqual(array_at_end('GGGGGGGGGGGG' + 'CCCTAA' * 500 + self.body, 500, 3000)[:3], (12, 2988, 'ok'))

    def test_junk_outside(self):
        s = 'TTAAAG' * 50 + 'CCCTAA' * 500 + self.body
        self.assertEqual(array_at_end(s, 500, 3000)[:3], (300, 3000, 'ok'))
        # junk ending in TTTAA reads on as TAACC into the array
        s = 'CTTTAA' * 50 + 'CCCTAA' * 500 + self.body
        self.assertEqual(array_at_end(s, 500, 3000)[:3], (297, 3003, 'ok'))

    def test_degenerate_telomere_inside(self):
        s = DEGENERATE + 'CCCTAA' * 500 + self.body
        self.assertEqual(array_at_end(s, 500, 3000)[0], 0)
        s = 'TTAAAG' * 50 + DEGENERATE + 'CCCTAA' * 500 + self.body
        self.assertEqual(array_at_end(s, 500, 3000)[0], 304)

    def test_beyond_max_outer(self):
        # the core window starts within max_outer + 100bp, but its outer edge is further in than max_outer
        s = 'TTAAAG' * 507 + 'CCCTAA' * 500 + self.body
        self.assertEqual(array_at_end(s, 500, 3000)[:3], (3042, 2958, 'beyond'))
        self.assertEqual(array_at_end(s, 500, 4000)[:3], (3042, 2958, 'ok'))
        # no core window that close
        self.assertEqual(array_at_end('TTAAAG' * 600 + 'CCCTAA' * 500 + self.body, 500, 3000)[:3], (0, 0, 'none'))

    def test_min_array(self):
        # 360bp of repeat: its last window is 60% telomeric, so its walk reaches 400 (min_array), while the array
        # ends with its last hexamer (the body has one by chance 13bp on, too far from the array to count)
        s = 'CCCTAA' * 60 + self.body
        self.assertEqual(array_at_end(s, 500, 3000)[:3], (0, 360, 'short'))
        self.assertEqual(array_at_end(s, 400, 3000)[:3], (0, 360, 'ok'))

    def test_inner_edge(self):
        """ the array ends with the last telomeric hexamer of its last window, not where the window does, in the
        subtelomere (HG02132#2's 6p array ends in CCCTCA, 11bp after the hexamer before it, then subtelomere, which
        used to be anchored too); a hexamer more than 12bp after the one before it is on its own, and doesn't count """
        sub = 'GCCGGCCCGCCCGCCCGGGTCTGACCTGAGGAGAACTGTGCTCC'
        s = 'CCCTAA' * 495 + 'CCCTCGCGGTACCCTCA' + sub + self.body
        self.assertEqual(array_at_end(s, 500, 3000)[:3], (0, 2987, 'ok'))
        s = 'CCCTAA' * 495 + 'CCCTCGCGGTAGGCCCTCA' + sub + self.body
        self.assertEqual(array_at_end(s, 500, 3000)[:3], (0, 2970, 'ok'))

    def test_open_end(self):
        """ an array whose walk runs out of the sequence may go on past it, whatever phase the sequence ends in (its
        inner edge, the end of its last hexamer, is only at the end of the sequence in one phase of six) """
        for k in range(6):
            s = ('CCCTAA' * 1000)[k:5000 + k]
            self.assertEqual(array_at_end(s, 500, 3000), (0, 5000 - (5000 + k) % 6, 'ok', True), k)
            # two windows of body after the array: it ends there
            self.assertFalse(array_at_end(s[:4800] + self.body, 500, 3000)[3], k)
        self.assertFalse(array_at_end(self.body, 500, 3000)[3])

    def test_perfect_run(self):
        self.assertEqual(perfect_run('CCCTAA' * 20 + 'CCCTGA' + 'CCCTAA' * 30 + 'A'), 180)
        self.assertEqual(perfect_run('CCCTAAA' * 100), 6)
        self.assertLess(perfect_run(self.body), 18)
        self.assertLess(perfect_run(degenerate_repeat(random.Random(1), 3000)), 72)

    def test_no_array(self):
        self.assertEqual(array_at_end(self.body, 500, 3000)[:3], (0, 0, 'none'))
        self.assertEqual(array_at_end('N' * 10000 + self.body, 500, 3000)[:3], (0, 0, 'none'))
        self.assertEqual(array_at_end('', 500, 3000)[:3], (0, 0, 'none'))
        self.assertEqual(array_at_end('CCCTAA', 500, 3000)[2], 'short')


class TestTrimPaf(unittest.TestCase):

    def test_matches_per_base_trim(self):
        """ the op-level trim equals the prototype's per-base one, on random cigars, both strands, any interval """
        rng = random.Random(7)
        for _ in range(3000):
            ops = [(rng.randint(1, 6), rng.choice('==XID')) for _ in range(rng.randint(1, 8))]
            if not any(op in '=X' for _, op in ops):
                ops.append((1, '='))
            qlen = sum(n for n, op in ops if op in '=XI')
            tlen = sum(n for n, op in ops if op in '=XD')
            qs, ts = rng.randint(0, 5), rng.randint(0, 5)
            t = ['q', str(qs + qlen + 5), str(qs), str(qs + qlen), rng.choice('+-'), 't', str(ts + tlen + 5), str(ts),
                 str(ts + tlen), '0', '0', '60', 'cg:Z:' + merge_ops(ops), 'tp:A:P']
            t[9] = str(sum(n for n, op in ops if op == '='))
            t[10] = str(sum(n for n, _ in ops))
            lo = rng.randint(qs - 2, qs + qlen + 1)
            hi = rng.randint(lo, qs + qlen + 2)
            expected = reference_trim(t, lo, hi)
            got = trim_paf(t, lo, hi)
            self.assertEqual(got, expected, (t, lo, hi))
            if expected is t:
                self.assertIs(got, t)

    def test_cigar_must_fit(self):
        t = ['q', '100', '10', '20', '+', 't', '100', '0', '10', '10', '10', '60', 'cg:Z:12=']
        with self.assertRaises(ValueError):
            trim_paf(t, 15, 100)
        with self.assertRaises(ValueError):
            trim_paf(t[:12], 15, 100)


class TestAnchorAlignment(unittest.TestCase):

    def test_reference_phase(self):
        """ an array in phase 2 for 22bp and then 2bp short, the rest in phase 4 (CHM13's chr5p array is like that):
        aligned as the others, at the offset most of its outer 600bp match, it starts with an insertion; as the
        reference's, at the offset its outermost bases are in, a repeat unit further in """
        unit = 'CCCTAA' * 200
        array = unit[2:24] + unit[26:1000]
        r, start, ops, used = anchor_alignment(array, 1200)
        self.assertEqual((r, start, ops[0]), (4, 4, (2, 'I')))
        r, start, ops, used = anchor_alignment(array, 1200, reference=True)
        self.assertEqual((r, start, ops[0][1]), (2, 8, '='))
        self.assertGreaterEqual(ops[0][0], 22)
        self.assertEqual(sum(n for n, op in ops if op in 'XI'), 0)

    def test_reference_outer_mismatch(self):
        """ a reference array whose outermost base is not telomeric starts with that base as a mismatch, not an
        insertion, in the phase of the rest, one repeat unit in """
        self.assertEqual(anchor_alignment('G' + 'CCCTAA' * 100, 1200, reference=True), (5, 11, [(1, 'X'), (600, '=')],
                                                                                          601))
        self.assertEqual(anchor_alignment('CCCTAA' * 100, 1200, reference=True), (0, 6, [(600, '=')], 600))

    def test_pin_outer_base(self):
        target = 'CCCTAA' * 10
        # a leading insertion: the array's first bases are paired with the target bases just outside
        self.assertEqual(pin_outer_base('GCCCTAA', target, 6, [(1, 'I'), (6, '=')]), (5, [(1, 'X'), (6, '=')]))
        self.assertEqual(pin_outer_base('ACCCTAA', target, 6, [(1, 'I'), (6, '=')]), (5, [(1, '='), (6, '=')]))
        # a leading deletion moves the start in
        self.assertEqual(pin_outer_base('CTAA', target, 6, [(2, 'D'), (4, '=')]), (8, [(4, '=')]))
        self.assertEqual(pin_outer_base('GCCTAA', target, 6, [(1, 'D'), (1, 'I'), (5, '=')]),
                         (6, [(1, 'X'), (5, '=')]))
        # nothing to do, nothing aligned, or no target outside
        self.assertEqual(pin_outer_base('CCCTAA', target, 6, [(6, '=')]), (6, [(6, '=')]))
        self.assertEqual(pin_outer_base('CC', target, 6, [(2, 'I')]), (6, [(2, 'I')]))
        self.assertEqual(pin_outer_base('GGCCCTAA', target, 1, [(2, 'I'), (6, '=')]), (1, [(2, 'I'), (6, '=')]))


class TestFastaEnds(unittest.TestCase):

    def setUp(self):
        self.work_dir = tempfile.mkdtemp(prefix='telomere_anchorsTest_')

    def tearDown(self):
        shutil.rmtree(self.work_dir)

    def test_ends(self):
        """ wrapped and unwrapped lines, a '>' in a description, an empty record, CRLF, no final line break, text
        before the first header; headers and lines cut by every block boundary """
        rng = random.Random(5)
        a, c, d = random_dna(rng, 1000), random_dna(rng, 7), random_dna(rng, 333)
        text = ('text before any header\n>a desc>x\n' + '\n'.join(a[i:i + 77] for i in range(0, 1000, 77)) +
                '\n>b\n>c\n' + c + '\r\n>d\n' + d[:100] + '\n' + d[100:] + '\n>e\nACGT')
        records = [('a', a), ('b', ''), ('c', c), ('d', d), ('e', 'ACGT')]
        for gz in (False, True):
            path = os.path.join(self.work_dir, 'x.fa')
            with (gzip.open if gz else open)(path, 'wb') as f:
                f.write(text.encode())
            for block_size in (1, 2, 5, 64, 1 << 20):
                self.assertEqual(list(fasta_ends(path, 50, block_size=block_size)),
                                 [(n, len(s), s[:50].encode(), s[-50:].encode()) for n, s in records], (gz, block_size))

    def test_header_at_end_without_line_break(self):
        path = os.path.join(self.work_dir, 'x.fa')
        with open(path, 'w') as f:
            f.write('>a\nACGT\n>b')
        self.assertEqual(list(fasta_ends(path, 10, block_size=3)), [('a', 4, b'ACGT', b'ACGT'), ('b', 0, b'', b'')])


class TestAnchorTelomeres(AnchorTestCase):

    def test_plain_arms(self):
        contigs = self.plain_contigs()
        anchors = self.run_anchors(contigs)[0]
        # every haplotype end, none of the reference's (it has no arrays)
        self.assertEqual(sorted((a['event'], a['end'], a['arm'], a['strand']) for a in anchors),
                         sorted((e, end, arm, '+') for e in ['H{}#1'.format(i) for i in range(6)]
                                for end, arm in (('5p', 'p'), ('3p', 'q'))))
        # T = longest anchored array + 1000, rounded up to a multiple of 6
        self.assertEqual(set(a['target_length'] for a in anchors if a['arm'] == 'p'), {7002})
        self.assertEqual(set(a['target_length'] for a in anchors if a['arm'] == 'q'), {5202})
        for a in anchors:
            self.assertEqual(a['outer_bp'], 0)
            if a['event'] != 'H4#1':   # canonical arrays align exactly, from the target's outer end
                self.assertEqual((a['identity'], a['offset']), (1.0, 0))
                self.assertEqual(a['t_end'] - a['t_start'], a['array_bp'])
        self.assertIn('12 telomere anchors (6 p to 7002bp, 6 q to 5202bp)', anchor_summary(anchors))

    def test_reference_anchor(self):
        """ the reference's arrays are anchored a repeat unit further in, in the phase of their outermost bases (its p
        array starts like CHM13's chr5p, see TestAnchorAlignment), and mirror when its contig is reverse complemented;
        the haplotypes' anchors are the same as without a reference """
        unit = 'CCCTAA' * 600
        ref = 'id=REF|' + CHROM
        haps = self.plain_contigs()[:6]

        def run(tag, reverse=False, **kwargs):
            contigs = haps + [('REF', ref) + self.chrom.contig(ref, unit[2:24] + unit[26:2426], 'CCCTAA' * 300,
                                                               reverse=reverse)]
            anchors, _, _, out_paf, _ = self.run_anchors(contigs, tag=tag, **kwargs)
            return anchors, self.paf_fields(out_paf)
        anchors, lines = run('r', reference='REF')
        ref_anchors = {a['arm']: a for a in anchors if a['event'] == 'REF'}
        self.assertEqual((ref_anchors['p']['offset'], ref_anchors['p']['t_start']), (2, 8))
        self.assertEqual((ref_anchors['q']['offset'], ref_anchors['q']['target_length'] - ref_anchors['q']['t_end']),
                         (0, 6))
        self.assertEqual(ref_anchors['q']['identity'], 1.0)
        # as a haplotype, its p array would start at offset 4 with an insertion
        plain_anchors, plain_lines = run('h')
        self.assertEqual([(a['offset'], a['t_start']) for a in plain_anchors
                          if a['event'] == 'REF' and a['arm'] == 'p'], [(4, 4)])
        is_hap_anchor = lambda t: '|tel_' in t[5] and not t[0].startswith('id=REF|')
        self.assertEqual([t for t in lines if is_hap_anchor(t)], [t for t in plain_lines if is_hap_anchor(t)])
        # the reverse complemented reference's anchor lines mirror the forward ones
        _, mirrored = run('m', reverse=True, reference='REF')
        is_ref_anchor = lambda t: '|tel_' in t[5] and t[0] == ref
        self.assertEqual(sorted(t for t in mirrored if is_ref_anchor(t)),
                         sorted(mirror_line(t) for t in lines if is_ref_anchor(t)))

    def test_trimmed_lines(self):
        """ the s1 and s3 lines reach into the arrays and are cut back to their inner edges; the line within the p
        array is dropped; the reference's lines are untouched """
        _, _, _, out_paf, _ = self.run_anchors(self.plain_contigs())
        lines = [t for t in self.paf_fields(out_paf) if '|tel_' not in t[5]]
        self.assertEqual(len(lines), 6 * 3 + 3)
        for t in lines:
            if t[0].startswith('id=H'):
                # the arrays are 600bp into s1 and the last 600 of s3: trimmed lines start and end on the body (H4's
                # variant array too, though it is not a whole number of 100bp windows, as it ends with a hexamer)
                if t[5].endswith('|s1'):
                    self.assertEqual(int(t[7]), 600)
                if t[5].endswith('|s3'):
                    self.assertEqual(int(t[8]), int(t[6]) - 600)

    def test_mirrored(self):
        """ every haplotype contig reverse complemented: every line of theirs is the exact mirror of the forward
        run's (the reference's are the same) """
        base = self.paf_fields(self.run_anchors(self.plain_contigs(), tag='f')[3])
        anchors, _, _, out_paf, _ = self.run_anchors(self.plain_contigs(lambda i: True), tag='r')
        conv = lambda t: mirror_line(t) if t[0].startswith('id=H') else t
        mirrored = self.paf_fields(out_paf)
        is_anchor = lambda t: '|tel_' in t[5]
        self.assertEqual([t for t in mirrored if not is_anchor(t)], [conv(t) for t in base if not is_anchor(t)])
        self.assertEqual(sorted(t for t in mirrored if is_anchor(t)), sorted(conv(t) for t in base if is_anchor(t)))
        self.assertEqual(set((a['end'], a['orientation'], a['arm'], a['strand']) for a in anchors),
                         {('5p', '-', 'q', '-'), ('3p', '-', 'p', '-')})

    def test_mixed_orientations(self):
        base = self.paf_fields(self.run_anchors(self.plain_contigs(), tag='f')[3])
        anchors, _, _, out_paf, _ = self.run_anchors(self.plain_contigs(lambda i: i % 2 == 1), tag='m')
        flipped = set('id=H{}#1|{}'.format(i, CHROM) for i in (1, 3, 5))
        conv = lambda t: mirror_line(t) if t[0] in flipped else t
        mixed = self.paf_fields(out_paf)
        is_anchor = lambda t: '|tel_' in t[5]
        self.assertEqual([t for t in mixed if not is_anchor(t)], [conv(t) for t in base if not is_anchor(t)])
        self.assertEqual(sorted(t for t in mixed if is_anchor(t)), sorted(conv(t) for t in base if is_anchor(t)))
        self.assertEqual(len(anchors), 12)

    def test_orientation_by_majority(self):
        """ the strand holding most aligned bases orients the contig: an inverted piece doesn't flip it """
        name = 'id=H0#1|' + CHROM
        seq, lines = self.chrom.contig(name, 'CCCTAA' * 300, '', reverse=True)
        # a 300bp piece of the reverse complemented contig aligned forward, to a segment of its own
        x = len(seq) // 2
        piece = ('id={}|s4'.format(MG), seq[x:x + 300])
        extra = paf_line(name, seq, x, x + 300, piece[0], piece[1], 0, 300)
        anchors = self.run_anchors([('H0#1', name, seq, lines + [extra])], mg_records=self.chrom.segments + [piece])[0]
        self.assertEqual([(a['end'], a['orientation'], a['arm']) for a in anchors], [('3p', '-', 'p')])

    def test_outer_junk(self):
        """ the prototype's junk cases: junk (CTTTAA, TTAAAG, TTTAAC repeats) outside the array is left out, degenerate
        telomere kept in, an edge within 11bp moved to the end, and an array past max_outer not anchored; each contig
        also reverse complemented, whose anchors must mirror the forward ones """
        contigs, expect = [], {}
        for i, p5 in enumerate(JUNK_CASES):
            p3 = JUNK_CASES[(i + 7) % len(JUNK_CASES)]
            o5, o3 = junk_sequence(p5, i), junk_sequence(p3, i + 1)
            for reverse in (False, True):
                event = 'J{}{}#1'.format(i, 'm' if reverse else '')
                name = 'id={}|{}'.format(event, CHROM)
                seq, lines = self.chrom.contig(name, 'CCCTAA' * 500, 'CCCTAA' * 500, o5, o3, reverse=reverse)
                contigs.append((event, name, seq, lines))
                e5 = expected_edge(p5, (o5 + 'CCCTAA' * 500 + self.chrom.body[:2000]).upper())
                e3 = expected_edge(p3, (o3 + 'CCCTAA' * 500 + rc(self.chrom.body[-2000:])).upper())
                expect[(event, '5p')], expect[(event, '3p')] = (e3, e5) if reverse else (e5, e3)
        anchors, _, _, out_paf, _ = self.run_anchors(contigs)
        got = {(a['event'], a['end']): a['outer_bp'] for a in anchors}
        self.assertEqual({k: got.get(k) for k in expect}, expect)
        # junk of 3003bp or more: 4 cases, each at both ends of a contig and its reverse complement
        self.assertEqual((len(anchors), sum(1 for e in expect.values() if e is None)), (88, 16))
        # the reverse complemented copy's anchor lines mirror the forward contig's
        lines = [t for t in self.paf_fields(out_paf) if '|tel_' in t[5]]
        fwd = [t for t in lines if not t[0].split('|')[0].endswith('m#1')]
        mir = set(tuple(t) for t in lines if t[0].split('|')[0].endswith('m#1'))
        self.assertEqual(set(tuple(mirror_line([t[0].replace('#1|', 'm#1|')] + t[1:])) for t in fwd), mir)

    def test_min_array(self):
        name = 'id=H0#1|' + CHROM
        seq, lines = self.chrom.contig(name, 'CCCTAA' * 60, 'CCCTAA' * 300)
        contigs = [('H0#1', name, seq, lines)]
        self.assertEqual([a['arm'] for a in self.run_anchors(contigs, tag='a')[0]], ['q'])
        self.assertEqual(sorted(a['arm'] for a in self.run_anchors(contigs, tag='b', min_array=400)[0]), ['p', 'q'])

    def test_max_outer(self):
        name = 'id=H0#1|' + CHROM
        seq, lines = self.chrom.contig(name, 'CCCTAA' * 500, 'CCCTAA' * 500, outer5='TTAAAG' * 600)
        contigs = [('H0#1', name, seq, lines)]
        self.assertEqual([a['arm'] for a in self.run_anchors(contigs, tag='a')[0]], ['q'])
        anchors = self.run_anchors(contigs, tag='b', max_outer=4000)[0]
        self.assertEqual(sorted((a['arm'], a['outer_bp']) for a in anchors), [('p', 3600), ('q', 0)])

    def test_scan_limit(self):
        """ an array longer than the bp read at a contig end (H2's 6000bp p array, with 5000bp read) is anchored as far
        as it was read, and logged; H5's 4800bp p array, followed by two windows of body within what was read, isn't """
        with mock.patch.object(telomere_anchors, 'SCAN', 2000):
            with self.assertLogs(telomere_anchors.logger, level='INFO') as log:
                anchors = self.run_anchors(self.plain_contigs())[0]
        self.assertTrue(any('1 arrays reach the 5000bp scan limit, such as id=H2#1|' + CHROM in l for l in log.output),
                        log.output)
        self.assertEqual([a['array_bp'] for a in anchors if a['event'] == 'H2#1' and a['arm'] == 'p'], [4998])

    def test_gzip(self):
        """ gzipped fastas (named .fa, as in cactus_cons) give the same outputs """
        _, _, _, paf_a, mg_a = self.run_anchors(self.plain_contigs(), tag='a')
        _, fasta_paths, _, paf_b, mg_b = self.run_anchors(self.plain_contigs(), gz=True, tag='b')
        with open(fasta_paths[MG], 'rb') as f:
            self.assertEqual(f.read(2), b'\x1f\x8b')
        self.assertSameFile(paf_a, paf_b)
        self.assertSameFile(mg_a, mg_b)

    def test_lowercase(self):
        """ soft-masked sequence gives the same outputs """
        _, _, _, paf_a, mg_a = self.run_anchors(self.plain_contigs(), tag='a')
        lower = [(e, n, s.lower(), l) for e, n, s, l in self.plain_contigs()]
        _, _, _, paf_b, mg_b = self.run_anchors(lower, tag='b')
        self.assertSameFile(paf_a, paf_b)
        self.assertSameFile(mg_a, mg_b)

    def test_threads(self):
        """ the process pool gives the same outputs as one thread """
        _, _, _, paf_a, mg_a = self.run_anchors(self.plain_contigs(), tag='a')
        with mock.patch.object(telomere_anchors, 'WORKER_BYTES', 1):
            _, _, _, paf_b, mg_b = self.run_anchors(self.plain_contigs(), tag='b', threads=3)
        self.assertSameFile(paf_a, paf_b)
        self.assertSameFile(mg_a, mg_b)

    def test_pool_size(self):
        """ a process per GB of fasta, at most 16 and at most threads (which can be a float, as toil's job.cores),
        and none for a small input, such as a small chromosome with many genomes """
        with mock.patch.object(telomere_anchors, 'ProcessPoolExecutor') as pool:
            self.assertIsInstance(telomere_anchors.process_pool(64, 466, 7 << 20), telomere_anchors.InlineExecutor)
            telomere_anchors.process_pool(64, 466, 116 << 30)
            telomere_anchors.process_pool(64.0, 466, 3 << 30)
            telomere_anchors.process_pool(4, 2, 10 << 30)
        self.assertEqual([call.kwargs['max_workers'] for call in pool.call_args_list], [16, 3, 2])

    def test_interstitial_repeat(self):
        """ a fragment of an assembly that starts at an interstitial telomeric repeat (ITS) mid-chromosome, as one that
        breaks next to it would: an ITS is degenerate, so it is not anchored (anchoring it to the telomere would make a
        false deletion from the telomere to the middle of the chromosome), while the fragment's real telomere is.  A
        plant's (CCCTAAA)n is not anchored either """
        rng = random.Random(4)
        body = self.chrom.body
        (n1, s1), (n2, s2), (n3, s3) = self.chrom.segments
        its = degenerate_repeat(rng, 900)
        self.assertEqual(array_at_end(its + body, 500, 3000)[2], 'ok')
        name = 'id=F#1|frag'
        frag = its + body[6000:] + rc('CCCTAA' * 300)
        lines = [paf_line(name, frag, 900, 2900, n2, s2, 2000, 4000),
                 paf_line(name, frag, 2900, 6900, n3, s3, 0, 4000)]
        plant = 'id=P#1|' + CHROM
        seq, plant_lines = self.chrom.contig(plant, 'CCCTAAA' * 300, 'CCCTAAA' * 300)
        contigs = self.plain_contigs()[:1] + [('F#1', name, frag, lines), ('P#1', plant, seq, plant_lines)]
        with self.assertLogs(telomere_anchors.logger, level='INFO') as log:
            anchors = self.run_anchors(contigs)[0]
        self.assertEqual(sorted((a['contig'], a['end']) for a in anchors),
                         sorted([(contigs[0][1], '5p'), (contigs[0][1], '3p'), (name, '3p')]))
        self.assertIn('degenerate_array 3', log.output[-1])
        self.assertTrue(any('3 arrays without 72bp of perfect repeat' in l for l in log.output), log.output)
        # with no floor, as before, the ITS would be anchored to the p telomere
        anchors = self.run_anchors(contigs, tag='b', min_run=0)[0]
        self.assertIn((name, '5p', 'p'), [(a['contig'], a['end'], a['arm']) for a in anchors])

    def test_reference_contigs(self):
        """ the rc:Z: tags of a contig's lines must put most of it on the chromosome: a job of several reference
        contigs (the --noSplit job, named all, or the chrOther bin) anchors nothing rather than give every telomere
        the same target per arm, and in a chromosome's job a contig mostly on another reference contig, or whose lines
        have no rc:Z: tags, is not anchored """
        contigs = self.plain_contigs()
        fasta_paths, paf_path = write_subproblem(self.work_dir, contigs, self.chrom.segments)
        for chrom in ('all', 'chrOther'):
            out_paf, out_mg = os.path.join(self.work_dir, chrom + '.paf'), os.path.join(self.work_dir, chrom + '.fa')
            with self.assertLogs(telomere_anchors.logger, level='INFO') as log:
                anchors = anchor_telomeres(chrom, fasta_paths, paf_path, out_paf, out_mg, minigraph_event=MG)
            self.assertEqual(anchors, [])
            self.assertSameFile(paf_path, out_paf)
            self.assertSameFile(fasta_paths[MG], out_mg)
            self.assertIn('other_chromosome 12', '\n'.join(log.output))
            self.assertIn('no contig with an array is mostly on reference contig ' + chrom, '\n'.join(log.output))

        def with_rc(lines, ref_contigs):
            return [[x for x in t if not x.startswith('rc:Z:')] + (['rc:Z:' + r] if r else [])
                    for t, r in zip(lines, ref_contigs)]
        # H0: 2/3 on chrU; H1: no tags; H2: 1/3 on chrU, still anchored
        contigs[0] = contigs[0][:3] + (with_rc(contigs[0][3], ['chrU', 'chrU', CHROM]),)
        contigs[1] = contigs[1][:3] + (with_rc(contigs[1][3], [None, None, None]),)
        contigs[2] = contigs[2][:3] + (with_rc(contigs[2][3], [CHROM, 'chrU', CHROM, CHROM]),)
        anchors = self.run_anchors(contigs, tag='b')[0]
        self.assertEqual(sorted(set(a['event'] for a in anchors)), ['H{}#1'.format(i) for i in range(2, 6)])

    def test_fuzz(self):
        """ random subproblems of short, noisy, junk-capped and soft-masked contigs, with random lines on either
        strand, a few on another reference contig: every output verifies """
        rng = random.Random(3)
        n_anchors = 0
        for case in range(80):
            contigs, segments = [], []
            for c in range(rng.randint(1, 4)):
                name = 'id=G{}#1|ctg{}'.format(c, case)
                junk = [junk_sequence(rng.choice(JUNK_CASES), rng.randint(0, 5)) if rng.random() < 0.3 else ''
                        for _ in range(2)]
                arrays = [noisy_telomere(rng, rng.choice([0, 300, 600, 1500, 5000]), rng.choice([0, 0.03, 0.1, 0.3]))
                          if rng.random() < 0.7 else '' for _ in range(2)]
                seq = junk[0] + arrays[0] + random_dna(rng, rng.choice([0, 50, 500, 3000, 20000])) + rc(arrays[1]) + \
                    rc(junk[1])
                if not seq:
                    continue
                if rng.random() < 0.3:
                    seq = seq.lower()
                lines = []
                for k in range(rng.randint(0, 5)):
                    q0 = rng.randrange(len(seq))
                    q1 = rng.randint(q0 + 1, min(len(seq), q0 + 6000))
                    strand = '+' if rng.random() < 0.7 else '-'
                    segment = list(seq[q0:q1].upper() if strand == '+' else rc(seq[q0:q1].upper()))
                    for _ in range(len(segment) // 200):   # a few edits, so cigars have X, I and D
                        i = rng.randrange(len(segment))
                        segment[i] = rng.choice('ACGT') if rng.random() < 0.6 else ''
                    segments.append(('id={}|s{}_{}'.format(MG, c, k), ''.join(segment) or 'A'))
                    lines.append(paf_line(name, seq, q0, q1, segments[-1][0], segments[-1][1], 0,
                                          len(segments[-1][1]), strand, CHROM if rng.random() < 0.9 else 'chrU'))
                contigs.append(('G{}#1'.format(c), name, seq, lines))
            segments = segments or [('id={}|s0'.format(MG), 'ACGT')]
            anchors = self.run_anchors(contigs, tag=str(case), mg_records=segments)[0]
            n_anchors += len(anchors)
        self.assertGreater(n_anchors, 20)

    def test_empty_paf(self):
        contigs = [(e, n, s, []) for e, n, s, _ in self.plain_contigs()]
        anchors, fasta_paths, paf_path, out_paf, out_mg = self.run_anchors(contigs)
        self.assertEqual(anchors, [])
        self.assertSameFile(paf_path, out_paf)
        self.assertSameFile(fasta_paths[MG], out_mg)

    def test_nothing_to_anchor(self):
        contigs = []
        for i in range(3):
            name = 'id=H{}#1|{}'.format(i, CHROM)
            contigs.append(('H{}#1'.format(i), name) + self.chrom.contig(name))
        for gz in (False, True):
            anchors, fasta_paths, paf_path, out_paf, out_mg = self.run_anchors(contigs, gz=gz, tag=str(gz))
            self.assertEqual(anchors, [])
            self.assertSameFile(paf_path, out_paf)
            self.assertSameFile(fasta_paths[MG], out_mg, gunzip_a=True)

    def test_production_like_inputs(self):
        """ a reference with N-run ends, and a fragmented assembly: tiny contigs, contigs without PAF lines, a short
        contig with an array, one made of telomere (two arrays that meet are not anchored) """
        contigs = list(self.plain_contigs()[:2])
        ref = 'id=GRCh38|' + CHROM
        seq, lines = self.chrom.contig(ref, outer5='N' * 10000, outer3='N' * 10000)
        contigs.append(('GRCh38', ref, seq, lines))
        rng = random.Random(9)
        n2 = self.chrom.segments[1]
        frag = lambda i: 'id=F#1|ctg{}'.format(i)
        contigs.append(('F#1', frag(0), random_dna(rng, 50), []))
        contigs.append(('F#1', frag(1), 'CCCTAA' * 50, []))
        contigs.append(('F#1', frag(2), 'CCCTAA' * 500 + random_dna(rng, 2000), []))   # can't be oriented
        body = n2[1][1000:3000]
        seq = 'CCCTAA' * 200 + body   # shorter than the 103kb scan: both ends read the whole contig
        contigs.append(('F#1', frag(3), seq, [paf_line(frag(3), seq, 1200, len(seq), n2[0], n2[1], 1000, 3000)]))
        # 100bp windows alternating between the strands: each end's walk spans the whole contig
        seq = ('CCCTAA' * 16 + 'CCCT' + 'TTAGGG' * 16 + 'TTAG') * 10
        contigs.append(('F#1', frag(4), seq, [paf_line(frag(4), seq, 0, 100, n2[0], n2[1], 0, 100)]))
        for i in range(5, 40):
            seq = random_dna(rng, rng.randint(200, 3000))
            contigs.append(('F#1', frag(i), seq, []))
        with self.assertLogs(telomere_anchors.logger, level='INFO') as log:
            anchors = self.run_anchors(contigs)[0]
        self.assertEqual(sorted((a['contig'], a['end']) for a in anchors),
                         sorted([('id=H0#1|' + CHROM, '5p'), ('id=H0#1|' + CHROM, '3p'), ('id=H1#1|' + CHROM, '5p'),
                                 ('id=H1#1|' + CHROM, '3p'), (frag(3), '5p')]))
        summary = [l for l in log.output if 'anchored 3 p ends and 2 q ends' in l]
        self.assertEqual(len(summary), 1, log.output)
        self.assertIn('no_paf_lines 1', summary[0])
        self.assertIn('arrays_overlap 2', summary[0])
        self.assertIn('short_array', summary[0])

    def test_odd_input(self):
        """ skipped and counted, never a crash: unparsable and blank PAF lines, a PAF length that disagrees with the
        fasta, PAF queries not in the fastas, a name in two fastas, an unreadable fasta, a line whose cigar doesn't
        fit, and a last line without a line break """
        contigs = self.plain_contigs()[:5]
        bad_line, mismatch, dup, untrimmable = (c[1] for c in contigs[:4])
        contigs[0] = contigs[0][:3] + (contigs[0][3] + [[bad_line, '5', '1']],)
        contigs[1] = contigs[1][:3] + ([[t[0], str(int(t[1]) + 1)] + t[2:] for t in contigs[1][3]],)
        contigs.append(('D#1', dup, 'CCCTAA' * 100, []))
        t = list(contigs[3][3][0])
        t[12] = 'cg:Z:{}='.format(int(t[3]) - int(t[2]) + 1)
        contigs[3] = contigs[3][:3] + ([t] + contigs[3][3][1:],)
        contigs.append(('U#1', 'id=U#1|' + CHROM) + self.chrom.contig('id=U#1|' + CHROM, 'CCCTAA' * 300))
        fasta_paths, paf_path = write_subproblem(self.work_dir, contigs, self.chrom.segments)
        os.remove(fasta_paths['U#1'])
        with open(paf_path) as f:
            text = f.read()
        with open(paf_path, 'w') as f:
            f.write('\n' + text +
                    'id=X#1|missing\t100\t0\t50\t+\tid={}|s2\t4000\t0\t50\t50\t50\t60\tcg:Z:50='.format(MG))
        out_paf, out_mg = os.path.join(self.work_dir, 'out.paf'), os.path.join(self.work_dir, 'out.mg.fa')
        with self.assertLogs(telomere_anchors.logger, level='INFO') as log:
            anchors = anchor_telomeres(CHROM, fasta_paths, paf_path, out_paf, out_mg, minigraph_event=MG)
        self.assertEqual(sorted(set(a['contig'] for a in anchors)), [untrimmable, contigs[4][1]])
        summary = log.output[-1]
        for reason in ('bad_paf_line 2', 'length_mismatch 2', 'duplicate_name 3'):
            self.assertIn(reason, summary)
        self.assertTrue(any('could not be trimmed' in l for l in log.output))
        self.assertTrue(any('not in the fastas' in l for l in log.output))
        self.assertTrue(any('U.1.fa' in l for l in log.output))
        out = open(out_paf).read().split('\n')
        self.assertEqual(out[0], '')
        self.assertIn(text.split('\n')[0], out)   # the bad contig's first line, unchanged
        self.assertIn('\t'.join(t), out)            # the untrimmable line, unchanged
        self.assertTrue(out[-2].endswith('ct:i:1') and out[-1] == '')
        self.assertIn('id=X#1|missing\t100', '\n'.join(out))

    def test_already_anchored(self):
        """ inputs that already hold the targets are left as they are """
        _, fasta_paths, _, out_paf, out_mg = self.run_anchors(self.plain_contigs())
        fasta_paths[MG] = out_mg
        again_paf, again_mg = os.path.join(self.work_dir, 'again.paf'), os.path.join(self.work_dir, 'again.fa')
        self.assertEqual(anchor_telomeres(CHROM, fasta_paths, out_paf, again_paf, again_mg, minigraph_event=MG), [])
        self.assertSameFile(out_paf, again_paf)
        self.assertSameFile(out_mg, again_mg)


class TestVerifier(AnchorTestCase):
    """ the verifier catches what it is there to catch (the prototype's negative controls) """

    def junk_contigs(self):
        contigs = []
        for i, outer in enumerate(['TTAAAG' * 50, 'CTTTAA' * 505, '']):
            name = 'id=J{}#1|{}'.format(i, CHROM)
            contigs.append(('J{}#1'.format(i), name) + self.chrom.contig(name, 'CCCTAA' * 500, 'CCCTAA' * 400, outer))
        return contigs

    def test_edge_ignoring_junk(self):
        with mock.patch.object(telomere_anchors, 'outer_edge', lambda s, c: 0):
            args = self.run_anchors(self.junk_contigs(), verify=False)[1:]
        with self.assertRaisesRegex(AssertionError, 'not a telomere/junk boundary'):
            verify_outputs(*args)

    def test_inner_edge_in_subtelomere(self):
        # H4's variant array ends 37bp before its last window does
        with mock.patch.object(telomere_anchors, 'inner_edge', lambda s, x: x):
            args = self.run_anchors(self.plain_contigs(), verify=False)[1:]
        with self.assertRaisesRegex(AssertionError, 'is not the end of a telomeric hexamer'):
            verify_outputs(*args)

    def test_reference_not_inset(self):
        """ a reference anchored like the haplotypes is caught """
        args = self.run_anchors(self.plain_contigs(), verify=False)[1:]
        verify_outputs(*args)
        with self.assertRaisesRegex(AssertionError, 'not pinned 6-11bp'):
            verify_outputs(*args, reference='H0#1')

    def test_past_max_outer(self):
        args = self.run_anchors(self.junk_contigs(), verify=False, max_outer=3100)[1:]
        verify_outputs(*args, max_outer=3100)
        with self.assertRaisesRegex(AssertionError, 'further inward than max_outer'):
            verify_outputs(*args, max_outer=3000)

    def test_wrong_cigar(self):
        fasta_paths, paf_path, out_paf, out_mg = self.run_anchors(self.plain_contigs())[1:]
        with open(out_paf) as f:
            lines = f.read().split('\n')
        i = next(i for i, l in enumerate(lines) if 'id=H0#1' in l and '|tel_{}p'.format(CHROM) in l)
        t = lines[i].split('\t')
        t[7], t[8] = str(int(t[7]) + 1), str(int(t[8]) + 1)   # off by one on the target: '=' runs mismatch
        lines[i] = '\t'.join(t)
        with open(out_paf, 'w') as f:
            f.write('\n'.join(lines))
        with self.assertRaisesRegex(AssertionError, 'is not a match|not pinned'):
            verify_outputs(fasta_paths, paf_path, out_paf, out_mg)

    def test_untrimmed_line(self):
        fasta_paths, paf_path, out_paf, out_mg = self.run_anchors(self.plain_contigs())[1:]
        with open(out_paf) as f:
            anchor_lines = [l for l in f if '|tel_' in l]
        with open(out_paf, 'w') as f:   # the input lines unchanged, then the anchors
            f.writelines(open(paf_path).readlines() + anchor_lines)
        with self.assertRaisesRegex(AssertionError, 'is not the expected'):
            verify_outputs(fasta_paths, paf_path, out_paf, out_mg)


if __name__ == '__main__':
    unittest.main()
