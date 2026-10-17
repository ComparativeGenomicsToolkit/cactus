#!/usr/bin/env python3

"""
Train a scoring model for abPOA with cPecan's pair HMM, from the genomes being aligned: the model the pangenome
pipeline (<graphmap lastTrain>) gives cactus-align, and minigraph's alignment penalties (lastTrainMap).

Random pieces of a few genomes are mapped to the reference with minimap2, and cPecan's five-state pair HMM (match;
a short and a long gap in each sequence) is fitted to those alignments by expectation maximisation:
cPecanRealign --outputExpectations gives the expected transitions and emissions in a band around each alignment,
and they are summed, symmetrised and normalised, from cPecan's default model, a few times over.  Pieces long
enough to hold long indels let EM fit the long-gap state too, which becomes abPOA's second gap piece.

This replaces LAST's last-train, which is GPL-licensed.  On 285 abPOA windows of human chr10 (STR sites and other
indels, GRCh38 and 12 haplotypes), the two halve the mismatches and gap runs of abPOA's untrained scores and are
within 1.5% of each other.

Copyright (C) 2026 by Glenn Hickey

Released under the MIT license, see LICENSE.txt
"""

import gzip
import os
import random
import re
from concurrent.futures import ThreadPoolExecutor

from cactus.shared.common import cactus_call

MATCH, SHORT_GAP_X, SHORT_GAP_Y, LONG_GAP_X, LONG_GAP_Y = range(5)   # cPecan's five states, in its order
# cPecan's built-in long-gap state (stateMachine5_construct), which training keeps when told not to fit it
# (train_long_gap): alignments too short to hold long indels turn it into a second short gap
LONG_GAP_OPEN = 0.001821479941473
LONG_GAP_EXTEND = 0.99656342579062
# the same states with the two sequences exchanged
XY_SWAP = [MATCH, SHORT_GAP_Y, SHORT_GAP_X, LONG_GAP_Y, LONG_GAP_X]
COMPLEMENT = bytes.maketrans(b'ACGTacgtNn', b'TGCAtgcaNn')
CIGAR_OPS = re.compile(r'(\d+)([MIDNX=])')


def open_fasta(path):
    """ a fasta file, gzipped (as an ancestor's is) or not, opened for reading bytes """
    with open(path, 'rb') as f:
        gzipped = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rb') if gzipped else open(path, 'rb')


def fasta_records(path, wanted=None):
    """ (name, sequence) of each record, or only of those named in wanted, the name being the header up to the
    first space (as minimap2, and so the PAF, have it) """
    name, chunks = None, []
    with open_fasta(path) as f:
        for line in f:
            if line.startswith(b'>'):
                if name is not None and (wanted is None or name in wanted):
                    yield name, b''.join(chunks)
                name, chunks = line[1:].split()[0].decode(), []
            elif wanted is None or name in wanted:
                chunks.append(line.rstrip())
    if name is not None and (wanted is None or name in wanted):
        yield name, b''.join(chunks)


def fasta_lengths(path):
    lengths, name, n = {}, None, 0
    with open_fasta(path) as f:
        for line in f:
            if line.startswith(b'>'):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split()[0].decode(), 0
            else:
                n += len(line.rstrip())
    if name is not None:
        lengths[name] = n
    return lengths


def sample_pieces(fasta_path, out_path, pieces, length, seed=1, max_n=0.05, prefix='trainPiece'):
    """ write up to pieces random, non-overlapping pieces of length bases of the genome at fasta_path, drawn
    uniformly over its sequences of at least that length, with no more than max_n of them N, as a fasta of
    sequences named <prefix><i>; and no more of them than would cover half the genome, so that a small genome
    is not mostly sampled twice over.  Returns how many were written. """
    lengths = fasta_lengths(fasta_path)
    usable = sorted(name for name, n in lengths.items() if n >= length)
    total = sum(lengths[name] for name in usable)
    pieces = min(pieces, total // (2 * length))
    if pieces <= 0:
        return 0
    rng = random.Random(seed)
    starts = {}
    for _ in range(pieces * 20):
        if sum(len(v) for v in starts.values()) >= pieces * 1.5:
            break
        x = rng.randrange(total)
        for name in usable:
            if x < lengths[name]:
                break
            x -= lengths[name]
        start = min(x, lengths[name] - length)
        if all(abs(start - s) >= length for s in starts.get(name, [])):
            starts.setdefault(name, []).append(start)
    written = 0
    with open(out_path, 'wb') as out:
        for name, seq in fasta_records(fasta_path, set(starts)):
            for start in sorted(starts[name]):
                piece = seq[start:start + length]
                if written < pieces and piece.count(b'N') + piece.count(b'n') <= max_n * length:
                    out.write(b'>%s%d\n%s\n' % (prefix.encode(), written, piece))
                    written += 1
    return written


def paf_to_cpecan_batches(paf_path, sample_path, genome_path, out_dir, batches):
    """ the alignments in paf_path between the sample (sample_path) and the genome (genome_path) as input for
    cPecanRealign, spread round-robin over batches pairs of files out_dir/<i>.fa and out_dir/<i>.cigar.  Each
    alignment's two spans are written as their own sequences, the genome's forward (X) and the sample's in the
    orientation that aligns to it (Y), with a cigar over them, so cPecanRealign sees every alignment on the
    forward strands.  Soft-masking is ignored, as it is by the rest of the pangenome pipeline (most inputs have
    none, and vg drops it), and abPOA aligns repeats like anything else.  Leaving out alignments mostly
    soft-masked, as last-train does, kept 323 of 1091 on a CHM13 chr10 reference (53% soft-masked), and with them
    went the long indels the long-gap state is fitted to.  Returns
    {'alignments', 'kept', 'aligned_bp', 'identity'}, the identity being the median over aligned bases, from the
    cigars' = and X (None for a mapper whose cigars have only M, which cPecanRealign takes all the same). """
    sample = dict(fasta_records(sample_path))
    records = []
    wanted = {}
    with open_fasta(paf_path) as f:
        for line in f:
            p = line.decode().rstrip('\n').split('\t')
            if p[0] in sample:
                q_name, q_start, q_end, strand, t_name, t_start, t_end = p[0], int(p[2]), int(p[3]), p[4], p[5], int(p[7]), int(p[8])
                sample_is_query = True
            elif p[5] in sample:
                q_name, q_start, q_end, strand, t_name, t_start, t_end = p[5], int(p[7]), int(p[8]), p[4], p[0], int(p[2]), int(p[3])
                sample_is_query = False
            else:
                continue
            cigar = next((x[5:] for x in p[12:] if x.startswith('cg:Z:')), None)
            if cigar is None:
                raise RuntimeError('PAF record without a cigar in {}: {}'.format(paf_path, line[:200]))
            records.append((q_name, q_start, q_end, strand, t_name, t_start, t_end, cigar, sample_is_query))
            wanted.setdefault(t_name, []).append(len(records) - 1)
    # the genome's spans, reading it once
    spans = {}
    for name, seq in fasta_records(genome_path, set(wanted)):
        for i in wanted[name]:
            spans[i] = seq[records[i][5]:records[i][6]]

    os.makedirs(out_dir, exist_ok=True)
    files = [(open(os.path.join(out_dir, '{}.fa'.format(b)), 'wb'), open(os.path.join(out_dir, '{}.cigar'.format(b)), 'w'))
             for b in range(batches)]
    stats = {'alignments': len(records), 'kept': 0, 'aligned_bp': 0}
    identities = []
    for i, (q_name, q_start, q_end, strand, t_name, t_start, t_end, cigar, sample_is_query) in enumerate(records):
        x, y = spans.get(i), sample[q_name][q_start:q_end]
        if x is None or len(x) != t_end - t_start:
            raise RuntimeError('{} has no span {}:{}-{} that {} aligns to'.format(genome_path, t_name, t_start, t_end, paf_path))
        if strand == '-':
            y = y.translate(COMPLEMENT)[::-1]
        # in the PAF's cigar I takes from its query and D from its target; cPecanRealign's sonLib cigar has D take
        # from X (the genome here) and I from Y (the sample)
        ops, aligned, matches, mismatches = [], 0, 0, 0
        for n, op in CIGAR_OPS.findall(cigar):
            n = int(n)
            if op in 'M=X':
                c = 'M'
                aligned += n
                matches += n if op == '=' else 0
                mismatches += n if op == 'X' else 0
            elif op == 'I':
                c = 'I' if sample_is_query else 'D'
            else:
                c = 'D' if sample_is_query else 'I'
            if ops and ops[-1][0] == c:
                ops[-1][1] += n
            else:
                ops.append([c, n])
        # the cigar runs along the PAF's target forward; when that is the sample, on the - strand the genome runs
        # backwards along it, so forward along the genome is the cigar reversed
        if strand == '-' and not sample_is_query:
            ops.reverse()
        stats['kept'] += 1
        stats['aligned_bp'] += aligned
        if matches + mismatches:
            identities.append((matches / (matches + mismatches), matches + mismatches))
        fa, cig = files[stats['kept'] % batches]
        fa.write(b'>x%d\n%s\n>y%d\n%s\n' % (i, x, i, y))
        # a sonLib cigar names its second sequence (Y) first
        cig.write('cigar: y{} 0 {} + x{} 0 {} + 0 {}\n'.format(i, len(y), i, len(x), ' '.join('{} {}'.format(c, n) for c, n in ops)))
    for fa, cig in files:
        fa.close()
        cig.close()
    identities.sort()
    half, seen = sum(weight for _, weight in identities) / 2, 0
    stats['identity'] = None
    for identity, weight in identities:
        seen += weight
        if seen >= half:
            stats['identity'] = identity
            break
    return stats


def read_expectations(path):
    """ (transitions, emissions, likelihood) from a cPecanRealign --outputExpectations file of a five-state model """
    with open(path) as f:
        t = [float(v) for v in f.readline().split()]
        e = [float(v) for v in f.readline().split()]
    if not t or int(t[0]) != 0 or len(t) != 27 or len(e) != 80:
        raise RuntimeError('{} is not the expectations of a five-state model'.format(path))
    return t[1:26], e, t[26]


def write_hmm(path, transitions, emissions):
    """ a five-state model as hmm_loadFromFile reads it, at full precision """
    with open(path, 'w') as f:
        f.write('0\t' + '\t'.join(repr(float(v)) for v in transitions) + '\t0.0\n')
        f.write('\t'.join(repr(float(v)) for v in emissions) + '\n')


def maximise(transitions, emissions, train_long_gap=False):
    """ the five-state model that the expected counts make, with the two sequences made alike (as pecan aligns
    them either way round) and the two strands alike: transitions[from * 5 + to], emissions[state * 16 + x * 4 + y],
    bases in the order A C G T.  The gap states share one base distribution, the letters they emit.  Unless
    train_long_gap, the long-gap state's transitions are held at cPecan's (LONG_GAP_OPEN, LONG_GAP_EXTEND), and the
    rest of the match state's share out what they leave in proportion to their counts: the most likely model with
    the long gaps fixed.  Train them only from alignments that have long gaps in them, as whole-haplotype ones do. """
    E = [[[max(emissions[s * 16 + x * 4 + y], 0.0) + 1e-6 for y in range(4)] for x in range(4)] for s in range(5)]
    T = [[max(transitions[i * 5 + j], 0.0) + 1e-6 for j in range(5)] for i in range(5)]
    T = [[(T[i][j] + T[XY_SWAP[i]][XY_SWAP[j]]) / 2 for j in range(5)] for i in range(5)]
    T = [[v / sum(row) for v in row] for row in T]
    if train_long_gap:
        return [v for row in T for v in row], gap_and_match_emissions(E)
    free = T[MATCH][MATCH] + T[MATCH][SHORT_GAP_X] + T[MATCH][SHORT_GAP_Y]
    for j in (MATCH, SHORT_GAP_X, SHORT_GAP_Y):
        T[MATCH][j] *= (1 - 2 * LONG_GAP_OPEN) / free
    T[MATCH][LONG_GAP_X] = T[MATCH][LONG_GAP_Y] = LONG_GAP_OPEN
    # as cPecan's own: extend or return to match, and a switch to the other sequence's long gap (which its recursion
    # does not use) of the open probability
    for g, other in ((LONG_GAP_X, LONG_GAP_Y), (LONG_GAP_Y, LONG_GAP_X)):
        T[g] = [0.0] * 5
        T[g][g], T[g][MATCH], T[g][other] = LONG_GAP_EXTEND, 1 - LONG_GAP_EXTEND, LONG_GAP_OPEN
    return [v for row in T for v in row], gap_and_match_emissions(E)


def gap_and_match_emissions(E):
    """ maximise's emissions: the match state's symmetrised over the two sequences and the two strands, and one base
    distribution for every gap state """
    # A C G T reversed is the complement of each: (x, y) and (3 - x, 3 - y) are the same pair on the other strand
    M = [[(E[MATCH][x][y] + E[MATCH][y][x] + E[MATCH][3 - x][3 - y] + E[MATCH][3 - y][3 - x]) / 4 for y in range(4)]
         for x in range(4)]
    total = sum(map(sum, M))
    M = [[v / total for v in row] for row in M]
    gap = [0.0] * 4
    for s in (SHORT_GAP_X, LONG_GAP_X):          # an X gap state emits x
        for x in range(4):
            gap[x] += sum(E[s][x])
    for s in (SHORT_GAP_Y, LONG_GAP_Y):          # a Y gap state emits y
        for y in range(4):
            gap[y] += sum(E[s][x][y] for x in range(4))
    gap = [(gap[b] + gap[3 - b]) / 2 for b in range(4)]
    gap = [v / sum(gap) for v in gap]
    emissions = []
    for s in range(5):
        for x in range(4):
            for y in range(4):
                emissions.append(M[x][y] if s == MATCH else gap[x] * gap[y])
    return emissions


def train_pair_hmm(batch_dir, batches, iterations, threads, train_long_gap=False):
    """ fit cPecan's five-state model to the alignments of paf_to_cpecan_batches by EM, iterations times from its
    default model, cPecanRealign running on threads batches at once, the long-gap state trained too with
    train_long_gap (see maximise).  Returns (transitions, emissions) and the match identity after each iteration.
    (cPecanRealign also reports a likelihood, but adds each matrix's once per diagonal of its band, so it is not
    the likelihood EM maximises and need not rise: the parameters settling is what to watch.) """
    hmm_path = None
    identities = []
    def expectations(b, iteration):
        cigar = os.path.join(batch_dir, '{}.cigar'.format(b))
        if os.path.getsize(cigar) == 0:
            return None
        out = os.path.join(batch_dir, '{}.{}.exp'.format(b, iteration))
        cactus_call(parameters=['cPecanRealign', os.path.join(batch_dir, '{}.fa'.format(b)),
                                '--outputExpectations={}'.format(out)] + (['--loadHmm={}'.format(hmm_path)] if hmm_path else []),
                    infile=cigar, outfile=out + '.stdout', rt_log_cmd=False)
        return read_expectations(out)
    for iteration in range(iterations):
        transitions, emissions = [0.0] * 25, [0.0] * 80
        with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
            for r in pool.map(lambda b: expectations(b, iteration), range(batches)):
                if r is not None:
                    transitions = [a + b for a, b in zip(transitions, r[0])]
                    emissions = [a + b for a, b in zip(emissions, r[1])]
        if sum(emissions) == 0:
            raise RuntimeError('cPecanRealign counted nothing in the alignments in {}'.format(batch_dir))
        transitions, emissions = maximise(transitions, emissions, train_long_gap)
        hmm_path = os.path.join(batch_dir, 'model.{}.hmm'.format(iteration))
        write_hmm(hmm_path, transitions, emissions)
        identities.append(sum(emissions[MATCH * 16 + x * 5] for x in range(4)))
    return transitions, emissions, identities


def combined_gap(transitions):
    """ the model's two gap states (in one sequence) as one affine gap, (open, extend) probabilities: open their total,
    extension the one that keeps their mean gap length -- for when EM has made both of them short gaps """
    T = lambda i, j: transitions[i * 5 + j]
    ds, es, dl, el = T(MATCH, SHORT_GAP_X), T(SHORT_GAP_X, SHORT_GAP_X), T(MATCH, LONG_GAP_X), T(LONG_GAP_X, LONG_GAP_X)
    mean_length = (ds / (1 - es) + dl / (1 - el)) / (ds + dl)
    return ds + dl, 1 - 1 / mean_length
