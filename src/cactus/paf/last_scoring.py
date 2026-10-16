#!/usr/bin/env python3

"""
Methods for training, reading and applying the scoring models of cactus-pangenome's abPOA alignments: trained with
cPecan's pair HMM (cactus.paf.pair_hmm_training), or given as a .train file (cactus's own, or last-train's)

Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey

Released under the MIT license, see LICENSE.txt
"""

from toil.job import Job
from toil.statsAndLogging import logger
from toil.lib.bioio import getLogLevelString
from sonLib.bioio import newickTreeParser
from toil.realtimeLogger import RealtimeLogger
import os
import math
from cactus.shared.common import cactus_call, getOptionalAttrib
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.common import findRequiredNode


def parse_train_file(train_file_path):
    """ read the .train file (cactus's, write_pair_hmm_train_file, or last-train's) into a dict
    subsitions are in, ex dict['A']['C'] and gaps are in
    dict['GAP-OPEN'] and dict['GAP-EXTEND'], and a trained second gap piece (cactus's
    #cactus -a2/-b2 lines) in dict['GAP-OPEN-2'] and dict['GAP-EXTEND-2'], and the single affine gap
    for aligners that take one piece (#cactus -as/-bs) in dict['GAP-OPEN-SINGLE'] and dict['GAP-EXTEND-SINGLE']
    to keep things simple, only symmetric matrices are accepted"""
    score_dict = {}
    with open(train_file_path, 'r') as train_file:
        for line in train_file:
            if line.startswith('#last -a') or line.startswith('#last -A'):
                key = 'GAP-OPEN'
                val = int(line.split()[-1].strip())
                if key in score_dict and score_dict[key] != val:
                    raise RuntimeError('Asymmetric gap score detected in {}: please use --gapsym with last-train'.format(
                        train_file_path))
                score_dict[key] = val
            elif line.startswith('#last -b') or line.startswith('#last -B'):
                key = 'GAP-EXTEND'
                val = int(line.split()[-1].strip())
                if key in score_dict and score_dict[key] != val:
                    raise RuntimeError('Asymmetric gap extend score detected in {}: please use --gapsym with last-train'.format(
                        train_file_path))
                score_dict[key] = val
            elif line.startswith('#cactus -a2') or line.startswith('#cactus -b2'):
                # the second piece of abPOA's convex gap, trained (see write_pair_hmm_train_file)
                score_dict['GAP-OPEN-2' if line.startswith('#cactus -a2') else 'GAP-EXTEND-2'] = int(line.split()[-1])
            elif line.startswith('#cactus -as') or line.startswith('#cactus -bs'):
                # the model's gaps as a single affine piece, for minigraph and minipoa
                score_dict['GAP-OPEN-SINGLE' if line.startswith('#cactus -as') else 'GAP-EXTEND-SINGLE'] = int(line.split()[-1])
            elif line[0] in ['A', 'C', 'G', 'T']:
                try:
                    row_toks = line.strip().split()
                    assert len(row_toks) == 5
                    key = line[0]
                    assert key not in score_dict
                    row_dict = { 'A' : int(row_toks[1]),
                                 'C' : int(row_toks[2]),
                                 'G' : int(row_toks[3]),
                                 'T' : int(row_toks[4]) }
                    score_dict[key] = row_dict
                except:
                    raise RuntimeError('Error parsing score matrix from {}'.format(train_file_path))
        for key in ['GAP-OPEN', 'GAP-EXTEND', 'A', 'C', 'G', 'T']:
            if key not in score_dict:
                raise RuntimeError('Information for {} not parsed from {}'.format(key, train_file_path))
        # apply_long_gap() needs an extension of at least 1 to scale by, and a model with a negative gap
        # open or a match score of 0 or less is not one the aligners can use.  Raising here makes
        # last_train() fall back to the defaults (or a borrowed model) rather than hand such a model on.
        if score_dict['GAP-EXTEND'] < 1 or score_dict['GAP-OPEN'] < 0 or min(score_dict[b][b] for b in 'ACGT') <= 0:
            raise RuntimeError('Unusable scoring model in {}: gap open {}, gap extend {}, matches {}'.format(
                train_file_path, score_dict['GAP-OPEN'], score_dict['GAP-EXTEND'], [score_dict[b][b] for b in 'ACGT']))
        for i in ['A', 'C', 'G', 'T']:
            rci = { 'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A' }[i]
            for j in ['A', 'C', 'G', 'T']:
                if score_dict[i][j] != score_dict[j][i]:
                    raise RuntimeError('Asymmetric score matrix detected in {}: please use --matsym with last-train'.format(
                        train_file_path))
                rcj = { 'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A' }[j]
                if score_dict[i][j] != score_dict[rci][rcj]:
                    raise RuntimeError('Reverse complement assymetry detected in score matrix in {}: please use --revsym  with last-train'.format(
                        train_file_path))

    return score_dict

def apply_long_gap(score_dict, open_factor, extend_factor):
    """ make a long gap open that's open_factor more expensive to open, but extend_factor cheaper to extend

    last-train's scores are small integers (its gap extension is typically 1), so the long-gap
    extension can only be made extend_factor times cheaper after scaling the whole model up by
    extend_factor.  The substitution scores and both gap costs are scaled together, which leaves the
    trained model unchanged, and the long-gap extension is then exactly the trained one.  (The
    matrix used to be scaled only when the trained extension was below extend_factor while the gaps
    were always scaled, so a model with a larger extension reached the aligner with its gaps
    extend_factor times too dear relative to its substitutions.) """
    assert open_factor > 1 and extend_factor >= 1
    for i in ['A', 'C', 'G', 'T']:
        for j in ['A', 'C', 'G', 'T']:
            score_dict[i][j] *= extend_factor
    score_dict['GAP-OPEN'] *= extend_factor
    score_dict['GAP-EXTEND'] *= extend_factor

    score_dict['GAP-OPEN-2'] = score_dict['GAP-OPEN'] * open_factor
    score_dict['GAP-EXTEND-2'] = max(1, score_dict['GAP-EXTEND'] // extend_factor)
    return reduce_scores(score_dict)

def reduce_scores(score_dict):
    """ divide every score of a model (substitutions and both gap pieces) by their greatest common
    divisor, and return it.  The aligners' choices, ties included, are unchanged by a common positive
    factor, so the alignments are identical with smaller numbers.  apply_long_gap() scales by
    extend_factor (30) to make room for a long-gap extension 30 times cheaper, and a learned extension
    that shares a factor with it leaves every score that many times larger than needed: a learned
    extension of 2 reached abPOA as 60 (GapOpen 1290; GapExtend 60; GapOpen2 3870; GapExtend2 2), more
    than the abPOA cactus pinned before could hold in its 16-bit DP (the HPRC oct9 CHM13 chr2 model
    aborted with "Error in cg_backtrack", and misaligned windows that did not abort).  Divided by 2 it
    is 645, 30, 1935, 1.  Extensions of 1 (all but one model so far) are left alone. """
    bases = ['A', 'C', 'G', 'T']
    gaps = ['GAP-OPEN', 'GAP-EXTEND', 'GAP-OPEN-2', 'GAP-EXTEND-2'] + \
        [k for k in ('GAP-OPEN-SINGLE', 'GAP-EXTEND-SINGLE') if k in score_dict]
    d = 0
    for v in [score_dict[i][j] for i in bases for j in bases] + [score_dict[k] for k in gaps]:
        d = math.gcd(d, abs(v))
    if d > 1:
        for i in bases:
            for j in bases:
                score_dict[i][j] //= d
        for k in gaps:
            score_dict[k] //= d
    return max(d, 1)

def abpoa_exact_depth(matrix, gap_open1, gap_ext1, gap_open2, gap_ext2, window):
    """ the deepest graph (rows on a path) for which abPOA can prove, up front, that a window-long
    sequence aligns exactly in 32 bits, in its worst case (AVX-512 lanes, a band narrower than the
    window).  This mirrors simd_abpoa_score_range() in cactus's abPOA: deeper graphs still align
    exactly, but abPOA has to check every score as it goes, and aborts with a message rather than
    give a wrong alignment if the scores ever outgrow 32 bits.  Trained models allow over 100
    window-lengths of depth; much less means the scores are far larger than they need to be.
    Path scores are left out, as cactus does not use them. """
    M, X = max(0, max(matrix)), max(0, -min(matrix))
    convex = gap_open1 > 0 and gap_open2 > 0  # abPOA's affine and linear modes use one piece
    oe = [gap_open1 + gap_ext1] + ([gap_open2 + gap_ext2] if convex else [])
    OE, EM, lanes = max(oe), max(gap_ext1, gap_ext2) if convex else gap_ext1, 16
    floor = -2**31 + (window + 2) * X + OE + lanes // 2 * EM
    ztop = floor + (window + 1) * M + M + 2 * X + 2 * OE
    # every real score is at least -max(X, min(oe)) * (window + depth)
    return max(0, (-ztop - 1) // max(1, X, min(oe)) - window)

def apply_scores_to_config(score_dict, config_xml, label=None):
    """ load the score dict into the config.  A model with only one gap piece (last-train's) gets a second, long-gap
    piece made from it with the <poa> TrainedGap*2Factor attributes; one with its own (cactus's) is used as it is.
    note that they are disabled with 0, untouched with None.
    label (the chromosome, say) only goes into the log, so a model can be traced back to its job """

    bar_node = findRequiredNode(config_xml, "bar")
    poa_node = findRequiredNode(bar_node, "poa")

    # todo: abPOA doesn't seem to be stable unless a reasonable long gap is set, ie something
    # that's more expensive to open and cheaper to extend.  We set this here using a factor
    # that's applied to the trained gap parameters
    long_gap_open_factor = int(poa_node.attrib['partialOrderAlignmentTrainedGapOpen2Factor'])
    long_gap_extend_factor = int(poa_node.attrib['partialOrderAlignmentTrainedGapExtension2Factor'])
    # the open factor used to be passed the extend factor, so TrainedGapOpen2Factor was read and
    # then discarded.  That went unnoticed while both defaulted to 3; the defaults are now 3 for the
    # open and 30 for the extension.
    trained_ext = score_dict['GAP-EXTEND']
    # a model with a second gap piece of its own is at abPOA's scale already (write_pair_hmm_train_file)
    two_pieces = 'GAP-OPEN-2' in score_dict
    divisor = reduce_scores(score_dict) if two_pieces else apply_long_gap(score_dict, long_gap_open_factor, long_gap_extend_factor)
    if divisor > 1 and not two_pieces:
        RealtimeLogger.info('Trained scoring model{} learned a gap extension of {}: its scores, multiplied by {} for the '
                            'long-gap piece, are all divided by their common factor {} (the alignments are the same, '
                            'the numbers smaller)'.format(' for ' + label if label else '', trained_ext,
                                                          long_gap_extend_factor, divisor))
    
    poa_node.attrib['partialOrderAlignmentGapOpenPenalty1'] = str(score_dict['GAP-OPEN'])
    poa_node.attrib['partialOrderAlignmentGapExtensionPenalty1'] = str(score_dict['GAP-EXTEND'])
    poa_node.attrib['partialOrderAlignmentGapOpenPenalty2'] = str(score_dict['GAP-OPEN-2'])
    poa_node.attrib['partialOrderAlignmentGapExtensionPenalty2'] = str(score_dict['GAP-EXTEND-2'])

    mismatch_scores = []
    match_scores = []
    for i in ['A', 'C', 'G', 'T']:
        for j in ['A', 'C', 'G', 'T']:
            if i != j:
                mismatch_scores.append(score_dict[i][j])
            else:
                match_scores.append(score_dict[i][j])

    # use min / max (match / mismatch) scores for N -- don't really have a better idea
    max_mismatch = max(mismatch_scores)
    min_match = min(match_scores)
            
    score_string = ''
    for i in ['A', 'C', 'G', 'T']:
        for j in ['A', 'C', 'G', 'T']:
            score_string += str(score_dict[i][j]) + ' '
        score_string += str(max_mismatch) + ' '
    for i in range(4):
        score_string += str(max_mismatch) + ' '
    # todo: do we want option to put a mismatch here
    # it seems like for pangenomes in particular we probably don't want N alignment
    score_string += str(min_match)
            
    poa_node.attrib['partialOrderAlignmentSubMatrix'] = score_string

    window = int(poa_node.attrib.get('partialOrderAlignmentWindow', 10000))
    depth = abpoa_exact_depth([int(v) for v in score_string.split()], score_dict['GAP-OPEN'], score_dict['GAP-EXTEND'],
                              score_dict['GAP-OPEN-2'], score_dict['GAP-EXTEND-2'], window)
    if depth < 50 * window:
        RealtimeLogger.warning('Trained scoring model{} is very large (gap open {}, gap extend {}, matrix {}): abPOA can only '
                               'guarantee exact {} bp alignments to graphs up to {} deep, and will check deeper ones as it goes, '
                               'aborting if they overflow.  Consider --scoresFile or lastTrain="0"'.format(
                                   ' for ' + label if label else '', score_dict['GAP-OPEN'], score_dict['GAP-EXTEND'],
                                   score_string, window, depth))

    # minipoa gets the learned gaps as a single affine piece: last-train fits a single affine gap model,
    # and a cactus model has its gaps as one piece as well (GAP-OPEN-SINGLE, both of its gap states as one).
    # The second piece abPOA gets is not handed on: last-train's is made from the first, as apply_long_gap
    # makes one, because abPOA is unstable without one.
    # minipoa's own second piece is turned off too: its shipped values are priced against the
    # shipped first piece, not the learned one.  minipoaSubMatrix is left alone: empty
    # means inherit <poa>'s matrix, which is the learned one we just wrote, so the trained
    # substitution scores reach minipoa too.
    minipoa_node = bar_node.find("minipoa")
    if minipoa_node is not None:
        minipoa_node.attrib['minipoaGapOpenPenalty'] = str(score_dict.get('GAP-OPEN-SINGLE', score_dict['GAP-OPEN']))
        minipoa_node.attrib['minipoaGapExtensionPenalty'] = str(score_dict.get('GAP-EXTEND-SINGLE', score_dict['GAP-EXTEND']))
        minipoa_node.attrib['minipoaGapOpenPenalty2'] = '0'
        minipoa_node.attrib['minipoaGapExtensionPenalty2'] = '0'
        RealtimeLogger.info("Overriding minipoa scores with trained values: GapOpen {}; GapExtend {} (single affine, as trained)".format(
            minipoa_node.attrib['minipoaGapOpenPenalty'],
            minipoa_node.attrib['minipoaGapExtensionPenalty']))

    # the line keeps its old prefix, which logs are grepped for; the label and scaling go on the end
    if two_pieces:
        scaling = 'both gap pieces as trained'
    else:
        scaling = 'learned gap extension {}; all scores multiplied by {}'.format(trained_ext, long_gap_extend_factor)
    RealtimeLogger.info("Overriding abPOA scores with trained values: GapOpen {}; GapExtend {}; GapOpen2 {}; GapExtend2 {}; SubMatrix {} ({}{}{})".format(
        poa_node.attrib['partialOrderAlignmentGapOpenPenalty1'],
        poa_node.attrib['partialOrderAlignmentGapExtensionPenalty1'],
        poa_node.attrib['partialOrderAlignmentGapOpenPenalty2'],
        poa_node.attrib['partialOrderAlignmentGapExtensionPenalty2'],
        poa_node.attrib['partialOrderAlignmentSubMatrix'], label + ': ' if label else '', scaling,
        ', then divided by their common factor {}'.format(divisor) if divisor > 1 else ''))

def minigraph_wfa_penalties(score_dict):
    """ convert a trained model into minigraph's base-alignment penalties, as (x, o, e) for
    its --wfa-pen option.

    minigraph fills the gaps between its chain anchors with WFA, which scores a match as 0, so only
    one match score M (the mean of the diagonal) and one mismatch score X (the mean of the rest)
    survive: there is no room for the transition/transversion split.  The Smith-Waterman to WFA
    transform is then x = 2(M + X), o = 2*open, e = 2*extend + M (doubled to stay integral).
    WFA's long-gap piece would have no room next to so cheap an extension, so --wfa-pen gets just one piece:
    last-train's, or a cactus model's gaps as a single piece (GAP-OPEN-SINGLE).  A cactus model's first piece
    is its short-gap state alone, 2bp long on average, whose dear extension, as the e everything is divided
    by, took human models from 8-9,11,1 down to 6,7,1, which in an HPRC 30-way left 0.7% more SNV errors.

    WFA works through every score up to the alignment's, so its cost grows with the size of the
    penalties: they are scaled so that e, much the smallest, is 1.  Human models come out at about
    8,11,1 against minigraph's default of 4,4,2 (plus a 15,1 long-gap piece), which is to say that
    mismatches and gap opens cost a bit more and gap extension an order of magnitude less.

    score_dict is as parse_train_file() returns it, and is not modified """
    bases = ['A', 'C', 'G', 'T']
    match = sum(score_dict[a][a] for a in bases) / 4.
    mismatch = -sum(score_dict[a][b] for a in bases for b in bases if a != b) / 12.
    x = 2. * (match + mismatch)
    gap_open = score_dict.get('GAP-OPEN-SINGLE', score_dict['GAP-OPEN'])
    gap_extend = score_dict.get('GAP-EXTEND-SINGLE', score_dict['GAP-EXTEND'])
    o = 2. * gap_open
    e = 2. * gap_extend + match
    if x <= 0 or o < 0 or e <= 0:
        raise RuntimeError('Scoring model does not convert to minigraph penalties: match {} mismatch {} gap open {} extend {}'.format(
            match, -mismatch, gap_open, gap_extend))
    return max(1, round(x / e)), round(o / e), 1

# lastz's default scores, HOXD70, and the base frequencies its scale is measured with.  A trained model is written
# at the same scale (lambda), which is the one abPOA's shipped scores are in
HOXD70 = [[91, -114, -31, -123], [-114, 100, -125, -31], [-31, -125, 100, -114], [-123, -31, -114, 91]]
HOXD70_FREQS = [0.26585, 0.23415, 0.23415, 0.26585]

def implied_lambda(matrix, freqs):
    """ the lambda with sum_ij f_i f_j exp(lambda s_ij) = 1: the scale of a score matrix, in nats per unit """
    def excess(lam):
        return sum(freqs[i] * freqs[j] * math.exp(lam * matrix[i][j]) for i in range(4) for j in range(4)) - 1
    lo, hi = 1e-6, 1.0
    if not (excess(lo) < 0 < excess(hi)):
        raise RuntimeError('score matrix {} has no positive lambda'.format(matrix))
    for _ in range(200):
        mid = (lo + hi) / 2
        if excess(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

def write_pair_hmm_train_file(path, transitions, emissions, comment, open2_factor=3, extend2_factor=30, long_gap_trained=True,
                              min_long_gap_mean=10.0):
    """ a pair HMM (pair_hmm_training.train_pair_hmm) as a .train file that parse_train_file reads, at abPOA's scale --
    HOXD70's units, as the shipped <poa> scores are -- rather than last-train's coarse integers, whose gap extension
    can only be 1 or 2 (and which last-train caps at exp(-1/t), which held 17 of the 25 human chromosomes' models of
    the oct9 HPRC build at one value).  Substitutions are log-odds (Durbin et al., ch. 4, end effects ignored):
    ln(p_ab / (f_a f_b)) + ln(m), m the match-continue probability, and a gap state of open and extend probabilities
    delta and epsilon costs d + (k - 1) e for k bases, d = -ln(delta (1 - epsilon) / m) and e = -ln(epsilon), which
    as an open plus k extensions is open = d - e, extend = e.  The two pieces of abPOA's convex gap are written as
    #last -a/-b and #cactus -a2/-b2, which apply_scores_to_config uses as they are, and the gaps as one affine
    piece (both states, pair_hmm_training.combined_gap) as #cactus -as/-bs, which minigraph's penalties
    (minigraph_wfa_penalties) and minipoa use.

    When the long-gap state was trained and came out a real long gap -- dearer to open than the short one, cheaper to
    extend, and with a mean length of at least min_long_gap_mean bases -- the two states are the two pieces.  On
    human chr10 haplotypes, 10kb pieces gave long states of 23-40bp, and those that EM made a second short gap had
    2.4-4.8bp; one of those (1073+22 against 909+48) cost abPOA 20% more events than the factors' piece.  Otherwise
    (it was held at cPecan's default, which is a prior for divergent sequence, or EM made it a second short gap, as
    2kb pieces of human haplotypes do: they hold next to no long indels) the first piece is all the trained gaps as
    one (pair_hmm_training.combined_gap; the short state's alone when the long one was held), and the second is made
    from it as apply_long_gap makes one, open2_factor times dearer to open and extend2_factor times cheaper to
    extend (<poa partialOrderAlignmentTrainedGapOpen2Factor> and ...Extension2Factor). """
    from cactus.paf.pair_hmm_training import MATCH, SHORT_GAP_X, LONG_GAP_X, combined_gap
    lam = implied_lambda(HOXD70, HOXD70_FREQS)
    P = [[emissions[MATCH * 16 + x * 4 + y] for y in range(4)] for x in range(4)]
    f = [sum(row) for row in P]
    m = transitions[MATCH * 5 + MATCH]
    matrix = [[int(round((math.log(P[x][y] / (f[x] * f[y])) + math.log(m)) / lam)) for y in range(4)] for x in range(4)]
    def piece(open_prob, extend_prob):
        d, e = -math.log(open_prob * (1 - extend_prob) / m), -math.log(extend_prob)
        return int(round((d - e) / lam)), max(1, int(round(e / lam)))
    a1, b1 = piece(transitions[MATCH * 5 + SHORT_GAP_X], transitions[SHORT_GAP_X * 5 + SHORT_GAP_X])
    a2, b2 = piece(transitions[MATCH * 5 + LONG_GAP_X], transitions[LONG_GAP_X * 5 + LONG_GAP_X])
    long_mean = 1 / (1 - transitions[LONG_GAP_X * 5 + LONG_GAP_X])
    if not (long_gap_trained and a2 > a1 and b2 < b1 and long_mean >= min_long_gap_mean):
        if long_gap_trained:
            a1, b1 = piece(*combined_gap(transitions))
        a2, b2 = a1 * open2_factor, max(1, b1 // extend2_factor)    # as apply_long_gap divides
    # both gap states as one affine piece, for the aligners that take only one (minigraph's WFA, minipoa): the short
    # state alone is 2bp on average, and as minigraph's gap it priced mismatches and gap opens too low.  Held at
    # cPecan's default, the long state is not the data's, so then it is the first piece
    s1, t1 = piece(*combined_gap(transitions)) if long_gap_trained else (a1, b1)
    with open(path, 'w') as out:
        out.write('# cactus pair HMM model (cPecan EM), at HOXD70\'s scale: {}\n'.format(comment))
        out.write('# substitution percent identity: {:.6}\n'.format(100 * sum(P[x][x] for x in range(4))))
        out.write('#last -a {}\n#last -A {}\n#last -b {}\n#last -B {}\n'.format(a1, a1, b1, b1))
        out.write('#cactus -a2 {}\n#cactus -b2 {}\n'.format(a2, b2))
        out.write('#cactus -as {}\n#cactus -bs {}\n'.format(s1, t1))
        out.write('# score matrix (query letters = columns, reference letters = rows):\n')
        out.write('       A      C      G      T\n')
        for base, row in zip('ACGT', matrix):
            out.write('{} {}\n'.format(base, ' '.join('{:6d}'.format(v) for v in row)))

def train_pair_hmm_file(ref_path, sample_path, train_file, work_dir, iterations=5, threads=1, train_long_gap=True,
                        open2_factor=3, extend2_factor=30, preset='asm5', min_mapq=60):
    """ train a model for aligning to the reference at ref_path and write it as train_file
    (write_pair_hmm_train_file): the pieces at sample_path (pair_hmm_training.sample_pieces) mapped to the reference
    with minimap2 (-x preset, primary alignments of mapping quality min_mapq or more), and cPecan's pair HMM fitted
    to the alignments by EM
    (pair_hmm_training.train_pair_hmm), training its long-gap state too with train_long_gap.  Returns
    (paf_to_cpecan_batches's statistics, the identity by iteration), or raises if there is nothing to train on. """
    from cactus.paf.pair_hmm_training import paf_to_cpecan_batches, train_pair_hmm
    paf_path = os.path.join(work_dir, 'sample.paf')
    cactus_call(parameters=['minimap2', '-c', '--eqx', '--secondary=no', '-x', preset, '-t', str(threads),
                            ref_path, sample_path], outfile=paf_path)
    batches = 4 * threads
    batch_dir = os.path.join(work_dir, 'batches')
    stats = paf_to_cpecan_batches(paf_path, sample_path, ref_path, batch_dir, batches, min_mapq=min_mapq)
    if stats['kept'] == 0:
        raise RuntimeError('none of the {} alignments of the sample is usable'.format(stats['alignments']))
    transitions, emissions, identities = train_pair_hmm(batch_dir, batches, iterations, threads, train_long_gap)
    write_pair_hmm_train_file(train_file, transitions, emissions, '{} alignments, {} bp'.format(
        stats['kept'], stats['aligned_bp']), open2_factor, extend2_factor, train_long_gap)
    return stats, identities

def last_train(job, config, seq_order, seq_id_map, ref_name=None):
    """ train the scoring model for aligning to the reference with cPecan's pair HMM (train_pair_hmm_file): up to
    <graphmap trainPartners> partner genomes (pick_train_partners), <graphmap trainPieces> pieces of
    <graphmap trainPieceLength> from them in all, the long-gap state trained too unless <graphmap trainLongGap="0">.
    The name is last-train's, which this replaces.  Returns the .train file's id, or None when there is nothing to
    train on or the fit fails.

    ref_name names the reference genome when it is not seq_order[0], as when extending an existing
    minigraph where the order holds only the genomes being added """

    if ref_name is None:
        assert len(seq_order) > 1
        ref_name = seq_order[0]
    elif ref_name not in seq_order:
        seq_order = [ref_name] + list(seq_order)

    name1 = ref_name

    # short circuit if ref sequence is too small to have a hope of training        
    if seq_id_map[name1].size < 500000:
        RealtimeLogger.warning('Input fasta for {} too small to train scoring model on.  Will fall back to defaults'.format(os.path.basename(name1)))
        return None

    graphmap_node = findRequiredNode(config, 'graphmap')
    poa_node = findRequiredNode(findRequiredNode(config, 'bar'), 'poa')
    partners = pick_train_partners(name1, seq_order, seq_id_map,
                                   getOptionalAttrib(graphmap_node, 'trainPartners', typeFn=int, default=4))

    # short circuit if we can't find anything to train on
    if not partners:
        RealtimeLogger.warning('Unable to find sequence to train scoring model on for {}.  Will fall back to defaults'.format(os.path.basename(name1)))
        return None
    names = ' '.join(os.path.basename(p) for p in partners)

    from cactus.paf.pair_hmm_training import sample_pieces
    work_dir = job.fileStore.getLocalTempDir()
    train_file = os.path.join(work_dir, '{}.train'.format(os.path.basename(name1)))
    pieces = getOptionalAttrib(graphmap_node, 'trainPieces', typeFn=int, default=1000)
    piece_length = getOptionalAttrib(graphmap_node, 'trainPieceLength', typeFn=int, default=10000)

    # a model is an optimization, not a requirement: whatever goes wrong here (too little alignment to fit, or a
    # model parse_train_file won't accept), the alignment falls back to the default scores, or to another
    # chromosome's model in batch mode
    try:
        # the pieces, shared between the partners, each read in turn and let go once sampled, so that only one
        # partner is on disk at a time
        sample_path = os.path.join(work_dir, 'sample.fa')
        sampled = 0
        with open(sample_path, 'wb') as sample:
            for i, partner in enumerate(partners):
                part_path = os.path.join(work_dir, 'sample.{}.fa'.format(i))
                sampled += sample_pieces(job.fileStore.readGlobalFile(seq_id_map[partner]), part_path,
                                         max(1, pieces // len(partners)), piece_length, seed=i + 1,
                                         prefix='trainPiece{}_'.format(i))
                job.fileStore.deleteLocalFile(seq_id_map[partner])
                with open(part_path, 'rb') as part:
                    sample.write(part.read())
        if sampled == 0:
            raise RuntimeError('the partners have no {}bp pieces to sample'.format(piece_length))
        stats, identities = train_pair_hmm_file(
            job.fileStore.readGlobalFile(seq_id_map[name1]), sample_path, train_file, work_dir,
            iterations=getOptionalAttrib(graphmap_node, 'trainIterations', typeFn=int, default=5),
            threads=max(1, int(job.cores)),
            train_long_gap=getOptionalAttrib(graphmap_node, 'trainLongGap', typeFn=bool, default=True),
            min_mapq=getOptionalAttrib(graphmap_node, 'trainMinMapQ', typeFn=int, default=60),
            open2_factor=int(poa_node.attrib['partialOrderAlignmentTrainedGapOpen2Factor']),
            extend2_factor=int(poa_node.attrib['partialOrderAlignmentTrainedGapExtension2Factor']))
        score_dict = parse_train_file(train_file)
    except Exception as e:
        RealtimeLogger.warning('Training scoring model for {} against {} failed, so it will fall back to defaults: {}'.format(
            os.path.basename(name1), names, e))
        return None

    RealtimeLogger.info('Trained scoring model for {} against {} (cPecan pair HMM, {} pieces, {} alignments, {} bp; EM '
                        'identity by iteration {}): gap {}+{}k, long gap {}+{}k, matrix {}'.format(
        os.path.basename(name1), names, sampled, stats['kept'], stats['aligned_bp'],
        ' '.join('{:.2f}'.format(100 * i) for i in identities), score_dict['GAP-OPEN'], score_dict['GAP-EXTEND'],
        score_dict['GAP-OPEN-2'], score_dict['GAP-EXTEND-2'],
        ' '.join(str(score_dict[a][b]) for a in 'ACGT' for b in 'ACGT')))
    return job.fileStore.writeGlobalFile(train_file)

def pick_train_partners(name1, seq_order, seq_id_map, count, **kwargs):
    """ up to count genomes to train name1's model against: pick_train_partner's candidates nearest its pick, the middle
    of them, so that the model is not one partner's (in the oct9 HPRC build the same chromosome's model moved from
    run to run with its partner) """
    candidates = pick_train_partner(name1, seq_order, seq_id_map, all_candidates=True, **kwargs)
    if not candidates:
        return []
    middle = len(candidates) // 2
    nearest = sorted(range(len(candidates)), key=lambda i: (abs(i - middle), i))[:count]
    return [candidates[i] for i in sorted(nearest)]

def pick_train_partner(name1, seq_order, seq_id_map, min_size=500000, min_ref_frac=0.5, size_band=2.0, all_candidates=False):
    """ choose the genome to train name1's model against (or, with all_candidates, the candidates, in order).

    The candidates are the genomes big enough to hold a meaningful amount of alignment: over
    min_size and over min_ref_frac of name1.  Of those, only the ones within size_band of their
    median size are kept, so that a fragmented or oversized (contaminated, or whole-genome instead
    of one chromosome) assembly doesn't set the model for everything.  The pick is then the middle
    of these in seq_order, which is by mash distance to the reference when minigraph sorts it: the
    furthest genome is by construction the outlier, and trained models barely depend on the partner
    anyway (gap open 42-45, extend 1 across the HPRC partners tried with last-train), so a typical
    one is the safest bet. """
    ref_size = seq_id_map[name1].size
    candidates = [seq for seq in seq_order if seq != name1 and seq in seq_id_map and
                  seq_id_map[seq].size > min_size and float(seq_id_map[seq].size) / float(ref_size) > min_ref_frac]
    if not candidates:
        return [] if all_candidates else None
    sizes = sorted(seq_id_map[seq].size for seq in candidates)
    median_size = sizes[len(sizes) // 2]
    typical = [seq for seq in candidates if median_size / size_band <= seq_id_map[seq].size <= median_size * size_band]
    if typical:
        candidates = typical
    return candidates if all_candidates else candidates[len(candidates) // 2]

def last_train_enabled(options, config_node):
    """ training is switched on and off by graphmap's lastTrain attribute.  --lastTrain is
    deprecated and only turns it on; --scoresFile turns it off since its model would be ignored """
    graphmap_node = findRequiredNode(config_node, "graphmap")
    if getattr(options, 'lastTrain', False):
        logger.warning('--lastTrain is deprecated: training is now on by default, and is toggled with the lastTrain '
                       'attribute of <graphmap> in the config')
        graphmap_node.attrib['lastTrain'] = '1'
    if getattr(options, 'scoresFile', None):
        return False
    return getOptionalAttrib(graphmap_node, 'lastTrain', typeFn=bool, default=False)
    

    
