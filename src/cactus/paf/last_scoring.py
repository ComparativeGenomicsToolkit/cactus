#!/usr/bin/env python3

"""
Methods for reading and computing pairwise alignment scoring matrices using last

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
    """ read the .train file (from last-train) into a dict
    subsitions are in, ex dict['A']['C'] and gaps are in
    dict['GAP-OPEN'] and dict['GAP-EXTEND']
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
    gaps = ['GAP-OPEN', 'GAP-EXTEND', 'GAP-OPEN-2', 'GAP-EXTEND-2']
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
    """ load the score dict into the config.  since last won't train long gaps,
    should be specified by the parameters.  note that they are disabled with 0, untouched with None.
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
    divisor = apply_long_gap(score_dict, long_gap_open_factor, long_gap_extend_factor)
    if divisor > 1:
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

    # minipoa gets the learned gaps verbatim, as a single affine piece: last-train fits a single
    # affine gap model, and the GAP-OPEN-2/GAP-EXTEND-2 pair above is a synthesised second piece
    # that exists only because abPOA is unstable without one, so it would be wrong to hand it on.
    # minipoa's own second piece is turned off too: its shipped values are priced against the
    # shipped first piece, not the learned one.  minipoaSubMatrix is left alone: empty
    # means inherit <poa>'s matrix, which is the learned one we just wrote, so the trained
    # substitution scores reach minipoa too.
    minipoa_node = bar_node.find("minipoa")
    if minipoa_node is not None:
        minipoa_node.attrib['minipoaGapOpenPenalty'] = str(score_dict['GAP-OPEN'])
        minipoa_node.attrib['minipoaGapExtensionPenalty'] = str(score_dict['GAP-EXTEND'])
        minipoa_node.attrib['minipoaGapOpenPenalty2'] = '0'
        minipoa_node.attrib['minipoaGapExtensionPenalty2'] = '0'
        RealtimeLogger.info("Overriding minipoa scores with trained values: GapOpen {}; GapExtend {} (single affine, as trained)".format(
            minipoa_node.attrib['minipoaGapOpenPenalty'],
            minipoa_node.attrib['minipoaGapExtensionPenalty']))

    # the line keeps its old prefix, which logs are grepped for; the label and scaling go on the end
    RealtimeLogger.info("Overriding abPOA scores with trained values: GapOpen {}; GapExtend {}; GapOpen2 {}; GapExtend2 {}; SubMatrix {} ({}learned gap extension {}; all scores multiplied by {}{})".format(
        poa_node.attrib['partialOrderAlignmentGapOpenPenalty1'],
        poa_node.attrib['partialOrderAlignmentGapExtensionPenalty1'],
        poa_node.attrib['partialOrderAlignmentGapOpenPenalty2'],
        poa_node.attrib['partialOrderAlignmentGapExtensionPenalty2'],
        poa_node.attrib['partialOrderAlignmentSubMatrix'], label + ': ' if label else '', trained_ext, long_gap_extend_factor,
        ', then divided by their common factor {}'.format(divisor) if divisor > 1 else ''))

def minigraph_wfa_penalties(score_dict):
    """ convert a last-train model into minigraph's base-alignment penalties, as (x, o, e) for
    its --wfa-pen option.

    minigraph fills the gaps between its chain anchors with WFA, which scores a match as 0, so only
    one match score M (the mean of the diagonal) and one mismatch score X (the mean of the rest)
    survive: there is no room for the transition/transversion split.  The Smith-Waterman to WFA
    transform is then x = 2(M + X), o = 2*open, e = 2*extend + M (doubled to stay integral).
    last-train fits a single affine gap, and with gap extension already cheap next to M (1 against
    ~6.5 on human), a cheaper long-gap piece would have no room either, so --wfa-pen gets just the
    one piece.

    WFA works through every score up to the alignment's, so its cost grows with the size of the
    penalties: they are scaled so that e, much the smallest, is 1.  Human models come out at about
    8,11,1 against minigraph's default of 4,4,2 (plus a 15,1 long-gap piece), which is to say that
    mismatches and gap opens cost a bit more and gap extension an order of magnitude less.

    score_dict is as parse_train_file() returns it, and is not modified """
    bases = ['A', 'C', 'G', 'T']
    match = sum(score_dict[a][a] for a in bases) / 4.
    mismatch = -sum(score_dict[a][b] for a in bases for b in bases if a != b) / 12.
    x = 2. * (match + mismatch)
    o = 2. * score_dict['GAP-OPEN']
    e = 2. * score_dict['GAP-EXTEND'] + match
    if x <= 0 or o < 0 or e <= 0:
        raise RuntimeError('Scoring model does not convert to minigraph penalties: match {} mismatch {} gap open {} extend {}'.format(
            match, -mismatch, score_dict['GAP-OPEN'], score_dict['GAP-EXTEND']))
    return max(1, round(x / e)), round(o / e), 1

def last_train(job, config, seq_order, seq_id_map, ref_name=None):
    """ run last_train on a pair of fasta files, using the first as the database.

    ref_name names the database genome when it is not seq_order[0], as when extending an existing
    minigraph where the order holds only the genomes being added """

    if ref_name is None:
        assert len(seq_order) > 1
        ref_name = seq_order[0]
    elif ref_name not in seq_order:
        seq_order = [ref_name] + list(seq_order)

    name1 = ref_name
    name2 = None

    # short circuit if ref sequence is too small to have a hope of training        
    if seq_id_map[name1].size < 500000:
        RealtimeLogger.warning('Input fasta for {} too small to train scoring model on.  Will fall back to defaults'.format(os.path.basename(name1)))
        return None
    
    name2 = pick_train_partner(name1, seq_order, seq_id_map)

    # short circuit if we can't find anything to train on
    if name2 is None:
       RealtimeLogger.warning('Unable to find sequence to train scoring model on for {}.  Will fall back to defaults'.format(os.path.basename(name1)))
       return None

    fa1_id, fa2_id = seq_id_map[name1], seq_id_map[name2]
    work_dir = job.fileStore.getLocalTempDir()
    fa1_path = os.path.join(work_dir, name1 + '.fa')
    job.fileStore.readGlobalFile(fa1_id, fa1_path)
    fa2_path = os.path.join(work_dir, name2 + '.fa')
    job.fileStore.readGlobalFile(fa2_id, fa2_path)

    # note: there are some specific options for distant genomes that should be
    # incorporated if/when this ever gets used in progressive cactus
    train_cmd = ['last-train', '--revsym', '--matsym', '--gapsym',
                 '-P', str(job.cores), name1 + '_db', fa2_path]
    train_file = os.path.join(work_dir, '{}_{}.train'.format(name1, name2))

    # a model is an optimization, not a requirement: whatever goes wrong here (last-train failing
    # to converge on too little alignment, or writing something parse_train_file won't accept),
    # the alignment falls back to the default scores, or to another chromosome's model in batch mode
    try:
        cactus_call(parameters=['lastdb', name1 + '_db', fa1_path, '-P', str(job.cores)])
        cactus_call(parameters=train_cmd, outfile=train_file)
        parse_train_file(train_file)
    except Exception as e:
        RealtimeLogger.warning('Training scoring model for {} against {} failed, so it will fall back to defaults: {}'.format(
            os.path.basename(name1), os.path.basename(name2), e))
        return None

    RealtimeLogger.info('Trained scoring model for {} against {}'.format(os.path.basename(name1), os.path.basename(name2)))
    return job.fileStore.writeGlobalFile(train_file)

def pick_train_partner(name1, seq_order, seq_id_map, min_size=500000, min_ref_frac=0.5, size_band=2.0):
    """ choose the genome to train name1's model against.

    The candidates are the genomes big enough to hold a meaningful amount of alignment: over
    min_size and over min_ref_frac of name1.  Of those, only the ones within size_band of their
    median size are kept, so that a fragmented or oversized (contaminated, or whole-genome instead
    of one chromosome) assembly doesn't set the model for everything.  The pick is then the middle
    of these in seq_order, which is by mash distance to the reference when minigraph sorts it: the
    furthest genome is by construction the outlier, and trained models barely depend on the partner
    anyway (gap open 42-45, extend 1 across the HPRC partners tried), so a typical one is the
    safest bet. """
    ref_size = seq_id_map[name1].size
    candidates = [seq for seq in seq_order if seq != name1 and seq in seq_id_map and
                  seq_id_map[seq].size > min_size and float(seq_id_map[seq].size) / float(ref_size) > min_ref_frac]
    if not candidates:
        return None
    sizes = sorted(seq_id_map[seq].size for seq in candidates)
    median_size = sizes[len(sizes) // 2]
    typical = [seq for seq in candidates if median_size / size_band <= seq_id_map[seq].size <= median_size * size_band]
    if typical:
        candidates = typical
    return candidates[len(candidates) // 2]

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
    

    
