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
    """ make a long gap open that's open_factor more expensive to open, but extend_factor cheaper to extend """
    assert open_factor > 1 and extend_factor >= 1
    # multiply everything else so we can make a smaller big gap extend
    if score_dict['GAP-EXTEND'] < extend_factor:
        for i in ['A', 'C', 'G', 'T']:
            for j in ['A', 'C', 'G', 'T']:
                score_dict[i][j] *= extend_factor
    score_dict['GAP-OPEN'] *= extend_factor
    score_dict['GAP-EXTEND'] *= extend_factor

    score_dict['GAP-OPEN-2'] = score_dict['GAP-OPEN'] * open_factor
    score_dict['GAP-EXTEND-2'] = max(1, int(score_dict['GAP-EXTEND'] / extend_factor))            

def apply_scores_to_config(score_dict, config_xml):
    """ load the score dict into the config.  since last won't train long gaps,
    should be specified by the parameters.  note that they are disabled with 0, untouched with None"""

    bar_node = findRequiredNode(config_xml, "bar")
    poa_node = findRequiredNode(bar_node, "poa")

    # todo: abPOA doesn't seem to be stable unless a reasonable long gap is set, ie something
    # that's more expensive to open and cheaper to extend.  We set this here using a factor
    # that's applied to the trained gap parameters
    long_gap_open_factor = int(poa_node.attrib['partialOrderAlignmentTrainedGapOpen2Factor'])
    long_gap_extend_factor = int(poa_node.attrib['partialOrderAlignmentTrainedGapExtension2Factor'])
    apply_long_gap(score_dict, long_gap_extend_factor, long_gap_extend_factor)
    
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

    RealtimeLogger.info("Overriding abPOA scores with trained values: GapOpen {}; GapExtend {}; GapOpen2 {}; GapExtend2 {}; SubMatrix {}".format(
        poa_node.attrib['partialOrderAlignmentGapOpenPenalty1'],
        poa_node.attrib['partialOrderAlignmentGapExtensionPenalty1'],
        poa_node.attrib['partialOrderAlignmentGapOpenPenalty2'],
        poa_node.attrib['partialOrderAlignmentGapExtensionPenalty2'],
        poa_node.attrib['partialOrderAlignmentSubMatrix']))

def last_train(job, config, seq_order, seq_id_map):
    """ run last_train on a pair of fasta files, using the first as the database """

    assert len(seq_order) > 1

    name1 = seq_order[0]
    name2 = None

    # short circuit if ref sequence is too small to have a hope of training        
    if seq_id_map[name1].size < 500000:
        RealtimeLogger.warning('Input fasta for {} too small to train scoring model on.  Will fall back to defaults'.format(os.path.basename(name1)))
        return None
    
    # determine sequence to compare to compare
    # it's the furthest in the order that's both greater than 500k and 50% of the size of the first
    rev_order = [seq for seq in reversed(seq_order)]
    for seq in rev_order[1:]:
        if seq != name1 and seq_id_map[seq].size > 500000 and (float(seq_id_map[seq].size) / float(seq_id_map[name1].size) > 0.5):
            name2 = seq
            break

    # short circuit if we can't find anything to train on
    if name2 is None:
       RealtimeLogger.warning('Unable to find sequence to train scoring model on for {}.  Will fall back to defaults'.format(os.path.basename(name1)))
       return None
                
    # sometimes the s
    fa1_id, fa2_id = seq_id_map[name1], seq_id_map[name2]
    work_dir = job.fileStore.getLocalTempDir()
    fa1_path = os.path.join(work_dir, name1 + '.fa')
    job.fileStore.readGlobalFile(fa1_id, fa1_path)
    fa2_path = os.path.join(work_dir, name2 + '.fa')
    job.fileStore.readGlobalFile(fa2_id, fa2_path)
        
    # make the database
    cactus_call(parameters=['lastdb', name1 + '_db', fa1_path, '-P', str(job.cores)])

    # do the training
    # note: there are some specific options for distant genomes that should be
    # incorporated if/when this ever gets used in progressive cactus
    train_cmd = ['last-train', '--revsym', '--matsym', '--gapsym',
                 '-P', str(job.cores), name1 + '_db', fa2_path]

    train_file = os.path.join(work_dir, '{}_{}.train'.format(name1, name2))
    cactus_call(parameters=train_cmd, outfile=train_file)

    return job.fileStore.writeGlobalFile(train_file)
    

    
