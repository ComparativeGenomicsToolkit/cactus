/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * poaAligner.h
 *
 *  Created on: 9 Oct 2020
 *      Author: glennhickey
 */

#ifndef POA_ALIGNER_H_
#define POA_ALIGNER_H_

#include "sonLib.h"
#include "multipleAligner.h"

/**
 * Drop-in replacement for makeAlignment() (in multipleAligner.h) that is based on an external POA library.
 * For larger numbers of sequences, it can be orders of magnitude faster than Pecan.  But... the maximum
 * length needs to be bounded or it'll just crash out right away trying to allocated a big matrix.   
 */
MultipleAlignment *makePartialOrderAlignment(StateMachine *sM, stList *seqFrags,
                                             float matchGamma,
                                             PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters);




/**
 * Convert the gapped matrix structure that abpoa returns into list of aligned pairs.  The algorithm used is
 * - For each column in the MSA
 *   - scan until the first non-gapped entry, mark its sequence as the anchor
 *   - keep scanning, reporting a pairwise alignment of each subsequent sequence and the anchor if non-gapped
 * So the number of pairs here is O(N * L) (or linear in the size of the MSA matrix)
 * The all pairs option enumerates all possible pairs of aligned positions. If this is zero then only a linear number
 * of alignment pairs are constructed per column.
 */
stList *poaMatrixToAlignedPairs(uint8_t** msaSeq, int numSeqs, int msaWidth, int score, stList* seqFrags, bool allPairs);

#endif
