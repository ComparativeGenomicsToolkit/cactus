/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * flowerAligner.h
 *
 *  Created on: 2 Jul 2010
 *      Author: benedictpaten
 */

#ifndef FLOWER_ALIGNER_H_
#define FLOWER_ALIGNER_H_

#include "pairwiseAligner.h"

/*
 * Constructs an alignment for the flower by constructing an alignment for each end
 * then filtering the alignments against each other so each position is a member of only one
 * end alignment. Spanning trees controls the number of pairwise alignments used
 * to construct the alignment, maxSequenceLength is the maximum length of a sequence to consider in the end alignment.
 * Model parameters is the parameters of the pairwise alignment model.
 */
stSortedSet *makeFlowerAlignment(StateMachine *sM, Flower *flower, int64_t spanningTrees,
        int64_t maxSequenceLength, bool useProgressiveMerging, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters, bool pruneOutStubAlignments);

/*
 * As above, but including alignments from disk.
 */
stSortedSet *makeFlowerAlignment3(StateMachine *sM, Flower *flower, stList *listOfEndAlignmentFiles, int64_t spanningTrees,
        int64_t maxSequenceLength, bool useProgressiveMerging, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters, bool pruneOutStubAlignments);

/*
 * Returns an end, if exists, that has cap involved in every adjacency, else returns null.
 * In case of tie, choose the end that would align forward along refEventName
 * If the dominant end is still flipped along the reference, set outIsFlipped (if it isn't NULL) to true
 */
End *getDominantEnd(Flower *flower, Name refEventName, bool* outIsFlipped);

/*
 * Ascertain which ends should be aligned separately.
 */
stSortedSet *getEndsToAlignSeparately(Flower *flower, int64_t maxSequenceLength, int64_t largeEndSize);

/*
 * The total number of unaligned bases in adjacencies incident with the end.
 */
int64_t getTotalAdjacencyLength(End *end);

#endif /* NETALIGNER_H_ */
