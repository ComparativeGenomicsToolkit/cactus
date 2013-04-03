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
stSortedSet *makeFlowerAlignment(Flower *flower, int32_t spanningTrees,
        int32_t maxSequenceLength, int32_t maximumNumberOfSequencesBeforeSwitchingToFast, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters, bool pruneOutStubAlignments);

/*
 * As above, but including alignments from disk.
 */
stSortedSet *makeFlowerAlignment3(Flower *flower, stList *listOfEndAlignmentFiles, int32_t spanningTrees,
        int32_t maxSequenceLength, int32_t maximumNumberOfSequencesBeforeSwitchingToFast, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters, bool pruneOutStubAlignments);

void writeEndAlignmentToDisk(End *end, stSortedSet *endAlignment, FILE *fileHandle);

stSortedSet *loadEndAlignmentFromDisk(FILE *fileHandle, End **end);

stSortedSet *getEndsToAlignSeparately(Flower *flower, int32_t maxSequenceLength, int32_t largeEndSize);

#endif /* NETALIGNER_H_ */
