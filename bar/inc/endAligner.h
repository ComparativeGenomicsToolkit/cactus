/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * endAligner.h
 *
 *  Created on: 1 Jul 2010
 *      Author: benedictpaten
 */

#ifndef ENDALIGNER_H_
#define ENDALIGNER_H_

#include "sonLib.h"
#include "cactus.h"
#include "pairwiseAligner.h"

typedef struct _AlignedPair {
    int64_t subsequenceIdentifier;
    int64_t position;
    bool strand;
    int64_t score;
    struct _AlignedPair *reverse;
} AlignedPair;

/*
 * Constructs the an aligned pair.
 */
AlignedPair *alignedPair_construct(int64_t subsequenceIdentifier1, int64_t position1, bool strand1,
        int64_t subsequenceIdentifier2, int64_t position2, bool strand2, int64_t score1, int64_t score2);

/*
 * Destruct the aligned pair.
 */
void alignedPair_destruct(AlignedPair *alignedPair);

/*
 * Compares two aligned pairs.
 */
int alignedPair_cmpFn(const AlignedPair *alignedPair1, const AlignedPair *alignedPair2);

/*
 * Creates a global alignment (as a set of aligned pairs) of the sequences from the end,
 * the pairs returned are ordered according
 * to the alignerPair comparison function.
 */
stSortedSet *makeEndAlignment(StateMachine *sM, End *end, int64_t spanningTrees, int64_t maxSequenceLength,
                              bool useProgressiveMerging, float gapGamma,
                              PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters,
                              int64_t poaWindow);

/*
 * Writes an end alignment to the given file.
 */
void writeEndAlignmentToDisk(End *end, stSortedSet *endAlignment, FILE *fileHandle);

/*
 * Loads an end alignment from the given file.
 */
stSortedSet *loadEndAlignmentFromDisk(Flower *flower, FILE *fileHandle, End **end);


#endif /* ENDALIGNER_H_ */
