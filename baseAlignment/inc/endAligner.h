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
    Name sequence;
    int32_t position;
    bool strand;
    int32_t score;
    struct _AlignedPair *reverse;
} AlignedPair;

/*
 * Constructs the an aligned pair.
 */
AlignedPair *alignedPair_construct(Name sequence1, int32_t position1, bool strand1,
                                   Name sequence2, int32_t position2, bool strand2, int32_t score);

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
stSortedSet *makeEndAlignment(End *end, int32_t spanningTrees, int32_t maxSequenceLength, float gapGamma, bool useBanding,
        PairwiseAlignmentBandingParameters *pairwiseAlignmentBandingParameters);

#endif /* ENDALIGNER_H_ */
