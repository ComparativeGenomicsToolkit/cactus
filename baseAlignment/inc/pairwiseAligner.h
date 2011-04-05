/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * pairwiseAligner.h
 *
 *  Created on: 5 Jul 2010
 *      Author: benedictpaten
 */

#ifndef PAIRWISEALIGNER_H_
#define PAIRWISEALIGNER_H_

//Constant that gives the integer value equal to probability 1. Integer probability zero is always 0.
#define PAIR_ALIGNMENT_PROB_1 10000000

/*
 * Gets the set of posterior match probabilities under a simple HMM model of alignment for two DNA sequences.
 */
stList *getAlignedPairs(const char *string1, const char *string2);

typedef  struct _pairwiseAlignmentBandingParameters {
    int32_t maxBandingSize; //No matrix biggger than this number squared will be computed
    int32_t minBandingSize; //Any matrix bigger than this number squared will be broken apart with banding
    int32_t minBandingConstraintDistance; //The minimum size of a dp matrix between banding constraints.
    int32_t minTraceBackDiag; //The x+y diagonal to leave between the cut point and the place we choose new cutpoints.
    int32_t minTraceGapDiags; //The distance to leave between a cutpoint and the traceback
    int32_t constraintDiagonalTrim; //Amount to remove from a diagonal to be considered for a banding constraint
    bool alignAmbiguityCharacters; //Align sequences that are not 'actgACTG' as ambiguity characters, else force them to be unaligned.
} PairwiseAlignmentBandingParameters;

PairwiseAlignmentBandingParameters *pairwiseAlignmentBandingParameters_construct();

void pairwiseAlignmentBandingParameters_destruct(PairwiseAlignmentBandingParameters *p);

/*
 * A heuristic, banded form of getAlignedPairs.
 */
stList *getAlignedPairs_Fast(const char *sX, const char *sY, PairwiseAlignmentBandingParameters *p);

#endif /* PAIRWISEALIGNER_H_ */
