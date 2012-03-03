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

typedef struct _pairwiseAlignmentBandingParameters {
    int32_t threshold; //Minimum posterior probability of a match to be added to the output
    int32_t minDiagsBetweenTraceBack; //Minimum x+y diagonals to leave between doing traceback.
    int32_t traceBackDiagonals; //Number of diagonals to leave between trace back diagonal
    int32_t diagonalExpansion; //The number of x-y diagonals to expand around an anchor point
    int32_t constraintDiagonalTrim; //Amount to remove from a diagonal to be considered for a banding constraint
    int64_t anchorMatrixBiggerThanThis; //Search for anchors on any matrix bigger than this
    int64_t repeatMaskMatrixBiggerThanThis; //Any matrix in the anchors bigger than this is searched for anchors using non-repeat masked sequences.
    int64_t splitMatrixBiggerThanThis; //Any matrix in the anchors bigger than this is split into two.
    bool alignAmbiguityCharacters;
} PairwiseAlignmentParameters;

PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters_construct();

void pairwiseAlignmentBandingParameters_destruct(PairwiseAlignmentParameters *p);

/*
 * Gets the set of posterior match probabilities under a simple HMM model of alignment for two DNA sequences.
 */
stList *getAlignedPairs(const char *string1, const char *string2, PairwiseAlignmentParameters *p);


#endif /* PAIRWISEALIGNER_H_ */
