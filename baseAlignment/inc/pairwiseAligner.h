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

/*
 * A heuristic, banded form of getAlignedPairs.
 */
stList *getAlignedPairs_Fast(const char *sX, const char *sY, int32_t bandingSize);

#endif /* PAIRWISEALIGNER_H_ */
