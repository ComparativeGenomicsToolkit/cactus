/*
 * chaining.h
 *
 *  Created on: 8 Mar 2012
 *      Author: benedictpaten
 */

#ifndef CHAINING_H_
#define CHAINING_H_

typedef struct _APair {
    int32_t x;
    int32_t y;
    int32_t fScore;
    int32_t bScore;
} APair;

typedef struct _APairArray {
    APair *aPairs;
    int32_t length;
} APairArray;

APairArray *aPairArray_construct(stList *blastPairs);

void aPairArray_destruct(APairArray *aPairArray);

void aPairArray_calculateForwardScores(APairArray *aPairArray);

void aPairArray_calculateBackwardScores(APairArray *aPairArray);

stList *filterToRemoveOverlap(stList *overlappingPairs);

/*
 * Gets back a sorted, monotonically increasing set of pairs.
 */
stList *getAnchorChain(stList *sortedBlastPairs, double alpha);

#endif /* CHAINING_H_ */
