/*
 * stPinchIterator.h
 *
 *  Created on: 21 Mar 2012
 *      Author: benedictpaten
 */

#ifndef ST_PINCH_ITERATOR_H_
#define ST_PINCH_ITERATOR_H_

#include "sonLib.h"
#include "stPinchGraphs.h"

typedef struct _stPinchIterator stPinchIterator;

/*
 * Get next alignment from iterator.
 */
stPinch *stPinchIterator_getNext(
        stPinchIterator *stPinchIterator);

/*
 * Reset the iterator, returning again to the beginning of the sequence.
 */
void stPinchIterator_reset(
        stPinchIterator *stPinchIterator);

/*
 * Cleanup the iterator
 */
void stPinchIterator_destruct(
        stPinchIterator *stPinchIterator);

/*
 * Get a pairwise alignment iterator from a file.
 */
stPinchIterator *stPinchIterator_constructFromFile(
        const char *alignmentFile);

/*
 * Get a pairwise alignment iterator from a list of alignments.
 * Does not cleanup the list or modify the list.
 */
stPinchIterator *stPinchIterator_constructFromList(
        stList *alignmentsList);

/*
 * Constructs iterator from aligned pairs.
 */
stPinchIterator *stPinchIterator_constructFromAlignedPairs(
        stSortedSet *alignedPairs, stPinch *(*getNextAlignedPairAlignment)(stSortedSetIterator *));

/*
 * Sets the amount to trim from the ends of each pinch in bases.
 */
void stPinchIterator_setTrim(stPinchIterator *pinchIterator, int64_t alignmentTrim);

#endif /* ST_PINCH_ITERATOR_H_ */
