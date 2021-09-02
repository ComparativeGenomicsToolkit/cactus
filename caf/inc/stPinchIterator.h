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

typedef struct _stPinchIterator {
    int64_t alignmentTrim;
    void *alignmentArg;
    stPinch *(*getNextAlignment)(void *, stPinch *);
    void *(*startAlignmentStack)(void *);
    void (*destructAlignmentArg)(void *);
} stPinchIterator;

/*
 * Get next alignment from iterator. pinchToFillOut is filled out and returned. A NULL return value indicates
 * there are no further pinches
 */
stPinch *stPinchIterator_getNext(stPinchIterator *stPinchIterator, stPinch *pinchToFillOut);

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
stPinchIterator *stPinchIterator_constructFromFile(const char *alignmentFile);

/*
 * Constructs iterator from aligned pairs.
 */
stPinchIterator *stPinchIterator_constructFromAlignedPairs(stSortedSet *alignedPairs,
                                                           stPinch *(*getNextAlignedPairAlignment)(stSortedSetIterator *, stPinch *));

/*
 * Sets the amount to trim from the ends of each pinch in bases.
 */
void stPinchIterator_setTrim(stPinchIterator *pinchIterator, int64_t alignmentTrim);

#endif /* ST_PINCH_ITERATOR_H_ */
