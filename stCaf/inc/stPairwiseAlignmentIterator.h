/*
 * pairwiseAlignmentIterator.h
 *
 *  Created on: 21 Mar 2012
 *      Author: benedictpaten
 */

#ifndef PAIRWISEALIGNMENTITERATOR_H_
#define PAIRWISEALIGNMENTITERATOR_H_

typedef struct _pairwiseAlignmentIterator PairwiseAlignmentIterator;

/*
 * Get next alignent from iterator.
 */
struct PairwiseAlignment *pairwiseAlignmentIterator_getNext(
        PairwiseAlignmentIterator *pairwiseAlignmentIterator);

/*
 * Reset the iterator, returning again to the beginning of the sequence.
 */
void pairwiseAlignmentIterator_reset(
        PairwiseAlignmentIterator *pairwiseAlignmentIterator);

/*
 * Cleanup the iterator
 */
void pairwiseAlignmentIterator_destruct(
        PairwiseAlignmentIterator *pairwiseAlignmentIterator);

/*
 * Cleanup an alignment produced by the sequence. DO NOT CALL destructPairwiseAlignment.
 */
void pairwiseAlignmentIterator_destructAlignment(PairwiseAlignmentIterator *pairwiseAlignmentIterator, struct PairwiseAlignment *pairwiseAlignment);

/*
 * Get a pairwise alignment iterator from a file.
 */
PairwiseAlignmentIterator *pairwiseAlignmentIterator_constructFromFile(
        const char *alignmentFile);

/*
 * Get a pairwise alignment iterator from a list of alignments.
 * Does not cleanup the list or modify the list.
 */
PairwiseAlignmentIterator *pairwiseAlignmentIterator_constructFromList(
        stList *alignmentsList);

/*
 * Constructs iterator from aligned pairs.
 */
PairwiseAlignmentIterator *pairwiseAlignmentIterator_constructFromAlignedPairs(
        stSortedSet *alignedPairs, struct PairwiseAlignment *(*getNextAlignedPairAlignment)(stSortedSetIterator *));

#endif /* PAIRWISEALIGNMENTITERATOR_H_ */
