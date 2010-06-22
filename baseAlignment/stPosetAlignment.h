/*
 * stPosetAlignment.h
 *
 *  Created on: 21-Jun-2010
 *      Author: benedictpaten
 */

#include "sonLib.h"

#ifndef STPOSETALIGNMENT_H_
#define STPOSETALIGNMENT_H_

/*
 * The anonmyous typedef declaration of the alignment structure.
 */
typedef struct _stPosetAlignment stPosetAlignment;

/*
 * Constructs a poset alignment
 */
stPosetAlignment *stPosetAlignment_construct(void);

/*
 * Destructs the poset alignment.
 */
void stPosetAlignment_destruct(stPosetAlignment *posetAlignment);

/*
 * Adds a sequence to the poset alignment. Length >= 0. All sequences must be added
 * before you add any pairwise alignments to the alignment structure.
 */
void stPosetAlignment_addSequence(stPosetAlignment *posetAlignment, int32_t length);

/*
 * Returns the number of sequences in the alignment.
 */
int32_t stPosetAlignment_getSequenceNumber(stPosetAlignment *posetAlignment);

/*
 * Returns the length of the sequence for the given sequence.
 */
int32_t stPosetAlignment_getSequenceLength(stPosetAlignment *posetAlignment, int32_t sequence);

/*
 * Adds a pairwise alignment, can only be added once you've added all the sequences to the alignment. It will
 * create an exception if you try adding a self alignment to the structure.
 */
void stPosetAlignment_addPairwiseAlignment(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t sequence2);

/*
 * Returns the number of pairwise alignments the structure contains.
 */
int32_t stPosetAlignment_getNumberPairwiseAlignments(stPosetAlignment *posetAlignment);

/*
 * Gets a sorted set of the pairwise alignments that it contains.
 */
stSortedSet *stPosetAlignment_getPairwiseAlignments(stPosetAlignment *posetAlignment);

/*
 * Returns non-zero if the pair can be added to the alignment, but does not actually add the pair to the alignment.
 */
bool stPosetAlignment_isPossible(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2);

/*
 * List stPosetAlignment_isPossible, but if the return value is non-zero the pair is added to the alignment.
 */
bool stPosetAlignment_add(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2);

#endif /* STPOSETALIGNMENT_H_ */
