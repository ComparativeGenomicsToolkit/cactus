/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_ATOM_H_
#define CACTUS_ATOM_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic block functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs the block, but not its ends.
 */
Block *block_construct(int64_t length, Flower *flower);

/*
 * Returns string name of the block.
 */
Name block_getName(Block *block);

/*
 * Returns a non zero if the block is oriented positively.
 * The orientation is arbitrary (it is not explicitly with respect to anything else), but is consistent.
 */
bool block_getOrientation(Block *block);

/*
 * Gets a positively oriented block .
 */
Block *block_getPositiveOrientation(Block *block);

/*
 * Returns a reversed view of the block (in the opposite orientation).
 */
Block *block_getReverse(Block *block);

/*
 * Returns the length in bases of the block.
 */
int64_t block_getLength(Block *block);

/*
 * Gets the flower the block is part of.
 */
Flower *block_getFlower(Block *block);

/*
 * Gets the 5 end of the block.
 */
End *block_get5End(Block *block);

/*
 * Gets the 3 end of the block.
 */
End *block_get3End(Block *block);

/*
 * Returns the number of instances (including any internal instances), the block contains.
 */
int64_t block_getInstanceNumber(Block *block);

/*
 * Gets the segment using its instance name as a key. Instance name is m of full name n.m.
 */
Segment *block_getInstance(Block *block, Name instanceName);

/*
 * Gets the first segment in the list.
 */
Segment *block_getFirst(Block *block);

/*
 * Gets the root cap of the block, if it is set, or returns NULL;
 */
Segment *block_getRootInstance(Block *block);

/*
 * Sets the root cap of the block. Will throw an error if the segment
 * is not part of the end, or already has a parent. Will set the root instances of left and
 * right ends.
 */
void block_setRootInstance(Block *block, Segment *segment);

/*
 * Gets an iterator to iterate over the segments.
 */
Block_InstanceIterator *block_getInstanceIterator(Block *block);

/*
 * Gets the next segment, or NULL if none remaining.
 */
Segment *block_getNext(Block_InstanceIterator *iterator);

/*
 * Gets the previous segment, or NULL if none remaining.
 */
Segment *block_getPrevious(Block_InstanceIterator *iterator);

/*
 * Duplicates the iterator.
 */
Block_InstanceIterator *block_copyInstanceIterator(Block_InstanceIterator *iterator);

/*
 * Destructs the iterator - should always be coupled with the iterator.
 */
void block_destructInstanceIterator(Block_InstanceIterator *block);

/*
 * Gets any chain associated with the block.
 */
Chain *block_getChain(Block *block);

/*
 * Get an arbitrary segment with whose event has the given name, else NULL.
 */
Segment *block_getSegmentForEvent(Block *block, Name eventName);

/*
 * Splits an block into two. The split point is equal to the length of the left block. This value must be less than the length
 * of the complete block and greater than zero, so that both sides of the block have greater than zero length.
 *
 * The original block and all its instances are destroyed by this operation. The new blocks are pointed to by the left and right block pointers.
 */
void block_split(Block *block, int64_t splitPoint, Block **leftBlock, Block **rightBlock);

/*
 * Checks (amongst other things) the following:
 * Checks the reverse is the mirror of the block.
 * Checks the two ends are block ends and there properties are consistent with the block.
 * That the block has non-zero length.
 * For each segment calls segment_check.
 */
void block_check(Block *block);

/*
 * Makes a newick string representation of the segments in the block.
 * Does not include any branch lengths currents.
 * Creates an assert error if the segments are not in a tree.
 * Returns null if the block contains no segments.
 * The names of the blocks are their names converted to strings.
 * Includes unary events in tree.
 */
char *block_makeNewickString(Block *block, bool includeInternalNames, bool includeUnaryEvents);

/*
 * Pushes the block into the higher level flower. Will not work if the block to be promoted is nested
 * within a chain of the higher level flower or is contained in a chain at the lower level (i.e. it must
 * be a trivial chain).
 */
void block_promote(Block *block);

/*
 * Returns non-zero if the ends of the block are not in a link.
 */
bool block_isTrivialChain(Block *block);

#endif
