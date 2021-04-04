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
 * Gets an iterator to iterate over the segments.
 */
Block_InstanceIterator *block_getInstanceIterator(Block *block);

/*
 * Gets the next segment, or NULL if none remaining.
 */
Segment *block_getNext(Block_InstanceIterator *iterator);

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
 * Checks (amongst other things) the following:
 * Checks the reverse is the mirror of the block.
 * Checks the two ends are block ends and there properties are consistent with the block.
 * That the block has non-zero length.
 * For each segment calls segment_check.
 */
void block_check(Block *block);

/*
 * Returns non-zero if the ends of the block are not in a link.
 */
bool block_isTrivialChain(Block *block);

#endif
