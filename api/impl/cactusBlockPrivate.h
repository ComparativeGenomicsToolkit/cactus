/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_ATOM_PRIVATE_H_
#define CACTUS_ATOM_PRIVATE_H_

#include "cactusGlobals.h"

struct _block_instanceIterator {
	Segment *segment;
	Block *block;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Block functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs the block, but not its ends.
 */
Block *block_construct2(Name name, int64_t length, End *leftEnd, End *rightEnd, Flower *flower);

/*
 * Destructs the block and all segments it contains.
 */
void block_destruct(Block *block);

/*
 * Adds in the instance to the block.
 */
void block_addInstance(Block *block, Segment *segment);

/*
 * Removes the instance from the block.
 */
void block_removeInstance(Block *block, Segment *segment);

#endif
