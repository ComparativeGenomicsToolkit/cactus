/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_ATOM_PRIVATE_H_
#define CACTUS_ATOM_PRIVATE_H_

#include "cactusGlobals.h"


typedef struct _blockContents {
	Name name;
	stSortedSet *segments;
	int64_t length;
	Flower *flower;
} BlockContents;

struct _block_instanceIterator {
	stSortedSetIterator *iterator;
	Block *block;
};

struct _block {
	BlockContents *blockContents;
	bool orientation;
	End *leftEnd;
	Block *rBlock;
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

/*
 * Write a binary representation of the block to the write function.
 */
void block_writeBinaryRepresentation(Block *block, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a flower into memory from a binary representation of the flower.
 */
Block *block_loadFromBinaryRepresentation(void **binaryString, Flower *flower);

/*
 * Sets the flower associated with the block.
 */
void block_setFlower(Block *block, Flower *flower);

#endif
