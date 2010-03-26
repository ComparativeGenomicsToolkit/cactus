#ifndef CACTUS_ATOM_PRIVATE_H_
#define CACTUS_ATOM_PRIVATE_H_

#include "cactusGlobals.h"

typedef struct _blockContents {
	Name name;
	struct avl_table *segments;
	int32_t length;
	Net *net;
} BlockContents;

struct _block_instanceIterator {
	struct avl_traverser *iterator;
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
Block *block_construct2(Name name, int32_t length, End *leftEnd, End *rightEnd, Net *net);

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
 * Loads a net into memory from a binary representation of the net.
 */
Block *block_loadFromBinaryRepresentation(void **binaryString, Net *net);

/*
 * Get a static instance (from the heap) with the name set.
 */
Block *block_getStaticNameWrapper(Name name);

#endif
