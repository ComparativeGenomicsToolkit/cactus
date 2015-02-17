/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_END_PRIVATE_H_
#define CACTUS_END_PRIVATE_H_

#include "cactusGlobals.h"

typedef struct _endContents {
	Cap *rootInstance;
	bool isStub;
	bool isAttached;
	Name name;
	Block *attachedBlock;
	stSortedSet *caps;
	Group *group;
	Flower *flower;
} EndContents;

struct _end_instanceIterator {
	stSortedSetIterator *iterator;
	End *end;
};

struct _end {
	EndContents *endContents;
	bool orientation;
	bool side;
	End *rEnd;
};


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//End functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs the end, but not any attached block.
 */
End *end_construct3(Name name, int64_t isStub, int64_t isAttached, int64_t side, Flower *flower);

/*
 * Destructs the end and any contained caps.
 */
void end_destruct(End *end);

/*
 * Sets the attached block.
 */
void end_setBlock(End *end, Block *block);

/*
 * Adds the cap to the end.
 */
void end_addInstance(End *end, Cap *cap);

/*
 * Removes the instance from the end.
 */
void end_removeInstance(End *end, Cap *cap);

/*
 * Write a binary representation of the end to the write function.
 */
void end_writeBinaryRepresentation(End *end, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a flower into memory from a binary representation of the flower.
 */
End *end_loadFromBinaryRepresentation(void **binaryString, Flower *flower);

/*
 * Hash key for an end, uses the name of the end to hash.. hence
 * the key doesn't care about the orientation.
 */
uint64_t end_hashKey(const void *o);

/*
 * Hash equals key, equal only if the two ends have the same name and orientation.
 */
int end_hashEqualsKey(const void *o, const void *o2);

/*
 * Sets the flower associated with the end.
 */
void end_setFlower(End *end, Flower *flower);


#endif
