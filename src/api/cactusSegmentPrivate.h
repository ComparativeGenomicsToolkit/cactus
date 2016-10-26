/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_ATOM_INSTANCE_PRIVATE_H_
#define CACTUS_ATOM_INSTANCE_PRIVATE_H_

#include "cactusGlobals.h"

struct _segment {
	Cap *_5Cap;
	Segment *rInstance;
	Name name;
	Block *block;
};


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private segment functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs segment with the two caps, which must both have the same instance name.
 */
Segment *segment_construct3(Name name, Block *block,
		Cap *_5Cap, Cap *_3Cap);

/*
 * Destruct the segment, does not destruct ends.
 */
void segment_destruct(Segment *segment);

/*
 * Write a binary representation of the segment to the write function.
 */
void segment_writeBinaryRepresentation(Segment *segment, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a flower into memory from a binary representation of the flower.
 */
Segment *segment_loadFromBinaryRepresentation(void **binaryString, Block *block);

#endif
