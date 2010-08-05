#ifndef CACTUS_SEQUENCE_PRIVATE_H_
#define CACTUS_SEQUENCE_PRIVATE_H_

#include "cactusGlobals.h"

struct _sequence {
	MetaSequence *metaSequence;
	Flower *flower;
};


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Write a binary representation of the sequence to the write function.
 */
void sequence_writeBinaryRepresentation(Sequence *sequence, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a sequence into memory from a binary representation of the sequence.
 */
Sequence *sequence_loadFromBinaryRepresentation(void **binaryString, Flower *flower);

/*
 * Get a static instance (from the heap) with the name set.
 */
Sequence *sequence_getStaticNameWrapper(Name name);

#endif
