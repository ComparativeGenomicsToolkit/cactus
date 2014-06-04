/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_META_SEQUENCE_PRIVATE_H_
#define CACTUS_META_SEQUENCE_PRIVATE_H_

#include "cactusGlobals.h"

struct _metaSequence {
	Name name;
	Name stringName;
	int64_t start;
	int64_t length;
	Name eventName;
	CactusDisk *cactusDisk;
	char *header;
	bool isTrivialSequence; //This flag is used to indicate if a sequence is trivial.
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private meta sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a meta sequence using an existing reference to a sequence in the sequence file.
 */
MetaSequence *metaSequence_construct2(Name name, int64_t start, int64_t length, Name stringName, const char *header,
		Name eventName, bool isTrivialSequence, CactusDisk *cactusDisk);

/*
 * Destructs a meta sequence.
 */
void metaSequence_destruct(MetaSequence *metaSequence);

/*
 * Creates a binary representation of the eventTree, returned as a char string.
 */
void metaSequence_writeBinaryRepresentation(MetaSequence *metaSequence, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a eventTree into memory from a binary representation of the eventTree.
 */
MetaSequence *metaSequence_loadFromBinaryRepresentation(void **binaryString, CactusDisk *cactusDisk);

#endif
