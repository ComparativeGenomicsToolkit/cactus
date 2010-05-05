#ifndef CACTUS_META_SEQUENCE_PRIVATE_H_
#define CACTUS_META_SEQUENCE_PRIVATE_H_

#include "cactusGlobals.h"

struct _metaSequence {
	Name name;
	int64_t fileOffset;
	int32_t start;
	int32_t length;
	Name eventName;
	NetDisk *netDisk;
	char *header;
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
MetaSequence *metaSequence_construct2(Name name, int32_t start, int32_t length, int64_t fileOffset, const char *header,
		Name eventName, NetDisk *netDisk);

/*
 * Destructs a meta sequence.
 */
void metaSequence_destruct(MetaSequence *metaSequence);

/*
 * Gets the file offset location of the string backing the metasequence.
 */
int64_t metaSequence_getFileOffset(MetaSequence *metaSequence);

/*
 * Creates a binary representation of the eventTree, returned as a char string.
 */
void metaSequence_writeBinaryRepresentation(MetaSequence *metaSequence, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a eventTree into memory from a binary representation of the eventTree.
 */
MetaSequence *metaSequence_loadFromBinaryRepresentation(void **binaryString, NetDisk *netDisk);

#endif
