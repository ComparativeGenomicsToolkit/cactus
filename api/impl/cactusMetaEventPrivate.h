#ifndef CACTUS_META_EVENT_PRIVATE_H_
#define CACTUS_META_EVENT_PRIVATE_H_

#include "cactusGlobals.h"

struct _metaEvent {
	Name name;
	char *header;
	CactusDisk *cactusDisk;
};


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta event private functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a meta event, which contains all the essential info for a event.
 */
MetaEvent *metaEvent_construct2(Name name, const char *header,
		CactusDisk *cactusDisk);

/*
 * Destructs a meta event.
 */
void metaEvent_destruct(MetaEvent *metaEvent);

/*
 * Creates a binary representation of the eventTree, returned as a char string.
 */
void metaEvent_writeBinaryRepresentation(MetaEvent *metaEvent, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a eventTree into memory from a binary representation of the eventTree.
 */
MetaEvent *metaEvent_loadFromBinaryRepresentation(void **binaryString, CactusDisk *cactusDisk);

#endif
