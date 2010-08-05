#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

MetaEvent *metaEvent_construct(const char *header,
		CactusDisk *cactusDisk) {
	return metaEvent_construct2(cactusDisk_getUniqueID(cactusDisk), header, cactusDisk);
}

MetaEvent *metaEvent_construct2(Name name, const char *header,
		CactusDisk *cactusDisk) {
	MetaEvent *metaEvent;
	metaEvent = st_malloc(sizeof(MetaEvent));
	metaEvent->name = name;
	metaEvent->header = stString_copy(header != NULL ? header : "");
	metaEvent->cactusDisk = cactusDisk;
	cactusDisk_addMetaEvent(cactusDisk, metaEvent);
	return metaEvent;
}

Name metaEvent_getName(MetaEvent *metaEvent) {
	return metaEvent->name;
}

const char *metaEvent_getHeader(MetaEvent *metaEvent) {
	return metaEvent->header;
}

void metaEvent_destruct(MetaEvent *metaEvent) {
	cactusDisk_unloadMetaEvent(metaEvent->cactusDisk, metaEvent);
	free(metaEvent->header);
	free(metaEvent);
}

/*
 * Serialisation functions.
 */

void metaEvent_writeBinaryRepresentation(MetaEvent *metaEvent, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_META_EVENT, writeFn);
	binaryRepresentation_writeName(metaEvent_getName(metaEvent), writeFn);
	binaryRepresentation_writeString(metaEvent_getHeader(metaEvent), writeFn);
}

MetaEvent *metaEvent_loadFromBinaryRepresentation(void **binaryString, CactusDisk *cactusDisk) {
	MetaEvent *metaEvent;
	Name name;
	char *header;

	metaEvent = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_META_EVENT) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		header = binaryRepresentation_getString(binaryString);
		metaEvent = metaEvent_construct2(name, header, cactusDisk);
		free(header);
	}
	return metaEvent;
}
