#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

MetaEvent *metaEvent_construct(const char *header,
		NetDisk *netDisk) {
	return metaEvent_construct2(netDisk_getUniqueID(netDisk), header, netDisk);
}

MetaEvent *metaEvent_construct2(Name name, const char *header,
		NetDisk *netDisk) {
	MetaEvent *metaEvent;
	metaEvent = mallocLocal(sizeof(MetaEvent));
	metaEvent->name = name;
	metaEvent->header = stringCopy(header != NULL ? header : "");
	metaEvent->netDisk = netDisk;
	netDisk_addMetaEvent(netDisk, metaEvent);
	return metaEvent;
}

Name metaEvent_getName(MetaEvent *metaEvent) {
	return metaEvent->name;
}

const char *metaEvent_getHeader(MetaEvent *metaEvent) {
	return metaEvent->header;
}

void metaEvent_destruct(MetaEvent *metaEvent) {
	netDisk_unloadMetaEvent(metaEvent->netDisk, metaEvent);
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

MetaEvent *metaEvent_loadFromBinaryRepresentation(void **binaryString, NetDisk *netDisk) {
	MetaEvent *metaEvent;
	Name name;
	char *header;

	metaEvent = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_META_EVENT) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		header = binaryRepresentation_getString(binaryString);
		metaEvent = metaEvent_construct2(name, header, netDisk);
		free(header);
	}
	return metaEvent;
}
