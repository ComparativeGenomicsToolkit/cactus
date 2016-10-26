/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

MetaSequence *metaSequence_construct2(Name name, int64_t start,
		int64_t length, Name stringName, const char *header,
		Name eventName, bool isTrivialSequence, CactusDisk *cactusDisk) {
	MetaSequence *metaSequence;

	metaSequence = st_malloc(sizeof(MetaSequence));
	metaSequence->name = name;
	assert(length >= 0);
	metaSequence->start = start;
	metaSequence->length = length;
	metaSequence->stringName = stringName;
	metaSequence->eventName = eventName;
	metaSequence->cactusDisk = cactusDisk;
	metaSequence->header = stString_copy(header != NULL ? header : "");
	metaSequence->isTrivialSequence = isTrivialSequence;

	cactusDisk_addMetaSequence(cactusDisk, metaSequence);
	return metaSequence;
}

MetaSequence *metaSequence_construct3(int64_t start, int64_t length,
        const char *string, const char *header, Name eventName,
        bool isTrivialSequence, CactusDisk *cactusDisk) {
    Name name;
    assert(strlen(string) == length);
    name = cactusDisk_addString(cactusDisk, string);
    return metaSequence_construct2(cactusDisk_getUniqueID(cactusDisk), start, length,
            name, header, eventName, isTrivialSequence, cactusDisk);
}

MetaSequence *metaSequence_construct(int64_t start, int64_t length,
		const char *string, const char *header, Name eventName, CactusDisk *cactusDisk) {
	return metaSequence_construct3(start, length, string, header, eventName, 0, cactusDisk);
}

void metaSequence_destruct(MetaSequence *metaSequence) {
	cactusDisk_removeMetaSequence(metaSequence->cactusDisk, metaSequence);
	free(metaSequence->header);
	free(metaSequence);
}

Name metaSequence_getName(MetaSequence *metaSequence) {
	return metaSequence->name;
}

int64_t metaSequence_getStart(MetaSequence *metaSequence) {
	return metaSequence->start;
}

int64_t metaSequence_getLength(MetaSequence *metaSequence) {
	return metaSequence->length;
}

Name metaSequence_getEventName(MetaSequence *metaSequence) {
	return metaSequence->eventName;
}

char *metaSequence_getString(MetaSequence *metaSequence, int64_t start, int64_t length, int64_t strand) {
	assert(start >= metaSequence_getStart(metaSequence));
	assert(length >= 0);
	assert(start + length <= metaSequence_getStart(metaSequence) + metaSequence_getLength(metaSequence));
	return cactusDisk_getString(metaSequence->cactusDisk, metaSequence->stringName, start - metaSequence_getStart(metaSequence), length, strand, metaSequence->length);
}

const char *metaSequence_getHeader(MetaSequence *metaSequence) {
	return metaSequence->header;
}


bool metaSequence_isTrivialSequence(MetaSequence *metaSequence) {
    return metaSequence->isTrivialSequence;
}

void metaSequence_setHeader(MetaSequence *metaSequence,
                            char *newHeader) {
	free(metaSequence->header);
	metaSequence->header = newHeader;
}

/*
 * Serialisation functions.
 */

void metaSequence_writeBinaryRepresentation(MetaSequence *metaSequence,
		void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_META_SEQUENCE, writeFn);
	binaryRepresentation_writeName(metaSequence_getName(metaSequence), writeFn);
	binaryRepresentation_writeInteger(metaSequence_getStart(metaSequence), writeFn);
	binaryRepresentation_writeInteger(metaSequence_getLength(metaSequence), writeFn);
	binaryRepresentation_writeName(metaSequence_getEventName(metaSequence), writeFn);
	binaryRepresentation_writeName(metaSequence->stringName, writeFn);
	binaryRepresentation_writeString(metaSequence_getHeader(metaSequence), writeFn);
	binaryRepresentation_writeBool(metaSequence_isTrivialSequence(metaSequence), writeFn);
}

MetaSequence *metaSequence_loadFromBinaryRepresentation(void **binaryString,
		CactusDisk *cactusDisk) {
	MetaSequence *metaSequence;
	Name name;
	int64_t start;
	int64_t length;
	Name stringName;
	Name eventName;
	char *header;

	metaSequence = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_META_SEQUENCE) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		start = binaryRepresentation_getInteger(binaryString);
		length = binaryRepresentation_getInteger(binaryString);
		eventName = binaryRepresentation_getName(binaryString);
		stringName = binaryRepresentation_getName(binaryString);
		header = binaryRepresentation_getString(binaryString);
		bool isTrivialSequence = binaryRepresentation_getBool(binaryString);
		metaSequence = metaSequence_construct2(name, start, length,
				stringName, header, eventName, isTrivialSequence, cactusDisk);
		free(header);
	}
	return metaSequence;
}

