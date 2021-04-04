/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Sequence *sequence_construct2(Name name, int64_t start,
		int64_t length, Name stringName, const char *header,
		Event *event, bool isTrivialSequence, CactusDisk *cactusDisk) {
	Sequence *sequence;

	sequence = st_malloc(sizeof(Sequence));
	sequence->name = name;
	assert(length >= 0);
	sequence->start = start;
	sequence->length = length;
	sequence->stringName = stringName;
	sequence->event = event;
	sequence->cactusDisk = cactusDisk;
	sequence->header = stString_copy(header != NULL ? header : "");
	sequence->isTrivialSequence = isTrivialSequence;

	cactusDisk_addSequence(cactusDisk, sequence);
	return sequence;
}

Sequence *sequence_construct3(int64_t start, int64_t length,
        const char *string, const char *header, Event *event,
        bool isTrivialSequence, CactusDisk *cactusDisk) {
    Name name;
    assert(strlen(string) == length);
    name = cactusDisk_addString(cactusDisk, string);
    return sequence_construct2(cactusDisk_getUniqueID(cactusDisk), start, length,
            name, header, event, isTrivialSequence, cactusDisk);
}

Sequence *sequence_construct(int64_t start, int64_t length,
		const char *string, const char *header, Event *event, CactusDisk *cactusDisk) {
	return sequence_construct3(start, length, string, header, event, 0, cactusDisk);
}

void sequence_destruct(Sequence *sequence) {
	cactusDisk_removeSequence(sequence->cactusDisk, sequence);
	free(sequence->header);
	free(sequence);
}

Name sequence_getName(Sequence *sequence) {
	return sequence->name;
}

int64_t sequence_getStart(Sequence *sequence) {
	return sequence->start;
}

int64_t sequence_getLength(Sequence *sequence) {
	return sequence->length;
}

Event *sequence_getEvent(Sequence *sequence) {
	return sequence->event;
}

char *sequence_getString(Sequence *sequence, int64_t start, int64_t length, int64_t strand) {
	assert(start >= sequence_getStart(sequence));
	assert(length >= 0);
	assert(start + length <= sequence_getStart(sequence) + sequence_getLength(sequence));
	return cactusDisk_getString(sequence->cactusDisk, sequence->stringName, start - sequence_getStart(sequence), length, strand, sequence->length);
}

const char *sequence_getHeader(Sequence *sequence) {
	return sequence->header;
}


bool sequence_isTrivialSequence(Sequence *sequence) {
    return sequence->isTrivialSequence;
}

void sequence_setHeader(Sequence *sequence,
                            char *newHeader) {
	free(sequence->header);
	sequence->header = newHeader;
}
