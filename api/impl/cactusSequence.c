/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Sequence *sequence_construct(MetaSequence *metaSequence, Flower *flower) {
    flower_addSequence(flower, metaSequence);
    return metaSequence;
	/*Sequence *sequence;
	sequence = st_malloc(sizeof(Sequence));
	sequence->metaSequence = metaSequence;
	sequence->flower = flower;
	flower_addSequence(flower, sequence);
	return sequence;*/
}

void sequence_destruct(Sequence *sequence) {
    assert(0);
	//flower_removeSequence(sequence_getFlower(sequence), sequence);
	//free(sequence);
}

MetaSequence *sequence_getMetaSequence(Sequence *sequence) {
    return sequence;
	//return sequence->metaSequence;
}

int64_t sequence_getStart(Sequence *sequence) {
	return metaSequence_getStart(sequence);
}

int64_t sequence_getLength(Sequence *sequence) {
	return metaSequence_getLength(sequence);
}

Name sequence_getName(Sequence *sequence) {
	return metaSequence_getName(sequence);
}

Event *sequence_getEvent(Sequence *sequence) {
	return metaSequence_getEvent(sequence); //NULL; //eventTree_getEvent(flower_getEventTree(sequence_getFlower(sequence)), metaSequence_getEventName(sequence->metaSequence));
}

/*Flower *sequence_getFlower(Sequence *sequence) {
    assert(0);
	return NULL; //sequence->flower;
}*/

char *sequence_getString(Sequence *sequence, int64_t start, int64_t length, bool strand) {
	return metaSequence_getString(sequence, start, length, strand);
}

const char *sequence_getHeader(Sequence *sequence) {
	return metaSequence_getHeader(sequence);
}

void sequence_check(Sequence *sequence) {
	/*Flower *flower = sequence_getFlower(sequence);
	cactusCheck(flower_getSequence(flower, sequence_getName(sequence)) == sequence); //properly connected to the flower..

	Group *parentGroup = flower_getParentGroup(flower);
	if(parentGroup != NULL) {
		Flower *parentFlower = group_getFlower(parentGroup);
		Sequence *parentSequence = flower_getSequence(parentFlower, sequence_getName(sequence));
		if(parentSequence != NULL) {
			cactusCheck(event_getName(sequence_getEvent(sequence)) == event_getName(sequence_getEvent(parentSequence)));
		}
	}*/
}

/*
 * Serialisation functions.
 */

void sequence_writeBinaryRepresentation(Sequence *sequence, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_SEQUENCE, writeFn);
	binaryRepresentation_writeName(sequence_getName(sequence), writeFn);
}

Sequence *sequence_loadFromBinaryRepresentation(void **binaryString, Flower *flower) {
	Sequence *sequence;

	sequence = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_SEQUENCE) {
		binaryRepresentation_popNextElementType(binaryString);
		sequence = sequence_construct(cactusDisk_getMetaSequence(flower_getCactusDisk(flower), binaryRepresentation_getName(binaryString)), flower);
	}
	return sequence;
}
