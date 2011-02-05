/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic pseudo-chromosome functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

PseudoChromosome *pseudoChromosome_construct(Reference *reference,
		End *_5End, End *_3End) {
	return pseudoChromosome_construct2(cactusDisk_getUniqueID(flower_getCactusDisk(reference_getFlower(reference))),
			reference, _5End, _3End);
}

Name pseudoChromosome_getName(PseudoChromosome *pseudoChromosome) {
	return pseudoChromosome->name;
}

End *pseudoChromosome_get5End(PseudoChromosome *pseudoChromosome) {
	return pseudoChromosome->_5End;
}

End *pseudoChromosome_get3End(PseudoChromosome *pseudoChromosome) {
	return pseudoChromosome->_3End;
}

Reference *pseudoChromosome_getReference(PseudoChromosome *pseudoChromosome) {
	return pseudoChromosome->reference;
}

int32_t pseudoChromosome_getPseudoAdjacencyNumber(PseudoChromosome *pseudoChromosome) {
	return stList_length(pseudoChromosome->pseudoAdjacencies);
}

PseudoAdjacency *pseudoChromosome_getPseudoAdjacencyByIndex(PseudoChromosome *pseudoChromosome, int32_t i) {
    return stList_get(pseudoChromosome->pseudoAdjacencies, i);
}

PseudoChromsome_PseudoAdjacencyIterator *pseudoChromosome_getPseudoAdjacencyIterator(PseudoChromosome *pseudoChromosome) {
	return stList_getIterator(pseudoChromosome->pseudoAdjacencies);
}

PseudoAdjacency *pseudoChromosome_getNextPseudoAdjacency(PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator) {
	return stList_getNext(pseudoAdjacencyIterator);
}

PseudoAdjacency *pseudoChromosome_getPreviousPseudoAdjacency(PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator) {
	return stList_getPrevious(pseudoAdjacencyIterator);
}

PseudoChromsome_PseudoAdjacencyIterator *pseudoChromosome_copyPseudoChromosomeIterator(PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator) {
	return stList_copyIterator(pseudoAdjacencyIterator);
}

void pseudoChromosome_destructPseudoAdjacencyIterator(PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator) {
	stList_destructIterator(pseudoAdjacencyIterator);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

PseudoChromosome *pseudoChromosome_construct2(Name name, Reference *reference,
		End *_5End, End *_3End) {
	PseudoChromosome *pseudoChromosome = st_malloc(sizeof(PseudoChromosome));
	assert(name != NULL_NAME);
	assert(reference != NULL);
	assert(_5End != NULL);
	assert(_3End != NULL);

	pseudoChromosome->pseudoAdjacencies = stList_construct3(0, (void (*)(void *))pseudoAdjacency_destruct);
	pseudoChromosome->_5End = end_getPositiveOrientation(_5End); //everything is on the positive orientation.
	pseudoChromosome->_3End = end_getPositiveOrientation(_3End);
	pseudoChromosome->reference = reference;
	pseudoChromosome->name = name;
	reference_addPseudoChromosome(reference, pseudoChromosome);
	return pseudoChromosome;
}

void pseudoChromosome_destruct(PseudoChromosome *pseudoChromosome) {
	reference_removePseudoChromosome(pseudoChromosome_getReference(pseudoChromosome), pseudoChromosome);
	stList_destruct(pseudoChromosome->pseudoAdjacencies);
	free(pseudoChromosome);
}

void pseudoChromosome_addPseudoAdjacency(PseudoChromosome *pseudoChromosome, PseudoAdjacency *pseudoAdjacency) {
	assert(pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome) == pseudoAdjacency_getIndex(pseudoAdjacency));
	stList_append(pseudoChromosome->pseudoAdjacencies, pseudoAdjacency);
}

static PseudoChromosome pseudoChromosome_getStaticNameWrapperP;
PseudoChromosome *pseudoChromosome_getStaticNameWrapper(Name name) {
	pseudoChromosome_getStaticNameWrapperP.name = name;
	return &pseudoChromosome_getStaticNameWrapperP;
}

void pseudoChromosome_writeBinaryRepresentation(PseudoChromosome *pseudoChromosome, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	PseudoChromsome_PseudoAdjacencyIterator *iterator;
	PseudoAdjacency *pseudoAdjacency;
	binaryRepresentation_writeElementType(CODE_PSEUDO_CHROMOSOME, writeFn);
	binaryRepresentation_writeName(pseudoChromosome_getName(pseudoChromosome), writeFn);
	binaryRepresentation_writeName(end_getName(pseudoChromosome_get5End(pseudoChromosome)), writeFn);
	binaryRepresentation_writeName(end_getName(pseudoChromosome_get3End(pseudoChromosome)), writeFn);
	iterator = pseudoChromosome_getPseudoAdjacencyIterator(pseudoChromosome);
	while((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(iterator)) != NULL) {
		pseudoAdjacency_writeBinaryRepresentation(pseudoAdjacency, writeFn);
	}
	pseudoChromosome_destructPseudoAdjacencyIterator(iterator);
	binaryRepresentation_writeElementType(CODE_PSEUDO_CHROMOSOME, writeFn);
}

PseudoChromosome *pseudoChromosome_loadFromBinaryRepresentation(void **binaryString, Reference *reference) {
	PseudoChromosome *pseudoChromosome = NULL;
	End *_5End, *_3End;
	Name name;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_PSEUDO_CHROMOSOME) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		_5End = flower_getEnd(reference_getFlower(reference), binaryRepresentation_getName(binaryString));
		_3End = flower_getEnd(reference_getFlower(reference), binaryRepresentation_getName(binaryString));
		pseudoChromosome = pseudoChromosome_construct2(name, reference, _5End, _3End);
		while(pseudoAdjacency_loadFromBinaryRepresentation(binaryString, pseudoChromosome));
		assert(binaryRepresentation_peekNextElementType(*binaryString) == CODE_PSEUDO_CHROMOSOME);
		binaryRepresentation_popNextElementType(binaryString);
	}
	return pseudoChromosome;
}
