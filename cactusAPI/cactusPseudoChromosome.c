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
	return pseudoChromosome_construct2(netDisk_getUniqueID(net_getNetDisk(reference_getNet(reference))),
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
	return sortedSet_getLength(pseudoChromosome->pseudoAdjacencies);
}

PseudoAdjacency *pseudoChromosome_getPseudoAdjacency(PseudoChromosome *pseudoChromosome, Name name) {
	PseudoAdjacency *pseudoAdjacency;
	pseudoAdjacency = pseudoAdjacency_getStaticNameWrapper(name);
	return sortedSet_find(pseudoChromosome->pseudoAdjacencies, pseudoAdjacency);
}

PseudoAdjacency *pseudoChromosome_getFirst(PseudoChromosome *pseudoChromosome) {
	return sortedSet_getFirst(pseudoChromosome->pseudoAdjacencies);
}

PseudoAdjacency *pseudoChromosome_getLast(PseudoChromosome *pseudoChromosome) {
	return sortedSet_getLast(pseudoChromosome->pseudoAdjacencies);
}

PseudoChromsome_PseudoAdjacencyIterator *pseudoChromosome_getPseudoAdjacencyIterator(PseudoChromosome *pseudoChromosome) {
	return iterator_construct(pseudoChromosome->pseudoAdjacencies);
}

PseudoAdjacency *pseudoChromosome_getNextPseudoAdjacency(PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator) {
	return iterator_getNext(pseudoAdjacencyIterator);
}

PseudoAdjacency *pseudoChromosome_getPreviousPseudoAdjacency(PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator) {
	return iterator_getPrevious(pseudoAdjacencyIterator);
}

PseudoChromsome_PseudoAdjacencyIterator *pseudoChromosome_copyPseudoChromosomeIterator(PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator) {
	return iterator_copy(pseudoAdjacencyIterator);
}

void pseudoChromosome_destructPseudoAdjacencyIterator(PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator) {
	iterator_destruct(pseudoAdjacencyIterator);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static int32_t pseudoChromosome_constructP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(pseudoAdjacency_getName((PseudoAdjacency *)o1),
							   pseudoAdjacency_getName((PseudoAdjacency *)o2));
}

PseudoChromosome *pseudoChromosome_construct2(Name name, Reference *reference,
		End *_5End, End *_3End) {
	PseudoChromosome *pseudoChromosome = malloc(sizeof(PseudoChromosome));
	assert(name != NULL_NAME);
	assert(reference != NULL);
	assert(_5End != NULL);
	assert(_3End != NULL);


	pseudoChromosome->pseudoAdjacencies = sortedSet_construct(pseudoChromosome_constructP);
	pseudoChromosome->_5End = end_getPositiveOrientation(_5End); //everything is on the positive orientation.
	pseudoChromosome->_3End = end_getPositiveOrientation(_3End);
	pseudoChromosome->reference = reference;
	pseudoChromosome->name = name;
	reference_addPseudoChromosome(reference, pseudoChromosome);
	return pseudoChromosome;
}

void pseudoChromosome_destruct(PseudoChromosome *pseudoChromosome) {
	reference_removePseudoChromosome(pseudoChromosome_getReference(pseudoChromosome), pseudoChromosome);
	PseudoAdjacency *pseudoAdjacency;
	while((pseudoAdjacency = pseudoChromosome_getFirst(pseudoChromosome)) != NULL) {
		pseudoAdjacency_destruct(pseudoAdjacency);
	}
	free(pseudoChromosome);
}

void pseudoChromosome_addPseudoAdjacency(PseudoChromosome *pseudoChromosome, PseudoAdjacency *pseudoAdjacency) {
	assert(sortedSet_find(pseudoChromosome->pseudoAdjacencies, pseudoAdjacency) == NULL);
	sortedSet_insert(pseudoChromosome->pseudoAdjacencies, pseudoAdjacency);
}

void pseudoChromosome_removePseudoAdjacency(PseudoChromosome *pseudoChromosome, PseudoAdjacency *pseudoAdjacency) {
	assert(sortedSet_find(pseudoChromosome->pseudoAdjacencies, pseudoAdjacency) != NULL);
	sortedSet_delete(pseudoChromosome->pseudoAdjacencies, pseudoAdjacency);
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
	binaryRepresentation_writeInteger(pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome), writeFn);
	iterator = pseudoChromosome_getPseudoAdjacencyIterator(pseudoChromosome);
	while((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(iterator)) != NULL) {
		pseudoAdjacency_writeBinaryRepresentation(pseudoAdjacency, writeFn);
	}
	pseudoChromosome_destructPseudoAdjacencyIterator(iterator);
}

PseudoChromosome *pseudoChromosome_loadFromBinaryRepresentation(void **binaryString, Reference *reference) {
	PseudoChromosome *pseudoChromosome = NULL;
	int32_t pseudoAdjacencyNumber;
	End *_5End, *_3End;
	Name name;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_PSEUDO_CHROMOSOME) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		_5End = net_getEnd(reference_getNet(reference), binaryRepresentation_getName(binaryString));
		_3End = net_getEnd(reference_getNet(reference), binaryRepresentation_getName(binaryString));
		pseudoChromosome = pseudoChromosome_construct2(name, reference, _5End, _3End);
		pseudoAdjacencyNumber = binaryRepresentation_getInteger(binaryString);
		while(pseudoAdjacencyNumber-- > 0) {
			pseudoAdjacency_loadFromBinaryRepresentation(binaryString, pseudoChromosome);
		}
	}
	return pseudoChromosome;
}
