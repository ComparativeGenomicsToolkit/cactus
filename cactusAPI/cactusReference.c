#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic reference functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Reference *reference_construct(Net *net) {
	return reference_construct2(netDisk_getUniqueID(net_getNetDisk(net)), net);
}

Name reference_getName(Reference *reference) {
	return reference->name;
}

Net *reference_getNet(Reference *reference) {
	return reference->net;
}

int32_t reference_getPseudoChromosomeNumber(Reference *reference) {
	return sortedSet_getLength(reference->pseudoChromosomes);
}

PseudoChromosome *reference_getPseudoChromosome(Reference *reference, Name name) {
	PseudoChromosome *pseudoChromosome;
	pseudoChromosome = pseudoChromosome_getStaticNameWrapper(name);
	return sortedSet_find(reference->pseudoChromosomes, pseudoChromosome);
}

PseudoChromosome *reference_getFirst(Reference *reference) {
	return sortedSet_getFirst(reference->pseudoChromosomes);
}

Reference_PseudoChromosomeIterator *reference_getPseudoChromosomeIterator(Reference *reference) {
	return iterator_construct(reference->pseudoChromosomes);
}

PseudoChromosome *reference_getNextPseudoChromosome(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
	return iterator_getNext(pseudoChromosomeIterator);
}

PseudoChromosome *reference_getPreviousPseudoChromosome(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
	return iterator_getPrevious(pseudoChromosomeIterator);
}

Reference_PseudoChromosomeIterator *reference_copyPseudoChromosomeIterator(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
	return iterator_copy(pseudoChromosomeIterator);
}

void reference_destructPseudoChromosomeIterator(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
	iterator_destruct(pseudoChromosomeIterator);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static int32_t reference_constructP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(pseudoChromosome_getName((PseudoChromosome *)o1),
			pseudoChromosome_getName((PseudoChromosome *)o2));
}

Reference *reference_construct2(Name name, Net *net) {
	Reference *reference = malloc(sizeof(Reference));
	//Setup the basic structure - a sorted set of pseudo-chromosomes.
	reference->pseudoChromosomes = sortedSet_construct(reference_constructP);
	//Link the reference and net.
	reference->net = net;
	net_addReference(net, reference);
	return reference;
}

void reference_destruct(Reference *reference) {
	PseudoChromosome *pseudoChromosome;
	net_removeReference(reference_getNet(reference), reference);
	while((pseudoChromosome = reference_getFirst(reference)) != NULL) {
			pseudoChromosome_destruct(pseudoChromosome);
	}
	sortedSet_destruct(reference->pseudoChromosomes, NULL);
	free(reference);
}

void reference_addPseudoChromosome(Reference *reference, PseudoChromosome *pseudoChromosome) {
	assert(sortedSet_find(reference->pseudoChromosomes, pseudoChromosome) == NULL);
	sortedSet_insert(reference->pseudoChromosomes, pseudoChromosome);
}

void reference_removePseudoChromosome(Reference *reference, PseudoChromosome *pseudoChromosome) {
	assert(sortedSet_find(reference->pseudoChromosomes, pseudoChromosome) != NULL);
	sortedSet_delete(reference->pseudoChromosomes, pseudoChromosome);
}

Reference reference_getStaticNameWrapperP;
Reference *reference_getStaticNameWrapper(Name name) {
	reference_getStaticNameWrapperP.name = name;
	return &reference_getStaticNameWrapperP;
}

void reference_writeBinaryRepresentation(Reference *reference, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	Reference_PseudoChromosomeIterator *iterator;
	PseudoChromosome *pseudoChromosome;
	binaryRepresentation_writeElementType(CODE_REFERENCE, writeFn);
	binaryRepresentation_writeName(reference_getName(reference), writeFn);
	binaryRepresentation_writeInteger(reference_getPseudoChromosomeNumber(reference), writeFn);

	iterator = reference_getPseudoChromosomeIterator(reference);
	while((pseudoChromosome = reference_getNextPseudoChromosome(iterator)) != NULL) {
		pseudoChromosome_writeBinaryRepresentation(pseudoChromosome, writeFn);
	}
	pseudoChromosome_destructPseudoAdjacencyIterator(iterator);
}

Reference *reference_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	Reference *reference = NULL;
	int32_t pseudoChromosomeNumber;
	Name name;

	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_REFERENCE) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		reference = reference_construct2(name, net);
		pseudoChromosomeNumber = binaryRepresentation_getInteger(binaryString);
		while(pseudoChromosomeNumber-- > 0) {
			pseudoChromosome_getName(pseudoChromosome_loadFromBinaryRepresentation(binaryString, reference));
		}
	}
	return reference;
}
