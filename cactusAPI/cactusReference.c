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

Hash *reference_getEndToPseudoAdjacencyHash(Reference *reference) {
	Hash *hash = hash_construct3(end_hashKey, end_hashEqualsKey, NULL, NULL);
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator =
		reference_getPseudoChromosomeIterator(reference);
	PseudoChromosome *pseudoChromosome;
	while((pseudoChromosome = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator =
			pseudoChromosome_getPseudoAdjacencyIterator(pseudoChromosome);
		PseudoAdjacency *pseudoAdjacency;
		while((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(pseudoAdjacencyIterator)) != NULL) {
			hash_insert(hash, pseudoAdjacency_get5End(pseudoAdjacency), pseudoAdjacency);
			hash_insert(hash, pseudoAdjacency_get3End(pseudoAdjacency), pseudoAdjacency);
		}
		pseudoChromosome_destructPseudoAdjacencyIterator(pseudoAdjacencyIterator);
	}
	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
	return hash;
}

void reference_check(Reference *reference) {
	Net *net = reference_getNet(reference);
	Hash *endsToPseudoAdjacencies = reference_getEndToPseudoAdjacencyHash(reference);

	//Going ends --> pseudo adjacencies.
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(end_isAttached(end) || end_isBlockEnd(end)) {
			PseudoAdjacency *pseudoAdjacency = hash_search(endsToPseudoAdjacencies, end);
			assert(pseudoAdjacency != NULL);
			assert(pseudoAdjacency_get5End(pseudoAdjacency) == end || pseudoAdjacency_get3End(pseudoAdjacency) == end);
		}
		else {
			assert(hash_search(endsToPseudoAdjacencies, end) == NULL); //check free stub end is not in the pseudo chromosomes..
		}
	}
	net_destructEndIterator(endIterator);

	//Going pseudo-adjacencies --> ends.
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator = reference_getPseudoChromosomeIterator(reference);
	PseudoChromosome *pseudoChromosome;
	int32_t i = 0;
	while((pseudoChromosome = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		//Here we check the structure of the pseudo chromosome also..
		assert(pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome) > 0); //must be at least one adjacency
		assert(pseudoChromosome_get5End(pseudoChromosome) == pseudoAdjacency_get5End(pseudoChromosome_getFirst(pseudoChromosome))); //check the 5 end matches the 5 end of the first pseudo adjacency.
		assert(pseudoChromosome_get3End(pseudoChromosome) == pseudoAdjacency_get3End(pseudoChromosome_getLast(pseudoChromosome))); //check the 5 end matches the 5 end of the first pseudo adjacency.

		PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator = pseudoChromosome_getPseudoAdjacencyIterator(pseudoChromosome);
		PseudoAdjacency *pseudoAdjacency, *previousPseudoAdjacency = NULL;
		while((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(pseudoAdjacencyIterator)) != NULL) {
			End *_5End = pseudoAdjacency_get5End(pseudoAdjacency);
			End *_3End = pseudoAdjacency_get3End(pseudoAdjacency);
			assert(_5End == end_getPositiveOrientation(_5End)); //check they are positive orientation
			assert(_3End == end_getPositiveOrientation(_3End)); //check they are positive orientation
			assert(hash_search(endsToPseudoAdjacencies, _5End) == pseudoAdjacency); //check these are represented in the hash.. so that the mapping is unique.
			assert(hash_search(endsToPseudoAdjacencies, _3End) == pseudoAdjacency);
			assert(end_getGroup(_5End) != NULL); //check the groups are the same for both sides of the adjacency.
			assert(end_getGroup(_5End) == end_getGroup(_3End));
			i++;
			assert(group_getLink(end_getGroup(_5End)) == group_getLink(end_getGroup(_3End))); //check if there is a link, they are in the same link.
			if(previousPseudoAdjacency != NULL) { //check the adjacency spans the block...
				assert(pseudoAdjacency_get3End(previousPseudoAdjacency) == end_getOtherBlockEnd(_5End));
				assert(_5End == end_getOtherBlockEnd(pseudoAdjacency_get3End(previousPseudoAdjacency))); //do it the other way, just for fun.
			}
			previousPseudoAdjacency = pseudoAdjacency;
		}
		pseudoChromosome_destructPseudoAdjacencyIterator(pseudoAdjacencyIterator);
	}
	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
	assert(i*2 == hash_size(endsToPseudoAdjacencies));

	//Cleanup
	hash_destruct(endsToPseudoAdjacencies);
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

void reference_setNet(Reference *reference, Net *net) {
	net_removeReference(reference_getNet(reference), reference);
	reference->net = net;
	net_addReference(net, reference);
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
