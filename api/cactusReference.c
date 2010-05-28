#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic reference functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static int reference_constructP(const void *o1, const void *o2) {
	return netMisc_nameCompare(pseudoChromosome_getName((PseudoChromosome *)o1),
			pseudoChromosome_getName((PseudoChromosome *)o2));
}

Reference *reference_construct(Net *net) {
	Reference *reference = st_malloc(sizeof(Reference));
	//Setup the basic structure - a sorted set of pseudo-chromosomes.
	reference->pseudoChromosomes = st_sortedSet_construct3(reference_constructP, NULL);
	//Link the reference and net.
	reference->net = net;
	net_setReference(net, reference);
	return reference;
}

Net *reference_getNet(Reference *reference) {
	return reference->net;
}

int32_t reference_getPseudoChromosomeNumber(Reference *reference) {
	return st_sortedSet_getLength(reference->pseudoChromosomes);
}

PseudoChromosome *reference_getPseudoChromosome(Reference *reference, Name name) {
	PseudoChromosome *pseudoChromosome;
	pseudoChromosome = pseudoChromosome_getStaticNameWrapper(name);
	return st_sortedSet_find(reference->pseudoChromosomes, pseudoChromosome);
}

PseudoChromosome *reference_getFirst(Reference *reference) {
	return st_sortedSet_getFirst(reference->pseudoChromosomes);
}

Reference_PseudoChromosomeIterator *reference_getPseudoChromosomeIterator(Reference *reference) {
	return st_sortedSet_getIterator(reference->pseudoChromosomes);
}

PseudoChromosome *reference_getNextPseudoChromosome(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
	return st_sortedSet_getNext(pseudoChromosomeIterator);
}

PseudoChromosome *reference_getPreviousPseudoChromosome(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
	return st_sortedSet_getPrevious(pseudoChromosomeIterator);
}

Reference_PseudoChromosomeIterator *reference_copyPseudoChromosomeIterator(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
	return st_sortedSet_copyIterator(pseudoChromosomeIterator);
}

void reference_destructPseudoChromosomeIterator(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
	st_sortedSet_destructIterator(pseudoChromosomeIterator);
}

st_Hash *reference_getEndToPseudoAdjacencyHash(Reference *reference) {
	st_Hash *hash = stHash_construct3(end_hashKey, end_hashEqualsKey, NULL, NULL);
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator =
		reference_getPseudoChromosomeIterator(reference);
	PseudoChromosome *pseudoChromosome;
	while((pseudoChromosome = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator =
			pseudoChromosome_getPseudoAdjacencyIterator(pseudoChromosome);
		PseudoAdjacency *pseudoAdjacency;
		while((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(pseudoAdjacencyIterator)) != NULL) {
			st_hash_insert(hash, pseudoAdjacency_get5End(pseudoAdjacency), pseudoAdjacency);
			st_hash_insert(hash, pseudoAdjacency_get3End(pseudoAdjacency), pseudoAdjacency);
		}
		pseudoChromosome_destructPseudoAdjacencyIterator(pseudoAdjacencyIterator);
	}
	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
	return hash;
}

void reference_check(Reference *reference) {
	Net *net = reference_getNet(reference);
	st_Hash *endsToPseudoAdjacencies = reference_getEndToPseudoAdjacencyHash(reference);

	//Going ends --> pseudo adjacencies.
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(end_isAttached(end) || end_isBlockEnd(end)) {
			PseudoAdjacency *pseudoAdjacency = st_hash_search(endsToPseudoAdjacencies, end);
			assert(pseudoAdjacency != NULL);
			assert(pseudoAdjacency_get5End(pseudoAdjacency) == end || pseudoAdjacency_get3End(pseudoAdjacency) == end);
		}
		else {
			assert(st_hash_search(endsToPseudoAdjacencies, end) == NULL); //check free stub end is not in the pseudo chromosomes..
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
			assert(st_hash_search(endsToPseudoAdjacencies, _5End) == pseudoAdjacency); //check these are represented in the hash.. so that the mapping is unique.
			assert(st_hash_search(endsToPseudoAdjacencies, _3End) == pseudoAdjacency);
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
	assert(i*2 == st_hash_size(endsToPseudoAdjacencies));

	//Cleanup
	st_hash_destruct(endsToPseudoAdjacencies);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void reference_destruct(Reference *reference) {
	net_removeReference(reference_getNet(reference), reference);
	PseudoChromosome *pseudoChromosome;
	while((pseudoChromosome = reference_getFirst(reference)) != NULL) {
			pseudoChromosome_destruct(pseudoChromosome);
	}
	st_sortedSet_destruct(reference->pseudoChromosomes);
	free(reference);
}

void reference_addPseudoChromosome(Reference *reference, PseudoChromosome *pseudoChromosome) {
	assert(st_sortedSet_find(reference->pseudoChromosomes, pseudoChromosome) == NULL);
	st_sortedSet_insert(reference->pseudoChromosomes, pseudoChromosome);
}

void reference_removePseudoChromosome(Reference *reference, PseudoChromosome *pseudoChromosome) {
	assert(st_sortedSet_find(reference->pseudoChromosomes, pseudoChromosome) != NULL);
	st_sortedSet_delete(reference->pseudoChromosomes, pseudoChromosome);
}

void reference_writeBinaryRepresentation(Reference *reference, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	Reference_PseudoChromosomeIterator *iterator;
	PseudoChromosome *pseudoChromosome;
	binaryRepresentation_writeElementType(CODE_REFERENCE, writeFn);
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

	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_REFERENCE) {
		binaryRepresentation_popNextElementType(binaryString);
		reference = reference_construct(net);
		pseudoChromosomeNumber = binaryRepresentation_getInteger(binaryString);
		while(pseudoChromosomeNumber-- > 0) {
			pseudoChromosome_loadFromBinaryRepresentation(binaryString, reference);
		}
	}
	return reference;
}
