/*
 * pseudoAdjacencyBuilder.c
 *
 *  Created on: 1 Apr 2010
 *      Author: benedictpaten
 *
 * This source file contains the functions which choose the pseudo adjacencies
 * in a pseudo chromosome. This is the vital bit!
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "reference.h"
#include "commonC.h"

typedef struct _adjacencyPair {
	End *end1;
	End *end2;
} AdjacencyPair;

AdjacencyPair *adjacencyPair_construct(End *end1, End *end2) {
	/*
	 * Constructs an adjacency pair. (an un-ordered pair of adjacencies).
	 */
	end1 = end_getPositiveOrientation(end1);
	end2 = end_getPositiveOrientation(end2);
	AdjacencyPair *adjacencyPair = malloc(sizeof(AdjacencyPair));
	if(netMisc_nameCompare(end_getName(end1), end_getName(end2) < 1)) {
		adjacencyPair->end2 = end1;
		adjacencyPair->end1 = end2;
	}
	else {
		adjacencyPair->end1 = end1;
		adjacencyPair->end2 = end2;
	}
	return adjacencyPair;
}

void adjacencyPair_destruct(AdjacencyPair *adjacencyPair) {
	/*
	 * Destructs an adjacency pair.
	 */
	free(adjacencyPair);
}

End *adjacencyPair_getEnd1(AdjacencyPair *adjacencyPair) {
	/*
	 * Gets the first end in the adjacency.
	 */
	return adjacencyPair->end1;
}

End *adjacencyPair_getEnd2(AdjacencyPair *adjacencyPair) {
	/*
	 * Gets the second end in the adjacency.
	 */
	return adjacencyPair->end2;
}

double adjacencyPair_getStrengthOfAdjacencyPair(AdjacencyPair *adjacencyPair) {
	/*
	 * Returns the strength of the adjacency between the two ends.
	 */
	End *end1 = adjacencyPair_getEnd1(adjacencyPair);
	End *end2 = adjacencyPair_getEnd2(adjacencyPair);
	int32_t i = end_getInstanceNumber(end1) + end_getInstanceNumber(end2);
	int32_t j = 0;
	//Walk through all the adjacencies and give increment j by 2 if
	//the adjacency is between the two ends.
	End_InstanceIterator *iterator = end_getInstanceIterator(end1);
	Cap *cap;
	while((cap = end_getNext(iterator)) != NULL) {
		Cap *cap2 = cap_getAdjacency(cap);
		assert(cap2 != NULL);
		if(end_getPositiveOrientation(cap_getEnd(cap2)) == end_getPositiveOrientation(end2)) {
			j += 2;
		}
	}
	end_destructInstanceIterator(iterator);
	return ((double)j) / i;
}

int adjacencyPair_cmpFnByStrength(AdjacencyPair **adjacencyPair1, AdjacencyPair **adjacencyPair2) {
	/*
	 * Compares two adjacencies and scores them according to amount of support.
	 */
	int32_t i = adjacencyPair_getStrengthOfAdjacencyPair(*adjacencyPair1);
	int32_t j = adjacencyPair_getStrengthOfAdjacencyPair(*adjacencyPair2);
	return i - j;
}

uint32_t adjacencyPair_hashKey(AdjacencyPair *adjacencyPair) {
	/*
	 * Hash function for adjacency pairs.
	 */
	return end_getName(adjacencyPair_getEnd1(adjacencyPair));
}

int32_t adjacencyPair_hashEqual(AdjacencyPair *adjacencyPair1, AdjacencyPair *adjacencyPair2) {
	/*
	 * Hash equals function for adjacency pairs.
	 */
	return adjacencyPair_getEnd1(adjacencyPair1) == adjacencyPair_getEnd1(adjacencyPair2) &&
			adjacencyPair_getEnd2(adjacencyPair1) == adjacencyPair_getEnd2(adjacencyPair2);
}

void addAdjacencyPairToHash(struct hashtable *adjacenciesHash, AdjacencyPair *adjacencyPair) {
	/*
	 * Adds the adjacency pair to the hash.
	 */
	hashtable_insert(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair), adjacencyPair);
	hashtable_insert(adjacenciesHash, adjacencyPair_getEnd2(adjacencyPair), adjacencyPair);
}

void removeAdjacency(struct hashtable *adjacencies, AdjacencyPair *adjacencyPair) {
	/*
	 * Removes the adjacency pair from the hash and destroys it.
	 */
	hashtable_remove(adjacencies, adjacencyPair_getEnd1(adjacencyPair), 0);
	hashtable_remove(adjacencies, adjacencyPair_getEnd2(adjacencyPair), 0);
	adjacencyPair_destruct(adjacencyPair);
}

struct List *makeListOfAdjacencyPairs(Net *net) {
	/*
	 * Get a list of all adjacencies between ends, excluding free stubs, in the net.
	 */
	End *end1;
	struct List *adjacencies = constructEmptyList(0, (void (*)(void *))adjacencyPair_destruct);
	Net_EndIterator *endIterator = net_getEndIterator(net);
	struct hashtable *adjacenciesHash = create_hashtable(0,
			(uint32_t (*)(void *))adjacencyPair_hashKey,
			(int32_t (*)(void *, void *))adjacencyPair_hashEqual, NULL, NULL);
	while((end1 = net_getNextEnd(endIterator)) != NULL) {
		End_InstanceIterator *capIterator = end_getInstanceIterator(end1);
		Cap *cap;
		while((cap = end_getNext(capIterator)) != NULL) {
			Cap *cap2 = cap_getAdjacency(cap);
			assert(cap2 != NULL);
			End *end2 = cap_getEnd(cap2);
			AdjacencyPair *adjacencyPair = adjacencyPair_construct(end1, end2);
			if(hashtable_search(adjacenciesHash, adjacencyPair) == NULL) {
				hashtable_insert(adjacenciesHash, adjacencyPair, adjacencyPair);
				listAppend(adjacencies, adjacencyPair);
			}
			else {
				adjacencyPair_destruct(adjacencyPair);
			}
		}
		end_destructInstanceIterator(capIterator);
	}
	net_destructEndIterator(endIterator);
	hashtable_destroy(adjacenciesHash, 0, 0);
	return adjacencies;
}

struct hashtable *choosePairing(struct List *adjacencies, Net *net) {
	/*
	 * Greedily picks the adjacencies from the list such that each end has one adjacency.
	 */
	struct hashtable *adjacenciesHash =
			create_hashtable(1, hashtable_key, hashtable_equalKey, NULL, NULL);
	int32_t i;
	for(i=0; i<adjacencies->length; i++) {
		AdjacencyPair *adjacencyPair;
		if(hashtable_search(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair)) == NULL &&
		   hashtable_search(adjacenciesHash, adjacencyPair_getEnd2(adjacencyPair)) == NULL) {
			addAdjacencyPairToHash(adjacenciesHash, adjacencyPair);
		}
		else {
			adjacencyPair_destruct(adjacencyPair);
		}
	}
	return adjacenciesHash;
}

void addPseudoPseudoAdjacencies(struct List *adjacencies, struct hashtable *adjacenciesHash, Net *net) {
	/*
	 * Adds pseudo-pseudo adjacencies to ends that have no valid adjacency pair in the hash. Added adjacencies
	 * are put in the hash and the list.
	 */
	End *end1;
	End *end2 = NULL;
	Net_EndIterator *endIterator = net_getEndIterator(net);
	while((end1 = net_getNextEnd(endIterator)) != NULL) {
		if(hashtable_search(adjacenciesHash, end1) == NULL) {
			if(end2 == NULL) {
				end2 = end1;
			}
			else {
				AdjacencyPair *adjacencyPair = adjacencyPair_construct(end1, end2);
				addAdjacencyPairToHash(adjacenciesHash, adjacencyPair);
				listAppend(adjacencies, adjacencyPair);
				end2 = NULL;
			}
		}
	}
}

void extendComponent(End *end, struct List *component, struct hashtable *componentsHash, struct hashtable *adjacencies) {
	assert(hashtable_search(componentsHash, end) == NULL);
	listAppend(component, end);
	hashtable_insert(componentsHash, end, component);
	AdjacencyPair *adjacencyPair =  hashtable_search(adjacencies, end);
	if(adjacencyPair != NULL) {
		extendComponent(adjacencyPair_getEnd1(adjacencyPair) != end ?
				adjacencyPair_getEnd1(adjacencyPair) : adjacencyPair_getEnd1(adjacencyPair),
				component, componentsHash, adjacencies);
	}
	if(end_isBlockEnd(end)) {
		extendComponent(block_getRightEnd(end_getBlock(end)), component, componentsHash, adjacencies);
	}
}

struct List *getConnectedComponents(struct hashtable *adjacencies, Net *net) {
	/*
	 * Gets a list of connected components of ends linked by adjacencies and blocks.
	 */
	struct List *components = constructEmptyList(0, NULL);
	struct hashtable *componentsHash =
				create_hashtable(1, hashtable_key, hashtable_equalKey, NULL, NULL);
	End *end;
	Net_EndIterator *endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		struct List *component = hashtable_search(componentsHash, end);
		if(component == NULL) {
			component = constructEmptyList(0, NULL);
			listAppend(components, component);
			extendComponent(end, component, componentsHash, adjacencies);
		}
	}
	hashtable_destroy(componentsHash, 0, 0);
	return components;
}

void getAttachedStubEndsInComponentP(End **end1, End **end2, End *end3) {
	if(end_isStubEnd(end3)) {
		assert(end_isAttached(end3));
		if(*end1 == NULL) {
			assert(*end2 == NULL);
			*end2 = end3;
		}
		else {
			*end1 = end3;
		}
	}
}

void getAttachedStubEndsInComponent(struct List *component, End **end1, End **end2) {
	/*
	 * Gets the one or two ends are in the component and assigns them, first to end1, then
	 * to end2. If there isn't an end, sets it null.
	 */
	int32_t i;
	*end1 = NULL;
	*end2 = NULL;
	for(i=0; i<component->length; i++) {
		End *end3 = component->list[i];
		if(end_isStubEnd(end3)) {
			assert(end_isAttached(end3));
			if(*end1 == NULL) {
				assert(*end2 == NULL);
				*end2 = end3;
			}
			else {
				*end1 = end3;
			}
		}
	}
}

void switchAdjacencies(AdjacencyPair *adjacencyPair1, AdjacencyPair *adjacencyPair2, struct hashtable *adjacencies,
		void **adjacencyPair3, void **adjacencyPair4) {
	/*
	 * Removes two adjacencies and creates two more (in random configuration of the two possible configs).
	 */
	*adjacencyPair3 = adjacencyPair_construct(adjacencyPair_getEnd1(adjacencyPair1), adjacencyPair_getEnd1(adjacencyPair2));
	*adjacencyPair4 = adjacencyPair_construct(adjacencyPair_getEnd2(adjacencyPair1), adjacencyPair_getEnd2(adjacencyPair2));
	addAdjacencyPairToHash(adjacencies, *adjacencyPair3);
	addAdjacencyPairToHash(adjacencies, *adjacencyPair4);
	adjacencyPair_destruct(adjacencyPair1);
	adjacencyPair_destruct(adjacencyPair2);
}

struct List *mergeCycle(struct List *component1, struct List *component2, int32_t index1, int32_t index2,
		struct hashtable *adjacencies) {
	/*
	 * Merges together two cycles or a cycle and a contig.
	 */
	struct List *mergedComponent = constructEmptyList(0, NULL);
	switchAdjacencies(component1->list[index1], component2->list[index2], adjacencies,
			&(component1->list[index1]), &(component2->list[index2]));
	listAppendArray(mergedComponent, component1->list, component1->length);
	listAppendArray(mergedComponent, component2->list, component2->length);
	destructList(component1);
	destructList(component2);
	return mergedComponent;
}

void mergeCycles(struct List *components, struct hashtable *adjacencies, Net *net) {
	/*
	 * Merges cycles not containing an attached stub end into components containing attached stub ends,
	 * returning a list of components such that each component contains two attached stub ends.
	 * Updates adjacencies as we go.
	 */
	int32_t i;
	End *end1, *end2;
	struct List *cycles = constructEmptyList(0, NULL);
	struct List *contigs = constructEmptyList(0, NULL);
	struct List *component;
	while(components->length) {
		component = components->list[--components->length];
		getAttachedStubEndsInComponent(component, &end1, &end2);
		if(end1 != NULL) {
			assert(end2 != NULL);
			listAppend(contigs, component);
		}
		else {
			listAppend(cycles, component);
		}
	}
	for(i=0; i<cycles->length; i++) {
		contigs->list[0] = mergeCycle(contigs->list[0], cycles->list[i], 0, 0, adjacencies);
	}
	listAppendArray(components, contigs->list, contigs->length);
	destructList(contigs);
	destructList(cycles);
}

AdjacencyPair *pickAdjacencyPairToBreak(struct List *component, struct hashtable *adjacencies) {
	/*
	 *
	 */
	assert(component->length > 0);
	End *end = component->list[0];
	return hashtable_search(adjacencies, end);
}

struct List *correctAttachedStubEndPairing(struct List *components, struct hashtable *adjacencies, Net *net, Reference *reference) {
	/*
	 * If a-b and c-d are two components respectively containing stub ends a and b and c and d and the correct pairing
	 * is a-c b-d then two adjacencies are broken and two are created to make this happen.
	 */
	struct hashtable *correctEndPairing = create_hashtable(1, hashtable_key, hashtable_equalKey, NULL, NULL);
	PseudoChromosome *pseudoChromosome;
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator = reference_getPseudoChromosomeIterator(reference);
	while((pseudoChromosome = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		AdjacencyPair *adjacencyPair = adjacencyPair_construct(pseudoChromosome_get5End(pseudoChromosome), pseudoChromosome_get3End(pseudoChromosome));
	    addAdjacencyPairToHash(correctEndPairing, adjacencyPair);
	}
	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);

	//Break apart contigs with mispaired ends.
	int32_t i, j;
	End *end1, *end2;
	for(i=0; i<components->length; i++) { //Get mispaired contigs, breaking them up.
		struct List *component = components->list[i];
		getAttachedStubEndsInComponent(components->list[i], &end1, &end2);
		assert(end1 != NULL);
		assert(end2 != NULL);
		AdjacencyPair *adjacencyPair1 = hashtable_search(correctEndPairing, end_getPositiveOrientation(end1)); //defensive
		AdjacencyPair *adjacencyPair2 = hashtable_search(correctEndPairing, end_getPositiveOrientation(end2));
		if(adjacencyPair1 != adjacencyPair2) {
			//Break an adjacency
			removeAdjacency(adjacencies, pickAdjacencyPairToBreak(component, adjacencies));
		}
	}

	//Recalculate the end components.
	components = getConnectedComponents(adjacencies, net);

	//Now go and reconnect the right components.
	for(i=1; i<components->length; i++) { //Get mispaired contigs, breaking them up.
		getAttachedStubEndsInComponent(components->list[i], &end1, &end2);
		assert(end1 != NULL);
		if(end2 == NULL) { //get the other component
			AdjacencyPair *adjacencyPair1 = hashtable_search(correctEndPairing, end_getPositiveOrientation(end1));
			for(j=0; j<i; j++) {
				End *end3, *end4;
				getAttachedStubEndsInComponent(components->list[j], &end3, &end4);
				if(end4 == NULL) {
					assert(end3 != NULL);
					AdjacencyPair *adjacencyPair2 = hashtable_search(correctEndPairing, end_getPositiveOrientation(end3));
					if(adjacencyPair1 == adjacencyPair2) {
						addAdjacencyPairToHash(adjacencies, adjacencyPair_construct(getFreeEnd(components->list[i]), getFreeEnd(components->list[j])));
					}
				}
			}
		}
	}
	hashtable_destroy(correctEndPairing, 0, 1);

	return getConnectedComponents(adjacencies, net);
}

void fillInPseudoAdjacenciesP(End *end1, PseudoChromosome *pseudoChromosome, struct hashtable *adjacencies) {
	AdjacencyPair *adjacencyPair = hashtable_search(adjacencies, end1);
	assert(adjacencyPair != NULL);
	//construct the new pseudo-adjacency!
	End *end2 = end1 == adjacencyPair_getEnd1(adjacencyPair) ? adjacencyPair_getEnd2(adjacencyPair) : adjacencyPair_getEnd1(adjacencyPair);
	pseudoAdjacency_construct(end1, end2, pseudoChromosome);
	if(pseudoChromosome_get3End(pseudoChromosome) != end2) {
		assert(end_isBlockEnd(end2));
		fillInPseudoAdjacenciesP(block_getRightEnd(end_getBlock(end2)), pseudoChromosome, adjacencies);
	}
}

static void fillInPseudoAdjacencies(Net *net, Reference *reference, struct hashtable *adjacencies) {
	/*
	 * Walks through the list of pseudo-chromosomes, then
	 * iterates through the adjacencies populating each psuedo chromosome with adjacencies.
	 */
	PseudoChromosome *pseudoChromosome;
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator = reference_getPseudoChromosomeIterator(reference);
	while((pseudoChromosome = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		//Now we fill in the pseudo-adjacencies
		End *end = pseudoChromosome_get5End(pseudoChromosome);
		fillInPseudoAdjacenciesP(end, pseudoChromosome, adjacencies);
	}
	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
}

void makePseudoAdjacencies(Net *net, Reference *reference) {
	/*
	 * Makes the pseudo adjacencies for each pseudo chromosome in reference.
	 */
	struct List *adjacencies = makeListOfAdjacencyPairs(net);
	qsort(adjacencies->list, adjacencies->length, sizeof(void *), (int (*)(const void *v, const void *))adjacencyPair_cmpFnByStrength);
	struct hashtable *adjacenciesHash = choosePairing(adjacencies, net);
	addPseudoPseudoAdjacencies(adjacencies, adjacenciesHash, net);
	struct List *components = getConnectedComponents(adjacenciesHash, net);
	mergeCycles(components, adjacenciesHash, net);
	correctAttachedStubEndPairing(components, adjacenciesHash, net);
	fillInPseudoAdjacencies(net, reference, adjacenciesHash);
	//Cleanup
	destructList(components);
	hashtable_destroy(adjacenciesHash, 0, 0);
	destructList(adjacencies); //this also cleans up all the actual adjacency structs.
}
