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
	/*
	 * Struct to represent a pair of 'adjacent' ends.
	 */
	End *end1;
	End *end2;
} AdjacencyPair;

static End *adjacencyPair_getEnd1(AdjacencyPair *adjacencyPair) {
	/*
	 * Gets the first end in the adjacency.
	 */
	return adjacencyPair->end1;
}

static End *adjacencyPair_getEnd2(AdjacencyPair *adjacencyPair) {
	/*
	 * Gets the second end in the adjacency.
	 */
	return adjacencyPair->end2;
}

static AdjacencyPair *adjacencyPair_construct(End *end1, End *end2) {
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
	assert(netMisc_nameCompare(end_getName(adjacencyPair_getEnd1(adjacencyPair)), end_getName(adjacencyPair_getEnd2(adjacencyPair))) >= 0);
	return adjacencyPair;
}

static void adjacencyPair_destruct(AdjacencyPair *adjacencyPair) {
	/*
	 * Destructs an adjacency pair.
	 */
	free(adjacencyPair);
}

static double adjacencyPair_getStrengthOfAdjacencyPair(AdjacencyPair *adjacencyPair) {
	/*
	 * Returns the strength of the adjacency between the two ends.
	 */
	End *end1 = adjacencyPair_getEnd1(adjacencyPair);
	End *end2 = adjacencyPair_getEnd2(adjacencyPair);
	int32_t i = end_getInstanceNumber(end1) + end_getInstanceNumber(end2);
	if(i == 0) { //avoid divide by zero.
		return 0;
	}
	int32_t j = 0;
	//Walk through all the adjacencies and give increment j by 2 if
	//the adjacency is between the two ends.
	End_InstanceIterator *iterator = end_getInstanceIterator(end1);
	Cap *cap;
	while((cap = end_getNext(iterator)) != NULL) {
		Cap *cap2 = cap_getAdjacency(cap);
		assert(cap2 != NULL);
		if(end_getPositiveOrientation(cap_getEnd(cap2)) == end2) {
			j += 2;
		}
	}
	end_destructInstanceIterator(iterator);
	return ((double)j) / i;
}

static int adjacencyPair_cmpFnByStrength(AdjacencyPair **adjacencyPair1, AdjacencyPair **adjacencyPair2) {
	/*
	 * Compares two adjacencies such that the adjacency with stronger support is greater.
	 */
	int32_t i = adjacencyPair_getStrengthOfAdjacencyPair(*adjacencyPair1);
	int32_t j = adjacencyPair_getStrengthOfAdjacencyPair(*adjacencyPair2);
	return i - j;
}

static uint32_t adjacencyPair_hashKey(AdjacencyPair *adjacencyPair) {
	/*
	 * Hash function for adjacency pairs.
	 */
	return end_getName(adjacencyPair_getEnd1(adjacencyPair));
}

static int32_t adjacencyPair_hashEqual(AdjacencyPair *adjacencyPair1, AdjacencyPair *adjacencyPair2) {
	/*
	 * Hash equals function for adjacency pairs.
	 */
	return adjacencyPair_getEnd1(adjacencyPair1) == adjacencyPair_getEnd1(adjacencyPair2) &&
			adjacencyPair_getEnd2(adjacencyPair1) == adjacencyPair_getEnd2(adjacencyPair2);
}

static struct List *makeListOfAdjacencyPairs(Net *net) {
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
			if(hashtable_search(adjacenciesHash, adjacencyPair) == NULL) { //
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

static void addAdjacencyPairToHash(struct hashtable *adjacenciesHash, AdjacencyPair *adjacencyPair) {
	/*
	 * Adds the adjacency pair to the hash.
	 */
	hashtable_insert(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair), adjacencyPair);
	hashtable_insert(adjacenciesHash, adjacencyPair_getEnd2(adjacencyPair), adjacencyPair);
}

static void removeAdjacency(struct hashtable *adjacencies, AdjacencyPair *adjacencyPair) {
	/*
	 * Removes the adjacency pair from the hash and destroys it.
	 */
	hashtable_remove(adjacencies, adjacencyPair_getEnd1(adjacencyPair), 0);
	hashtable_remove(adjacencies, adjacencyPair_getEnd2(adjacencyPair), 0);
	adjacencyPair_destruct(adjacencyPair);
}

static struct hashtable *choosePairing(struct List *adjacencies, Net *net) {
	/*
	 * Greedily picks the adjacencies from the list such that each end has one adjacency.
	 * Destroys the input list in the process.
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
	destructList(adjacencies);
	return adjacenciesHash;
}

static void addPseudoPseudoAdjacencies(struct hashtable *adjacenciesHash, Net *net) {
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
				end2 = NULL;
			}
		}
	}
	assert(end2 == NULL); //there must be an even number of pairs.
}

static void extendComponent(End *end, struct List *component, struct hashtable *componentsHash, struct hashtable *adjacencies) {
	/*
	 * Sub-function of get connected components, extends the component.
	 */
	assert(hashtable_search(componentsHash, end) == NULL);
	listAppend(component, end);
	hashtable_insert(componentsHash, end, component);
	AdjacencyPair *adjacencyPair = hashtable_search(adjacencies, end);
	if(adjacencyPair != NULL) {
		extendComponent(adjacencyPair_getEnd1(adjacencyPair) != end ?
				adjacencyPair_getEnd1(adjacencyPair) : adjacencyPair_getEnd2(adjacencyPair),
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
	struct List *components = constructEmptyList(0, (void (*)(void *))destructList);
	struct hashtable *componentsHash =
				create_hashtable(1, hashtable_key, hashtable_equalKey, NULL, NULL);
	End *end;
	Net_EndIterator *endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) { //iterates over the positive oriented ends.
		struct List *component = hashtable_search(componentsHash, end);
		if(component == NULL) {
			component = constructEmptyList(0, NULL);
			listAppend(components, component);
			extendComponent(end, component, componentsHash, adjacencies);
		}
	}
	net_destructEndIterator(endIterator);
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
	 * Gets the one or two stub ends in the component and assigns them, first to end1, then
	 * to end2. If there isn't an end, sets it null.
	 */
	int32_t i;
	*end1 = NULL;
	*end2 = NULL;
	for(i=0; i<component->length; i++) {
		End *end3 = component->list[i];
		if(end_isStubEnd(end3)) {
			assert(end_isAttached(end3));
			assert(*end2 == NULL);
			if(*end1 == NULL) {
				*end1 = end3;
			}
			else {
				*end2 = end3;
			}
		}
	}
}

void switchAdjacencies(AdjacencyPair *adjacencyPair1, AdjacencyPair *adjacencyPair2, struct hashtable *adjacencies,
		AdjacencyPair **adjacencyPair3, AdjacencyPair **adjacencyPair4) {
	/*
	 * Removes two adjacencies and destroys them and creates two more (in random configuration of the two possible configs).
	 */
	*adjacencyPair3 = adjacencyPair_construct(adjacencyPair_getEnd1(adjacencyPair1), adjacencyPair_getEnd1(adjacencyPair2));
	*adjacencyPair4 = adjacencyPair_construct(adjacencyPair_getEnd2(adjacencyPair1), adjacencyPair_getEnd2(adjacencyPair2));
	removeAdjacency(adjacencies, adjacencyPair1);
	removeAdjacency(adjacencies, adjacencyPair2);
	addAdjacencyPairToHash(adjacencies, *adjacencyPair3);
	addAdjacencyPairToHash(adjacencies, *adjacencyPair4);
}

void splitIntoContigsAndCycles(struct List *components, struct List **contigs, struct List **cycles) {
	/*
	 * Divides the list of components into rings and components terminated by a pair of stubs.
	 */
	int32_t i;
	End *end1, *end2;
	*cycles = constructEmptyList(0, NULL);
	*contigs = constructEmptyList(0, NULL);
	struct List *component;
	for(i=0; i<components->length; i++) {
		component = components->list[--components->length];
		getAttachedStubEndsInComponent(component, &end1, &end2);
		if(end1 != NULL) {
			assert(end2 != NULL);
			listAppend(*contigs, component);
		}
		else {
			listAppend(*cycles, component);
		}
	}
}

AdjacencyPair *getAdjacencyPair(struct List *component, struct hashtable *adjacencies) {
	/*
	 * Gets an adjacency pair from a component.
	 */
	assert(component->length > 0);
	AdjacencyPair *adjacencyPair = hashtable_search(adjacencies, component->list[0]);
	assert(adjacencyPair != NULL);
	return adjacencyPair;
}

void mergeCycles(struct hashtable *adjacencies, Net *net) {
	/*
	 * Merges cycles not containing an attached stub end into components containing attached stub ends,
	 * Updates adjacencies as we go and destroys the list of components.
	 */
	int32_t i;
	struct List *contigs, *cycles;
	AdjacencyPair *adjacencyPair1, *adjacencyPair2;
	AdjacencyPair *adjacencyPair3, *adjacencyPair4;

	//Get the comopnents to merge.
	struct List *components = getConnectedComponents(adjacencies, net);

	//Get the cycles and contigs
	splitIntoContigsAndCycles(components, &contigs, &cycles);

	//Get the adjacency pair in a cycle with which to merge.
	assert(contigs->length > 0);
	adjacencyPair1 = getAdjacencyPair(contigs->list[0], adjacencies);

	for(i=0; i<cycles->length; i++) { //Now merge the cycles
		adjacencyPair2 = getAdjacencyPair(cycles->list[i], adjacencies); //the adjacency to break
		switchAdjacencies(adjacencyPair1, adjacencyPair2, adjacencies, &adjacencyPair3, &adjacencyPair4);
		adjacencyPair1 = adjacencyPair3; //now ensure we have new adjacency to break for the next round
	}
	destructList(contigs);
	destructList(cycles);
	destructList(components);
}

static struct hashtable *getCorrectEndPairing(Reference *reference) {
	/*
	 * Constructs a hash set of adjacency pairs, each pair being the 5 and 3 prime
	 * ends of a pseudo chromosome in the reference.
	 */
	struct hashtable *correctEndPairing = create_hashtable(0, (uint32_t (*)(void *))adjacencyPair_hashKey,
																  (int32_t (*)(void *, void *))adjacencyPair_hashEqual, (void (*)(void *))adjacencyPair_destruct, NULL);
	PseudoChromosome *pseudoChromosome;
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator = reference_getPseudoChromosomeIterator(reference);
	while((pseudoChromosome = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		AdjacencyPair *adjacencyPair = adjacencyPair_construct(pseudoChromosome_get5End(pseudoChromosome), pseudoChromosome_get3End(pseudoChromosome));
		addAdjacencyPairToHash(correctEndPairing, adjacencyPair);
	}
	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
	return correctEndPairing;
}

void breakApartMisPairedContigs(struct hashtable *correctEndPairing,
		struct hashtable *adjacencies, Net *net) {
	/*
	 * Breaks apart mispaired contigs.
	 */
	int32_t i;
	End *end1, *end2;
	struct List *components = getConnectedComponents(adjacencies, net);
	for(i=0; i<components->length; i++) { //Get mispaired contigs, breaking them up.
		struct List *component = components->list[i];
		getAttachedStubEndsInComponent(components->list[i], &end1, &end2);
		assert(end1 != NULL);
		assert(end2 != NULL);
		AdjacencyPair *adjacencyPair1 = hashtable_search(correctEndPairing, end1); //defensive
		AdjacencyPair *adjacencyPair2 = hashtable_search(correctEndPairing, end2);
		if(adjacencyPair1 != adjacencyPair2) {
			//Break an adjacency
			removeAdjacency(adjacencies, getAdjacencyPair(component, adjacencies));
		}
	}
	destructList(components);
}

End *getFreeEnd(struct List *component, struct hashtable *adjacencies) {
	/*
	 * Gets the end in the component with no adjacency.
	 */
	int32_t i;
	End *end = NULL;
	for(i=0; i<component->length; i++) {
		End *end2 = component->list[i];
		if(hashtable_search(adjacencies, end2) == NULL) {
			assert(end == NULL);
			end = end2;
		}
	}
	assert(end != NULL);
	return end;
}

void pairBrokenContigs(struct hashtable *correctEndPairing,
		struct hashtable *adjacencies, Net *net) {
	/*
	 * Pairs broken contigs together so that each has two stubs, paired
	 * according to the correct end pairing.
	 */
	//Now go and reconnect the right components.
	int32_t i, j;
	End *end1, *end2;
	struct List *components = getConnectedComponents(adjacencies, net);
	for(i=1; i<components->length; i++) { //Get mispaired contigs, breaking them up.
		getAttachedStubEndsInComponent(components->list[i], &end1, &end2);
		assert(end1 != NULL);
		if(end2 == NULL) { //get the other component
			AdjacencyPair *adjacencyPair1 = hashtable_search(correctEndPairing, end1);
			for(j=0; j<i; j++) {
				End *end3, *end4;
				getAttachedStubEndsInComponent(components->list[j], &end3, &end4);
				if(end4 == NULL) {
					assert(end3 != NULL);
					AdjacencyPair *adjacencyPair2 = hashtable_search(correctEndPairing, end3);
					if(adjacencyPair1 == adjacencyPair2) {
						addAdjacencyPairToHash(adjacencies,
								adjacencyPair_construct(getFreeEnd(components->list[i], adjacencies),
								getFreeEnd(components->list[j], adjacencies)));
					}
				}
			}
		}
	}
	destructList(components);
}

void correctAttachedStubEndPairing(struct hashtable *adjacencies, Net *net, Reference *reference) {
	/*
	 * If a-b and c-d are two components respectively containing stub ends a and b and c and d and the correct pairing
	 * is a-c b-d then two adjacencies are broken and two are created to make this happen.
	 * Destructs the list of components as we go.
	 */
	//Get the correct end pairings of the pseudo-telomeres in a hash of pseudo adjacencies.
	struct hashtable *correctEndPairing = getCorrectEndPairing(reference);
	breakApartMisPairedContigs(correctEndPairing, adjacencies, net);
	pairBrokenContigs(correctEndPairing, adjacencies, net);
	hashtable_destroy(correctEndPairing, 0, 1);
}

static void fillInPseudoAdjacenciesP(End *end1, PseudoChromosome *pseudoChromosome, struct hashtable *adjacencies) {
	/*
	 * Traverses the connected component of the contig constructing the pseudo-adjacencies in the pseudo-chromosome.
	 */
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
	addPseudoPseudoAdjacencies(adjacenciesHash, net);
	mergeCycles(adjacenciesHash, net);
	correctAttachedStubEndPairing(adjacenciesHash, net, reference);
	fillInPseudoAdjacencies(net, reference, adjacenciesHash);
	//Cleanup
	hashtable_destroy(adjacenciesHash, 0, 0);
	destructList(adjacencies); //this also cleans up all the actual adjacency structs.
}
