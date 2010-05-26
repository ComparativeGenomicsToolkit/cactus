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
	assert(end1 != end2);
	assert(end_getName(end1) != end_getName(end2));
	AdjacencyPair *adjacencyPair = mallocLocal(sizeof(AdjacencyPair));
	if(netMisc_nameCompare(end_getName(end1), end_getName(end2)) <= 0) {
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

End *adjacencyPair_getOtherEnd(AdjacencyPair *adjacencyPair, End *end) {
	/*
	 * Gets the other end in the adjacency pair.
	 */
	assert(end_getOrientation(end));
	assert(adjacencyPair_getEnd1(adjacencyPair) == end || adjacencyPair_getEnd2(adjacencyPair) == end);
	return adjacencyPair_getEnd1(adjacencyPair) == end ? adjacencyPair_getEnd2(adjacencyPair) : adjacencyPair_getEnd1(adjacencyPair);
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
		if(cap2 != NULL && end_getPositiveOrientation(cap_getEnd(cap2)) == end2) {
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
	double i = adjacencyPair_getStrengthOfAdjacencyPair(*adjacencyPair1);
	double j = adjacencyPair_getStrengthOfAdjacencyPair(*adjacencyPair2);
	return i > j ? 1 : (i < j ? -1 : 0);
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
	struct List *adjacencies = constructEmptyList(0, NULL);
	Net_EndIterator *endIterator = net_getEndIterator(net);
	st_Hash *adjacenciesHash = st_hash_construct3(
			(uint32_t (*)(void *))adjacencyPair_hashKey,
			(int32_t (*)(void *, void *))adjacencyPair_hashEqual, NULL, NULL);
	while((end1 = net_getNextEnd(endIterator)) != NULL) {
		if(end_isBlockEnd(end1) || end_isAttached(end1)) {
			assert(end_getPositiveOrientation(end1) == end1); //they are always in the positive orientation, right?
			End_InstanceIterator *capIterator = end_getInstanceIterator(end1);
			Cap *cap;
			while((cap = end_getNext(capIterator)) != NULL) {
				Cap *cap2 = cap_getAdjacency(cap); //we allow both internal inferred and leaf adjacencies.
				if(cap2 != NULL) {
					End *end2 = end_getPositiveOrientation(cap_getEnd(cap2));
					if(end1 != end2 && (end_isBlockEnd(end2) || end_isAttached(end2))) { //we don't allow free stubs or adjacency pairs which create self loops, as we can traverse a node twice!
						AdjacencyPair *adjacencyPair = adjacencyPair_construct(end1, end2);
						if(st_hash_search(adjacenciesHash, adjacencyPair) == NULL) { //
							st_hash_insert(adjacenciesHash, adjacencyPair, adjacencyPair);
							listAppend(adjacencies, adjacencyPair);
						}
						else {
							adjacencyPair_destruct(adjacencyPair);
						}
					}
				}
			}
			end_destructInstanceIterator(capIterator);
		}
	}
	net_destructEndIterator(endIterator);
	st_hash_destruct(adjacenciesHash);
	return adjacencies;
}

static void addAdjacencyPairToHash(st_Hash *adjacenciesHash, AdjacencyPair *adjacencyPair) {
	/*
	 * Adds the adjacency pair to the hash.
	 */
	st_hash_insert(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair), adjacencyPair);
	st_hash_insert(adjacenciesHash, adjacencyPair_getEnd2(adjacencyPair), adjacencyPair);
}

static void removeAndDestructAdjacency(st_Hash *adjacencies, AdjacencyPair *adjacencyPair) {
	/*
	 * Removes the adjacency pair from the hash and destroys it.
	 */
	assert(st_hash_remove(adjacencies, adjacencyPair_getEnd1(adjacencyPair)) == adjacencyPair);
	assert(st_hash_remove(adjacencies, adjacencyPair_getEnd2(adjacencyPair)) == adjacencyPair);
	adjacencyPair_destruct(adjacencyPair);
}

static st_Hash *choosePairing(struct List *adjacencies, Net *net) {
	/*
	 * Greedily picks the adjacencies from the list such that each end has one adjacency.
	 * Destroys the input list in the process.
	 */
	st_Hash *adjacenciesHash = st_hash_construct();
#ifdef BEN_DEBUG
	double strength = INT32_MAX;
#endif
	while(adjacencies->length > 0) {
		AdjacencyPair *adjacencyPair = adjacencies->list[--adjacencies->length];
#ifdef BEN_DEBUG
		double d = adjacencyPair_getStrengthOfAdjacencyPair(adjacencyPair);
		assert(d <= strength);
		strength = d;
#endif
		if(st_hash_search(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair)) == NULL &&
		   st_hash_search(adjacenciesHash, adjacencyPair_getEnd2(adjacencyPair)) == NULL) {
			addAdjacencyPairToHash(adjacenciesHash, adjacencyPair);
		}
		else {
			adjacencyPair_destruct(adjacencyPair);
		}
	}
	assert(adjacencies->length == 0);
	destructList(adjacencies);
	return adjacenciesHash;
}

void adjacenciesHash_cleanUp(Net *net, st_Hash *adjacencies) {
	/*
	 * Frees the adjacencies pairs in the adjacencies hash safely.
	 */
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		AdjacencyPair *adjacencyPair = st_hash_search(adjacencies, end);
		if(adjacencyPair != NULL) {
			removeAndDestructAdjacency(adjacencies, adjacencyPair);
		}
	}
	net_destructEndIterator(endIterator);
	assert(st_hash_size(adjacencies) == 0);
	st_hash_destruct(adjacencies);
}

static void addPseudoPseudoAdjacencies(st_Hash *adjacenciesHash, Net *net) {
	/*
	 * Adds pseudo-pseudo adjacencies to ends that have no valid adjacency pair in the hash. Added adjacencies
	 * are put in the hash and the list.
	 */
	End *end1;
	End *end2 = NULL;
	Net_EndIterator *endIterator = net_getEndIterator(net);
	while((end1 = net_getNextEnd(endIterator)) != NULL) {
		if(end_isBlockEnd(end1) || end_isAttached(end1)) { //not interested in free stubs.
			assert(end_getPositiveOrientation(end1) == end1); //should be positive orientation for this to work.
			if(st_hash_search(adjacenciesHash, end1) == NULL) {
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
	}
	assert(end2 == NULL); //there must be an even number of pairs.
}

static void extendComponent(End *end, struct List *component, st_Hash *componentsHash, st_Hash *adjacencies) {
	/*
	 * Sub-function of get connected components, extends the component.
	 */
	assert(end_getOrientation(end));
	if(st_hash_search(componentsHash, end) == NULL) {
		listAppend(component, end);
		st_hash_insert(componentsHash, end, component);
		AdjacencyPair *adjacencyPair = st_hash_search(adjacencies, end);
		if(adjacencyPair != NULL) {
			extendComponent(adjacencyPair_getOtherEnd(adjacencyPair, end), component, componentsHash, adjacencies);
		}
		if(end_isBlockEnd(end)) {
			extendComponent(end_getOtherBlockEnd(end), component, componentsHash, adjacencies);
		}
	}
}

struct List *getConnectedComponents(st_Hash *adjacencies, Net *net) {
	/*
	 * Gets a list of connected components of ends linked by adjacencies and blocks.
	 */
	struct List *components = constructEmptyList(0, (void (*)(void *))destructList);
	st_Hash *componentsHash = st_hash_construct();
	End *end;
	Net_EndIterator *endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) { //iterates over the positive oriented ends.
		if(end_isBlockEnd(end) || end_isAttached(end)) { //ignore free stubs
			struct List *component = st_hash_search(componentsHash, end);
			if(component == NULL) {
				component = constructEmptyList(0, NULL);
				listAppend(components, component);
				extendComponent(end, component, componentsHash, adjacencies);
			}
		}
	}
	net_destructEndIterator(endIterator);
	st_hash_destruct(componentsHash);
	return components;
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

void switchAdjacencies(AdjacencyPair *adjacencyPair1, AdjacencyPair *adjacencyPair2, st_Hash *adjacencies,
		AdjacencyPair **adjacencyPair3, AdjacencyPair **adjacencyPair4) {
	/*
	 * Removes two adjacencies and destroys them and creates two more (in random configuration of the two possible configs).
	 */
	*adjacencyPair3 = adjacencyPair_construct(adjacencyPair_getEnd1(adjacencyPair1), adjacencyPair_getEnd1(adjacencyPair2));
	*adjacencyPair4 = adjacencyPair_construct(adjacencyPair_getEnd2(adjacencyPair1), adjacencyPair_getEnd2(adjacencyPair2));
	removeAndDestructAdjacency(adjacencies, adjacencyPair1);
	removeAndDestructAdjacency(adjacencies, adjacencyPair2);
	addAdjacencyPairToHash(adjacencies, *adjacencyPair3);
	addAdjacencyPairToHash(adjacencies, *adjacencyPair4);
}

void splitIntoContigsAndCycles(struct List *components, struct List **contigs, struct List **cycles) {
	/*
	 * Divides the list of components into rings and components terminated by a pair of stubs.
	 */
	End *end1, *end2;
	*cycles = constructEmptyList(0, NULL);
	*contigs = constructEmptyList(0, NULL);
	struct List *component;
	while(components->length > 0) {
		component = components->list[--components->length];
		getAttachedStubEndsInComponent(component, &end1, &end2);
		if(end1 != NULL) {
			assert(end2 != NULL);
			listAppend(*contigs, component);
		}
		else {
			assert(end2 == NULL);
			listAppend(*cycles, component);
		}
	}
}

AdjacencyPair *getAdjacencyPair(struct List *component, st_Hash *adjacencies) {
	/*
	 * Gets an adjacency pair from a component.
	 */
	assert(component->length > 0);
	AdjacencyPair *adjacencyPair = st_hash_search(adjacencies, component->list[0]);
	assert(adjacencyPair != NULL);
	return adjacencyPair;
}

void mergeCycles(st_Hash *adjacencies, Net *net) {
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

static st_Hash *getCorrectEndPairing(Reference *reference) {
	/*
	 * Constructs a hash set of adjacency pairs, each pair being the 5 and 3 prime
	 * ends of a pseudo chromosome in the reference.
	 */
	st_Hash *correctEndPairing = st_hash_construct();
	PseudoChromosome *pseudoChromosome;
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator = reference_getPseudoChromosomeIterator(reference);
	while((pseudoChromosome = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		AdjacencyPair *adjacencyPair = adjacencyPair_construct(pseudoChromosome_get5End(pseudoChromosome), pseudoChromosome_get3End(pseudoChromosome));
		addAdjacencyPairToHash(correctEndPairing, adjacencyPair);
	}
	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
	return correctEndPairing;
}

void breakApartMisPairedContigs(st_Hash *correctEndPairing,
		st_Hash *adjacencies, Net *net) {
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
		AdjacencyPair *adjacencyPair1 = st_hash_search(correctEndPairing, end1); //defensive
		AdjacencyPair *adjacencyPair2 = st_hash_search(correctEndPairing, end2);
		if(adjacencyPair1 != adjacencyPair2) {
			//Break an adjacency
			removeAndDestructAdjacency(adjacencies, getAdjacencyPair(component, adjacencies));
		}
	}
	destructList(components);
}

End *getFreeEnd(struct List *component, st_Hash *adjacencies) {
	/*
	 * Gets the end in the component with no adjacency.
	 */
	int32_t i;
	End *end = NULL;
	for(i=0; i<component->length; i++) {
		End *end2 = component->list[i];
		if(st_hash_search(adjacencies, end2) == NULL) {
			assert(end == NULL);
			end = end2;
		}
	}
	assert(end != NULL);
	return end;
}

void pairBrokenContigs(st_Hash *correctEndPairing,
		st_Hash *adjacencies, Net *net) {
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
			AdjacencyPair *adjacencyPair1 = st_hash_search(correctEndPairing, end1);
			for(j=0; j<i; j++) {
				End *end3, *end4;
				getAttachedStubEndsInComponent(components->list[j], &end3, &end4);
				assert(end3 != NULL);
				if(end4 == NULL) {
					AdjacencyPair *adjacencyPair2 = st_hash_search(correctEndPairing, end3);
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

void correctAttachedStubEndPairing(st_Hash *adjacencies, Net *net, Reference *reference) {
	/*
	 * If a-b and c-d are two components respectively containing stub ends a and b and c and d and the correct pairing
	 * is a-c b-d then two adjacencies are broken and two are created to make this happen.
	 * Destructs the list of components as we go.
	 */
	//Get the correct end pairings of the pseudo-telomeres in a hash of pseudo adjacencies.
	st_Hash *correctEndPairing = getCorrectEndPairing(reference);
	breakApartMisPairedContigs(correctEndPairing, adjacencies, net);
	pairBrokenContigs(correctEndPairing, adjacencies, net);
	adjacenciesHash_cleanUp(net, correctEndPairing);
}

static void fillInPseudoAdjacenciesP(End *end1, PseudoChromosome *pseudoChromosome, st_Hash *adjacencies) {
	/*
	 * Traverses the connected component of the contig constructing the pseudo-adjacencies in the pseudo-chromosome.
	 */
	AdjacencyPair *adjacencyPair = st_hash_search(adjacencies, end1);
	assert(adjacencyPair != NULL);
	//construct the new pseudo-adjacency!
	End *end2 = adjacencyPair_getOtherEnd(adjacencyPair, end1);
	Group *group1 = end_getGroup(end1);
	Group *group2 = end_getGroup(end2);
	if(group1 != group2) {
		Block *pseudoBlock = block_construct(1, end_getNet(end1));
		End *end3 = block_get5End(pseudoBlock), *end4 = block_get3End(pseudoBlock);
		if(group_isLink(group1)) {
			link_split(group_getLink(group1));
		}
		if(group_isLink(group2)) {
			link_split(group_getLink(group2));
		}
		end_setGroup(end3, group1);
		end_setGroup(end4, group2);
		//Push the new ends into the children
		pushEndIntoChildNets(end3);
		pushEndIntoChildNets(end4);
		pseudoAdjacency_construct(end1, end3, pseudoChromosome);
		pseudoAdjacency_construct(end4, end2, pseudoChromosome);
		//Push block into children
	}
	else {
		pseudoAdjacency_construct(end1, end2, pseudoChromosome);
	}
	if(pseudoChromosome_get3End(pseudoChromosome) != end2) {
		assert(end_isBlockEnd(end2));
		fillInPseudoAdjacenciesP(end_getOtherBlockEnd(end2), pseudoChromosome, adjacencies);
	}
}

static void fillInPseudoAdjacencies(Net *net, Reference *reference, st_Hash *adjacencies) {
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
	st_Hash *adjacenciesHash = choosePairing(adjacencies, net);
	addPseudoPseudoAdjacencies(adjacenciesHash, net);
	mergeCycles(adjacenciesHash, net);
	correctAttachedStubEndPairing(adjacenciesHash, net, reference);
	fillInPseudoAdjacencies(net, reference, adjacenciesHash);
	adjacenciesHash_cleanUp(net, adjacenciesHash);
}
