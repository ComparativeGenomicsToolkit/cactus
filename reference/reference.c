/*
 * reference.c
 *
 *  Created on: 1 Apr 2010
 *      Author: benedictpaten
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

static struct List *getAttachedEnds(Net *net) {
	/*
	 * Get the top level attached ends.
	 */
	struct List *list = constructEmptyList(0, NULL);
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(end_isAttached(end) && end_getInstanceNumber(end) > 0) {
			listAppend(list, end);
		}
	}
	net_destructEndIterator(endIterator);
	return list;
}

static void makePseudoChromosomesFromPairs(struct List *ends, Reference *reference) {
	/*
	 * Make pseudo chromosomes from list of paired ends
	 */
	assert(ends->length >= 2); //the reconstruction must contain 2
	//or more attached top level ends to build a reference genome currently.
	assert(ends->length % 2 == 0); //we must have an even number of ends.
	int32_t i;
	for(i=0; i<=ends->length; i+=2) {
		End *end1 = ends->list[i];
		End *end2 = ends->list[i+1];
		pseudoChromosome_construct(reference, end1, end2);
	}
}

static void makePseudoChromosomes(Net *net, Reference *reference, int (*cmpFn)(End **, End **)) {
	/*
	 * Uses the above functions to construct an set of pairs of ends and then construct the pseudo-chromsomes,
	 * but without any pseudo-adjacencies.
	 */
	struct List *ends = getAttachedEnds(net);
	/* Sort the ends and so construct pairing of attached ends to construct pseudo chromosomes*/
	qsort(ends->list, ends->length, sizeof(void *),
				(int (*)(const void *v, const void *))cmpFn);
	makePseudoChromosomesFromPairs(ends, reference);
	destructList(ends);
}

static int makeTopLevelPseudoChromosomes_cmpEnds(End **end1, End **end2) {
	/*
	 * Sorts the attached ends according to the sequences they are connected to.
	 */
	assert(end_getInstanceNumber(*end1) == 1); //this is the top level
	assert(end_getInstanceNumber(*end2) == 1);
	Cap *cap1 = end_getFirst(*end1);
	Cap *cap2 = end_getFirst(*end2);
	int32_t i = netMisc_nameCompare(sequence_getName(cap_getSequence(cap1)), sequence_getName(cap_getSequence(cap2)));
	if(i == 0) {
		assert(cap_getSide(cap1) != cap_getSide(cap2));
		return cap_getSide(cap1) ? 1 : -1;
	}
	else {
		return i;
	}
}

void makeTopLevelPseudoChromosomes(Net *net, Reference *reference) {
	makePseudoChromosomes(net, reference, makeTopLevelPseudoChromosomes_cmpEnds);
}

static int makeIntermediateLevelPseudoChromosomes_cmpEnds(End **end1, End **end2) {
	/*
	 * Sorts the attached ends according to the pairing in the higher level reference.
	 */
	Net *net = end_getNet(*end1);
	Group *parentGroup = net_getParentGroup(net);
	assert(parentGroup != NULL);
	assert(net_getReferenceNumber(net) == 1);
	Reference *reference = net_getFirstReference(net);
	PseudoAdjacency *pseudoAdjacency1 = reference_getPseudoAjacencyByEnd(reference, net_getEnd(net, end_getName(*end1)));
	PseudoAdjacency *pseudoAdjacency2 = reference_getPseudoAjacencyByEnd(reference, net_getEnd(net, end_getName(*end2)));
	PseudoChromosome *pseudoChromosome1 = pseudoAdjacency_getPseudoChromosome(pseudoAdjacency1);
	PseudoChromosome *pseudoChromosome2 = pseudoAdjacency_getPseudoChromosome(pseudoAdjacency2);

	//Sort first by pseudo chromosome
	int32_t i = netMisc_nameCompare(pseudoChromosome_getName(pseudoChromosome1), pseudoChromosome_getName(pseudoChromosome2));
	if(i != 0) {
		return i;
	}
	//Then by pseudo-adjacency
	i = netMisc_nameCompare(pseudoAdjacency_getName(pseudoAdjacency1), pseudoAdjacency_getName(pseudoAdjacency2));
	if(i != 0) {
		return i;
	}
	//if in the same pseudo adjacency then
	//return the proper pair according to ordering of ends in pair.
	return end_getName(pseudoAdjacency_get5End(pseudoAdjacency1)) == end_getName(*end1) ? 1 : -1;
}

void makeIntermediateLevelPseudoChromosomes(Net *net, Reference *reference) {
	makePseudoChromosomes(net, reference, makeIntermediateLevelPseudoChromosomes_cmpEnds);
}

void mergeGroupsLinkedByPseudoAdjacencies(Net *net, Reference *reference) {
	PseudoChromosome *pseudoChromosome;
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator = reference_getPseudoChromosomeIterator(reference);
	while((pseudoChromosome = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		//iterate through all the pseudo chromosomes
		PseudoAdjacency *pseudoAdjacency;
		PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator = pseudoChromosome_getPseudoAdjacencyIterator(pseudoChromosome);
		//iterate through the pseudo adjacencies in the pseudo chromosomes.
		while((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(pseudoAdjacencyIterator)) != NULL) {
			End *_5End = pseudoAdjacency_get5End(pseudoAdjacency);
			End *_3End = pseudoAdjacency_get5End(pseudoAdjacency);
			//if the groups are distinct then we call the merge function..
			if(end_getGroup(_5End) != end_getGroup(_3End)) {
				group_mergeGroups(end_getGroup(_5End), end_getGroup(_3End));
			}
		}
		pseudoChromosome_destructPseudoChromosomeIterator(pseudoAdjacencyIterator);
	}

	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
}

void addReferenceToNet(Net *net) {
	//Function will currently only work if no reference has already been added.
	assert(net_getReferenceNumber(net) == 0);
	Reference *reference = reference_construct(net);

	if(net_getParentGroup(net) == NULL) {
		/*
		 * If this is the top level net then we will create the pseudo-chromosomes based
		 * upon the set of attached stubs.. (we will throw an error ?! if we don't have at least one pair of
		 * attached stubs). We do this in the order of the sequences that were passed to us.
		 */
		makeTopLevelPseudoChromosomes(net, reference);
	}
	else {
		/*
		 * Else this is not the top level net, and we must locate the pseudo-adjacencies in the parent
		 * reference, to establish the ends of pseudo-chromosomes.
		 * We do this in the order of the parent reference's pseudo-adjacencies.
		 */
		makeIntermediateLevelPseudoChromosomes(net, reference);
	}
	/*
	 * Having defined the ordered pseudo chromosomes, we fill in the pseudo adjacencies.
	 * For each pseudo-chromosome..
	 */
	makePseudoAdjacencies(net, reference);
	/*
	 * Now merge groups linked by novel pseudo-adjacencies.
	 */
	mergeGroupsLinkedByPseudoAdjacencies(net, reference);
}

