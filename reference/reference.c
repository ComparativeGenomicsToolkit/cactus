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
#include "commonC.h"

int32_t getFreeStubEndNumber(Group *group) {
	End *end;
	Group_EndIterator *endIterator = group_getEndIterator(group);
	int32_t i = 0;
	while((end = group_getNextEnd(endIterator)) != NULL) {
		if(end_isStubEnd(end) && end_isFree(end)) {
			i++;
		}
	}
	group_destructEndIterator(endIterator);
	return i;
}

Group *getSpareGroup(Net *net) {
	//First try and find group with an odd number of ends..
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group, *group2 = NULL;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		int32_t i = getFreeStubEndNumber(group);
		int32_t j = group_getEndNumber(group) - i; //the number of block ends and attached stub ends.
		assert(j >= 0);
		if(j % 2) {
			assert(!group_isLink(group));
			net_destructGroupIterator(groupIterator);
			return group;
		}
		else if(!group_isLink(group)){
			group2 = group;
		}
	}
	net_destructGroupIterator(groupIterator);

	//Else get a group without a link..
	if(group2 != NULL) {
		assert(!group_isLink(group2));
		return group2;
	}

	//Else all groups are link groups.. so get the first one and remove the link from the chain..
	assert(net_getGroupNumber(net) > 0);
	group = net_getFirstGroup(net);
	assert(group != NULL);
	assert(group_isLink(group));
	link_split(group_getLink(group));
	assert(!group_isLink(group));
	return group;
}

void pushEndIntoChildNets(End *end) {
	assert(end_isAttached(end) || end_isBlockEnd(end));
	Group *group = end_getGroup(end);
	assert(group != NULL);
	if(!group_isLeaf(group)) {
		Net *nestedNet = group_getNestedNet(group);
		assert(net_getEnd(nestedNet, end_getName(end)) == NULL);
		End *end2 = end_copyConstruct(end, nestedNet);
		end_setGroup(end2, getSpareGroup(nestedNet));
		assert(end_getGroup(end2) != NULL);
		assert(!group_isLink(end_getGroup(end2)));
		//Now call recursively
		pushEndIntoChildNets(end2);
	}
}

static struct List *getAttachedStubEnds(Net *net) {
	/*
	 * Get the top level attached ends.
	 */
	struct List *list = constructEmptyList(0, NULL);
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(end_isAttached(end)) {
			assert(end_isStubEnd(end));
			assert(!end_isFree(end));
			listAppend(list, end);
		}
	}
	net_destructEndIterator(endIterator);
	assert(list->length % 2 == 0);
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
	for(i=0; i<ends->length; i+=2) {
		End *end1 = ends->list[i];
		End *end2 = ends->list[i+1];
		assert(end1 != NULL);
		assert(end2 != NULL);
		assert(end_isStubEnd(end1));
		assert(end_isAttached(end1));
		assert(end_isStubEnd(end2));
		assert(end_isAttached(end2));
		pseudoChromosome_construct(reference, end1, end2);
	}
}

static void linkZeroSizeGroups(Net *net) {
	/*
	 * Adds block ends into groups with zero non-free stub ends,
	 * so that they will be properly included in the traversal.
	 */
	Group *group;
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(getFreeStubEndNumber(group) == group_getEndNumber(group)) { //We need to add some new blocks to the problem to link it into the problem..
			assert(group_isTangle(group));
			//Construct two pseudo blocks so that the group gets included.
			Block *block = block_construct(1, net);
			End *_5End = block_get5End(block);
			End *_3End = block_get3End(block);
			end_setGroup(_5End, group);
			end_setGroup(_3End, group);
			pushEndIntoChildNets(_5End);
			pushEndIntoChildNets(_3End);
		}
	}
	net_destructGroupIterator(groupIterator);
}

static End *constructTopLevelAttachedStub(Net *net) {
	/*
	 * Creates a top level attached stub and propogates it into the children.
	 */
	End *end = end_construct(1, net);
	end_setGroup(end, getSpareGroup(net));
	pushEndIntoChildNets(end);
	return end;
}

static void makePseudoChromosomes(Net *net, Reference *reference, int (*cmpFn)(End **, End **)) {
	/*
	 * Uses the above functions to construct a set of pairs of ends and then construct the pseudo-chromsomes,
	 * but without any pseudo-adjacencies.
	 */
	struct List *ends = getAttachedStubEnds(net);
	if(ends->length == 0) {
		assert(net_getParentGroup(net) == NULL); //If there are no attached ends / block ends in the problem we must be in the parent net or else our code to add them recursively has gone wrong
		listAppend(ends, constructTopLevelAttachedStub(net));
		listAppend(ends, constructTopLevelAttachedStub(net));
	}
	//Now think ahead to include all child groups in the ordering.
	linkZeroSizeGroups(net); //these ends will be block ends, not attached ends, so no need to add to ends list.
	/* Sort the ends and so construct pairing of attached ends to construct pseudo chromosomes*/
	qsort(ends->list, ends->length, sizeof(void *),
				(int (*)(const void *v, const void *))cmpFn);
	makePseudoChromosomesFromPairs(ends, reference);
	destructList(ends);
}

static Cap *makeTopLevelPseudoChromosomesP(End *end) {
	/*
	 * Gets the leaf cap from a top level attached stub. There can only
	 * be one such cap per top level attached stub, currently.
	 */
	Cap *cap = NULL, *cap2;
	End_InstanceIterator *iterator = end_getInstanceIterator(end);
	while((cap2 = end_getNext(iterator)) != NULL) {
		if(cap_getChildNumber(cap2) == 0) {
			assert(cap == NULL);
			cap = cap2;
		}
	}
	end_destructInstanceIterator(iterator);
	assert(cap != NULL);
	return cap;
}

static int makeTopLevelPseudoChromosomes_cmpEnds(End **end1, End **end2) {
	/*
	 * Sorts the attached ends according to the sequences they are connected to.
	 */
	Cap *cap1 = makeTopLevelPseudoChromosomesP(*end1);
	Cap *cap2 = makeTopLevelPseudoChromosomesP(*end2);
	Sequence *sequence1 = cap_getSequence(cap1);
	Sequence *sequence2 = cap_getSequence(cap2);
	assert(sequence1 != NULL);
	assert(sequence2 != NULL);
	int32_t i = netMisc_nameCompare(sequence_getName(sequence1), sequence_getName(sequence2));
	if(i == 0) {
		assert(cap_getSide(cap1) != cap_getSide(cap2));
		return cap_getSide(cap1) ? -1 : 1; //sort 5' to 3' (i.e. with the 5 prime with a lower index to the 3 prime side)
	}
	else {
		return i;
	}
}

void makeTopLevelPseudoChromosomes(Net *net, Reference *reference) {
	makePseudoChromosomes(net, reference, makeTopLevelPseudoChromosomes_cmpEnds);
}

static st_Hash *makeIntermediateLevelPseudoChromosomes_cmpEndsP = NULL;
static Net *makeIntermediateLevelPseudoChromosomes_parentNet = NULL;

static int makeIntermediateLevelPseudoChromosomes_cmpEnds(End **end1, End **end2) {
	/*
	 * Sorts the attached ends according to the pairing in the higher level reference.
	 */
	End *end3 = net_getEnd(makeIntermediateLevelPseudoChromosomes_parentNet, end_getName(*end1));
	End *end4 = net_getEnd(makeIntermediateLevelPseudoChromosomes_parentNet, end_getName(*end2));
	assert(end3 != NULL);
	assert(end4 != NULL);
	assert(end_getOrientation(end3));
	assert(end_getOrientation(end4));
	assert(end_isAttached(end3) || end_isBlockEnd(end3));
	assert(end_isAttached(end4) || end_isBlockEnd(end4));

	PseudoAdjacency *pseudoAdjacency1 = st_hash_search(makeIntermediateLevelPseudoChromosomes_cmpEndsP, end3);
	PseudoAdjacency *pseudoAdjacency2 = st_hash_search(makeIntermediateLevelPseudoChromosomes_cmpEndsP, end4);
	assert(pseudoAdjacency1 != NULL);
	assert(pseudoAdjacency2 != NULL);
	PseudoChromosome *pseudoChromosome1 = pseudoAdjacency_getPseudoChromosome(pseudoAdjacency1);
	PseudoChromosome *pseudoChromosome2 = pseudoAdjacency_getPseudoChromosome(pseudoAdjacency2);

	//Sort first by pseudo chromosome
	int32_t i = netMisc_nameCompare(pseudoChromosome_getName(pseudoChromosome1), pseudoChromosome_getName(pseudoChromosome2));
	if(i != 0) { //we order, so the chromosomes reflect the ordering of the parent adjacencies on the circle.
		return i;
	}
	assert(pseudoChromosome1 == pseudoChromosome2);
	//Then by pseudo-adjacency
	i = netMisc_nameCompare(pseudoAdjacency_getName(pseudoAdjacency1), pseudoAdjacency_getName(pseudoAdjacency2));
	if(i != 0) {
		return i;
	}
	assert(pseudoAdjacency1 == pseudoAdjacency2);
	//if in the same pseudo adjacency then
	//return the proper pair according to ordering of ends in pair.
	assert(end3 == pseudoAdjacency_get5End(pseudoAdjacency1) || end3 == pseudoAdjacency_get3End(pseudoAdjacency1));
	assert(end4 == pseudoAdjacency_get5End(pseudoAdjacency1) || end4 == pseudoAdjacency_get3End(pseudoAdjacency1));
	return pseudoAdjacency_get5End(pseudoAdjacency1) == end3 ? -1 : 1;
}

void makeIntermediateLevelPseudoChromosomes(Net *net, Reference *reference) {
	Group *parentGroup = net_getParentGroup(net);
	assert(parentGroup != NULL);
	Net *parentNet = group_getNet(parentGroup);
	Reference *parentReference = net_getReference(parentNet);
	assert(parentReference != NULL);

	makeIntermediateLevelPseudoChromosomes_parentNet = parentNet;
	makeIntermediateLevelPseudoChromosomes_cmpEndsP = reference_getEndToPseudoAdjacencyHash(parentReference);

	makePseudoChromosomes(net, reference, makeIntermediateLevelPseudoChromosomes_cmpEnds);

	st_hash_destruct(makeIntermediateLevelPseudoChromosomes_cmpEndsP);
}

void addReferenceToNet(Net *net) {
	Reference *reference = net_getReference(net);
	if(reference != NULL) {
		return; //we've already built it, so no need to do it again!
	}
	reference = reference_construct(net);

	if(net_getGroupNumber(net) == 0) { //In this case we have nothing to add, and no point in continuing.
		assert(net_getEndNumber(net) == 0);
		return;
	}

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
	 * Now check the reference created for goodness.
	 */
#ifdef BEN_DEBUG
	reference_check(reference);
#endif
}

