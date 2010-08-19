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
#include "balanceTangles.h"

int32_t getFreeStubEndNumber(Group *group) {
    End *end;
    Group_EndIterator *endIterator = group_getEndIterator(group);
    int32_t i = 0;
    while ((end = group_getNextEnd(endIterator)) != NULL) {
        if (end_isStubEnd(end) && end_isFree(end)) {
            i++;
        }
    }
    group_destructEndIterator(endIterator);
    return i;
}

Group *getSpareGroup(Flower *flower) {
    //First try and find group with an odd number of ends..
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group, *group2 = NULL;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        int32_t i = getFreeStubEndNumber(group);
        int32_t j = group_getEndNumber(group) - i; //the number of block ends and attached stub ends.
        assert(j >= 0);
        if (j % 2) {
            assert(!group_isLink(group));
            flower_destructGroupIterator(groupIterator);
            return group;
        } else if (!group_isLink(group)) {
            group2 = group;
        }
    }
    flower_destructGroupIterator(groupIterator);

    //Else get a group without a link..
    if (group2 != NULL) {
        assert(!group_isLink(group2));
        return group2;
    }

    //Else all groups are link groups.. so get the first one and remove the link from the chain..
    assert(flower_getGroupNumber(flower) > 0);
    group = flower_getFirstGroup(flower);
    assert(group != NULL);
    assert(group_isLink(group));
    link_split(group_getLink(group));
    assert(!group_isLink(group));
    return group;
}

void pushEndIntoChildFlowers(End *end) {
    assert(end_isAttached(end) || end_isBlockEnd(end));
    Group *group = end_getGroup(end);
    assert(group != NULL);
    if (!group_isLeaf(group)) {
        Flower *nestedFlower = group_getNestedFlower(group);
        assert(flower_getEnd(nestedFlower, end_getName(end)) == NULL);
        End *end2 = end_copyConstruct(end, nestedFlower);
        end_setGroup(end2, getSpareGroup(nestedFlower));
        assert(end_getGroup(end2) != NULL);
        assert(!group_isLink(end_getGroup(end2)));
        //Now call recursively
        pushEndIntoChildFlowers(end2);
    }
}

static struct List *getAttachedStubEnds(Flower *flower) {
    /*
     * Get the top level attached ends.
     */
    struct List *list = constructEmptyList(0, NULL);
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        if (end_isAttached(end)) {
            assert(end_isStubEnd(end));
            assert(!end_isFree(end));
            listAppend(list, end);
        }
    }
    flower_destructEndIterator(endIterator);
    assert(list->length % 2 == 0);
    return list;
}

static void makePseudoChromosomesFromPairs(struct List *ends,
        Reference *reference) {
    /*
     * Make pseudo chromosomes from list of paired ends
     */
    assert(ends->length >= 2); //the reconstruction must contain 2
    //or more attached top level ends to build a reference genome currently.
    assert(ends->length % 2 == 0); //we must have an even number of ends.
    int32_t i;
    for (i = 0; i < ends->length; i += 2) {
        End *end1 = ends->list[i];
        End *end2 = ends->list[i + 1];
        assert(end1 != NULL);
        assert(end2 != NULL);
        assert(end_isStubEnd(end1));
        assert(end_isAttached(end1));
        assert(end_isStubEnd(end2));
        assert(end_isAttached(end2));
        pseudoChromosome_construct(reference, end1, end2);
    }
}

static void linkZeroSizeGroups(Flower *flower) {
    /*
     * Adds block ends into groups with zero non-free stub ends,
     * so that they will be properly included in the traversal.
     */
    Group *group;
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        if (getFreeStubEndNumber(group) == group_getEndNumber(group)) { //We need to add some new blocks to the problem to link it into the problem..
            assert(group_isTangle(group));
            //Construct two pseudo blocks so that the group gets included.
            Block *block = block_construct(1, flower);
            End *_5End = block_get5End(block);
            End *_3End = block_get3End(block);
            end_setGroup(_5End, group);
            end_setGroup(_3End, group);
            pushEndIntoChildFlowers(_5End);
            pushEndIntoChildFlowers(_3End);
        }
    }
    flower_destructGroupIterator(groupIterator);
}

static End *constructTopLevelAttachedStub(Flower *flower) {
    /*
     * Creates a top level attached stub and propogates it into the children.
     */
    End *end = end_construct(1, flower);
    end_setGroup(end, getSpareGroup(flower));
    pushEndIntoChildFlowers(end);
    return end;
}

static void makePseudoChromosomes(Flower *flower, Reference *reference,
        int(*cmpFn)(End **, End **)) {
    /*
     * Uses the above functions to construct a set of pairs of ends and then construct the pseudo-chromsomes,
     * but without any pseudo-adjacencies.
     */
    struct List *ends = getAttachedStubEnds(flower);
    if (ends->length == 0) {
        assert(flower_getParentGroup(flower) == NULL); //If there are no attached ends / block ends in the problem we must be in the parent flower or else our code to add them recursively has gone wrong
        listAppend(ends, constructTopLevelAttachedStub(flower));
        listAppend(ends, constructTopLevelAttachedStub(flower));
    }
    //Now think ahead to include all child groups in the ordering.
    linkZeroSizeGroups(flower); //these ends will be block ends, not attached ends, so no need to add to ends list.
    /* Sort the ends and so construct pairing of attached ends to construct pseudo chromosomes*/
    qsort(ends->list, ends->length, sizeof(void *), (int(*)(const void *v,
            const void *)) cmpFn);
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
    while ((cap2 = end_getNext(iterator)) != NULL) {
        if (cap_getChildNumber(cap2) == 0) {
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
    int32_t i = cactusMisc_nameCompare(sequence_getName(sequence1),
            sequence_getName(sequence2));
    if (i == 0) {
        assert(cap_getSide(cap1) != cap_getSide(cap2));
        return cap_getSide(cap1) ? -1 : 1; //sort 5' to 3' (i.e. with the 5 prime with a lower index to the 3 prime side)
    } else {
        return i;
    }
}

void makeTopLevelPseudoChromosomes(Flower *flower, Reference *reference) {
    makePseudoChromosomes(flower, reference,
            makeTopLevelPseudoChromosomes_cmpEnds);
}

static stHash *makeIntermediateLevelPseudoChromosomes_cmpEndsP = NULL;
static Flower *makeIntermediateLevelPseudoChromosomes_parentFlower = NULL;

static int makeIntermediateLevelPseudoChromosomes_cmpEnds(End **end1,
        End **end2) {
    /*
     * Sorts the attached ends according to the pairing in the higher level reference.
     */
    End *end3 = flower_getEnd(
            makeIntermediateLevelPseudoChromosomes_parentFlower, end_getName(
                    *end1));
    End *end4 = flower_getEnd(
            makeIntermediateLevelPseudoChromosomes_parentFlower, end_getName(
                    *end2));
    assert(end3 != NULL);
    assert(end4 != NULL);
    assert(end_getOrientation(end3));
    assert(end_getOrientation(end4));
    assert(end_isAttached(end3) || end_isBlockEnd(end3));
    assert(end_isAttached(end4) || end_isBlockEnd(end4));

    PseudoAdjacency *pseudoAdjacency1 = stHash_search(
            makeIntermediateLevelPseudoChromosomes_cmpEndsP, end3);
    PseudoAdjacency *pseudoAdjacency2 = stHash_search(
            makeIntermediateLevelPseudoChromosomes_cmpEndsP, end4);
    assert(pseudoAdjacency1 != NULL);
    assert(pseudoAdjacency2 != NULL);
    PseudoChromosome *pseudoChromosome1 = pseudoAdjacency_getPseudoChromosome(
            pseudoAdjacency1);
    PseudoChromosome *pseudoChromosome2 = pseudoAdjacency_getPseudoChromosome(
            pseudoAdjacency2);

    //Sort first by pseudo chromosome
    int32_t i = cactusMisc_nameCompare(pseudoChromosome_getName(
            pseudoChromosome1), pseudoChromosome_getName(pseudoChromosome2));
    if (i != 0) { //we order, so the chromosomes reflect the ordering of the parent adjacencies on the circle.
        return i;
    }
    assert(pseudoChromosome1 == pseudoChromosome2);
    //Then by pseudo-adjacency
    i = cactusMisc_nameCompare(pseudoAdjacency_getName(pseudoAdjacency1),
            pseudoAdjacency_getName(pseudoAdjacency2));
    if (i != 0) {
        return i;
    }
    assert(pseudoAdjacency1 == pseudoAdjacency2);
    //if in the same pseudo adjacency then
    //return the proper pair according to ordering of ends in pair.
    assert(end3 == pseudoAdjacency_get5End(pseudoAdjacency1) || end3 == pseudoAdjacency_get3End(pseudoAdjacency1));
    assert(end4 == pseudoAdjacency_get5End(pseudoAdjacency1) || end4 == pseudoAdjacency_get3End(pseudoAdjacency1));
    return pseudoAdjacency_get5End(pseudoAdjacency1) == end3 ? -1 : 1;
}

void makeIntermediateLevelPseudoChromosomes(Flower *flower,
        Reference *reference) {
    Group *parentGroup = flower_getParentGroup(flower);
    assert(parentGroup != NULL);
    Flower *parentFlower = group_getFlower(parentGroup);
    Reference *parentReference = flower_getReference(parentFlower);
    assert(parentReference != NULL);

    makeIntermediateLevelPseudoChromosomes_parentFlower = parentFlower;
    makeIntermediateLevelPseudoChromosomes_cmpEndsP
            = reference_getEndToPseudoAdjacencyHash(parentReference);

    makePseudoChromosomes(flower, reference,
            makeIntermediateLevelPseudoChromosomes_cmpEnds);

    stHash_destruct(makeIntermediateLevelPseudoChromosomes_cmpEndsP);
}

void addReferenceToFlower(Flower *flower) {
#ifdef BEN_DEBUG
    flower_check(flower);
#endif
    //Ensure the cactus graph is balanced.
    balanceTangles(flower);

    Reference *reference = flower_getReference(flower);
    if (reference != NULL) {
        return; //we've already built it, so no need to do it again!
    }
    reference = reference_construct(flower);

    if (flower_getGroupNumber(flower) == 0) { //In this case we have nothing to add, and no point in continuing.
        assert(flower_getEndNumber(flower) == 0);
        return;
    }

    if (flower_getParentGroup(flower) == NULL) {
        /*
         * If this is the top level flower then we will create the pseudo-chromosomes based
         * upon the set of attached stubs.. (we will throw an error ?! if we don't have at least one pair of
         * attached stubs). We do this in the order of the sequences that were passed to us.
         */
        makeTopLevelPseudoChromosomes(flower, reference);
    } else {
        /*
         * Else this is not the top level flower, and we must locate the pseudo-adjacencies in the parent
         * reference, to establish the ends of pseudo-chromosomes.
         * We do this in the order of the parent reference's pseudo-adjacencies.
         */
        makeIntermediateLevelPseudoChromosomes(flower, reference);
    }
    /*
     * Having defined the ordered pseudo chromosomes, we fill in the pseudo adjacencies.
     * For each pseudo-chromosome..
     */
    makePseudoAdjacencies(flower, reference);

    /*
     * Now check the reference created for goodness.
     */
#ifdef BEN_DEBUG
    reference_check(reference);
    flower_check(flower);
#endif
}

/*
void addReferenceToFlower2(Flower *flower) {
#ifdef BEN_DEBUG
    flower_check(flower);
#endif

    Reference *reference = flower_getReference(flower);
    if (reference != NULL) {
        return; //we've already built it, so no need to do it again!
    }
    reference = reference_construct(flower);

    if (flower_getGroupNumber(flower) == 0) { //In this case we have nothing to add, and no point in continuing.
        assert(flower_getEndNumber(flower) == 0);
        return;
    }

    if (flower_isTerminal(flower)) {
        //We arbitrarily pair the ends.
        assert(flower_getAttachedStubEndNumber(flower) % 2 == 0);
        assert(flower_getGroupNumber(flower) == 1);
        End *end;
        Flower_EndIterator *endIt = flower_getEndIterator(flower);
        stList *list = stList_construct();
        while ((end = flower_getNextEnd(endIt)) != NULL) {
            assert(end_isStubEnd(end));
            if (end_isAttached(end)) {
                stList_append(list, end);
            }
        }
        assert(stList_length(list) > 0);
        assert(stList_length(list) % 2 == 0);
        for (int32_t i = 0; i < stList_length(list); i += 2) {
            pseudoChromosome_construct(reference, stList_get(list, i), stList_get(list, i + 1));
        }
        stList_destruct(list);
        flower_destructEndIterator(endIt);
    } else {
        //We must construct the pseudo chromosomes according to the nested nets.
        End *end;
        Flower_EndIterator *endIt = flower_getEndIterator(flower);
        while ((end = flower_getNextEnd(endIt)) != NULL) {
            if (end_isAttached(end) && end_isStubEnd(end)
                    && end_getPseudoAdjacency(end) == NULL) {
                stList *list = stList_construct();
                do {
                    Group *group = end_getGroup(end);
                    assert(!group_isLeaf(group));
                    Flower *nestedFlower = group_getNestedFlower(group);
                    End *nestedEnd = flower_getEnd(nestedFlower, end_getName(
                            end));
                    assert(nestedEnd != NULL);
                    PseudoAdjacency *nestedPseudoAdjacency =
                            end_getPseudoAdjacency(nestedEnd);
                    PseudoChromosome
                            *nestedPseudoChromosome =
                                    pseudoAdjacency_getPseudoChromosome(
                                            nestedPseudoAdjacency);
                    assert(pseudoChromosome_get5End(nestedPseudoChromosome) == nestedEnd || pseudoChromosome_get3End(nestedPseudoChromosome) == nestedEnd);
                    End *otherNestedEnd =
                            pseudoChromosome_get5End(nestedPseudoChromosome)
                                    == nestedEnd ? pseudoChromosome_get3End(
                                    nestedPseudoChromosome)
                                    : pseudoChromosome_get5End(
                                            nestedPseudoChromosome);
                    End *otherEnd = flower_getEnd(flower, end_getName(
                            otherNestedEnd));
                    assert(end_getPseudoAdjacency(otherEnd) == NULL);
                    assert(otherEnd != NULL);
                    stList_append(list, end);
                    stList_append(list, otherEnd);
                    end = end_getOtherBlockEnd(otherEnd);
                } while (end != NULL);
                PseudoChromosome *pseudoChromosome =
                        pseudoChromosome_construct(reference, stList_get(list,
                                0), stList_peek(list));
                for (int32_t i = 0; i < stList_length(list); i += 2) {
                    pseudoAdjacency_construct(stList_get(list, i), stList_get(
                            list, i + 1), pseudoChromosome);
                }
                stList_destruct(list);
            }
        }
        flower_destructEndIterator(endIt);
    }
}
*/
