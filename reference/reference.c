/*
 * reference.c
 *
 *  Created on: 1 Apr 2010
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "cactus.h"

#include "balanceTangles.h"
#include "pseudoChromosomeBuilder.h"
#include "pseudoAdjacencyBuilder.h"


void addReferenceToFlower(Flower *flower) {
#ifdef BEN_DEBUG
    flower_check(flower);
#endif
    //Ensure the cactus graph is balanced.
    balanceTangles(flower);
    return;

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
         * upon the set of attached stubs.
         * We do this in the order of the sequences that were passed to us.
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

Group *getSpareGroup(Flower *flower) {
    //First try and find group with an odd number of ends..
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group, *group2 = NULL;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        int32_t i = group_getFreeStubEndNumber(group);
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

static void linkZeroSizeGroups(Flower *flower) {
    *
     * Adds block ends into groups with zero non-free stub ends,
     * so that they will be properly included in the traversal.
     *
    Group *group;
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        if (group_getFreeStubEndNumber(group) == group_getEndNumber(group)) { //We need to add some new blocks to the problem to link it into the problem..
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
    *
     * Creates a top level attached stub and propogates it into the children.
     *
    End *end = end_construct(1, flower);
    end_setGroup(end, getSpareGroup(flower));
    pushEndIntoChildFlowers(end);
    return end;
}

*/
