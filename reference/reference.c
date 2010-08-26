/*
 * reference.c
 *
 *  Created on: 1 Apr 2010
 *      Author: benedictpaten
 *
 * This source file contains the top level algorithms for building the reference.
 */

#include "cactus.h"
#include "sonLib.h"
#include "adjacencyPairs.h"
#include "adjacencyPairsHash.h"
#include "mergeCycles.h"
#include "hyperChains.h"
#include "balanceTangles.h"
#include "pressedFlowers.h"

static void makeTerminalReference(Flower *flower, stHash *adjacenciesHash) {
    /*
     * Builds the reference for a terminal flower, given the adjacencies hash.
     */
#ifdef BEN_DEBUG
    flower_check(flower);
#endif
    assert(flower_isTerminal(flower));
    /*
     * For atomicity, if we find we have already built a reference we rebuild it..
     */
    if(flower_getReference(flower) != NULL) {
        reference_destruct(flower_getReference(flower));
    }

    Reference *reference = reference_construct(flower);
    End *end;
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    while((end = flower_getNextEnd(endIt)) != NULL) {
        assert(!end_isBlockEnd(end));
        if(end_isAttached(end)) {
            AdjacencyPair *adjacencyPair = stHash_search(adjacenciesHash, end);
            assert(adjacencyPair != NULL);
            if(adjacencyPair_getEnd1(adjacencyPair) == end) { //Just build the pseudo chromosome once..
                End *end2 = adjacencyPair_getEnd2(adjacencyPair);
                PseudoChromosome *pseudoChromosome = pseudoChromosome_construct(reference,
                        end, end2);
                //Now complete the pseudo chromosome by adding in the pseudo adjacency.
                pseudoAdjacency_construct(end, end2, pseudoChromosome);
            }
        }
    }
    flower_destructEndIterator(endIt);

#ifdef BEN_DEBUG
    //Check we have paired all the attached stub ends in the terminal flower.
    assert(reference_getPseudoChromosomeNumber(reference)*2 == flower_getAttachedStubEndNumber(flower));
    flower_check(flower);
#endif
}

static void balanceTanglesRecursively(Flower *flower) {
    /*
     * Balance every tangle descendant (see the pressed flowers method)
     * of the flower, thus allowing us to build a set of adjacency pairs
     * for all the terminal ends.
     */
    if(!flower_isTerminal(flower)) {
        balanceTangles(flower);
        Group *group;
        Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
        while((group = flower_getNextGroup(groupIt)) != NULL) {
            assert(!group_isLeaf(group));
            if(group_isTangle(group)) {
                balanceTanglesRecursively(group_getNestedFlower(group));
            }
        }
        flower_destructGroupIterator(groupIt);
    }
}

void constructReference_topDownPhase(Flower *flower) {
#ifdef BEN_DEBUG
    flower_check(flower);
#endif
    if(flower_getParentGroup(flower) != NULL && group_isTangle(flower_getParentGroup(flower))) {
        /*
         * In this case we've must have already built the child tangle terminal nets.
         */
        return;
    }
    /*
     * Balance the tangle descendants of the flower and the flower itself.
     */
    balanceTanglesRecursively(flower);
    /*
     * Construct the initial pairing of the terminal ends.
     */
    stHash *adjacenciesHash = adjacencyHash_constructInitialPairs(flower);
    /*
     * Find the set of 'hyper chains' linking the ends in the terminal problems.
     */
    stHash *hyperChains = constructHyperChains(flower);
    /*
     * Get rid of cyclic components from our initial pairing.
     */
    mergeCycles(adjacenciesHash, hyperChains);
    /*
     * Create the references in for the terminal flowers in the set of
     * pressed flowers for the flower.
     */
    stList *pressedFlowers = getListOfPressedFlowers(flower);
    while(stList_length(pressedFlowers) > 0) {
        makeTerminalReference(stList_pop(pressedFlowers), adjacenciesHash);
    }
    /*
     * Finally cleanup.
     */
    stList_destruct(pressedFlowers);
    adjacencyHash_destruct(adjacenciesHash, flower);
    stHash_destruct(hyperChains);
}

void constructReference_bottomUpPhase(Flower *flower) {
#ifdef BEN_DEBUG
    flower_check(flower);
#endif

    if (flower_isTerminal(flower)) {
        /*
         * We build the reference for the terminal problems in the top down phase,
         * so no need to do anything.
         */
        assert(flower_getReference(flower) != NULL);
    } else {
        /*
         * For atomicity check if a flower does not already have a reference
         * and destory it if it does.
         */
        if(flower_getReference(flower) != NULL) {
            reference_destruct(flower_getReference(flower));
        }

        Reference *reference = reference_construct(flower);
        /*
         * We construct the pseudo chromosomes according to the nested nets.
         */
        End *end;
        Flower_EndIterator *endIt = flower_getEndIterator(flower);
        while ((end = flower_getNextEnd(endIt)) != NULL) {
            if (end_isAttached(end) && end_isStubEnd(end)
                    && end_getPseudoAdjacency(end) == NULL) {
                stList *list = stList_construct();
                do { //This loop constructs a pseudo chromosome..
                    Group *group = end_getGroup(end);
                    assert(!group_isLeaf(group));
                    Flower *nestedFlower = group_getNestedFlower(group);
                    assert(flower_getReference(nestedFlower) != NULL);
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
                /*
                 * Now we can finally build the pseudo chromosome.
                 */
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
#ifdef BEN_DEBUG
    reference_check(flower_getReference(flower));
    flower_check(flower);
#endif
}
