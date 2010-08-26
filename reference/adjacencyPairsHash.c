/*
 * adjacencyPairsHash.c
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "cactus.h"
#include "adjacencyPairs.h"
#include "pressedFlowers.h"

void adjacencyHash_add(stHash *adjacenciesHash, AdjacencyPair *adjacencyPair) {
#ifdef BEN_DEBUG
    assert(stHash_search(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair)) == NULL);
    assert(stHash_search(adjacenciesHash, adjacencyPair_getEnd2(adjacencyPair)) == NULL);
#endif
    stHash_insert(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair),
            adjacencyPair);
    stHash_insert(adjacenciesHash, adjacencyPair_getEnd2(adjacencyPair),
            adjacencyPair);
}

void adjacencyHash_remove(stHash *adjacencies,
        AdjacencyPair *adjacencyPair) {
#ifdef BEN_DEBUG
    assert(stHash_search(adjacencies, adjacencyPair_getEnd1(adjacencyPair)) == adjacencyPair);
    assert(stHash_search(adjacencies, adjacencyPair_getEnd2(adjacencyPair)) == adjacencyPair);
#endif
    stHash_remove(adjacencies, adjacencyPair_getEnd1(adjacencyPair));
    stHash_remove(adjacencies, adjacencyPair_getEnd2(adjacencyPair));
}

void adjacencyHash_destruct(stHash *adjacencies, Flower *flower) {
    stHashIterator *endIterator = stHash_getIterator(adjacencies);
    End *end;
    stList *adjacencyPairs = stList_construct3(0, (void (*)(void *))adjacencyPair_destruct);
    while ((end = stHash_getNext(endIterator)) != NULL) {
        AdjacencyPair *adjacencyPair = stHash_search(adjacencies, end);
        if(adjacencyPair_getEnd1(adjacencyPair) == end) {
            stList_append(adjacencyPairs, adjacencyPair);
        }
    }
    stHash_destructIterator(endIterator);
    stList_destruct(adjacencyPairs);
    stHash_destruct(adjacencies);
}

static void makeListOfAdjacencyPairsP(Flower *flower, stList *adjacencies) {
    /*
     * Get a list of all adjacencies between ends in the flower.
     */
    End *end1;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    stHash *adjacenciesHash = stHash_construct3(
            (uint32_t(*)(const void *)) adjacencyPair_hashKey, (int(*)(
                    const void *, const void *)) adjacencyPair_hashEqual, NULL,
            NULL);
    while ((end1 = flower_getNextEnd(endIterator)) != NULL) {
        assert(!end_isBlockEnd(end1));
        if (end_isAttached(end1)) {
            assert(end_getPositiveOrientation(end1) == end1); //they are always in the positive orientation, right?
            End_InstanceIterator *capIterator = end_getInstanceIterator(end1);
            Cap *cap;
            while ((cap = end_getNext(capIterator)) != NULL) {
                Cap *cap2 = cap_getAdjacency(cap); //we allow both internal inferred and leaf adjacencies.
                if (cap2 != NULL) {
                    End *end2 = end_getPositiveOrientation(cap_getEnd(cap2));
                    assert(!end_isBlockEnd(end2));
                    if (end1 != end2 && end_isAttached(end2)) { //we don't allow free stubs or adjacency pairs which create self loops, as we can traverse a node twice!
                        AdjacencyPair *adjacencyPair = adjacencyPair_construct(
                                end1, end2);
                        if (stHash_search(adjacenciesHash, adjacencyPair)
                                == NULL) { //
                            stHash_insert(adjacenciesHash, adjacencyPair,
                                    adjacencyPair);
                            stList_append(adjacencies, adjacencyPair);
                        } else {
                            adjacencyPair_destruct(adjacencyPair);
                        }
                    }
                }
            }
            end_destructInstanceIterator(capIterator);
        }
    }
    flower_destructEndIterator(endIterator);
    stHash_destruct(adjacenciesHash);
}

static stList *makeListOfAdjacencyPairs(Flower *flower) {
    /*
     * Makes a list of possible adjacency pairs between ends in the effected flowers of a flower.
     * Each adjacency pair will be created if there exists an actual adjacency between the pair of ends.
     */
    stList *pressedFlowers = getListOfPressedFlowers(flower);
    stList *adjacencies = stList_construct();
    while (stList_length(pressedFlowers) > 0) {
        makeListOfAdjacencyPairsP(stList_pop(pressedFlowers), adjacencies);
    }
    stList_destruct(pressedFlowers);
    return adjacencies;
}

static stHash *choosePairing(stList *adjacencies) {
    /*
     * Greedily picks the adjacencies from the list such that each end has one adjacency.
     * Destroys the input list in the process.
     */
    stHash *adjacenciesHash = stHash_construct();
#ifdef BEN_DEBUG
    double strength = INT32_MAX;
#endif
    while (stList_length(adjacencies) > 0) {
        AdjacencyPair *adjacencyPair = stList_pop(adjacencies);
#ifdef BEN_DEBUG
        double d = adjacencyPair_getStrengthOfAdjacencyPair(adjacencyPair);
        assert(d <= strength);
        strength = d;
#endif
        if (stHash_search(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair))
                == NULL && stHash_search(adjacenciesHash,
                adjacencyPair_getEnd2(adjacencyPair)) == NULL) {
            adjacencyHash_add(adjacenciesHash, adjacencyPair);
        } else {
            adjacencyPair_destruct(adjacencyPair);
        }
    }
    assert(stList_length(adjacencies) == 0);
    stList_destruct(adjacencies);
    return adjacenciesHash;
}

static void addPseudoPseudoAdjacenciesP(stHash *adjacenciesHash, Flower *flower) {
    /*
     * Adds pseudo-pseudo adjacencies to ends that have no valid adjacency pair in the hash. Added adjacencies
     * are put in the hash and the list.
     */
    End *end1;
    End *end2 = NULL;
    assert(flower_isTerminal(flower));
    assert(flower_getGroupNumber(flower) == 1);
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);

    while ((end1 = flower_getNextEnd(endIterator))) {
        assert(!end_isBlockEnd(end1));
        assert(end_getOrientation(end1)); //should be positive orientation for this to work.
        if (end_isAttached(end1)) { //not interested in free stubs.
            if (stHash_search(adjacenciesHash, end1) == NULL) {
                if (end2 == NULL) {
                    end2 = end1;
                } else {
                    AdjacencyPair *adjacencyPair = adjacencyPair_construct(
                            end1, end2);
                    adjacencyHash_add(adjacenciesHash, adjacencyPair);
                    end2 = NULL;
                }
            }
        }
    }
    flower_destructEndIterator(endIterator);
    assert(end2 == NULL); //there must be an even number of pairs in each group.
}

static void addPseudoPseudoAdjacencies(stHash *adjacenciesHash, Flower *flower) {
    stList *pressedFlowers = getListOfPressedFlowers(flower);
    while (stList_length(pressedFlowers) > 0) {
        addPseudoPseudoAdjacenciesP(adjacenciesHash,
                stList_pop(pressedFlowers));
    }
    stList_destruct(pressedFlowers);
}

stHash *adjacencyHash_constructInitialPairs(Flower *flower) {
    stList *adjacencies = makeListOfAdjacencyPairs(flower);
    stList_sort(adjacencies,
            (int(*)(const void *, const void *)) adjacencyPair_cmpFnByStrength);
    stHash *adjacenciesHash = choosePairing(adjacencies);
    addPseudoPseudoAdjacencies(adjacenciesHash, flower);
    return adjacenciesHash;
}
