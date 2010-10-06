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
#include "reference.h"
#include "matchingAlgorithms.h"

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
#ifdef BEN_DEBUG
        assert(!end_isBlockEnd(end1));
#endif
        if (end_isAttached(end1)) {
            assert(end_getPositiveOrientation(end1) == end1); //they are always in the positive orientation, right?
            End_InstanceIterator *capIterator = end_getInstanceIterator(end1);
            Cap *cap;
            while ((cap = end_getNext(capIterator)) != NULL) {
                Cap *cap2 = cap_getAdjacency(cap); //we allow both internal inferred and leaf adjacencies.
                if (cap2 != NULL) {
                    End *end2 = end_getPositiveOrientation(cap_getEnd(cap2));
#ifdef BEN_DEBUG
                    assert(!end_isBlockEnd(end2));
#endif
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

static stSortedSet *makeSetOfEnds(Flower *flower) {
    /*
     * Makes a set of ends.
     */
    stSortedSet *ends = stSortedSet_construct();
    stList *pressedFlowers = getListOfPressedFlowers(flower);
    while (stList_length(pressedFlowers) > 0) {
        End *end;
        Flower_EndIterator *endIterator = flower_getEndIterator(stList_pop(pressedFlowers));
        while ((end = flower_getNextEnd(endIterator)) != NULL) {
            assert(!end_isBlockEnd(end));
            if (end_isAttached(end)) {
                stSortedSet_insert(ends, end);
            }
        }
        flower_destructEndIterator(endIterator);
    }
    stList_destruct(pressedFlowers);
    return ends;
}

static void addPseudoPseudoAdjacenciesP(stHash *adjacenciesHash, Flower *flower) {
    /*
     * Adds pseudo-pseudo adjacencies to ends that have no valid adjacency pair in the hash. Added adjacencies
     * are put in the hash and the list.
     */
    End *end1;
    End *end2 = NULL;
#ifdef BEN_DEBUG
    assert(flower_isTerminal(flower));
    assert(flower_getAttachedStubEndNumber(flower) % 2 == 0);
    assert(flower_getGroupNumber(flower) == 1);
#endif
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

static stHash *getEndsToIntsHash(stSortedSet *ends) {
    /*
     * Creates a hash of ends to integers.
     */
    stSortedSetIterator *it = stSortedSet_getIterator(ends);
    End *end;
    stHash *endsToInts = stHash_construct2(NULL, (void (*)(void *))stIntTuple_destruct);
    int32_t counter = 0;
    while((end = stSortedSet_getNext(it)) != NULL) {
        stIntTuple *intTuple = stIntTuple_construct(1, counter++);
        stHash_insert(endsToInts, end, intTuple);
    }
    stSortedSet_destructIterator(it);
    return endsToInts;
}

static stHash *getTupleEdgesToAdjacencyEdgesHash(stList *adjacencies, stHash *endsToInts) {
    /*
     * Creates a hash of tuple edges to adjacency edges.
     */
    stHash *tupleEdgesToAdjacencyEdges = stHash_construct3((uint32_t (*)(const void *))stIntTuple_hashKey,
            (int (*)(const void *, const void *))stIntTuple_equalsFn, (void (*)(void *))stIntTuple_destruct, NULL);
    for(int32_t i=0; i<stList_length(adjacencies); i++) {
        AdjacencyPair *adjacencyPair = stList_get(adjacencies, i);
#ifdef BEN_DEBUG
        assert(adjacencyPair != NULL);
#endif
        int32_t j = stIntTuple_getPosition(stHash_search(endsToInts, adjacencyPair_getEnd1(adjacencyPair)), 0);
        int32_t k = stIntTuple_getPosition(stHash_search(endsToInts, adjacencyPair_getEnd2(adjacencyPair)), 0);
        int32_t weight = adjacencyPair_getStrengthOfAdjacencyPair(adjacencyPair);
        //We add it in both directions
        stHash_insert(tupleEdgesToAdjacencyEdges, stIntTuple_construct(3, j, k, weight), adjacencyPair);
    }
    return tupleEdgesToAdjacencyEdges;
}

static stList *convertTupleEdgesToAdjacencies(stList *matching,
        stHash *tupleEdgesToAdjacencyEdges) {
    /*
     * Convert the tuple adjacencies into adjacencies.
     */
    stList *chosenAdjacencyPairs = stList_construct();
    for(int32_t i=0; i<stList_length(matching); i++) {
        stIntTuple *edge = stList_get(matching, i);
        AdjacencyPair *adjacencyPair = stHash_search(tupleEdgesToAdjacencyEdges, edge);
#ifdef BEN_DEBUG
        assert(adjacencyPair != NULL);
#endif
        stList_append(chosenAdjacencyPairs, adjacencyPair);
    }
    return chosenAdjacencyPairs;
}

static stHash *makeAdjacencyHashForMatching(stList *matching) {
    /*
     * Makes a hash of ends to adjacencies edges for a matching.
     */
    stHash *adjacenciesHash = stHash_construct();
    for(int32_t i=0; i<stList_length(matching); i++) {
        AdjacencyPair *adjacencyPair = stList_get(matching, i);
        adjacencyHash_add(adjacenciesHash, adjacencyPair);
    }
    return adjacenciesHash;
}

stHash *adjacencyHash_constructInitialPairs(Flower *flower, MatchingAlgorithm referenceAlgorithm) {
    /*
     * Get the initial list of adjacency pairs.
     */
    stSortedSet *ends = makeSetOfEnds(flower);
    stList *adjacencies = makeListOfAdjacencyPairs(flower);

    /*
     * Each end is assigned an integer, in random order starting from zero.
     */
    stHash *endsToInts = getEndsToIntsHash(ends);
    int32_t nodeNumber = stSortedSet_size(ends);
#ifdef BEN_DEBUG
    assert(nodeNumber % 2 == 0);
#endif
    stHash *tupleEdgesToAdjacencies = getTupleEdgesToAdjacencyEdgesHash(adjacencies, endsToInts);
#ifdef BEN_DEBUG
    assert(stHash_size(tupleEdgesToAdjacencies) == stList_length(adjacencies));
#endif
    stList *tupleEdges = stHash_getKeys(tupleEdgesToAdjacencies);
#ifdef BEN_DEBUG
    assert(stList_length(tupleEdges) == stHash_size(tupleEdgesToAdjacencies));
#endif
    /*
     * Create the matching.
     */
    stList *matchingTuples = NULL;
    if(referenceAlgorithm == greedy || nodeNumber > 500) { //We use greedy it the problem is too big
        matchingTuples = chooseMatching_greedy(tupleEdges, nodeNumber);
    }
    else if (referenceAlgorithm == blossom5) {
        matchingTuples = chooseMatching_blossom5(tupleEdges, nodeNumber);
    }
    else if (referenceAlgorithm == maxCardinality) {
        matchingTuples = chooseMatching_maximumCardinalityMatching(tupleEdges, nodeNumber);
    }
    else if (referenceAlgorithm == maxWeight) {
        matchingTuples = chooseMatching_maximumWeightMatching(tupleEdges, nodeNumber);
    }
    else {
        stThrowNew(REFERENCE_BUILDING_EXCEPTION, "Unrecognised matching algorithm: %i", referenceAlgorithm);
    }

    /*
     * Convert the edges in the tuple back to the adjacency edges.
     */
    stList *matchingAdjacencies = convertTupleEdgesToAdjacencies(matchingTuples, tupleEdgesToAdjacencies);
    assert(stList_length(matchingTuples) == stList_length(matchingAdjacencies));

    /*
     * We construct the hash for the chosen pairs.
     */
    stHash *adjacenciesHash = makeAdjacencyHashForMatching(matchingAdjacencies);
    assert(stHash_size(adjacenciesHash) == stList_length(matchingAdjacencies)*2);

    /*
     * Pseudo-pseudo adjacencies
     */
    addPseudoPseudoAdjacencies(adjacenciesHash, flower);
    assert(stHash_size(adjacenciesHash) == nodeNumber); //the matching is now perfect

    /*
     * We cleanup.
     */
    //Get rid of the temp files..
    stSortedSet_destruct(ends);
    stHash_destruct(endsToInts);
    stHash_destruct(tupleEdgesToAdjacencies);
    stList_destruct(tupleEdges);
    stList_destruct(matchingTuples);
    stList_destruct(matchingAdjacencies);

    //Cleanup the adjacency pairs we didn't use..
    while(stList_length(adjacencies) > 0) {
        AdjacencyPair *adjacencyPair = stList_pop(adjacencies);
        AdjacencyPair *chosenAdjacencyPair = stHash_search(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair));
        assert(chosenAdjacencyPair != NULL);
        if(adjacencyPair != chosenAdjacencyPair) {
            adjacencyPair_destruct(adjacencyPair);
        }
    }
    stList_destruct(adjacencies);

    return adjacenciesHash;
}
