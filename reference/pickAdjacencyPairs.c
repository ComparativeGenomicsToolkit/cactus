
#include "sonLib.h"
#include "cactus.h"
#include "adjacencyPairs.h"
#include "adjacencyPairsHash.h"

static stList *makeListOfAdjacencyPairs(Flower *flower) {
    /*
     * Get a list of all adjacencies between ends, excluding free stubs, in the flower.
     */
    End *end1;
    stList *adjacencies = stList_construct();
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    stHash *adjacenciesHash = stHash_construct3(
            (uint32_t(*)(const void *)) adjacencyPair_hashKey, (int(*)(
                    const void *, const void *)) adjacencyPair_hashEqual, NULL,
            NULL);
    while ((end1 = flower_getNextEnd(endIterator)) != NULL) {
        if (end_isBlockEnd(end1) || end_isAttached(end1)) {
            assert(end_getPositiveOrientation(end1) == end1); //they are always in the positive orientation, right?
            End_InstanceIterator *capIterator = end_getInstanceIterator(end1);
            Cap *cap;
            while ((cap = end_getNext(capIterator)) != NULL) {
                Cap *cap2 = cap_getAdjacency(cap); //we allow both internal inferred and leaf adjacencies.
                if (cap2 != NULL) {
                    End *end2 = end_getPositiveOrientation(cap_getEnd(cap2));
                    if (end1 != end2 && (end_isBlockEnd(end2)
                            || end_isAttached(end2))) { //we don't allow free stubs or adjacency pairs which create self loops, as we can traverse a node twice!
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
    return adjacencies;
}

static stHash *choosePairing(stList *adjacencies, Flower *flower) {
    /*
     * Greedily picks the adjacencies from the list such that each end has one adjacency.
     * Destroys the input list in the process.
     */
    assert(flower != NULL); //we may need this parameter later..
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

static void addPseudoPseudoAdjacencies(stHash *adjacenciesHash, Flower *flower) {
    /*
     * Adds pseudo-pseudo adjacencies to ends that have no valid adjacency pair in the hash. Added adjacencies
     * are put in the hash and the list.
     */
    End *end1;
    End *end2 = NULL;
    Flower_EndIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
#ifdef BEN_DEBUG
    	int32_t i = group_getEndNumber(group) - group_getFreeStubEndNumber(group);
    	assert(i > 0);
    	assert(i % 2 == 0);
#endif
    	Group_EndIterator *endIterator = group_getEndIterator(group);
    	while((end1 = group_getNextEnd(endIterator))) {
			if (end_isBlockEnd(end1) || end_isAttached(end1)) { //not interested in free stubs.
				assert(end_getPositiveOrientation(end1) == end1); //should be positive orientation for this to work.
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
    	group_destructEndIterator(endIterator);
    	assert(end2 == NULL); //there must be an even number of pairs in each group.
    }
    flower_destructGroupIterator(groupIterator);
}

stHash *pickAdjacencyPairs(Flower *flower, Reference *reference) {
    stList *adjacencies = makeListOfAdjacencyPairs(flower);
    stList_sort(adjacencies, (int (*)(const void *, const void *))adjacencyPair_cmpFnByStrength);
    stHash *adjacenciesHash = choosePairing(adjacencies, flower);
    addPseudoPseudoAdjacencies(adjacenciesHash, flower);
#ifdef BEN_DEBUG
    assert(stHash_size(adjacenciesHash) == flower_getEndNumber(flower) - flower_getFreeStubEndNumber(flower));
#endif
    return adjacenciesHash;
}
