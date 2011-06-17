#include "cactus.h"
#include "sonLib.h"
#include "cycleConstrainedMatchingAlgorithms.h"

static stList *getExtraAttachedStubsFromParent(Flower *flower) {
    /*
     * Copy any stubs not present in the flower from the parent group.
     * At the end there will be an even number of attached stubs in the problem.
     * Returns the list of new ends, which need to be assigned a group.
     */
    Group *parentGroup = flower_getParentGroup(flower);
    if (parentGroup != NULL) {
        Group *parentEndIt = group_getEndIterator(parentGroup);
        End *parentEnd;
        while ((parentEnd = group_getNextEnd(parentEndIt)) != NULL) {
            if (end_isAttached(end) || end_isBlockEnd(end)) {
                End *end = flower_getEnd(flower, end_getName(parentEnd));
                if (end == NULL) { //We have found an end that needs to be pushed into the child.
                    stList_append(newEnds, end_copyConstruct(parentEnd, flower)); //At this point it has no associated group;
                }
            }
        }
        group_destructEndIterator(parentEndIt);
    }
    assert(flower_getAttachedStubEndNumber(flower) % 2 == 0);
    assert(flower_getAttachedStubEndNumber(flower) > 0);
    return newEnds;
}

static void iterateOverTangleEnds(Flower *flower,
        void(*fn)(End *end, void *extraArg), void *extraArg) {
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_isTangle(group)) {
            Group_EndIterator *endIt = group_getEndIterator(group);
            End *end;
            while ((end = group_getNextEnd(endIt)) != NULL) {
                fn(end, extraArg);
            }
            group_destructEndIterator(endIt);
        }
    }
    flower_destructGroupIterator(groupIt);
}


static void getMapsOfNonFreeStubEndsToIntegersP(End *end, void *extraArgs) {
    stHash *endsToInts = extraArgs[0];
    int32_t *endCount = extraArgs[1];
    stHash_insert(end, stIntTuple_construct(1, (*endCount)++));
}

stHash *getMapsOfNonFreeStubEndsToIntegers(Group *group) {
    stHash *endsToInts = stHash_construct2(NULL, (void (*)(void *))stIntTuple_destruct);
    int32_t endCount = 0;
    void *args = { endsToInts, &endCount };
    iterateOverTangleEnds(flower, getExtraAttachedStubsFromParentP, newEnds);
    return endsToInts;
}

static stHash *getIntsToEnds(stHash *endToInts) {
    stHash *intToEnds = stHash_construct3((int (*)(const void *))stIntTuple_hashKey, (int (*)(const void *, const void *))stIntTuple_hashKey, NULL, NULL);



    return intToEnds;
}

End *getAdjacentEndFromParent(End *end) {
    /*
     * Get the adjacent stub end by looking at the reference adjacency in the parent.
     */
    Flower *flower = end_getFlower(end);
    Group *parentGroup = flower_getParentGroup(flower);
    assert(parentGroup != NULL);
    End *parentEnd = group_getEnd(parentGroup, end_getName(end));
    assert(parentEnd != NULL);
    //Now trace the parent end..
    End *adjacentParentEnd = NULL;
    assert(adjacentParentEnd != NULL);
    End *adjacentEnd = flower_getEnd(flower, end_getName(adjacentParentEnd));
    assert(adjacentEnd != NULL);
    return adjacentEnd;
}

stList *getStubEdgesFromParent(Flower *flower, stHash *endsToInts) {
    /*
     * For each attached stub in the flower, get the end in the parent group, and add its
     * adjacency to the group.
     */
    stList *stubEdges = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end = NULL;
    stSortedSet *endsSeen = stSortedSet_construct();
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isAttached(end) && end_isStubEnd(end) && stSortedSet_search(
                endsSeen, end) == NULL) {
            assert(end_getGroup(end) != NULL);
            assert(group_isTangle(end_getGroup(end)));
            End *adjacentEnd = getAdjacentEndFromParent(end);
            assert(adjacentEnd != NULL);
            assert(end_getGroup(adjacentEnd) != NULL);
            assert(group_isTangle(end_getGroup(adjacentEnd)));
            assert(stSortedSet_search(endsSeen, adjacentEnd) == NULL);
            assert(adjacentEnd != end);
            stSortedSet_insert(endsSeen, end);
            stSortedSet_insert(endsSeen, adjacentEnd);
            stIntTuple *endInt = stHash_search(endsToInts, end);
            assert(endInt != NULL);
            stIntTuple *adjacenctEndInt =
                    stHash_search(endsToInts, adjacentEnd);
            assert(adjacentEndInt != NULL);
            stList_append(
                    stubEdges,
                    stIntTuple_construct(2, stIntTuple_getPosition(endInt, 0),
                            stIntTuple_getPosition(adjacenctEndInt, 0)));
        }
    }
    flower_destructEndIterator(endIt);
    stSortedSet_destruct(endsSeen);
    return stubEdges;
}

void getArbitraryStubEdgesP(stSortedSet *nodesSet, int32_t node) {
    stIntTuple *i = stIntTuple_construct(1, node);
    if(stSortedSet_search(nodesSet, i) != NULL) {
        stSortedSet_remove(nodesSet, i);
    }
    stIntTuple_destruct(i);
}

stList *getArbitraryStubEdges(stHash *endsToInts, stList *chainEdges) {
    stList *nodes = stHash_getValues(endsToInts);
    stSortedSet *nodesSet = stList_getSortedSet(endInts, (int (*)(const void *, const void *))stIntTuple_cmpFn);

    for(int32_t i=0; i<stList_length(chainEdges); i++) {
        stIntTuple *edge = stList_get(chainEdges, i);
        getArbitraryStubEdgesP(nodesSet, stIntTuple_getPosition(edge, 0));
        getArbitraryStubEdgesP(nodesSet, stIntTuple_getPosition(edge, 1));
    }

    assert(stSortedSet_size(nodesSet) > 0);
    assert(stSortedSet_size(nodesSet) % 2 == 0);
    stList_destruct(nodes);
    nodes = stSortedSet_getList(nodesSet);
    stList *stubEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int32_t i=0; i<stList_length(nodes); i+=2) {
        stList_append(stubEdges,
                stIntTuple_construct(2, stIntTuple_getPosition(stList_get(nodes, i), 0),
                stIntTuple_getPosition(stList_get(nodes, i+1), 0)));
    }
    return stubEdges;
}

End *getAdjacentEndArbitrarily(End *end, stHash *endsSeen) {
    /*
     * Return another stub end, randomly.
     */
    Flower *flower = end_getFlower(end);
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end2;
    while ((end2 = flower_getNextEnd(endIt)) != NULL) {
        if (stHash_search(endsSeen, end2) == NULL) {
            flower_destructEndIterator(endIt);
            return end2;
        }
    }
    flower_destructEndIterator(endIt);
    return NULL;
}

static stList *getStubEdges(Flower *flower, stHash *endsToInts) {
    /*
     * Get the stub edges for the flower.
     */
    return flower_getParentGroup(flower) != NULL ? getStubEdges2(flower,
            endsToInts, getAdjacentEndFromParent) : getStubEdges2(flower,
            endsToInts, getAdjacentEndArbitrarily);
}

stIntTuple *getEdge(End *end1, End *end2, stHash *endsToInts) {
    return stIntTuple_construct(2,
            stIntTuple_getPosition(stHash_search(endToInts, end1), 0),
            stIntTuple_getPosition(stHash_search(endToInts, end2), 0));
}

static stList *getChainEdges(Flower *flower, stHash *endsToInts) {
    /*
     * Get the chain edges for the flower.
     */
    stList *chainEdges = stList_construct();
    Flower_ChainIterator *chainIt = flower_getChainIterator(flower);
    Chain *chain;
    while ((chain = flower_getNextChain(chainIt)) != NULL) {
        End *_5End = link_get5End(chain_getFirst(chain));
        End *_3End = link_get3End(chain_getLast(chain));
        if (end_isBlockEnd(_5End) && end_isBlockEnd(_3End)) {
            End *end1 = end_getOtherBlockEnd(_5End);
            End *end1 = end_getOtherBlockEnd(_3End);
            assert(stHash_search(endToInsts, end1) != NULL);
            assert(stHash_search(endToInsts, end2) != NULL);
            stList_append(chainEdges, getEdge(end1, end2, endToInts));
        }
    }
    flower_destructChainIterator(chainIt);
    return chainEdges;
}

void getAdjacencyEdgesP(End *end, void *args) {
    stHash *endToInts = args[0];
    stList *adjacencyEdges = args[1];
    End_InstanceIterator *instanceIt = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(instanceIt)) != NULL) {
        cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
        if (cap_getSide(cap)) { //Use this side property to avoid adding edge twice.
            stList_append(adjacencyEdges,
                    getEdge(end, cap_getEnd(cap_getAdjacency(cap)), endToInts));
        }
    }
    end_destructInstanceIterator(instanceIt);
}

static stList *getAdjacencyEdges(Flower *flower, stHash *endsToInts) {
    /*
     * Get the adjacency edges for the flower.
     */
    stList *adjacencyEdges = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    void *args = { endsToInts, adjacencyEdges };
    iterateOverTangleEnds(flower, getAdjacencyEdgesP, args);
    stList_sort(adjacencyEdges,
            (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    stList *weightedAdjacencyEdges = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    int32_t weight = 0;
    stIntTuple *edge = NULL;
    while (stList_length(adjacencyEdges) > 0) {
        stIntTuple *edge2 = stList_pop(adjacencyEdges);
        if (edge != NULL && stIntTuple_cmpFn(edge, edge2) == 0) {
            weight++;
        } else {
            if (edge != NULL) {
                stList_append(
                        weightedAdjacencyEdges,
                        stIntTuple_construct(3,
                                stIntTuple_getPosition(edge, 0),
                                stIntTuple_getPosition(edge, 1), weight));
            }
            edge = edge2;
            weight = 1;
        }
    }
    if (edge != NULL) {
        stList_append(
                weightedAdjacencyEdges,
                stIntTuple_construct(3, stIntTuple_getPosition(edge, 0),
                        stIntTuple_getPosition(edge, 1), weight));
    }
    return adjacencyEdges;
}

End *getEndFromInt(stHash *intToEnds, int32_t node) {
    stIntTuple *i = stIntTuple_construct(1, node);
    End *end = stHash_search(intToEnds, i);
    assert(end != NULL);
    stIntTuple_destruct(i);
    return end;
}

static void addBridgeBlocks(Flower *flower, stList *chosenAdjacencyEdges, stHash *intsToEnds) {
    /*
     * For edge adjacency edge that bridges between to groups add a block with an end in each.
     */
    for(int32_t i=0; i<stList_length(chosenAdjacencyEdges); i++) {
        stIntTuple *edge = stList_get(chosenAdjacencyEdges, i);
        End *end1 = getEndFromInt(intsToEnds, stIntTuple_getPosition(edge, 0));
        End *end2 = getEndFromInt(intsToEnds, stIntTuple_getPosition(edge, 1));
        if(end_getGroup(end1) != end_getGroup(end2)) { //We need to add a bridge block.
            Block *block = block_construct(1, flower);
            end_setGroup(block_get5End(block), end1_getGroup())
        }
    }
}

static void addAdjacenciesAndSegments(Flower *flower,
        stList *chosenAdjacencyEdges) {
    /*
     * For each chosen adjacency edge ensure there exists a cap in each incident end (adding segments as needed),
     * and a node.
     */
}

void topDown(Flower *flower, Name eventName,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber)) {
    /*
     * Get any extra ends to balance the group from the parent problem
     */
    stList *newEnds = getExtraAttachedStubsFromParent(flower);
    int32_t nodeNumber = group_getAttachedStubEndNumber(flower)
            + group_getBlockEndNumber(flower);
    assert(nodeNumber % 2 == 0);
    assert(nodeNumber > 0);

    bool hasParent = flower_getParentGroup(flower) != NULL;

    /*
     * Create a map to identify the ends by integers.
     */
    stHash *endsToInts = getMapsOfNonFreeStubEndsToIntegers(flower);
    stHash *intsToEnds = getIntsToEnds(endToInts);

    /*
     * Get the chain edges.
     */
    stList *chainEdges = getChainEdges(group);

    /*
     * Get the stub edges.
     */
    stList *stubEdges = hasParent ? getStubEdgesFromParent(flower, endsToInts) : getArbitraryStubEdges(endsToInts, chainEdges);

    /*
     * Get the adjacency edges.
     */
    stList *adjacencyEdges = getAdjacencyEdges(group);

    /*
     * Calculate the matching
     */
    stList *chosenAdjacencyEdges = chooseMatching(nodeNumber, adjacencyEdges,
            stubEdges, chainEdges, !hasParent, matchingAlgorithm);

    /*
     * Add in any extra blocks to bridge between groups.
     */
    addBridgeBlocks(flower, chosenAdjacencyEdges, intsToEnds);

    /*
     * Get link adjacencies.
     */
    stList *linkEdges = getLinkEdges(flower);
    stList_appendAll(chosenAdjacencyEdges, linkEdges);

    /*
     * Add the reference genome into flower
     */
    addAdjacenciesAndSegments(flower, chosenAdjacencyEdges);

    /*
     * Ensure the newly created ends have a group.
     */
    assignGroups(newEnds);


    /*
     * Cleanup
     */
    stList_destruct(newEnds);
}
