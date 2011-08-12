#include "cactus.h"
#include "sonLib.h"
#include "cactusCycleConstrainedMatchingAlgorithms.h"
#include "cactusMatchingAlgorithms.h"
#include "perfectMatching.h"
#include "shared.h"
#include "checkEdges.h"

const char *REFERENCE_BUILDING_EXCEPTION = "REFERENCE_BUILDING_EXCEPTION";

////////////////////////////////////
////////////////////////////////////
//Get the reference event
////////////////////////////////////
////////////////////////////////////

static Event *getReferenceEvent(Flower *flower,
        const char *referenceEventHeader) {
    /*
     * Create the reference event, with the given header, for the flower.
     */
    EventTree *eventTree = flower_getEventTree(flower);
    Event *referenceEvent = eventTree_getEventByHeader(eventTree,
            referenceEventHeader);
    if (referenceEvent == NULL) {
        Group *parentGroup = flower_getParentGroup(flower);
        if (parentGroup == NULL) {
            //We are the root, so make a new event..
            return event_construct3(referenceEventHeader, INT32_MAX,
                    eventTree_getRootEvent(eventTree), eventTree);
        } else {
            Event *parentEvent = eventTree_getEventByHeader(
                    flower_getEventTree(group_getFlower(parentGroup)),
                    referenceEventHeader);
            assert(parentEvent != NULL);
            Event *event = event_construct(event_getName(parentEvent),
                    referenceEventHeader, INT32_MAX,
                    eventTree_getRootEvent(eventTree), eventTree);
            assert(event_getName(event) == event_getName(parentEvent));
            return event;
        }
    }
    return referenceEvent;
}

////////////////////////////////////
////////////////////////////////////
//Get ends from parent
////////////////////////////////////
////////////////////////////////////

static stList *getExtraAttachedStubsFromParent(Flower *flower) {
    /*
     * Copy any stubs not present in the flower from the parent group.
     * At the end there will be an even number of attached stubs in the problem.
     * Returns the list of new ends, which need to be assigned a group.
     */
    Group *parentGroup = flower_getParentGroup(flower);
    stList *newEnds = stList_construct();
    if (parentGroup != NULL) {
        Group_EndIterator *parentEndIt = group_getEndIterator(parentGroup);
        End *parentEnd;
        while ((parentEnd = group_getNextEnd(parentEndIt)) != NULL) {
            if (end_isAttached(parentEnd) || end_isBlockEnd(parentEnd)) {
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

////////////////////////////////////
////////////////////////////////////
//Functions for creating a map between integers representing the nodes of the graph and
//the ends in the flower
////////////////////////////////////
////////////////////////////////////

static stHash *getMapOfTangleEndsToNodes(Flower *flower) {
    /*
     * Iterates over new ends (from the parent, without a group, and those ends in tangles).
     */
    stHash *endsToNodes = stHash_construct2(NULL,
            (void(*)(void *)) stIntTuple_destruct);
    int32_t endCount = 0;
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        Group *group = end_getGroup(end);
        if (group != NULL) {
            if (group_isTangle(group)) {
                if (end_isAttached(end) || end_isBlockEnd(end)) {
                    stHash_insert(endsToNodes, end,
                            stIntTuple_construct(1, endCount++));
                }
            } else {
                assert(group_isLink(group));
            }
        } else {
            assert(end_isStubEnd(end));
            assert(end_isAttached(end));
            stHash_insert(endsToNodes, end, stIntTuple_construct(1, endCount++));
        }
    }
    flower_destructEndIterator(endIt);
    return endsToNodes;
}

////////////////////////////////////
////////////////////////////////////
//Chain edge weight functions
////////////////////////////////////
////////////////////////////////////

static Cap *traceAdjacency(Cap *cap, stSortedSet *activeEnds) {
    while (1) {
        cap = cap_getAdjacency(cap);
        assert(cap != NULL);
        End *end = end_getPositiveOrientation(cap_getEnd(cap));
        if (stSortedSet_search(activeEnds, end) != NULL) {
            return cap;
        }
        if (end_isStubEnd(cap_getEnd(cap))) {
            assert(end_isFree(cap_getEnd(cap)));
            return NULL;
        }
        cap = cap_getOtherSegmentCap(cap);
        assert(cap != NULL);
    }
}

static int getTotalDistinctConnections(End *end, stSortedSet *activeEnds) {
    stSortedSet *endsSeen = stSortedSet_construct();
    End_InstanceIterator *instanceIt = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(instanceIt)) != NULL) {
        Cap *adjacentCap = cap_getAdjacency(cap) != NULL ? traceAdjacency(cap,
                activeEnds) : NULL;
        if (adjacentCap != NULL) {
            End *adjacentEnd = end_getPositiveOrientation(cap_getEnd(adjacentCap));
            assert(end_getOrientation(adjacentEnd));
            assert(stSortedSet_search(activeEnds, adjacentEnd) != NULL);
            stSortedSet_insert(endsSeen, adjacentEnd);
        }
    }
    end_destructInstanceIterator(instanceIt);
    int32_t distinctConnections = stSortedSet_size(endsSeen);
    stSortedSet_destruct(endsSeen);
    return distinctConnections > 0 ? distinctConnections : 1;
}

static int getTotalBlockLength(Chain *chain) {
    int32_t i = 0;
    Link *link = chain_getFirst(chain);
    while (link != NULL) {
        End *end = link_get3End(link);
        if (end_isBlockEnd(end)) {
            i += block_getLength(end_getBlock(end));
        }
        else {
            assert(link_getPreviousLink(link) == NULL);
        }
        Link *nLink = link_getNextLink(link);
        if (nLink == NULL) {
            end = link_get5End(link);
            if (end_isBlockEnd(end)) {
                i += block_getLength(end_getBlock(end));
            }
        }
        else {
            assert(end_isBlockEnd(link_get5End(link)));
            assert(end_isBlockEnd(link_get3End(nLink)));
            assert(block_getPositiveOrientation(end_getBlock(link_get5End(link))) == block_getPositiveOrientation(end_getBlock(link_get3End(nLink))));
        }
        link = nLink;
    }
    return i;
}

static int getTotalNumberOfSequencesP(End *_3End, End *_5End) {
    return (end_getInstanceNumber(_3End) + end_getInstanceNumber(_5End)) / 2;
}

static int getTotalNumberOfSequences(Chain *chain) {
    return getTotalNumberOfSequencesP(link_get3End(chain_getFirst(chain)),
            link_get5End(chain_getLast(chain)));
}

static double getNonTrivialChainWeight(End *end1, End *end2, Chain *chain,
        int32_t code, stSortedSet *activeEnds) {
    double i = code % 2 == 0 ? chain_getAverageInstanceBaseLength(chain)
            : getTotalBlockLength(chain);
    if (code >= 4) {
        i *= getTotalNumberOfSequences(chain);
    }
    if (code == 2 || code == 3 || code == 6 || code == 7) {
        i /= getTotalDistinctConnections(end1, activeEnds)
                * getTotalDistinctConnections(end2, activeEnds);
    }
    return i;
}

static double getTrivialChainWeight(End *end1, End *end2, Block *block,
        int32_t code, stSortedSet *activeEnds) {
    double i = block_getLength(block);
    if (code >= 4) {
        i *= getTotalNumberOfSequencesP(block_get5End(block),
                block_get3End(block));
    }
    if (code == 2 || code == 3 || code == 6 || code == 7) {
        i /= getTotalDistinctConnections(end1, activeEnds)
                * getTotalDistinctConnections(end2, activeEnds);
    }
    assert(i >= 0);
    return i;
}

/*
 * Chain length function:
 * 0: Avg. instance length
 * 1: Total block length
 * 2: Avg. instance length / total number of connections
 * 3: Total block length / total number of connections
 * 4: Avg. instance length * number of sequences
 * 5. Total block length  * number of sequences
 * 6. Avg. instance length * number of sequences / total number of connections
 * 7. Total block length * number of sequences / total number of connections
 */
double getChainWeight(End *end1, End *end2, int32_t code,
        stSortedSet *activeEnds) {
    assert(code >= 0 && code <= 7);
    assert(end_isBlockEnd(end1));
    assert(end_isBlockEnd(end2));
    assert(group_isTangle(end_getGroup(end1)));
    assert(group_isTangle(end_getGroup(end2)));
    if (end_getPositiveOrientation(end_getOtherBlockEnd(end1)) == end2) {
        assert(end_getPositiveOrientation(end_getOtherBlockEnd(end2)) == end1);
        return getTrivialChainWeight(end1, end2, end_getBlock(end1), code,
                activeEnds);
    } else {
        Link *link = group_getLink(end_getGroup(end_getOtherBlockEnd(end1)));
        assert(link != NULL);
        return getNonTrivialChainWeight(end1, end2, link_getChain(link), code,
                activeEnds);
    }
}

////////////////////////////////////
////////////////////////////////////
//Chain edges
////////////////////////////////////
////////////////////////////////////

End *getEndFromNode(stHash *nodesToEnds, int32_t node) {
    /*
     * Get the end for the given node.
     */
    stIntTuple *i = stIntTuple_construct(1, node);
    End *end = stHash_search(nodesToEnds, i);
    assert(end != NULL);
    assert(end_getOrientation(end));
    stIntTuple_destruct(i);
    return end;
}

static stIntTuple *getEdge2(End *end1, End *end2, stHash *endsToNodes) {
    assert(stHash_search(endsToNodes, end1) != NULL);
    assert(stHash_search(endsToNodes, end2) != NULL);
    return constructEdge(
            stIntTuple_getPosition(stHash_search(endsToNodes, end1), 0),
            stIntTuple_getPosition(stHash_search(endsToNodes, end2), 0));
}

static void getNonTrivialChainEdges(Flower *flower, stHash *endsToNodes,
        stList *chainEdges) {
    /*
     * Get the non-trivial chain edges for the flower.
     */
    Flower_ChainIterator *chainIt = flower_getChainIterator(flower);
    Chain *chain;
    while ((chain = flower_getNextChain(chainIt)) != NULL) {
        End *_5End = link_get3End(chain_getFirst(chain));
        End *_3End = link_get5End(chain_getLast(chain));
        if (end_isBlockEnd(_5End) && end_isBlockEnd(_3End)) {
            End *end1 = end_getOtherBlockEnd(_5End);
            End *end2 = end_getOtherBlockEnd(_3End);

            assert(end_getGroup(end1) != NULL);
            assert(end_getGroup(end2) != NULL);
            assert(group_isTangle(end_getGroup(end1)));
            assert(group_isTangle(end_getGroup(end2)));

            assert(stHash_search(endsToNodes, end1) != NULL);
            assert(stHash_search(endsToNodes, end2) != NULL);

            stList_append(chainEdges, getEdge2(end1, end2, endsToNodes));
        }
    }
    flower_destructChainIterator(chainIt);
}

static void getTrivialChainEdges(Flower *flower, stHash *endsToNodes,
        stList *chainEdges) {
    /*
     * Get the trivial chain edges for the flower.
     */
    Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIt)) != NULL) {
        End *_5End = block_get5End(block);
        End *_3End = block_get3End(block);
        assert(end_getGroup(_5End) != NULL);
        assert(end_getGroup(_3End) != NULL);
        if (group_isTangle(end_getGroup(_5End)) && group_isTangle(
                end_getGroup(_3End))) {
            stList_append(chainEdges, getEdge2(_5End, _3End, endsToNodes));
        }
    }
    flower_destructBlockIterator(blockIt);
}

static stList *getChainEdges(Flower *flower, stHash *endsToNodes) {
    /*
     * Get the trivial chain edges for the flower.
     */
    stList *chainEdges = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    getNonTrivialChainEdges(flower, endsToNodes, chainEdges);
    getTrivialChainEdges(flower, endsToNodes, chainEdges);
    return chainEdges;
}

static stHash *scoreChainEdges(stList *chainEdges, stSortedSet *activeEnds,
        stHash *nodesToEnds, int32_t code) {
    stHash *chainEdgeScores = stHash_construct3(
            (uint32_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn, NULL,
            (void(*)(void *)) stDoubleTuple_destruct);
    for (int32_t i = 0; i < stList_length(chainEdges); i++) {
        stIntTuple *chainEdge = stList_get(chainEdges, i);
        int32_t node1 = stIntTuple_getPosition(chainEdge, 0);
        int32_t node2 = stIntTuple_getPosition(chainEdge, 1);
        End *end1 = getEndFromNode(nodesToEnds, node1);
        End *end2 = getEndFromNode(nodesToEnds, node2);
        assert(end1 != NULL);
        assert(end2 != NULL);
        stDoubleTuple *score = stDoubleTuple_construct(1,
                getChainWeight(end1, end2, code, activeEnds));
        assert(stHash_search(chainEdgeScores, chainEdge) == NULL);
        stHash_insert(chainEdgeScores, chainEdge, score);
        assert(stHash_search(chainEdgeScores, chainEdge) != NULL);
    }
    return chainEdgeScores;
}

////////////////////////////////////
////////////////////////////////////
//Stub edges
////////////////////////////////////
////////////////////////////////////

static Cap *getCapWithEvent(End *end, Name eventName) {
    /*
     * Get the first encountered cap in the end with the given event.
     */
    Cap *cap;
    End_InstanceIterator *instanceIt = end_getInstanceIterator(end);
    while ((cap = end_getNext(instanceIt)) != NULL) {
        if (event_getName(cap_getEvent(cap)) == eventName) {
            end_destructInstanceIterator(instanceIt);
            return cap;
        }
    }
    end_destructInstanceIterator(instanceIt);
    return NULL;
}

static End *getAdjacentEndFromParent(End *end, Event *referenceEvent) {
    /*
     * Get the adjacent stub end by looking at the reference adjacency in the parent.
     */
    Flower *flower = end_getFlower(end);
    Group *parentGroup = flower_getParentGroup(flower);
    assert(parentGroup != NULL);
    End *parentEnd = group_getEnd(parentGroup, end_getName(end));
    assert(parentEnd != NULL);
    //Now trace the parent end..
    Cap *cap = getCapWithEvent(parentEnd, event_getName(referenceEvent));
    assert(cap != NULL);
    Cap *adjacentCap = cap_getAdjacency(cap);
    assert(adjacentCap != NULL);
    End *adjacentParentEnd = cap_getEnd(adjacentCap);
    End *adjacentEnd = flower_getEnd(flower, end_getName(adjacentParentEnd));
    assert(adjacentEnd != NULL);
    return adjacentEnd;
}

static stList *getStubEdgesFromParent(Flower *flower, stHash *endsToNodes,
        Event *referenceEvent) {
    /*
     * For each attached stub in the flower, get the end in the parent group, and add its
     * adjacency to the group.
     */
    stList *stubEdges = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    stList *ends = stHash_getKeys(endsToNodes);
    stSortedSet *endsSeen = stSortedSet_construct();
    for (int32_t i = 0; i < stList_length(ends); i++) {
        End *end = stList_get(ends, i);
        assert(end_getOrientation(end));
        if (end_isStubEnd(end) && stSortedSet_search(endsSeen, end) == NULL) {
            End *adjacentEnd = getAdjacentEndFromParent(end, referenceEvent);
            assert(end_getOrientation(adjacentEnd));

            stSortedSet_insert(endsSeen, end);
            stSortedSet_insert(endsSeen, adjacentEnd);
            stList_append(stubEdges, getEdge2(end, adjacentEnd, endsToNodes));
        }
    }
    stSortedSet_destruct(endsSeen);
    stList_destruct(ends);
    return stubEdges;
}

////////////////////////////////////
////////////////////////////////////
//Adjacency edges
////////////////////////////////////
////////////////////////////////////

static void getAdjacencyEdgesP(End *end, stSortedSet *activeEnds,
        stHash *endsToNodes, stList *adjacencyEdges) {
    /*
     * This function adds adjacencies
     */
    End_InstanceIterator *instanceIt = end_getInstanceIterator(end);
    Cap *cap;
    assert(stHash_search(endsToNodes, end) != NULL);
    assert(end_isAttached(end) || end_isBlockEnd(end));
    assert(end_getOrientation(end));
    while ((cap = end_getNext(instanceIt)) != NULL) {
        if (cap_getSequence(cap) != NULL) {
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            if (cap_getSide(cap)) { //Use this side property to avoid adding edges twice.
                Cap *adjacentCap = traceAdjacency(cap, activeEnds);
                End *adjacentEnd;
                if (adjacentCap != NULL) {
                    assert(traceAdjacency(adjacentCap, activeEnds) == cap);
                    adjacentEnd = end_getPositiveOrientation(
                            cap_getEnd(adjacentCap));
                    if (adjacentEnd != end) {
                        assert(
                                end_isAttached(adjacentEnd) || end_isBlockEnd(
                                        adjacentEnd));
                        assert(end_getGroup(end) != NULL);
                        assert(end_getGroup(adjacentEnd) != NULL);
                        assert(group_isTangle(end_getGroup(end)));
                        assert(group_isTangle(end_getGroup(adjacentEnd)));
                        assert(stHash_search(endsToNodes, adjacentEnd) != NULL);
                        stList_append(adjacencyEdges,
                                getEdge2(end, adjacentEnd, endsToNodes));
                    }
                }
            }
        }
    }
    end_destructInstanceIterator(instanceIt);
}

static stList *getWeightedAdjacencyEdges(stList *adjacencyEdges) {
    /*
     * Iterate over the list, building weighted edges according to the number of duplicates.
     */

    /*
     * Sort the original adjacencies.
     */
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
            stIntTuple_destruct(edge2);
        } else {
            if (edge != NULL) {
                assert(weight >= 1);
                stList_append(
                        weightedAdjacencyEdges,
                        constructWeightedEdge(stIntTuple_getPosition(edge, 0),
                                stIntTuple_getPosition(edge, 1), weight));
                stIntTuple_destruct(edge);
            }
            edge = edge2;
            weight = 1;
        }
    }
    if (edge != NULL) {
        assert(weight >= 1);
        stList_append(
                weightedAdjacencyEdges,
                constructWeightedEdge(stIntTuple_getPosition(edge, 0),
                        stIntTuple_getPosition(edge, 1), weight));
        stIntTuple_destruct(edge);
    }
    return weightedAdjacencyEdges;
}

static stList *getAdjacencyEdges(Flower *flower, stSortedSet *activeEnds,
        stHash *endsToNodes) {
    /*
     * Get the adjacency edges for the flower.
     */

    /*
     * First build a list of adjacencies from the tangle ends.
     */
    stList *adjacencyEdges = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    stSortedSetIterator *endIt = stSortedSet_getIterator(activeEnds);
    End *end;
    while ((end = stSortedSet_getNext(endIt)) != NULL) {
        getAdjacencyEdgesP(end, activeEnds, endsToNodes, adjacencyEdges);
    }
    stSortedSet_destructIterator(endIt);

    /*
     * Convert to weighted adjacency edges.
     */
    stList *weightedAdjacencyEdges = getWeightedAdjacencyEdges(adjacencyEdges);
    stList_destruct(adjacencyEdges);

    return weightedAdjacencyEdges;
}

////////////////////////////////////
////////////////////////////////////
//Functions to add adjacencies and segments, given the chosen edges
////////////////////////////////////
////////////////////////////////////

static Cap *makeCapWithEvent(End *end, Event *referenceEvent) {
    /*
     * Returns a cap in the end that is part of the given reference event.
     */
    Cap *cap = getCapWithEvent(end, event_getName(referenceEvent));
    if (cap == NULL) {
        if (end_isBlockEnd(end)) {
            Block *block = end_getBlock(end);
            Segment *segment = segment_construct(block, referenceEvent);
            if (block_getRootInstance(block) != NULL) {
                segment_makeParentAndChild(block_getRootInstance(block),
                        segment);
            }
            cap = getCapWithEvent(end, event_getName(referenceEvent));
        } else { //Look in parent problem for event
            assert(end_isAttached(end));
            Group *parentGroup = flower_getParentGroup(end_getFlower(end));
            if (parentGroup != NULL) {
                End *parentEnd = group_getEnd(parentGroup, end_getName(end));
                assert(parentEnd != NULL);
                Cap *parentCap = getCapWithEvent(parentEnd,
                        event_getName(referenceEvent));
                assert(parentCap != NULL);
                cap = cap_copyConstruct(end, parentCap);
            } else {
                cap = cap_construct(end, referenceEvent);
            }
            if (end_getRootInstance(end) != NULL) {
                cap_makeParentAndChild(end_getRootInstance(end), cap);
            }
        }
    }
    assert(cap != NULL);
    assert(getCapWithEvent(end, event_getName(referenceEvent)) != NULL);
    return cap;
}

static void addAdjacenciesAndSegmentsP(End *end1, End *end2,
        Event *referenceEvent) {
    Cap *cap1 = makeCapWithEvent(end1, referenceEvent);
    Cap *cap2 = makeCapWithEvent(end2, referenceEvent);
    assert(cap_getAdjacency(cap1) == NULL);
    assert(cap_getAdjacency(cap2) == NULL);
    cap_makeAdjacent(cap1, cap2);
}

static void addLinkAdjacenciesAndSegments(Flower *flower, Event *referenceEvent) {
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_isLink(group)) {
            Link *link = group_getLink(group);
            addAdjacenciesAndSegmentsP(link_get3End(link), link_get5End(link),
                    referenceEvent);
        }
    }
    flower_destructGroupIterator(groupIt);
}

static void addTangleAdjacenciesAndSegments(Flower *flower,
        stList *chosenAdjacencyEdges, stHash *nodesToEnds,
        Event *referenceEvent) {
    /*
     * For each chosen adjacency edge ensure there exists a reference cap in each incident end (adding segments as needed),
     * then create the an adjacency for each cap. Where adjacencies link ends in different groups we add extra blocks to 'bridge' the
     * groups.
     */
    for (int32_t i = 0; i < stList_length(chosenAdjacencyEdges); i++) {
        stIntTuple *edge = stList_get(chosenAdjacencyEdges, i);
        End *end1 =
                getEndFromNode(nodesToEnds, stIntTuple_getPosition(edge, 0));
        End *end2 =
                getEndFromNode(nodesToEnds, stIntTuple_getPosition(edge, 1));
        assert(end1 != end2);
        if (end_getGroup(end1) != NULL && end_getGroup(end2) != NULL
                && end_getGroup(end1) != end_getGroup(end2)) {
            /*
             * We build a 'bridging block'.
             */
            Block *block = block_construct(1, flower);
            if (flower_builtTrees(flower)) { //Add a root segment
                Event *event = eventTree_getRootEvent(
                        flower_getEventTree(flower));
                assert(event != NULL);
                block_setRootInstance(block, segment_construct(block, event));
            }

            end_setGroup(block_get5End(block), end_getGroup(end1));
            assert(group_isTangle(end_getGroup(end1)));
            addAdjacenciesAndSegmentsP(end1, block_get5End(block),
                    referenceEvent);

            end_setGroup(block_get3End(block), end_getGroup(end2));
            assert(group_isTangle(end_getGroup(end2)));
            addAdjacenciesAndSegmentsP(end2, block_get3End(block),
                    referenceEvent);
        } else {
            addAdjacenciesAndSegmentsP(end1, end2, referenceEvent);
        }
    }
}

////////////////////////////////////
////////////////////////////////////
//Put the new ends into groups.
////////////////////////////////////
////////////////////////////////////

static void assignGroups(stList *newEnds, Flower *flower, Event *referenceEvent) {
    /*
     * Put ends into groups.
     */
    for (int32_t i = 0; i < stList_length(newEnds); i++) {
        End *end = stList_get(newEnds, i);
        if (end_getGroup(end) == NULL) { //It may already have been assigned.
            Cap *cap = getCapWithEvent(end, event_getName(referenceEvent));
            assert(cap != NULL);
            Cap *adjacentCap = cap_getAdjacency(cap);
            assert(adjacentCap != NULL);
            End *adjacentEnd = cap_getEnd(adjacentCap);
            Group *group = end_getGroup(adjacentEnd);
            if (group == NULL) {
                /*
                 * Two new ends link to one another, hence we need to find an existing tangle group.
                 */
                Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
                while ((group = flower_getNextGroup(groupIt)) != NULL) {
                    if (group_isTangle(group)) {
                        break;
                    }
                }
                assert(group_isTangle(group));
                flower_destructGroupIterator(groupIt);
                end_setGroup(adjacentEnd, group);
            }
            end_setGroup(end, group);
        }
    }
}

////////////////////////////////////
////////////////////////////////////
//Main function
////////////////////////////////////
////////////////////////////////////

stSortedSet *getEndsFromNodes(stSortedSet *nodes, stHash *nodesToEnds) {
    stSortedSet *ends = stSortedSet_construct();
    stSortedSetIterator *it = stSortedSet_getIterator(nodes);
    stIntTuple *node;
    while ((node = stSortedSet_getNext(it)) != NULL) {
        End *end = stHash_search(nodesToEnds, node);
        assert(end != NULL);
        stSortedSet_insert(ends, end);
    }
    stSortedSet_destructIterator(it);
    return ends;
}

void makeEdgesAClique(stList *edges, stSortedSet *nodes, int32_t defaultWeight) {
    /*
     * Adds edges to the list to make the set of edges a clique.
     */
    stHash *nodesToEdges = getNodesToEdgesHash(edges);
    stSortedSetIterator *nodeIt = stSortedSet_getIterator(nodes);
    stIntTuple *node;
    while ((node = stSortedSet_getNext(nodeIt)) != NULL) {
        int32_t i = stIntTuple_getPosition(node, 0);
        stSortedSetIterator *nodeIt2 = stSortedSet_getIteratorFrom(nodes, node);
        stIntTuple *node2 = stSortedSet_getNext(nodeIt2);
        assert(node == node2);
        while ((node2 = stSortedSet_getNext(nodeIt2)) != NULL) {
            assert(node2 != node);
            int32_t j = stIntTuple_getPosition(node2, 0);
            assert(i < j);
            if (getEdgeForNodes(i, j, nodesToEdges) == NULL) {
                stList_append(edges, constructWeightedEdge(i, j, defaultWeight));
            }
        }
        stSortedSet_destructIterator(nodeIt2);
    }
    stSortedSet_destructIterator(nodeIt);
    stHash_destruct(nodesToEdges);
}

stSortedSet *getActiveNodes(stList *chainEdges, stHash *nodesToEnds) {
    stList *nodes = stHash_getKeys(nodesToEnds);
    stSortedSet *chainNodes = getNodeSetOfEdges(chainEdges);

    stSortedSet *activeNodes = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(nodes); i++) {
        stIntTuple *node = stList_get(nodes, i);
        if (stSortedSet_search(chainNodes, node) == NULL) {
            addNodeToSet(activeNodes, stIntTuple_getPosition(node, 0));
        }
    }

    stSortedSet_destruct(chainNodes);
    stList_destruct(nodes);

    return activeNodes;
}

stSortedSet *getActiveEnds(stSortedSet *activeNodes, stHash *nodesToEnds) {
    stSortedSetIterator *it = stSortedSet_getIterator(activeNodes);
    stSortedSet *activeEnds = stSortedSet_construct();
    stIntTuple *node;
    while ((node = stSortedSet_getNext(it)) != NULL) {
        End *end = stHash_search(nodesToEnds, node);
        assert(end != NULL);
        assert(stSortedSet_search(activeEnds, end) == NULL);
        stSortedSet_insert(activeEnds, end);
    }
    stSortedSet_destructIterator(it);
    return activeEnds;
}

stList *getStubEdges(Flower *flower, stHash *endsToNodes,
        stSortedSet *activeEnds, stSortedSet *activeNodes,
        Event *referenceEvent,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber)) {
    /*
     * Gets the stub edges for a problem
     */
    if (flower_getParentGroup(flower) != NULL) {
        /*
         * Copy them from the parent.
         */
        return getStubEdgesFromParent(flower, endsToNodes, referenceEvent);
    }
    /*
     * Create a matching for the parent stub edges.
     */
    stList *allAdjacencyEdges = getAdjacencyEdges(flower, activeEnds,
            endsToNodes);
    checkEdges(allAdjacencyEdges, activeNodes, 1, 0);
    makeEdgesAClique(allAdjacencyEdges, activeNodes, 0);
    stList *chosenAdjacencyEdges = getPerfectMatching(activeNodes,
            allAdjacencyEdges, matchingAlgorithm);
    stList *stubEdges = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(chosenAdjacencyEdges); i++) {
        stIntTuple *adjacencyEdge = stList_get(chosenAdjacencyEdges, i);
        stList_append(
                stubEdges,
                constructEdge(stIntTuple_getPosition(adjacencyEdge, 0),
                        stIntTuple_getPosition(adjacencyEdge, 1)));
    }
    stList_destruct(chosenAdjacencyEdges);
    stList_destruct(allAdjacencyEdges);
    return stubEdges;
}

stHash *getActiveChainEdgesP2 = NULL;
int getActiveChainEdgesP(const void *o, const void *o2) {
    stDoubleTuple *score1 = stHash_search(getActiveChainEdgesP2, (void *) o);
    stDoubleTuple *score2 = stHash_search(getActiveChainEdgesP2, (void *) o2);
    assert(score1 != NULL);
    assert(score2 != NULL);
    assert(stDoubleTuple_getPosition(score1, 0) >= 0);
    assert(stDoubleTuple_getPosition(score2, 0) >= 0);
    return stDoubleTuple_cmpFn(score1, score2);
}

stList *getActiveChainEdges(stList *chainEdges, stHash *nodesToEnds,
        stSortedSet *activeEnds, stSortedSet *activeNodes,
        int32_t maxNumberOfChainsToSolvePerRound, int32_t chainWeightCode) {
    /*
     * Get the set of chain edges to add to the matching.
     */

    /*
     * Update the scores of the chain edges, given the new set of active nodes.
     */
    getActiveChainEdgesP2 = scoreChainEdges(chainEdges, activeEnds,
            nodesToEnds, chainWeightCode);
    stList_sort(chainEdges, getActiveChainEdgesP);

    for (int32_t i = 0; i < stList_length(chainEdges) - 1; i++) {
        assert(
                stDoubleTuple_getPosition(
                        stHash_search(getActiveChainEdgesP2,
                                stList_get(chainEdges, i)), 0)
                        <= stDoubleTuple_getPosition(
                                stHash_search(getActiveChainEdgesP2,
                                        stList_get(chainEdges, i + 1)), 0));
    }

    /*
     * Get the chain edges and update the adjacency edges and active nodes as we go.
     */
    stList *activeChainEdges = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    double pChainScore = INT32_MAX;

    while (stList_length(chainEdges) > 0 && stList_length(activeChainEdges)
            < maxNumberOfChainsToSolvePerRound) {
        stIntTuple *chainEdge = stList_pop(chainEdges);
        stList_append(activeChainEdges, chainEdge);

        /*
         * Update the set of active nodes and ends.
         */
        int32_t node1 = stIntTuple_getPosition(chainEdge, 0);
        int32_t node2 = stIntTuple_getPosition(chainEdge, 1);
        addNodeToSet(activeNodes, node1);
        addNodeToSet(activeNodes, node2);
        End *end1 = getEndFromNode(nodesToEnds, node1);
        assert(end1 != NULL);
        End *end2 = getEndFromNode(nodesToEnds, node2);
        assert(end2 != NULL);
        assert(stSortedSet_search(activeEnds, end1) == NULL);
        assert(stSortedSet_search(activeEnds, end2) == NULL);
        stSortedSet_insert(activeEnds, end1);
        stSortedSet_insert(activeEnds, end2);

        stDoubleTuple *score = stHash_search(getActiveChainEdgesP2, chainEdge);
        assert(score != NULL);
        assert(
                pChainScore == INT32_MAX || stDoubleTuple_getPosition(score, 0)
                        <= pChainScore);
        pChainScore = stDoubleTuple_getPosition(score, 0);
    }
    stHash_destruct(getActiveChainEdgesP2);

    return activeChainEdges;
}

void buildReferenceTopDown(Flower *flower, const char *referenceEventHeader,
        int32_t maxNumberOfChainsToSolvePerRound,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber),
        int32_t chainWeightCode, bool recalculateMatchingEachCycle) {
    /*
     * Implements the following pseudo code.
     *
     * S = initial constraints
     * C = chains
     * G = adjacency graph
     * X = integer > 0
     * def layeredMatching(C, S, G, X):
     *      while C not empty:
     *          C' = largest(C, X) #Get largest X chains
     *          N = nodes(C \cup S) #Get nodes at ends of constraint and chain edges
     *          A = adjacencies(N, G) #Get adjacencies between N, as traced through G
     *          S = matching(N, C', S, A) #Calculate the best matching, respecting set of constraints S
     *          C = C \ C' #Update the list of chains left to process
     *          return S #All edges ultimately become stubs
     */

    /*
     * Get the reference event
     */
    Event *referenceEvent = getReferenceEvent(flower, referenceEventHeader);

    /*
     * Get any extra ends to balance the group from the parent problem.
     */
    stList *newEnds = getExtraAttachedStubsFromParent(flower);

    /*
     * Create a map to identify the ends by integers.
     */
    stHash *endsToNodes = getMapOfTangleEndsToNodes(flower);
    stHash *nodesToEnds = stHash_invert(endsToNodes,
            (uint32_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn, NULL,
            NULL);
    int32_t nodeNumber = stHash_size(endsToNodes);
    assert(nodeNumber % 2 == 0);

    /*
     * Get the chain edges.
     */
    stList *chainEdges = getChainEdges(flower, endsToNodes);

    /*
     * Get the active nodes..
     */
    stSortedSet *activeNodes = getActiveNodes(chainEdges, nodesToEnds);
    stSortedSet *activeEnds = getActiveEnds(activeNodes, nodesToEnds);

    /*
     * Get the stub edges and chosen edges.
     */
    stList * stubEdges = getStubEdges(flower, endsToNodes, activeEnds,
            activeNodes, referenceEvent, matchingAlgorithm);

    /*
     * Check the edges and nodes before starting to calculate the matching.
     */
    st_logDebug(
            "Starting to build the reference for flower %lli with %i stubs and %i ends\n",
            flower_getName(flower), stSortedSet_size(activeNodes),
            stHash_size(nodesToEnds));

    assert(stSortedSet_size(activeEnds) == stSortedSet_size(activeNodes));
    assert(stList_length(stubEdges) * 2 == stSortedSet_size(activeNodes));
    for (int32_t i = 0; i < stList_length(stubEdges); i++) {
        stIntTuple *edge = stList_get(stubEdges, i);
        assert(nodeInSet(activeNodes, stIntTuple_getPosition(edge, 0)));
        assert(nodeInSet(activeNodes, stIntTuple_getPosition(edge, 1)));
    }
    assert(
            2 * (stList_length(stubEdges) + stList_length(chainEdges))
                    == nodeNumber);

    while (stList_length(chainEdges) > 0) {
        /*
         * Get the chain edges to consider.
         */
        stList *activeChainEdges = getActiveChainEdges(chainEdges, nodesToEnds,
                activeEnds, activeNodes, maxNumberOfChainsToSolvePerRound,
                chainWeightCode);

        /*
         * Get the adjacency edges
         */
        stList *nonZeroWeightAdjacencyEdges = getAdjacencyEdges(flower,
                activeEnds, endsToNodes);
        stList_setDestructor(nonZeroWeightAdjacencyEdges, NULL);
        stList *allAdjacencyEdges = stList_copy(nonZeroWeightAdjacencyEdges,
                (void(*)(void *)) stIntTuple_destruct);
        makeEdgesAClique(allAdjacencyEdges, activeNodes, 0);
        stSortedSet *allAdjacencyEdgesSet = stList_getSortedSet(
                allAdjacencyEdges,
                (int(*)(const void *, const void *)) stIntTuple_cmpFn);

        stList *chosenAdjacencyEdges;
        if (recalculateMatchingEachCycle) {
            /*
             * Recalculate the matching
             */
            chosenAdjacencyEdges = getPerfectMatching(activeNodes,
                    allAdjacencyEdges, matchingAlgorithm);
        } else {
            /*
             * Add the chain edges to the matching
             */
            chosenAdjacencyEdges = stList_construct();
            stList *stubAndChainEdges = stList_copy(stubEdges, NULL);
            stList_appendAll(stubAndChainEdges, activeChainEdges);
            for (int32_t i = 0; i < stList_length(stubAndChainEdges); i++) {
                stIntTuple *edge = stList_get(stubAndChainEdges, i);
                stIntTuple *adjacencyEdge = getWeightedEdgeFromSet(
                        stIntTuple_getPosition(edge, 0),
                        stIntTuple_getPosition(edge, 1), allAdjacencyEdgesSet);
                assert(adjacencyEdge != NULL);
                stList_append(chosenAdjacencyEdges, adjacencyEdge);
            }
            stList_destruct(stubAndChainEdges);
        }
        assert(
                stList_length(chosenAdjacencyEdges) * 2 == stSortedSet_size(
                        activeNodes));

        st_logDebug(
                "Going to build a reference matching for the flower %lli with %i nodes, %i adjacencies with cardinality %i and weight %i,  %i stub edges and %i chain edges, %i attached stubs, %i free stubs and %i block ends,  %i ends total and %i chains\n",
                flower_getName(flower), stSortedSet_size(activeNodes),
                stList_length(allAdjacencyEdges),
                matchingCardinality(allAdjacencyEdges),
                matchingWeight(allAdjacencyEdges), stList_length(stubEdges),
                stList_length(chainEdges),
                flower_getAttachedStubEndNumber(flower),
                flower_getFreeStubEndNumber(flower),
                flower_getBlockEndNumber(flower), flower_getEndNumber(flower),
                flower_getChainNumber(flower));

        /*
         * Make it obey the cyclic constraints..
         */
        stList *updatedChosenAdjacencyEdges =
                makeMatchingObeyCyclicConstraints(activeNodes,
                        chosenAdjacencyEdges, allAdjacencyEdgesSet,
                        nonZeroWeightAdjacencyEdges, stubEdges,
                        activeChainEdges, recalculateMatchingEachCycle);
        stList_destruct(chosenAdjacencyEdges);
        chosenAdjacencyEdges = updatedChosenAdjacencyEdges;

        /*
         * Update the stub edges.
         */
        stList_destruct(stubEdges);
        stubEdges = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
        while (stList_length(chosenAdjacencyEdges) > 0) {
            stIntTuple *edge = stList_pop(chosenAdjacencyEdges);
            stList_append(
                    stubEdges,
                    constructEdge(stIntTuple_getPosition(edge, 0),
                            stIntTuple_getPosition(edge, 1)));
        }

        /*
         * Cleanup
         */
        stList_destruct(chosenAdjacencyEdges);
        stList_destruct(nonZeroWeightAdjacencyEdges);
        stSortedSet_destruct(allAdjacencyEdgesSet);
        stList_destruct(allAdjacencyEdges);
        stList_destruct(activeChainEdges);
    }

    /*
     * Check the matching we have.
     */
    checkEdges(stubEdges, activeNodes, 1, 0);
    assert(stList_length(stubEdges) * 2 == nodeNumber);
    assert(stSortedSet_size(activeEnds) == stSortedSet_size(activeNodes));
    assert(stSortedSet_size(activeEnds) == nodeNumber);

    /*
     * Add the reference genome into flower
     */
    addLinkAdjacenciesAndSegments(flower, referenceEvent);
    addTangleAdjacenciesAndSegments(flower, stubEdges, nodesToEnds,
            referenceEvent);

    /*
     * Ensure the newly created ends have a group.
     */
    assignGroups(newEnds, flower, referenceEvent);

    /*
     * Cleanup
     */
    stSortedSet_destruct(activeEnds);
    stSortedSet_destruct(activeNodes);
    stList_destruct(newEnds);
    stHash_destruct(endsToNodes);
    stHash_destruct(nodesToEnds);
    stList_destruct(chainEdges);
    stList_destruct(stubEdges);
}

