#include "cactus.h"
#include "sonLib.h"
#include "checkEdges.h"
#include "adjacencyProblem.h"
#include "perfectMatching.h"
#include "checkEdges.h"
#include "cactusMatchingAlgorithms.h"
#include <math.h>

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
    if(flower_getBlockEndNumber(flower) > 0) {
        assert(flower_getAttachedStubEndNumber(flower) > 0);
    }
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
//Adjacency edges
////////////////////////////////////
////////////////////////////////////

static stList *calculateZP(Cap *cap, stHash *endsToNodes) {
    /*
     * Get the list of caps that represent the ends of the chains and stubs within a sequence.
     */
    assert(!cap_getSide(cap));
    assert(end_isStubEnd(end_getPositiveOrientation(cap_getEnd(cap))));
    stList *caps = stList_construct();
    bool b = 0;
    while (1) {
        End *end = end_getPositiveOrientation(cap_getEnd(cap));
        if (stHash_search(endsToNodes, end) != NULL) {
            assert(!cap_getSide(cap));
            if(stList_length(caps) > 0) {
                assert(b);
            }
            b = 0;
            stList_append(caps, cap);
        }
        cap = cap_getAdjacency(cap);
        assert(cap != NULL);
        end = end_getPositiveOrientation(cap_getEnd(cap));
        if (stHash_search(endsToNodes, end) != NULL) {
            assert(cap_getSide(cap));
            if(stList_length(caps) > 0) {
                assert(!b);
            }
            b = 1;
            stList_append(caps, cap);
        }
        if (end_isStubEnd(end)) {
            return caps;
        }
        assert(cap != cap_getOtherSegmentCap(cap));
        cap = cap_getOtherSegmentCap(cap);
        assert(cap != NULL);
    }
}

static Cap *calculateZP4(Cap *cap, stHash *endsToNodes) {
    if(cap_getOtherSegmentCap(cap) == NULL) {
        return NULL;
    }
    while (1) {
        cap = cap_getOtherSegmentCap(cap);
        assert(cap != NULL);
        End *end = end_getPositiveOrientation(cap_getEnd(cap));
        if (stHash_search(endsToNodes, end) != NULL) {
            return cap;
        }
        cap = cap_getAdjacency(cap);
        assert(cap != NULL);
        end = end_getPositiveOrientation(cap_getEnd(cap));
        assert(stHash_search(endsToNodes, end) == NULL);
        if (end_isStubEnd(end)) {
            if(!end_isFree(end)) {
                assert(!flower_hasParentGroup(end_getFlower(end)));
            }
            return NULL;
        }
    }
}

static int32_t calculateZP2(Cap *cap, stHash *endsToNodes) {
    /*
     * Calculate the length of a segment that can be traversed from a cap,
     * before hitting the end of the sequence of one of the other ends in the set
     * endsToNodes.
     */
    assert(cap_getStrand(cap));
    Sequence *sequence = cap_getSequence(cap);
    assert(sequence != NULL);
    Cap *otherCap = calculateZP4(cap, endsToNodes);
    int32_t capLength;
    if (otherCap == NULL) {
        capLength = cap_getSide(cap) ? sequence_getLength(sequence) + sequence_getStart(sequence)
                - cap_getCoordinate(cap) : cap_getCoordinate(cap)
                - sequence_getStart(sequence) + 1;
    } else {
        capLength = cap_getSide(cap) ? cap_getCoordinate(otherCap)
                - cap_getCoordinate(cap) + 1 : cap_getCoordinate(cap)
                - cap_getCoordinate(otherCap) + 1;
    }
    assert(capLength >= 0);
    if(capLength == 0) { //Give stubs at the end of sequences some weight.
        assert(end_isStubEnd(cap_getEnd(cap)));
        capLength = 1;
    }
    return capLength;
}

double *calculateZ(Flower *flower, stHash *endsToNodes, double theta) {
    /*
     * Calculate the zScores between all ends.
     */
    int32_t nodeNumber = stHash_size(endsToNodes);
    assert(nodeNumber % 2 == 0);
    double *z = st_calloc(nodeNumber * nodeNumber, sizeof(double));
    for(int32_t i=0; i<nodeNumber*nodeNumber; i++) { //Setting to zero for paranoid reasons.
        z[i] = 0.0;
    }
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end)) {
            End_InstanceIterator *capIt = end_getInstanceIterator(end);
            Cap *cap;
            while ((cap = end_getNext(capIt)) != NULL) {
                cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                if (!cap_getSide(cap) && cap_getSequence(cap) != NULL) {
                    stList *caps = calculateZP(cap, endsToNodes);

                    /*
                     * Calculate the lengths of the sequences following the 3 caps, for efficiency.
                     */
                    int32_t *capSizes = st_malloc(
                            sizeof(int32_t) * stList_length(caps));
                    for (int32_t i = 0; i < stList_length(caps); i++) {
                        Cap *cap = stList_get(caps, i);
                        capSizes[i] = calculateZP2(cap, endsToNodes);
                    }

                    /*
                     * Iterate through all pairs of 5' and 3' caps to calculate additions to scores.
                     */
                    for (int32_t i = (stList_length(caps) > 0 && cap_getSide(stList_get(caps, 0))) ? 1 : 0;
                            i < stList_length(caps); i+=2) {
                        Cap *_3Cap = stList_get(caps, i);
                        assert(!cap_getSide(_3Cap));
                        int32_t _3CapSize = capSizes[i];
                        int32_t _3Node = stIntTuple_getPosition(
                                stHash_search(
                                        endsToNodes,
                                        end_getPositiveOrientation(
                                                cap_getEnd(_3Cap))), 0);
                        assert(_3Node >= 0);
                        assert(_3Node < nodeNumber);
                        for (int32_t j = i+1; j < stList_length(caps); j+=2) {
                            Cap *_5Cap = stList_get(caps, j);
                            assert(cap_getSide(_5Cap));
                            int32_t _5Node = stIntTuple_getPosition(
                                    stHash_search(
                                            endsToNodes,
                                            end_getPositiveOrientation(
                                                    cap_getEnd(_5Cap))), 0);
                            int32_t _5CapSize = capSizes[j];
                            assert(_5Node >= 0);
                            assert(_5Node < nodeNumber);
                            int32_t diff = cap_getCoordinate(_5Cap)
                                    - cap_getCoordinate(_3Cap);
                            assert(diff >= 1);
                            double score = calculateZScore(_5CapSize, _3CapSize, diff, theta);
                            assert(score >= -0.0001);
                            if(score <= 0.0) {
                                score = 1e-10; //Make slightly non-zero.
                            }
                            assert(score > 0.0);
                            z[_5Node * nodeNumber + _3Node] += score;
                            z[_3Node * nodeNumber + _5Node] = z[_5Node * nodeNumber + _3Node];
                            assert(z[_5Node * nodeNumber + _3Node] == z[_3Node * nodeNumber + _5Node]);
                            assert(z[_5Node * nodeNumber + _3Node] > 0.0);
                        }
                    }
                    stList_destruct(caps);
                    free(capSizes);
                }
            }
            end_destructInstanceIterator(capIt);
        }
    }
    flower_destructEndIterator(endIt);
    for(int32_t i=0; i<nodeNumber*nodeNumber; i++) {
        assert(z[i] >= 0.0);
    }

    return z;
}

////////////////////////////////////
////////////////////////////////////
//Chain edges
////////////////////////////////////
////////////////////////////////////

static End *getEndFromNode(stHash *nodesToEnds, int32_t node) {
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

static stList *getStubEdges(Flower *flower, stHash *endsToNodes,
        Event *referenceEvent,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber), stList *chainEdges) {
    int32_t nodeNumber = stHash_size(endsToNodes);

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
     * Get stub nodes
     */
    stSortedSet *chainNodeSet = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, (void (*)(void *))stIntTuple_destruct);
    stList *stubNodes = stList_construct();
    for(int32_t i=0; i<stList_length(chainEdges); i++) {
        stIntTuple *chainEdge = stList_get(chainEdges, i);
        assert(stIntTuple_length(chainEdge) == 2);
        stIntTuple *node1 = stIntTuple_construct(1, stIntTuple_getPosition(chainEdge, 0));
        stIntTuple *node2 = stIntTuple_construct(1, stIntTuple_getPosition(chainEdge, 1));
        assert(stSortedSet_search(chainNodeSet, node1) == NULL);
        stSortedSet_insert(chainNodeSet, node1);
        assert(stSortedSet_search(chainNodeSet, node2) == NULL);
        stSortedSet_insert(chainNodeSet, node2);
    }
    stHashIterator *it = stHash_getIterator(endsToNodes);
    End *end;
    while((end = stHash_getNext(it)) != NULL) {
        stIntTuple *node = stHash_search(endsToNodes, end);
        if(stSortedSet_search(chainNodeSet, node) == NULL) {
            stList_append(stubNodes, node);
        }
    }
    stHash_destructIterator(it);
    stSortedSet_destruct(chainNodeSet);

    /*
     * Make Z for the stubs using a theta of 0.0
     */
    double *z = calculateZ(flower, endsToNodes, 0.0);

    st_logDebug("Building a matching for %i stub nodes in the top level problem from %i total stubs of which %i attached , %i total ends, %i chains, %i blocks %i groups and %i sequences\n", stList_length(stubNodes), flower_getStubEndNumber(flower), flower_getAttachedStubEndNumber(flower), flower_getEndNumber(flower), flower_getChainNumber(flower), flower_getBlockNumber(flower), flower_getGroupNumber(flower), flower_getSequenceNumber(flower));

    /*
     * Create a matching for the parent stub edges.
     */
    stList *adjacencyEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int32_t i=0; i<stList_length(stubNodes); i++) {
        int32_t node1 = stIntTuple_getPosition(stList_get(stubNodes, i), 0);
        for(int32_t j=i+1; j<stList_length(stubNodes); j++) {
            int32_t node2 = stIntTuple_getPosition(stList_get(stubNodes, j), 0);
            assert(z[node1 * nodeNumber + node2] >= 0);
            double score = round(z[node1 * nodeNumber + node2]);
            assert(score >= 0);
            int32_t score2 = score > INT32_MAX ? INT32_MAX : score;
            assert(score2 >= 0);
            stList_append(adjacencyEdges, constructWeightedEdge(node1, node2, score2));
        }
    }
    stSortedSet *stubNodesSet = stList_getSortedSet(stubNodes, (int (*)(const void *, const void *))stIntTuple_cmpFn);

    checkEdges(adjacencyEdges, stubNodesSet, 1, 0);

    stList *chosenAdjacencyEdges = getPerfectMatching(stubNodesSet,
            adjacencyEdges, matchingAlgorithm);

    stList *stubEdges = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(chosenAdjacencyEdges); i++) {
        stIntTuple *adjacencyEdge = stList_get(chosenAdjacencyEdges, i);
        assert(stIntTuple_length(adjacencyEdge) == 3);
        stList_append(
                stubEdges,
                constructEdge(stIntTuple_getPosition(adjacencyEdge, 0),
                        stIntTuple_getPosition(adjacencyEdge, 1)));
    }

    stList_destruct(chosenAdjacencyEdges);
    stList_destruct(adjacencyEdges);
    stSortedSet_destruct(stubNodesSet);
    stList_destruct(stubNodes);
    free(z);

    return stubEdges;
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
    assert(getCapWithEvent(end, event_getName(referenceEvent)) == cap);
    return cap;
}

static void addAdjacenciesAndSegmentsP(End *end1, End *end2,
        Event *referenceEvent) {
    Cap *cap1 = makeCapWithEvent(end1, referenceEvent);
    Cap *cap2 = makeCapWithEvent(end2, referenceEvent);
    assert(cap_getAdjacency(cap1) == NULL);
    assert(cap_getAdjacency(cap2) == NULL);
    assert(cap1 != cap2);
    assert(cap_getPositiveOrientation(cap1) != cap_getPositiveOrientation(cap2));
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

void buildReferenceTopDown(Flower *flower, const char *referenceEventHeader,
        int32_t permutations,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber),
        double (*temperature)(double),
        double theta) {
    /*
     * Implements a greedy algorithm and gibbs sampler to find a solution to the adjacency problem for a net.
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
     * Calculate z function
     */
    double *z = calculateZ(flower, endsToNodes, theta);

    /*
     * Get the chain edges.
     */
    stList *chainEdges = getChainEdges(flower, endsToNodes);

    /*
     * Get the stub edges and chosen edges.
     */
    stList * stubEdges = getStubEdges(flower, endsToNodes,
            referenceEvent, matchingAlgorithm, chainEdges);

    /*
     * Check the edges and nodes before starting to calculate the matching.
     */
    st_logDebug(
            "Starting to build the reference for flower %lli with %i stubs and %i ends\n",
            flower_getName(flower), stList_length(stubEdges), stList_length(chainEdges));

    assert(
            2 * (stList_length(stubEdges) + stList_length(chainEdges))
                    == nodeNumber);

    double maxPossibleScore = calculateMaxZ(nodeNumber, z);
    double totalScoreAfterGreedy = 0.0;
    stList *reference = makeReferenceGreedily(stubEdges, chainEdges, z, nodeNumber, &totalScoreAfterGreedy);
    double totalScoreAfterGreedy2 = calculateZScoreOfReference(reference, nodeNumber, z);
    st_logDebug("The score of the initial solution is %f (recalculated %f) out of a max possible %f\n", totalScoreAfterGreedy, totalScoreAfterGreedy2, maxPossibleScore);
    //assert(totalScoreAfterGreedy <= totalScoreAfterGreedy2 + 0.01);
    //assert(totalScoreAfterGreedy2 <= totalScoreAfterGreedy + 0.01);


    //gibbsSamplingWithSimulatedAnnealing(reference, chainEdges, z, permutations, temperature, 0);
    //double totalScoreAfterSimulatedAnnealing = calculateZScoreOfReference(reference, nodeNumber, z);
    //st_logDebug("The score of the sampled solution is %f after %i rounds of permutation out of a max possible %f\n", totalScoreAfterSimulatedAnnealing, permutations, maxPossibleScore);

    gibbsSamplingWithSimulatedAnnealing(reference, chainEdges, z, permutations, NULL, 1);
    double totalScoreAfterGreedySampling = calculateZScoreOfReference(reference, nodeNumber, z);
    st_logDebug("The score of the final solution is %f after %i rounds of greedy permutation out of a max possible %f\n", totalScoreAfterGreedySampling, permutations, maxPossibleScore);
    //assert(totalScoreAfterGreedySampling + 0.01 >= totalScoreAfterGreedy); //totalScoreAfterSimulatedAnnealing);

    stList *chosenEdges = convertReferenceToAdjacencyEdges(reference);

    /*
     * Check the matching we have.
     */
    assert(stList_length(chosenEdges) * 2 == nodeNumber);
    stList *nodes = stHash_getValues(endsToNodes);
    assert(nodeNumber == stList_length(nodes));
    stSortedSet *nodesSet = stList_getSortedSet(nodes, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    assert(nodeNumber == stSortedSet_size(nodesSet));
    checkEdges(chosenEdges, nodesSet, 1, 0);
    stSortedSet_destruct(nodesSet);
    stList_destruct(nodes);

    /*
     * Add the reference genome into flower
     */
    addLinkAdjacenciesAndSegments(flower, referenceEvent);
    addTangleAdjacenciesAndSegments(flower, chosenEdges, nodesToEnds,
            referenceEvent);

    /*
     * Ensure the newly created ends have a group.
     */
    assignGroups(newEnds, flower, referenceEvent);

    /*
     * Cleanup
     */
    free(z);
    stList_destruct(newEnds);
    stList_destruct(stubEdges);
    stHash_destruct(endsToNodes);
    stHash_destruct(nodesToEnds);
    stList_destruct(chainEdges);
    stList_destruct(chosenEdges);
    stList_destruct(reference);
}
