#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

static void stCaf_constructEmptyPinchGraphP(End *end, stSegment *segment, bool orientation, stHash *endsToBlocks) {
    end = end_getPositiveOrientation(end);
    stBlock *block = stHash_search(endsToBlocks, end);
    block = block == NULL ? stBlock_construct3(segment, orientation) : stBlock_pinch2(block, segment, orientation);
    stHash_insert(endsToBlocks, end, block);
}

stThreadSet *stCaf_constructEmptyPinchGraph(Flower *flower) {
    stThreadSet *threadSet = stThreadSet_construct();
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    stHash *endsToBlocks = stHash_construct();
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        assert(!end_isBlockEnd(end));
        End_InstanceIterator *capIt = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(capIt)) != NULL) {
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            if (!cap_getSide(cap)) {
                Cap *adjacentCap = cap_getAdjacency(cap);
                assert(adjacentCap != NULL);
                assert(cap_getSide(adjacentCap));
                assert(cap_getCoordinate(cap) < cap_getCoordinate(adjacentCap));
                stThread *thread = stThreadSet_addThread(threadSet, cap_getName(cap), cap_getCoordinate(cap),
                        cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) + 1);
                stSegment *_5PrimeSegment = stThread_getFirst(thread);
                stSegment_split(_5PrimeSegment, cap_getCoordinate(cap));
                stCaf_constructEmptyPinchGraphP(end, _5PrimeSegment, 1, endsToBlocks);
                stSegment *_3PrimeSegment = stSegment_get3Prime(_5PrimeSegment);
                stSegment_split(_3PrimeSegment, cap_getCoordinate(adjacentCap));
                stCaf_constructEmptyPinchGraphP(end, _3PrimeSegment, 0, endsToBlocks);
            }
        }
        end_destructInstanceIterator(capIt);
    }
    flower_destructEndIterator(endIt);
    stHash_destruct(endsToBlocks);
    return threadSet;
}

void stCaf_addAlignmentsToPinchGraph(stThreadSet *threadSet, stPinch *(*pinchGenerator)()) {
    stPinch *pinch;
    while ((pinch = pinchGenerator()) != NULL) {
        stThread *thread1 = stThreadSet_getThread(threadSet, pinch->name1);
        stThread *thread2 = stThreadSet_getThread(threadSet, pinch->name2);
        assert(thread1 != NULL && thread2 != NULL);
        stThread_pinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand);
    }
    stThreadSet_joinTrivialBoundaries(threadSet);
}

void attachThreadToDeadEndComponentP(stBlock *block, bool orientation, stList *deadEndAdjacencyComponent, stSortedSet *adjacencyComponents,
        stHash *edgeEndsToAdjacencyComponents) {
    assert(block != NULL);
    stEnd staticEnd = stEnd_constructStatic(block, orientation);
    stList *component = stHash_remove(edgeEndsToAdjacencyComponents, &staticEnd);
    assert(component != NULL);
    assert(stList_length(component) == 1);
    stEnd *end = stList_get(component, 0);
    stList_destruct(component);
    stList_append(deadEndAdjacencyComponent, end);
    stHash_insert(edgeEndsToAdjacencyComponents, end, deadEndAdjacencyComponent);
}

void attachThreadToDeadEndComponent(stThread *thread, stList *deadEndAdjacencyComponent, stSortedSet *adjacencyComponents,
        stHash *edgeEndsToAdjacencyComponents) {
    attachThreadToDeadEndComponentP(stSegment_getBlock(stThread_getFirst(thread)), 1, deadEndAdjacencyComponent, adjacencyComponents,
            edgeEndsToAdjacencyComponents);
    attachThreadToDeadEndComponentP(stSegment_getBlock(stThread_getLast(thread)), 0, deadEndAdjacencyComponent, adjacencyComponents,
            edgeEndsToAdjacencyComponents);
}

void stCaf_constructDeadEndComponentP(End *end, stThreadSet *threadSet, stList *deadEndAdjacencyComponent,
        stSortedSet *adjacencyComponents, stHash *edgeEndsToAdjacencyComponents) {
    End_InstanceIterator *capIt = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(capIt)) != NULL) {
        cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
        if (!cap_getSide(cap)) {
            end_destructInstanceIterator(capIt);
            stThread *thread = stThreadSet_getThread(threadSet, cap_getName(cap));
            assert(thread != NULL);
            attachThreadToDeadEndComponent(thread, deadEndAdjacencyComponent, adjacencyComponents, edgeEndsToAdjacencyComponents);
            return;
        }
    }
    end_destructInstanceIterator(capIt);
    assert(0);
}

stList *stCaf_constructDeadEndComponent(Flower *flower, stThreadSet *threadSet, stSortedSet *adjacencyComponents,
        stHash *edgeEndsToAdjacencyComponents) {
    stList *deadEndAdjacencyComponent = stList_construct();
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        assert(!end_isBlockEnd(end));
        if (end_isAttached(end)) {
            stCaf_constructDeadEndComponentP(end, threadSet, deadEndAdjacencyComponent, adjacencyComponents, edgeEndsToAdjacencyComponents);
        }
    }
    flower_destructEndIterator(endIt);
    return deadEndAdjacencyComponent;
}

bool threadIsAttachedToDeadEndComponent(stThread *thread, stList *deadEndComponent, stHash *edgeEndsToAdjacencyComponents) {
    stBlock *block = stSegment_getBlock(stThread_getFirst(thread));
    assert(block != NULL);
    stEnd staticEnd = stEnd_constructStatic(block, 1);
    stList *adjacencyComponent = stHash_search(edgeEndsToAdjacencyComponents, &staticEnd);
    assert(adjacencyComponent != NULL);
    return adjacencyComponent == deadEndComponent;
}

bool threadComponentIsAttachedToDeadEndComponent(stList *threadComponent, stList *deadEndComponent, stHash *edgeEndsToAdjacencyComponents) {
    for (int32_t i = 0; i < stList_length(threadComponent); i++) {
        if (threadIsAttachedToDeadEndComponent(stList_get(threadComponent, i), deadEndComponent, edgeEndsToAdjacencyComponents)) {
            return 1;
        }
    }
    return 0;
}

void attachThreadComponentToDeadEndComponent(stList *threadComponent, stList *deadEndComponent, stSortedSet *adjacencyComponents,
        stHash *edgeEndsToAdjacencyComponents) {
    stThread *longestThread = NULL;
    int64_t maxLength = 0;
    for (int32_t i = 0; i < stList_length(threadComponent); i++) {
        stThread *thread = stList_get(threadComponent, i);
        if (stThread_getLength(thread) > maxLength) {
            longestThread = thread;
            maxLength = stThread_getLength(thread);
        }
    }
    assert(longestThread != NULL);
    attachThreadToDeadEndComponent(longestThread, deadEndComponent, adjacencyComponents, edgeEndsToAdjacencyComponents);
}

void stCaf_attachUnattachedThreadComponents(Flower *flower, stThreadSet *threadSet, stList *deadEndComponent,
        stSortedSet *adjacencyComponents, stHash *edgeEndsToAdjacencyComponents, bool markEndsAttached) {
    stSortedSet *threadComponents = stThreadSet_getThreadComponents(threadSet);
    if (stSortedSet_size(threadComponents) > 1) {
        stSortedSetIterator *threadIt = stSortedSet_getIterator(threadComponents);
        stList *threadComponent;
        while((threadComponent = stSortedSet_getNext(threadIt)) != NULL) {
            if (!threadComponentIsAttachedToDeadEndComponent(threadComponent, deadEndComponent, edgeEndsToAdjacencyComponents)) {
                attachThreadComponentToDeadEndComponent(threadComponent, deadEndComponent, adjacencyComponents,
                        edgeEndsToAdjacencyComponents);
            }
        }
        stSortedSet_destructIterator(threadIt);
    }
    stSortedSet_destruct(threadComponents);
}

stCactusNode *getCactusNode(stEnd *end, stHash *edgeEndsToAdjacencyComponents, stHash *adjacencyComponentsToCactusNodes) {
    stList *adjacencyComponent = stHash_search(edgeEndsToAdjacencyComponents, end);
    assert(adjacencyComponent != NULL);
    stCactusNode *cactusNode = stHash_search(adjacencyComponentsToCactusNodes, adjacencyComponent);
    assert(cactusNode != NULL);
    return cactusNode;
}

static void *mergeNodeObjects(void *a, void *b) {
    stList *groups1 = a;
    stList *groups2 = b;
    assert(groups1 != groups2);
    stList_appendAll(groups1, groups2);
    stList_setDestructor(groups2, NULL);
    stList_destruct(groups2);
    return groups1;
}

static void *makeNodeObject(stList *adjacencyComponent) {
    stList *groups = stList_construct();
    stList_append(groups, adjacencyComponent);
    return groups;
}

stCactusGraph *stCaf_constructCactusGraph(stSortedSet *adjacencyComponents, stList *deadEndComponent,
        stHash *edgeEndsToAdjacencyComponents, stCactusNode **startCactusNode) {
    stCactusGraph *cactusGraph = stCactusGraph_construct();
    stHash *adjacencyComponentsToCactusNodes = stHash_construct();
    //Make the nodes
    *startCactusNode = stCactusNode_construct(cactusGraph, makeNodeObject(deadEndComponent));
    stHash_insert(adjacencyComponentsToCactusNodes, deadEndComponent, *startCactusNode);
    stSortedSetIterator *it = stSortedSet_getIterator(adjacencyComponents);
    stList *adjacencyComponent;
    while ((adjacencyComponent = stSortedSet_getNext(it)) != NULL) {
        if (adjacencyComponent != deadEndComponent) {
            stHash_insert(adjacencyComponentsToCactusNodes, adjacencyComponent,
                    stCactusNode_construct(cactusGraph, makeNodeObject(adjacencyComponent)));
        }
    }
    stSortedSet_destructIterator(it);

    //Make a map of edge ends to themselves for memory lookup
    stHash *endsToNonStaticEnds = stHash_construct();
    stHashIterator *edgeEndIt = stHash_getIterator(edgeEndsToAdjacencyComponents);
    stEnd *end;
    while ((end = stHash_getNext(edgeEndIt)) != NULL) {
        stHash_insert(endsToNonStaticEnds, end, end);
    }
    stHash_destructIterator(edgeEndIt);
    //Make the edges
    edgeEndIt = stHash_getIterator(edgeEndsToAdjacencyComponents);
    while ((end = stHash_getNext(edgeEndIt)) != NULL) {
        stCactusNode *cactusNode1 = getCactusNode(end, edgeEndsToAdjacencyComponents, adjacencyComponentsToCactusNodes);
        stEnd end2 = stEnd_constructStatic(stEnd_getBlock(end), !stEnd_getOrientation(end));
        stEnd *otherEnd = stHash_search(endsToNonStaticEnds, &end2);
        assert(otherEnd != NULL);
        stCactusNode *cactusNode2 = getCactusNode(otherEnd, edgeEndsToAdjacencyComponents, adjacencyComponentsToCactusNodes);
        stCactusEdgeEnd_construct(cactusGraph, cactusNode1, cactusNode2, end, otherEnd);
    }
    stHash_destructIterator(edgeEndIt);
    stHash_destruct(endsToNonStaticEnds);
    stHash_destruct(adjacencyComponentsToCactusNodes);

    //Run the cactus-ifying functions
    stCactusGraph_collapseToCactus(cactusGraph, mergeNodeObjects, *startCactusNode);
    stCactusGraph_markCycles(cactusGraph, *startCactusNode);
    stCactusGraph_collapseBridges(cactusGraph, *startCactusNode, mergeNodeObjects);

    return cactusGraph;
}

static void makeBlockP(stEnd *pinchBlockEnd, End *end, stHash *pinchBlockEndsToEnds) {
    assert(stHash_search(pinchBlockEndsToEnds, pinchBlockEnd) == NULL);
    stHash_insert(pinchBlockEndsToEnds, stEnd_construct(stEnd_getBlock(pinchBlockEnd), stEnd_getOrientation(pinchBlockEnd)), end);
}

static void makeBlock(stCactusEdgeEnd *cactusEdgeEnd, Flower *parentFlower, Flower *flower, stHash *pinchBlockEndsToEnds) {
    stEnd *pinchBlockEnd = stCactusEdgeEnd_getObject(cactusEdgeEnd);
    assert(pinchBlockEnd != NULL);
    stBlock *pinchBlock = stEnd_getBlock(pinchBlockEnd);
    Block *block = block_construct(stBlock_getLength(pinchBlock), flower);
    stSegment *pinchSegment;
    stBlockIt pinchSegmentIt = stBlock_getSegmentIterator(pinchBlock);
    while ((pinchSegment = stBlockIt_getNext(&pinchSegmentIt))) {
        Cap *cap = flower_getCap(parentFlower, stSegment_getName(pinchSegment));
        assert(cap != NULL);
        Sequence *sequence = cap_getSequence(cap);
        segment_construct2(
                stEnd_getOrientation(pinchBlockEnd) ^ stSegment_getBlockOrientation(pinchSegment) ? block : block_getReverse(block),
                stSegment_getStart(pinchSegment), 1, sequence);
    }
    makeBlockP(pinchBlockEnd, block_get5End(block), pinchBlockEndsToEnds);
    stEnd *otherPinchBlockEnd = stCactusEdgeEnd_getObject(stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd));
    makeBlockP(otherPinchBlockEnd, block_get3End(block), pinchBlockEndsToEnds);
}

End *convertPinchBlockEndToEnd(stEnd *pinchBlockEnd, stHash *pinchBlockEndsToEnds, Flower *flower) {
    End *end = stHash_search(pinchBlockEndsToEnds, pinchBlockEnd);
    if (end == NULL) {
        return NULL;
    }
    End *end2 = flower_getEnd(flower, end_getName(end));
    if (end2 == NULL) {
        assert(end_isFree(end));
        assert(end_isStubEnd(end));
        //Copy the end down the hierarchy
        Group *parentGroup = flower_getParentGroup(flower);
        end2 = convertPinchBlockEndToEnd(pinchBlockEnd, pinchBlockEndsToEnds, group_getFlower(parentGroup));
        assert(end2 != NULL);
        if (end_getGroup(end2) == NULL) {
            end_setGroup(end2, parentGroup);
        }
        assert(parentGroup == end_getGroup(end2));
        end2 = end_copyConstruct(end2, flower);
        assert(end2 != NULL);
        assert(end_getFlower(end2) == flower);
    }
    return end2;
}

End *convertCactusEdgeEndToEnd(stCactusEdgeEnd *cactusEdgeEnd, stHash *pinchBlockEndsToEnds, Flower *flower) {
    return convertPinchBlockEndToEnd(stCactusEdgeEnd_getObject(cactusEdgeEnd), pinchBlockEndsToEnds, flower);
}

void getPinchBlockEndsToEndsHashPP(stBlock *pinchBlock, bool orientation, End *end, stHash *pinchBlockEndsToEnds) {
    stEnd pinchEnd = stEnd_constructStatic(pinchBlock, orientation);
    if (stHash_search(pinchBlockEndsToEnds, &pinchEnd) == NULL) {
        stHash_insert(pinchBlockEndsToEnds, stEnd_construct(pinchBlock, orientation), end);
    }
}

void getPinchBlockEndsToEndsHashP(stBlock *pinchBlock, Cap *cap, stHash *pinchBlockEndsToEnds) {
    assert(cap != NULL);
    End *end = cap_getEnd(cap);
    assert(end != NULL);
    getPinchBlockEndsToEndsHashPP(pinchBlock, 1, end, pinchBlockEndsToEnds);
    getPinchBlockEndsToEndsHashPP(pinchBlock, 0, end, pinchBlockEndsToEnds);
}

stHash *getPinchBlockEndsToEndsHash(stThreadSet *threadSet, Flower *parentFlower) {
    stHash *pinchBlockEndsToEnds = stHash_construct3(stEnd_hashFn, stEnd_equalsFn, (void(*)(void *)) stEnd_destruct, NULL);
    stThreadIt threadIt = stThreadSet_getIterator(threadSet);
    stThread *thread;
    while ((thread = stThreadIt_getNext(&threadIt))) {
        Cap *cap = flower_getCap(parentFlower, stThread_getName(thread));
        getPinchBlockEndsToEndsHashP(stSegment_getBlock(stThread_getFirst(thread)), cap, pinchBlockEndsToEnds);
        getPinchBlockEndsToEndsHashP(stSegment_getBlock(stThread_getLast(thread)), cap_getAdjacency(cap), pinchBlockEndsToEnds);
    }
    return pinchBlockEndsToEnds;
}

void makeChain(stCactusEdgeEnd *cactusEdgeEnd, Flower *flower, stHash *pinchBlockEndsToEnds, Flower *parentFlower, stList *stack) {
    cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd);
    if (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd)) { //We have a chain
        Chain *chain = chain_construct(flower);
        do {
            stCactusEdgeEnd *linkedCactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd);
            if (convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchBlockEndsToEnds, flower) == NULL) { //Make subsequent block
                makeBlock(linkedCactusEdgeEnd, parentFlower, flower, pinchBlockEndsToEnds);
            }
            assert(stCactusEdgeEnd_getNode(cactusEdgeEnd) == stCactusEdgeEnd_getNode(linkedCactusEdgeEnd));
            Group *group = group_construct2(flower);
            End *end1 = convertCactusEdgeEndToEnd(cactusEdgeEnd, pinchBlockEndsToEnds, flower);
            End *end2 = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchBlockEndsToEnds, flower);
            end_setGroup(end1, group);
            end_setGroup(end2, group);
            link_construct(end1, end2, group, chain);
            stList_append(stack, stCactusEdgeEnd_getNode(cactusEdgeEnd));
            stList_append(stack, group_makeNestedFlower(group));
            cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(linkedCactusEdgeEnd);
        } while (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd));
    }
}

void makeChains(stCactusNode *cactusNode, Flower *flower, stHash *pinchBlockEndsToEnds, Flower *parentFlower, stList *stack) {
    stCactusNode_edgeEndIt cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
    stCactusEdgeEnd *cactusEdgeEnd;
    while ((cactusEdgeEnd = stCactusNode_edgeEndIt_getNext(&cactusEdgeEndIt))) {
        if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) { //We have some sort of chain
            End *end = convertCactusEdgeEndToEnd(cactusEdgeEnd, pinchBlockEndsToEnds, flower);
            stCactusEdgeEnd *linkedCactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd), *startEnd = NULL;
            if (end != NULL) {
                End *end2;
                if ((end2 = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchBlockEndsToEnds, flower)) != NULL) {
                    assert(end_getSide(end) != end_getSide(end2));
                }
                startEnd = end_getSide(end) ? cactusEdgeEnd : linkedCactusEdgeEnd;
            } else {
                end = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchBlockEndsToEnds, flower);
                if (end != NULL) {
                    startEnd = end_getSide(end) ? linkedCactusEdgeEnd : cactusEdgeEnd;
                } else {
                    makeBlock(linkedCactusEdgeEnd, parentFlower, flower, pinchBlockEndsToEnds);
                    startEnd = linkedCactusEdgeEnd;
                }
            }
            assert(startEnd != NULL);
            makeChain(startEnd, flower, pinchBlockEndsToEnds, parentFlower, stack);
        }
    }
}

void makeTangles(stCactusNode *cactusNode, Flower *flower, stHash *pinchBlockEndsToEnds, stList *deadEndComponent) {
    stList *adjacencyComponents = stCactusNode_getObject(cactusNode);
    for (int32_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stList *adjacencyComponent = stList_get(adjacencyComponents, i);
        if (adjacencyComponent != deadEndComponent) {
            Group *group = group_construct2(flower);
            for (int32_t j = 0; j < stList_length(adjacencyComponent); j++) {
                stEnd *pinchBlockEnd = stList_get(adjacencyComponent, j);
                assert(pinchBlockEnd != NULL);
                End *end = convertPinchBlockEndToEnd(pinchBlockEnd, pinchBlockEndsToEnds, flower);
                assert(end != NULL);
                end_setGroup(end, group);
            }
        }
    }
}

static int addAdjacenciesToLeafCapsP(Cap *cap1, Cap *cap2) {
    assert(cap_getStrand(cap1) && cap_getStrand(cap2));
    Sequence *sequence1 = cap_getSequence(cap1);
    Sequence *sequence2 = cap_getSequence(cap2);
    int32_t i = cactusMisc_nameCompare(sequence_getName(sequence1), sequence_getName(sequence2));
    if (i == 0) {
        int32_t j = cap_getCoordinate(cap1);
        int32_t k = cap_getCoordinate(cap2);
        i = j - k;
        if (i == 0) {
            assert(cap_getSegment(cap1) == cap_getSegment(cap2));
            j = cap_getSide(cap1);
            k = cap_getSide(cap2);
            assert((j && !k) || (!j && k));
            i = j ? -1 : 1;
        }
    }
    return i;
}

static void addAdjacenciesToLeafCaps(Flower *flower) {
    /*
     * Build a list of caps, then sort them.
     */
    stList *list = stList_construct();
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(instanceIterator)) != NULL) {
            if (!cap_isInternal(cap)) {
                if (!cap_getStrand(cap)) {
                    cap = cap_getReverse(cap);
                }
                stList_append(list, cap);
            }
        }
        end_destructInstanceIterator(instanceIterator);
    }
    flower_destructEndIterator(endIterator);
    assert(stList_length(list) % 2 == 0);

    /*
     * Sort the caps.
     */
    stList_sort(list, (int (*)(const void *, const void *))addAdjacenciesToLeafCapsP);

    /*
     * Now make the adjacencies.
     */
    for (int32_t i = 1; i < stList_length(list); i += 2) {
        Cap *cap = stList_get(list, i-1);
        Cap *cap2 = stList_get(list, i);
        cap_makeAdjacent(cap, cap2);
    }

    /*
     * Clean up.
     */
    stList_destruct(list);
}

static void addAdjacenciesToEnds(Flower *flower) {
    /*
     * Add the adjacencies between leaf caps of a flower. Does not touch internal instances.
     * All ends must be in the flower before this function is called.
     */
    addAdjacenciesToLeafCaps(flower);
    Group_EndIterator *adjacencyIterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(adjacencyIterator)) != NULL) {
        if (group_getNestedFlower(group) != NULL) {
            addAdjacenciesToEnds(group_getNestedFlower(group));
        }
    }
    flower_destructGroupIterator(adjacencyIterator);
}

void stCaf_convertCactusGraphToFlowers(stThreadSet *threadSet, stCactusNode *startCactusNode, Flower *parentFlower,
        stList *deadEndComponent) {
    stList *stack = stList_construct();
    stList_append(stack, startCactusNode);
    stList_append(stack, parentFlower);
    stHash *pinchBlockEndsToEnds = getPinchBlockEndsToEndsHash(threadSet, parentFlower);
    while (stList_length(stack) > 0) {
        Flower *flower = stList_pop(stack);
        stCactusNode *cactusNode = stList_pop(stack);
        makeChains(cactusNode, flower, pinchBlockEndsToEnds, parentFlower, stack);
        makeTangles(cactusNode, flower, pinchBlockEndsToEnds, deadEndComponent);
    }
    stHash_destruct(pinchBlockEndsToEnds);
    stList_destruct(stack);
    addAdjacenciesToEnds(parentFlower);
}

void stCore(Flower *flower, stPinch *(*pinchGenerator)()) {
    //Create empty pinch graph from flower
    stThreadSet *threadSet = stCaf_constructEmptyPinchGraph(flower);

    //Add alignments to pinch graph
    stCaf_addAlignmentsToPinchGraph(threadSet, pinchGenerator);

    //Get adjacency components
    stHash *edgeEndsToAdjacencyComponents;
    stSortedSet *adjacencyComponents = stList_convertToSortedSet(
            stThreadSet_getAdjacencyComponents2(threadSet, &edgeEndsToAdjacencyComponents));

    //Merge together dead end component
    stList *deadEndComponent = stCaf_constructDeadEndComponent(flower, threadSet, adjacencyComponents, edgeEndsToAdjacencyComponents);

    //Join unattached components of graph by dead ends to dead end component, and make other ends 'attached' if necessary
    stCaf_attachUnattachedThreadComponents(flower, threadSet, deadEndComponent, adjacencyComponents, edgeEndsToAdjacencyComponents, 1);

    //Create cactus
    stCactusNode *startCactusNode;
    stCactusGraph *cactusGraph = stCaf_constructCactusGraph(adjacencyComponents, deadEndComponent, edgeEndsToAdjacencyComponents,
            &startCactusNode);

    //Convert cactus graph/pinch graph to API
    stCaf_convertCactusGraphToFlowers(threadSet, startCactusNode, flower, deadEndComponent);

    //Cleanup
    stCactusGraph_destruct(cactusGraph);
    stHash_destruct(edgeEndsToAdjacencyComponents);
    stSortedSet_destruct(adjacencyComponents);
    stThreadSet_destruct(threadSet);
}
