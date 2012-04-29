#include "sonLib.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "cactus.h"
#include "stCaf.h"

stSortedSet *stList_convertToSortedSet(stList *list) {
    stSortedSet *set = stList_getSortedSet(list, NULL);
    stSortedSet_setDestructor(set, stList_getDestructor(list));
    stList_setDestructor(list, NULL);
    stList_destruct(list);
    return set;
}

static void makeBlock(End *end, stSegment *segment, bool orientation, stHash *endsToBlocks) {
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
        while ((cap = end_getNext(end)) != NULL) {
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
                makeBlock(end, _5PrimeSegment, 1, endsToBlocks);
                stSegment *_3PrimeSegment = stSegment_get3Prime(_5PrimeSegment);
                stSegment_split(_3PrimeSegment, cap_getCoordinate(adjacentCap));
                makeBlock(end, _3PrimeSegment, 0, endsToBlocks);
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
    while ((cap = end_getNext(end)) != NULL) {
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
            stCaf_constructDeadEndComponentP(end, deadEndAdjacencyComponent, adjacencyComponents, edgeEndsToAdjacencyComponents);
        }
    }
    flower_destructEndIterator(endIt);
    return deadEndAdjacencyComponent;
}

bool threadIsAttachedToDeadEndComponent(stThread *thread, stList *deadEndComponent, stHash *edgeEndsToAdjacencyComponents) {
    stThread *thread = stList_get(threadComponent, j);
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
    stList *threadComponents = stThreadSet_getThreadComponents(threadSet);
    if (stList_length(threadComponents) > 1) {
        for (int32_t i = 0; i < stList_length(threadComponents); i++) {
            stList *threadComponent = stList_get(threadComponents, i);
            if (!threadComponentIsAttachedToDeadEndComponent(threadComponent)) {
                attachThreadComponentToDeadEndComponent(threadComponent, deadEndComponent, adjacencyComponents,
                        edgeEndsToAdjacencyComponents);
            }
        }
    }
    stList_destruct(threadComponents);
}

stList *getCactusNode(stEnd *end, stHash *edgeEndsToAdjacencyComponents, stHash *adjacencyComponentsToCactusNodes) {
    stList *adjacencyComponent = stHash_search(edgeEndsToAdjacencyComponents, end);
    assert(adjacencyComponent != NULL);
    stCactusNode *cactusNode = stHash_search(adjacencyComponentsToCactusNodes, adjacencyComponent);
    assert(cactusNode != NULL);
    return cactusNode;
}

static void *mergeNodeObjects(const void *a, const void *b) {
    stList *groups1 = a;
    stList *groups2 = a;
    assert(groups1 != group2);
    stList_appendAll(groups1, groups2);
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
        stEnd *otherEnd = stHash_search(endsToNonStaticEnds, &stEnd_constructStatic(stEnd_getBlock(end), !stEnd_getOrientation(end)));
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

void makeBlock(stCactusEdgeEnd *cactusEdgeEnd, Flower *flower, stHash *edgeEndsToEnds) {

}

void setAdjacencies();

End *convertEdgeEndToEnd(stCactusEdgeEnd *cactusEdgeEnd) {
    return stHash_search(edgeEndsToEnds, stCactusEdgeEnd_getObject(cactusEdgeEnd));
}

void stCaf_convertCactusGraphToFlowers(stCactusGraph *cactusGraph, stCactusNode *startCactusNode, Flower *flower, stList *deadEndComponent) {
    stList *stack = stList_construct();
    stList_append(stack, startCactusNode);
    stList_append(stack, flower);
    stHash *edgeEndsToEnds = stHash_construct();
    while (stList_length(stack) > 0) {
        flower = stList_pop(stack);
        stCactusNode *cactusNode = stList_pop(stack);
        //Make the chains
        stCactusNode_edgeEndIt cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
        stCactusEdgeEnd *cactusEdgeEnd;
        while ((cactusEdgeEnd = stCactusNode_edgeEndIt_getNext(&cactusEdgeEndIt))) {
            if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) { //Iterate around the chain
                if (stHash_search(edgeEndsToEnds, stCactusEdgeEnd_getObject(cactusEdgeEnd)) == NULL) { //Make first block
                    makeBlock(cactusEdgeEnd, flower, edgeEndsToEnds);
                }
                cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd);
                if (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd)) { //We have a chain
                    Chain *chain = chain_construct(flower);
                    do {
                        stCactusEdgeEnd *linkedCactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd);
                        if (stHash_search(edgeEndsToEnds, stCactusEdgeEnd_getObject(linkedCactusEdgeEnd)) == NULL) { //Make subsequent block
                            makeBlock(linkedCactusEdgeEnd, flower, edgeEndsToEnds);
                        }
                        Group *group = group_construct2(flower);
                        End *end1 = stHash_search(edgeEndsToEnds, stCactusEdgeEnd_getObject(cactusEdgeEnd));
                        End *end2 = stHash_search(edgeEndsToEnds, stCactusEdgeEnd_getObject(linkedCactusEdgeEnd));
                        end_setGroup(end1, group);
                        end_setGroup(end2, group);
                        link_construct(end1, end2, group, chain);
                        stList_append(stack, stCactusEdgeEnd_getNode(cactusEdgeEnd));
                        stList_append(stack, group_makeNestedFlower(group));
                        setAdjacencies(group);
                        cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(linkedCactusEdgeEnd);
                    } while(!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd));
                }
            }
        }
        //Now construct the tangle groups
        stList *adjacencyComponents = stCactusNode_getObject(cactusNode);
        for (int32_t i = 0; i < stList_length(adjacencyComponents); i++) {
            stList *adjacencyComponent = stList_get(adjacencyComponents, i);
            if (adjacencyComponent != deadEndComponent) {
                assert(stList_length(adjacencyComponent) != 2);
                Group *group = group_construct2(flower);
                for (int32_t j = 0; j < stList_length(adjacencyComponent); j++) {
                    stEdgeEnd *edgeEnd = stList_get(adjacencyComponent, j);
                    End *end = stHash_search(edgeEndsToEnds, edgeEnd);
                    assert(edgeEnd != NULL);
                    end_setGroup(end, group);
                }
                //Now set the adjacencies
                setAdjacencies(group);
            }
        }
    }
    stList_destruct(stack);
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
    stCaf_attachUnattachedThreadComponents(flower, deadEndComponent, adjacencyComponents, edgeEndsToAdjacencyComponents, 1);

    //Create cactus
    stCactusNode *startCactusNode;
    stCactusGraph *cactusGraph = stCaf_constructCactusGraph(adjacencyComponents, deadEndComponent, edgeEndsToAdjacencyComponents,
            &startCactusNode);

    //Convert cactus graph/pinch graph to API
    stCaf_convertCactusGraphToFlowers(cactusGraph, startCactusNode, flower);

    //Cleanup
    stCactusGraph_destruct(cactusGraph);
    stHash_destruct(edgeEndsToAdjacencyComponents);
    stSortedSet_destruct(adjacencyComponents);
    stThreadSet_destruct(threadSet);
}
