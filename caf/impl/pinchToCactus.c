#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

///////////////////////////////////////////////////////////////////////////
// Construct dead end component
///////////////////////////////////////////////////////////////////////////

static void attachPinchBlockEndToAnotherComponent(stPinchBlock *pinchBlock, bool orientation, stList *anotherComponent,
        stHash *pinchEndsToAdjacencyComponents) {
    assert(pinchBlock != NULL);
    assert(stPinchBlock_getLength(pinchBlock) == 1);
    stPinchEnd staticPinchEnd = stPinchEnd_constructStatic(pinchBlock, orientation);
    stList *component = stHash_remove(pinchEndsToAdjacencyComponents, &staticPinchEnd);
    assert(component != NULL);
    assert(stList_length(component) == 1);
    stPinchEnd *pinchEnd = stList_get(component, 0);
    assert(stPinchEnd_equalsFn(&staticPinchEnd, pinchEnd));
    stList_setDestructor(component, NULL);
    stList_destruct(component);
    stList_append(anotherComponent, pinchEnd);
    stHash_insert(pinchEndsToAdjacencyComponents, pinchEnd, anotherComponent);
}

static stList *stCaf_constructDeadEndComponent(Flower *flower, stPinchThreadSet *threadSet, stHash *pinchEndsToAdjacencyComponents) {
    /*
     * Locates the ends of all the attached ends and merges together their 'dead end' components to create a single
     * 'dead end' component, as described in the JCB cactus paper.
     */
    //For each block end at the end of a thread, attach to dead end component if associated end if attached
    stList *deadEndAdjacencyComponent = stList_construct3(0, (void(*)(void *)) stPinchEnd_destruct);
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *pinchThread;
    while ((pinchThread = stPinchThreadSetIt_getNext(&threadIt))) {
        Cap *cap = flower_getCap(flower, stPinchThread_getName(pinchThread));
        assert(cap != NULL);
        End *end1 = cap_getEnd(cap), *end2 = cap_getEnd(cap_getAdjacency(cap));
        assert(end1 != NULL && end2 != NULL);
        if (end_isAttached(end1)) {
            stPinchSegment *pinchSegment = stPinchThread_getFirst(pinchThread);
            stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
            if (stPinchBlock_getFirst(pinchBlock) == pinchSegment) { //We only want to do this once
                attachPinchBlockEndToAnotherComponent(pinchBlock, stPinchSegment_getBlockOrientation(pinchSegment),
                        deadEndAdjacencyComponent, pinchEndsToAdjacencyComponents);
            }
        }
        if (end_isAttached(end2)) {
            stPinchSegment *pinchSegment = stPinchThread_getLast(pinchThread);
            stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
            if (stPinchBlock_getFirst(pinchBlock) == pinchSegment) { //And only once for the other end
                attachPinchBlockEndToAnotherComponent(pinchBlock, !stPinchSegment_getBlockOrientation(pinchSegment),
                        deadEndAdjacencyComponent, pinchEndsToAdjacencyComponents);
            }
        }
    }
    return deadEndAdjacencyComponent;
}

///////////////////////////////////////////////////////////////////////////
// Attach unatttached thread components
///////////////////////////////////////////////////////////////////////////

static bool threadIsAttachedToDeadEndComponent5Prime(stPinchThread *thread, stList *deadEndComponent, stHash *pinchEndsToAdjacencyComponents) {
    stPinchSegment *pinchSegment = stPinchThread_getFirst(thread);
    stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
    assert(pinchBlock != NULL);
    stPinchEnd staticPinchEnd = stPinchEnd_constructStatic(pinchBlock, stPinchSegment_getBlockOrientation(pinchSegment));
    stList *adjacencyComponent = stHash_search(pinchEndsToAdjacencyComponents, &staticPinchEnd);
    assert(adjacencyComponent != NULL);
    return adjacencyComponent == deadEndComponent;
}

static bool threadIsAttachedToDeadEndComponent3Prime(stPinchThread *thread, stList *deadEndComponent, stHash *pinchEndsToAdjacencyComponents) {
    stPinchSegment *pinchSegment = stPinchThread_getLast(thread);
    stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
    assert(pinchBlock != NULL);
    stPinchEnd staticPinchEnd = stPinchEnd_constructStatic(pinchBlock, !stPinchSegment_getBlockOrientation(pinchSegment));
    stList *adjacencyComponent = stHash_search(pinchEndsToAdjacencyComponents, &staticPinchEnd);
    assert(adjacencyComponent != NULL);
    return adjacencyComponent == deadEndComponent;
}

static bool threadComponentIsAttachedToDeadEndComponent(stList *threadComponent, stList *deadEndComponent,
        stHash *pinchEndsToAdjacencyComponents) {
    for (int32_t i = 0; i < stList_length(threadComponent); i++) {
        stPinchThread *thread = stList_get(threadComponent, i);
        if (threadIsAttachedToDeadEndComponent5Prime(thread, deadEndComponent, pinchEndsToAdjacencyComponents) ||
                threadIsAttachedToDeadEndComponent3Prime(thread, deadEndComponent, pinchEndsToAdjacencyComponents)) {
            return 1;
        }
    }
    return 0;
}

static void attachThreadToDeadEndComponent(stPinchThread *thread, stList *deadEndAdjacencyComponent, stHash *pinchEndsToAdjacencyComponents) {
    stPinchSegment *segment = stPinchThread_getFirst(thread);
    attachPinchBlockEndToAnotherComponent(stPinchSegment_getBlock(segment), stPinchSegment_getBlockOrientation(segment),
            deadEndAdjacencyComponent, pinchEndsToAdjacencyComponents);
    segment = stPinchThread_getLast(thread);
    attachPinchBlockEndToAnotherComponent(stPinchSegment_getBlock(segment), !stPinchSegment_getBlockOrientation(segment),
            deadEndAdjacencyComponent, pinchEndsToAdjacencyComponents);
}

static void attachThreadComponentToDeadEndComponent(stList *threadComponent, stList *deadEndComponent,
        stHash *pinchEndsToAdjacencyComponents, bool markEndsAttached, Flower *flower) {
    stPinchThread *longestPinchThread = NULL;
    int64_t maxLength = 0;
    for (int32_t i = 0; i < stList_length(threadComponent); i++) {
        stPinchThread *pinchThread = stList_get(threadComponent, i);
        if (stPinchThread_getLength(pinchThread) > maxLength) {
            longestPinchThread = pinchThread;
            maxLength = stPinchThread_getLength(pinchThread);
        }
    }
    assert(longestPinchThread != NULL);
    attachThreadToDeadEndComponent(longestPinchThread, deadEndComponent, pinchEndsToAdjacencyComponents);
    if (markEndsAttached) { //Get the ends and attach them
        Cap *cap = flower_getCap(flower, stPinchSegment_getName(stPinchThread_getFirst(longestPinchThread))); //The following three lines isolates the sequence associated with a segment.
        assert(cap != NULL);
        end_makeAttached(cap_getEnd(cap));
        end_makeAttached(cap_getEnd(cap_getAdjacency(cap)));
    }
}

static void stCaf_attachUnattachedThreadComponents(Flower *flower, stPinchThreadSet *threadSet, stList *deadEndComponent,
        stHash *pinchEndsToAdjacencyComponents, bool markEndsAttached) {
    /*
     * Locates threads components which have no dead ends part of the dead end component, and then
     * connects them, picking the longest thread to attach them.
     */
    stSortedSet *threadComponents = stPinchThreadSet_getThreadComponents(threadSet);
    assert(stSortedSet_size(threadComponents) > 0);
    stSortedSetIterator *threadIt = stSortedSet_getIterator(threadComponents);
    stList *threadComponent;
    while ((threadComponent = stSortedSet_getNext(threadIt)) != NULL) {
        if (!threadComponentIsAttachedToDeadEndComponent(threadComponent, deadEndComponent,
                pinchEndsToAdjacencyComponents)) {
            attachThreadComponentToDeadEndComponent(threadComponent, deadEndComponent, pinchEndsToAdjacencyComponents,
                    markEndsAttached, flower);
        }
    }
    stSortedSet_destructIterator(threadIt);
    stSortedSet_destruct(threadComponents);
}

///////////////////////////////////////////////////////////////////////////
// Create a cactus graph from a pinch graph
///////////////////////////////////////////////////////////////////////////

static stCactusNode *getCactusNode(stPinchEnd *pinchEnd, stHash *pinchEndsToAdjacencyComponents, stHash *adjacencyComponentsToCactusNodes) {
    stList *adjacencyComponent = stHash_search(pinchEndsToAdjacencyComponents, pinchEnd);
    assert(adjacencyComponent != NULL);
    stCactusNode *cactusNode = stHash_search(adjacencyComponentsToCactusNodes, adjacencyComponent);
    assert(cactusNode != NULL);
    return cactusNode;
}

void *stCaf_mergeNodeObjects(void *a, void *b) {
    stList *adjacencyComponents1 = a;
    stList *adjacencyComponents2 = b;
    assert(adjacencyComponents1 != adjacencyComponents2);
    stList_appendAll(adjacencyComponents1, adjacencyComponents2);
    stList_setDestructor(adjacencyComponents2, NULL);
    stList_destruct(adjacencyComponents2);
    return adjacencyComponents1;
}

static void *makeNodeObject(stList *adjacencyComponent) {
    stList *adjacencyComponents = stList_construct3(0, (void(*)(void *)) stList_destruct);
    stList_append(adjacencyComponents, adjacencyComponent);
    return adjacencyComponents;
}

static bool isDeadEndStubComponent(stList *adjacencyComponent, stPinchEnd *pinchEnd) {
    if(stList_length(adjacencyComponent) != 1) {
        return 0;
    }
    stPinchSegment *pinchSegment = stPinchBlock_getFirst(stPinchEnd_getBlock(pinchEnd));
    return (stPinchEnd_traverse5Prime(stPinchEnd_getOrientation(pinchEnd), pinchSegment) ? stPinchSegment_get5Prime(pinchSegment) :
            stPinchSegment_get3Prime(pinchSegment)) == NULL;
}

static stCactusGraph *stCaf_constructCactusGraph(stList *deadEndComponent, stHash *pinchEndsToAdjacencyComponents, stCactusNode **startCactusNode) {
    /*
     * Constructs a cactus graph from a set of pinch graph components, including the dead end component. Returns a cactus
     * graph, and assigns 'startCactusNode' to the cactus node containing the dead end component.
     */
    stCactusGraph *cactusGraph = stCactusGraph_construct2((void(*)(void *)) stList_destruct, NULL);
    stHash *adjacencyComponentsToCactusNodes = stHash_construct();
    stHash *pinchEndsToNonStaticPinchEnds = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, NULL, NULL);

    //Make the nodes
    *startCactusNode = stCactusNode_construct(cactusGraph, makeNodeObject(deadEndComponent));
    stHash_insert(adjacencyComponentsToCactusNodes, deadEndComponent, *startCactusNode);
    stHashIterator *pinchEndIt = stHash_getIterator(pinchEndsToAdjacencyComponents);
    stPinchEnd *pinchEnd;
    while ((pinchEnd = stHash_getNext(pinchEndIt)) != NULL) {
        //Make a map of edge ends to themselves for memory lookup
        stHash_insert(pinchEndsToNonStaticPinchEnds, pinchEnd, pinchEnd);
        stList *adjacencyComponent = stHash_search(pinchEndsToAdjacencyComponents, pinchEnd);
        assert(adjacencyComponent != NULL);
        if (stHash_search(adjacencyComponentsToCactusNodes, adjacencyComponent) == NULL) {
            if (isDeadEndStubComponent(adjacencyComponent, pinchEnd)) { //Going to be a bridge to nowhere, so we join it - this ensures
                //that all dead end nodes of free stubs end up in the same node as their non-dead end counterparts.
                assert(pinchEnd == stList_get(adjacencyComponent, 0));
                stPinchEnd otherPinchEnd = stPinchEnd_constructStatic(stPinchEnd_getBlock(pinchEnd), !stPinchEnd_getOrientation(pinchEnd));
                stList *otherAdjacencyComponent = stHash_search(pinchEndsToAdjacencyComponents, &otherPinchEnd);
                assert(otherAdjacencyComponent != NULL && adjacencyComponent != otherAdjacencyComponent);
                stCactusNode *cactusNode = stHash_search(adjacencyComponentsToCactusNodes, otherAdjacencyComponent);
                if (cactusNode == NULL) {
                    cactusNode = stCactusNode_construct(cactusGraph, makeNodeObject(otherAdjacencyComponent));
                    stHash_insert(adjacencyComponentsToCactusNodes, otherAdjacencyComponent, cactusNode);
                }
                stList *adjacencyComponents = stCactusNode_getObject(cactusNode);
                stList_append(adjacencyComponents, adjacencyComponent);
                stHash_insert(adjacencyComponentsToCactusNodes, adjacencyComponent, cactusNode);
            } else {
                stHash_insert(adjacencyComponentsToCactusNodes, adjacencyComponent,
                        stCactusNode_construct(cactusGraph, makeNodeObject(adjacencyComponent)));
            }
        }
    }
    stHash_destructIterator(pinchEndIt);

    //Make the edges
    pinchEndIt = stHash_getIterator(pinchEndsToAdjacencyComponents);
    while ((pinchEnd = stHash_getNext(pinchEndIt)) != NULL) {
        if (stPinchEnd_getOrientation(pinchEnd)) { //Assure we make the edge only once
            assert(stPinchEnd_getBlock(pinchEnd) != NULL);
            stPinchEnd pinchEnd2Static = stPinchEnd_constructStatic(stPinchEnd_getBlock(pinchEnd), 0);
            stPinchEnd *pinchEnd2 = stHash_search(pinchEndsToNonStaticPinchEnds, &pinchEnd2Static);
            assert(pinchEnd2 != NULL);
            assert(pinchEnd != pinchEnd2);
            stCactusNode *cactusNode1 = getCactusNode(pinchEnd, pinchEndsToAdjacencyComponents, adjacencyComponentsToCactusNodes);
            stCactusNode *cactusNode2 = getCactusNode(pinchEnd2, pinchEndsToAdjacencyComponents, adjacencyComponentsToCactusNodes);
            stCactusEdgeEnd_construct(cactusGraph, cactusNode1, cactusNode2, pinchEnd, pinchEnd2);
        }
    }
    stHash_destructIterator(pinchEndIt);
    stHash_destruct(pinchEndsToNonStaticPinchEnds);
    stHash_destruct(adjacencyComponentsToCactusNodes);

    //Run the cactus-ifying functions
    stCactusGraph_collapseToCactus(cactusGraph, stCaf_mergeNodeObjects, *startCactusNode);
    stCactusGraph_collapseBridges(cactusGraph, *startCactusNode, stCaf_mergeNodeObjects);

    return cactusGraph;
}

///////////////////////////////////////////////////////////////////////////
// Function that draws together above functions to generate a cactus graph from a pinch graph.
///////////////////////////////////////////////////////////////////////////

stCactusGraph *stCaf_getCactusGraphForThreadSet(Flower *flower, stPinchThreadSet *threadSet, stCactusNode **startCactusNode,
        stList **deadEndComponent, bool attachEndsInFlower) {
    //Get adjacency components
    stHash *pinchEndsToAdjacencyComponents;
    stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents2(threadSet, &pinchEndsToAdjacencyComponents);
    stList_setDestructor(adjacencyComponents, NULL);
    stList_destruct(adjacencyComponents);

    //Merge together dead end component
    *deadEndComponent = stCaf_constructDeadEndComponent(flower, threadSet, pinchEndsToAdjacencyComponents);

    //Join unattached components of graph by dead ends to dead end component, and make other ends 'attached' if necessary
    stCaf_attachUnattachedThreadComponents(flower, threadSet, *deadEndComponent, pinchEndsToAdjacencyComponents, attachEndsInFlower);

    //Create cactus
    stCactusGraph *cactusGraph = stCaf_constructCactusGraph(*deadEndComponent, pinchEndsToAdjacencyComponents, startCactusNode);

    //Cleanup (the memory is owned by the cactus graph, so this does not break anything)
    stHash_destruct(pinchEndsToAdjacencyComponents);

    return cactusGraph;
}
