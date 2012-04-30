#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

///////////////////////////////////////////////////////////////////////////
// Construct a pinch graph from an empty cactus
///////////////////////////////////////////////////////////////////////////

static void stCaf_constructEmptyPinchGraphP(End *end, stPinchSegment *segment, bool orientation, stHash *endsToBlocks) {
    assert(stPinchSegment_getLength(segment) == 1);
    end = end_getPositiveOrientation(end);
    stPinchBlock *block = stHash_search(endsToBlocks, end);
    if (block == NULL) {
        stHash_insert(endsToBlocks, end, stPinchBlock_construct3(segment, orientation));
    } else {
        stPinchBlock_pinch2(block, segment, orientation);
    }
}

stPinchThreadSet *stCaf_constructEmptyPinchGraph(Flower *flower) {
    stPinchThreadSet *threadSet = stPinchThreadSet_construct();
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
                assert(cap_getSide(adjacentCap));
                assert(cap_getCoordinate(cap) < cap_getCoordinate(adjacentCap));
                stPinchThread *thread = stPinchThreadSet_addThread(threadSet, cap_getName(cap), cap_getCoordinate(cap),
                        cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) + 1);
                stPinchThread_split(thread, cap_getCoordinate(cap));
                stPinchThread_split(thread, cap_getCoordinate(adjacentCap) - 1);
                stPinchSegment *_5PrimeSegment = stPinchThread_getFirst(thread);
                stPinchSegment *_3PrimeSegment = stPinchThread_getLast(thread);
                assert(stPinchSegment_getStart(_5PrimeSegment) == cap_getCoordinate(cap));
                assert(stPinchSegment_getStart(_3PrimeSegment) == cap_getCoordinate(adjacentCap));
                stCaf_constructEmptyPinchGraphP(end, _5PrimeSegment, 1, endsToBlocks);
                stCaf_constructEmptyPinchGraphP(cap_getEnd(adjacentCap), _3PrimeSegment, 0, endsToBlocks);
            }
        }
        end_destructInstanceIterator(capIt);
    }
    flower_destructEndIterator(endIt);
    stHash_destruct(endsToBlocks);
    return threadSet;
}

///////////////////////////////////////////////////////////////////////////
// Add alignments to pinch graph
///////////////////////////////////////////////////////////////////////////

static void stCaf_ensureEndsAreDistinct(stPinchThreadSet *threadSet) {
    /*
     * Ensures the blocks at the ends of threads are distinct.
     */
    stPinchThread *thread;
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
        stPinchThread_split(thread, stPinchThread_getStart(thread));
        stPinchThread_split(thread, stPinchThread_getStart(thread) + stPinchThread_getLength(thread) - 1);
    }
}

void stCaf_addAlignmentsToPinchGraph(stPinchThreadSet *threadSet, stPinch *(*pinchGenerator)()) {
    stPinch *pinch;
    while ((pinch = pinchGenerator()) != NULL) {
        stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
        stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
        assert(thread1 != NULL && thread2 != NULL);
        assert(stPinchThread_getStart(thread1) < pinch->start1);
        assert(stPinchThread_getStart(thread1) + stPinchThread_getLength(thread1) > pinch->start1 + pinch->length);
        assert(stPinchThread_getStart(thread2) < pinch->start2);
        assert(stPinchThread_getStart(thread2) + stPinchThread_getLength(thread2) > pinch->start2 + pinch->length);
        stPinchThread_pinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand);
    }
    stPinchThreadSet_joinTrivialBoundaries(threadSet);
    stCaf_ensureEndsAreDistinct(threadSet); //Ensure that end blocks are not joined at there ends.
}

///////////////////////////////////////////////////////////////////////////
// Construct dead end component
///////////////////////////////////////////////////////////////////////////

static void attachThreadToDeadEndComponentP(stPinchBlock *pinchBlock, bool orientation, stList *deadEndAdjacencyComponent,
        stSortedSet *adjacencyComponents, stHash *pinchEndsToAdjacencyComponents) {
    assert(pinchBlock != NULL);
    assert(stPinchBlock_getLength(pinchBlock) == 1);
    stPinchEnd staticPinchEnd = stPinchEnd_constructStatic(pinchBlock, orientation);
    stList *component = stHash_remove(pinchEndsToAdjacencyComponents, &staticPinchEnd);
    assert(component != NULL);
    assert(stList_length(component) == 1);
    stPinchEnd *pinchEnd = stList_get(component, 0);
    stList_destruct(component);
    stList_append(deadEndAdjacencyComponent, pinchEnd);
    stHash_insert(pinchEndsToAdjacencyComponents, pinchEnd, deadEndAdjacencyComponent);
}

static void attachThreadToDeadEndComponent(stPinchThread *thread, stList *deadEndAdjacencyComponent, stSortedSet *adjacencyComponents,
        stHash *pinchEndsToAdjacencyComponents) {
    attachThreadToDeadEndComponentP(stPinchSegment_getBlock(stPinchThread_getFirst(thread)), 1, deadEndAdjacencyComponent,
            adjacencyComponents, pinchEndsToAdjacencyComponents);
    attachThreadToDeadEndComponentP(stPinchSegment_getBlock(stPinchThread_getLast(thread)), 0, deadEndAdjacencyComponent,
            adjacencyComponents, pinchEndsToAdjacencyComponents);
}

static void stCaf_constructDeadEndComponentP(End *end, stPinchThreadSet *threadSet, stList *deadEndAdjacencyComponent,
        stSortedSet *adjacencyComponents, stHash *pinchEndsToAdjacencyComponents) {
    End_InstanceIterator *capIt = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(capIt)) != NULL) {
        cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
        if (!cap_getSide(cap)) {
            end_destructInstanceIterator(capIt);
            stPinchThread *thread = stPinchThreadSet_getThread(threadSet, cap_getName(cap));
            assert(thread != NULL);
            attachThreadToDeadEndComponent(thread, deadEndAdjacencyComponent, adjacencyComponents, pinchEndsToAdjacencyComponents);
            return;
        }
    }
    end_destructInstanceIterator(capIt);
}

stList *stCaf_constructDeadEndComponent(Flower *flower, stPinchThreadSet *threadSet, stSortedSet *adjacencyComponents,
        stHash *pinchEndsToAdjacencyComponents) {
    stList *deadEndAdjacencyComponent = stList_construct();
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        assert(!end_isBlockEnd(end));
        if (end_isAttached(end)) {
            stCaf_constructDeadEndComponentP(end, threadSet, deadEndAdjacencyComponent, adjacencyComponents, pinchEndsToAdjacencyComponents);
        }
    }
    flower_destructEndIterator(endIt);
    return deadEndAdjacencyComponent;
}

///////////////////////////////////////////////////////////////////////////
// Attach unatttached thread components
///////////////////////////////////////////////////////////////////////////

static bool threadIsAttachedToDeadEndComponent(stPinchThread *thread, stList *deadEndComponent, stHash *pinchEndsToAdjacencyComponents) {
    stPinchBlock *pinchBlock = stPinchSegment_getBlock(stPinchThread_getFirst(thread));
    assert(pinchBlock != NULL);
    stPinchEnd staticPinchEnd = stPinchEnd_constructStatic(pinchBlock, 1);
    stList *adjacencyComponent = stHash_search(pinchEndsToAdjacencyComponents, &staticPinchEnd);
    assert(adjacencyComponent != NULL);
    return adjacencyComponent == deadEndComponent;
}

static bool threadComponentIsAttachedToDeadEndComponent(stList *threadComponent, stList *deadEndComponent,
        stHash *pinchEndsToAdjacencyComponents) {
    for (int32_t i = 0; i < stList_length(threadComponent); i++) {
        if (threadIsAttachedToDeadEndComponent(stList_get(threadComponent, i), deadEndComponent, pinchEndsToAdjacencyComponents)) {
            return 1;
        }
    }
    return 0;
}

static void attachThreadComponentToDeadEndComponent(stList *threadComponent, stList *deadEndComponent, stSortedSet *adjacencyComponents,
        stHash *pinchEndsToAdjacencyComponents) {
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
    attachThreadToDeadEndComponent(longestPinchThread, deadEndComponent, adjacencyComponents, pinchEndsToAdjacencyComponents);
}

void stCaf_attachUnattachedThreadComponents(Flower *flower, stPinchThreadSet *threadSet, stList *deadEndComponent,
        stSortedSet *adjacencyComponents, stHash *pinchEndsToAdjacencyComponents, bool markEndsAttached) {
    stSortedSet *threadComponents = stPinchThreadSet_getThreadComponents(threadSet);
    assert(stSortedSet_size(threadComponents) > 0);
    if (stSortedSet_size(threadComponents) > 1) {
        stSortedSetIterator *threadIt = stSortedSet_getIterator(threadComponents);
        stList *threadComponent;
        while ((threadComponent = stSortedSet_getNext(threadIt)) != NULL) {
            if (!threadComponentIsAttachedToDeadEndComponent(threadComponent, deadEndComponent, pinchEndsToAdjacencyComponents)) {
                attachThreadComponentToDeadEndComponent(threadComponent, deadEndComponent, adjacencyComponents,
                        pinchEndsToAdjacencyComponents);
            }
        }
        stSortedSet_destructIterator(threadIt);
    }
    stSortedSet_destruct(threadComponents);
    //#ifdef BEN_DEBUG
    threadComponents = stPinchThreadSet_getThreadComponents(threadSet);
    assert(stSortedSet_size(threadComponents) == 1);
    stSortedSet_destruct(threadComponents);
    //#endif
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

static void *mergeNodeObjects(void *a, void *b) {
    stList *adjacencyComponents1 = a;
    stList *adjacencyComponents2 = b;
    assert(adjacencyComponents1 != adjacencyComponents2);
    stList_appendAll(adjacencyComponents1, adjacencyComponents2);
    stList_setDestructor(adjacencyComponents2, NULL);
    stList_destruct(adjacencyComponents2);
    return adjacencyComponents1;
}

static void *makeNodeObject(stList *adjacencyComponent) {
    stList *adjacencyComponents = stList_construct();
    stList_append(adjacencyComponents, adjacencyComponent);
    return adjacencyComponents;
}

stCactusGraph *stCaf_constructCactusGraph(stSortedSet *adjacencyComponents, stList *deadEndComponent,
        stHash *pinchEndsToAdjacencyComponents, stCactusNode **startCactusNode) {
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
    stHash *pinchEndsToNonStaticPinchEnds = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, NULL, NULL);
    stHashIterator *pinchEndIt = stHash_getIterator(pinchEndsToAdjacencyComponents);
    stPinchEnd *pinchEnd;
    while ((pinchEnd = stHash_getNext(pinchEndIt)) != NULL) {
        stHash_insert(pinchEndsToNonStaticPinchEnds, pinchEnd, pinchEnd);
    }
    stHash_destructIterator(pinchEndIt);

    //Make the edges
    pinchEndIt = stHash_getIterator(pinchEndsToAdjacencyComponents);
    while ((pinchEnd = stHash_getNext(pinchEndIt)) != NULL) {
        if (stPinchEnd_getOrientation(pinchEnd)) {  //Assure we make the edge only once
            stPinchEnd pinchEnd2Static = stPinchEnd_constructStatic(stPinchEnd_getBlock(pinchEnd), !stPinchEnd_getOrientation(pinchEnd));
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
    stCactusGraph_collapseToCactus(cactusGraph, mergeNodeObjects, *startCactusNode);
    stCactusGraph_collapseBridges(cactusGraph, *startCactusNode, mergeNodeObjects);

    return cactusGraph;
}

///////////////////////////////////////////////////////////////////////////
// Convert the complete cactus graph/pinch graph into filled out set of flowers
///////////////////////////////////////////////////////////////////////////

//Functions used to build hash between pinchEnds and flower ends.

static void getPinchBlockEndsToEndsHashPP(stPinchBlock *pinchBlock, bool orientation, End *end, stHash *pinchEndsToEnds) {
    stPinchEnd pinchEnd = stPinchEnd_constructStatic(pinchBlock, orientation);
    if (stHash_search(pinchEndsToEnds, &pinchEnd) == NULL) {
        stHash_insert(pinchEndsToEnds, stPinchEnd_construct(pinchBlock, orientation), end);
    }
    else {
        assert(stHash_search(pinchEndsToEnds, &pinchEnd) == end);
    }
}

static void getPinchBlockEndsToEndsHashP(stPinchSegment *pinchSegment, bool endOrientation, Cap *cap, stHash *pinchEndsToEnds) {
    stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
    assert(pinchBlock != NULL);
    assert(cap != NULL);
    End *end = end_getPositiveOrientation(cap_getEnd(cap));
    assert(end != NULL);
    assert(!end_isBlockEnd(end));
    getPinchBlockEndsToEndsHashPP(pinchBlock, endOrientation, end, pinchEndsToEnds);
    getPinchBlockEndsToEndsHashPP(pinchBlock, !endOrientation, end_getReverse(end), pinchEndsToEnds);
}

static stHash *getPinchEndsToEndsHash(stPinchThreadSet *threadSet, Flower *parentFlower) {
    stHash *pinchEndsToEnds = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void(*)(void *)) stPinchEnd_destruct, NULL);
    stPinchThreadSetIt pinchThreadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *pinchThread;
    while ((pinchThread = stPinchThreadSetIt_getNext(&pinchThreadIt))) {
        Cap *cap = flower_getCap(parentFlower, stPinchThread_getName(pinchThread));
        assert(cap != NULL);
        stPinchSegment *pinchSegment = stPinchThread_getFirst(pinchThread);
        getPinchBlockEndsToEndsHashP(pinchSegment, !stPinchSegment_getBlockOrientation(pinchSegment), cap, pinchEndsToEnds);
        pinchSegment = stPinchThread_getLast(pinchThread);
        getPinchBlockEndsToEndsHashP(pinchSegment, stPinchSegment_getBlockOrientation(pinchSegment), cap_getAdjacency(cap), pinchEndsToEnds);
    }
    return pinchEndsToEnds;
}

//Functions for going from cactus/pinch ends to flower ends and updating flower structure as necessary

static End *convertPinchBlockEndToEnd(stPinchEnd *pinchEnd, stHash *pinchEndsToEnds, Flower *flower) {
    End *end = stHash_search(pinchEndsToEnds, pinchEnd);
    if (end == NULL) { //Happens if pinch end represents end of a block in flower that has not yet been defined.
        return NULL;
    }
    End *end2 = flower_getEnd(flower, end_getName(end));
    if (end2 == NULL) { //Happens if is free stub end not yet defined in given flower - but must be present somewhere in the hierarchy.
        assert(end_isFree(end));
        assert(end_isStubEnd(end));
        //Copy the end down the hierarchy
        Group *parentGroup = flower_getParentGroup(flower);
        assert(parentGroup != NULL); //Can not be at the top of the hierarchy, else would be defined.
        end2 = convertPinchBlockEndToEnd(pinchEnd, pinchEndsToEnds, group_getFlower(parentGroup));
        assert(end2 != NULL);
        assert(end_getGroup(end2) == NULL); //This happens because group of the end is only defined at the given flower.
        end_setGroup(end2, parentGroup);
        end2 = end_copyConstruct(end2, flower);
        assert(end2 != NULL);
        assert(end_getFlower(end2) == flower);
    }
    return end2;
}

static End *convertCactusEdgeEndToEnd(stCactusEdgeEnd *cactusEdgeEnd, stHash *pinchEndsToEnds, Flower *flower) {
    return convertPinchBlockEndToEnd(stCactusEdgeEnd_getObject(cactusEdgeEnd), pinchEndsToEnds, flower);
}

//Functions to create blocks

static void makeBlockP(stPinchEnd *pinchEnd, End *end, stHash *pinchEndsToEnds) {
    assert(stHash_search(pinchEndsToEnds, pinchEnd) == NULL);
    stHash_insert(pinchEndsToEnds, stPinchEnd_construct(stPinchEnd_getBlock(pinchEnd), stPinchEnd_getOrientation(pinchEnd)), end);
}

static void makeBlock(stCactusEdgeEnd *cactusEdgeEnd, Flower *parentFlower, Flower *flower, stHash *pinchEndsToEnds) {
    stPinchEnd *pinchEnd = stCactusEdgeEnd_getObject(cactusEdgeEnd);
    assert(pinchEnd != NULL);
    stPinchBlock *pinchBlock = stPinchEnd_getBlock(pinchEnd);
    Block *block = block_construct(stPinchBlock_getLength(pinchBlock), flower);
    stPinchSegment *pinchSegment;
    stPinchBlockIt pinchSegmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    while ((pinchSegment = stPinchBlockIt_getNext(&pinchSegmentIt))) {
        Cap *cap = flower_getCap(parentFlower, stPinchSegment_getName(pinchSegment)); //The following three lines isolates the sequence associated with a segment.
        assert(cap != NULL);
        Sequence *sequence = cap_getSequence(cap);
        assert(sequence != NULL);
        segment_construct2(
                stPinchEnd_getOrientation(pinchEnd) ^ stPinchSegment_getBlockOrientation(pinchSegment) ? block_getReverse(block) : block,
                stPinchSegment_getStart(pinchSegment), 1, sequence);
    }
    makeBlockP(pinchEnd, block_get5End(block), pinchEndsToEnds);
    stPinchEnd *otherPinchBlockEnd = stCactusEdgeEnd_getObject(stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd));
    makeBlockP(otherPinchBlockEnd, block_get3End(block), pinchEndsToEnds);
}

//Functions to generate the chains of a flower

static void makeChain(stCactusEdgeEnd *cactusEdgeEnd, Flower *flower, stHash *pinchEndsToEnds, Flower *parentFlower, stList *stack) {
    cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd);
    if (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd)) { //We have a non-trivial chain
        Chain *chain = chain_construct(flower);
        do {
            stCactusEdgeEnd *linkedCactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd);
            if (convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchEndsToEnds, flower) == NULL) { //Make subsequent block
                makeBlock(linkedCactusEdgeEnd, parentFlower, flower, pinchEndsToEnds);
            }
            assert(stCactusEdgeEnd_getNode(cactusEdgeEnd) == stCactusEdgeEnd_getNode(linkedCactusEdgeEnd));
            Group *group = group_construct2(flower);
            End *end1 = convertCactusEdgeEndToEnd(cactusEdgeEnd, pinchEndsToEnds, flower);
            End *end2 = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchEndsToEnds, flower);
            assert(end1 != NULL);
            assert(end2 != NULL);
            assert(end_getOrientation(end1));
            assert(end_getOrientation(end2));
            assert(!end_getSide(end1));
            assert(end_getSide(end2));
            end_setGroup(end1, group);
            end_setGroup(end2, group);
            link_construct(end1, end2, group, chain);
            stList_append(stack, stCactusEdgeEnd_getNode(cactusEdgeEnd));
            stList_append(stack, group_makeNestedFlower(group));
            cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(linkedCactusEdgeEnd);
        } while (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd));
    }
}

static void makeChains(stCactusNode *cactusNode, Flower *flower, stHash *pinchEndsToEnds, Flower *parentFlower, stList *stack) {
    stCactusNodeEdgeEndIt cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
    stCactusEdgeEnd *cactusEdgeEnd;
    while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt))) {
        if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) { //We have some sort of chain
            End *end = convertCactusEdgeEndToEnd(cactusEdgeEnd, pinchEndsToEnds, flower);
            stCactusEdgeEnd *linkedCactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd), *startCactusEdgeEnd = NULL;
            if (end != NULL) {
//#ifdef BEN_DEBUG
                End *end2;
                if ((end2 = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchEndsToEnds, flower)) != NULL) {
                    assert(end_getSide(end) != end_getSide(end2));
                }
//#endif
                startCactusEdgeEnd = end_getSide(end) ? cactusEdgeEnd : linkedCactusEdgeEnd;
            } else {
                end = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchEndsToEnds, flower);
                if (end != NULL) {
                    startCactusEdgeEnd = end_getSide(end) ? linkedCactusEdgeEnd : cactusEdgeEnd;
                } else {
                    makeBlock(linkedCactusEdgeEnd, parentFlower, flower, pinchEndsToEnds);
                    startCactusEdgeEnd = linkedCactusEdgeEnd;
                }
            }
            assert(startCactusEdgeEnd != NULL);
            makeChain(startCactusEdgeEnd, flower, pinchEndsToEnds, parentFlower, stack);
        }
    }
}

//Functions to make tangles.

static void makeTangles(stCactusNode *cactusNode, Flower *flower, stHash *pinchEndsToEnds, stList *deadEndComponent) {
    stList *adjacencyComponents = stCactusNode_getObject(cactusNode);
    for (int32_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stList *adjacencyComponent = stList_get(adjacencyComponents, i);
        if (adjacencyComponent != deadEndComponent) {
            Group *group = group_construct2(flower);
            for (int32_t j = 0; j < stList_length(adjacencyComponent); j++) {
                stPinchEnd *pinchEnd = stList_get(adjacencyComponent, j);
                assert(pinchEnd != NULL);
                End *end = convertPinchBlockEndToEnd(pinchEnd, pinchEndsToEnds, flower);
                assert(end != NULL);
                end_setGroup(end, group);
            }
        }
    }
}

//Main function

void stCaf_convertCactusGraphToFlowers(stPinchThreadSet *threadSet, stCactusNode *startCactusNode, Flower *parentFlower,
        stList *deadEndComponent) {
    stList *stack = stList_construct();
    stList_append(stack, startCactusNode);
    stList_append(stack, parentFlower);
    stHash *pinchEndsToEnds = getPinchEndsToEndsHash(threadSet, parentFlower);
    while (stList_length(stack) > 0) {
        Flower *flower = stList_pop(stack);
        stCactusNode *cactusNode = stList_pop(stack);
        makeChains(cactusNode, flower, pinchEndsToEnds, parentFlower, stack);
        makeTangles(cactusNode, flower, pinchEndsToEnds, deadEndComponent);
    }
    stHash_destruct(pinchEndsToEnds);
    stList_destruct(stack);
    stCaf_addAdjacencies(parentFlower);
}

///////////////////////////////////////////////////////////////////////////
// Functions for actually filling out cactus
///////////////////////////////////////////////////////////////////////////

void stCore(Flower *flower, stPinch *(*pinchGenerator)()) {
    //Create empty pinch graph from flower
    stPinchThreadSet *threadSet = stCaf_constructEmptyPinchGraph(flower);

    //Add alignments to pinch graph
    stCaf_addAlignmentsToPinchGraph(threadSet, pinchGenerator);

    //Get adjacency components
    stHash *pinchEndsToAdjacencyComponents;
    stSortedSet *adjacencyComponents = stList_convertToSortedSet(
            stPinchThreadSet_getAdjacencyComponents2(threadSet, &pinchEndsToAdjacencyComponents));

    //Merge together dead end component
    stList *deadEndComponent = stCaf_constructDeadEndComponent(flower, threadSet, adjacencyComponents, pinchEndsToAdjacencyComponents);

    //Join unattached components of graph by dead ends to dead end component, and make other ends 'attached' if necessary
    stCaf_attachUnattachedThreadComponents(flower, threadSet, deadEndComponent, adjacencyComponents, pinchEndsToAdjacencyComponents, 1);

    //Create cactus
    stCactusNode *startCactusNode;
    stCactusGraph *cactusGraph = stCaf_constructCactusGraph(adjacencyComponents, deadEndComponent, pinchEndsToAdjacencyComponents,
            &startCactusNode);

    //Convert cactus graph/pinch graph to API
    stCaf_convertCactusGraphToFlowers(threadSet, startCactusNode, flower, deadEndComponent);

    //Cleanup
    stCactusGraph_destruct(cactusGraph);
    stHash_destruct(pinchEndsToAdjacencyComponents);
    stSortedSet_destruct(adjacencyComponents);
    stPinchThreadSet_destruct(threadSet);
}
