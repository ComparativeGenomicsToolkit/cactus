#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
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
        assert(stPinchThread_getLength(thread) >= 2);
        stPinchThread_split(thread, stPinchThread_getStart(thread) + stPinchThread_getLength(thread) - 2);
    }
}

void stCaf_addAlignmentsToPinchGraph(stPinchThreadSet *threadSet, stPinchIterator *pinchIterator) {
    stPinchIterator_reset(pinchIterator);
    stPinch *pinch;
    while ((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
        stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
        stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
        assert(thread1 != NULL && thread2 != NULL);
        assert(stPinchThread_getStart(thread1) < pinch->start1);
        assert(stPinchThread_getStart(thread1) + stPinchThread_getLength(thread1) > pinch->start1 + pinch->length);
        assert(stPinchThread_getStart(thread2) < pinch->start2);
        assert(stPinchThread_getStart(thread2) + stPinchThread_getLength(thread2) > pinch->start2 + pinch->length);
        stPinchThread_pinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand);
        stPinchIterator_destructAlignment(pinchIterator, pinch);
    }
    stPinchThreadSet_joinTrivialBoundaries(threadSet);
    stCaf_ensureEndsAreDistinct(threadSet); //Ensure that end blocks are not joined at there ends.
}

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

stList *stCaf_constructDeadEndComponent(Flower *flower, stPinchThreadSet *threadSet,
        stHash *pinchEndsToAdjacencyComponents) {
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

static bool threadIsAttachedToDeadEndComponent(stPinchThread *thread, stList *deadEndComponent, stHash *pinchEndsToAdjacencyComponents) {
    stPinchSegment *pinchSegment = stPinchThread_getFirst(thread);
    stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
    assert(pinchBlock != NULL);
    stPinchEnd staticPinchEnd = stPinchEnd_constructStatic(pinchBlock, stPinchSegment_getBlockOrientation(pinchSegment));
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

static void attachThreadToDeadEndComponent(stPinchThread *thread, stList *deadEndAdjacencyComponent,
        stHash *pinchEndsToAdjacencyComponents) {
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

void stCaf_attachUnattachedThreadComponents(Flower *flower, stPinchThreadSet *threadSet, stList *deadEndComponent,
       stHash *pinchEndsToAdjacencyComponents, bool markEndsAttached) {
    stSortedSet *threadComponents = stPinchThreadSet_getThreadComponents(threadSet);
    assert(stSortedSet_size(threadComponents) > 0);
    if (stSortedSet_size(threadComponents) > 1) {
        stSortedSetIterator *threadIt = stSortedSet_getIterator(threadComponents);
        stList *threadComponent;
        while ((threadComponent = stSortedSet_getNext(threadIt)) != NULL) {
            if (!threadComponentIsAttachedToDeadEndComponent(threadComponent, deadEndComponent, pinchEndsToAdjacencyComponents)) {
                attachThreadComponentToDeadEndComponent(threadComponent, deadEndComponent,
                        pinchEndsToAdjacencyComponents, markEndsAttached, flower);
            }
        }
        stSortedSet_destructIterator(threadIt);
    }
    stSortedSet_destruct(threadComponents);
#ifdef BEN_DEBUG
    threadComponents = stPinchThreadSet_getThreadComponents(threadSet);
    assert(stSortedSet_size(threadComponents) == 1);
    stSortedSet_destruct(threadComponents);
#endif
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
    stList *adjacencyComponents = stList_construct3(0, (void (*)(void *))stList_destruct);
    stList_append(adjacencyComponents, adjacencyComponent);
    return adjacencyComponents;
}

stCactusGraph *stCaf_constructCactusGraph(stList *deadEndComponent,
        stHash *pinchEndsToAdjacencyComponents, stCactusNode **startCactusNode) {
    stCactusGraph *cactusGraph = stCactusGraph_construct2((void(*)(void *)) stList_destruct, NULL);
    stHash *adjacencyComponentsToCactusNodes = stHash_construct();
	stHash *pinchEndsToNonStaticPinchEnds = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, NULL, NULL);
    
    //Make the nodes
    *startCactusNode = stCactusNode_construct(cactusGraph, makeNodeObject(deadEndComponent));
    stHash_insert(adjacencyComponentsToCactusNodes, deadEndComponent, *startCactusNode);
    stHashIterator *pinchEndIt = stHash_getIterator(pinchEndsToAdjacencyComponents);
    stPinchEnd *pinchEnd;
    while((pinchEnd = stHash_getNext(pinchEndIt)) != NULL) {
    	//Make a map of edge ends to themselves for memory lookup
    	stHash_insert(pinchEndsToNonStaticPinchEnds, pinchEnd, pinchEnd); 
        stList *adjacencyComponent = stHash_search(pinchEndsToAdjacencyComponents, pinchEnd);
        assert(adjacencyComponent != NULL);
        if (stHash_search(adjacencyComponentsToCactusNodes, adjacencyComponent) == NULL) {
            if (stList_length(adjacencyComponent) == 1) { //Going to be a bridge to nowhere, so we join it - this ensures
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
    stCactusGraph_collapseToCactus(cactusGraph, mergeNodeObjects, *startCactusNode);
    stCactusGraph_collapseBridges(cactusGraph, *startCactusNode, mergeNodeObjects);

    return cactusGraph;
}

///////////////////////////////////////////////////////////////////////////
// Function that draws together above functions to generate a cactus graph from a pinch graph.
///////////////////////////////////////////////////////////////////////////

stCactusGraph *stCaf_getCactusGraphForThreadSet(Flower *flower, stPinchThreadSet *threadSet, stCactusNode **startCactusNode, stList **deadEndComponent) {
    //Get adjacency components
    stHash *pinchEndsToAdjacencyComponents;
    stList *adjacencyComponents =
            stPinchThreadSet_getAdjacencyComponents2(threadSet, &pinchEndsToAdjacencyComponents);
    stList_setDestructor(adjacencyComponents, NULL);
    stList_destruct(adjacencyComponents);

    //Merge together dead end component
    *deadEndComponent = stCaf_constructDeadEndComponent(flower, threadSet, pinchEndsToAdjacencyComponents);

    //Join unattached components of graph by dead ends to dead end component, and make other ends 'attached' if necessary
    stCaf_attachUnattachedThreadComponents(flower, threadSet, *deadEndComponent, pinchEndsToAdjacencyComponents, 1);

    //Create cactus
    stCactusGraph *cactusGraph = stCaf_constructCactusGraph(*deadEndComponent, pinchEndsToAdjacencyComponents,
            startCactusNode);

    //Cleanup (the memory is owned by the cactus graph, so this does not break anything)
    stHash_destruct(pinchEndsToAdjacencyComponents);

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
    } else {
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
    assert(end_getOrientation(end));
    assert(!end_getOrientation(end_getReverse(end)));
    getPinchBlockEndsToEndsHashPP(pinchBlock, endOrientation, end_getReverse(end), pinchEndsToEnds);
    getPinchBlockEndsToEndsHashPP(pinchBlock, !endOrientation, end, pinchEndsToEnds);
}

static stHash *getPinchEndsToEndsHash(stPinchThreadSet *threadSet, Flower *parentFlower) {
    stHash *pinchEndsToEnds = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void(*)(void *)) stPinchEnd_destruct, NULL);
    stPinchThreadSetIt pinchThreadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *pinchThread;
    while ((pinchThread = stPinchThreadSetIt_getNext(&pinchThreadIt))) {
        Cap *cap = flower_getCap(parentFlower, stPinchThread_getName(pinchThread));
        assert(cap != NULL);
        stPinchSegment *pinchSegment = stPinchThread_getFirst(pinchThread);
        getPinchBlockEndsToEndsHashP(pinchSegment, stPinchSegment_getBlockOrientation(pinchSegment), cap, pinchEndsToEnds);
        pinchSegment = stPinchThread_getLast(pinchThread);
        getPinchBlockEndsToEndsHashP(pinchSegment, !stPinchSegment_getBlockOrientation(pinchSegment), cap_getAdjacency(cap),
                pinchEndsToEnds);
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
        end2 = end_copyConstruct(end_getPositiveOrientation(end2), flower);
        assert(end2 != NULL);
        assert(end_getFlower(end2) == flower);
    }
    assert(end_getOrientation(end2));
    return end_getOrientation(end) ? end2 : end_getReverse(end2);
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
        Cap *parentCap = flower_getCap(parentFlower, stPinchSegment_getName(pinchSegment)); //The following three lines isolates the sequence associated with a segment.
        assert(parentCap != NULL);
        Sequence *parentSequence = cap_getSequence(parentCap);
        assert(parentSequence != NULL);
        Sequence *sequence = flower_getSequence(flower, sequence_getName(parentSequence));
        if (sequence == NULL) {
            sequence = sequence_construct(cactusDisk_getMetaSequence(flower_getCactusDisk(flower), sequence_getName(sequence)), flower);
        }
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
            assert(end_isBlockEnd(end1) || end_isAttached(end1));
            assert(end_isBlockEnd(end2) || end_isAttached(end2));
            end_setGroup(end1, group);
            end_setGroup(end2, group);
            link_construct(end1, end2, group, chain);
            //Make a nested group
            Flower *nestedFlower = group_makeEmptyNestedFlower(group);
            end_copyConstruct(end1, nestedFlower);
            end_copyConstruct(end2, nestedFlower);
            assert(flower_getGroupNumber(nestedFlower) == 0);
            //Fill out stack
            stList_append(stack, stCactusEdgeEnd_getNode(cactusEdgeEnd));
            stList_append(stack, nestedFlower);
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
            assert(linkedCactusEdgeEnd != NULL);
            if (end != NULL) {
#ifdef BEN_DEBUG
                End *end2;
                if ((end2 = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchEndsToEnds, flower)) != NULL) {
                    assert(end_getSide(end) != end_getSide(end2));
                }
#endif
                startCactusEdgeEnd = end_getSide(end) ? cactusEdgeEnd : linkedCactusEdgeEnd;
            } else {
                end = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchEndsToEnds, flower);
                if (end != NULL) {
                    if (end_getSide(end)) {
                        startCactusEdgeEnd = linkedCactusEdgeEnd;
                    } else {
                        makeBlock(cactusEdgeEnd, parentFlower, flower, pinchEndsToEnds);
                        startCactusEdgeEnd = cactusEdgeEnd;
                    }
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
            if (stList_length(adjacencyComponent) == 1) { //Deal with components for dead ends of free stubs
                End *end = convertPinchBlockEndToEnd(stList_get(adjacencyComponent, 0), pinchEndsToEnds, flower);
                assert(end != NULL);
                if (!end_getOrientation(end)) {
                    continue;
                }
            }
            Group *group = group_construct2(flower);
            for (int32_t j = 0; j < stList_length(adjacencyComponent); j++) {
                End *end = convertPinchBlockEndToEnd(stList_get(adjacencyComponent, j), pinchEndsToEnds, flower);
                assert(end != NULL);
                assert(end_getOrientation(end));
                assert(end_getGroup(end) == NULL);
                end_setGroup(end, group);
            }
        }
    }
}

//Sets the 'built-blocks flag' for all the flowers in the subtree, including the given flower.

static void setBlocksBuilt(Flower *flower) {
    //#ifdef BEN_DEBUG
    assert(!flower_builtBlocks(flower));
    //#endif
    flower_setBuiltBlocks(flower, 1);
    Flower_GroupIterator *iterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(iterator)) != NULL) {
        if (!group_isLeaf(group)) {
            setBlocksBuilt(group_getNestedFlower(group));
        }
    }
    flower_destructGroupIterator(iterator);
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
        assert(flower_getAttachedStubEndNumber(flower) > 0);
        stCactusNode *cactusNode = stList_pop(stack);
        makeChains(cactusNode, flower, pinchEndsToEnds, parentFlower, stack);
        makeTangles(cactusNode, flower, pinchEndsToEnds, deadEndComponent);
    }
    stHash_destruct(pinchEndsToEnds);
    stList_destruct(stack);
    stCaf_addAdjacencies(parentFlower);
    setBlocksBuilt(parentFlower);
}

///////////////////////////////////////////////////////////////////////////
// Functions for actually filling out cactus
///////////////////////////////////////////////////////////////////////////

static void initialiseFlowerForFillingOut(Flower *flower) {
#ifdef BEN_DEBUG
    assert(!flower_builtBlocks(flower)); //We can't do this if we've already built blocks for the flower!.
    flower_check(flower);
    assert(flower_isTerminal(flower));
    assert(flower_getGroupNumber(flower) == 1);
    assert(group_isLeaf(flower_getFirstGroup(flower))); //this should be true by the previous assert
    //Destruct any chain
    assert(flower_getChainNumber(flower) <= 1);
#endif
    if (flower_getChainNumber(flower) == 1) {
        Chain *chain = flower_getFirstChain(flower);
        chain_destruct(chain);
    }
    group_destruct(flower_getFirstGroup(flower));
}

/*void stThreadSet_trimAlignments(stThreadSet *threadSet, int32_t minimumBlockLength) {

}

void stThreadSet_filterAlignments(stThreadSet *threadSet, bool (*blockFilterFn)(stPinchBlock *)) {

}

stList *stCaf_getBlocksInChainsLessThanGivenLength(stCactusGraph *cactusGraph, int32_t minimumChainLength) {

}*/

void stCaf_core(Flower *flower, stPinchIterator *pinchIterator) {
    //Setup the empty flower that will be filled out
    initialiseFlowerForFillingOut(flower);

    //Create empty pinch graph from flower
    stPinchThreadSet *threadSet = stCaf_constructEmptyPinchGraph(flower);

    //Add alignments to pinch graph
    stCaf_addAlignmentsToPinchGraph(threadSet, pinchIterator);

    stCactusNode *startCactusNode;
    stList *deadEndComponent;
    stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent);

    //Convert cactus graph/pinch graph to API
    stCaf_convertCactusGraphToFlowers(threadSet, startCactusNode, flower, deadEndComponent);

    //Cleanup
    stCactusGraph_destruct(cactusGraph);
    stPinchThreadSet_destruct(threadSet);

#ifdef BEN_DEBUG
    flower_checkRecursive(flower);
#endif
}
