#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

///////////////////////////////////////////////////////////////////////////
// Core functions for melting
///////////////////////////////////////////////////////////////////////////

static bool isThreadEnd(stPinchBlock *pinchBlock) {
    stPinchSegment *pinchSegment = stPinchBlock_getFirst(pinchBlock);
    bool threadEnd = pinchSegment != NULL && (stPinchSegment_get3Prime(pinchSegment) == NULL || stPinchSegment_get5Prime(pinchSegment)
            == NULL);
    if (threadEnd) {
        assert(stPinchBlock_getLength(pinchBlock) == 1);
    }
    return threadEnd;
}

static void processChain(stCactusEdgeEnd *cactusEdgeEnd, void(*edgeEndFn)(stPinchBlock *, void *), void *extraArg, bool recursive) {
    while (1) {
        stPinchEnd *pinchEnd = stCactusEdgeEnd_getObject(cactusEdgeEnd);
        assert(pinchEnd != NULL);
        stPinchBlock *pinchBlock = stPinchEnd_getBlock(pinchEnd);
        assert(pinchBlock != NULL);
        edgeEndFn(pinchBlock, extraArg);
        assert(stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd)) == cactusEdgeEnd);
        cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd);
        if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd)) {
            break;
        }
        if (recursive) {
            stCactusNode *node = stCactusEdgeEnd_getNode(cactusEdgeEnd);
            stCactusNodeEdgeEndIt it = stCactusNode_getEdgeEndIt(node);
            stCactusEdgeEnd *cactusEdgeEnd2;
            while ((cactusEdgeEnd2 = stCactusNodeEdgeEndIt_getNext(&it)) != NULL) {
                if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd2) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd2)) {
                    processChain(cactusEdgeEnd2, edgeEndFn, extraArg, 1);
                }
            }
        }
        assert(stCactusEdgeEnd_getLink(stCactusEdgeEnd_getLink(cactusEdgeEnd)) == cactusEdgeEnd);
        cactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd);
    }
}

static void addBlock(stPinchBlock *block, void *extraArg) {
    if (!isThreadEnd(block)) {
        stList_append(extraArg, block);
    }
}

static void addChainBlocksToBlocksToDelete(stCactusEdgeEnd *cactusEdgeEnd, stList *blocksToDelete) {
    processChain(cactusEdgeEnd, addBlock, blocksToDelete, 0);
}

static void addBlockRecursive(stPinchBlock *block, void *extraArg) {
    if (!isThreadEnd(block)) {
        stSet_insert(extraArg, block);
    }
}

static void addChainBlocksToBlocksToDeleteRecursive(stCactusEdgeEnd *cactusEdgeEnd, stSet *blocksToDelete) {
    processChain(cactusEdgeEnd, addBlockRecursive, blocksToDelete, 1);
}

static void addLength(stPinchBlock *block, void *extraArg) {
    *((int64_t *) extraArg) += stPinchBlock_getLength(block);
}

static int64_t getChainLength(stCactusEdgeEnd *cactusEdgeEnd) {
    int64_t length = 0;
    processChain(cactusEdgeEnd, addLength, &length, 0);
    return length;
}

static stList *stCaf_getBlocksInChainsLessThanGivenLength(stCactusGraph *cactusGraph, int64_t minimumChainLength) {
    stList *blocksToDelete = stList_construct3(0, (void(*)(void *)) stPinchBlock_destruct);
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(cactusGraph);
    stCactusNode *cactusNode;
    while ((cactusNode = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stCactusNodeEdgeEndIt cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
        stCactusEdgeEnd *cactusEdgeEnd;
        while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt)) != NULL) {
            if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) {
                if (getChainLength(cactusEdgeEnd) < minimumChainLength) {
                    addChainBlocksToBlocksToDelete(cactusEdgeEnd, blocksToDelete);
                }
            }
        }
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
    return blocksToDelete;
}

bool isBubbleChain(stCactusEdgeEnd *cactusEdgeEnd, stHash *pinchEndsToAdjacencyComponents, bool blockFilterfn(stPinchBlock *)) {
    stCactusEdgeEnd *cactusEdgeEnd2 = stCactusEdgeEnd_getLink(cactusEdgeEnd);
    assert(cactusEdgeEnd != cactusEdgeEnd2);
    assert(stCactusEdgeEnd_getNode(cactusEdgeEnd) == stCactusEdgeEnd_getNode(cactusEdgeEnd2));

    if (!blockFilterfn(stPinchEnd_getBlock(stCactusEdgeEnd_getObject(cactusEdgeEnd))) && !blockFilterfn(
            stPinchEnd_getBlock(stCactusEdgeEnd_getObject(cactusEdgeEnd2)))) {
        return 0;
    }

    stPinchEnd *pinchEnd1 = stCactusEdgeEnd_getObject(cactusEdgeEnd);
    stList *adjacencyComponent = stHash_search(pinchEndsToAdjacencyComponents, pinchEnd1);
    assert(adjacencyComponent != NULL);
    stPinchEnd *pinchEnd2 = stCactusEdgeEnd_getObject(cactusEdgeEnd2);
    stList *adjacencyComponent2 = stHash_search(pinchEndsToAdjacencyComponents, pinchEnd2);
    assert(adjacencyComponent2 != NULL);
    if(adjacencyComponent == adjacencyComponent2) {
        return 1;
    }

    //How many ends are we connected to.
    return stPinchEnd_getNumberOfConnectedPinchEnds(pinchEnd1) == 1 || stPinchEnd_getNumberOfConnectedPinchEnds(pinchEnd2) == 1;
}

static stSet *stCaf_getBlocksInBubbleChains(stPinchThreadSet *threadSet, stCactusGraph *cactusGraph, bool blockFilterfn(stPinchBlock *)) {
    stHash *pinchEndsToAdjacencyComponents;
    stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents2(threadSet, &pinchEndsToAdjacencyComponents);
    stSet *blocksToDelete = stSet_construct2((void(*)(void *)) stPinchBlock_destruct);
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(cactusGraph);
    stCactusNode *cactusNode;
    while ((cactusNode = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stCactusNodeEdgeEndIt cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
        stCactusEdgeEnd *cactusEdgeEnd;
        while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt)) != NULL) {
            if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) {
                if (isBubbleChain(cactusEdgeEnd, pinchEndsToAdjacencyComponents, blockFilterfn)) {
                    addChainBlocksToBlocksToDeleteRecursive(cactusEdgeEnd, blocksToDelete);
                }
            }
        }
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
    //cleanup
    stHash_destruct(pinchEndsToAdjacencyComponents);
    stList_destruct(adjacencyComponents);

    return blocksToDelete;
}

static void trimAlignments(stPinchThreadSet *threadSet, int64_t blockEndTrim) {
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *block = stPinchThreadSetBlockIt_getNext(&blockIt);
    while (block != NULL) {
        stPinchBlock *block2 = stPinchThreadSetBlockIt_getNext(&blockIt);
        if (!isThreadEnd(block)) {
            stPinchBlock_trim(block, blockEndTrim);
        }
        block = block2;
    }
}

static void filterAlignments(stPinchThreadSet *threadSet, bool(*blockFilterFn)(stPinchBlock *)) {
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *block = stPinchThreadSetBlockIt_getNext(&blockIt);
    while (block != NULL) {
        stPinchBlock *block2 = stPinchThreadSetBlockIt_getNext(&blockIt);
        if (!isThreadEnd(block) && blockFilterFn(block)) {
            stPinchBlock_destruct(block);
        }
        block = block2;
    }
}

static bool blockFilterByDegree(stPinchBlock *pinchBlock) {
    if (stPinchBlock_getDegree(pinchBlock) < 2) {
        return 1;
    }
    return 0;
}

void stCaf_melt(Flower *flower, stPinchThreadSet *threadSet, bool blockFilterfn(stPinchBlock *), int64_t blockEndTrim,
        int64_t minimumChainLength) {
    //First trim
    if (blockEndTrim > 0) {
        trimAlignments(threadSet, blockEndTrim);
    }
    //Then filter blocks
    filterAlignments(threadSet, blockFilterByDegree);

    if (blockFilterfn != NULL) {
        //Then filter blocks
        filterAlignments(threadSet, blockFilterByDegree);
        stCactusNode *startCactusNode;
        stList *deadEndComponent;
        stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, INT64_MAX,
                0.0);
        stSet *blocksToDelete = stCaf_getBlocksInBubbleChains(threadSet, cactusGraph, blockFilterfn);
        //Cleanup cactus
        stCactusGraph_destruct(cactusGraph);
        (void)blocksToDelete;
        stSet_destruct(blocksToDelete); //This will destroy the blocks
    }

    //Now apply the minimum chain length filter
    if (minimumChainLength > 1) {
        stCactusNode *startCactusNode;
        stList *deadEndComponent;
        stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, INT64_MAX,
                0.0);
        stList *blocksToDelete = stCaf_getBlocksInChainsLessThanGivenLength(cactusGraph, minimumChainLength);
        //Cleanup cactus
        stCactusGraph_destruct(cactusGraph);
        stList_destruct(blocksToDelete); //This will destroy the blocks
    }
    //Now heal up the trivial boundaries
    stCaf_joinTrivialBoundaries(threadSet);
}

///////////////////////////////////////////////////////////////////////////
// Functions for calculating required species/tree coverage
///////////////////////////////////////////////////////////////////////////

void stCaf_calculateRequiredFractionsOfSpecies(Flower *flower, float requiredIngroupFraction, float requiredOutgroupFraction,
        float requiredAllFraction, int64_t *requiredIngroupSpecies, int64_t *requiredOutgroupSpecies, int64_t *requiredAllSpecies) {
    if (requiredIngroupFraction <= 0.0 && requiredOutgroupFraction <= 0.0 && requiredAllFraction <= 0.0) {
        *requiredAllSpecies = 0;
        *requiredOutgroupSpecies = 0;
        *requiredAllSpecies = 0;
    }
    EventTree *eventTree = flower_getEventTree(flower);
    Event *event;
    int64_t outgroupEventNumber = 0;
    int64_t ingroupEventNumber = 0;
    EventTree_Iterator *eventIt = eventTree_getIterator(eventTree);
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        if (event_getChildNumber(event) == 0) {
            if (event_isOutgroup(event)) {
                outgroupEventNumber++;
            } else {
                ingroupEventNumber++;
            }
        }
    }
    eventTree_destructIterator(eventIt);
    *requiredOutgroupSpecies = 0.5 + outgroupEventNumber * requiredOutgroupFraction;
    *requiredIngroupSpecies = 0.5 + ingroupEventNumber * requiredIngroupFraction;
    if (*requiredIngroupSpecies == 0) {
        *requiredIngroupSpecies = 1;
    }
    *requiredAllSpecies = 0.5 + (ingroupEventNumber + outgroupEventNumber) * requiredAllFraction;
}

Event *getEvent(stPinchSegment *segment, Flower *flower) {
    Event *event = cap_getEvent(flower_getCap(flower, stPinchSegment_getName(segment)));
    assert(event != NULL);
    return event;
}

bool stCaf_containsRequiredSpecies(stPinchBlock *pinchBlock, Flower *flower, int64_t requiredIngroupSpecies,
        int64_t requiredOutgroupSpecies, int64_t requiredAllSpecies) {
    if (requiredIngroupSpecies <= 0 && requiredOutgroupSpecies <= 0 && requiredAllSpecies <= 0) {
        return 1;
    }
    int64_t outgroupSequences = 0;
    int64_t ingroupSequences = 0;
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        Event *event = getEvent(segment, flower);
        if (event_isOutgroup(event)) {
            outgroupSequences++;
        } else {
            ingroupSequences++;
        }
    }
    return ingroupSequences >= requiredIngroupSpecies && outgroupSequences >= requiredOutgroupSpecies && outgroupSequences
            + ingroupSequences >= requiredAllSpecies;
}

static bool stCaf_containsMultipleCopiesOfSpecies(stPinchBlock *pinchBlock, Flower *flower, bool(*acceptableEvent)(Event *)) {
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    stPinchSegment *segment;
    stHash *seen = stHash_construct();
    while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        Event *event = getEvent(segment, flower);
        if (acceptableEvent(event)) {
            if (stHash_search(seen, event) != NULL) {
                stHash_destruct(seen);
                return 1;
            }
            stHash_insert(seen, event, event);
        }
    }
    stHash_destruct(seen);
    return 0;
}

static bool returnTrue(Event *event) {
    return 1;
}

bool stCaf_containsMultipleCopiesOfAnySpecies(stPinchBlock *pinchBlock, Flower *flower) {
    return stCaf_containsMultipleCopiesOfSpecies(pinchBlock, flower, returnTrue);
}

static bool isIngroup(Event *event) {
    return !event_isOutgroup(event);
}

bool stCaf_containsMultipleCopiesOfIngroupSpecies(stPinchBlock *pinchBlock, Flower *flower) {
    return stCaf_containsMultipleCopiesOfSpecies(pinchBlock, flower, isIngroup);
}

bool stCaf_containsMultipleCopiesOfOutgroupSpecies(stPinchBlock *pinchBlock, Flower *flower) {
    return stCaf_containsMultipleCopiesOfSpecies(pinchBlock, flower, event_isOutgroup);
}

bool stCaf_treeCoverage(stPinchBlock *pinchBlock, Flower *flower) {
    EventTree *eventTree = flower_getEventTree(flower);
    Event *commonAncestorEvent = NULL;
    stPinchSegment *segment;
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    while ((segment = stPinchBlockIt_getNext(&segmentIt))) {
        Event *event = getEvent(segment, flower);
        commonAncestorEvent = commonAncestorEvent == NULL ? event : eventTree_getCommonAncestor(event, commonAncestorEvent);
    }
    assert(commonAncestorEvent != NULL);
    float treeCoverage = 0.0;
    stHash *hash = stHash_construct();

    segmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    while ((segment = stPinchBlockIt_getNext(&segmentIt))) {
        Event *event = getEvent(segment, flower);
        while (event != commonAncestorEvent && stHash_search(hash, event) == NULL) {
            treeCoverage += event_getBranchLength(event);
            stHash_insert(hash, event, event);
            event = event_getParent(event);
        }
    }

    float wholeTreeCoverage = event_getSubTreeBranchLength(event_getChild(eventTree_getRootEvent(eventTree), 0));
    assert(wholeTreeCoverage >= 0.0);
    if (wholeTreeCoverage <= 0.0) { //deal with case all leaf branches are not empty.
        return 0.0;
    }
    treeCoverage /= wholeTreeCoverage;
    assert(treeCoverage >= -0.001);
    assert(treeCoverage <= 1.0001);
    return treeCoverage;
}
