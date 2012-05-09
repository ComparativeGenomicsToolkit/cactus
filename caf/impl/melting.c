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
    bool threadEnd = pinchSegment != NULL && (stPinchSegment_get3Prime(pinchSegment) == NULL || stPinchSegment_get5Prime(pinchSegment) == NULL);
    if(threadEnd) {
        assert(stPinchBlock_getLength(pinchBlock) == 1);
    }
    return threadEnd;
}

static void processChain(stCactusEdgeEnd *cactusEdgeEnd, void(*edgeEndFn)(stPinchBlock *, void *), void *extraArg) {
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
        assert(stCactusEdgeEnd_getLink(stCactusEdgeEnd_getLink(cactusEdgeEnd)) == cactusEdgeEnd);
        cactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd);
    }
}

static void addBlock(stPinchBlock *block, void *extraArg) {
    if(!isThreadEnd(block)) {
        stList_append(extraArg, block);
    }
}

static void addChainBlocksToBlocksToDelete(stCactusEdgeEnd *cactusEdgeEnd, stList *blocksToDelete) {
    processChain(cactusEdgeEnd, addBlock, blocksToDelete);
}

static void addLength(stPinchBlock *block, void *extraArg) {
    *((int64_t *) extraArg) += stPinchBlock_getLength(block);
}

static int64_t getChainLength(stCactusEdgeEnd *cactusEdgeEnd) {
    int64_t length = 0;
    processChain(cactusEdgeEnd, addLength, &length);
    return length;
}

static stList *stCaf_getBlocksInChainsLessThanGivenLength(stCactusGraph *cactusGraph, int32_t minimumChainLength) {
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

static void trimAlignments(stPinchThreadSet *threadSet, int32_t blockEndTrim) {
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

void stCaf_melt(Flower *flower, stPinchThreadSet *threadSet, bool blockFilterfn(stPinchBlock *), int32_t blockEndTrim,
        int64_t minimumChainLength) {
    //First trim
    if (blockEndTrim > 0) {
        trimAlignments(threadSet, blockEndTrim);
    }
    //Then filter blocks
    if (blockFilterfn != NULL) {
        filterAlignments(threadSet, blockFilterfn);
    }
    //Now apply the minimum chain length filter
    if (minimumChainLength > 1) {
        stCactusNode *startCactusNode;
        stList *deadEndComponent;
        stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0);
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
        float requiredAllFraction, int32_t *requiredIngroupSpecies,  int32_t *requiredOutgroupSpecies, int32_t *requiredAllSpecies) {
    if (requiredIngroupFraction <= 0.0 && requiredOutgroupFraction <= 0.0 && requiredAllFraction <= 0.0) {
        *requiredAllSpecies = 0;
        *requiredOutgroupSpecies = 0;
        *requiredAllSpecies = 0;
    }
    EventTree *eventTree = flower_getEventTree(flower);
    Event *event;
    int32_t outgroupEventNumber = 0;
    int32_t ingroupEventNumber = 0;
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
    *requiredOutgroupSpecies = outgroupEventNumber * requiredOutgroupFraction;
    *requiredIngroupSpecies = ingroupEventNumber * requiredIngroupFraction;
    *requiredAllSpecies = (ingroupEventNumber + outgroupEventNumber) * requiredAllFraction;
}

static Event *getEvent(stPinchSegment *segment, Flower *flower) {
    Event *event = cap_getEvent(flower_getCap(flower, stPinchSegment_getName(segment)));
    assert(event != NULL);
    return event;
}

bool stCaf_containsRequiredSpecies(stPinchBlock *pinchBlock, Flower *flower, int32_t requiredIngroupSpecies,
        int32_t requiredOutgroupSpecies, int32_t requiredAllSpecies) {
    if (requiredIngroupSpecies <= 0 && requiredOutgroupSpecies <= 0 && requiredAllSpecies <= 0) {
        return 1;
    }
    int32_t outgroupSequences = 0;
    int32_t ingroupSequences = 0;
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
