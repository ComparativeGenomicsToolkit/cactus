#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

///////////////////////////////////////////////////////////////////////////
// Core functions for melting
///////////////////////////////////////////////////////////////////////////

static void processChain(stCactusEdgeEnd *cactusEdgeEnd, void (*edgeEndFn)(stPinchBlock *, void *), void *extraArg) {
    while (1) {
        stPinchEnd *pinchEnd = stCactusEdgeEnd_getObject(cactusEdgeEnd);
        assert(pinchEnd != NULL);
        stPinchBlock *pinchBlock = stPinchEnd_getBlock(pinchEnd);
        assert(pinchBlock != NULL);
        edgeEndFn(pinchBlock, extraArg);
        cactusEdgeEnd = stCactusEdgeEnd_getNextEdgeEnd(cactusEdgeEnd);
        if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd)) {
            break;
        }
        cactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd);
    }
}

static void addBlock(stPinchBlock *block, void *extraArg) {
    stList_append(extraArg, block);
}

static void addChainBlocksToBlocksToDelete(stCactusEdgeEnd *cactusEdgeEnd, stList *blocksToDelete) {
    processChain(cactusEdgeEnd, addBlock, blocksToDelete);
}

static void addLength(stPinchBlock *block, void *extraArg) {
    *((int64_t *)extraArg) += stPinchBlock_getLength(block);
}

static int64_t getChainLength(stCactusEdgeEnd *cactusEdgeEnd) {
    int64_t length;
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


static bool isThreadEnd(stPinchBlock *pinchBlock) {
    stPinchSegment *pinchSegment = stPinchBlock_getFirst(pinchBlock);
    return pinchSegment != NULL &&
    (stPinchSegment_get3Prime(pinchSegment) == NULL || stPinchSegment_get5Prime(pinchSegment) == NULL);
}

static void trimAlignments(stPinchThreadSet *threadSet, int32_t blockEndTrim) {
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *block = stPinchThreadSetBlockIt_getNext(&blockIt);
    while (block != NULL) {
        stPinchBlock *block2 = stPinchThreadSetBlockIt_getNext(&blockIt);
        if(!isThreadEnd(block)) {
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
}

///////////////////////////////////////////////////////////////////////////
// Functions for calculating required species
///////////////////////////////////////////////////////////////////////////

bool stCaf_containsRequiredSpecies(stPinchBlock *pinchBlock, Flower *flower, int32_t requiredIngroupSpecies,
        int32_t requiredOutgroupSpecies, int32_t requiredAllSpecies) {
    int32_t outgroupSequences = 0;
    int32_t ingroupSequences = 0;
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        Event *event = cap_getEvent(flower_getCap(flower, stPinchSegment_getName(segment)));
        assert(event != NULL);
        if (event_isOutgroup(event)) {
            outgroupSequences++;
        } else {
            ingroupSequences++;
        }
    }
    return ingroupSequences >= requiredIngroupSpecies && outgroupSequences >= requiredOutgroupSpecies
            && outgroupSequences + ingroupSequences >= requiredAllSpecies;
}

void stCaf_calculateRequiredFractionsOfSpecies(Flower *flower,
        float requiredIngroupFraction, float requiredOutgroupFraction, float requiredAllFraction,
        int32_t *requiredOutgroupSpecies, int32_t *requiredIngroupSpecies, int32_t *requiredAllSpecies) {
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
