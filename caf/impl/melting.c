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
} //Adding dummy comment

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

static void addLength(stPinchBlock *block, void *extraArg) {
    *((int64_t *) extraArg) += stPinchBlock_getLength(block);
}

static int64_t getChainLength(stCactusEdgeEnd *cactusEdgeEnd) {
    int64_t length = 0;
    processChain(cactusEdgeEnd, addLength, &length, 0);
    return length;
}

stList *stCaf_getBlocksInChainsLessThanGivenLength(stCactusGraph *cactusGraph, int64_t minimumChainLength) {
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

static int64_t pathLength(stList *blocks) {
    int64_t accum = 0;
    for (int64_t i = 0; i < stList_length(blocks); i++) {
        accum += stPinchBlock_getLength(stList_get(blocks, i));
    }
    return accum;
}

void stCaf_undoChainsSmallerThanThis_preserveNonUndoableChains(stOnlineCactus *cactus, stPinchThreadSet *threadSet,
                                                               stList *pinches, stPinchUndo *undo,
                                                               int64_t minimumChainLength) {
    stList *worstPath = stOnlineCactus_getGloballyWorstMaximalChainOrBridgePath(cactus);
    while (pathLength(worstPath) < minimumChainLength) {
        // Need to undo stuff.
        // Go through all the blocks we pinched and find the one with the lowest chain (or maximal bridge-path) length.
        int64_t worstUndoablePathScore = INT64_MAX;
        stPinchBlock *worstUndoableBlock = NULL;
        for (int64_t i = 0; i < stList_length(pinches); i++) {
            stPinch *pinch = stList_get(pinches, i);
            stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
            stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
            stPinchSegment *segment = stPinchThread_getSegment(thread1, pinch->start1);
            int64_t pinchEnd = pinch->start1 + pinch->length;
            bool finishedThread1 = false;
            while (true) {
                stPinchBlock *block = stPinchSegment_getBlock(segment);
                int64_t nil; // for ignoring value of offset & length
                if (block != NULL && stPinchUndo_findOffsetForBlock(undo, threadSet, block,
                                                                    &nil, &nil)) {
                    stList *chainOrBridgePath = stOnlineCactus_getMaximalChainOrBridgePath(cactus, block);
                    int64_t pathScore = pathLength(chainOrBridgePath);
                    if (pathScore < worstUndoablePathScore || (pathScore == worstUndoablePathScore && stPinchBlock_getLength(block) < stPinchBlock_getLength(worstUndoableBlock))) {
                        worstUndoablePathScore = pathScore;
                        worstUndoableBlock = block;
                    }
                }
                segment = stPinchSegment_get3Prime(segment);
                if ((segment == NULL || (stPinchSegment_getThread(segment) == thread1  && stPinchSegment_getStart(segment) >= pinchEnd)) && !finishedThread1) {
                    // Switch to thread2 (or the second section of thread1). (We have to check both
                    // threads, in case one has collapsed in on
                    // itself.)
                    segment = stPinchThread_getSegment(thread2, pinch->start2);
                    pinchEnd = pinch->start2 + pinch->length;
                    finishedThread1 = true;
                } else if (segment == NULL || stPinchSegment_getStart(segment) >= pinchEnd) {
                    assert(segment == NULL || stPinchSegment_getThread(segment) == thread2);
                    break;
                }
            }
        }
        if (worstUndoableBlock == NULL) {
            st_errAbort("Couldn't find an undoable block, and there are still short chains");
        }

        // Undo this block.
        int64_t undoOffset;
        int64_t undoLength;
        stPinchUndo_findOffsetForBlock(undo, threadSet, worstUndoableBlock, &undoOffset, &undoLength);
        stPinchThreadSet_partiallyUndoPinch(threadSet, undo, undoOffset, undoLength);

        worstPath = stOnlineCactus_getGloballyWorstMaximalChainOrBridgePath(cactus);
    }
}

void stCaf_undoChainsSmallerThanThis_onlyUndo(stOnlineCactus *cactus, stPinchThreadSet *threadSet,
                                              stList *pinches, stPinchUndo *undo,
                                              int64_t minimumChainLength) {
    for (;;) {
        // Go through all the blocks we pinched and find the one with the lowest chain (or maximal bridge-path) length.
        int64_t worstUndoablePathScore = INT64_MAX;
        stPinchBlock *worstUndoableBlock = NULL;
        for (int64_t i = 0; i < stList_length(pinches); i++) {
            stPinch *pinch = stList_get(pinches, i);
            stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
            stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
            stPinchSegment *segment = stPinchThread_getSegment(thread1, pinch->start1);
            int64_t pinchEnd = pinch->start1 + pinch->length;
            bool finishedThread1 = false;
            while (true) {
                stPinchBlock *block = stPinchSegment_getBlock(segment);
                int64_t nil; // for ignoring value of offset & length
                if (block != NULL && stPinchUndo_findOffsetForBlock(undo, threadSet, block,
                                                                    &nil, &nil)) {
                    stList *chainOrBridgePath = stOnlineCactus_getMaximalChainOrBridgePath(cactus, block);
                    int64_t pathScore = pathLength(chainOrBridgePath);
                    if (pathScore < worstUndoablePathScore || (pathScore == worstUndoablePathScore && stPinchBlock_getLength(block) < stPinchBlock_getLength(worstUndoableBlock))) {
                        worstUndoablePathScore = pathScore;
                        worstUndoableBlock = block;
                    }
                }
                segment = stPinchSegment_get3Prime(segment);
                if ((segment == NULL || (stPinchSegment_getThread(segment) == thread1  && stPinchSegment_getStart(segment) >= pinchEnd)) && !finishedThread1) {
                    // Switch to thread2 (or the second section of thread1). (We have to check both
                    // threads, in case one has collapsed in on
                    // itself.)
                    segment = stPinchThread_getSegment(thread2, pinch->start2);
                    pinchEnd = pinch->start2 + pinch->length;
                    finishedThread1 = true;
                } else if (segment == NULL || stPinchSegment_getStart(segment) >= pinchEnd) {
                    assert(segment == NULL || stPinchSegment_getThread(segment) == thread2);
                    break;
                }
            }
        }

        if (worstUndoablePathScore < minimumChainLength) {
            // Undo this block.
            int64_t undoOffset;
            int64_t undoLength;
            stPinchUndo_findOffsetForBlock(undo, threadSet, worstUndoableBlock, &undoOffset, &undoLength);
            stPinchThreadSet_partiallyUndoPinch(threadSet, undo, undoOffset, undoLength);
        } else {
            break;
        }
    }
}

void stCaf_undoChainsSmallerThanThis_removeNonUndoableChains(stOnlineCactus *cactus, stPinchThreadSet *threadSet,
                                                             stList *pinches, stPinchUndo *undo,
                                                             int64_t minimumChainLength) {
    stList *worstPath = stOnlineCactus_getGloballyWorstMaximalChainOrBridgePath(cactus);
    while (pathLength(worstPath) < minimumChainLength) {
        // Need to undo stuff.
        // Go through all the blocks we pinched and find the one with the lowest chain (or maximal bridge-path) length.
        bool foundBlock = false;
        for (int64_t i = 0; i < stList_length(worstPath); i++) {
            stPinchBlock *block = stList_get(worstPath, i);
            // Undo this block if possible.
            int64_t undoOffset;
            int64_t undoLength;

            if (stPinchUndo_findOffsetForBlock(undo, threadSet, block, &undoOffset, &undoLength)) {
                stPinchThreadSet_partiallyUndoPinch(threadSet, undo, undoOffset, undoLength);
                foundBlock = true;
                break;
            }
        }
        if (!foundBlock) {
            // Need to copy the list, as it will be invalidated after a block destruct.
            stList *worstPathCopy = stList_construct();
            stList_appendAll(worstPathCopy, worstPath);
            for (int64_t i = 0; i < stList_length(worstPathCopy); i++) {
                stPinchBlock *block = stList_get(worstPathCopy, i);
                stPinchBlock_destruct(block);
            }
        }
        worstPath = stOnlineCactus_getGloballyWorstMaximalChainOrBridgePath(cactus);
    }
}

void stCaf_undoChainsSmallerThanThis_onlyRemove(stOnlineCactus *cactus, stPinchThreadSet *threadSet,
                                                stList *pinches, stPinchUndo *undo,
                                                int64_t minimumChainLength) {
    stList *worstPath = stOnlineCactus_getGloballyWorstMaximalChainOrBridgePath(cactus);
    while (pathLength(worstPath) < minimumChainLength) {
        // Need to copy the list, as it will be invalidated after a block destruct.
        stList *worstPathCopy = stList_construct();
        stList_appendAll(worstPathCopy, worstPath);
        for (int64_t i = 0; i < stList_length(worstPathCopy); i++) {
            stPinchBlock *block = stList_get(worstPathCopy, i);
            stPinchBlock_destruct(block);
        }
        worstPath = stOnlineCactus_getGloballyWorstMaximalChainOrBridgePath(cactus);
    }
}

void stCaf_melt(Flower *flower, stPinchThreadSet *threadSet, bool blockFilterfn(stPinchBlock *), int64_t blockEndTrim,
        int64_t minimumChainLength, bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds) {
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
        stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, INT64_MAX,
                0.0, breakChainsAtReverseTandems, maximumMedianSpacingBetweenLinkedEnds);
        stList *blocksToDelete = stCaf_getBlocksInChainsLessThanGivenLength(cactusGraph, minimumChainLength);

        printf("A melting round is destroying %" PRIi64 " blocks with an average degree "
               "of %lf from chains with length less than %" PRIi64 ". Total aligned bases"
               " lost: %" PRIu64 "\n",
               stList_length(blocksToDelete), stCaf_averageBlockDegree(blocksToDelete),
               minimumChainLength, stCaf_totalAlignedBases(blocksToDelete));

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

Event *getEvent(stPinchSegment *segment, Flower *flower) {
    Event *event = cap_getEvent(flower_getCap(flower, stPinchSegment_getName(segment)));
    assert(event != NULL);
    return event;
}

bool stCaf_containsRequiredSpecies(stPinchBlock *pinchBlock, Flower *flower, int64_t minimumIngroupDegree,
        int64_t minimumOutgroupDegree, int64_t requiredAllSpecies) {
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
    return ingroupSequences >= minimumIngroupDegree && outgroupSequences >= minimumOutgroupDegree && outgroupSequences
            + ingroupSequences >= requiredAllSpecies;
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

///////////////////////////////////////////////////////////////////////////
// Misc. functions
///////////////////////////////////////////////////////////////////////////

double stCaf_averageBlockDegree(stList *blocks) {
    if (stList_length(blocks) == 0) {
        return 0.0;
    }
    uint64_t total = 0;
    for (int64_t i = 0; i < stList_length(blocks); i++) {
        total += stPinchBlock_getDegree(stList_get(blocks, i));
    }
    return ((double) total) / stList_length(blocks);
}

uint64_t stCaf_totalAlignedBases(stList *blocks) {
    uint64_t ret = 0;
    for (int64_t i = 0; i < stList_length(blocks); i++) {
        stPinchBlock *block = stList_get(blocks, i);
        ret += stPinchBlock_getDegree(block) * stPinchBlock_getLength(block);
    }
    return ret;
}
