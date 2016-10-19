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

static bool isTelomere(stPinchEnd *end, stSet *deadEndComponent) {
    stPinchSegment *segment = stPinchBlock_getFirst(end->block);
    bool atEndOfThread = stPinchThread_getFirst(stPinchSegment_getThread(segment)) == segment || stPinchThread_getLast(stPinchSegment_getThread(segment)) == segment;
    bool inDeadEndComponent = stSet_search(deadEndComponent, end);
    return atEndOfThread || inDeadEndComponent;
}

static bool endSetContainsTelomere(stSet *endSet, stSet *deadEndComponent) {
    stSetIterator *it = stSet_getIterator(endSet);
    bool containsTelomere = false;
    stPinchEnd *end;
    while ((end = stSet_getNext(it)) != NULL) {
        if (isTelomere(end, deadEndComponent)) {
            containsTelomere = true;
            break;
        }
    }
    stSet_destructIterator(it);
    return containsTelomere;
}

static bool endsDoNotHaveSameThreadComposition(stPinchEnd *end1, stPinchEnd *end2) {
    if (stPinchBlock_getDegree(end1->block) != stPinchBlock_getDegree(end2->block)) {
        return true;
    }
    stPinchBlockIt it1 = stPinchBlock_getSegmentIterator(end1->block);
    stSet *threads1 = stSet_construct();
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&it1)) != NULL) {
        stSet_insert(threads1, stPinchSegment_getThread(segment));
    }

    stSet *threads2 = stSet_construct();
    stPinchBlockIt it2 = stPinchBlock_getSegmentIterator(end2->block);
    while ((segment = stPinchBlockIt_getNext(&it2)) != NULL) {
        stSet_insert(threads2, stPinchSegment_getThread(segment));
    }

    bool sameThreadComposition = true;
    stSet *intersection = stSet_getIntersection(threads1, threads2);
    if (stSet_size(intersection) != stSet_size(threads1) || stSet_size(intersection) != stSet_size(threads2)) {
        sameThreadComposition = false;
    }

    stSet_destruct(threads1);
    stSet_destruct(threads2);
    stSet_destruct(intersection);
    return !sameThreadComposition;
}

static bool chainConnectsToTelomere(stCactusEdgeEnd *chainEnd, stSet *deadEndComponent) {
    stPinchEnd *end1 = stCactusEdgeEnd_getObject(chainEnd);
    stPinchEnd *end2 = stCactusEdgeEnd_getObject(stCactusEdgeEnd_getLink(chainEnd));

    if (endsDoNotHaveSameThreadComposition(end1, end2)) {
        // One or more of the threads ran into a stub end and
        // appeared/disappeared partway through the chain
        return true;
    }

    stSet *connectedEnds1 = stPinchEnd_getConnectedPinchEnds(end1);
    stSet *connectedEnds2 = stPinchEnd_getConnectedPinchEnds(end2);

    bool connectedToTelomere = false;
    if (endSetContainsTelomere(connectedEnds1, deadEndComponent) ||
        endSetContainsTelomere(connectedEnds2, deadEndComponent)) {
        // Connected to one or more attached ends or stub ends.
        connectedToTelomere = true;
    }

    stSet_destruct(connectedEnds1);
    stSet_destruct(connectedEnds2);
    return connectedToTelomere;
}

// Determine whether the chain is recoverable (i.e. will bar phase be
// expected to pick it back up?).
static bool chainIsRecoverable(stCactusEdgeEnd *chainEnd, stSet *deadEndComponent) {
    stPinchEnd *end1 = stCactusEdgeEnd_getObject(chainEnd);
    stPinchEnd *end2 = stCactusEdgeEnd_getObject(stCactusEdgeEnd_getLink(chainEnd));

    stSet *connectedEnds1 = stPinchEnd_getConnectedPinchEnds(end1);
    stSet *connectedEnds2 = stPinchEnd_getConnectedPinchEnds(end2);

    stSet *sharedEnds = stSet_getIntersection(connectedEnds1, connectedEnds2);

    bool recoverable = true;
    if (isTelomere(end1, deadEndComponent) || isTelomere(end2, deadEndComponent)) {
        // Chain containing only the telomere/stub end
        recoverable =  false;
    } else if (stSet_size(sharedEnds) != 0) {
        // The two ends link to the same end.
        recoverable = false;
    } else if (stSet_size(connectedEnds1) != 1 && stSet_size(connectedEnds2) != 1) {
        // Both ends link to more than one end.
        recoverable = false;
    } else if (stSet_search(connectedEnds1, end2)) {
        // A duplication (link connecting the two child chain ends).
        assert(stSet_search(connectedEnds2, end1));
        recoverable = false;
    }

    stSet_destruct(sharedEnds);
    stSet_destruct(connectedEnds1);
    stSet_destruct(connectedEnds2);
    return recoverable;
}

// Abstracts out getting the only corresponding chain end from a set of pinch ends of size 1.
static stCactusEdgeEnd *getChainEndFromSingletonSet(stSet *ends,
                                                    stHash *pinchEndToChainEnd) {
    assert(stSet_size(ends) == 1);
    stSetIterator *it = stSet_getIterator(ends);
    stPinchEnd *connectedPinchEnd = stSet_getNext(it);
    stSet_destructIterator(it);

    stCactusEdgeEnd *chainEnd = stHash_search(pinchEndToChainEnd, connectedPinchEnd);
    assert(chainEnd != NULL);
    if (!stCactusEdgeEnd_getLinkOrientation(chainEnd)) {
        chainEnd = stCactusEdgeEnd_getLink(chainEnd);
    }
    return chainEnd;
}

// Mark down which chain(s) this (recoverable) chain is recoverable given.
static void markRecoverableAdjacencies(stCactusEdgeEnd *recoverableChainEnd,
                                       stHash *pinchEndToChainEnd,
                                       stHash *chainToRecoverableAdjacencies) {
    stPinchEnd *end1 = stCactusEdgeEnd_getObject(recoverableChainEnd);
    stPinchEnd *end2 = stCactusEdgeEnd_getObject(stCactusEdgeEnd_getLink(recoverableChainEnd));

    stSet *connectedEnds1 = stPinchEnd_getConnectedPinchEnds(end1);
    stSet *connectedEnds2 = stPinchEnd_getConnectedPinchEnds(end2);

    stList *recoverableAdjacencies = stList_construct();
    // We can safely assume there are no shared ends since the chain
    // is known to be recoverable. So all we have to check for is that
    // there is only one connected end. If so, this chain is
    // recoverable given the other.
    if (stSet_size(connectedEnds1) == 1) {
        stCactusEdgeEnd *connectedChainEnd = getChainEndFromSingletonSet(connectedEnds1,
                                                                         pinchEndToChainEnd);
        stList_append(recoverableAdjacencies, connectedChainEnd);
    }

    if (stSet_size(connectedEnds2) == 1) {
        stCactusEdgeEnd *connectedChainEnd = getChainEndFromSingletonSet(connectedEnds2,
                                                                         pinchEndToChainEnd);
        stList_append(recoverableAdjacencies, connectedChainEnd);
    }

    stHash_insert(chainToRecoverableAdjacencies, recoverableChainEnd, recoverableAdjacencies);

    stSet_destruct(connectedEnds1);
    stSet_destruct(connectedEnds2);
}

/*
 * Get a mapping from pinch ends to the canonical chain end for their chain.
 */

static stHash *getPinchEndToChainEndHash(stCactusGraph *cactusGraph) {
    stHash *pinchEndToChainEnd = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, NULL, NULL);
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(cactusGraph);
    stCactusNode *cactusNode;
    while ((cactusNode = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stCactusNodeEdgeEndIt cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
        stCactusEdgeEnd *cactusEdgeEnd;
        while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt)) != NULL) {
            if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) {
                // This is the canonical end for an unvisited chain. We
                // iterate over all the ends in the chain, mapping them to
                // this canonical end.
                stCactusEdgeEnd *chainEnd = cactusEdgeEnd;
                stCactusEdgeEnd *curEnd = cactusEdgeEnd;
                do {
                    stPinchEnd *pinchEnd = stCactusEdgeEnd_getObject(curEnd);
                    stHash_insert(pinchEndToChainEnd, pinchEnd, chainEnd);
                    if (stCactusEdgeEnd_getLinkOrientation(curEnd)) {
                        curEnd = stCactusEdgeEnd_getLink(curEnd);
                    } else {
                        curEnd = stCactusEdgeEnd_getOtherEdgeEnd(curEnd);
                    }
                } while (curEnd != chainEnd);
            }
        }
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
    return pinchEndToChainEnd;
}

// For a given cactus node, recurse through all nodes below it and
// find recoverable chains below them. Then find recoverable chains
// below the current node given its parent chain.
static void getRecoverableChains_R(stCactusNode *cactusNode, stCactusEdgeEnd *parentChain, stSet *deadEndComponent, Flower *flower, bool (*recoverabilityFilter)(stCactusEdgeEnd *, Flower *), stHash *pinchEndToChainEnd, stSet *recoverableChains, stList *telomereAdjacentChains, stHash *chainToRecoverableAdjacencies) {
    stCactusNodeEdgeEndIt cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
    stCactusEdgeEnd *cactusEdgeEnd;
    while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt)) != NULL) {
        if ((parentChain == NULL
             || (cactusEdgeEnd != parentChain
                 && cactusEdgeEnd != stCactusEdgeEnd_getLink(parentChain)))
            && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)
            && stCactusEdgeEnd_getOtherNode(cactusEdgeEnd) != cactusNode) {
            // Found a new chain below this node.
            assert(stCactusEdgeEnd_isChainEnd(cactusEdgeEnd));
            getRecoverableChains_R(stCactusEdgeEnd_getOtherNode(cactusEdgeEnd),
                                   stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd),
                                   deadEndComponent,
                                   flower,
                                   recoverabilityFilter,
                                   pinchEndToChainEnd,
                                   recoverableChains,
                                   telomereAdjacentChains,
                                   chainToRecoverableAdjacencies);
        }
    }

    if (parentChain != NULL) {
        // Visit the next node on this chain (unless it's where we started).
        stCactusEdgeEnd *nextEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd_getLink(parentChain));
        if (!stCactusEdgeEnd_isChainEnd(nextEdgeEnd)) {
            getRecoverableChains_R(stCactusEdgeEnd_getNode(nextEdgeEnd), nextEdgeEnd, deadEndComponent, flower, recoverabilityFilter, pinchEndToChainEnd, recoverableChains, telomereAdjacentChains, chainToRecoverableAdjacencies);
        }
    }

    cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
    while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt)) != NULL) {
        if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) {
            if ((recoverabilityFilter == NULL || recoverabilityFilter(cactusEdgeEnd, flower)) && chainIsRecoverable(cactusEdgeEnd, deadEndComponent)) {
                stSet_insert(recoverableChains, cactusEdgeEnd);
                markRecoverableAdjacencies(cactusEdgeEnd, pinchEndToChainEnd, chainToRecoverableAdjacencies);
                if (chainConnectsToTelomere(cactusEdgeEnd, deadEndComponent)) {
                    stList_append(telomereAdjacentChains, cactusEdgeEnd);
                }
            }
        }
    }
}

static stList *getRecoverableChains(stCactusGraph *cactusGraph, stCactusNode *startCactusNode, stSet *deadEndComponent, Flower *flower, bool (*recoverabilityFilter)(stCactusEdgeEnd *, Flower *)) {
    stHash *pinchEndToChainEnd = getPinchEndToChainEndHash(cactusGraph);

    stSet *recoverableChainSet = stSet_construct();
    stList *telomereAdjacentChains = stList_construct();
    stHash *chainToRecoverableAdjacencies = stHash_construct2(NULL, (void (*)(void *)) stList_destruct);
    getRecoverableChains_R(startCactusNode, NULL, deadEndComponent, flower, recoverabilityFilter, pinchEndToChainEnd, recoverableChainSet, telomereAdjacentChains, chainToRecoverableAdjacencies);

    // Remove anchors that are connected to telomeres and are not
    // transitively connected to an unrecoverable chain. This ensures
    // that we don't lose alignment by deeming all chains recoverable
    // and not keeping any anchors to recover them.
    for (int64_t i = 0; i < stList_length(telomereAdjacentChains); i++) {
        stCactusEdgeEnd *telomereAdjacentChain = stList_get(telomereAdjacentChains, i);
        stCactusEdgeEnd *curChain = telomereAdjacentChain;
        stCactusEdgeEnd *prevChain = NULL;
        bool neededAsAnchor = false;
        while (stSet_search(recoverableChainSet, curChain)) {
            stList *recoverableAdjacencies = stHash_search(chainToRecoverableAdjacencies, curChain);
            assert(stList_length(recoverableAdjacencies) > 0);
            assert(stList_length(recoverableAdjacencies) <= 2);
            bool foundValidAdjacency = false;
            for (int64_t j = 0; j < stList_length(recoverableAdjacencies); j++) {
                stCactusEdgeEnd *recoverableAdjacency = stList_get(recoverableAdjacencies, j);
                stPinchEnd *adjacencyEnd1 = stCactusEdgeEnd_getObject(recoverableAdjacency);
                stPinchEnd *adjacencyEnd2 = stCactusEdgeEnd_getObject(stCactusEdgeEnd_getLink(recoverableAdjacency));
                if (recoverableAdjacency != prevChain &&
                    !isTelomere(adjacencyEnd1, deadEndComponent) &&
                    !isTelomere(adjacencyEnd2, deadEndComponent)) {
                    prevChain = curChain;
                    curChain = recoverableAdjacency;
                    foundValidAdjacency = true;
                    break;
                }
            }
            if (!foundValidAdjacency) {
                neededAsAnchor = true;
                break;
            }
        }
        if (neededAsAnchor) {
            stSet_remove(recoverableChainSet, telomereAdjacentChain);
        }
    }
    stList_destruct(telomereAdjacentChains);
    stHash_destruct(chainToRecoverableAdjacencies);

    // Convert the recoverable chains set into a list.
    stList *recoverableChains = stList_construct();
    stSetIterator *it = stSet_getIterator(recoverableChainSet);
    stCactusEdgeEnd *chainEnd;
    while ((chainEnd = stSet_getNext(it)) != NULL) {
        stList_append(recoverableChains, chainEnd);
    }
    stSet_destructIterator(it);
    stSet_destruct(recoverableChainSet);
    stHash_destruct(pinchEndToChainEnd);
    return recoverableChains;
}

static int64_t numColumns(stList *blocks) {
    int64_t total = 0;
    for (int64_t i = 0; i < stList_length(blocks); i++) {
        stPinchBlock *block = stList_get(blocks, i);
        total += stPinchBlock_getLength(block);
    }
    return total;
}

static int64_t totalAlignedBases(stList *blocks) {
    int64_t total = 0;
    for (int64_t i = 0; i < stList_length(blocks); i++) {
        stPinchBlock *block = stList_get(blocks, i);
        total += stPinchBlock_getLength(block) * stPinchBlock_getDegree(block);
    }
    return total;
}

void stCaf_meltRecoverableChains(Flower *flower, stPinchThreadSet *threadSet, bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds, bool (*recoverabilityFilter)(stCactusEdgeEnd *, Flower *), int64_t maxNumIterations, int64_t maxRecoverableChainLength) {
    while (maxNumIterations-- > 0) {
        stCactusNode *startCactusNode;
        stList *deadEndComponent;
        // FIXME: We shouldn't really have to rebuild the cactus graph
        // every time. Instead we should just be able to do multiple
        // iterations over the same graph, keeping track of the chains
        // whose underlying blocks we've already deleted.
        stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, 0,
                                                                      0.0, breakChainsAtReverseTandems, maximumMedianSpacingBetweenLinkedEnds);

        // Construct a queryable set of stub ends.
        stSet *deadEndComponentSet = stSet_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, NULL);
        for (int64_t i = 0; i < stList_length(deadEndComponent); i++) {
            stSet_insert(deadEndComponentSet, stList_get(deadEndComponent, i));
        }

        stList *recoverableChains = getRecoverableChains(cactusGraph, startCactusNode, deadEndComponentSet, flower, recoverabilityFilter);

        stList *blocksToDelete = stList_construct3(0, (void(*)(void *)) stPinchBlock_destruct);
        for (int64_t i = 0; i < stList_length(recoverableChains); i++) {
            stCactusEdgeEnd *chainEnd = stList_get(recoverableChains, i);
            if (getChainLength(chainEnd) <= maxRecoverableChainLength) {
                addChainBlocksToBlocksToDelete(chainEnd, blocksToDelete);
            }
        }
        int64_t numRecoverableBlocks = stList_length(blocksToDelete);
        printf("Destroying %" PRIi64 " recoverable blocks\n", numRecoverableBlocks);
        printf("The blocks covered %" PRIi64 " columns for a total of %" PRIi64 " aligned bases\n", numColumns(blocksToDelete), totalAlignedBases(blocksToDelete));
        stList_destruct(recoverableChains);
        stList_destruct(blocksToDelete);

        stCactusGraph_destruct(cactusGraph);
        stSet_destruct(deadEndComponentSet);

        if (numRecoverableBlocks == 0) {
            // We didn't delete anything this round; we can safely
            // stop since we haven't changed the graph at all.
            break;
        }
    }
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
