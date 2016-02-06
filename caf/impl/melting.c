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
        //Cleanup cactus
        stCactusGraph_destruct(cactusGraph);
        stList_destruct(blocksToDelete); //This will destroy the blocks
    }
    //Now heal up the trivial boundaries
    stCaf_joinTrivialBoundaries(threadSet);
}

Flower *debugFlower;

static char *getEndStr(stPinchEnd *end) {
    stPinchBlockIt it = stPinchBlock_getSegmentIterator(end->block);
    stList *segmentNames = stList_construct3(0, free);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&it)) != NULL) {
        Cap *cap = flower_getCap(debugFlower, stPinchSegment_getName(segment));
        const char *header = sequence_getHeader(cap_getSequence(cap));
        const char *genome = event_getHeader(cap_getEvent(cap));
        stList_append(segmentNames, stString_print("%s.%s|%" PRIi64 "-%" PRIi64, genome, header, stPinchSegment_getStart(segment), stPinchSegment_getStart(segment) + stPinchSegment_getLength(segment)));
    }
    char *ret = stString_join2(",", segmentNames);
    stList_destruct(segmentNames);
    return ret;
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

// Determine whether the chain is recoverable (i.e. will bar phase be
// expected to pick it back up?).
static bool chainIsRecoverable(stCactusEdgeEnd *chainEnd, stSet *deadEndComponent) {
    stPinchEnd *end1 = stCactusEdgeEnd_getObject(chainEnd);
    stPinchEnd *end2 = stCactusEdgeEnd_getObject(stCactusEdgeEnd_getLink(chainEnd));

    stSet *connectedEnds1 = stPinchEnd_getConnectedPinchEnds(end1);
    stSet *connectedEnds2 = stPinchEnd_getConnectedPinchEnds(end2);

    printf("c1: %s c2: %s\n", getEndStr(end1), getEndStr(end2));
    printf("ends connected to c1:\n");
    stSetIterator *it = stSet_getIterator(connectedEnds1);
    stPinchEnd *end;
    while ((end = stSet_getNext(it)) != NULL) {
        printf("%s\n", getEndStr(end));
    }
    stSet_destructIterator(it);
    printf("ends connected to c2:\n");
    it = stSet_getIterator(connectedEnds2);
    while ((end = stSet_getNext(it)) != NULL) {
        printf("%s\n", getEndStr(end));
    }
    stSet_destructIterator(it);

    if (isTelomere(end1, deadEndComponent) || isTelomere(end2, deadEndComponent)) {
        // Chain containing only the telomere/stub end
        printf("telomere\n");
        return false;
    }

    if (endsDoNotHaveSameThreadComposition(end1, end2)) {
        // One or more of the threads ran into a stub end and
        // appeared/disappeared partway through the chain
        printf("stubbed chain\n");
        return false;
    }


    stSet *sharedEnds = stSet_getIntersection(connectedEnds1, connectedEnds2);

    bool recoverable = true;
    if (stSet_size(sharedEnds) != 0) {
        // The two ends link to the same end.
        recoverable = false;
    } else if (stSet_size(connectedEnds1) != 1 && stSet_size(connectedEnds2) != 1) {
        // Both ends link to more than one end.
        recoverable = false;
    } else if (endSetContainsTelomere(connectedEnds1, deadEndComponent) || endSetContainsTelomere(connectedEnds2, deadEndComponent)) {
        // Connected to one or more attached ends or stub ends.
        recoverable = false;
    } else if (stSet_search(connectedEnds1, end2)) {
        // A duplication (link connecting the two child chain ends).
        assert(stSet_search(connectedEnds2, end1));
        recoverable = false;
    }

    stSet_destruct(sharedEnds);
    stSet_destruct(connectedEnds1);
    stSet_destruct(connectedEnds2);
    printf("recoverable: %d\n", recoverable);
    return recoverable;
}

// FIXME: remove
static int64_t numNodesVisited = 0;
static int64_t numChainsVisited = 0;
static int64_t numRecoverableChains = 0;

// For a given cactus node, recurse through all nodes below it and
// find recoverable chains below them. Then find recoverable chains
// below the current node given its parent chain.
static void getRecoverableChains_R(stCactusNode *cactusNode, stCactusEdgeEnd *parentChain, stSet *deadEndComponent, bool (*chainFilter)(stCactusEdgeEnd *), stList *recoverableChains) {
    numNodesVisited++;
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
                                   chainFilter,
                                   recoverableChains);
        }
    }

    if (parentChain != NULL) {
        // Visit the next node on this chain (unless it's where we started).
        stCactusEdgeEnd *nextEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd_getLink(parentChain));
        if (!stCactusEdgeEnd_isChainEnd(nextEdgeEnd)) {
            getRecoverableChains_R(stCactusEdgeEnd_getNode(nextEdgeEnd), nextEdgeEnd, deadEndComponent, chainFilter, recoverableChains);
        }
    }

    cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
    while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt)) != NULL) {
        if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) {
            if ((chainFilter == NULL || chainFilter(cactusEdgeEnd)) && chainIsRecoverable(cactusEdgeEnd, deadEndComponent)) {
                numRecoverableChains++;
                stList_append(recoverableChains, cactusEdgeEnd);
            }
            numChainsVisited++;
        }
    }
}

static stList *getRecoverableChains(stCactusNode *startCactusNode, stSet *deadEndComponent, bool (*chainFilter)(stCactusEdgeEnd *)) {
    stList *ret = stList_construct();
    getRecoverableChains_R(startCactusNode, NULL, deadEndComponent, chainFilter, ret);
    printf("Visited %" PRIi64 " cactus nodes while getting recoverable chains\n", numNodesVisited);
    printf("Found %" PRIi64 " / %" PRIi64 " recoverable chains\n", numRecoverableChains, numChainsVisited);
    return ret;
}

void stCaf_meltRecoverableChains(Flower *flower, stPinchThreadSet *threadSet, bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds, bool (*chainFilter)(stCactusEdgeEnd *)) {
    debugFlower = flower;
    stCactusNode *startCactusNode;
    stList *deadEndComponent;
    stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, INT64_MAX,
                                                                  0.0, breakChainsAtReverseTandems, maximumMedianSpacingBetweenLinkedEnds);

    // Construct a queryable set of stub ends.
    stSet *deadEndComponentSet = stSet_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, NULL);
    for (int64_t i = 0; i < stList_length(deadEndComponent); i++) {
        stSet_insert(deadEndComponentSet, stList_get(deadEndComponent, i));
    }

    stList *recoverableChains = getRecoverableChains(startCactusNode, deadEndComponentSet, chainFilter);
    stList *blocksToDelete = stList_construct3(0, (void(*)(void *)) stPinchBlock_destruct);
    for (int64_t i = 0; i < stList_length(recoverableChains); i++) {
        stCactusEdgeEnd *chainEnd = stList_get(recoverableChains, i);
        addChainBlocksToBlocksToDelete(chainEnd, blocksToDelete);
    }
    printf("Destroying %" PRIi64 " recoverable blocks\n", stList_length(blocksToDelete));
    stList_destruct(recoverableChains);
    stList_destruct(blocksToDelete);

    // FIXME: remove
    stCactusGraphNodeIt *it = stCactusGraphNodeIterator_construct(cactusGraph);
    stCactusNode *node;
    int64_t numCactusNodes = 0, numChains = 0;
    while ((node = stCactusGraphNodeIterator_getNext(it)) != NULL) {
        numCactusNodes++;
        numChains += stCactusNode_getChainNumber(node);
    }
    printf("There were actually %" PRIi64 " nodes and %" PRIi64 " chains in the graph\n", numCactusNodes, numChains);
    stCactusGraphNodeIterator_destruct(it);

    stCactusGraph_destruct(cactusGraph);
    stSet_destruct(deadEndComponentSet);
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
