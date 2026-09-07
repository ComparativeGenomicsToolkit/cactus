#include <stdlib.h>
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

static void filterAlignments(stPinchThreadSet *threadSet, bool(*blockFilterFn)(stPinchBlock *, void *extraArg),
                             void *extraArg) {
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *block = stPinchThreadSetBlockIt_getNext(&blockIt);
    while (block != NULL) {
        stPinchBlock *block2 = stPinchThreadSetBlockIt_getNext(&blockIt);
        if (!isThreadEnd(block) && blockFilterFn(block, extraArg)) {
            stPinchBlock_destruct(block);
        }
        block = block2;
    }
}

int64_t stCaf_melt(Flower *flower, stPinchThreadSet *threadSet, bool blockFilterfn(stPinchBlock *, void *extraArg),
                void *extraArg, int64_t blockEndTrim, int64_t minimumChainLength,
                bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds) {
    double t = stCaf_now(), trimTime = 0.0, filterTime = 0.0, graphTime = 0.0, scanTime = 0.0, deleteTime = 0.0;
    int64_t blocksDestroyed = 0;

    //First trim
    if (blockEndTrim > 0) {
        trimAlignments(threadSet, blockEndTrim);
        trimTime = stCaf_now() - t;
        t = stCaf_now();
    }

    //Then filter blocks
    if (blockFilterfn != NULL) {
        filterAlignments(threadSet, blockFilterfn, extraArg);
        filterTime = stCaf_now() - t;
        t = stCaf_now();
    }

    //Now apply the minimum chain length filter
    if (minimumChainLength > 1) {
        stCactusNode *startCactusNode;
        stPinchComponent *deadEndComponent;
        stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, INT64_MAX,
                0.0, breakChainsAtReverseTandems, maximumMedianSpacingBetweenLinkedEnds);
        graphTime = stCaf_now() - t;
        t = stCaf_now();
        stList *blocksToDelete = stCaf_getBlocksInChainsLessThanGivenLength(cactusGraph, minimumChainLength);
        scanTime = stCaf_now() - t;
        t = stCaf_now();
        blocksDestroyed = stList_length(blocksToDelete);

        st_logInfo("A melting round is destroying %" PRIi64 " blocks with an average degree "
               "of %lf from chains with length less than %" PRIi64 ". Total aligned bases"
               " lost: %" PRIu64 "\n",
               stList_length(blocksToDelete), stCaf_averageBlockDegree(blocksToDelete),
               minimumChainLength, stCaf_totalAlignedBases(blocksToDelete));

        //Cleanup cactus
        stCaf_destructCactusGraph(cactusGraph, threadSet);
        stList_destruct(blocksToDelete); //This will destroy the blocks
        deleteTime = stCaf_now() - t;
        t = stCaf_now();
    }
    //Now heal up the trivial boundaries
    stCaf_joinTrivialBoundaries(threadSet);
    st_logInfo("caf-timing: melt minChain=%" PRIi64 " trim %.3fs filter %.3fs graph %.3fs scan %.3fs delete %.3fs join %.3fs destroyed %" PRIi64 "\n",
               minimumChainLength, trimTime, filterTime, graphTime, scanTime, deleteTime, stCaf_now() - t, blocksDestroyed);
    return blocksDestroyed;
}

///////////////////////////////////////////////////////////////////////////
// Melting by chain length alone, with a join confined to what the deletions touched
///////////////////////////////////////////////////////////////////////////

/*
 * A deleted block's segment with the coordinates it had, so that a segment already absorbed into a
 * merged run (and freed) can be recognised and skipped without being read.
 */
typedef struct _deletedSegment {
    stPinchSegment *segment;
    int64_t threadIndex, start;
} DeletedSegment;

static int deletedSegment_cmp(const void *a, const void *b) {
    const DeletedSegment *x = a, *y = b;
    if (x->threadIndex != y->threadIndex) {
        return x->threadIndex < y->threadIndex ? -1 : 1;
    }
    return x->start < y->start ? -1 : (x->start > y->start ? 1 : 0);
}

/*
 * Destroys the blocks and merges every run of block-less segments the deletions created or extended,
 * leaving the graph as the per-thread pass of stPinchThreadSet_joinTrivialBoundaries would: each run
 * merged into its leftmost segment. Only the deleted segments and their runs are touched.
 */
static void destroyBlocksJoiningRuns(stList *blocksToDelete) {
    DeletedSegment stackSegments[64];
    for (int64_t i = 0; i < stList_length(blocksToDelete); i++) {
        stPinchBlock *block = stList_get(blocksToDelete, i);
        int64_t degree = stPinchBlock_getDegree(block);
        DeletedSegment *segments = degree <= 64 ? stackSegments : st_malloc(degree * sizeof(DeletedSegment));
        int64_t n = 0;
        stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(block);
        stPinchSegment *segment;
        while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
            segments[n].segment = segment;
            segments[n].threadIndex = stPinchThread_getIndex(stPinchSegment_getThread(segment));
            segments[n].start = stPinchSegment_getStart(segment);
            n++;
        }
        assert(n == degree);
        // In thread order, so that a segment absorbed into the run of an earlier one is recognised by
        // lying inside that run's extent
        qsort(segments, n, sizeof(DeletedSegment), deletedSegment_cmp);
        stPinchBlock_destruct(block);
        int64_t survivorThread = -1, survivorEnd = 0;
        for (int64_t j = 0; j < n; j++) {
            if (segments[j].threadIndex == survivorThread && segments[j].start < survivorEnd) {
                continue; //absorbed into the previous run, and freed
            }
            stPinchSegment *survivor = stPinchSegment_joinTrivialBoundaries(segments[j].segment);
            survivorThread = segments[j].threadIndex;
            survivorEnd = stPinchSegment_getStart(survivor) + stPinchSegment_getLength(survivor);
        }
        if (segments != stackSegments) {
            free(segments);
        }
    }
    stList_setDestructor(blocksToDelete, NULL); //the blocks are gone
}

/*
 * The block-end pass of stPinchThreadSet_joinTrivialBoundaries, restricted to the blocks at the ends
 * of threads. After a melt those are the only blocks whose boundaries can be trivial: a boundary is
 * trivial only where two blocks meet with no segment between them, deleting blocks only ever puts
 * block-less segments between blocks, and the graph was at a join fixpoint before the melt except for
 * the thread-end splits stCaf_ensureEndsAreDistinct made after the previous join. Each block is
 * visited once, when its head segment is reached, as the full pass would visit it.
 */
static int64_t joinTrivialBoundariesAtThreadEnds(stPinchThreadSet *threadSet) {
    int64_t joins = 0;
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *thread;
    while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
        stPinchSegment *first = stPinchThread_getFirst(thread);
        stPinchBlock *block = stPinchSegment_getBlock(first);
        if (block != NULL && stPinchBlock_getFirst(block) == first) {
            joins += stPinchBlock_joinTrivialBoundaries(block);
        }
        stPinchSegment *last = stPinchThread_getLast(thread); //read after the join above, which may have changed it
        block = stPinchSegment_getBlock(last);
        if (block != NULL && stPinchBlock_getFirst(block) == last) {
            joins += stPinchBlock_joinTrivialBoundaries(block);
        }
    }
    return joins;
}

static int64_t totalSegmentCount(stPinchThreadSet *threadSet) {
    int64_t count = 0;
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *thread;
    while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
        count += stPinchThread_getSegmentCount(thread);
    }
    return count;
}

int64_t stCaf_meltChains(Flower *flower, stPinchThreadSet *threadSet, int64_t minimumChainLength,
                bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds) {
    assert(minimumChainLength > 1);
    double t = stCaf_now();
    stCactusNode *startCactusNode;
    stPinchComponent *deadEndComponent;
    stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, INT64_MAX,
            0.0, breakChainsAtReverseTandems, maximumMedianSpacingBetweenLinkedEnds);
    double graphTime = stCaf_now() - t;
    t = stCaf_now();
    stList *blocksToDelete = stCaf_getBlocksInChainsLessThanGivenLength(cactusGraph, minimumChainLength);
    double scanTime = stCaf_now() - t;
    t = stCaf_now();
    int64_t blocksDestroyed = stList_length(blocksToDelete);

    st_logInfo("A melting round is destroying %" PRIi64 " blocks with an average degree "
           "of %lf from chains with length less than %" PRIi64 ". Total aligned bases"
           " lost: %" PRIu64 "\n",
           stList_length(blocksToDelete), stCaf_averageBlockDegree(blocksToDelete),
           minimumChainLength, stCaf_totalAlignedBases(blocksToDelete));

    stCaf_destructCactusGraph(cactusGraph, threadSet);
    destroyBlocksJoiningRuns(blocksToDelete);
    stList_destruct(blocksToDelete);
    double deleteTime = stCaf_now() - t;
    t = stCaf_now();
    int64_t joins = joinTrivialBoundariesAtThreadEnds(threadSet);
    stCaf_ensureEndsAreDistinct(threadSet);
    double joinTime = stCaf_now() - t;
    st_logInfo("caf-timing: melt minChain=%" PRIi64 " trim 0.000s filter 0.000s graph %.3fs scan %.3fs delete %.3fs join %.3fs destroyed %" PRIi64 " end-joins %" PRIi64 "\n",
               minimumChainLength, graphTime, scanTime, deleteTime, joinTime, blocksDestroyed, joins);

    if (getenv("CACTUS_CAF_CHECK_JOIN") != NULL) {
        // The full join must now find nothing to do
        int64_t segmentsBefore = totalSegmentCount(threadSet);
        int64_t changes = stPinchThreadSet_joinTrivialBoundaries(threadSet);
        stCaf_ensureEndsAreDistinct(threadSet);
        int64_t segmentsAfter = totalSegmentCount(threadSet);
        if (changes != 0 || segmentsAfter != segmentsBefore) {
            st_errAbort("CACTUS_CAF_CHECK_JOIN: the full join after the melt with minimum chain length %" PRIi64
                        " made %" PRIi64 " changes and took the segment count from %" PRIi64 " to %" PRIi64 "\n",
                        minimumChainLength, changes, segmentsBefore, segmentsAfter);
        }
        st_logInfo("caf-timing: check-join ok\n");
    }
    return blocksDestroyed;
}

static bool isTelomere(stPinchEnd *end, stPinchComponent *deadEndComponent) {
    stPinchSegment *segment = stPinchBlock_getFirst(end->block);
    bool atEndOfThread = stPinchThread_getFirst(stPinchSegment_getThread(segment)) == segment || stPinchThread_getLast(stPinchSegment_getThread(segment)) == segment;
    bool inDeadEndComponent = stPinchEnd_getComponent(end) == deadEndComponent;
    return atEndOfThread || inDeadEndComponent;
}

static bool endSetContainsTelomere(stSet *endSet, stPinchComponent *deadEndComponent) {
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

static bool chainConnectsToTelomere(stCactusEdgeEnd *chainEnd, stPinchComponent *deadEndComponent) {
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
static bool chainIsRecoverable(stCactusEdgeEnd *chainEnd, stPinchComponent *deadEndComponent) {
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
static stCactusEdgeEnd *getChainEndFromSingletonSet(stSet *ends) {
    assert(stSet_size(ends) == 1);
    stSetIterator *it = stSet_getIterator(ends);
    stPinchEnd *connectedPinchEnd = stSet_getNext(it);
    stSet_destructIterator(it);

    stCactusEdgeEnd *chainEnd = stPinchEnd_getData(connectedPinchEnd);
    assert(chainEnd != NULL);
    if (!stCactusEdgeEnd_getLinkOrientation(chainEnd)) {
        chainEnd = stCactusEdgeEnd_getLink(chainEnd);
    }
    return chainEnd;
}

// Mark down which chain(s) this (recoverable) chain is recoverable given.
static void markRecoverableAdjacencies(stCactusEdgeEnd *recoverableChainEnd,
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
        stCactusEdgeEnd *connectedChainEnd = getChainEndFromSingletonSet(connectedEnds1);
        stList_append(recoverableAdjacencies, connectedChainEnd);
    }

    if (stSet_size(connectedEnds2) == 1) {
        stCactusEdgeEnd *connectedChainEnd = getChainEndFromSingletonSet(connectedEnds2);
        stList_append(recoverableAdjacencies, connectedChainEnd);
    }

    stHash_insert(chainToRecoverableAdjacencies, recoverableChainEnd, recoverableAdjacencies);

    stSet_destruct(connectedEnds1);
    stSet_destruct(connectedEnds2);
}

/*
 * Record on each pinch end (in its data slot) the canonical chain end for its chain.
 */

static void setPinchEndToChainEnd(stCactusGraph *cactusGraph, stPinchThreadSet *threadSet) {
    stPinchThreadSet_clearEndData(threadSet); //the slots held the cactus nodes while the graph was built
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
                    stPinchEnd_setData(pinchEnd, chainEnd);
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
}

// For a given cactus node, recurse through all nodes below it and
// find recoverable chains below them. Then find recoverable chains
// below the current node given its parent chain.
static void getRecoverableChains_R(stCactusNode *cactusNode, stCactusEdgeEnd *parentChain, stPinchComponent *deadEndComponent, Flower *flower, bool (*recoverabilityFilter)(stCactusEdgeEnd *, Flower *), stSet *recoverableChains, stList *telomereAdjacentChains, stHash *chainToRecoverableAdjacencies) {
    while(1) {
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
                                       recoverableChains,
                                       telomereAdjacentChains,
                                       chainToRecoverableAdjacencies);
            }
        }

        cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
        while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt)) != NULL) {
            if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) {
                if ((recoverabilityFilter == NULL || recoverabilityFilter(cactusEdgeEnd, flower)) &&
                    chainIsRecoverable(cactusEdgeEnd, deadEndComponent)) {
                    stSet_insert(recoverableChains, cactusEdgeEnd);
                    markRecoverableAdjacencies(cactusEdgeEnd, chainToRecoverableAdjacencies);
                    if (chainConnectsToTelomere(cactusEdgeEnd, deadEndComponent)) {
                        stList_append(telomereAdjacentChains, cactusEdgeEnd);
                    }
                }
            }
        }

        if (parentChain == NULL) { // If we're not traversing a chain, we're done
            break;
        }
        // Visit the next node on this chain (unless it's where we started).
        stCactusEdgeEnd *nextEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd_getLink(parentChain));
        if (stCactusEdgeEnd_isChainEnd(nextEdgeEnd)) { // If at the end, stop
            break;
        }
        // Now get the next link in the chain
        assert(cactusNode != stCactusEdgeEnd_getNode(nextEdgeEnd));
        cactusNode = stCactusEdgeEnd_getNode(nextEdgeEnd);
        parentChain = nextEdgeEnd;
    }
}

static stList *getRecoverableChains(stCactusGraph *cactusGraph, stCactusNode *startCactusNode, stPinchComponent *deadEndComponent, Flower *flower, bool (*recoverabilityFilter)(stCactusEdgeEnd *, Flower *), stPinchThreadSet *threadSet) {
    setPinchEndToChainEnd(cactusGraph, threadSet);

    stSet *recoverableChainSet = stSet_construct();
    stList *telomereAdjacentChains = stList_construct();
    stHash *chainToRecoverableAdjacencies = stHash_construct2(NULL, (void (*)(void *)) stList_destruct);
    getRecoverableChains_R(startCactusNode, NULL, deadEndComponent, flower, recoverabilityFilter, recoverableChainSet, telomereAdjacentChains, chainToRecoverableAdjacencies);

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

int64_t stCaf_meltRecoverableChains(Flower *flower, stPinchThreadSet *threadSet, bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds, bool (*recoverabilityFilter)(stCactusEdgeEnd *, Flower *), int64_t maxNumIterations, int64_t maxRecoverableChainLength) {
    int64_t iteration = 0, totalBlocksDestroyed = 0;
    while (maxNumIterations-- > 0) {
        double t = stCaf_now();
        stCactusNode *startCactusNode;
        stPinchComponent *deadEndComponent;
        // FIXME: We shouldn't really have to rebuild the cactus graph
        // every time. Instead we should just be able to do multiple
        // iterations over the same graph, keeping track of the chains
        // whose underlying blocks we've already deleted.
        stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, 0,
                                                                      0.0, breakChainsAtReverseTandems, maximumMedianSpacingBetweenLinkedEnds);
        double graphTime = stCaf_now() - t;
        t = stCaf_now();
        double endSetTime = 0.0; //dead end membership is now read off the block end records

        stList *recoverableChains = getRecoverableChains(cactusGraph, startCactusNode, deadEndComponent, flower, recoverabilityFilter, threadSet);

        stList *blocksToDelete = stList_construct3(0, (void(*)(void *)) stPinchBlock_destruct);
        for (int64_t i = 0; i < stList_length(recoverableChains); i++) {
            stCactusEdgeEnd *chainEnd = stList_get(recoverableChains, i);
            if (getChainLength(chainEnd) <= maxRecoverableChainLength) {
                addChainBlocksToBlocksToDelete(chainEnd, blocksToDelete);
            }
        }
        double findTime = stCaf_now() - t;
        t = stCaf_now();
        int64_t numRecoverableBlocks = stList_length(blocksToDelete);
        st_logInfo("Destroying %" PRIi64 " recoverable blocks\n", numRecoverableBlocks);
        st_logInfo("The blocks covered %" PRIi64 " columns for a total of %" PRIi64 " aligned bases\n", numColumns(blocksToDelete), totalAlignedBases(blocksToDelete));
        stList_destruct(recoverableChains);
        stList_destruct(blocksToDelete);

        stCaf_destructCactusGraph(cactusGraph, threadSet);
        st_logInfo("caf-timing: recoverable iter %" PRIi64 " graph %.3fs endset %.3fs find %.3fs delete %.3fs destroyed %" PRIi64 "\n",
                   iteration++, graphTime, endSetTime, findTime, stCaf_now() - t, numRecoverableBlocks);
        totalBlocksDestroyed += numRecoverableBlocks;

        if (numRecoverableBlocks == 0) {
            // We didn't delete anything this round; we can safely
            // stop since we haven't changed the graph at all.
            break;
        }
    }
    return totalBlocksDestroyed;
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
