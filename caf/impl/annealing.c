#include "sonLib.h"
#include "cactus.h"
#include "pairwiseAlignment.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stCactusGraphs.h"
#include "stOnlineCactus.h"
#include "stCaf.h"

///////////////////////////////////////////////////////////////////////////
// Code to safely join all the trivial boundaries in the pinch graph, while
// respecting end blocks.
///////////////////////////////////////////////////////////////////////////

static void stCaf_ensureEndsAreDistinct(stPinchThreadSet *threadSet) {
    /*
     * Ensures the blocks at the ends of threads are distinct.
     */
    stPinchThread *thread;
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
        stPinchThread_split(thread, stPinchThread_getStart(thread));
        assert(stPinchThread_getLength(thread) > 1);
        stPinchThread_split(thread, stPinchThread_getStart(thread) + stPinchThread_getLength(thread) - 2);
    }
}

void stCaf_joinTrivialBoundaries(stPinchThreadSet *threadSet) {
    stPinchThreadSet_joinTrivialBoundaries(threadSet);
    stCaf_ensureEndsAreDistinct(threadSet);
}

///////////////////////////////////////////////////////////////////////////
// Basic annealing function
///////////////////////////////////////////////////////////////////////////

void stCaf_anneal2(stPinchThreadSet *threadSet, stPinch *(*pinchIterator)(void *), void *extraArg) {
    stPinch *pinch;
    while ((pinch = pinchIterator(extraArg)) != NULL) {
        stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
        stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
        assert(thread1 != NULL && thread2 != NULL);
        stPinchThread_pinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand);
    }
}

static void stCaf_annealWithFilter2(stPinchThreadSet *threadSet, stPinch *(*pinchIterator)(void *), void *extraArg, bool (*filterFn)(stPinchSegment *, stPinchSegment *)) {
    stPinch *pinch;
    while ((pinch = pinchIterator(extraArg)) != NULL) {
        stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
        stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
        assert(thread1 != NULL && thread2 != NULL);
        stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand, filterFn);
    }
}

void stCaf_anneal(stPinchThreadSet *threadSet, stPinchIterator *pinchIterator, bool (*filterFn)(stPinchSegment *, stPinchSegment *)) {
    stPinchIterator_reset(pinchIterator);
    if(filterFn != NULL) {
        stCaf_annealWithFilter2(threadSet, (stPinch *(*)(void *)) stPinchIterator_getNext, pinchIterator, filterFn);
    }
    else {
        stCaf_anneal2(threadSet, (stPinch *(*)(void *)) stPinchIterator_getNext, pinchIterator);
    }
    stCaf_joinTrivialBoundaries(threadSet);
}

double averageBlockDegree(stList *blocks) {
    if (stList_length(blocks) == 0) {
        return 0.0;
    }
    uint64_t total = 0;
    for (int64_t i = 0; i < stList_length(blocks); i++) {
        total += stPinchBlock_getDegree(stList_get(blocks, i));
    }
    return ((double) total)/stList_length(blocks);
}

void stCaf_annealPreventingSmallChains(Flower *flower, stPinchThreadSet *threadSet,
                                       stOnlineCactus *cactus,
                                       const char *alignmentsFile,
                                       stList *alignmentsList,
                                       int64_t alignmentTrim,
                                       bool (*filterFn)(stPinchSegment *, stPinchSegment *),
                                       stList *minimumChainLengths,
                                       stCaf_meltingMethod meltingMethod) {
    stListIterator *listIt = NULL;
    FILE *alignments = NULL;
    struct PairwiseAlignment *alignment;
    if (alignmentsFile != NULL) {
        alignments = fopen(alignmentsFile, "r");
        alignment = cigarRead(alignments);
    } else {
        listIt = stList_getIterator(alignmentsList);
        alignment = stList_getNext(listIt);
    }
    size_t numAlignments = 0;
    size_t numPinches = 0;
    while (alignment != NULL) {
        numAlignments++;
        // TODO: make this prettier
        stList *pairwiseAlignments = stList_construct();
        stList_append(pairwiseAlignments, alignment);
        stPinchIterator *pinchIterator = stPinchIterator_constructFromList(pairwiseAlignments);
        stPinchIterator_setTrim(pinchIterator, alignmentTrim);
        stList *pinches = stList_construct3(0, (void (*)(void *)) stPinch_destruct);
        stPinch *pinch;
        while ((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
            stList_append(pinches, stPinch_construct(pinch->name1, pinch->name2,
                                                     pinch->start1, pinch->start2,
                                                     pinch->length, pinch->strand));
        }
        stPinchIterator_destruct(pinchIterator);
        stList_destruct(pairwiseAlignments);
        numPinches += stList_length(pinches);
        stPinchUndo *undo = stPinchThreadSet_prepareGappedUndo(threadSet, pinches);
        stPinchThread *thread1 = NULL;
        stPinchThread *thread2 = NULL;
        for (int64_t i = 0; i < stList_length(pinches); i++) {
            pinch = stList_get(pinches, i);
            thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
            thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);

            // Apply the pinch to the graph.
            assert(thread1 != NULL && thread2 != NULL);
            if (filterFn == NULL) {
                stPinchThread_pinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand);
            } else {
                stPinchThread_filterPinch(thread1, thread2, pinch->start1, pinch->start2, pinch->length, pinch->strand, filterFn);
            }
        }

        for (int64_t i = 0; i < stList_length(minimumChainLengths); i++) {
            int64_t minimumChainLength = stIntTuple_get(stList_get(minimumChainLengths, i), 0);
            if (meltingMethod == PRESERVE_NON_UNDOABLE_CHAINS) {
                stCaf_undoChainsSmallerThanThis_preserveNonUndoableChains(cactus, threadSet, pinches, undo, minimumChainLength);
            } else if (meltingMethod == REMOVE_NON_UNDOABLE_CHAINS) {
                stCaf_undoChainsSmallerThanThis_removeNonUndoableChains(cactus, threadSet, pinches, undo, minimumChainLength);
            } else if (meltingMethod == ONLY_UNDO) {
                stCaf_undoChainsSmallerThanThis_onlyUndo(cactus, threadSet, pinches, undo, minimumChainLength);
            } else if (meltingMethod == ONLY_REMOVE) {
                stCaf_undoChainsSmallerThanThis_onlyRemove(cactus, threadSet, pinches, undo, minimumChainLength);
            }
        }
        stPinchUndo_destruct(undo);
        stList_destruct(pinches);
        if (alignmentsFile != NULL) {
            alignment = cigarRead(alignments);
        } else {
            alignment = stList_getNext(listIt);
        }
    }
    printf("Added %" PRIi64 " gapped alignments containing %" PRIi64 " pinches to the graph.\n", numAlignments, numPinches);
}

///////////////////////////////////////////////////////////////////////////
// Annealing function that ignores homologies between bases not in the same adjacency component.
///////////////////////////////////////////////////////////////////////////

static int64_t getIntersectionLength(int64_t start1, int64_t start2, stPinchInterval *pinchInterval1,
        stPinchInterval *pinchInterval2) {
    int64_t length1 = pinchInterval1->length + pinchInterval1->start - start1;
    int64_t length2 = pinchInterval2->length + pinchInterval2->start - start2;
    assert(length1 > 0 && length2 > 0);
    return length1 > length2 ? length2 : length1;
}

static int64_t getIntersectionLengthReverse(int64_t start1, int64_t end2, stPinchInterval *pinchInterval1,
        stPinchInterval *pinchInterval2) {
    int64_t length1 = pinchInterval1->length + pinchInterval1->start - start1;
    int64_t length2 = end2 - pinchInterval2->start + 1;
    assert(length1 > 0 && length2 > 0);
    return length1 > length2 ? length2 : length1;
}

static stPinchInterval *updatePinchInterval(int64_t start, stPinchInterval *pinchInterval,
        stSortedSet *adjacencyComponentIntervals) {
    return start < pinchInterval->start + pinchInterval->length ? pinchInterval : stPinchIntervals_getInterval(
            adjacencyComponentIntervals, pinchInterval->name, start);
}

static stPinchInterval *updatePinchIntervalReverse(int64_t end, stPinchInterval *pinchInterval,
        stSortedSet *adjacencyComponentIntervals) {
    return end >= pinchInterval->start ? pinchInterval : stPinchIntervals_getInterval(adjacencyComponentIntervals,
            pinchInterval->name, end);
}

static int64_t min(int64_t i, int64_t j) {
    return i < j ? i : j;
}

static void alignSameComponents(stPinch *pinch, stPinchThreadSet *threadSet, stSortedSet *adjacencyComponentIntervals, bool (*filterFn)(stPinchSegment *, stPinchSegment *)) {
    stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
    stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
    assert(thread1 != NULL && thread2 != NULL);
    stPinchInterval *pinchInterval1 = stPinchIntervals_getInterval(adjacencyComponentIntervals, pinch->name1,
            pinch->start1);
    int64_t offset = 0;
    if (pinch->strand) { //A bit redundant code wise, but fast.
        stPinchInterval *pinchInterval2 = stPinchIntervals_getInterval(adjacencyComponentIntervals, pinch->name2,
                pinch->start2);
        while (offset < pinch->length) {
            assert(pinchInterval1 != NULL && pinchInterval2 != NULL);
            int64_t length = min(getIntersectionLength(pinch->start1 + offset, pinch->start2 + offset, pinchInterval1,
                    pinchInterval2), pinch->length - offset);
            if (stPinchInterval_getLabel(pinchInterval1) == stPinchInterval_getLabel(pinchInterval2)) {
                if(filterFn != NULL) {
                    stPinchThread_filterPinch(thread1, thread2, pinch->start1 + offset, pinch->start2 + offset, length, 1, filterFn);
                }
                else {
                    stPinchThread_pinch(thread1, thread2, pinch->start1 + offset, pinch->start2 + offset, length, 1);
                }
            }
            offset += length;
            pinchInterval1 = updatePinchInterval(pinch->start1 + offset, pinchInterval1, adjacencyComponentIntervals);
            pinchInterval2 = updatePinchInterval(pinch->start2 + offset, pinchInterval2, adjacencyComponentIntervals);
        }
    } else {
        int64_t end2 = pinch->start2 + pinch->length - 1;
        stPinchInterval *pinchInterval2 = stPinchIntervals_getInterval(adjacencyComponentIntervals, pinch->name2, end2);
        while (offset < pinch->length) {
            assert(pinchInterval1 != NULL && pinchInterval2 != NULL);
            int64_t length = min(getIntersectionLengthReverse(pinch->start1 + offset, end2 - offset, pinchInterval1,
                    pinchInterval2), pinch->length - offset);
            if (stPinchInterval_getLabel(pinchInterval1) == stPinchInterval_getLabel(pinchInterval2)) {
                if(filterFn != NULL) {
                    stPinchThread_filterPinch(thread1, thread2, pinch->start1 + offset, end2 - offset - length + 1, length, 0, filterFn);
                }
                else {
                    stPinchThread_pinch(thread1, thread2, pinch->start1 + offset, end2 - offset - length + 1, length, 0);
                }
            }
            offset += length;
            pinchInterval1 = updatePinchInterval(pinch->start1 + offset, pinchInterval1, adjacencyComponentIntervals);
            pinchInterval2 = updatePinchIntervalReverse(end2 - offset, pinchInterval2, adjacencyComponentIntervals);
        }
    }
}

static stSortedSet *getAdjacencyComponentIntervals(stPinchThreadSet *threadSet, stList **adjacencyComponents) {
    stHash *pinchEndsToAdjacencyComponents;
    *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents2(threadSet, &pinchEndsToAdjacencyComponents);
    stSortedSet *adjacencyComponentIntervals = stPinchThreadSet_getLabelIntervals(threadSet,
            pinchEndsToAdjacencyComponents);
    stHash_destruct(pinchEndsToAdjacencyComponents);
    return adjacencyComponentIntervals;
}

void stCaf_annealBetweenAdjacencyComponents2(stPinchThreadSet *threadSet, stPinch *(*pinchIterator)(void *),
        void *extraArg, bool (*filterFn)(stPinchSegment *, stPinchSegment *)) {
    //Get the adjacency component intervals
    stList *adjacencyComponents;
    stSortedSet *adjacencyComponentIntervals = getAdjacencyComponentIntervals(threadSet, &adjacencyComponents);
    //Now do the actual alignments.
    stPinch *pinch;
    while ((pinch = pinchIterator(extraArg)) != NULL) {
        alignSameComponents(pinch, threadSet, adjacencyComponentIntervals, filterFn);
    }
    stSortedSet_destruct(adjacencyComponentIntervals);
    stList_destruct(adjacencyComponents);
}

void stCaf_annealBetweenAdjacencyComponents(stPinchThreadSet *threadSet, stPinchIterator *pinchIterator, bool (*filterFn)(stPinchSegment *, stPinchSegment *)) {
    stPinchIterator_reset(pinchIterator);
    stCaf_annealBetweenAdjacencyComponents2(threadSet, (stPinch *(*)(void *)) stPinchIterator_getNext, pinchIterator, filterFn);
    stCaf_joinTrivialBoundaries(threadSet);
}
