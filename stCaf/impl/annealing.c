#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

///////////////////////////////////////////////////////////////////////////
// Ensure that flower ends are in blocks without other sequence.
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

///////////////////////////////////////////////////////////////////////////
// Basic annealing function
///////////////////////////////////////////////////////////////////////////

void stCaf_anneal(stPinchThreadSet *threadSet, stPinchIterator *pinchIterator) {
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
    stCaf_ensureEndsAreDistinct(threadSet);
}

///////////////////////////////////////////////////////////////////////////
// Annealing function that ignores homologies between bases not in the same adjacency component.
///////////////////////////////////////////////////////////////////////////

static int64_t getIntersectionLength(int64_t start1, int64_t start2, stPinchInterval *pinchInterval1,
        stPinchInterval *pinchInterval2) {
    int64_t length1 = pinchInterval1->length + pinchInterval1->start - start1;
    int64_t length2 = pinchInterval2->length + pinchInterval2->start - start2;
    assert(length1 > 0 && length2 > 0);
    return length1 > length2 ? length1 : length2;
}

static int64_t getIntersectionLengthReverse(int64_t start1, int64_t end2, stPinchInterval *pinchInterval1,
        stPinchInterval *pinchInterval2) {
    int64_t length1 = pinchInterval1->length + pinchInterval1->start - start1;
    int64_t length2 = end2 - pinchInterval2->start + 1;
    assert(length1 > 0 && length2 > 0);
    return length1 > length2 ? length1 : length2;
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

static void alignSameComponents(stPinch *pinch, stPinchThreadSet *threadSet, stSortedSet *adjacencyComponentIntervals) {
    stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
    stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);
    assert(thread1 != NULL && thread2 != NULL);
    stPinchInterval *pinchInterval1 = stPinchIntervals_getInterval(adjacencyComponentIntervals, pinch->name1,
            pinch->start1);
    assert(pinchInterval1 != NULL);
    int64_t offset = 0;
    if (pinch->strand) { //A but redundant, but fast.
        stPinchInterval *pinchInterval2 = stPinchIntervals_getInterval(adjacencyComponentIntervals, pinch->name2,
                pinch->start2);
        assert(pinchInterval2 != NULL);
        while (offset < pinch->length) {
            int64_t length = getIntersectionLength(pinch->start1 + offset, pinch->start2 + offset, pinchInterval1,
                    pinchInterval2);
            if (stPinchInterval_getLabel(pinchInterval1) == stPinchInterval_getLabel(pinchInterval2)) {
                stPinchThread_pinch(thread1, thread2, pinch->start1 + offset, pinch->start2 + offset, length, 1);
            }
            offset += length;
            pinchInterval1 = updatePinchInterval(pinch->start1 + offset, pinchInterval1, adjacencyComponentIntervals);
            pinchInterval2 = updatePinchInterval(pinch->start2 + offset, pinchInterval2, adjacencyComponentIntervals);
        }
    } else {
        stPinchInterval *pinchInterval2 = stPinchIntervals_getInterval(adjacencyComponentIntervals, pinch->name2,
                pinch->start2 + pinch->length - 1);
        assert(pinchInterval2 != NULL);
        int64_t end2 = pinch->start2 + pinch->length - 1;
        while (offset < pinch->length) {
            int64_t length = getIntersectionLengthReverse(pinch->start1 + offset, end2 - offset, pinchInterval1,
                    pinchInterval2);
            if (stPinchInterval_getLabel(pinchInterval1) == stPinchInterval_getLabel(pinchInterval2)) {
                stPinchThread_pinch(thread1, thread2, pinch->start1 + offset, end2 - offset - length + 1, length, 0);
            }
            offset += length;
            pinchInterval1 = updatePinchInterval(pinch->start1 + offset, pinchInterval1, adjacencyComponentIntervals);
            pinchInterval2 = updatePinchIntervalReverse(end2 - offset, pinchInterval2, adjacencyComponentIntervals);
        }
    }
}

void stCaf_annealBetweenAdjacencyComponents(stPinchThreadSet *threadSet, stPinchIterator *pinchIterator) {
    //Get adjacency components
    stHash *pinchEndsToAdjacencyComponents;
    stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents2(threadSet, &pinchEndsToAdjacencyComponents);
    stList_setDestructor(adjacencyComponents, NULL);
    stList_destruct(adjacencyComponents);
    //Get the adjacency component intervals
    stSortedSet *adjacencyComponentIntervals = stPinchThreadSet_getLabelIntervals(threadSet,
            pinchEndsToAdjacencyComponents);
    //Now do the actual alignments.
    stPinchIterator_reset(pinchIterator);
    stPinch *pinch;
    while ((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
        alignSameComponents(pinch, threadSet, adjacencyComponentIntervals);
    }
    stPinchThreadSet_joinTrivialBoundaries(threadSet);
    stCaf_ensureEndsAreDistinct(threadSet);
    stSortedSet_destruct(adjacencyComponentIntervals);
}
