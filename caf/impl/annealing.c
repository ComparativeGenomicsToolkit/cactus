#include "sonLib.h"
#include "cactus.h"
#include "pairwiseAlignment.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stCactusGraphs.h"
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

static int uint64_cmp(const uint64_t *x, const uint64_t *y) {
    if (*x < *y) {
        return -1;
    } else if (*x == *y) {
        return 0;
    } else {
        return 1;
    }
}

// Print a set of statistics (avg, median, max, min) for degree and
// support percentage in the pinch graph.
void printThreadSetStatistics(stPinchThreadSet *threadSet, Flower *flower, FILE *f)
{
    // Naively finds the max, median, and min by sorting: the lists
    // will have "only" millions of elements, so they should fit
    // comfortably into tens of MB of memory.

    uint64_t numBlocks = stPinchThreadSet_getTotalBlockNumber(threadSet);

    uint64_t *blockDegrees = malloc(numBlocks * sizeof(uint64_t));
    double totalDegree = 0.0;

    uint64_t totalAlignedBases = 0;

    stPinchThreadSetBlockIt it = stPinchThreadSet_getBlockIt(threadSet);
    uint64_t i = 0;
    stPinchBlock *block;
    while ((block = stPinchThreadSetBlockIt_getNext(&it)) != NULL) {
        blockDegrees[i] = stPinchBlock_getDegree(block);
        totalDegree += stPinchBlock_getDegree(block);

        totalAlignedBases += stPinchBlock_getLength(block) * stPinchBlock_getDegree(block);

        i++;
    }

    printf("There were %" PRIu64 " blocks in the sequence graph, representing %" PRIi64
           " total aligned bases\n", numBlocks, totalAlignedBases);

    qsort(blockDegrees, numBlocks, sizeof(uint64_t),
          (int (*)(const void *, const void *)) uint64_cmp);
    printf("Block degree stats: min %" PRIu64 ", avg %lf, median %" PRIu64 ", max %" PRIu64 "\n",
           blockDegrees[0], totalDegree/numBlocks, blockDegrees[(numBlocks - 1) / 2],
           blockDegrees[numBlocks - 1]);
    free(blockDegrees);
}

void stCaf_annealPreventingSmallChains(Flower *flower, stPinchThreadSet *threadSet,
                                       const char *alignmentsFile,
                                       stList *alignmentsList,
                                       int64_t alignmentTrim,
                                       bool (*filterFn)(stPinchSegment *, stPinchSegment *),
                                       stList *minimumChainLengths,
                                       int64_t numAlignmentsPerBatch,
                                       double maxRedundantFraction) {
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
    stList *pinches = stList_construct3(0, (void (*)(void *)) stPinch_destruct);
    while (alignment != NULL) {
        numAlignments++;

        // Get the set of pinches from the next gapped alignment.
        // TODO: make this prettier
        stList *pairwiseAlignments = stList_construct();
        stList_append(pairwiseAlignments, alignment);
        stPinchIterator *pinchIterator = stPinchIterator_constructFromList(pairwiseAlignments);
        stPinchIterator_setTrim(pinchIterator, alignmentTrim);
        stPinch *pinch;
        while ((pinch = stPinchIterator_getNext(pinchIterator)) != NULL) {
            stList_append(pinches, stPinch_construct(pinch->name1, pinch->name2,
                                                     pinch->start1, pinch->start2,
                                                     pinch->length, pinch->strand));
        }
        stPinchIterator_destruct(pinchIterator);
        stList_destruct(pairwiseAlignments);
        numPinches += stList_length(pinches);
        if (alignmentsFile != NULL) {
            alignment = cigarRead(alignments);
        } else {
            alignment = stList_getNext(listIt);
        }

        // First, we anneal (adding alignments to the pinch graph).
        stPinchThread *thread1 = NULL;
        stPinchThread *thread2 = NULL;

        // Count the redundant pairs in the batch, so we can enforce the % redundant threshold.
        int64_t redundantPairs = 0;
        int64_t totalPairs = 0;
        for (int64_t i = 0; i < stList_length(pinches); i++) {
            pinch = stList_get(pinches, i);
            thread1 = stPinchThreadSet_getThread(threadSet, pinch->name1);
            thread2 = stPinchThreadSet_getThread(threadSet, pinch->name2);

            // Record the number of pairs in this pinch that are already aligned (redundant pairs).
            for (int64_t j = 0; j < pinch->length; j++) {
                int64_t pos1 = pinch->start1 + j;
                int64_t pos2 = pinch->strand ? pinch->start2 + j : pinch->start2 + pinch->length - 1 - j;
                stPinchSegment *segment1 = stPinchThread_getSegment(thread1, pos1);
                stPinchSegment *segment2 = stPinchThread_getSegment(thread2, pos2);
                if (stPinchSegment_getBlock(segment1) != NULL && stPinchSegment_getBlock(segment1) == stPinchSegment_getBlock(segment2)) {
                    int64_t offset1 = pos1 - stPinchSegment_getStart(segment1);
                    int64_t offset2 = pos2 - stPinchSegment_getStart(segment2);
                    if ((pinch->strand && offset1 == offset2) || (!pinch->strand && stPinchBlock_getLength(stPinchSegment_getBlock(segment1)) - 1 - offset2 == offset1)) {
                        redundantPairs++;
                    }
                }
            }
            totalPairs += pinch->length;
        }

        double redundantFraction = ((double) redundantPairs) / totalPairs;
        if (redundantFraction > maxRedundantFraction) {
            // This batch must be rejected, since it doesn't pass the redundancy filter.
            stList_destruct(pinches);
            pinches = stList_construct3(0, (void (*)(void *)) stPinch_destruct);
            continue;
        }

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
        // If the batch is done, then start the melting for the batch.
        if (numAlignments % numAlignmentsPerBatch == 0) {
            // Melt, enforcing each minimum chain length in succession.
            printf("Online melting round after %" PRIi64 " alignments. Thread set statistics before melting:\n", numAlignments);
            printThreadSetStatistics(threadSet, flower, stdout);
            for (int64_t i = 0; i < stList_length(minimumChainLengths); i++) {
                int64_t minimumChainLength = stIntTuple_get(stList_get(minimumChainLengths, i), 0);
                stCaf_melt(flower, threadSet, NULL, 0, minimumChainLength, 0, INT64_MAX);
            }
            printf("Thread set statistics after online melting:\n");
            printThreadSetStatistics(threadSet, flower, stdout);
        }

        stList_destruct(pinches);
        pinches = stList_construct3(0, (void (*)(void *)) stPinch_destruct);
    }
    printf("Added %" PRIi64 " gapped alignments containing %" PRIi64 " pinches to the graph.\n", numAlignments, numPinches);
}
