/*
 * stPinchGraphs.c
 *
 *  Created on: 11 Apr 2012
 *      Author: benedictpaten
 */

//Basic data structures

#include <stdlib.h>
#include "sonLib.h"
#include "stPinchGraphs.h"

const char *ST_PINCH_GRAPH_EXCEPTION_ID = "ST_PINCH_GRAPH_EXCEPTION_ID";

struct _stThreadSet {
    stList *threads;
    stHash *threadsHash;
};

struct _stThread {
    int64_t name;
    int64_t start;
    int64_t length;
    stSortedSet *segments;
};

struct _stSegment {
    stThread *thread;
    int64_t start;
    stSegment *pSegment;
    stSegment *nSegment;
    stBlock *block;
    bool blockOrientation;
    stSegment *nBlockSegment;
};

struct _stBlock {
    uint32_t degree;
    stSegment *headSegment;
    stSegment *tailSegment;
};

//Blocks

static void connectBlockToSegment(stSegment *segment, bool orientation,
        stBlock *block, stSegment *nBlockSegment) {
    segment->block = block;
    segment->blockOrientation = orientation;
    segment->nBlockSegment = nBlockSegment;
}

stBlock *stBlock_construct2(stSegment *segment1) {
    stBlock *block = st_malloc(sizeof(stBlock));
    block->headSegment = segment1;
    block->tailSegment = segment1;
    connectBlockToSegment(segment1, 1, block, NULL);
    block->degree = 1;
    return block;
}

stBlock *stBlock_construct(stSegment *segment1, bool orientation1,
        stSegment *segment2, bool orientation2) {
    if (stSegment_getLength(segment1) != stSegment_getLength(segment2)) {
        stThrowNew(
                ST_COMPRESSION_EXCEPTION_ID,
                "Two segments that are being pinched have different lengths: %i %i\n",
                stSegment_getLength(segment1), stSegment_getLength(segment2));
    }
    stBlock *block = st_malloc(sizeof(stBlock));
    block->headSegment = segment1;
    block->tailSegment = segment2;
    connectBlockToSegment(segment1, orientation1, block, segment2);
    connectBlockToSegment(segment2, orientation2, block, NULL);
    block->degree = 2;
    return block;
}

void stBlock_destruct(stBlock *block) {
    stBlockIt blockIt = stBlock_getSegmentIterator(block);
    stSegment *segment = stBlockIt_getNext(&blockIt);
    while (segment != NULL) {
        stSegment *nSegment = stBlockIt_getNext(&blockIt);
        connectBlockToSegment(segment, 0, NULL, NULL);
        segment = nSegment;
    }
    free(block);
}

stBlock *stBlock_pinch(stBlock *block1, stBlock *block2, bool orientation) {
    if (block1 == block2) {
        return block1; //Already joined
    }
    if (stBlock_getDegree(block1) > stBlock_getDegree(block2)) { //Avoid merging large blocks into small blocks
        return stBlock_pinch(block2, block1, orientation);
    }
    if (stBlock_getLength(block1) != stBlock_getLength(block2)) {
        stThrowNew(
                ST_PINCH_GRAPH_EXCEPTION_ID,
                "Two segments that are being pinched have different lengths: %i %i\n",
                stBlock_getLength(block1), stBlock_getLength(block2));
    }
    stBlockIt blockIt = stBlock_getSegmentIterator(block2);
    stSegment *segment = stBlockIt_getNext(&blockIt);
    while (segment != NULL) {
        stSegment *nSegment = stBlockIt_getNext(&blockIt);
        bool segmentOrientation = stSegment_getBlockOrientation(segment);
        stBlock_pinch2(
                block1,
                segment,
                (segmentOrientation && orientation) || (!segmentOrientation
                        && !orientation));
        segment = nSegment;
    }
    free(block2);
    return block1;
}

stBlock *stBlock_pinch2(stBlock *block, stSegment *segment, bool orientation) {
    assert(block->tailSegment != NULL);
    assert(block->tailSegment->nBlockSegment == NULL);
    block->tailSegment->nBlockSegment = segment;
    connectBlockToSegment(segment, orientation, block, NULL);
    block->tailSegment = segment;
    block->degree++;
    return block;
}

stBlockIt stBlock_getSegmentIterator(stBlock *block) {
    stBlockIt blockIt;
    blockIt.segment = block->headSegment;
    return blockIt;
}

stSegment *stBlockIt_getNext(stBlockIt *blockIt) {
    stSegment *segment = blockIt->segment;
    if (segment != NULL) {
        blockIt->segment = segment->nBlockSegment;
    }
    return segment;
}

uint32_t stBlock_getDegree(stBlock *block) {
    return block->degree;
}

stSegment *stBlock_getFirst(stBlock *block) {
    assert(block->headSegment != NULL);
    return block->headSegment;
}

int64_t stBlock_getLength(stBlock *block) {
    return stSegment_getLength(stBlock_getFirst(block));
}

//Segments

int64_t stSegment_getStart(stSegment *segment) {
    return segment->start;
}

int64_t stSegment_getLength(stSegment *segment) {
    assert(segment->nSegment != NULL);
    return segment->nSegment->start - segment->start;
}

stBlock *stSegment_getBlock(stSegment *segment) {
    return segment->block;
}

bool stSegment_getBlockOrientation(stSegment *segment) {
    return segment->blockOrientation;
}

stSegment *stSegment_get5Prime(stSegment *segment) {
    return segment->pSegment;
}

stSegment *stSegment_get3Prime(stSegment *segment) {
    return segment->nSegment->nSegment != NULL ? segment->nSegment : NULL;
}

int64_t stSegment_getName(stSegment *segment) {
    return stThread_getName(segment->thread);
}

stThread *stSegment_getThread(stSegment *segment) {
    return segment->thread;
}

//Private segment functions

void stSegment_destruct(stSegment *segment) {
    if (stSegment_getBlock(segment) != NULL) {
        stBlock_destruct(stSegment_getBlock(segment));
    }
    free(segment);
}

int stSegment_compareFunction(const stSegment *segment1,
        const stSegment *segment2) {
    return segment1->start < segment2->start ? -1 : (segment1->start
            > segment2->start ? 1 : 0);
}

static stSegment *stSegment_construct(int64_t start, stThread *thread) {
    stSegment *segment = st_calloc(1, sizeof(stSegment));
    segment->start = start;
    segment->thread = thread;
    return segment;
}

static stSegment *stSegment_splitP(stSegment *segment, int64_t leftBlockLength) {
    stSegment *nSegment = segment->nSegment;
    assert(nSegment != NULL);
    stSegment *rightSegment = stSegment_construct(
            stSegment_getStart(segment) + leftBlockLength, segment->thread);
    segment->nSegment = rightSegment;
    rightSegment->pSegment = segment;
    rightSegment->nSegment = nSegment;
    nSegment->pSegment = rightSegment;
    stSortedSet_insert(segment->thread->segments, rightSegment);
    return rightSegment;
}

void stSegment_split(stSegment *segment, int64_t leftSideOfSplitPoint) {
    if (leftSideOfSplitPoint == stSegment_getStart(segment)
            + stSegment_getLength(segment) - 1) { //There is already a break
        return;
    }
    int64_t leftSegmentLength = leftSideOfSplitPoint - stSegment_getStart(
            segment) + 1;
    assert(leftSegmentLength > 0);
    stBlock *block;
    if ((block = stSegment_getBlock(segment)) != NULL) {
        int64_t rightSegmentLength = stBlock_getLength(block)
                - leftSegmentLength;
        assert(rightSegmentLength > 0);
        if (!stSegment_getBlockOrientation(segment)) {
            int64_t i = rightSegmentLength;
            rightSegmentLength = leftSegmentLength;
            leftSegmentLength = i;
        }
        stBlockIt blockIt = stBlock_getSegmentIterator(block);
        segment = stBlockIt_getNext(&blockIt);
        assert(segment != NULL);
        stSegment *pSegment = NULL;
        stBlock *block2;
        if (stSegment_getBlockOrientation(segment)) {
            stSegment *segment2 = stSegment_splitP(segment, leftSegmentLength);
            block2 = stBlock_construct2(segment2);
            pSegment = segment;
        } else {
            stSegment *segment2 = stSegment_splitP(segment, rightSegmentLength);
            block->headSegment = segment2;
            connectBlockToSegment(segment2, 0, block, segment->nBlockSegment);
            if (segment2->nBlockSegment == NULL) {
                block->tailSegment = segment2;
            }
            block2 = stBlock_construct2(segment);
            segment->blockOrientation = 0; //This gets sets positive by default.
            pSegment = segment2;
        }
        while ((segment = stBlockIt_getNext(&blockIt)) != NULL) {
            if (stSegment_getBlockOrientation(segment)) {
                stSegment *segment2 = stSegment_splitP(segment,
                        leftSegmentLength);
                stBlock_pinch2(block2, segment2, 1);
                pSegment = segment;
            } else {
                stSegment *segment2 = stSegment_splitP(segment,
                        rightSegmentLength);
                pSegment->nBlockSegment = segment2;
                connectBlockToSegment(segment2, 0, block,
                        segment->nBlockSegment);
                if (segment2->nBlockSegment == NULL) {
                    block->tailSegment = segment2;
                }
                stBlock_pinch2(block2, segment, 0);
                pSegment = segment2;
            }
        }
    } else {
        stSegment_splitP(segment, leftSegmentLength);
    }
}

//Thread

int64_t stThread_getName(stThread *thread) {
    return thread->name;
}

int64_t stThread_getStart(stThread *thread) {
    return thread->start;
}

int64_t stThread_getLength(stThread *thread) {
    return thread->length;
}

stSegment *stThread_getSegment(stThread *thread, int64_t coordinate) {
    static stSegment segment;
    segment.start = coordinate;
    stSegment *segment2 = stSortedSet_searchLessThanOrEqual(thread->segments,
            &segment);
    if (segment2 == NULL) {
        return NULL;
    }
    assert(stSegment_getStart(segment2) <= coordinate);
    if (stSegment_getStart(segment2) + stSegment_getLength(segment2)
            <= coordinate) {
        return NULL;
    }
    return segment2;
}

stSegment *stThread_getFirst(stThread *thread) {
    return stSortedSet_getFirst(thread->segments);
}

void stThread_split(stThread *thread, int64_t leftSideOfSplitPoint) {
    stSegment *segment = stThread_getSegment(thread, leftSideOfSplitPoint);
    if (segment == NULL) {
        return;
    }
    stSegment_split(segment, leftSideOfSplitPoint);
}

void stThread_joinTrivialBoundaries(stThread *thread) {
    stSegment *segment = stThread_getFirst(thread);
    do {
        if (stSegment_getBlock(segment) == NULL) {
            while (1) {
                stSegment *nSegment = stSegment_get3Prime(segment);
                if (nSegment != NULL) {
                    stBlock *nBlock = stSegment_getBlock(nSegment);
                    if (nBlock == NULL) {
                        //Trivial join
                        segment->nSegment = nSegment->nSegment;
                        assert(nSegment->nSegment != NULL);
                        nSegment->nSegment->pSegment = segment;
                        stSortedSet_remove(thread->segments, nSegment);
                        stSegment_destruct(nSegment);
                        continue;
                    }
                }
                break;
            }
        }
    } while ((segment = stSegment_get3Prime(segment)) != NULL);
}

stSegment *stThread_pinchP(stThread *thread, int64_t start) {
    stSegment *segment1 = stThread_getSegment(thread, start);
    assert(segment1 != NULL);
    if (stSegment_getStart(segment1) != start) {
        stSegment_split(segment1, start - 1);
        segment1 = stSegment_get3Prime(segment1);
    }
    return segment1;
}

stSegment *stThread_pinchTrim(stSegment *segment, bool strand, int32_t length) {
    assert(length > 0);
    if (stSegment_getLength(segment) <= length) {
        return segment;
    }
    if (strand) {
        stSegment_split(segment, stSegment_getStart(segment) + length - 1);
        return segment;
    }
    stSegment_split(
            segment,
            stSegment_getStart(segment) + stSegment_getLength(segment) - 1
                    - length);
    return stSegment_get3Prime(segment);
}

void stThread_pinch(stThread *thread1, stThread *thread2, int64_t start1,
        int64_t start2, int64_t length, bool strand2) {
    if (length == 0) {
        return;
    }
    stSegment *segment1 = stThread_pinchP(thread1, start1), *segment2;
    if (strand2) {
        segment2 = stThread_pinchP(thread2, start2);
    } else {
        segment2 = stThread_getSegment(thread2, start2 + length - 1);
        stSegment_split(segment2, start2 + length - 1);
    }
    while (length > 0) {
        if (segment1 == segment2) {
            if (strand2) {
                return; //This is a trivial alignment
            }
            if (stSegment_getLength(segment1) > 1) { //Split the block in two
                segment2 = stThread_pinchTrim(segment2, 0,
                        stSegment_getLength(segment1) / 2);
            }
        }
        do {
            int64_t i = stSegment_getLength(segment1);
            segment2 = stThread_pinchTrim(segment2, strand2,
                    length > i ? i : length);
            i = stSegment_getLength(segment2);
            stThread_pinchTrim(segment1, 1, length > i ? i : length);
        } while (stSegment_getLength(segment1) != stSegment_getLength(segment2));
        stBlock *block1, *block2;
        if ((block1 = stSegment_getBlock(segment1)) == NULL) {
            block1 = stBlock_construct2(segment1);
        }
        if ((block2 = stSegment_getBlock(segment2)) == NULL) {
            block2 = stBlock_construct2(segment2);
        }
        bool bO1 = stSegment_getBlockOrientation(segment1);
        bool bO2 = stSegment_getBlockOrientation(segment2);
        bool alignmentOrientation = ((bO1 == bO2) && strand2) || ((bO1 != bO2)
                && !strand2);
        if (block1 == block2) {
            if (stSegment_getLength(segment1) > 1 && !alignmentOrientation) {
                segment2 = stThread_pinchTrim(segment2, strand2,
                        stSegment_getLength(segment2) / 2);
                continue;
            }
        }
        block1 = stBlock_pinch(block1, block2, alignmentOrientation);
        length -= stSegment_getLength(segment1);
        segment1 = stSegment_get3Prime(segment1);
        segment2 = strand2 ? stSegment_get3Prime(segment2)
                : stSegment_get5Prime(segment2);
    }
}

//Private functions

static stThread *stThread_construct(int64_t name, int64_t start, int64_t length) {
    stThread *thread = st_malloc(sizeof(stThread));
    thread->name = name;
    thread->start = start;
    thread->length = length;
    thread->segments = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stSegment_compareFunction,
            (void(*)(void *)) stSegment_destruct);
    stSegment *segment = stSegment_construct(start, thread);
    stSegment *terminatorSegment = stSegment_construct(start + length, thread);
    segment->nSegment = terminatorSegment;
    terminatorSegment->pSegment = segment;
    stSortedSet_insert(thread->segments, segment);
    return thread;
}

static void stThread_destruct(stThread *thread) {
    stSortedSet_destruct(thread->segments);
    free(thread);
}

static uint32_t stThread_hashKey(const stThread *thread) {
    return thread->name;
}

static int stThread_equals(const stThread *thread1, const stThread *thread2) {
    return thread1->name == thread2->name;
}

//Thread set

stThreadSet *stThreadSet_construct() {
    stThreadSet *threadSet = st_malloc(sizeof(stThreadSet));
    threadSet->threads = stList_construct3(0,
            (void(*)(void *)) stThread_destruct);
    threadSet->threadsHash = stHash_construct3(
            (uint32_t(*)(const void *)) stThread_hashKey,
            (int(*)(const void *, const void *)) stThread_equals, NULL, NULL);
    return threadSet;
}

void stThreadSet_destruct(stThreadSet *threadSet) {
    stList_destruct(threadSet->threads);
    stHash_destruct(threadSet->threadsHash);
    free(threadSet);
}

stThread *stThreadSet_addThread(stThreadSet *threadSet, int64_t name,
        int64_t start, int64_t length) {
    stThread *thread = stThread_construct(name, start, length);
    assert(stThreadSet_getThread(threadSet, name) == NULL);
    stHash_insert(threadSet->threadsHash, thread, thread);
    stList_append(threadSet->threads, thread);
    return thread;
}

stThread *stThreadSet_getThread(stThreadSet *threadSet, int64_t name) {
    static stThread thread;
    thread.name = name;
    return stHash_search(threadSet->threadsHash, &thread);
}

int32_t stThreadSet_getSize(stThreadSet *threadSet) {
    return stList_length(threadSet->threads);
}

stThreadIt stThreadSet_getIterator(stThreadSet *threadSet) {
    stThreadIt threadIt;
    threadIt.threadSet = threadSet;
    threadIt.index = 0;
    return threadIt;
}

stThread *stThreadIt_getNext(stThreadIt *threadIt) {
    if (threadIt->index < stThreadSet_getSize(threadIt->threadSet)) {
        return stList_get(threadIt->threadSet->threads, threadIt->index++);
    }
    return NULL;
}

stSegment *getNextBlockSegment(stSegment *segment) {
    while (1) {
        if (segment == NULL) {
            return NULL;
        }
        stBlock *block;
        if ((block = stSegment_getBlock(segment)) != NULL && stBlock_getFirst(
                block) == segment) {
            return segment;
        }
        segment = stSegment_get3Prime(segment);
    }
}

void stThreadSet_joinTrivialBoundaries(stThreadSet *threadSet) {
    stThreadIt threadIt = stThreadSet_getIterator(threadSet);
    stThread *thread;
    while ((thread = stThreadIt_getNext(&threadIt)) != NULL) {
        stThread_joinTrivialBoundaries(thread);
        stSegment *segment = getNextBlockSegment(stThread_getFirst(thread));
        while (segment != NULL) {
            stEnd end;
            end.block = stSegment_getBlock(segment);
            end.orientation = 0;
            if (stEnd_boundaryIsTrivial(end)) {
                stEnd_joinTrivialBoundary(end);
            }
            end.orientation = 1;
            if (stEnd_boundaryIsTrivial(end)) {
                stEnd_joinTrivialBoundary(end);
            }
            segment = getNextBlockSegment(stSegment_get3Prime(segment));
        }
    }
}

stSegment *stThreadSet_getSegment(stThreadSet *threadSet, int64_t name,
        int64_t coordinate) {
    stThread *thread = stThreadSet_getThread(threadSet, name);
    if (thread == NULL) {
        return NULL;
    }
    return stThread_getSegment(thread, coordinate);
}

//convenience functions

stThreadSetSegmentIt stThreadSet_getSegmentIt(stThreadSet *threadSet) {
    stThreadSetSegmentIt segmentIt;
    segmentIt.threadIt = stThreadSet_getIterator(threadSet);
    segmentIt.segment = NULL;
    return segmentIt;
}

stSegment *stThreadSetSegmentIt_getNext(stThreadSetSegmentIt *segmentIt) {
    while (segmentIt->segment == NULL) {
        stThread *thread = stThreadIt_getNext(&segmentIt->threadIt);
        if (thread == NULL) {
            return NULL;
        }
        segmentIt->segment = stThread_getFirst(thread);
    }
    stSegment *segment = segmentIt->segment;
    segmentIt->segment = stSegment_get3Prime(segment);
    return segment;
}

stThreadSetBlockIt stThreadSet_getBlockIt(stThreadSet *threadSet) {
    stThreadSetBlockIt blockIt;
    blockIt.segmentIt = stThreadSet_getSegmentIt(threadSet);
    return blockIt;
}

stBlock *stThreadSetBlockIt_getNext(stThreadSetBlockIt *blockIt) {
    while (1) {
        stSegment *segment =
                stThreadSetSegmentIt_getNext(&(blockIt->segmentIt));
        if (segment == NULL) {
            return NULL;
        }
        stBlock *block;
        if ((block = stSegment_getBlock(segment)) != NULL && stBlock_getFirst(
                block) == segment) {
            return block;
        }
    }
}

void stThreadSet_getAdjacencyComponentsP2(stHash *endsToAdjacencyComponents,
        stList *adjacencyComponent, stEnd *end) {
    stList *stack = stList_construct();
    stList_append(adjacencyComponent, end);
    stHash_insert(endsToAdjacencyComponents, end, adjacencyComponent);
    stList_append(stack, end);
    while (stList_length(stack) > 0) {
        end = stList_pop(stack);
        stBlockIt blockIt = stBlock_getSegmentIterator(end->block);
        stSegment *segment;
        stEnd end2;
        while ((segment = stBlockIt_getNext(&blockIt)) != NULL) {
            bool _5PrimeTraversal = stEnd_traverse5Prime(end->orientation,
                    segment);
            while (1) {
                segment = _5PrimeTraversal ? stSegment_get5Prime(segment)
                        : stSegment_get3Prime(segment);
                if (segment == NULL) {
                    break;
                }
                stBlock *block = stSegment_getBlock(segment);
                if (block != NULL) {
                    end2.block = block;
                    end2.orientation = stEnd_endOrientation(_5PrimeTraversal,
                            segment);
                    if (stHash_search(endsToAdjacencyComponents, &end2) == NULL) {
                        stEnd *end3 = stEnd_construct(end2.block,
                                end2.orientation);
                        stList_append(adjacencyComponent, end3);
                        stHash_insert(endsToAdjacencyComponents, end3,
                                adjacencyComponent);
                        stList_append(stack, end3);
                    }
                    break;
                }
            }
        }
    }
    stList_destruct(stack);
}

void stThreadSet_getAdjacencyComponentsP(stHash *endsToAdjacencyComponents,
        stList *adjacencyComponents, stBlock *block, bool orientation) {
    stEnd end;
    end.orientation = orientation;
    end.block = block;
    stList *adjacencyComponent = stHash_search(endsToAdjacencyComponents, &end);
    if (adjacencyComponent == NULL) {
        adjacencyComponent = stList_construct3(0,
                (void(*)(void *)) stEnd_destruct);
        stList_append(adjacencyComponents, adjacencyComponent);
        stThreadSet_getAdjacencyComponentsP2(endsToAdjacencyComponents,
                adjacencyComponent, stEnd_construct(block, orientation));
    }
}

stList *stThreadSet_getAdjacencyComponents(stThreadSet *threadSet) {
    //stSortedSet *seen = stSortedSet_construct();
    stHash *endsToAdjacencyComponents = stHash_construct3(stEnd_hashFn,
            stEnd_equalsFn, NULL, NULL);
    stList *adjacencyComponents = stList_construct3(0,
            (void(*)(void *)) stList_destruct);
    stThreadSetBlockIt blockIt = stThreadSet_getBlockIt(threadSet);
    stBlock *block;
    while ((block = stThreadSetBlockIt_getNext(&blockIt)) != NULL) {
        stThreadSet_getAdjacencyComponentsP(endsToAdjacencyComponents,
                adjacencyComponents, block, 0);
        stThreadSet_getAdjacencyComponentsP(endsToAdjacencyComponents,
                adjacencyComponents, block, 1);
    }
    stHash_destruct(endsToAdjacencyComponents);
    return adjacencyComponents;
}

//stEnd

//Block ends

stEnd *stEnd_construct(stBlock *block, bool orientation) {
    stEnd *end = st_malloc(sizeof(stEnd));
    end->block = block;
    end->orientation = orientation;
    return end;
}

void stEnd_destruct(stEnd *end) {
    free(end);
}

stBlock *stEnd_getBlock(stEnd *end) {
    return end->block;
}

bool stEnd_getOrientation(stEnd *end) {
    return end->orientation;
}

int stEnd_equalsFn(const void *a, const void *b) {
    const stEnd *end1 = a, *end2 = b;
    return end1->block == end2->block && end1->orientation == end2->orientation;
}

uint32_t stEnd_hashFn(const void *a) {
    const stEnd *end1 = a;
    int64_t i = (int64_t) end1->block;
    return (uint32_t) i;
}

bool stEnd_traverse5Prime(bool endOrientation, stSegment *segment) {
    return !(stSegment_getBlockOrientation(segment) ^ endOrientation);
}

bool stEnd_endOrientation(bool _5PrimeTraversal, stSegment *segment) {
    return _5PrimeTraversal ^ stSegment_getBlockOrientation(segment);
}

bool stEnd_boundaryIsTrivial(stEnd end) {
    stBlockIt segmentIt = stBlock_getSegmentIterator(end.block);
    stSegment *segment = stBlockIt_getNext(&segmentIt);
    bool _5PrimeTraversal = stEnd_traverse5Prime(end.orientation, segment);
    segment = _5PrimeTraversal ? stSegment_get5Prime(segment)
            : stSegment_get3Prime(segment);
    stBlock *block;
    if (segment == NULL || (block = stSegment_getBlock(segment)) == NULL
            || block == end.block || stBlock_getDegree(block)
            != stBlock_getDegree(end.block)) {
        return 0;
    }
    bool endOrientation = stEnd_endOrientation(_5PrimeTraversal, segment);
    while ((segment = stBlockIt_getNext(&segmentIt)) != NULL) {
        _5PrimeTraversal = stEnd_traverse5Prime(end.orientation, segment);
        segment = _5PrimeTraversal ? stSegment_get5Prime(segment)
                : stSegment_get3Prime(segment);
        if (segment == NULL) {
            return 0;
        }
        stBlock *block2 = stSegment_getBlock(segment);
        if (block2 == NULL || block != block2 || stEnd_endOrientation(
                _5PrimeTraversal, segment) != endOrientation) {
            return 0;
        }
    }
    return 1;
}

static void merge3Prime(stSegment *segment) {
    stSegment *nSegment = segment->nSegment;
    assert(nSegment != NULL && nSegment != segment);
    stSortedSet_remove(segment->thread->segments, nSegment);
    assert(nSegment->block == NULL);
    assert(nSegment->nSegment != NULL);
    segment->nSegment = nSegment->nSegment;
    nSegment->nSegment->pSegment = segment;
    stSegment_destruct(nSegment);
}

static void merge5Prime(stSegment *segment) {
    stSegment *pSegment = segment->pSegment;
    assert(pSegment != NULL && pSegment != segment);
    stSortedSet_remove(segment->thread->segments, pSegment);
    assert(pSegment->block == NULL);
    segment->pSegment = pSegment->pSegment;
    if (pSegment->pSegment != NULL) {
        pSegment->pSegment->nSegment = segment;
    }
    assert(pSegment->start < segment->start);
    segment->start = pSegment->start;
    stSegment_destruct(pSegment);
}

void stEnd_joinTrivialBoundary(stEnd end) {
    stSegment *segment = stBlock_getFirst(end.block);
    assert(segment != NULL);
    bool _5PrimeTraversal = stEnd_traverse5Prime(end.orientation, segment);
    segment = _5PrimeTraversal ? stSegment_get5Prime(segment)
            : stSegment_get3Prime(segment);
    assert(
            segment != NULL && stSegment_getBlock(segment) != NULL
                    && stSegment_getBlock(segment) != end.block);
    stBlock_destruct(stSegment_getBlock(segment)); //get rid of the old block
    stBlockIt segmentIt = stBlock_getSegmentIterator(end.block);
    while ((segment = stBlockIt_getNext(&segmentIt)) != NULL) {
        bool _5PrimeTraversal = stEnd_traverse5Prime(end.orientation, segment);
        if (_5PrimeTraversal) {
            merge5Prime(segment);
        } else {
            merge3Prime(segment);
        }
    }
}

//stPinch

void stPinch_fillOut(stPinch *pinch, int64_t name1, int64_t name2,
        int64_t start1, int64_t start2, int64_t length, bool strand) {
    pinch->name1 = name1;
    pinch->name2 = name2;
    pinch->start1 = start1;
    pinch->start2 = start2;
    pinch->length = length;
    pinch->strand = strand;
}

stPinch *stPinch_construct(int64_t name1, int64_t name2, int64_t start1,
        int64_t start2, int64_t length, bool strand) {
    stPinch *pinch = st_malloc(sizeof(stPinch));
    stPinch_fillOut(pinch, name1, name2, start1, start2, length, strand);
    return pinch;
}

void stPinch_destruct(stPinch *pinch) {
    free(pinch);
}

