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

struct _stPinchThreadSet {
    stList *threads;
    stHash *threadsHash;
};

struct _stPinchThread {
    int64_t name;
    int64_t start;
    int64_t length;
    stSortedSet *segments;
};

struct _stPinchSegment {
    stPinchThread *thread;
    int64_t start;
    stPinchSegment *pSegment;
    stPinchSegment *nSegment;
    stPinchBlock *block;
    bool blockOrientation;
    stPinchSegment *nBlockSegment;
};

struct _stPinchBlock {
    uint32_t degree;
    stPinchSegment *headSegment;
    stPinchSegment *tailSegment;
};

//Blocks

static void connectBlockToSegment(stPinchSegment *segment, bool orientation, stPinchBlock *block, stPinchSegment *nBlockSegment) {
    segment->block = block;
    segment->blockOrientation = orientation;
    segment->nBlockSegment = nBlockSegment;
}

stPinchBlock *stPinchBlock_construct3(stPinchSegment *segment, bool orientation) {
    stPinchBlock *block = st_malloc(sizeof(stPinchBlock));
    block->headSegment = segment;
    block->tailSegment = segment;
    connectBlockToSegment(segment, orientation, block, NULL);
    block->degree = 1;
    return block;
}

stPinchBlock *stPinchBlock_construct2(stPinchSegment *segment) {
    return stPinchBlock_construct3(segment, 1);
}

stPinchBlock *stPinchBlock_construct(stPinchSegment *segment1, bool orientation1, stPinchSegment *segment2, bool orientation2) {
    if (stPinchSegment_getLength(segment1) != stPinchSegment_getLength(segment2)) {
        stThrowNew(ST_COMPRESSION_EXCEPTION_ID, "Two segments that are being pinched have different lengths: %i %i\n",
                stPinchSegment_getLength(segment1), stPinchSegment_getLength(segment2));
    }
    stPinchBlock *block = st_malloc(sizeof(stPinchBlock));
    block->headSegment = segment1;
    block->tailSegment = segment2;
    connectBlockToSegment(segment1, orientation1, block, segment2);
    connectBlockToSegment(segment2, orientation2, block, NULL);
    block->degree = 2;
    return block;
}

void stPinchBlock_destruct(stPinchBlock *block) {
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment = stPinchBlockIt_getNext(&blockIt);
    while (segment != NULL) {
        stPinchSegment *nSegment = stPinchBlockIt_getNext(&blockIt);
        connectBlockToSegment(segment, 0, NULL, NULL);
        segment = nSegment;
    }
    free(block);
}

stPinchBlock *stPinchBlock_pinch(stPinchBlock *block1, stPinchBlock *block2, bool orientation) {
    if (block1 == block2) {
        return block1; //Already joined
    }
    if (stPinchBlock_getDegree(block1) > stPinchBlock_getDegree(block2)) { //Avoid merging large blocks into small blocks
        return stPinchBlock_pinch(block2, block1, orientation);
    }
    if (stPinchBlock_getLength(block1) != stPinchBlock_getLength(block2)) {
        stThrowNew(ST_PINCH_GRAPH_EXCEPTION_ID, "Two segments that are being pinched have different lengths: %i %i\n",
                stPinchBlock_getLength(block1), stPinchBlock_getLength(block2));
    }
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block2);
    stPinchSegment *segment = stPinchBlockIt_getNext(&blockIt);
    while (segment != NULL) {
        stPinchSegment *nSegment = stPinchBlockIt_getNext(&blockIt);
        bool segmentOrientation = stPinchSegment_getBlockOrientation(segment);
        stPinchBlock_pinch2(block1, segment, (segmentOrientation && orientation) || (!segmentOrientation && !orientation));
        segment = nSegment;
    }
    free(block2);
    return block1;
}

stPinchBlock *stPinchBlock_pinch2(stPinchBlock *block, stPinchSegment *segment, bool orientation) {
    assert(block->tailSegment != NULL);
    assert(block->tailSegment->nBlockSegment == NULL);
    block->tailSegment->nBlockSegment = segment;
    connectBlockToSegment(segment, orientation, block, NULL);
    block->tailSegment = segment;
    block->degree++;
    return block;
}

stPinchBlockIt stPinchBlock_getSegmentIterator(stPinchBlock *block) {
    stPinchBlockIt blockIt;
    blockIt.segment = block->headSegment;
    return blockIt;
}

stPinchSegment *stPinchBlockIt_getNext(stPinchBlockIt *blockIt) {
    stPinchSegment *segment = blockIt->segment;
    if (segment != NULL) {
        blockIt->segment = segment->nBlockSegment;
    }
    return segment;
}

uint32_t stPinchBlock_getDegree(stPinchBlock *block) {
    return block->degree;
}

stPinchSegment *stPinchBlock_getFirst(stPinchBlock *block) {
    assert(block->headSegment != NULL);
    return block->headSegment;
}

int64_t stPinchBlock_getLength(stPinchBlock *block) {
    return stPinchSegment_getLength(stPinchBlock_getFirst(block));
}

void stPinchBlock_trim(stPinchBlock *block, int32_t blockEndTrim) {
    if (blockEndTrim <= 0) {
        return;
    }
    if (stPinchBlock_getLength(block) > 2 * blockEndTrim) {
        stPinchSegment *segment = stPinchBlock_getFirst(block);
        stPinchSegment_split(segment, stPinchSegment_getStart(segment) + blockEndTrim - 1);
        block = stPinchSegment_getBlock(segment);
        assert(block != NULL);
        stPinchBlock_destruct(block);
        segment = stPinchSegment_get3Prime(segment);
        assert(segment != NULL);
        assert(stPinchSegment_getBlock(segment) != NULL);
        stPinchSegment_split(segment, stPinchSegment_getStart(segment) + stPinchSegment_getLength(segment) - 1 - blockEndTrim);
        segment = stPinchSegment_get3Prime(segment);
        assert(segment != NULL);
        block = stPinchSegment_getBlock(segment);
        assert(block != NULL);
        stPinchBlock_destruct(block);
    } else { //Two short, so we just destroy it
        stPinchBlock_destruct(block);
    }
}

//Segments

int64_t stPinchSegment_getStart(stPinchSegment *segment) {
    return segment->start;
}

int64_t stPinchSegment_getLength(stPinchSegment *segment) {
    assert(segment->nSegment != NULL);
    return segment->nSegment->start - segment->start;
}

stPinchBlock *stPinchSegment_getBlock(stPinchSegment *segment) {
    return segment->block;
}

bool stPinchSegment_getBlockOrientation(stPinchSegment *segment) {
    return segment->blockOrientation;
}

stPinchSegment *stPinchSegment_get5Prime(stPinchSegment *segment) {
    return segment->pSegment;
}

stPinchSegment *stPinchSegment_get3Prime(stPinchSegment *segment) {
    return segment->nSegment->nSegment != NULL ? segment->nSegment : NULL;
}

int64_t stPinchSegment_getName(stPinchSegment *segment) {
    return stPinchThread_getName(segment->thread);
}

stPinchThread *stPinchSegment_getThread(stPinchSegment *segment) {
    return segment->thread;
}

//Private segment functions

void stPinchSegment_destruct(stPinchSegment *segment) {
    if (stPinchSegment_getBlock(segment) != NULL) {
        stPinchBlock_destruct(stPinchSegment_getBlock(segment));
    }
    free(segment);
}

int stPinchSegment_compareFunction(const stPinchSegment *segment1, const stPinchSegment *segment2) {
    return segment1->start < segment2->start ? -1 : (segment1->start > segment2->start ? 1 : 0);
}

static stPinchSegment *stPinchSegment_construct(int64_t start, stPinchThread *thread) {
    stPinchSegment *segment = st_calloc(1, sizeof(stPinchSegment));
    segment->start = start;
    segment->thread = thread;
    return segment;
}

static stPinchSegment *stPinchSegment_splitP(stPinchSegment *segment, int64_t leftBlockLength) {
    stPinchSegment *nSegment = segment->nSegment;
    assert(nSegment != NULL);
    stPinchSegment *rightSegment = stPinchSegment_construct(stPinchSegment_getStart(segment) + leftBlockLength, segment->thread);
    segment->nSegment = rightSegment;
    rightSegment->pSegment = segment;
    rightSegment->nSegment = nSegment;
    nSegment->pSegment = rightSegment;
    stSortedSet_insert(segment->thread->segments, rightSegment);
    return rightSegment;
}

void stPinchSegment_split(stPinchSegment *segment, int64_t leftSideOfSplitPoint) {
    if (leftSideOfSplitPoint == stPinchSegment_getStart(segment) + stPinchSegment_getLength(segment) - 1) { //There is already a break
        return;
    }
    int64_t leftSegmentLength = leftSideOfSplitPoint - stPinchSegment_getStart(segment) + 1;
    assert(leftSegmentLength > 0);
    stPinchBlock *block;
    if ((block = stPinchSegment_getBlock(segment)) != NULL) {
        int64_t rightSegmentLength = stPinchBlock_getLength(block) - leftSegmentLength;
        assert(rightSegmentLength > 0);
        if (!stPinchSegment_getBlockOrientation(segment)) {
            int64_t i = rightSegmentLength;
            rightSegmentLength = leftSegmentLength;
            leftSegmentLength = i;
        }
        stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
        segment = stPinchBlockIt_getNext(&blockIt);
        assert(segment != NULL);
        stPinchSegment *pSegment = NULL;
        stPinchBlock *block2;
        if (stPinchSegment_getBlockOrientation(segment)) {
            stPinchSegment *segment2 = stPinchSegment_splitP(segment, leftSegmentLength);
            block2 = stPinchBlock_construct2(segment2);
            pSegment = segment;
        } else {
            stPinchSegment *segment2 = stPinchSegment_splitP(segment, rightSegmentLength);
            block->headSegment = segment2;
            connectBlockToSegment(segment2, 0, block, segment->nBlockSegment);
            if (segment2->nBlockSegment == NULL) {
                block->tailSegment = segment2;
            }
            block2 = stPinchBlock_construct2(segment);
            segment->blockOrientation = 0; //This gets sets positive by default.
            pSegment = segment2;
        }
        while ((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
            if (stPinchSegment_getBlockOrientation(segment)) {
                stPinchSegment *segment2 = stPinchSegment_splitP(segment, leftSegmentLength);
                stPinchBlock_pinch2(block2, segment2, 1);
                pSegment = segment;
            } else {
                stPinchSegment *segment2 = stPinchSegment_splitP(segment, rightSegmentLength);
                pSegment->nBlockSegment = segment2;
                connectBlockToSegment(segment2, 0, block, segment->nBlockSegment);
                if (segment2->nBlockSegment == NULL) {
                    block->tailSegment = segment2;
                }
                stPinchBlock_pinch2(block2, segment, 0);
                pSegment = segment2;
            }
        }
    } else {
        stPinchSegment_splitP(segment, leftSegmentLength);
    }
}

//Thread

int64_t stPinchThread_getName(stPinchThread *thread) {
    return thread->name;
}

int64_t stPinchThread_getStart(stPinchThread *thread) {
    return thread->start;
}

int64_t stPinchThread_getLength(stPinchThread *thread) {
    return thread->length;
}

stPinchSegment *stPinchThread_getSegment(stPinchThread *thread, int64_t coordinate) {
    static stPinchSegment segment;
    segment.start = coordinate;
    stPinchSegment *segment2 = stSortedSet_searchLessThanOrEqual(thread->segments, &segment);
    if (segment2 == NULL) {
        return NULL;
    }
    assert(stPinchSegment_getStart(segment2) <= coordinate);
    if (stPinchSegment_getStart(segment2) + stPinchSegment_getLength(segment2) <= coordinate) {
        return NULL;
    }
    return segment2;
}

stPinchSegment *stPinchThread_getFirst(stPinchThread *thread) {
    return stSortedSet_getFirst(thread->segments);
}

stPinchSegment *stPinchThread_getLast(stPinchThread *thread) {
    return stSortedSet_getLast(thread->segments);
}

void stPinchThread_split(stPinchThread *thread, int64_t leftSideOfSplitPoint) {
    stPinchSegment *segment = stPinchThread_getSegment(thread, leftSideOfSplitPoint);
    if (segment == NULL) {
        return;
    }
    stPinchSegment_split(segment, leftSideOfSplitPoint);
}

void stPinchThread_joinTrivialBoundaries(stPinchThread *thread) {
    stPinchSegment *segment = stPinchThread_getFirst(thread);
    do {
        if (stPinchSegment_getBlock(segment) == NULL) {
            while (1) {
                stPinchSegment *nSegment = stPinchSegment_get3Prime(segment);
                if (nSegment != NULL) {
                    stPinchBlock *nBlock = stPinchSegment_getBlock(nSegment);
                    if (nBlock == NULL) {
                        //Trivial join
                        segment->nSegment = nSegment->nSegment;
                        assert(nSegment->nSegment != NULL);
                        nSegment->nSegment->pSegment = segment;
                        stSortedSet_remove(thread->segments, nSegment);
                        stPinchSegment_destruct(nSegment);
                        continue;
                    }
                }
                break;
            }
        }
    } while ((segment = stPinchSegment_get3Prime(segment)) != NULL);
}

stPinchSegment *stPinchThread_pinchP(stPinchThread *thread, int64_t start) {
    stPinchSegment *segment1 = stPinchThread_getSegment(thread, start);
    assert(segment1 != NULL);
    if (stPinchSegment_getStart(segment1) != start) {
        stPinchSegment_split(segment1, start - 1);
        segment1 = stPinchSegment_get3Prime(segment1);
    }
    return segment1;
}

stPinchSegment *stPinchThread_pinchTrim(stPinchSegment *segment, bool strand, int32_t length) {
    assert(length > 0);
    if (stPinchSegment_getLength(segment) <= length) {
        return segment;
    }
    if (strand) {
        stPinchSegment_split(segment, stPinchSegment_getStart(segment) + length - 1);
        return segment;
    }
    stPinchSegment_split(segment, stPinchSegment_getStart(segment) + stPinchSegment_getLength(segment) - 1 - length);
    return stPinchSegment_get3Prime(segment);
}

void stPinchThread_pinch(stPinchThread *thread1, stPinchThread *thread2, int64_t start1, int64_t start2, int64_t length, bool strand2) {
    assert(length >= 0);
    if (length == 0) {
        return;
    }
    assert(stPinchThread_getStart(thread1) <= start1);
    assert(stPinchThread_getStart(thread1) + stPinchThread_getLength(thread1) >= start1 + length);
    assert(stPinchThread_getStart(thread2) <= start2);
    assert(stPinchThread_getStart(thread2) + stPinchThread_getLength(thread2) >= start2 + length);
    stPinchSegment *segment1 = stPinchThread_pinchP(thread1, start1), *segment2;
    if (strand2) {
        segment2 = stPinchThread_pinchP(thread2, start2);
    } else {
        segment2 = stPinchThread_getSegment(thread2, start2 + length - 1);
        stPinchSegment_split(segment2, start2 + length - 1);
    }
    while (length > 0) {
        if (segment1 == segment2) {
            if (strand2) {
                return; //This is a trivial alignment
            }
            if (stPinchSegment_getLength(segment1) > 1) { //Split the block in two
                segment2 = stPinchThread_pinchTrim(segment2, 0, stPinchSegment_getLength(segment1) / 2);
            }
        }
        do {
            int64_t i = stPinchSegment_getLength(segment1);
            segment2 = stPinchThread_pinchTrim(segment2, strand2, length > i ? i : length);
            i = stPinchSegment_getLength(segment2);
            stPinchThread_pinchTrim(segment1, 1, length > i ? i : length);
        } while (stPinchSegment_getLength(segment1) != stPinchSegment_getLength(segment2));
        stPinchBlock *block1, *block2;
        if ((block1 = stPinchSegment_getBlock(segment1)) == NULL) {
            block1 = stPinchBlock_construct2(segment1);
        }
        if ((block2 = stPinchSegment_getBlock(segment2)) == NULL) {
            block2 = stPinchBlock_construct2(segment2);
        }
        bool bO1 = stPinchSegment_getBlockOrientation(segment1);
        bool bO2 = stPinchSegment_getBlockOrientation(segment2);
        bool alignmentOrientation = ((bO1 == bO2) && strand2) || ((bO1 != bO2) && !strand2);
        if (block1 == block2) {
            if (stPinchSegment_getLength(segment1) > 1 && !alignmentOrientation) {
                segment2 = stPinchThread_pinchTrim(segment2, strand2, stPinchSegment_getLength(segment2) / 2);
                continue;
            }
        }
        block1 = stPinchBlock_pinch(block1, block2, alignmentOrientation);
        length -= stPinchSegment_getLength(segment1);
        segment1 = stPinchSegment_get3Prime(segment1);
        segment2 = strand2 ? stPinchSegment_get3Prime(segment2) : stPinchSegment_get5Prime(segment2);
    }
}

//Private functions

static stPinchThread *stPinchThread_construct(int64_t name, int64_t start, int64_t length) {
    stPinchThread *thread = st_malloc(sizeof(stPinchThread));
    thread->name = name;
    thread->start = start;
    thread->length = length;
    thread->segments = stSortedSet_construct3((int(*)(const void *, const void *)) stPinchSegment_compareFunction,
            (void(*)(void *)) stPinchSegment_destruct);
    stPinchSegment *segment = stPinchSegment_construct(start, thread);
    stPinchSegment *terminatorSegment = stPinchSegment_construct(start + length, thread);
    segment->nSegment = terminatorSegment;
    terminatorSegment->pSegment = segment;
    stSortedSet_insert(thread->segments, segment);
    return thread;
}

static void stPinchThread_destruct(stPinchThread *thread) {
    stPinchSegment *segment = stPinchThread_getLast(thread);
    free(segment->nSegment);
    stSortedSet_destruct(thread->segments);
    free(thread);
}

static uint32_t stPinchThread_hashKey(const stPinchThread *thread) {
    return thread->name;
}

static int stPinchThread_equals(const stPinchThread *thread1, const stPinchThread *thread2) {
    return thread1->name == thread2->name;
}

//Thread set

stPinchThreadSet *stPinchThreadSet_construct() {
    stPinchThreadSet *threadSet = st_malloc(sizeof(stPinchThreadSet));
    threadSet->threads = stList_construct3(0, (void(*)(void *)) stPinchThread_destruct);
    threadSet->threadsHash = stHash_construct3((uint32_t(*)(const void *)) stPinchThread_hashKey,
            (int(*)(const void *, const void *)) stPinchThread_equals, NULL, NULL);
    return threadSet;
}

void stPinchThreadSet_destruct(stPinchThreadSet *threadSet) {
    stList_destruct(threadSet->threads);
    stHash_destruct(threadSet->threadsHash);
    free(threadSet);
}

stPinchThread *stPinchThreadSet_addThread(stPinchThreadSet *threadSet, int64_t name, int64_t start, int64_t length) {
    stPinchThread *thread = stPinchThread_construct(name, start, length);
    assert(stPinchThreadSet_getThread(threadSet, name) == NULL);
    stHash_insert(threadSet->threadsHash, thread, thread);
    stList_append(threadSet->threads, thread);
    return thread;
}

stPinchThread *stPinchThreadSet_getThread(stPinchThreadSet *threadSet, int64_t name) {
    static stPinchThread thread;
    thread.name = name;
    return stHash_search(threadSet->threadsHash, &thread);
}

int32_t stPinchThreadSet_getSize(stPinchThreadSet *threadSet) {
    return stList_length(threadSet->threads);
}

stPinchThreadSetIt stPinchThreadSet_getIt(stPinchThreadSet *threadSet) {
    stPinchThreadSetIt threadIt;
    threadIt.threadSet = threadSet;
    threadIt.index = 0;
    return threadIt;
}

stPinchThread *stPinchThreadSetIt_getNext(stPinchThreadSetIt *threadIt) {
    if (threadIt->index < stPinchThreadSet_getSize(threadIt->threadSet)) {
        return stList_get(threadIt->threadSet->threads, threadIt->index++);
    }
    return NULL;
}

void stPinchThreadSet_joinTrivialBoundaries(stPinchThreadSet *threadSet) {
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *thread;
    while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
        stPinchThread_joinTrivialBoundaries(thread);
    }
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *block;
    while ((block = stPinchThreadSetBlockIt_getNext(&blockIt))) {
        stPinchEnd end = stPinchEnd_constructStatic(block, 0);
        if (stPinchEnd_boundaryIsTrivial(end)) {
            stPinchEnd_joinTrivialBoundary(end);
        }
        end.orientation = 1;
        if (stPinchEnd_boundaryIsTrivial(end)) {
            stPinchEnd_joinTrivialBoundary(end);
        }
    }
}

stPinchSegment *stPinchThreadSet_getSegment(stPinchThreadSet *threadSet, int64_t name, int64_t coordinate) {
    stPinchThread *thread = stPinchThreadSet_getThread(threadSet, name);
    if (thread == NULL) {
        return NULL;
    }
    return stPinchThread_getSegment(thread, coordinate);
}

//convenience functions

stPinchThreadSetSegmentIt stPinchThreadSet_getSegmentIt(stPinchThreadSet *threadSet) {
    stPinchThreadSetSegmentIt segmentIt;
    segmentIt.threadIt = stPinchThreadSet_getIt(threadSet);
    segmentIt.segment = NULL;
    return segmentIt;
}

stPinchSegment *stPinchThreadSetSegmentIt_getNext(stPinchThreadSetSegmentIt *segmentIt) {
    if (segmentIt->segment != NULL) {
        segmentIt->segment = stPinchSegment_get3Prime(segmentIt->segment);
    }
    while (segmentIt->segment == NULL) {
        stPinchThread *thread = stPinchThreadSetIt_getNext(&segmentIt->threadIt);
        if (thread == NULL) {
            return NULL;
        }
        segmentIt->segment = stPinchThread_getFirst(thread);
    }
    return segmentIt->segment;
}

stPinchThreadSetBlockIt stPinchThreadSet_getBlockIt(stPinchThreadSet *threadSet) {
    stPinchThreadSetBlockIt blockIt;
    blockIt.segmentIt = stPinchThreadSet_getSegmentIt(threadSet);
    return blockIt;
}

stPinchBlock *stPinchThreadSetBlockIt_getNext(stPinchThreadSetBlockIt *blockIt) {
    while (1) {
        stPinchSegment *segment = stPinchThreadSetSegmentIt_getNext(&(blockIt->segmentIt));
        if (segment == NULL) {
            return NULL;
        }
        stPinchBlock *block;
        if ((block = stPinchSegment_getBlock(segment)) != NULL && stPinchBlock_getFirst(block) == segment) {
            return block;
        }
    }
}

int32_t stPinchThreadSet_getTotalBlockNumber(stPinchThreadSet *threadSet) {
    int32_t blockCount = 0;
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    while (stPinchThreadSetBlockIt_getNext(&blockIt) != NULL) {
        blockCount++;
    }
    return blockCount;
}

void stPinchThreadSet_getAdjacencyComponentsP2(stHash *endsToAdjacencyComponents, stList *adjacencyComponent, stPinchEnd *end) {
    stList *stack = stList_construct();
    stList_append(adjacencyComponent, end);
    stHash_insert(endsToAdjacencyComponents, end, adjacencyComponent);
    stList_append(stack, end);
    while (stList_length(stack) > 0) {
        end = stList_pop(stack);
        stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(end->block);
        stPinchSegment *segment;
        while ((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
            bool _5PrimeTraversal = stPinchEnd_traverse5Prime(end->orientation, segment);
            while (1) {
                segment = _5PrimeTraversal ? stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);
                if (segment == NULL) {
                    break;
                }
                stPinchBlock *block = stPinchSegment_getBlock(segment);
                if (block != NULL) {
                    stPinchEnd end2 = stPinchEnd_constructStatic(block, stPinchEnd_endOrientation(_5PrimeTraversal, segment));
                    if (stHash_search(endsToAdjacencyComponents, &end2) == NULL) {
                        stPinchEnd *end3 = stPinchEnd_construct(end2.block, end2.orientation);
                        stList_append(adjacencyComponent, end3);
                        stHash_insert(endsToAdjacencyComponents, end3, adjacencyComponent);
                        stList_append(stack, end3);
                    }
                    break;
                }
            }
        }
    }
    stList_destruct(stack);
}

void stPinchThreadSet_getAdjacencyComponentsP(stHash *endsToAdjacencyComponents, stList *adjacencyComponents, stPinchBlock *block,
        bool orientation) {
    stPinchEnd end = stPinchEnd_constructStatic(block, orientation);
    stList *adjacencyComponent = stHash_search(endsToAdjacencyComponents, &end);
    if (adjacencyComponent == NULL) {
        adjacencyComponent = stList_construct3(0, (void(*)(void *)) stPinchEnd_destruct);
        stList_append(adjacencyComponents, adjacencyComponent);
        stPinchThreadSet_getAdjacencyComponentsP2(endsToAdjacencyComponents, adjacencyComponent, stPinchEnd_construct(block, orientation));
    }
}

stList *stPinchThreadSet_getAdjacencyComponents2(stPinchThreadSet *threadSet, stHash **endsToAdjacencyComponents) {
    *endsToAdjacencyComponents = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, NULL, NULL);
    stList *adjacencyComponents = stList_construct3(0, (void(*)(void *)) stList_destruct);
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *block;
    while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
        stPinchThreadSet_getAdjacencyComponentsP(*endsToAdjacencyComponents, adjacencyComponents, block, 0);
        stPinchThreadSet_getAdjacencyComponentsP(*endsToAdjacencyComponents, adjacencyComponents, block, 1);
    }
    return adjacencyComponents;
}

stList *stPinchThreadSet_getAdjacencyComponents(stPinchThreadSet *threadSet) {
    stHash *endsToAdjacencyComponents;
    stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents2(threadSet, &endsToAdjacencyComponents);
    stHash_destruct(endsToAdjacencyComponents);
    return adjacencyComponents;
}

stSortedSet *stPinchThreadSet_getThreadComponents(stPinchThreadSet *threadSet) {
    stHash *threadToComponentHash = stHash_construct();
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *thread;
    //Make a component for each thread.
    while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
        stList *threadComponent = stList_construct();
        stList_append(threadComponent, thread);
        stHash_insert(threadToComponentHash, thread, threadComponent);
    }
    //Now join components progressively according to blocks
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *block;
    while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
        stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(block);
        stPinchSegment *segment = stPinchBlockIt_getNext(&segmentIt);
        assert(segment != NULL);
        stList *threadComponent1 = stHash_search(threadToComponentHash, stPinchSegment_getThread(segment));
        assert(threadComponent1 != NULL);
        while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
            stList *threadComponent2 = stHash_search(threadToComponentHash, stPinchSegment_getThread(segment));
            assert(threadComponent2 != NULL);
            if (threadComponent1 != threadComponent2) {
                stList_appendAll(threadComponent1, threadComponent2);
                for (int32_t i = 0; i < stList_length(threadComponent2); i++) {
                    stPinchThread *thread = stList_get(threadComponent2, i);
                    stHash_insert(threadToComponentHash, thread, threadComponent1);
                }
                stList_destruct(threadComponent2);
            }
        }
    }
    //Get a list of the components
    stList *values = stHash_getValues(threadToComponentHash);
    stSortedSet *threadComponentsSet = stList_getSortedSet(values, NULL);
    stSortedSet_setDestructor(threadComponentsSet, (void(*)(void *)) stList_destruct);
    //Cleanup
    stList_destruct(values);
    stHash_destruct(threadToComponentHash);
    return threadComponentsSet;
}

//stPinchEnd

//Block ends

void stPinchEnd_fillOut(stPinchEnd *end, stPinchBlock *block, bool orientation) {
    end->block = block;
    end->orientation = orientation;
}

stPinchEnd *stPinchEnd_construct(stPinchBlock *block, bool orientation) {
    stPinchEnd *end = st_malloc(sizeof(stPinchEnd));
    stPinchEnd_fillOut(end, block, orientation);
    return end;
}

stPinchEnd stPinchEnd_constructStatic(stPinchBlock *block, bool orientation) {
    stPinchEnd end;
    stPinchEnd_fillOut(&end, block, orientation);
    return end;
}

void stPinchEnd_destruct(stPinchEnd *end) {
    free(end);
}

stPinchBlock *stPinchEnd_getBlock(stPinchEnd *end) {
    return end->block;
}

bool stPinchEnd_getOrientation(stPinchEnd *end) {
    return end->orientation;
}

int stPinchEnd_equalsFn(const void *a, const void *b) {
    const stPinchEnd *end1 = a, *end2 = b;
    return end1->block == end2->block && end1->orientation == end2->orientation;
}

uint32_t stPinchEnd_hashFn(const void *a) {
    const stPinchEnd *end1 = a;
    int64_t i = (int64_t) end1->block;
    return (uint32_t) i;
}

bool stPinchEnd_traverse5Prime(bool endOrientation, stPinchSegment *segment) {
    return !(stPinchSegment_getBlockOrientation(segment) ^ endOrientation);
}

bool stPinchEnd_endOrientation(bool _5PrimeTraversal, stPinchSegment *segment) {
    return _5PrimeTraversal ^ stPinchSegment_getBlockOrientation(segment);
}

bool stPinchEnd_boundaryIsTrivial(stPinchEnd end) {
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(end.block);
    stPinchSegment *segment = stPinchBlockIt_getNext(&segmentIt);
    bool _5PrimeTraversal = stPinchEnd_traverse5Prime(end.orientation, segment);
    segment = _5PrimeTraversal ? stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);
    stPinchBlock *block;
    if (segment == NULL || (block = stPinchSegment_getBlock(segment)) == NULL || block == end.block || stPinchBlock_getDegree(block)
            != stPinchBlock_getDegree(end.block)) {
        return 0;
    }
    bool endOrientation = stPinchEnd_endOrientation(_5PrimeTraversal, segment);
    while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        _5PrimeTraversal = stPinchEnd_traverse5Prime(end.orientation, segment);
        segment = _5PrimeTraversal ? stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);
        if (segment == NULL) {
            return 0;
        }
        stPinchBlock *block2 = stPinchSegment_getBlock(segment);
        if (block2 == NULL || block != block2 || stPinchEnd_endOrientation(_5PrimeTraversal, segment) != endOrientation) {
            return 0;
        }
    }
    return 1;
}

static void merge3Prime(stPinchSegment *segment) {
    stPinchSegment *nSegment = segment->nSegment;
    assert(nSegment != NULL && nSegment != segment);
    stSortedSet_remove(segment->thread->segments, nSegment);
    assert(nSegment->block == NULL);
    assert(nSegment->nSegment != NULL);
    segment->nSegment = nSegment->nSegment;
    nSegment->nSegment->pSegment = segment;
    stPinchSegment_destruct(nSegment);
}

static void merge5Prime(stPinchSegment *segment) {
    stPinchSegment *pSegment = segment->pSegment;
    assert(pSegment != NULL && pSegment != segment);
    stSortedSet_remove(segment->thread->segments, pSegment);
    assert(pSegment->block == NULL);
    segment->pSegment = pSegment->pSegment;
    if (pSegment->pSegment != NULL) {
        pSegment->pSegment->nSegment = segment;
    }
    assert(pSegment->start < segment->start);
    segment->start = pSegment->start;
    stPinchSegment_destruct(pSegment);
}

void stPinchEnd_joinTrivialBoundary(stPinchEnd end) {
    stPinchSegment *segment = stPinchBlock_getFirst(end.block);
    assert(segment != NULL);
    bool _5PrimeTraversal = stPinchEnd_traverse5Prime(end.orientation, segment);
    segment = _5PrimeTraversal ? stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);
    assert(segment != NULL && stPinchSegment_getBlock(segment) != NULL && stPinchSegment_getBlock(segment) != end.block);
    stPinchBlock_destruct(stPinchSegment_getBlock(segment)); //get rid of the old block
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(end.block);
    while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        bool _5PrimeTraversal = stPinchEnd_traverse5Prime(end.orientation, segment);
        if (_5PrimeTraversal) {
            merge5Prime(segment);
        } else {
            merge3Prime(segment);
        }
    }
}

//stPinch

void stPinch_fillOut(stPinch *pinch, int64_t name1, int64_t name2, int64_t start1, int64_t start2, int64_t length, bool strand) {
    pinch->name1 = name1;
    pinch->name2 = name2;
    pinch->start1 = start1;
    pinch->start2 = start2;
    pinch->length = length;
    pinch->strand = strand;
}

stPinch *stPinch_construct(int64_t name1, int64_t name2, int64_t start1, int64_t start2, int64_t length, bool strand) {
    stPinch *pinch = st_malloc(sizeof(stPinch));
    stPinch_fillOut(pinch, name1, name2, start1, start2, length, strand);
    return pinch;
}

stPinch stPinch_constructStatic(int64_t name1, int64_t name2, int64_t start1, int64_t start2, int64_t length, bool strand) {
    stPinch pinch;
    stPinch_fillOut(&pinch, name1, name2, start1, start2, length, strand);
    return pinch;
}

void stPinch_destruct(stPinch *pinch) {
    free(pinch);
}

//stPinchInterval

void stPinchInterval_fillOut(stPinchInterval *pinchInterval, int64_t name, int64_t start, int64_t length, void *label) {
    pinchInterval->name = name;
    pinchInterval->start = start;
    pinchInterval->length = length;
    pinchInterval->label = label;
}

stPinchInterval stPinchInterval_constructStatic(int64_t name, int64_t start, int64_t length, void *label) {
    stPinchInterval interval;
    stPinchInterval_fillOut(&interval, name, start, length, label);
    return interval;
}

stPinchInterval *stPinchInterval_construct(int64_t name, int64_t start, int64_t length, void *label) {
    stPinchInterval *interval = st_malloc(sizeof(stPinchInterval));
    stPinchInterval_fillOut(interval, name, start, length, label);
    return interval;
}

int64_t stPinchInterval_getName(stPinchInterval *pinchInterval) {
    return pinchInterval->name;
}

int64_t stPinchInterval_getStart(stPinchInterval *pinchInterval) {
    return pinchInterval->start;
}

int64_t stPinchInterval_getLength(stPinchInterval *pinchInterval) {
    return pinchInterval->length;
}

void *stPinchInterval_getLabel(stPinchInterval *pinchInterval) {
    return pinchInterval->label;
}

void stPinchInterval_destruct(stPinchInterval *pinchInterval) {
    free(pinchInterval);
}

void stPinchThreadSet_getLabelIntervalsP2(stPinchThread *thread, stSortedSet *pinchIntervals, int64_t start, void *label) {
    int64_t end = stPinchThread_getLength(thread) + stPinchThread_getStart(thread);
    if (start < end) {
        stSortedSet_insert(pinchIntervals, stPinchInterval_construct(stPinchThread_getName(thread), start, end - start, label));
    }
}

void stPinchThreadSet_getLabelIntervalsP(stPinchThread *thread, stHash *pinchEndsToLabels, stSortedSet *pinchIntervals) {
    stPinchSegment *segment = stPinchThread_getFirst(thread);
    if (segment == NULL) {
        return;
    }
    int64_t start = stPinchSegment_getStart(segment);
    void *label = NULL;
    do {
        stPinchBlock *block;
        while ((block = stPinchSegment_getBlock(segment)) == NULL) {
            segment = stPinchSegment_get3Prime(segment);
            if (segment == NULL) {
                stPinchThreadSet_getLabelIntervalsP2(thread, pinchIntervals, start, label);
                return;
            }
        }
        stPinchEnd pinchEnd = stPinchEnd_constructStatic(block, !stPinchSegment_getBlockOrientation(segment));
        void *label2 = stHash_search(pinchEndsToLabels, &pinchEnd);
        assert(label2 != NULL);
        if (label == NULL) {
            pinchEnd.orientation = !pinchEnd.orientation;
            label = stHash_search(pinchEndsToLabels, &pinchEnd);
            assert(label != NULL);
        }
        //#ifdef BEN_DEBUG
        else {
            pinchEnd.orientation = !pinchEnd.orientation;
            assert(label == stHash_search(pinchEndsToLabels, &pinchEnd));
        }
        //#endif
        if (label != label2) {
            int64_t end = stPinchSegment_getStart(segment) + stPinchSegment_getLength(segment) / 2;
            if (start < end) {
                stSortedSet_insert(pinchIntervals, stPinchInterval_construct(stPinchThread_getName(thread), start, end - start, label));
            }
            start = end;
            label = label2;
        }
        segment = stPinchSegment_get3Prime(segment);
    } while (segment != NULL);
    stPinchThreadSet_getLabelIntervalsP2(thread, pinchIntervals, start, label);
}

stSortedSet *stPinchThreadSet_getLabelIntervals(stPinchThreadSet *threadSet, stHash *pinchEndsToLabels) {
    stSortedSet *pinchIntervals = stSortedSet_construct3((int(*)(const void *, const void *)) stPinchInterval_compareFunction,
            (void(*)(void *)) stPinchInterval_destruct);
    stPinchThread *thread;
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    while ((thread = stPinchThreadSetIt_getNext(&threadIt))) {
        stPinchThreadSet_getLabelIntervalsP(thread, pinchEndsToLabels, pinchIntervals);
    }
    return pinchIntervals;
}

static int cmp64s(int64_t i, int64_t j) {
    return i > j ? 1 : (i < j ? -1 : 0);
}

int stPinchInterval_compareFunction(const stPinchInterval *interval1, const stPinchInterval *interval2) {
    int i = cmp64s(interval1->name, interval2->name);
    if (i != 0) {
        return i;
    }
    i = cmp64s(interval1->start, interval2->start);
    if (i != 0) {
        return i;
    }
    return cmp64s(interval2->length, interval1->length);
}

stPinchInterval *stPinchIntervals_getInterval(stSortedSet *pinchIntervals, int64_t name, int64_t position) {
    stPinchInterval interval;
    stPinchInterval_fillOut(&interval, name, position, 1, NULL);
    stPinchInterval *interval2 = stSortedSet_searchLessThanOrEqual(pinchIntervals, &interval);
    if (interval2 == NULL || stPinchInterval_getName(interval2) != name || stPinchInterval_getStart(interval2) + stPinchInterval_getLength(
            interval2) <= position) {
        return NULL;
    }
    assert(stPinchInterval_getStart(interval2) <= position);
    assert(stPinchInterval_getStart(interval2) + stPinchInterval_getLength(interval2) > position);
    return interval2;
}

//Random pinch graph generation, used for testing

static void getRandomPosition(stPinchThreadSet *threadSet, stPinchThread **thread, int64_t *position, bool *strand) {
    *thread = st_randomChoice(threadSet->threads);
    *position = st_randomInt(stPinchThread_getStart(*thread), stPinchThread_getStart(*thread) + stPinchThread_getLength(*thread));
    *strand = st_random() > 0.5;
}

stPinch stPinchThreadSet_getRandomPinch(stPinchThreadSet *threadSet) {
    stPinchThread *thread1, *thread2;
    int64_t start1, start2;
    bool strand1, strand2;
    getRandomPosition(threadSet, &thread1, &start1, &strand1);
    getRandomPosition(threadSet, &thread2, &start2, &strand2);
    int32_t i = stPinchThread_getStart(thread1) + stPinchThread_getLength(thread1) - start1;
    int32_t j = stPinchThread_getStart(thread2) + stPinchThread_getLength(thread2) - start2;
    assert(i >= 0 && j >= 0);
    i = i > j ? j : i;
    int64_t length = i == 0 ? 0 : st_randomInt(0, i);
    stPinch pinch = stPinch_constructStatic(stPinchThread_getName(thread1), stPinchThread_getName(thread2), start1, start2, length, strand1 == strand2);
    return pinch;
}

stPinchThreadSet *stPinchThreadSet_getRandomEmptyGraph() {
    stPinchThreadSet *threadSet = stPinchThreadSet_construct();
    int32_t randomThreadNumber = st_randomInt(2, 10);
    for (int32_t threadIndex = 0; threadIndex < randomThreadNumber; threadIndex++) {
        int32_t start = st_randomInt(1, 100);
        int32_t length = st_randomInt(1, 100);
        int32_t threadName = threadIndex + 4;
        stPinchThreadSet_addThread(threadSet, threadName, start, length);
    }
    return threadSet;
}

stPinchThreadSet *stPinchThreadSet_getRandomGraph() {
    stPinchThreadSet *threadSet = stPinchThreadSet_getRandomEmptyGraph();
    //Randomly push them together, updating both sets, and checking that set of alignments is what we expect
    while (st_random() > 0.05) {
        stPinch pinch = stPinchThreadSet_getRandomPinch(threadSet);
        stPinchThread_pinch(stPinchThreadSet_getThread(threadSet, pinch.name1), stPinchThreadSet_getThread(threadSet, pinch.name2),
                pinch.start1, pinch.start2, pinch.length, pinch.strand);
    }
    return threadSet;
}

