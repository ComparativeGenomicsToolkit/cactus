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

stBlock *stBlock_construct(stSegment *segment1, bool orientation1,
        stSegment *segment2, bool orientation2) {
    if(stSegment_getLength(segment1) != stSegment_getLength(segment2)) {
        stThrowNew(ST_COMPRESSION_EXCEPTION_ID, "Two segments that are being pinched have different lengths: %i %i\n", stSegment_getLength(segment1), stSegment_getLength(segment2));
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
    if(stBlock_getDegree(block1) > stBlock_getDegree(block2)) { //Avoid merging large blocks into small blocks
        return stBlock_pinch(block2, block1, orientation);
    }
    if(stBlock_getLength(block1) != stBlock_getLength(block2)) {
        stThrowNew(ST_COMPRESSION_EXCEPTION_ID, "Two segments that are being pinched have different lengths: %i %i\n", stBlock_getLength(block1), stBlock_getLength(block2));
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
    return block1;
}

stBlock *stBlock_pinch2(stBlock *block, stSegment *segment, bool orientation) {
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

static stSegment *stSegment_construct(int64_t start) {
    stSegment *segment = st_calloc(1, sizeof(stSegment));
    segment->start = start;
    return segment;
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

stSegment *stThread_get5Prime(stThread *thread, stSegment *segment) {
    return segment->pSegment;
}

stSegment *stThread_get3Prime(stThread *thread, stSegment *segment) {
    return segment->nSegment->nSegment != NULL ? segment->nSegment : NULL;
}

static stSegment *stThread_splitP(stSegment *segment, int64_t leftBlockLength) {
    stSegment *nSegment = segment->nSegment;
    assert(nSegment != NULL);
    stSegment *rightSegment = stSegment_construct(stSegment_getStart(segment) + leftBlockLength);
    segment->nSegment = rightSegment;
    rightSegment->pSegment = segment;
    rightSegment->nSegment = nSegment;
    nSegment->pSegment = rightSegment;
    //Thread
    stSortedSet_insert(thread->segments, rightSegment);
    return rightSegment;
}

void stThread_split(stThread *thread, int64_t leftSideOfSplitPoint) {
    stSegment *segment = stThread_getSegment(thread, leftSideOfSplitPoint);
    if (segment == NULL) {
        return;
    }
    if (leftSideOfSplitPoint == stSegment_getStart(segment)
            + stSegment_getLength(segment) - 1) { //There is already a break
        return;
    }
    int64_t leftBlockLength = leftSideOfSplitPoint - stSegment_getStart(segment) + 1;
    if(stSegment_getBlock(segment) != NULL) {
        stBlockIt blockIt = stBlock_getSegmentIterator(stSegment_getBlock(segment));
        segment = stBlockIt_getNext(&blockIt);
        stSegment *segment2 = stBlockIt_getNext(&blockIt);
        stBlock *rightBlock = stBlock_construct(stThread_splitP(segment, leftBlockLength), stSegment_getBlockOrientation(segment),
                stThread_splitP(segment2, leftBlockLength), stSegment_getBlockOrientation(segment2));
        while ((segment = stBlockIt_getNext(&blockIt)) != NULL) {
            stBlock_pinch2(rightBlock, stThread_splitP(segment, leftBlockLength), stSegment_getBlockOrientation(segment));
        }
    }
    else {
        stThread_splitP(segment, leftBlockLength);
    }
}

void stThread_joinTrivialBoundaries(stThread *thread) {
    stSegment *segment = stThread_getFirst(thread);
    do {
        stBlock *block = stSegment_getBlock(segment);
        while (1) {
            stSegment *nSegment = stThread_get3Prime(thread, segment);
            if (nSegment != NULL) {
                stBlock *nBlock = stSegment_getBlock(nSegment);
                if (block == NULL && nBlock == NULL) {
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
    } while ((segment = stThread_get3Prime(thread, segment)));
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
    stSegment *segment = stSegment_construct(start);
    stSegment *terminatorSegment = stSegment_construct(start + length);
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

void stThreadSet_joinTrivialBoundaries(stThreadSet *threadSet) {
    stThreadIt threadIt = stThreadSet_getIterator(threadSet);
    stThread *thread;
    while ((thread = stThreadIt_getNext(&threadIt)) != NULL) {
        stThread_joinTrivialBoundaries(thread);
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

