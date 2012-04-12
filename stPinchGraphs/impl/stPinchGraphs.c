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
    stSegment *headSegment;
    stSegment *tailSegment;
};

//Blocks

static void connectBlockToSegment(stSegment *segment, bool orientation,
        stBlock *block, stSegment *nSegment, stSegment *pSegment) {
    segment->block = block;
    segment->blockOrientation = orientation;
    segment->nSegment = nSegment;
    segment->pSegment = pSegment;
}

stBlock *stBlock_construct(stSegment *segment1, bool orientation1,
        stSegment *segment2, bool orientation2) {
    stBlock *block = st_malloc(sizeof(stBlock));
    block->headSegment = segment1;
    block->tailSegment = segment2;
    connectBlockToSegment(segment1, orientation1, block, NULL, segment2);
    connectBlockToSegment(segment2, orientation2, block, segment1, NULL);
    return block;
}

stBlock *stBlock_pinch(stBlock *block1, stBlock *block2, bool orientation) {
    stBlockIt blockIt = stBlock_getSegmentIterator(block2);
    stSegment *segment;
    while ((segment = stBlockIt_getNext(blockIt)) != NULL) {
        bool segmentOrientation = stSegment_getBlockOrientation(segment);
        stSegment_removeFromBlock(segment);
        stBlock_pinch2(
                block1,
                segment,
                (segmentOrientation && orientation) || (!segmentOrientation
                        && !orientation));
    }
    return block1;
}

stBlock *stBlock_pinch2(stBlock *block, stSegment *segment, bool orientation) {
    block->tailSegment->nSegment = segment;
    connectBlockToSegment(segment, orientation, block, block->tailSegment, NULL);
    block->tailSegment = segment;
    return block;
}

stBlockIt stBlock_getSegmentIterator(stBlock *block) {
    stBlockIt blockIt;
    blockIt.segment = block->headSegment;
    return blockIt;
}

stSegment *stBlockIt_getNext(stBlockIt blockIt) {
    stSegment *segment = blockIt.segment;
    if (segment != NULL) {
        blockIt.segment = segment->nSegment;
    }
    return segment;
}

//Segments

int64_t stSegment_getStart(stSegment *segment) {
    return segment->start;
}

int64_t stSegment_getLength(stSegment *segment) {
    return segment->nSegment->start - segment->start;
}

stBlock *stSegment_getBlock(stSegment *segment) {
    return segment->block;
}

bool stSegment_getBlockOrientation(stSegment *segment) {
    return segment->blockOrientation;
}

void stSegment_removeFromBlock(stSegment *segment) {
    stBlock *block = stSegment_getBlock(segment);
    assert(block != NULL);
    if (segment->pSegment != NULL) {
        segment->pSegment->nSegment = segment->nSegment;
    } else {
        assert(block->headSegment == segment);
        block->headSegment = segment->nSegment;
    }
    if (segment->nSegment != NULL) {
        segment->nSegment->pSegment = segment->pSegment;
    } else {
        assert(block->tailSegment == segment);
        block->tailSegment = segment->pSegment;
    }
    connectBlockToSegment(segment, 0, NULL, NULL, NULL);
    if (block->headSegment == block->tailSegment) {
        connectBlockToSegment(block->headSegment, 0, NULL, NULL, NULL);
        free(block);
    }
}

//Private segment functions

void stSegment_destruct(stSegment *segment) {
    if (stSegment_getBlock(segment) != NULL) {
        stSegment_removeFromBlock(segment);
    }
    free(segment);
}

int stSegment_compareFunction(const stSegment *segment1,
        const stSegment *segment2) {
    return segment1->start < segment2->start ? -1 : (segment1->start
            > segment2->start ? 1 : 0);
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
    if (stSegment_getStart(segment2) + stSegment_getLength(segment2) <= coordinate) {
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

static stSegment *stSegment_construct(int64_t start) {
    stSegment *segment = st_calloc(1, sizeof(stSegment));
    segment->start = start;
    return segment;
}

void stThread_split(stThread *thread, int64_t leftSideOfSplitPoint) {
    stSegment *segment = stThread_getSegment(thread, leftSideOfSplitPoint);
    if (segment == NULL) {
        return;
    }
    if (leftSideOfSplitPoint == stSegment_getLength(segment)
            + stSegment_getStart(segment) - 1) { //There is already a break
        return;
    }
    stSegment *nSegment = segment->nSegment;
    stSegment *rightSegment = stSegment_construct(leftSideOfSplitPoint + 1);
    segment->nBlockSegment = rightSegment;
    rightSegment->pSegment = segment;
    rightSegment->nSegment = nSegment;
    nSegment->pSegment = rightSegment;
    stSortedSet_insert(thread->segments, rightSegment);
}

void stThread_joinTrivialBoundaries(stThread *thread) {
    stSegment *segment = stThread_getFirst(thread);
    do {
        stBlock *block = stSegment_getBlock(segment);
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
            }
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

stThread *stThreadIt_getNext(stThreadIt threadIt) {
    if (threadIt.index < stThreadSet_getSize(threadIt.threadSet)) {
        stList_get(threadIt.threadSet->threads, threadIt.index++);
    }
    return NULL;
}

void stThreadSet_joinTrivialBoundaries(stThreadSet *threadSet) {
    stThreadIt threadIt = stThreadSet_getIterator(threadSet);
    stThread *thread;
    while ((thread = stThreadIt_getNext(threadIt)) != NULL) {
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

