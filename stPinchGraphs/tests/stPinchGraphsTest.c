/*
 * stPinchGraphsTest.c
 *
 *  Created on: 11 Apr 2012
 *      Author: benedictpaten
 *
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stPinchGraphs.h"

static stThreadSet *threadSet = NULL;
static int64_t name1 = 0, start1 = 1, length1 = ((int64_t) INT32_MAX) + 1;
static int64_t leftSplitPoint1 = 10, leftSplitPoint2 = 11;
static stThread *thread1;
static int64_t name2 = INT32_MAX, start2 = 4, length2 = 10;
static stThread *thread2;

static void teardown() {
    if (threadSet != NULL) {
        stThreadSet_destruct(threadSet);
        threadSet = NULL;
    }
}

static void setup() {
    teardown();
    threadSet = stThreadSet_construct();
    thread1 = stThreadSet_addThread(threadSet, name1, start1, length1);
    thread2 = stThreadSet_addThread(threadSet, name2, start2, length2);
    stThread_split(thread1, leftSplitPoint1);
    stThread_split(thread1, leftSplitPoint2);
}

static void testStThreadSet(CuTest *testCase) {
    setup();
    CuAssertIntEquals(testCase, 2, stThreadSet_getSize(threadSet));
    CuAssertPtrEquals(testCase, thread1, stThreadSet_getThread(threadSet, name1));
    CuAssertPtrEquals(testCase, thread2, stThreadSet_getThread(threadSet, name2));
    CuAssertPtrEquals(testCase, NULL, stThreadSet_getThread(threadSet, 5));
    stThreadIt threadIt = stThreadSet_getIterator(threadSet);
    CuAssertPtrEquals(testCase, thread1, stThreadIt_getNext(&threadIt));
    CuAssertPtrEquals(testCase, thread2, stThreadIt_getNext(&threadIt));
    CuAssertPtrEquals(testCase, NULL, stThreadIt_getNext(&threadIt));
    CuAssertPtrEquals(testCase, NULL, stThreadIt_getNext(&threadIt));

    //Test segment iterators
    stSortedSet *segmentSet = stSortedSet_construct();
    stThreadSetSegmentIt segmentIt = stThreadSet_getSegmentIt(threadSet);
    stSegment *segment;
    while ((segment = stThreadSetSegmentIt_getNext(&segmentIt)) != NULL) {
        CuAssertTrue(testCase, stSortedSet_search(segmentSet, segment) == NULL);
        stSortedSet_insert(segmentSet, segment);
    }
    int32_t segmentCount = 0;
    stThread *threads[] = { thread1, thread2 };
    for (int32_t i = 0; i < 2; i++) {
        stThread *thread = threads[i];
        stSegment *segment = stThread_getFirst(thread);
        while (segment != NULL) {
            segmentCount++;
            CuAssertTrue(testCase, stSortedSet_search(segmentSet, segment) != NULL);
            segment = stSegment_get3Prime(segment);
        }
    }
    CuAssertIntEquals(testCase, segmentCount, stSortedSet_size(segmentSet));

    //Now test block iterator
    stThreadSetBlockIt blockIt = stThreadSet_getBlockIt(threadSet);
    CuAssertPtrEquals(testCase, NULL, stThreadSetBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stThreadSetBlockIt_getNext(&blockIt));
    stBlock *block1 = stBlock_construct2(stThread_getFirst(thread1));
    stBlock *block2 = stBlock_construct2(stThread_getFirst(thread2));
    blockIt = stThreadSet_getBlockIt(threadSet);
    CuAssertPtrEquals(testCase, block1, stThreadSetBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, block2, stThreadSetBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stThreadSetBlockIt_getNext(&blockIt));

    teardown();
}

static void testStThreadAndSegment(CuTest *testCase) {
    setup();
    CuAssertIntEquals(testCase, name1, stThread_getName(thread1));
    CuAssertIntEquals(testCase, start1, stThread_getStart(thread1));
    CuAssertIntEquals(testCase, length1, stThread_getLength(thread1));
    CuAssertIntEquals(testCase, name2, stThread_getName(thread2));
    CuAssertIntEquals(testCase, start2, stThread_getStart(thread2));
    CuAssertIntEquals(testCase, length2, stThread_getLength(thread2));

    //Test the first thread.
    CuAssertPtrEquals(testCase, NULL, stThread_getSegment(thread2, start2 - 1));
    CuAssertPtrEquals(testCase, stThread_getFirst(thread2), stThread_getSegment(thread2, start2));
    stSegment *segment = stThread_getFirst(thread2);
    CuAssertTrue(testCase, segment != NULL);
    CuAssertIntEquals(testCase, start2, stSegment_getStart(segment));
    CuAssertIntEquals(testCase, length2, stSegment_getLength(segment));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment));

    //Test the second thread.
    CuAssertPtrEquals(testCase, NULL, stThread_getSegment(thread1, start1 - 1));
    CuAssertPtrEquals(testCase, stThread_getFirst(thread1), stThread_getSegment(thread1, start1));
    segment = stThread_getFirst(thread1);
    CuAssertTrue(testCase, segment != NULL);
    CuAssertPtrEquals(testCase, NULL, stSegment_get5Prime(segment));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment));
    CuAssertIntEquals(testCase, start1, stSegment_getStart(segment));
    CuAssertIntEquals(testCase, leftSplitPoint1 - start1 + 1, stSegment_getLength(segment));

    stSegment *segment2 = stSegment_get3Prime(segment);
    CuAssertTrue(testCase, segment2 != NULL);
    CuAssertPtrEquals(testCase, segment, stSegment_get5Prime(segment2));
    CuAssertIntEquals(testCase, leftSplitPoint1 + 1, stSegment_getStart(segment2));
    CuAssertIntEquals(testCase, leftSplitPoint2 - leftSplitPoint1, stSegment_getLength(segment2));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment));

    stSegment *segment3 = stSegment_get3Prime(segment2);
    CuAssertTrue(testCase, segment3 != NULL);
    CuAssertPtrEquals(testCase, segment2, stSegment_get5Prime(segment3));
    CuAssertPtrEquals(testCase, NULL, stSegment_get3Prime(segment3));
    CuAssertIntEquals(testCase, leftSplitPoint2 + 1, stSegment_getStart(segment3));
    CuAssertIntEquals(testCase, length1 + start1 - leftSplitPoint2 - 1, stSegment_getLength(segment3));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment));

    //Test stThread_getSegment
    CuAssertPtrEquals(testCase, segment, stThread_getSegment(thread1, start1));
    CuAssertPtrEquals(testCase, segment, stThread_getSegment(thread1, leftSplitPoint1));
    CuAssertPtrEquals(testCase, segment2, stThread_getSegment(thread1, leftSplitPoint1 + 1));
    CuAssertPtrEquals(testCase, segment2, stThread_getSegment(thread1, leftSplitPoint2));
    CuAssertPtrEquals(testCase, segment3, stThread_getSegment(thread1, leftSplitPoint2 + 1));
    CuAssertPtrEquals(testCase, segment3, stThread_getSegment(thread1, start1 + length1 - 1));
    CuAssertPtrEquals(testCase, NULL, stThread_getSegment(thread1, start1 + length1));

    //Test stThread_joinTrivialBoundaries
    stThread_joinTrivialBoundaries(thread1);
    //Now should be just one segment
    segment = stThread_getFirst(thread1);
    CuAssertTrue(testCase, segment != NULL);
    CuAssertPtrEquals(testCase, NULL, stSegment_get5Prime(segment));
    CuAssertIntEquals(testCase, start1, stSegment_getStart(segment));
    CuAssertIntEquals(testCase, length1, stSegment_getLength(segment));
    CuAssertPtrEquals(testCase, NULL, stSegment_get3Prime(segment));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment));
    CuAssertPtrEquals(testCase, segment, stThread_getSegment(thread1, start1));
    CuAssertPtrEquals(testCase, segment, stThread_getSegment(thread1, start1 + length1 - 1));
    CuAssertPtrEquals(testCase, NULL, stThread_getSegment(thread1, start1 - 1));
    CuAssertPtrEquals(testCase, NULL, stThread_getSegment(thread1, start1 + length1));

    teardown();
}

static void testStBlock_NoSplits(CuTest *testCase) {
    setup();
    static int64_t name3 = 5, start3 = 0, length3 = 20;
    stThread *thread3 =
            stThreadSet_addThread(threadSet, name3, start3, length3);
    stThread_split(thread3, 4);
    stThread_split(thread3, 9);
    stThread_split(thread3, 14);
    stSegment *segment1 = stThread_getFirst(thread3);
    stSegment *segment2 = stSegment_get3Prime(segment1);
    stSegment *segment3 = stSegment_get3Prime(segment2);
    stSegment *segment4 = stSegment_get3Prime(segment3);
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stSegment_getLength(segment2));
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stSegment_getLength(segment3));
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stSegment_getLength(segment4));
    CuAssertPtrEquals(testCase, NULL, stSegment_get3Prime(segment4));
    stBlock *block = stBlock_construct(segment1, 1, segment2, 0);
    CuAssertIntEquals(testCase, 2, stBlock_getDegree(block));
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stBlock_getLength(block));
    //get block
    CuAssertPtrEquals(testCase, segment1, stBlock_getFirst(block));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment2));
    //get block orientation
    CuAssertIntEquals(testCase, 1, stSegment_getBlockOrientation(segment1));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment2));
    //test iterator
    stBlockIt blockIt = stBlock_getSegmentIterator(block);
    CuAssertPtrEquals(testCase, segment1, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment2, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stBlockIt_getNext(&blockIt));

    //Now try pinching in a segment
    stBlock_pinch2(block, segment3, 1);
    CuAssertPtrEquals(testCase, segment1, stBlock_getFirst(block));
    CuAssertIntEquals(testCase, 3, stBlock_getDegree(block));
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stBlock_getLength(block));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment3));
    CuAssertIntEquals(testCase, 1, stSegment_getBlockOrientation(segment3));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment2));
    CuAssertIntEquals(testCase, 1, stSegment_getBlockOrientation(segment1));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment2));
    blockIt = stBlock_getSegmentIterator(block);
    CuAssertPtrEquals(testCase, segment1, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment2, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment3, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stBlockIt_getNext(&blockIt));

    //Now try removing from the block
    stBlock_destruct(stSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment2));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment3));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment1));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment2));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment3));

    //Now try merging two blocks and undoing them
    block = stBlock_pinch(
            stBlock_pinch(stBlock_construct2(segment1),
                    stBlock_construct2(segment2), 0),
            stBlock_construct(segment3, 0, segment4, 1), 0);
    CuAssertIntEquals(testCase, 4, stBlock_getDegree(block));
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stBlock_getLength(block));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment1));
    CuAssertIntEquals(testCase, 1, stSegment_getBlockOrientation(segment1));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment2));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment2));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment3));
    CuAssertIntEquals(testCase, 1, stSegment_getBlockOrientation(segment3));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment4));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment4));

    CuAssertPtrEquals(testCase, segment1, stBlock_getFirst(block));
    blockIt = stBlock_getSegmentIterator(block);
    CuAssertPtrEquals(testCase, segment1, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment2, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment3, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment4, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stBlockIt_getNext(&blockIt));

    //Test block destruct
    stBlock_destruct(block);
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment1));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment1));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment2));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment2));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment3));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment3));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment4));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment4));
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), 5);
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stSegment_getLength(segment2));
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stSegment_getLength(segment3));
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stSegment_getLength(segment4));

    //Now try merging two blocks of uneven length
    stThread_split(thread3, 0);
    segment1 = stThread_getFirst(thread3);
    segment2 = stSegment_get3Prime(segment1);
    stTry {
            stBlock_construct(segment1, 1, segment2, 1);
            CuAssertTrue(testCase, 0);
        }stCatch(ST_PINCH_GRAPH_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(ST_PINCH_GRAPH_EXCEPTION_ID));
            }stTryEnd

    //Now make a block with a single element
    block = stBlock_construct2(segment1);
    CuAssertIntEquals(testCase, 1, stBlock_getDegree(block));
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stBlock_getLength(block));
    //get block
    CuAssertPtrEquals(testCase, segment1, stBlock_getFirst(block));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment1));
    //get block orientation
    CuAssertIntEquals(testCase, 1, stSegment_getBlockOrientation(segment1));
    //test iterator
    blockIt = stBlock_getSegmentIterator(block);
    CuAssertPtrEquals(testCase, segment1, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stBlockIt_getNext(&blockIt));
    stBlock_destruct(block);
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment1));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment1));

    teardown();

}

static void testStBlock_Splits(CuTest *testCase) {
    /*
     * Tests splitting of segments that are aligned. Put in seperate function as draws together previous code.
     */
    setup();
    static int64_t name3 = 5, start3 = 0, length3 = 15;
    stThread *thread3 =
            stThreadSet_addThread(threadSet, name3, start3, length3);
    stThread_split(thread3, 4);
    stThread_split(thread3, 9);
    stSegment *segment1 = stThread_getFirst(thread3);
    stSegment *segment2 = stSegment_get3Prime(segment1);
    stSegment *segment3 = stSegment_get3Prime(segment2);
    CuAssertPtrEquals(testCase, NULL, stSegment_get3Prime(segment3));
    stBlock_pinch2(stBlock_construct(segment1, 1, segment2, 0), segment3, 0);

    stThread_split(thread3, 0);

    //Now traverse through thread and check all is okay
    segment1 = stThread_getFirst(thread3);
    CuAssertTrue(testCase, segment1 != NULL);
    CuAssertIntEquals(testCase, 1, stSegment_getLength(segment1));
    segment2 = stSegment_get3Prime(segment1);
    CuAssertTrue(testCase, segment2 != NULL);
    CuAssertIntEquals(testCase, 4, stSegment_getLength(segment2));
    segment3 = stSegment_get3Prime(segment2);
    CuAssertTrue(testCase, segment3 != NULL);
    CuAssertIntEquals(testCase, 1, stSegment_getLength(segment3));
    stSegment *segment4 = stSegment_get3Prime(segment3);
    CuAssertTrue(testCase, segment4 != NULL);
    CuAssertIntEquals(testCase, 4, stSegment_getLength(segment4));
    stSegment *segment5 = stSegment_get3Prime(segment4);
    CuAssertTrue(testCase, segment5 != NULL);
    CuAssertIntEquals(testCase, 1, stSegment_getLength(segment5));
    stSegment *segment6 = stSegment_get3Prime(segment5);
    CuAssertTrue(testCase, segment6 != NULL);
    CuAssertIntEquals(testCase, 4, stSegment_getLength(segment6));
    CuAssertPtrEquals(testCase, NULL, stSegment_get3Prime(segment6));
    //Check the thread.
    CuAssertPtrEquals(testCase, thread3, stSegment_getThread(segment1));
    CuAssertPtrEquals(testCase, thread3, stSegment_getThread(segment2));
    CuAssertPtrEquals(testCase, thread3, stSegment_getThread(segment3));
    CuAssertPtrEquals(testCase, thread3, stSegment_getThread(segment4));
    CuAssertPtrEquals(testCase, thread3, stSegment_getThread(segment5));
    CuAssertPtrEquals(testCase, thread3, stSegment_getThread(segment6));
    //Check the name
    CuAssertIntEquals(testCase, name3, stSegment_getName(segment1));
    CuAssertIntEquals(testCase, name3, stSegment_getName(segment2));
    CuAssertIntEquals(testCase, name3, stSegment_getName(segment3));
    CuAssertIntEquals(testCase, name3, stSegment_getName(segment4));
    CuAssertIntEquals(testCase, name3, stSegment_getName(segment5));
    CuAssertIntEquals(testCase, name3, stSegment_getName(segment6));

    //Now check blocks
    stBlock *block1 = stSegment_getBlock(segment1);
    CuAssertTrue(testCase, block1 != NULL);
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stBlock_getLength(block1));
    CuAssertIntEquals(testCase, 3, stBlock_getDegree(block1));
    stBlock *block2 = stSegment_getBlock(segment2);
    CuAssertTrue(testCase, block2 != NULL);
    CuAssertIntEquals(testCase, stSegment_getLength(segment2), stBlock_getLength(block2));
    CuAssertIntEquals(testCase, 3, stBlock_getDegree(block2));

    stBlockIt blockIt = stBlock_getSegmentIterator(block1);
    CuAssertPtrEquals(testCase, segment1, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment3, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment5, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, block1, stSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, block1, stSegment_getBlock(segment3));
    CuAssertPtrEquals(testCase, block1, stSegment_getBlock(segment5));

    blockIt = stBlock_getSegmentIterator(block2);
    CuAssertPtrEquals(testCase, segment2, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment4, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment6, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, block2, stSegment_getBlock(segment2));
    CuAssertPtrEquals(testCase, block2, stSegment_getBlock(segment4));
    CuAssertPtrEquals(testCase, block2, stSegment_getBlock(segment6));

    //Now do merge
    stThreadSet_joinTrivialBoundaries(threadSet);
    stThreadSetBlockIt threadBlockIt = stThreadSet_getBlockIt(threadSet);
    block1 = stThreadSetBlockIt_getNext(&threadBlockIt);
    CuAssertTrue(testCase, block1 != NULL);
    CuAssertPtrEquals(testCase, NULL, stThreadSetBlockIt_getNext(&threadBlockIt));
    blockIt = stBlock_getSegmentIterator(block1);
    CuAssertPtrEquals(testCase, segment1, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment3, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment5, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, block1, stSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, block1, stSegment_getBlock(segment3));
    CuAssertPtrEquals(testCase, block1, stSegment_getBlock(segment5));
    CuAssertIntEquals(testCase, 5, stBlock_getLength(block1));

    teardown();
}

static bool areAligned(stThread *thread1, int32_t base1, stThread *thread2,
        int32_t base2, bool orientation) {
    stSegment *segment1 = stThread_getSegment(thread1, base1);
    stSegment *segment2 = stThread_getSegment(thread2, base2);
    assert(segment1 != NULL && segment2 != NULL);
    stBlock *block1 = stSegment_getBlock(segment1);
    stBlock *block2 = stSegment_getBlock(segment2);
    if (block1 == NULL) {
        return 0;
    }
    if (block1 != block2) {
        return 0;
    }
    int32_t offset1 = base1 - stSegment_getStart(segment1);
    int32_t offset2 = base2 - stSegment_getStart(segment2);
    assert(base1 >= 0 && base2 >= 0);
    if (orientation) {
        return offset1 == offset2;
    }
    return stBlock_getLength(block1) - offset2 == offset1;
}

static void testStThread_pinchP(CuTest *testCase, int32_t segmentNumber,
        int64_t start, int64_t lengths[], int64_t blockDegrees[],
        stThread *thread) {
    stSegment *segment = stThread_getFirst(thread);
    int64_t i = start;
    for (int32_t j = 0; j < segmentNumber; j++) {
        CuAssertTrue(testCase, segment != NULL);
        CuAssertIntEquals(testCase, stThread_getName(thread), stSegment_getName(segment));
        CuAssertIntEquals(testCase, i, stSegment_getStart(segment));
        CuAssertIntEquals(testCase, lengths[j], stSegment_getLength(segment));
        //Block
        if (blockDegrees[j] == 1) {
            CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment));
        } else {
            stBlock *block = stSegment_getBlock(segment);
            CuAssertTrue(testCase, block != NULL);
            CuAssertIntEquals(testCase, blockDegrees[j], stBlock_getDegree(block));
        }

        i += stSegment_getLength(segment);
        if (j + 1 == segmentNumber) {
            CuAssertTrue(testCase, stSegment_get3Prime(segment) == NULL);
        } else {
            stSegment *segment2 = stSegment_get3Prime(segment);
            CuAssertTrue(testCase, segment2 != NULL);
            CuAssertPtrEquals(testCase, segment, stSegment_get5Prime(segment2));
            segment = segment2;
        }
    }
}

static void testStThread_pinch(CuTest *testCase) {
    setup();
    stThread_pinch(thread1, thread2, 5, 5, 8, 1);
    int64_t lengths1[] = { 4, 6, 1, 1, length1 - 12 };
    int64_t blockDegrees1[] = { 1, 2, 2, 2, 1 };
    testStThread_pinchP(testCase, 5, start1, lengths1, blockDegrees1, thread1);
    int64_t lengths2[] = { 1, 5, 1, 1, 1, 1 };
    int64_t blockDegrees2[] = { 1, 2, 2, 2, 2, 1 };
    testStThread_pinchP(testCase, 6, start2, lengths2, blockDegrees2, thread2);

    stThread_pinch(thread1, thread2, 4, 10, 4, 0);
    int64_t lengths1b[] = { 3, 1, 1, 1, 1, 2, 1, 1, 1, length1 - 12 };
    int64_t blockDegrees1b[] = { 1, 2, 4, 4, 4, 2, 4, 4, 4, 1 };
    testStThread_pinchP(testCase, 10, start1, lengths1b, blockDegrees1b,
            thread1);
    int64_t lengths2b[] = { 1, 5, 1, 1, 1, 1 };
    int64_t blockDegrees2b[] = { 1, 2, 4, 4, 4, 2 };
    testStThread_pinchP(testCase, 6, start2, lengths2b, blockDegrees2b, thread2);

    stThread_pinch(thread1, thread2, 4, 10, 4, 0); //Doing the same thing again should not affect the result
    testStThread_pinchP(testCase, 10, start1, lengths1b, blockDegrees1b,
            thread1);
    testStThread_pinchP(testCase, 6, start2, lengths2b, blockDegrees2b, thread2);
    //nor should a zero length pinch
    stThread_pinch(thread1, thread2, 5, 10, 0, 0);
    testStThread_pinchP(testCase, 10, start1, lengths1b, blockDegrees1b,
            thread1);
    testStThread_pinchP(testCase, 6, start2, lengths2b, blockDegrees2b, thread2);

    //Check a subset of the homology groups
    CuAssertTrue(testCase, areAligned(thread1, 4, thread2, 12, 0));

    CuAssertTrue(testCase, areAligned(thread1, 5, thread2, 5, 1));
    CuAssertTrue(testCase, areAligned(thread1, 5, thread2, 12, 0));
    CuAssertTrue(testCase, areAligned(thread1, 5, thread1, 11, 0));

    CuAssertTrue(testCase, areAligned(thread1, 6, thread2, 6, 1));
    CuAssertTrue(testCase, areAligned(thread1, 6, thread2, 11, 0));
    CuAssertTrue(testCase, areAligned(thread1, 6, thread1, 10, 0));

    CuAssertTrue(testCase, areAligned(thread1, 7, thread2, 7, 1));
    CuAssertTrue(testCase, areAligned(thread1, 7, thread2, 10, 0));
    CuAssertTrue(testCase, areAligned(thread1, 7, thread1, 10, 0));

    CuAssertTrue(testCase, areAligned(thread1, 8, thread2, 8, 1));

    CuAssertTrue(testCase, areAligned(thread1, 9, thread2, 10, 1));

    teardown();
}

//Functions that implement a very simple merging of alignment positions

static void addColumn(stHash *columns, int64_t name, int64_t position,
        int64_t strand) {
    assert(position > 0);
    stInt64Tuple *p = stInt64Tuple_construct(2, name,
            strand ? position : -position);
    stSortedSet *column = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stInt64Tuple_cmpFn, NULL);
    stSortedSet_insert(column, p);
    stHash_insert(columns, p, column);
}

static stSortedSet *getColumn(stHash *columns, int64_t name, int64_t position,
        int64_t strand) {
    stInt64Tuple *p = stInt64Tuple_construct(2, name, position);
    stSortedSet *column = stHash_search(columns, p);
    assert(column != NULL);
    stInt64Tuple_destruct(p);
    return column;
}

static void mergePositions(stHash *columns, int64_t name1, int64_t start1,
        bool strand1, int64_t name2, int64_t start2, bool strand2) {
    stSortedSet *column1 = getColumn(columns, name1, start1, strand1);
    stSortedSet *column2 = getColumn(columns, name2, start2, strand2);
    stSortedSet *mergedColumn = stSortedSet_getUnion(column1, column2);
    stSortedSetIterator *it = stSortedSet_getIterator(mergedColumn);
    stInt64Tuple *position;
    while ((position = stSortedSet_getNext(it)) != NULL) {
        assert(stHash_remove(columns, position) != NULL);
        stHash_insert(columns, position, mergedColumn);
    }
    stSortedSet_destructIterator(it);
    stSortedSet_destruct(column1);
    stSortedSet_destruct(column2);
}

static void mergePositionsSymmetric(stHash *columns, int64_t name1,
        int64_t start1, bool strand1, int64_t name2, int64_t start2,
        bool strand2) {
    mergePositions(columns, name1, start1, strand1, name2, start2, strand2);
    mergePositions(columns, name1, start1, !strand1, name2, start2, !strand2);
}

stHash *getUnalignedColumns(stThreadSet *threadSet) {
    stHash *columns = stHash_construct3((uint32_t (*)(const void *))stInt64Tuple_hashKey,
            (int (*)(const void *, const void *))stInt64Tuple_equalsFn, (void (*)(void *))stInt64Tuple_destruct, NULL);
    stThreadIt it = stThreadSet_getIterator(threadSet);
    stThread *thread;
    while ((thread = stThreadIt_getNext(&it))) {
        for (int32_t i = 0; i < stThread_getLength(thread); i++) {
            addColumn(columns, stThread_getName(thread),
                    stThread_getStart(thread) + i, 1);
            addColumn(columns, stThread_getName(thread),
                    stThread_getStart(thread) + i, 0);
        }
    }
    return columns;
}

static void decodePosition(stThreadSet *threadSet,
        stInt64Tuple *alignedPosition, stThread **thread, int64_t *position,
        bool *strand) {
    *thread = stThreadSet_getThread(threadSet, stInt64Tuple_getPosition(alignedPosition, 0));
    *position = stInt64Tuple_getPosition(alignedPosition, 1);
    assert(*position != 0);
    *strand = 1;
    if(*position < 0) {
        *strand = 0;
        *position *= -1;
    }
}

static void randomPosition(stList *threadList, stThread **thread,
        int64_t *position, bool *strand) {
    *thread = st_randomChoice(threadList);
    *position = st_randomInt(stThread_getStart(*thread), stThread_getStart(*thread) + stThread_getLength(*thread));
    *position = st_random() > 0.5;
}

static void randomPinch(stList *threadList, stThread **thread1,
        int64_t *start1, bool *strand1, stThread **thread2, int64_t *start2,
        bool *strand2, int64_t *length) {
    randomPosition(threadList, thread1, start1, strand1);
    randomPosition(threadList, thread2, start2, strand2);
    int32_t i = stThread_getStart(*thread1) + stThread_getLength(*thread1)
            - *start1;
    int32_t j = stThread_getStart(*thread2) + stThread_getLength(*thread2)
            - *start2;
    *length = st_randomInt(0, i > j ? j : i);
}

static stThreadSet *getRandomThreadSet() {
    stThreadSet *threadSet = stThreadSet_construct();
    int32_t randomThreadNumber = st_randomInt(2, 10);
    for (int32_t threadIndex = 0; threadIndex < randomThreadNumber; threadIndex++) {
        int32_t start = st_randomInt(1, 100);
        int32_t length = st_randomInt(0, 100);
        int32_t threadName = threadIndex + 4;
        stThreadSet_addThread(threadSet, threadName, start, length);
    }
    return threadSet;
}

static stList *getThreadList(stThreadSet *threadSet) {
    stList *threadList = stList_construct();
    stThreadIt threadIt = stThreadSet_getIterator(threadSet);
    stThread *thread;
    while((thread = stThreadIt_getNext(&threadIt)) != NULL) {
        stList_append(threadList, thread);
    }
    return threadList;
}

static stThreadSet *getRandomPinchGraph() {
    stThreadSet *threadSet = getRandomThreadSet();
    stList *threadList = getThreadList(threadSet);
    //Randomly push them together, updating both sets, and checking that set of alignments is what we expect
    while (st_random() > 0.01) {
        stThread *thread1, *thread2;
        int64_t start1, start2, length;
        bool strand1, strand2;
        randomPinch(threadList, &thread1, &start1, &strand1, &thread2, &start1,
                &strand2, &length);
        stThread_pinch(thread1, thread2, start1, start2, length,
                strand1 == strand2);
    }
    stList_destruct(threadList);
    return threadSet;
}

static void testStThread_pinch_randomTests(CuTest *testCase) {
    stThreadSet *threadSet = getRandomThreadSet();
    stList *threadList = getThreadList(threadSet);
    stHash *columns = getUnalignedColumns(threadSet);

    //Randomly push them together, updating both sets, and checking that set of alignments is what we expect
    while (st_random() > 0.01) {
        stThread *thread1, *thread2;
        int64_t start1, start2, length;
        bool strand1, strand2;
        randomPinch(threadList, &thread1, &start1, &strand1, &thread2, &start1,
                &strand2, &length);
        stThread_pinch(thread1, thread2, start1, start2, length,
                strand1 == strand2);
        //now do all the pushing together of the equivalence classes
        for (int32_t i = 0; i < length; i++) {
            mergePositionsSymmetric(columns, stThread_getName(thread1),
                    start1 + i, strand1, stThread_getName(thread2), start2 + i,
                    strand2);
        }
    }
    stList_destruct(threadList);

    //Check they are equivalent
    stHashIterator *hashIt = stHash_getIterator(columns);
    stInt64Tuple *key;
    while ((key = stHash_getNext(hashIt)) != NULL) {
        stThread *thread1, *thread2;
        int64_t position1, position2;
        bool strand1, strand2;
        stSortedSet *column = stHash_search(columns, key);
        stList *columnList = stSortedSet_getList(column);
        for (int32_t i = 0; i < stList_length(columnList); i++) {
            stInt64Tuple *alignedPosition1 = stList_get(columnList, i);
            decodePosition(threadSet, alignedPosition1, &thread1, &position1, &strand1);
            for (int32_t j = i + 1; j < stList_length(columnList); j++) {
                stInt64Tuple *alignedPosition2 = stList_get(columnList, j);
                decodePosition(threadSet, alignedPosition2, &thread2, &position2, &strand2);
                CuAssertTrue(testCase, areAligned(thread1, position1, thread2,
                                position2, strand1 == strand2));
            }
        }
    }
}

static void checkNotTrivialBoundary(CuTest *testCase, stSegment *segment2,
        stSegment *segment) {

}

static void testStThreadSet_joinTrivialBoundaries_randomTests(CuTest *testCase) {
    for (int32_t test = 0; test < 100; test++) {
        stThreadSet *threadSet = getRandomPinchGraph();
        stThreadSet_joinTrivialBoundaries(threadSet);
        stThreadSetSegmentIt segmentIt = stThreadSet_getSegmentIt(threadSet);
        stSegment *segment;
        while ((segment = stThreadSetSegmentIt_getNext(&segmentIt))) {
            //Check no trivial boundaries
            stSegment *segment2 = stSegment_get5Prime(segment);
            if (segment2 != NULL) {
                checkNotTrivialBoundary(testCase, segment2, segment);
            }
            segment2 = stSegment_get3Prime(segment);
            if (segment2 != NULL) {
                checkNotTrivialBoundary(testCase, segment, segment2);
            }
        }
        stThreadSet_destruct(threadSet);
    }
}

static void testStThreadSet_getAdjacencyComponents(CuTest *testCase) {
    setup();

    teardown();
}

static void testStThreadSet_getAdjacencyComponents_randomTests(CuTest *testCase) {
    for (int32_t test = 0; test < 100; test++) {
        stThreadSet *threadSet = getRandomPinchGraph();
        //stList *adjacencyComponents = stThreadSet_getAdjacencyComponents(
        //       threadSet);
        //Check all ends in one adjacency component

        //Check all connected nodes in same adjacency component
        stThreadSet_destruct(threadSet);
    }
}

CuSuite* stPinchGraphsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testStThreadSet);
    SUITE_ADD_TEST(suite, testStThreadAndSegment);
    SUITE_ADD_TEST(suite, testStBlock_NoSplits);
    SUITE_ADD_TEST(suite, testStBlock_Splits);
    SUITE_ADD_TEST(suite, testStThreadSet_getAdjacencyComponents);
    SUITE_ADD_TEST(suite, testStThread_pinch);
    SUITE_ADD_TEST(suite, testStThreadSet_joinTrivialBoundaries_randomTests);
    SUITE_ADD_TEST(suite, testStThreadSet_getAdjacencyComponents_randomTests);
    SUITE_ADD_TEST(suite, testStThread_pinch_randomTests);
    return suite;
}
