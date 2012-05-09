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

static stPinchThreadSet *threadSet = NULL;
static int64_t name1 = 0, start1 = 1, length1 = ((int64_t) INT32_MAX) + 1;
static int64_t leftSplitPoint1 = 10, leftSplitPoint2 = 11;
static stPinchThread *thread1;
static int64_t name2 = INT32_MAX, start2 = 4, length2 = 10;
static stPinchThread *thread2;

static void teardown() {
    if (threadSet != NULL) {
        stPinchThreadSet_destruct(threadSet);
        threadSet = NULL;
    }
}

static void setup() {
    teardown();
    threadSet = stPinchThreadSet_construct();
    thread1 = stPinchThreadSet_addThread(threadSet, name1, start1, length1);
    thread2 = stPinchThreadSet_addThread(threadSet, name2, start2, length2);
    stPinchThread_split(thread1, leftSplitPoint1);
    stPinchThread_split(thread1, leftSplitPoint2);
}

static void testStPinchThreadSet(CuTest *testCase) {
    setup();
    CuAssertIntEquals(testCase, 2, stPinchThreadSet_getSize(threadSet));
    CuAssertPtrEquals(testCase, thread1, stPinchThreadSet_getThread(threadSet, name1));
    CuAssertPtrEquals(testCase, thread2, stPinchThreadSet_getThread(threadSet, name2));
    CuAssertPtrEquals(testCase, NULL, stPinchThreadSet_getThread(threadSet, 5));
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    CuAssertPtrEquals(testCase, thread1, stPinchThreadSetIt_getNext(&threadIt));
    CuAssertPtrEquals(testCase, thread2, stPinchThreadSetIt_getNext(&threadIt));
    CuAssertPtrEquals(testCase, NULL, stPinchThreadSetIt_getNext(&threadIt));
    CuAssertPtrEquals(testCase, NULL, stPinchThreadSetIt_getNext(&threadIt));

    //Test segment iterators
    stSortedSet *segmentSet = stSortedSet_construct();
    stPinchThreadSetSegmentIt segmentIt = stPinchThreadSet_getSegmentIt(threadSet);
    stPinchSegment *segment;
    while ((segment = stPinchThreadSetSegmentIt_getNext(&segmentIt)) != NULL) {
        CuAssertTrue(testCase, stSortedSet_search(segmentSet, segment) == NULL);
        stSortedSet_insert(segmentSet, segment);
    }
    int32_t segmentCount = 0;
    stPinchThread *threads[] = { thread1, thread2 };
    for (int32_t i = 0; i < 2; i++) {
        stPinchThread *thread = threads[i];
        stPinchSegment *segment = stPinchThread_getFirst(thread);
        while (segment != NULL) {
            segmentCount++;
            CuAssertTrue(testCase, stSortedSet_search(segmentSet, segment) != NULL);
            segment = stPinchSegment_get3Prime(segment);
        }
    }
    CuAssertIntEquals(testCase, segmentCount, stSortedSet_size(segmentSet));
    stSortedSet_destruct(segmentSet);

    //Now test block iterator
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    CuAssertPtrEquals(testCase, NULL, stPinchThreadSetBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stPinchThreadSetBlockIt_getNext(&blockIt));
    stPinchBlock *block1 = stPinchBlock_construct2(stPinchThread_getFirst(thread1));
    stPinchBlock *block2 = stPinchBlock_construct2(stPinchThread_getFirst(thread2));
    blockIt = stPinchThreadSet_getBlockIt(threadSet);
    CuAssertPtrEquals(testCase, block1, stPinchThreadSetBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, block2, stPinchThreadSetBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stPinchThreadSetBlockIt_getNext(&blockIt));

    teardown();
}

static void testStPinchThreadAndSegment(CuTest *testCase) {
    setup();
    CuAssertIntEquals(testCase, name1, stPinchThread_getName(thread1));
    CuAssertIntEquals(testCase, start1, stPinchThread_getStart(thread1));
    CuAssertIntEquals(testCase, length1, stPinchThread_getLength(thread1));
    CuAssertIntEquals(testCase, name2, stPinchThread_getName(thread2));
    CuAssertIntEquals(testCase, start2, stPinchThread_getStart(thread2));
    CuAssertIntEquals(testCase, length2, stPinchThread_getLength(thread2));

    //Test the first thread.
    CuAssertPtrEquals(testCase, NULL, stPinchThread_getSegment(thread2, start2 - 1));
    CuAssertPtrEquals(testCase, stPinchThread_getFirst(thread2), stPinchThread_getSegment(thread2, start2));
    stPinchSegment *segment = stPinchThread_getFirst(thread2);
    CuAssertTrue(testCase, segment != NULL);
    CuAssertIntEquals(testCase, start2, stPinchSegment_getStart(segment));
    CuAssertIntEquals(testCase, length2, stPinchSegment_getLength(segment));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment));

    //Test the second thread.
    CuAssertPtrEquals(testCase, NULL, stPinchThread_getSegment(thread1, start1 - 1));
    CuAssertPtrEquals(testCase, stPinchThread_getFirst(thread1), stPinchThread_getSegment(thread1, start1));
    segment = stPinchThread_getFirst(thread1);
    CuAssertTrue(testCase, segment != NULL);
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_get5Prime(segment));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment));
    CuAssertIntEquals(testCase, start1, stPinchSegment_getStart(segment));
    CuAssertIntEquals(testCase, leftSplitPoint1 - start1 + 1, stPinchSegment_getLength(segment));

    stPinchSegment *segment2 = stPinchSegment_get3Prime(segment);
    CuAssertTrue(testCase, segment2 != NULL);
    CuAssertPtrEquals(testCase, segment, stPinchSegment_get5Prime(segment2));
    CuAssertIntEquals(testCase, leftSplitPoint1 + 1, stPinchSegment_getStart(segment2));
    CuAssertIntEquals(testCase, leftSplitPoint2 - leftSplitPoint1, stPinchSegment_getLength(segment2));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment));

    stPinchSegment *segment3 = stPinchSegment_get3Prime(segment2);
    CuAssertTrue(testCase, segment3 != NULL);
    CuAssertPtrEquals(testCase, segment2, stPinchSegment_get5Prime(segment3));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_get3Prime(segment3));
    CuAssertIntEquals(testCase, leftSplitPoint2 + 1, stPinchSegment_getStart(segment3));
    CuAssertIntEquals(testCase, length1 + start1 - leftSplitPoint2 - 1, stPinchSegment_getLength(segment3));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment));

    //Test stPinchThread_getSegment
    CuAssertPtrEquals(testCase, segment, stPinchThread_getSegment(thread1, start1));
    CuAssertPtrEquals(testCase, segment, stPinchThread_getSegment(thread1, leftSplitPoint1));
    CuAssertPtrEquals(testCase, segment2, stPinchThread_getSegment(thread1, leftSplitPoint1 + 1));
    CuAssertPtrEquals(testCase, segment2, stPinchThread_getSegment(thread1, leftSplitPoint2));
    CuAssertPtrEquals(testCase, segment3, stPinchThread_getSegment(thread1, leftSplitPoint2 + 1));
    CuAssertPtrEquals(testCase, segment3, stPinchThread_getSegment(thread1, start1 + length1 - 1));
    CuAssertPtrEquals(testCase, NULL, stPinchThread_getSegment(thread1, start1 + length1));

    //Test stPinchThread_joinTrivialBoundaries
    stPinchThread_joinTrivialBoundaries(thread1);
    //Now should be just one segment
    segment = stPinchThread_getFirst(thread1);
    CuAssertTrue(testCase, segment != NULL);
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_get5Prime(segment));
    CuAssertIntEquals(testCase, start1, stPinchSegment_getStart(segment));
    CuAssertIntEquals(testCase, length1, stPinchSegment_getLength(segment));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_get3Prime(segment));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment));
    CuAssertPtrEquals(testCase, segment, stPinchThread_getSegment(thread1, start1));
    CuAssertPtrEquals(testCase, segment, stPinchThread_getSegment(thread1, start1 + length1 - 1));
    CuAssertPtrEquals(testCase, NULL, stPinchThread_getSegment(thread1, start1 - 1));
    CuAssertPtrEquals(testCase, NULL, stPinchThread_getSegment(thread1, start1 + length1));

    //Test stPinchThreadSet_getTotalBlockNumber
    CuAssertIntEquals(testCase, 0, stPinchThreadSet_getTotalBlockNumber(threadSet));

    teardown();
}

static void testStPinchBlock_NoSplits(CuTest *testCase) {
    setup();
    static int64_t name3 = 5, start3 = 0, length3 = 20;
    stPinchThread *thread3 = stPinchThreadSet_addThread(threadSet, name3, start3, length3);
    stPinchThread_split(thread3, 4);
    stPinchThread_split(thread3, 9);
    stPinchThread_split(thread3, 14);
    stPinchSegment *segment1 = stPinchThread_getFirst(thread3);
    stPinchSegment *segment2 = stPinchSegment_get3Prime(segment1);
    stPinchSegment *segment3 = stPinchSegment_get3Prime(segment2);
    stPinchSegment *segment4 = stPinchSegment_get3Prime(segment3);
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchSegment_getLength(segment2));
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchSegment_getLength(segment3));
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchSegment_getLength(segment4));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_get3Prime(segment4));
    stPinchBlock *block = stPinchBlock_construct(segment1, 1, segment2, 0);
    CuAssertIntEquals(testCase, 2, stPinchBlock_getDegree(block));
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchBlock_getLength(block));
    //get block
    CuAssertPtrEquals(testCase, segment1, stPinchBlock_getFirst(block));
    CuAssertPtrEquals(testCase, block, stPinchSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, block, stPinchSegment_getBlock(segment2));
    //get block orientation
    CuAssertIntEquals(testCase, 1, stPinchSegment_getBlockOrientation(segment1));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment2));
    //test iterator
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    CuAssertPtrEquals(testCase, segment1, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment2, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stPinchBlockIt_getNext(&blockIt));

    //Now try pinching in a segment
    stPinchBlock_pinch2(block, segment3, 1);
    CuAssertPtrEquals(testCase, segment1, stPinchBlock_getFirst(block));
    CuAssertIntEquals(testCase, 3, stPinchBlock_getDegree(block));
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchBlock_getLength(block));
    CuAssertPtrEquals(testCase, block, stPinchSegment_getBlock(segment3));
    CuAssertIntEquals(testCase, 1, stPinchSegment_getBlockOrientation(segment3));
    CuAssertPtrEquals(testCase, block, stPinchSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, block, stPinchSegment_getBlock(segment2));
    CuAssertIntEquals(testCase, 1, stPinchSegment_getBlockOrientation(segment1));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment2));
    blockIt = stPinchBlock_getSegmentIterator(block);
    CuAssertPtrEquals(testCase, segment1, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment2, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment3, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stPinchBlockIt_getNext(&blockIt));

    //Now try removing from the block
    stPinchBlock_destruct(stPinchSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment2));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment3));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment1));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment2));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment3));

    //Now try merging two blocks and undoing them
    block = stPinchBlock_pinch(stPinchBlock_pinch(stPinchBlock_construct2(segment1), stPinchBlock_construct2(segment2), 0),
            stPinchBlock_construct(segment3, 0, segment4, 1), 0);
    CuAssertIntEquals(testCase, 4, stPinchBlock_getDegree(block));
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchBlock_getLength(block));
    CuAssertPtrEquals(testCase, block, stPinchSegment_getBlock(segment1));
    CuAssertIntEquals(testCase, 1, stPinchSegment_getBlockOrientation(segment1));
    CuAssertPtrEquals(testCase, block, stPinchSegment_getBlock(segment2));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment2));
    CuAssertPtrEquals(testCase, block, stPinchSegment_getBlock(segment3));
    CuAssertIntEquals(testCase, 1, stPinchSegment_getBlockOrientation(segment3));
    CuAssertPtrEquals(testCase, block, stPinchSegment_getBlock(segment4));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment4));

    CuAssertPtrEquals(testCase, segment1, stPinchBlock_getFirst(block));
    blockIt = stPinchBlock_getSegmentIterator(block);
    CuAssertPtrEquals(testCase, segment1, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment2, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment3, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment4, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stPinchBlockIt_getNext(&blockIt));

    //Test block destruct
    stPinchBlock_destruct(block);
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment1));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment1));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment2));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment2));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment3));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment3));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment4));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment4));
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), 5);
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchSegment_getLength(segment2));
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchSegment_getLength(segment3));
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchSegment_getLength(segment4));

    //Now try merging two blocks of uneven length
    stPinchThread_split(thread3, 0);
    segment1 = stPinchThread_getFirst(thread3);
    segment2 = stPinchSegment_get3Prime(segment1);
    stTry {
            stPinchBlock_construct(segment1, 1, segment2, 1);
            CuAssertTrue(testCase, 0);
        }stCatch(ST_PINCH_GRAPH_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(ST_PINCH_GRAPH_EXCEPTION_ID));
            }stTryEnd

    //Now make a block with a single element
    block = stPinchBlock_construct2(segment1);
    CuAssertIntEquals(testCase, 1, stPinchBlock_getDegree(block));
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchBlock_getLength(block));
    //get block
    CuAssertPtrEquals(testCase, segment1, stPinchBlock_getFirst(block));
    CuAssertPtrEquals(testCase, block, stPinchSegment_getBlock(segment1));
    //get block orientation
    CuAssertIntEquals(testCase, 1, stPinchSegment_getBlockOrientation(segment1));
    //test iterator
    blockIt = stPinchBlock_getSegmentIterator(block);
    CuAssertPtrEquals(testCase, segment1, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stPinchBlockIt_getNext(&blockIt));
    stPinchBlock_destruct(block);
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment1));
    CuAssertIntEquals(testCase, 0, stPinchSegment_getBlockOrientation(segment1));

    //Test stPinchThreadSet_getTotalBlockNumber
    CuAssertIntEquals(testCase, 0, stPinchThreadSet_getTotalBlockNumber(threadSet));

    teardown();

}

static void testStPinchBlock_Splits(CuTest *testCase) {
    /*
     * Tests splitting of segments that are aligned. Put in seperate function as draws together previous code.
     */
    setup();
    static int64_t name3 = 5, start3 = 0, length3 = 15;
    stPinchThread *thread3 = stPinchThreadSet_addThread(threadSet, name3, start3, length3);
    stPinchThread_split(thread3, 4);
    stPinchThread_split(thread3, 9);
    stPinchSegment *segment1 = stPinchThread_getFirst(thread3);
    stPinchSegment *segment2 = stPinchSegment_get3Prime(segment1);
    stPinchSegment *segment3 = stPinchSegment_get3Prime(segment2);
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_get3Prime(segment3));
    stPinchBlock_pinch2(stPinchBlock_construct(segment1, 1, segment2, 1), segment3, 1);

    stPinchThread_split(thread3, 0);

    //Now traverse through thread and check all is okay
    segment1 = stPinchThread_getFirst(thread3);
    CuAssertTrue(testCase, segment1 != NULL);
    CuAssertIntEquals(testCase, 1, stPinchSegment_getLength(segment1));
    segment2 = stPinchSegment_get3Prime(segment1);
    CuAssertTrue(testCase, segment2 != NULL);
    CuAssertIntEquals(testCase, 4, stPinchSegment_getLength(segment2));
    segment3 = stPinchSegment_get3Prime(segment2);
    CuAssertTrue(testCase, segment3 != NULL);
    CuAssertIntEquals(testCase, 1, stPinchSegment_getLength(segment3));
    stPinchSegment *segment4 = stPinchSegment_get3Prime(segment3);
    CuAssertTrue(testCase, segment4 != NULL);
    CuAssertIntEquals(testCase, 4, stPinchSegment_getLength(segment4));
    stPinchSegment *segment5 = stPinchSegment_get3Prime(segment4);
    CuAssertTrue(testCase, segment5 != NULL);
    CuAssertIntEquals(testCase, 1, stPinchSegment_getLength(segment5));
    stPinchSegment *segment6 = stPinchSegment_get3Prime(segment5);
    CuAssertTrue(testCase, segment6 != NULL);
    CuAssertIntEquals(testCase, 4, stPinchSegment_getLength(segment6));
    CuAssertPtrEquals(testCase, NULL, stPinchSegment_get3Prime(segment6));
    //Check the thread.
    CuAssertPtrEquals(testCase, thread3, stPinchSegment_getThread(segment1));
    CuAssertPtrEquals(testCase, thread3, stPinchSegment_getThread(segment2));
    CuAssertPtrEquals(testCase, thread3, stPinchSegment_getThread(segment3));
    CuAssertPtrEquals(testCase, thread3, stPinchSegment_getThread(segment4));
    CuAssertPtrEquals(testCase, thread3, stPinchSegment_getThread(segment5));
    CuAssertPtrEquals(testCase, thread3, stPinchSegment_getThread(segment6));
    //Check the name
    CuAssertIntEquals(testCase, name3, stPinchSegment_getName(segment1));
    CuAssertIntEquals(testCase, name3, stPinchSegment_getName(segment2));
    CuAssertIntEquals(testCase, name3, stPinchSegment_getName(segment3));
    CuAssertIntEquals(testCase, name3, stPinchSegment_getName(segment4));
    CuAssertIntEquals(testCase, name3, stPinchSegment_getName(segment5));
    CuAssertIntEquals(testCase, name3, stPinchSegment_getName(segment6));

    //Now check blocks
    stPinchBlock *block1 = stPinchSegment_getBlock(segment1);
    CuAssertTrue(testCase, block1 != NULL);
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment1), stPinchBlock_getLength(block1));
    CuAssertIntEquals(testCase, 3, stPinchBlock_getDegree(block1));
    stPinchBlock *block2 = stPinchSegment_getBlock(segment2);
    CuAssertTrue(testCase, block2 != NULL);
    CuAssertIntEquals(testCase, stPinchSegment_getLength(segment2), stPinchBlock_getLength(block2));
    CuAssertIntEquals(testCase, 3, stPinchBlock_getDegree(block2));

    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block1);
    CuAssertPtrEquals(testCase, segment1, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment3, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment5, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, block1, stPinchSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, block1, stPinchSegment_getBlock(segment3));
    CuAssertPtrEquals(testCase, block1, stPinchSegment_getBlock(segment5));

    blockIt = stPinchBlock_getSegmentIterator(block2);
    CuAssertPtrEquals(testCase, segment2, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment4, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment6, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, block2, stPinchSegment_getBlock(segment2));
    CuAssertPtrEquals(testCase, block2, stPinchSegment_getBlock(segment4));
    CuAssertPtrEquals(testCase, block2, stPinchSegment_getBlock(segment6));

    //Test block iterator
    stPinchThreadSetBlockIt threadBlockIt = stPinchThreadSet_getBlockIt(threadSet);
    CuAssertTrue(testCase, stPinchThreadSetBlockIt_getNext(&threadBlockIt) != NULL);
    CuAssertTrue(testCase, stPinchThreadSetBlockIt_getNext(&threadBlockIt) != NULL);
    CuAssertPtrEquals(testCase, NULL, stPinchThreadSetBlockIt_getNext(&threadBlockIt));

    //Test stPinchThreadSet_getTotalBlockNumber
    CuAssertIntEquals(testCase, 2, stPinchThreadSet_getTotalBlockNumber(threadSet));

    //Now do merge
    stPinchThreadSet_joinTrivialBoundaries(threadSet);
    threadBlockIt = stPinchThreadSet_getBlockIt(threadSet);
    block1 = stPinchThreadSetBlockIt_getNext(&threadBlockIt);
    CuAssertTrue(testCase, block1 != NULL);
    CuAssertPtrEquals(testCase, NULL, stPinchThreadSetBlockIt_getNext(&threadBlockIt));
    blockIt = stPinchBlock_getSegmentIterator(block1);
    CuAssertPtrEquals(testCase, segment1, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment3, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, segment5, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, NULL, stPinchBlockIt_getNext(&blockIt));
    CuAssertPtrEquals(testCase, block1, stPinchSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, block1, stPinchSegment_getBlock(segment3));
    CuAssertPtrEquals(testCase, block1, stPinchSegment_getBlock(segment5));
    CuAssertIntEquals(testCase, 5, stPinchBlock_getLength(block1));

    //Test stPinchThreadSet_getTotalBlockNumber
    CuAssertIntEquals(testCase, 1, stPinchThreadSet_getTotalBlockNumber(threadSet));

    teardown();
}

static bool areAligned(stPinchThread *thread1, int32_t base1, stPinchThread *thread2, int32_t base2, bool orientation) {
    stPinchSegment *segment1 = stPinchThread_getSegment(thread1, base1);
    stPinchSegment *segment2 = stPinchThread_getSegment(thread2, base2);
    assert(segment1 != NULL && segment2 != NULL);
    stPinchBlock *block1 = stPinchSegment_getBlock(segment1);
    stPinchBlock *block2 = stPinchSegment_getBlock(segment2);
    if (block1 == NULL) {
        return 0;
    }
    if (block1 != block2) {
        return 0;
    }
    int32_t offset1 = base1 - stPinchSegment_getStart(segment1);
    int32_t offset2 = base2 - stPinchSegment_getStart(segment2);
    assert(base1 >= 0 && base2 >= 0);
    if (orientation) {
        return offset1 == offset2;
    }
    return stPinchBlock_getLength(block1) - 1 - offset2 == offset1;
}

static void testStPinchThread_pinchP(CuTest *testCase, int32_t segmentNumber, int64_t start, int64_t lengths[], int64_t blockDegrees[],
        stPinchThread *thread) {
    stPinchSegment *segment = stPinchThread_getFirst(thread);
    int64_t i = start;
    for (int32_t j = 0; j < segmentNumber; j++) {
        CuAssertTrue(testCase, segment != NULL);
        CuAssertIntEquals(testCase, stPinchThread_getName(thread), stPinchSegment_getName(segment));
        CuAssertIntEquals(testCase, i, stPinchSegment_getStart(segment));
        CuAssertIntEquals(testCase, lengths[j], stPinchSegment_getLength(segment));
        //Block
        if (blockDegrees[j] == 1) {
            CuAssertPtrEquals(testCase, NULL, stPinchSegment_getBlock(segment));
        } else {
            stPinchBlock *block = stPinchSegment_getBlock(segment);
            CuAssertTrue(testCase, block != NULL);
            CuAssertIntEquals(testCase, blockDegrees[j], stPinchBlock_getDegree(block));
        }

        i += stPinchSegment_getLength(segment);
        if (j + 1 == segmentNumber) {
            CuAssertTrue(testCase, stPinchSegment_get3Prime(segment) == NULL);
        } else {
            stPinchSegment *segment2 = stPinchSegment_get3Prime(segment);
            CuAssertTrue(testCase, segment2 != NULL);
            CuAssertPtrEquals(testCase, segment, stPinchSegment_get5Prime(segment2));
            segment = segment2;
        }
    }
}

static void testStPinchThread_pinch(CuTest *testCase) {
    setup();
    stPinchThread_pinch(thread1, thread2, 5, 5, 8, 1);
    int64_t lengths1[] = { 4, 6, 1, 1, length1 - 12 };
    int64_t blockDegrees1[] = { 1, 2, 2, 2, 1 };
    testStPinchThread_pinchP(testCase, 5, start1, lengths1, blockDegrees1, thread1);
    st_logInfo("First thread, first pinch okay\n");
    int64_t lengths2[] = { 1, 6, 1, 1, 1 };
    int64_t blockDegrees2[] = { 1, 2, 2, 2, 1 };
    testStPinchThread_pinchP(testCase, 5, start2, lengths2, blockDegrees2, thread2);
    st_logInfo("Second thread, first pinch okay\n");

    stPinchThread_pinch(thread1, thread2, 4, 10, 4, 0);
    int64_t lengths1b[] = { 3, 1, 1, 1, 1, 2, 1, 1, 1, length1 - 12 };
    int64_t blockDegrees1b[] = { 1, 2, 4, 4, 4, 2, 4, 4, 4, 1 };
    testStPinchThread_pinchP(testCase, 10, start1, lengths1b, blockDegrees1b, thread1);
    st_logInfo("First thread, second pinch okay\n");
    int64_t lengths2b[] = { 1, 1, 1, 1, 2, 1, 1, 1, 1 };
    int64_t blockDegrees2b[] = { 1, 4, 4, 4, 2, 4, 4, 4, 2 };
    testStPinchThread_pinchP(testCase, 9, start2, lengths2b, blockDegrees2b, thread2);
    st_logInfo("Second thread, second pinch okay\n");

    stPinchThread_pinch(thread1, thread2, 4, 10, 4, 0); //Doing the same thing again should not affect the result
    testStPinchThread_pinchP(testCase, 10, start1, lengths1b, blockDegrees1b, thread1);
    testStPinchThread_pinchP(testCase, 9, start2, lengths2b, blockDegrees2b, thread2);
    st_logInfo("Third pinch okay\n");
    //nor should a zero length pinch
    stPinchThread_pinch(thread1, thread2, 5, 10, 0, 0);
    testStPinchThread_pinchP(testCase, 10, start1, lengths1b, blockDegrees1b, thread1);
    testStPinchThread_pinchP(testCase, 9, start2, lengths2b, blockDegrees2b, thread2);
    st_logInfo("Fourth pinch okay\n");

    //Check a subset of the homology groups
    CuAssertTrue(testCase, areAligned(thread1, 4, thread2, 13, 0));

    CuAssertTrue(testCase, areAligned(thread1, 5, thread2, 5, 1));
    CuAssertTrue(testCase, areAligned(thread1, 5, thread2, 12, 0));
    CuAssertTrue(testCase, areAligned(thread1, 5, thread1, 12, 0));

    CuAssertTrue(testCase, areAligned(thread1, 6, thread2, 6, 1));
    CuAssertTrue(testCase, areAligned(thread1, 6, thread2, 11, 0));
    CuAssertTrue(testCase, areAligned(thread1, 6, thread1, 11, 0));

    CuAssertTrue(testCase, areAligned(thread1, 7, thread2, 7, 1));
    CuAssertTrue(testCase, areAligned(thread1, 7, thread2, 10, 0));
    CuAssertTrue(testCase, areAligned(thread1, 7, thread1, 10, 0));

    CuAssertTrue(testCase, areAligned(thread1, 8, thread2, 8, 1));

    CuAssertTrue(testCase, areAligned(thread1, 9, thread2, 9, 1));

    teardown();
}

//Functions that implement a very simple merging of alignment positions

static void addColumn(stHash *columns, int64_t name, int64_t position, int64_t strand) {
    assert(position > 0);
    stInt64Tuple *p = stInt64Tuple_construct(2, name, (strand ? position : -position));
    stSortedSet *column = stSortedSet_construct3((int(*)(const void *, const void *)) stInt64Tuple_cmpFn, NULL);
    stSortedSet_insert(column, p);
    stHash_insert(columns, p, column);
}

static stSortedSet *getColumn(stHash *columns, int64_t name, int64_t position, int64_t strand) {
    assert(position != 0);
    stInt64Tuple *p = stInt64Tuple_construct(2, name, strand ? position : -position);
    stSortedSet *column = stHash_search(columns, p);
    assert(column != NULL);
    stInt64Tuple_destruct(p);
    return column;
}

static void mergePositions(stHash *columns, int64_t name1, int64_t start1, bool strand1, int64_t name2, int64_t start2, bool strand2) {
    stSortedSet *column1 = getColumn(columns, name1, start1, strand1);
    stSortedSet *column2 = getColumn(columns, name2, start2, strand2);
    stSortedSet *column2R = getColumn(columns, name2, start2, !strand2);
    if (column1 != column2 && column1 != column2R) {
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
}

static void mergePositionsSymmetric(stHash *columns, int64_t name1, int64_t start1, bool strand1, int64_t name2, int64_t start2,
        bool strand2) {
    mergePositions(columns, name1, start1, strand1, name2, start2, strand2);
    mergePositions(columns, name1, start1, !strand1, name2, start2, !strand2);
}

stHash *getUnalignedColumns(stPinchThreadSet *threadSet) {
    stHash *columns = stHash_construct3((uint32_t(*)(const void *)) stInt64Tuple_hashKey,
            (int(*)(const void *, const void *)) stInt64Tuple_equalsFn, (void(*)(void *)) stInt64Tuple_destruct, NULL);
    stPinchThreadSetIt it = stPinchThreadSet_getIt(threadSet);
    stPinchThread *thread;
    while ((thread = stPinchThreadSetIt_getNext(&it))) {
        for (int32_t i = 0; i < stPinchThread_getLength(thread); i++) {
            addColumn(columns, stPinchThread_getName(thread), stPinchThread_getStart(thread) + i, 1);
            addColumn(columns, stPinchThread_getName(thread), stPinchThread_getStart(thread) + i, 0);
        }
    }
    return columns;
}

static void decodePosition(stPinchThreadSet *threadSet, stInt64Tuple *alignedPosition, stPinchThread **thread, int64_t *position,
        bool *strand) {
    *thread = stPinchThreadSet_getThread(threadSet, stInt64Tuple_getPosition(alignedPosition, 0));
    assert(*thread != NULL);
    *position = stInt64Tuple_getPosition(alignedPosition, 1);
    assert(*position != 0);
    *strand = 1;
    if (*position < 0) {
        *strand = 0;
        *position *= -1;
    }
}

static void testStPinchThread_pinch_randomTests(CuTest *testCase) {
    for (int32_t test = 0; test < 100; test++) {
        st_logInfo("Starting random pinch test %i\n", test);
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomEmptyGraph();
        stHash *columns = getUnalignedColumns(threadSet);

        //Randomly push them together, updating both sets, and checking that set of alignments is what we expect
        double threshold = st_random();
        while (st_random() > threshold) {
            stPinch pinch = stPinchThreadSet_getRandomPinch(threadSet);
            stPinchThread_pinch(stPinchThreadSet_getThread(threadSet, pinch.name1), stPinchThreadSet_getThread(threadSet, pinch.name2),
                    pinch.start1, pinch.start2, pinch.length, pinch.strand);
            //now do all the pushing together of the equivalence classes
            for (int32_t i = 0; i < pinch.length; i++) {
                mergePositionsSymmetric(columns, pinch.name1, pinch.start1 + i, 1, pinch.name2,
                        pinch.strand ? pinch.start2 + i : pinch.start2 + pinch.length - 1 - i, pinch.strand);
            }
        }
        if (st_random() > 0.5) {
            stPinchThreadSet_joinTrivialBoundaries(threadSet); //Checks this function does not affect result
        }

        //Check they are equivalent
        stHashIterator *hashIt = stHash_getIterator(columns);
        stInt64Tuple *key;
        while ((key = stHash_getNext(hashIt)) != NULL) {
            stPinchThread *thread1, *thread2;
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
                    CuAssertTrue(testCase, areAligned(thread1, position1, thread2, position2, strand1 == strand2));
                }
            }
            stList_destruct(columnList);
        }
        stHash_destructIterator(hashIt);

        stPinchThreadSet_destruct(threadSet);
        stList *columnList = stHash_getValues(columns);
        stSortedSet *columnSet = stList_getSortedSet(columnList, NULL);
        stSortedSet_setDestructor(columnSet, (void(*)(void *)) stSortedSet_destruct);
        stList_destruct(columnList);
        stHash_destruct(columns);
        stSortedSet_destruct(columnSet);
    }
}

static void testStPinchThreadSet_joinTrivialBoundaries_randomTests(CuTest *testCase) {
    for (int32_t test = 0; test < 100; test++) {
        st_logInfo("Starting random trivial boundaries test %i\n", test);
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomGraph();
        stPinchThreadSet_joinTrivialBoundaries(threadSet);
        stPinchThreadSetSegmentIt segmentIt = stPinchThreadSet_getSegmentIt(threadSet);
        stPinchSegment *segment;
        while ((segment = stPinchThreadSetSegmentIt_getNext(&segmentIt))) {
            if (stPinchSegment_getBlock(segment) == NULL) {
                stPinchSegment *segment2 = stPinchSegment_get5Prime(segment);
                CuAssertTrue(testCase, segment2 == NULL || stPinchSegment_getBlock(segment2) != NULL);
                segment2 = stPinchSegment_get3Prime(segment);
                CuAssertTrue(testCase, segment2 == NULL || stPinchSegment_getBlock(segment2) != NULL);
            }
        }
        stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
        stPinchBlock *block;
        while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
            stPinchEnd end;
            end.block = block;
            end.orientation = 0;
            CuAssertTrue(testCase, !stPinchEnd_boundaryIsTrivial(end));
            end.orientation = 1;
            CuAssertTrue(testCase, !stPinchEnd_boundaryIsTrivial(end));
        }
        stPinchThreadSet_destruct(threadSet);
    }
}

static void testStPinchThreadSet_getAdjacencyComponents(CuTest *testCase) {
    //return;
    setup();
    //Quick check that it returns what we expect
    stPinchThread_pinch(thread1, thread2, 5, 5, 8, 1);
    stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents(threadSet);
    CuAssertIntEquals(testCase, 4, stList_length(adjacencyComponents));
    stList_destruct(adjacencyComponents);
    teardown();
}

static void testStPinchThreadSet_getAdjacencyComponents_randomTests(CuTest *testCase) {
    //return;
    for (int32_t test = 0; test < 100; test++) {
        st_logInfo("Starting random adjacency component test %i\n", test);
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomGraph();
        stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents(threadSet);
        //Check all ends in one adjacency component
        stHash *ends = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, NULL, NULL);
        for (int32_t i = 0; i < stList_length(adjacencyComponents); i++) {
            stList *adjacencyComponent = stList_get(adjacencyComponents, i);
            for (int32_t j = 0; j < stList_length(adjacencyComponent); j++) {
                stPinchEnd *end = stList_get(adjacencyComponent, j);
                CuAssertPtrEquals(testCase, NULL, stHash_search(ends, end));
                stHash_insert(ends, end, end);
            }
        }
        stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
        int32_t blockNumber = 0;
        while ((stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
            blockNumber++;
        }
        CuAssertIntEquals(testCase, 2 * blockNumber, stHash_size(ends));
        //Check all connected nodes in same adjacency component
        stPinchThreadSet_destruct(threadSet);
        stHash_destruct(ends);
        stList_destruct(adjacencyComponents);
    }
}

static void testStPinchThreadSet_getThreadComponents(CuTest *testCase) {
    for (int32_t test = 0; test < 100; test++) {
        st_logInfo("Starting random thread component test %i\n", test);
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomGraph();
        stSortedSet *threadComponents = stPinchThreadSet_getThreadComponents(threadSet);
        stSortedSetIterator *it = stSortedSet_getIterator(threadComponents);
        stList *threadComponent;
        stSortedSet *threads = stSortedSet_construct();
        while ((threadComponent = stSortedSet_getNext(it)) != NULL) {
            stSortedSet *threadComponentSet = stList_getSortedSet(threadComponent, NULL);
            for (int32_t i = 0; i < stList_length(threadComponent); i++) {
                stPinchThread *thread = stList_get(threadComponent, i);
                CuAssertTrue(testCase, stSortedSet_search(threads, thread) == NULL); //Check is unique;
                stSortedSet_insert(threads, thread);
                stPinchSegment *segment = stPinchThread_getFirst(thread);
                while (segment != NULL) {
                    stPinchBlock *block;
                    if ((block = stPinchSegment_getBlock(segment)) != NULL) {
                        stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(block);
                        stPinchSegment *segment2;
                        while ((segment2 = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
                            CuAssertTrue(testCase,
                                    stSortedSet_search(threadComponentSet, stPinchSegment_getThread(segment2)) != NULL);
                        }
                    }
                    segment = stPinchSegment_get3Prime(segment);
                }
            }
            stSortedSet_destruct(threadComponentSet);
        }
        st_logInfo("There were %i threads in %i components\n", stPinchThreadSet_getSize(threadSet), stSortedSet_size(threadComponents));
        CuAssertIntEquals(testCase, stPinchThreadSet_getSize(threadSet), stSortedSet_size(threads));
        stSortedSet_destructIterator(it);
        stSortedSet_destruct(threadComponents);
        stSortedSet_destruct(threads);
        stPinchThreadSet_destruct(threadSet);
    }
}

static void testStPinchThreadSet_trimAlignments_randomTests(CuTest *testCase) {
    //return;
    for (int32_t test = 0; test < 100; test++) {
        st_logInfo("Starting random block trim test %i\n", test);
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomGraph();
        int32_t trim = st_randomInt(0, 10);
        stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
        stPinchBlock *block = stPinchThreadSetBlockIt_getNext(&blockIt);
        while (block != NULL) {
            stPinchBlock *block2 = stPinchThreadSetBlockIt_getNext(&blockIt);
            stPinchBlock_trim(block, trim);
            block = block2;
        }
        stPinchThreadSet_destruct(threadSet);
    }
}

static void testStPinchInterval(CuTest *testCase) {
    setup();
    stPinchInterval *interval = stPinchInterval_construct(name1, start1, length1, testCase);
    CuAssertIntEquals(testCase, name1, stPinchInterval_getName(interval));
    CuAssertIntEquals(testCase, start1, stPinchInterval_getStart(interval));
    CuAssertIntEquals(testCase, length1, stPinchInterval_getLength(interval));
    CuAssertPtrEquals(testCase, testCase, stPinchInterval_getLabel(interval));
    stPinchInterval interval2 = stPinchInterval_constructStatic(name2, start2, length2, testCase);
    CuAssertIntEquals(testCase, name2, stPinchInterval_getName(&interval2));
    CuAssertIntEquals(testCase, start2, stPinchInterval_getStart(&interval2));
    CuAssertIntEquals(testCase, length2, stPinchInterval_getLength(&interval2));
    CuAssertPtrEquals(testCase, testCase, stPinchInterval_getLabel(&interval2));
    CuAssertIntEquals(testCase, 0, stPinchInterval_compareFunction(interval, interval));
    CuAssertIntEquals(testCase, 0, stPinchInterval_compareFunction(&interval2, &interval2));
    CuAssertIntEquals(testCase, -1, stPinchInterval_compareFunction(interval, &interval2));
    CuAssertIntEquals(testCase, 1, stPinchInterval_compareFunction(&interval2, interval));
    stSortedSet *intervals = stSortedSet_construct3((int(*)(const void *, const void *)) stPinchInterval_compareFunction, NULL);
    CuAssertPtrEquals(testCase, NULL, stSortedSet_search(intervals, interval));
    CuAssertPtrEquals(testCase, NULL, stSortedSet_search(intervals, &interval2));
    stSortedSet_insert(intervals, interval);
    stSortedSet_insert(intervals, &interval2);
    CuAssertPtrEquals(testCase, interval, stSortedSet_search(intervals, interval));
    CuAssertPtrEquals(testCase, &interval2, stSortedSet_search(intervals, &interval2));
    CuAssertPtrEquals(testCase, interval, stPinchIntervals_getInterval(intervals, name1, start1));
    CuAssertPtrEquals(testCase, interval, stPinchIntervals_getInterval(intervals, name1, start1 + length1 - 1));
    CuAssertPtrEquals(testCase, NULL, stPinchIntervals_getInterval(intervals, name1, start1 + length1));
    CuAssertPtrEquals(testCase, NULL, stPinchIntervals_getInterval(intervals, name1, start1 - 1));
    stSortedSet_destruct(intervals);
    stPinchInterval_destruct(interval);
    teardown();
}

static void testStPinchThreadSet_getLabelIntervals(CuTest *testCase) {
    setup();
    //Tests when there are no blocks in the problem
    stHash *pinchEndsToAdjacencyComponents;
    stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents2(threadSet, &pinchEndsToAdjacencyComponents);
    stSortedSet *intervals = stPinchThreadSet_getLabelIntervals(threadSet, pinchEndsToAdjacencyComponents);
    CuAssertIntEquals(testCase, 2, stSortedSet_size(intervals));
    stPinchInterval *interval = stPinchIntervals_getInterval(intervals, name1, start1);
    CuAssertTrue(testCase, interval != NULL);
    CuAssertIntEquals(testCase, name1, stPinchInterval_getName(interval));
    CuAssertIntEquals(testCase, start1, stPinchInterval_getStart(interval));
    CuAssertIntEquals(testCase, length1, stPinchInterval_getLength(interval));
    CuAssertPtrEquals(testCase, NULL, stPinchInterval_getLabel(interval));
    interval = stPinchIntervals_getInterval(intervals, name2, start2);
    CuAssertTrue(testCase, interval != NULL);
    CuAssertIntEquals(testCase, name2, stPinchInterval_getName(interval));
    CuAssertIntEquals(testCase, start2, stPinchInterval_getStart(interval));
    CuAssertIntEquals(testCase, length2, stPinchInterval_getLength(interval));
    CuAssertPtrEquals(testCase, NULL, stPinchInterval_getLabel(interval));
    stHash_destruct(pinchEndsToAdjacencyComponents);
    stSortedSet_destruct(intervals);
    stList_destruct(adjacencyComponents);
    teardown();
}

static stList *getAdjacencyComponentP(stPinchSegment *segment, int64_t position, stHash *pinchEndsToAdjacencyComponents) {
    stPinchBlock *block = stPinchSegment_getBlock(segment);
    assert(block != NULL);
    bool orientation = (position >= stPinchSegment_getLength(segment) / 2) ^ stPinchSegment_getBlockOrientation(segment);
    stPinchEnd pinchEnd = stPinchEnd_constructStatic(block, orientation);
    stList *adjacencyComponent = stHash_search(pinchEndsToAdjacencyComponents, &pinchEnd);
    assert(adjacencyComponent != NULL);
    return adjacencyComponent;
}

static stList *getAdjacencyComponent(stPinchSegment *segment, int64_t position, stHash *pinchEndsToAdjacencyComponents) {
    if (stPinchSegment_getBlock(segment) != NULL) {
        return getAdjacencyComponentP(segment, position, pinchEndsToAdjacencyComponents);
    }
    stPinchSegment *segment2 = stPinchSegment_get3Prime(segment);
    while (segment2 != NULL) {
        if (stPinchSegment_getBlock(segment2) != NULL) {
            return getAdjacencyComponentP(segment2, -1, pinchEndsToAdjacencyComponents);
        }
        segment2 = stPinchSegment_get3Prime(segment2);
    }
    segment2 = stPinchSegment_get5Prime(segment);
    while (segment2 != NULL) {
        if (stPinchSegment_getBlock(segment2) != NULL) {
            return getAdjacencyComponentP(segment2, INT32_MAX, pinchEndsToAdjacencyComponents);
        }
        segment2 = stPinchSegment_get5Prime(segment2);
    }
    return NULL;
}

static void testStPinchThreadSet_getLabelIntervals_randomTests(CuTest *testCase) {
    //return;
    for (int32_t test = 0; test < 100; test++) {
        st_logInfo("Starting random get label intervals test %i\n", test);
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomGraph();
        stHash *pinchEndsToAdjacencyComponents;
        stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents2(threadSet, &pinchEndsToAdjacencyComponents);
        stSortedSet *intervals = stPinchThreadSet_getLabelIntervals(threadSet, pinchEndsToAdjacencyComponents);
        //Check every base is in a label interval
        stPinchThreadSetSegmentIt segmentIt = stPinchThreadSet_getSegmentIt(threadSet);
        stPinchSegment *segment;
        while ((segment = stPinchThreadSetSegmentIt_getNext(&segmentIt)) != NULL) {
            for (int32_t i = 0; i < stPinchSegment_getLength(segment); i++) {
                stPinchInterval *interval = stPinchIntervals_getInterval(intervals, stPinchSegment_getName(segment),
                        stPinchSegment_getStart(segment) + i);
                CuAssertTrue(testCase, interval != NULL);
                CuAssertIntEquals(testCase, stPinchSegment_getName(segment), stPinchInterval_getName(interval));
                CuAssertTrue(testCase, stPinchSegment_getStart(segment) + i >= stPinchInterval_getStart(interval));
                CuAssertTrue(
                        testCase,
                        stPinchSegment_getStart(segment) + i < stPinchInterval_getStart(interval)
                        + stPinchInterval_getLength(interval));
                //Now check out stupid way of calculating the adjacency component is the same as the calculated by the label function.
                CuAssertPtrEquals(testCase, getAdjacencyComponent(segment, i, pinchEndsToAdjacencyComponents),
                        stPinchInterval_getLabel(interval));
            }
        }
        //Cleanup
        stSortedSet_destruct(intervals);
        stHash_destruct(pinchEndsToAdjacencyComponents);
        stList_destruct(adjacencyComponents);
        stPinchThreadSet_destruct(threadSet);
    }
}

CuSuite* stPinchGraphsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testStPinchThreadSet);
    SUITE_ADD_TEST(suite, testStPinchThreadAndSegment);
    SUITE_ADD_TEST(suite, testStPinchBlock_NoSplits);
    SUITE_ADD_TEST(suite, testStPinchBlock_Splits);
    SUITE_ADD_TEST(suite, testStPinchThread_pinch);
    SUITE_ADD_TEST(suite, testStPinchThread_pinch_randomTests);
    SUITE_ADD_TEST(suite, testStPinchThreadSet_getAdjacencyComponents);
    SUITE_ADD_TEST(suite, testStPinchThreadSet_getAdjacencyComponents_randomTests);
    SUITE_ADD_TEST(suite, testStPinchThreadSet_joinTrivialBoundaries_randomTests);
    SUITE_ADD_TEST(suite, testStPinchThreadSet_getThreadComponents);
    SUITE_ADD_TEST(suite, testStPinchThreadSet_trimAlignments_randomTests);
    SUITE_ADD_TEST(suite, testStPinchInterval);
    SUITE_ADD_TEST(suite, testStPinchThreadSet_getLabelIntervals);
    SUITE_ADD_TEST(suite, testStPinchThreadSet_getLabelIntervals_randomTests);

    return suite;
}
