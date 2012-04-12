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
    CuAssertPtrEquals(testCase, thread1, stThreadIt_getNext(threadIt));
    CuAssertPtrEquals(testCase, thread2, stThreadIt_getNext(threadIt));
    CuAssertPtrEquals(testCase, NULL, stThreadIt_getNext(threadIt));
    CuAssertPtrEquals(testCase, NULL, stThreadIt_getNext(threadIt));
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
    CuAssertPtrEquals(testCase, NULL, stThread_get5Prime(thread1, segment));
    CuAssertIntEquals(testCase, start1, stSegment_getStart(segment));
    CuAssertIntEquals(testCase, leftSplitPoint1 - start1, stSegment_getLength(segment));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment));

    stSegment *segment2 = stThread_get3Prime(thread1, segment);
    CuAssertTrue(testCase, segment2 != NULL);
    CuAssertPtrEquals(testCase, segment, stThread_get5Prime(thread1, segment2));
    CuAssertIntEquals(testCase, leftSplitPoint1+1, stSegment_getStart(segment2));
    CuAssertIntEquals(testCase, leftSplitPoint2 - leftSplitPoint1+1, stSegment_getLength(segment2));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment));

    stSegment *segment3 = stThread_get3Prime(thread1, segment2);
    CuAssertTrue(testCase, segment3 != NULL);
    CuAssertPtrEquals(testCase, segment2, stThread_get5Prime(thread1, segment3));
    CuAssertPtrEquals(testCase, NULL, stThread_get3Prime(thread1, segment3));
    CuAssertIntEquals(testCase, leftSplitPoint2+1, stSegment_getStart(segment3));
    CuAssertIntEquals(testCase, length1 - leftSplitPoint2+1 + start1, stSegment_getLength(segment3));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment));

    //Test stThread_getSegment
    CuAssertPtrEquals(testCase, segment, stThread_getSegment(thread1, start1));
    CuAssertPtrEquals(testCase, segment, stThread_getSegment(thread1, leftSplitPoint1-1));
    CuAssertPtrEquals(testCase, segment2, stThread_getSegment(thread1, leftSplitPoint1));
    CuAssertPtrEquals(testCase, segment2, stThread_getSegment(thread1, leftSplitPoint2-1));
    CuAssertPtrEquals(testCase, segment3, stThread_getSegment(thread1, leftSplitPoint2));
    CuAssertPtrEquals(testCase, segment3, stThread_getSegment(thread1, start1 + length1 - 1));
    CuAssertPtrEquals(testCase, NULL, stThread_getSegment(thread1, start1 + length1 - 1));

    //Test stThread_joinTrivialBoundaries
    stThread_joinTrivialBoundaries(thread1);
    //Now should be just one segment
    segment = stThread_getFirst(thread1);
    CuAssertTrue(testCase, segment != NULL);
    CuAssertPtrEquals(testCase, NULL, stThread_get5Prime(thread1, segment));
    CuAssertPtrEquals(testCase, NULL, stThread_get3Prime(thread1, segment));
    CuAssertIntEquals(testCase, start1, stSegment_getStart(segment));
    CuAssertIntEquals(testCase, length1, stSegment_getLength(segment));
    CuAssertPtrEquals(testCase, NULL, stSegment_getBlock(segment));
    CuAssertIntEquals(testCase, 0, stSegment_getBlockOrientation(segment));
    CuAssertPtrEquals(testCase, segment, stThread_getSegment(thread1, start1));
    CuAssertPtrEquals(testCase, segment, stThread_getSegment(thread1, start1+length1-1));
    CuAssertPtrEquals(testCase, NULL, stThread_getSegment(thread1, start1-1));
    CuAssertPtrEquals(testCase, NULL, stThread_getSegment(thread1, start1+length1));

    teardown();
}

static void testStBlock(CuTest *testCase) {
    setup();
    static int64_t name3 = 5, start3 = 0, length3 = 15;
    stThread *thread3 = stThreadSet_addThread(threadSet, name3, start3, length3);
    stThread_split(thread3, 4);
    stThread_split(thread3, 9);
    stSegment *segment1 = stThread_getFirst(stThread3);
    stSegment *segment2 = stThread_get3Prime(segment1);
    stSegment *segment3 = stThread_get3Prime(segment2);
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stSegment_getLength(segment2));
    CuAssertIntEquals(testCase, stSegment_getLength(segment1), stSegment_getLength(segment3));
    CuAssertPtrEquals(testCase, NULL, stThread_get3Prime(segment3));
    stBlock *block = stBlock_construct(segment1, 1, segment2, 0);
    //get block
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment1));
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment2));
    //get block orientation
    CuAssertIntEquals(testCase, 1, stSegment_getBlock(segment1));
    CuAssertIntEquals(testCase, 0, stSegment_getBlock(segment2));
    //test iterator
    stBlockIt blockIt = stBlock_getSegmentIterator(block);
    CuAssertPtrEquals(testCase, segment1, stBlock_getNext(blockIt));
    CuAssertPtrEquals(testCase, segment2, stBlock_getNext(blockIt));
    CuAssertPtrEquals(testCase, NULL, stBlock_getNext(blockIt));
    CuAssertPtrEquals(testCase, NULL, stBlock_getNext(blockIt));

    //Now try pinching in a segment
    stBlock_pinch2(block, segment3, 1);
    CuAssertPtrEquals(testCase, block, stSegment_getBlock(segment2));


    //Now try removing from the block

    teardown();
}

CuSuite* stPinchGraphsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testStThreadSet);
    SUITE_ADD_TEST(suite, testStThreadAndSegment);
    SUITE_ADD_TEST(suite, testStBlock);
    return suite;
}
