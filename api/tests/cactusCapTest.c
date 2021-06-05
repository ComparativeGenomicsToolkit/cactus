/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusEndsTestShared.h"

static bool nestedTest = 0;

void cactusCapTestSetup(CuTest* testCase) {
    if (!nestedTest) {
        cactusEndsTestSharedSetup(testCase->name);
    }
}

void cactusCapTestTeardown(CuTest* testCase) {
    if (!nestedTest) {
        cactusEndsTestSharedTeardown(testCase->name);
    }
}

void testCap_construct(CuTest* testCase) {
    nestedTest = 0;
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, rootCap != NULL);
    cactusCapTestTeardown(testCase);
}

void testCap_getName(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getName(rootCap) != NULL_NAME);
    CuAssertTrue(testCase, end_getInstance(end, cap_getName(rootCap)) == cap_getReverse(rootCap));
    CuAssertTrue(testCase, cap_getName(leaf2Cap) != NULL_NAME);
    CuAssertTrue(testCase, end_getInstance(end, cap_getName(leaf2Cap)) == leaf2Cap);

    cactusCapTestTeardown(testCase);
}

void testCap_getOrientation(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getOrientation(rootCap) == end_getOrientation(cap_getEnd(rootCap)));
    CuAssertTrue(testCase, cap_getOrientation(leaf1Cap) == end_getOrientation(cap_getEnd(leaf1Cap)));
    CuAssertTrue(testCase, cap_getOrientation(leaf2Cap) == end_getOrientation(cap_getEnd(leaf2Cap)));

    CuAssertTrue(testCase, cap_getOrientation(cap_getReverse(rootCap)) == end_getOrientation(end_getReverse(cap_getEnd(rootCap))));
    CuAssertTrue(testCase, cap_getOrientation(cap_getReverse(leaf1Cap)) == end_getOrientation(end_getReverse(cap_getEnd(leaf1Cap))));
    CuAssertTrue(testCase, cap_getOrientation(cap_getReverse(leaf2Cap)) == end_getOrientation(end_getReverse(cap_getEnd(leaf2Cap))));

    CuAssertTrue(testCase, cap_getOrientation(leaf1Cap) == cap_getOrientation(rootCap));
    CuAssertTrue(testCase, cap_getOrientation(leaf1Cap) != cap_getOrientation(leaf2Cap));

    cactusCapTestTeardown(testCase);
}

void testCap_getReverse(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getReverse(rootCap) != NULL);
    CuAssertTrue(testCase, cap_getReverse(cap_getReverse(rootCap)) == rootCap);
    cactusCapTestTeardown(testCase);
}

void testCap_getEvent(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getEvent(rootCap) == rootEvent);
    CuAssertTrue(testCase, cap_getEvent(cap_getReverse(rootCap)) == rootEvent);
    cactusCapTestTeardown(testCase);
}

void testCap_getEnd(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getEnd(rootCap) == end_getReverse(end));
    CuAssertTrue(testCase, cap_getEnd(cap_getReverse(rootCap)) == end);
    CuAssertTrue(testCase, cap_getEnd(leaf2Cap) == end);
    CuAssertTrue(testCase, cap_getEnd(cap_getReverse(leaf2Cap)) == end_getReverse(end));
    cactusCapTestTeardown(testCase);
}

void testCap_getSegment(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    Block *block = block_construct(2, flower);
    Segment *segment = segment_construct(block, rootEvent);
    CuAssertTrue(testCase, cap_getSegment(segment_get5Cap(segment)) == segment);
    CuAssertTrue(testCase, cap_getSegment(segment_get3Cap(segment)) == segment);
    CuAssertTrue(testCase, cap_getOrientation(segment_get5Cap(segment)) == segment_getOrientation(segment));
    CuAssertTrue(testCase, cap_getOrientation(segment_get3Cap(segment)) == segment_getOrientation(segment));
    CuAssertTrue(testCase, cap_getSegment(cap_getReverse(segment_get5Cap(segment))) == segment_getReverse(segment));
    CuAssertTrue(testCase, cap_getSegment(cap_getReverse(segment_get3Cap(segment))) == segment_getReverse(segment));
    cactusCapTestTeardown(testCase);
}

void testCap_getOtherSegmentCap(CuTest *testCase) {
    cactusCapTestSetup(testCase);

    Block *block = block_construct(3, flower);
    Segment *segment = segment_construct2(block, 2, 1, sequence);
    Cap *_5Cap = segment_get5Cap(segment);
    Cap *_3Cap = segment_get3Cap(segment);

    CuAssertTrue(testCase, cap_getOtherSegmentCap(leaf1Cap) == NULL);
    CuAssertTrue(testCase, cap_getOtherSegmentCap(cap_getReverse(leaf1Cap)) == NULL);

    CuAssertTrue(testCase, cap_getOtherSegmentCap(_5Cap) == _3Cap);
    CuAssertTrue(testCase, cap_getOtherSegmentCap(_3Cap) == _5Cap);

    CuAssertTrue(testCase, cap_getOtherSegmentCap(cap_getReverse(_5Cap)) == cap_getReverse(_3Cap));
    CuAssertTrue(testCase, cap_getOtherSegmentCap(cap_getReverse(_3Cap)) == cap_getReverse(_5Cap));

    cactusCapTestTeardown(testCase);
}

void testCap_segmentCoordinates(CuTest* testCase) {
    /*
     * Tests the coordinates of an segment and its 5 and 3 prime caps.
     */
    cactusCapTestSetup(testCase);

    Block *block = block_construct(3, flower);
    Segment *segment = segment_construct2(block, 2, 1, sequence);
    Cap *_5Cap = segment_get5Cap(segment);
    Cap *_3Cap = segment_get3Cap(segment);

    CuAssertTrue(testCase, cap_getSide(_5Cap));
    CuAssertTrue(testCase, !cap_getSide(_3Cap));
    CuAssertTrue(testCase, cap_getStrand(_5Cap));
    CuAssertTrue(testCase, cap_getStrand(_3Cap));
    CuAssertIntEquals(testCase, 2, cap_getCoordinate(_5Cap));
    CuAssertIntEquals(testCase, 4, cap_getCoordinate(_3Cap));
    CuAssertTrue(testCase, segment_getStrand(segment));
    CuAssertIntEquals(testCase, 2, segment_getStart(segment));
    CuAssertIntEquals(testCase, 3, segment_getLength(segment));

    CuAssertTrue(testCase, !cap_getSide(cap_getReverse(_5Cap)));
    CuAssertTrue(testCase, cap_getSide(cap_getReverse(_3Cap)));
    CuAssertTrue(testCase, !cap_getStrand(cap_getReverse(_5Cap)));
    CuAssertTrue(testCase, !cap_getStrand(cap_getReverse(_3Cap)));
    CuAssertIntEquals(testCase, 2, cap_getCoordinate(cap_getReverse(_5Cap)));
    CuAssertIntEquals(testCase, 4, cap_getCoordinate(cap_getReverse(_3Cap)));
    CuAssertTrue(testCase, !segment_getStrand(segment_getReverse(segment)));
    CuAssertIntEquals(testCase, 4, segment_getStart(segment_getReverse(segment)));
    CuAssertIntEquals(testCase, 3, segment_getLength(segment_getReverse(segment)));

    cactusCapTestTeardown(testCase);
}

void testCap_segmentCoordinatesReverseStrand(CuTest* testCase) {
    /*
     * Tests the coordinates of an segment and its 5 and 3 prime caps.
     */
    cactusCapTestSetup(testCase);

    Block *block = block_construct(3, flower);
    Segment *segment = segment_construct2(block, 2, 0, sequence);
    Cap *_5Cap = segment_get5Cap(segment);
    Cap *_3Cap = segment_get3Cap(segment);

    CuAssertTrue(testCase, cap_getSide(_5Cap));
    CuAssertTrue(testCase, !cap_getSide(_3Cap));
    CuAssertTrue(testCase, !cap_getStrand(_5Cap));
    CuAssertTrue(testCase, !cap_getStrand(_3Cap));
    CuAssertIntEquals(testCase, 4, cap_getCoordinate(_5Cap));
    CuAssertIntEquals(testCase, 2, cap_getCoordinate(_3Cap));
    CuAssertTrue(testCase, !segment_getStrand(segment));
    CuAssertIntEquals(testCase, 4, segment_getStart(segment));
    CuAssertIntEquals(testCase, 3, segment_getLength(segment));

    CuAssertTrue(testCase, !cap_getSide(cap_getReverse(_5Cap)));
    CuAssertTrue(testCase, cap_getSide(cap_getReverse(_3Cap)));
    CuAssertTrue(testCase, cap_getStrand(cap_getReverse(_5Cap)));
    CuAssertTrue(testCase, cap_getStrand(cap_getReverse(_3Cap)));
    CuAssertIntEquals(testCase, 4, cap_getCoordinate(cap_getReverse(_5Cap)));
    CuAssertIntEquals(testCase, 2, cap_getCoordinate(cap_getReverse(_3Cap)));
    CuAssertTrue(testCase, segment_getStrand(segment_getReverse(segment)));
    CuAssertIntEquals(testCase, 2, segment_getStart(segment_getReverse(segment)));
    CuAssertIntEquals(testCase, 3, segment_getLength(segment_getReverse(segment)));

    cactusCapTestTeardown(testCase);
}

void testCap_getCoordinate(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getCoordinate(rootCap) == INT64_MAX);
    CuAssertTrue(testCase, cap_getCoordinate(cap_getReverse(rootCap)) == INT64_MAX);
    CuAssertTrue(testCase, cap_getCoordinate(leaf1Cap) == 4);
    CuAssertTrue(testCase, cap_getCoordinate(cap_getReverse(leaf1Cap)) == 4);
    cactusCapTestTeardown(testCase);
}

void testCap_setCoordinate(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getCoordinate(rootCap) == INT64_MAX);
    CuAssertTrue(testCase, cap_getStrand(rootCap));
    CuAssertTrue(testCase, cap_getSequence(rootCap) == NULL);
    cap_setCoordinates(rootCap, 5, 0, NULL);
    CuAssertTrue(testCase, cap_getCoordinate(rootCap) == 5);
    CuAssertTrue(testCase, !cap_getStrand(rootCap));
    CuAssertTrue(testCase, cap_getSequence(rootCap) == NULL);
    cap_setCoordinates(rootCap, INT64_MAX, 1, NULL);
    CuAssertTrue(testCase, cap_getCoordinate(rootCap) == INT64_MAX);
    CuAssertTrue(testCase, cap_getStrand(rootCap));
    CuAssertTrue(testCase, cap_getSequence(rootCap) == NULL);
    cactusCapTestTeardown(testCase);
}

void testCap_getStrand(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getStrand(rootCap));
    CuAssertTrue(testCase, !cap_getStrand(cap_getReverse(rootCap)));
    CuAssertTrue(testCase, cap_getStrand(leaf1Cap));
    CuAssertTrue(testCase, !cap_getStrand(cap_getReverse(leaf1Cap)));
    CuAssertTrue(testCase, !cap_getStrand(leaf2Cap));
    CuAssertTrue(testCase, cap_getStrand(cap_getReverse(leaf2Cap)));
    cactusCapTestTeardown(testCase);
}

void testCap_getSide(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, !cap_getSide(rootCap));
    CuAssertTrue(testCase, cap_getSide(cap_getReverse(rootCap)));
    CuAssertTrue(testCase, !cap_getSide(leaf1Cap));
    CuAssertTrue(testCase, cap_getSide(cap_getReverse(leaf1Cap)));
    CuAssertTrue(testCase, cap_getSide(leaf2Cap));
    CuAssertTrue(testCase, !cap_getSide(cap_getReverse(leaf2Cap)));
    cactusCapTestTeardown(testCase);
}

void testCap_getSequence(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getSequence(rootCap) == NULL);
    CuAssertTrue(testCase, cap_getSequence(cap_getReverse(rootCap)) == NULL);
    CuAssertTrue(testCase, cap_getSequence(leaf1Cap) == sequence);
    CuAssertTrue(testCase, cap_getSequence(cap_getReverse(leaf1Cap)) == sequence);
    cactusCapTestTeardown(testCase);
}

void testCap_adjacent(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getAdjacency(leaf1Cap) == NULL);
    CuAssertTrue(testCase, cap_getAdjacency(leaf3Cap) == NULL);
    cap_makeAdjacent(leaf1Cap, leaf3Cap);
    CuAssertTrue(testCase, cap_getAdjacency(leaf1Cap) == cap_getReverse(leaf3Cap));
    CuAssertTrue(testCase, cap_getAdjacency(leaf3Cap) == cap_getReverse(leaf1Cap));
    CuAssertTrue(testCase, cap_getAdjacency(cap_getReverse(leaf1Cap)) == leaf3Cap);
    CuAssertTrue(testCase, cap_getAdjacency(cap_getReverse(leaf3Cap)) == leaf1Cap);
    cactusCapTestTeardown(testCase);
}

CuSuite* cactusCapTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testCap_getName);
    SUITE_ADD_TEST(suite, testCap_getOrientation);
    SUITE_ADD_TEST(suite, testCap_getReverse);
    SUITE_ADD_TEST(suite, testCap_getEvent);
    SUITE_ADD_TEST(suite, testCap_getEnd);
    SUITE_ADD_TEST(suite, testCap_getSegment);
    SUITE_ADD_TEST(suite, testCap_getOtherSegmentCap);
    SUITE_ADD_TEST(suite, testCap_segmentCoordinates);
    SUITE_ADD_TEST(suite, testCap_segmentCoordinatesReverseStrand);
    SUITE_ADD_TEST(suite, testCap_getCoordinate);
    SUITE_ADD_TEST(suite, testCap_setCoordinate);
    SUITE_ADD_TEST(suite, testCap_getStrand);
    SUITE_ADD_TEST(suite, testCap_getSide);
    SUITE_ADD_TEST(suite, testCap_getSequence);
    SUITE_ADD_TEST(suite, testCap_adjacent);
    SUITE_ADD_TEST(suite, testCap_construct);
    return suite;
}
