/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusBlocksTestShared.h"

static bool nestedTest = 0;

static void cactusSegmentTestTeardown(CuTest* testCase) {
	if(!nestedTest) {
		cactusBlocksTestSharedTeardown(testCase->name);
	}
}

static void cactusSegmentTestSetup(CuTest* testCase) {
	if(!nestedTest) {
		cactusBlocksTestSharedSetup(testCase->name);
	}
}

void testSegment_construct(CuTest* testCase) {
	nestedTest = 0;
	cactusSegmentTestSetup(testCase);
	CuAssertTrue(testCase, rootSegment != NULL);
	CuAssertTrue(testCase, leaf1Segment != NULL);
	cactusSegmentTestTeardown(testCase);
}

void testSegment_getBlock(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);
	CuAssertTrue(testCase, segment_getBlock(rootSegment) == block_getReverse(block));
	CuAssertTrue(testCase, segment_getBlock(leaf1Segment) == block);
	CuAssertTrue(testCase, segment_getBlock(leaf2Segment) == block_getReverse(block));

	CuAssertTrue(testCase, segment_getBlock(segment_getReverse(rootSegment)) == block);
	CuAssertTrue(testCase, segment_getBlock(segment_getReverse(leaf1Segment)) == block_getReverse(block));
	CuAssertTrue(testCase, segment_getBlock(segment_getReverse(leaf2Segment)) == block);

	cactusSegmentTestTeardown(testCase);
}

void testSegment_getName(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);

	CuAssertTrue(testCase, segment_getName(rootSegment) != NULL_NAME);
	CuAssertTrue(testCase, block_getInstance(block, segment_getName(rootSegment)) == segment_getReverse(rootSegment));
	CuAssertTrue(testCase, block_getInstance(block_getReverse(block), segment_getName(rootSegment)) == rootSegment);

	CuAssertTrue(testCase, segment_getName(leaf1Segment) != NULL_NAME);
	CuAssertTrue(testCase, block_getInstance(block, segment_getName(leaf1Segment)) == leaf1Segment);
	CuAssertTrue(testCase, block_getInstance(block_getReverse(block), segment_getName(leaf1Segment)) == segment_getReverse(leaf1Segment));

	CuAssertTrue(testCase, segment_getName(leaf2Segment) != NULL_NAME);
	CuAssertTrue(testCase, block_getInstance(block, segment_getName(leaf2Segment)) == segment_getReverse(leaf2Segment));
	CuAssertTrue(testCase, block_getInstance(block_getReverse(block), segment_getName(leaf2Segment)) == leaf2Segment);

	cactusSegmentTestTeardown(testCase);
}

void testSegment_getOrientation(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);

	CuAssertTrue(testCase, segment_getOrientation(rootSegment) == block_getOrientation(block_getReverse(block)));
	CuAssertTrue(testCase, segment_getOrientation(leaf1Segment) == block_getOrientation(block));
	CuAssertTrue(testCase, segment_getOrientation(leaf2Segment) == block_getOrientation(block_getReverse(block)));

	CuAssertTrue(testCase, segment_getOrientation(leaf1Segment) != segment_getOrientation(leaf2Segment));

	cactusSegmentTestTeardown(testCase);
}

void testSegment_getReverse(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);

	CuAssertTrue(testCase, segment_getReverse(rootSegment) != NULL);
	CuAssertTrue(testCase, segment_getReverse(segment_getReverse(rootSegment)) == rootSegment);

	CuAssertTrue(testCase, segment_getReverse(leaf1Segment) != NULL);
	CuAssertTrue(testCase, segment_getReverse(segment_getReverse(leaf1Segment)) == leaf1Segment);

	CuAssertTrue(testCase, segment_getReverse(leaf2Segment) != NULL);
	CuAssertTrue(testCase, segment_getReverse(segment_getReverse(leaf2Segment)) == leaf2Segment);

	cactusSegmentTestTeardown(testCase);
}

void testSegment_getEvent(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);

	CuAssertTrue(testCase, segment_getEvent(rootSegment) == rootEvent);
	CuAssertTrue(testCase, segment_getEvent(leaf1Segment) == leafEvent);
	CuAssertTrue(testCase, segment_getEvent(leaf2Segment) == leafEvent);

	cactusSegmentTestTeardown(testCase);
}

void testSegment_getStart(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);

	CuAssertTrue(testCase, segment_getStart(rootSegment) == INT64_MAX);
	CuAssertTrue(testCase, segment_getStart(segment_getReverse(rootSegment)) == INT64_MAX);

	CuAssertTrue(testCase, segment_getStart(leaf1Segment) == 2);
	CuAssertTrue(testCase, segment_getStart(segment_getReverse(leaf1Segment)) == 4);

	CuAssertTrue(testCase, segment_getStart(leaf2Segment) == 6);
	CuAssertTrue(testCase, segment_getStart(segment_getReverse(leaf2Segment)) == 4);

	cactusSegmentTestTeardown(testCase);
}

void testSegment_getStrand(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);

	CuAssertTrue(testCase, segment_getStrand(rootSegment));
	CuAssertTrue(testCase, segment_getStrand(rootSegment) != segment_getStrand(segment_getReverse(rootSegment)));
	CuAssertTrue(testCase, segment_getStrand(rootSegment) == cap_getStrand(segment_get5Cap(rootSegment)));
	CuAssertTrue(testCase, segment_getStrand(rootSegment) == cap_getStrand(segment_get3Cap(rootSegment)));
	CuAssertTrue(testCase, segment_getStrand(segment_getReverse(rootSegment)) == cap_getStrand(cap_getReverse(segment_get5Cap(rootSegment))));
	CuAssertTrue(testCase, segment_getStrand(segment_getReverse(rootSegment)) == cap_getStrand(cap_getReverse(segment_get3Cap(rootSegment))));

	CuAssertTrue(testCase, segment_getStrand(leaf1Segment));
	CuAssertTrue(testCase, segment_getStrand(leaf1Segment) != segment_getStrand(segment_getReverse(leaf1Segment)));
	CuAssertTrue(testCase, segment_getStrand(leaf1Segment) == cap_getStrand(segment_get5Cap(leaf1Segment)));
	CuAssertTrue(testCase, segment_getStrand(leaf1Segment) == cap_getStrand(segment_get3Cap(leaf1Segment)));
	CuAssertTrue(testCase, segment_getStrand(segment_getReverse(leaf1Segment)) == cap_getStrand(cap_getReverse(segment_get5Cap(leaf1Segment))));
	CuAssertTrue(testCase, segment_getStrand(segment_getReverse(leaf1Segment)) == cap_getStrand(cap_getReverse(segment_get3Cap(leaf1Segment))));

	CuAssertTrue(testCase, !segment_getStrand(leaf2Segment));
	CuAssertTrue(testCase, segment_getStrand(leaf2Segment) != segment_getStrand(segment_getReverse(leaf2Segment)));
	CuAssertTrue(testCase, segment_getStrand(leaf2Segment) == cap_getStrand(segment_get5Cap(leaf2Segment)));
	CuAssertTrue(testCase, segment_getStrand(leaf2Segment) == cap_getStrand(segment_get3Cap(leaf2Segment)));
	CuAssertTrue(testCase, segment_getStrand(segment_getReverse(leaf2Segment)) == cap_getStrand(cap_getReverse(segment_get5Cap(leaf2Segment))));
	CuAssertTrue(testCase, segment_getStrand(segment_getReverse(leaf2Segment)) == cap_getStrand(cap_getReverse(segment_get3Cap(leaf2Segment))));

	cactusSegmentTestTeardown(testCase);
}

void testSegment_getLength(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);

	CuAssertTrue(testCase, segment_getLength(rootSegment) == block_getLength(block));
	CuAssertTrue(testCase, segment_getLength(segment_getReverse(rootSegment)) == block_getLength(block));

	CuAssertTrue(testCase, segment_getLength(leaf1Segment) == block_getLength(block));
	CuAssertTrue(testCase, segment_getLength(segment_getReverse(leaf1Segment)) == block_getLength(block));

	CuAssertTrue(testCase, segment_getLength(leaf2Segment) == block_getLength(block));
	CuAssertTrue(testCase, segment_getLength(segment_getReverse(leaf2Segment)) == block_getLength(block));

	cactusSegmentTestTeardown(testCase);
}

void testSegment_getSequence(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);

	CuAssertTrue(testCase, segment_getSequence(rootSegment) == NULL);
	CuAssertTrue(testCase, segment_getSequence(segment_getReverse(rootSegment)) == NULL);

	CuAssertTrue(testCase, segment_getSequence(leaf1Segment) == sequence);
	CuAssertTrue(testCase, segment_getSequence(segment_getReverse(leaf1Segment)) == sequence);

	CuAssertTrue(testCase, segment_getSequence(leaf2Segment) == sequence);
	CuAssertTrue(testCase, segment_getSequence(segment_getReverse(leaf2Segment)) == sequence);

	cactusSegmentTestTeardown(testCase);
}

void testSegment_getString(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);

	CuAssertTrue(testCase, segment_getString(rootSegment) == NULL);
	CuAssertTrue(testCase, segment_getString(segment_getReverse(rootSegment)) == NULL);

	CuAssertStrEquals(testCase, "CTG", segment_getString(leaf1Segment));
	CuAssertStrEquals(testCase, "CAG", segment_getString(segment_getReverse(leaf1Segment)));

	CuAssertStrEquals(testCase, "GTC", segment_getString(leaf2Segment));
	CuAssertStrEquals(testCase, "GAC", segment_getString(segment_getReverse(leaf2Segment)));

	cactusSegmentTestTeardown(testCase);
}

void testSegment_get5And3End(CuTest* testCase) {
	cactusSegmentTestSetup(testCase);

	Cap *cap5 = segment_get5Cap(rootSegment);
	Cap *cap3 = segment_get3Cap(rootSegment);

	CuAssertTrue(testCase, cap_getSegment(cap5) == rootSegment);
	CuAssertTrue(testCase, cap_getSegment(cap3) == rootSegment);
	CuAssertTrue(testCase, cap_getOrientation(cap5) == segment_getOrientation(rootSegment));
	CuAssertTrue(testCase, cap_getOrientation(cap3) == segment_getOrientation(rootSegment));

	cap5 = segment_get5Cap(segment_getReverse(rootSegment));
	cap3 = segment_get3Cap(segment_getReverse(rootSegment));

	CuAssertTrue(testCase, cap_getSegment(cap5) == segment_getReverse(rootSegment));
	CuAssertTrue(testCase, cap_getSegment(cap3) == segment_getReverse(rootSegment));
	CuAssertTrue(testCase, cap_getOrientation(cap5) == segment_getOrientation(segment_getReverse(rootSegment)));
	CuAssertTrue(testCase, cap_getOrientation(cap3) == segment_getOrientation(segment_getReverse(rootSegment)));

	cactusSegmentTestTeardown(testCase);
}

CuSuite* cactusSegmentTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testSegment_getBlock);
	SUITE_ADD_TEST(suite, testSegment_getName);
	SUITE_ADD_TEST(suite, testSegment_getOrientation);
	SUITE_ADD_TEST(suite, testSegment_getReverse);
	SUITE_ADD_TEST(suite, testSegment_getEvent);
	SUITE_ADD_TEST(suite, testSegment_getStart);
	SUITE_ADD_TEST(suite, testSegment_getStrand);
	SUITE_ADD_TEST(suite, testSegment_getLength);
	SUITE_ADD_TEST(suite, testSegment_getSequence);
	SUITE_ADD_TEST(suite, testSegment_getString);
	SUITE_ADD_TEST(suite, testSegment_get5And3End);
	SUITE_ADD_TEST(suite, testSegment_construct);
	return suite;
}
