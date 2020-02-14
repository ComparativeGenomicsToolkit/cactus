/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "flowersShared.h"
#include "adjacencySequences.h"

static void testAdjacencySequence_1(CuTest *testCase) {
   setup(testCase);
   AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap1, INT64_MAX);
   CuAssertTrue(testCase, adjacencySequence->subsequenceIdentifier == cap_getName(cap1)); //sequence_getName(sequence1));
   CuAssertIntEquals(testCase, adjacencySequence->start, 1);
   CuAssertIntEquals(testCase, adjacencySequence->strand, 1);
   CuAssertIntEquals(testCase, adjacencySequence->length, 4);
   CuAssertStrEquals(testCase, "ACTG", adjacencySequence->string);
   adjacencySequence_destruct(adjacencySequence);
   teardown(testCase);
}

static void testAdjacencySequence_2(CuTest *testCase) {
   setup(testCase);
   AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap_getReverse(cap2), INT64_MAX);
   CuAssertTrue(testCase, adjacencySequence->subsequenceIdentifier == cap_getName(cap1)); //sequence_getName(sequence1));
   CuAssertIntEquals(testCase, adjacencySequence->start, 4);
   CuAssertIntEquals(testCase, adjacencySequence->strand, 0);
   CuAssertIntEquals(testCase, adjacencySequence->length, 4);
   CuAssertStrEquals(testCase, "CAGT", adjacencySequence->string);
   adjacencySequence_destruct(adjacencySequence);
   teardown(testCase);
}

static void testAdjacencySequence_3(CuTest *testCase) {
   setup(testCase);
   AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap1, 2);
   CuAssertTrue(testCase, adjacencySequence->subsequenceIdentifier ==  cap_getName(cap1)); //sequence_getName(sequence1));
   CuAssertIntEquals(testCase, adjacencySequence->start, 1);
   CuAssertIntEquals(testCase, adjacencySequence->strand, 1);
   CuAssertIntEquals(testCase, adjacencySequence->length, 2);
   CuAssertStrEquals(testCase, "AC", adjacencySequence->string);
   adjacencySequence_destruct(adjacencySequence);
   teardown(testCase);
}

static void testAdjacencySequence_4(CuTest *testCase) {
   setup(testCase);
   AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap_getReverse(cap2), 0);
   CuAssertTrue(testCase, adjacencySequence->subsequenceIdentifier ==  cap_getName(cap1)); //sequence_getName(sequence1));
   CuAssertIntEquals(testCase, adjacencySequence->start, 4);
   CuAssertIntEquals(testCase, adjacencySequence->strand, 0);
   CuAssertIntEquals(testCase, adjacencySequence->length, 0);
   CuAssertStrEquals(testCase, "", adjacencySequence->string);
   adjacencySequence_destruct(adjacencySequence);
   teardown(testCase);
}

static void testAdjacencySequence_5(CuTest *testCase) {
    setup(testCase);
   AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap7, INT64_MAX);
   CuAssertTrue(testCase, adjacencySequence->subsequenceIdentifier == cap_getName(cap8)); //sequence_getName(sequence2));
   CuAssertIntEquals(testCase, adjacencySequence->start, 6);
   CuAssertIntEquals(testCase, adjacencySequence->strand, 0);
   CuAssertIntEquals(testCase, adjacencySequence->length, 6);
   CuAssertStrEquals(testCase, "CCGGTT", adjacencySequence->string);
   adjacencySequence_destruct(adjacencySequence);
   teardown(testCase);
}

static void testAdjacencySequence_6(CuTest *testCase) {
    setup(testCase);
   AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap9, INT64_MAX);
   CuAssertTrue(testCase, adjacencySequence->subsequenceIdentifier == cap_getName(cap9)); //sequence_getName(sequence3));
   CuAssertIntEquals(testCase, adjacencySequence->start, 1);
   CuAssertIntEquals(testCase, adjacencySequence->strand, 1);
   CuAssertIntEquals(testCase, adjacencySequence->length, 4);
   CuAssertStrEquals(testCase, "CGGG", adjacencySequence->string);
   adjacencySequence_destruct(adjacencySequence);
   teardown(testCase);
}

static void testAdjacencySequence_7(CuTest *testCase) {
   setup(testCase);
   AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap11, INT64_MAX);
   CuAssertTrue(testCase, adjacencySequence->subsequenceIdentifier == cap_getName(cap11)); //sequence_getName(sequence4));
   CuAssertIntEquals(testCase, adjacencySequence->start, 2);
   CuAssertIntEquals(testCase, adjacencySequence->strand, 1);
   CuAssertIntEquals(testCase, adjacencySequence->length, 0);
   CuAssertStrEquals(testCase, "", adjacencySequence->string);
   adjacencySequence_destruct(adjacencySequence);
   teardown(testCase);
}

CuSuite* adjacencySequenceTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testAdjacencySequence_1);
    SUITE_ADD_TEST(suite, testAdjacencySequence_2);
    SUITE_ADD_TEST(suite, testAdjacencySequence_3);
    SUITE_ADD_TEST(suite, testAdjacencySequence_4);
    SUITE_ADD_TEST(suite, testAdjacencySequence_5);
    SUITE_ADD_TEST(suite, testAdjacencySequence_6);
    SUITE_ADD_TEST(suite, testAdjacencySequence_7);
    return suite;
}
