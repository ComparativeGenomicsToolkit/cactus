/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusReferenceTestShared.h"

static void testTeardown() {
    if (!nestedTest) {
        cactusReferenceTestSharedTeardown();
    }
}

static void testSetup() {
    if (!nestedTest) {
        cactusReferenceTestSharedSetup();
    }
}

static void testPseudoAdjacency_construct(CuTest* testCase) {
    testSetup();
    assert(testCase != NULL);
    //already tested by shared code.
    testTeardown();
}

static void testPseudoAdjacency_getName(CuTest* testCase) {
    testSetup();
    CuAssertTrue(testCase, pseudoAdjacency_getName(pseudoAdjacency1) != NULL_NAME);
    CuAssertTrue(testCase, pseudoAdjacency_getName(pseudoAdjacency2) != NULL_NAME);
    CuAssertTrue(testCase, pseudoAdjacency_getName(pseudoAdjacency1) != pseudoAdjacency_getName(pseudoAdjacency2));
    testTeardown();
}

static void testPseudoAdjacency_get5End(CuTest* testCase) {
    testSetup();
    CuAssertTrue(testCase, pseudoAdjacency_get5End(pseudoAdjacency1) == end1);
    CuAssertTrue(testCase, pseudoAdjacency_get5End(pseudoAdjacency2) == end3);
    CuAssertTrue(testCase, pseudoAdjacency_get5End(pseudoAdjacency3) == end5);
    CuAssertTrue(testCase, pseudoAdjacency_get5End(pseudoAdjacency4) == end7);

    CuAssertTrue(testCase, end_getPseudoAdjacency(end1) == pseudoAdjacency1);
    CuAssertTrue(testCase, end_getPseudoAdjacency(end3) == pseudoAdjacency2);
    CuAssertTrue(testCase, end_getPseudoAdjacency(end5) == pseudoAdjacency3);
    CuAssertTrue(testCase, end_getPseudoAdjacency(end7) == pseudoAdjacency4);
    testTeardown();
}

static void testPseudoAdjacency_get3End(CuTest* testCase) {
    testSetup();
    CuAssertTrue(testCase, pseudoAdjacency_get3End(pseudoAdjacency1) == end2);
    CuAssertTrue(testCase, pseudoAdjacency_get3End(pseudoAdjacency2) == end4);
    CuAssertTrue(testCase, pseudoAdjacency_get3End(pseudoAdjacency3) == end6);
    CuAssertTrue(testCase, pseudoAdjacency_get3End(pseudoAdjacency4) == end8);

    CuAssertTrue(testCase, end_getPseudoAdjacency(end2) == pseudoAdjacency1);
    CuAssertTrue(testCase, end_getPseudoAdjacency(end4) == pseudoAdjacency2);
    CuAssertTrue(testCase, end_getPseudoAdjacency(end6) == pseudoAdjacency3);
    CuAssertTrue(testCase, end_getPseudoAdjacency(end8) == pseudoAdjacency4);
    testTeardown();
}

static void testPseudoAdjacency_getPseudoChromosome(CuTest* testCase) {
    testSetup();
    CuAssertTrue(testCase, pseudoAdjacency_getPseudoChromosome(pseudoAdjacency1) == pseudoChromosome1);
    CuAssertTrue(testCase, pseudoAdjacency_getPseudoChromosome(pseudoAdjacency2) == pseudoChromosome1);
    CuAssertTrue(testCase, pseudoAdjacency_getPseudoChromosome(pseudoAdjacency3) == pseudoChromosome1);
    CuAssertTrue(testCase, pseudoAdjacency_getPseudoChromosome(pseudoAdjacency4) == pseudoChromosome2);
    testTeardown();
}

static void testPseudoAdjacency_getIndex(CuTest* testCase) {
    return;
    testSetup();
    CuAssertTrue(testCase, pseudoAdjacency_getIndex(pseudoAdjacency1) == 0);
    CuAssertTrue(testCase, pseudoAdjacency_getIndex(pseudoAdjacency2) == 1);
    CuAssertTrue(testCase, pseudoAdjacency_getIndex(pseudoAdjacency3) == 2);
    //CuAssertTrue(testCase, pseudoAdjacency_getIndex(pseudoAdjacency4) == 3);
    testTeardown();
}

CuSuite* cactusPseudoAdjacencyTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, testPseudoAdjacency_getName);
    SUITE_ADD_TEST(suite, testPseudoAdjacency_get5End);
    SUITE_ADD_TEST(suite, testPseudoAdjacency_get3End);
    SUITE_ADD_TEST(suite, testPseudoAdjacency_getPseudoChromosome);
    SUITE_ADD_TEST(suite, testPseudoAdjacency_getIndex);
    SUITE_ADD_TEST(suite, testPseudoAdjacency_construct);

    return suite;
}
