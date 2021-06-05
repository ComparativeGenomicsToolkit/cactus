/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusChainsTestShared.h"

static bool nestedTest = 0;

void cactusLinkTestSetup(CuTest* testCase) {
    if (!nestedTest) {
        cactusChainsSharedTestSetup(testCase->name);
    }
}

void cactusLinkTestTeardown(CuTest* testCase) {
    if (!nestedTest) {
        cactusChainsSharedTestTeardown(testCase->name);
    }
}

void testLink_construct(CuTest* testCase) {
    nestedTest = 0;
    cactusLinkTestSetup(testCase);
    CuAssertTrue(testCase, link1 != NULL);
    CuAssertTrue(testCase, link2 != NULL);
    cactusLinkTestTeardown(testCase);
}

void testLink_getNextLink(CuTest* testCase) {
    cactusLinkTestSetup(testCase);
    CuAssertTrue(testCase, link_getNextLink(link1) == link2);
    CuAssertTrue(testCase, link_getNextLink(link2) == NULL);
    cactusLinkTestTeardown(testCase);
}

void testLink_getGroup(CuTest* testCase) {
    cactusLinkTestSetup(testCase);
    CuAssertTrue(testCase, link_getGroup(link1) == group1);
    CuAssertTrue(testCase, link_getGroup(link2) == group2);
    cactusLinkTestTeardown(testCase);
}

void testLink_getLeft(CuTest* testCase) {
    cactusLinkTestSetup(testCase);
    CuAssertTrue(testCase, link_get3End(link1) == end1);
    CuAssertTrue(testCase, link_get3End(link2) == block_get3End(block));
    cactusLinkTestTeardown(testCase);
}

void testLink_getRight(CuTest* testCase) {
    cactusLinkTestSetup(testCase);
    CuAssertTrue(testCase, link_get5End(link1) == block_get5End(block));
    CuAssertTrue(testCase, link_get5End(link2) == end2);
    cactusLinkTestTeardown(testCase);
}

void testLink_getChain(CuTest* testCase) {
    cactusLinkTestSetup(testCase);
    CuAssertTrue(testCase, link_getChain(link1) == chain);
    CuAssertTrue(testCase, link_getChain(link2) == chain);
    cactusLinkTestTeardown(testCase);
}

void testLink_isTrivial(CuTest *testCase) {
    cactusLinkTestSetup(testCase);
    CuAssertTrue(testCase, !link_isTrivial(link1));
    CuAssertTrue(testCase, !link_isTrivial(link2));
    CuAssertTrue(testCase, link_isTrivial(link4));
    //CuAssertTrue(testCase, !link_isTrivial(link5));
    cactusLinkTestTeardown(testCase);
}

CuSuite* cactusLinkTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testLink_getNextLink);
    SUITE_ADD_TEST(suite, testLink_getGroup);
    SUITE_ADD_TEST(suite, testLink_getLeft);
    SUITE_ADD_TEST(suite, testLink_getRight);
    SUITE_ADD_TEST(suite, testLink_getChain);
    SUITE_ADD_TEST(suite, testLink_isTrivial);
    SUITE_ADD_TEST(suite, testLink_construct);
    return suite;
}
