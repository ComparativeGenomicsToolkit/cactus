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

void testLink_getPreviousLink(CuTest* testCase) {
    cactusLinkTestSetup(testCase);
    CuAssertTrue(testCase, link_getPreviousLink(link2) == link1);
    CuAssertTrue(testCase, link_getPreviousLink(link1) == NULL);
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

void testLink_split(CuTest *testCase) {
    cactusLinkTestSetup(testCase);
    CuAssertTrue(testCase, flower_getChainNumber(flower) == 3);
    CuAssertTrue(testCase, chain_getLength(chain) == 2);
    CuAssertTrue(testCase, group_getLink(group1) == link1);
    link_split(link1);
    CuAssertTrue(testCase, group_getLink(group1) == NULL);
    CuAssertTrue(testCase, flower_getChainNumber(flower) == 3);
    chain = flower_getFirstChain(flower);
    Flower_ChainIterator *chainIt = flower_getChainIterator(flower);
    while((chain = flower_getNextChain(chainIt)) != NULL && (chain == chain2 || chain == chain3));
    flower_destructChainIterator(chainIt);
    assert(chain != NULL);
    CuAssertTrue(testCase, chain_getLength(chain) == 1);
    Link *link3 = chain_getFirst(chain);
    CuAssertTrue(testCase, link_get3End(link3) == block_get3End(block));
    CuAssertTrue(testCase, link_get5End(link3) == end2);
    CuAssertTrue(testCase, group_getLink(group2) == link3);
    link_split(chain_getFirst(chain));
    CuAssertTrue(testCase, flower_getChainNumber(flower) == 2);
    CuAssertTrue(testCase, group_getLink(group1) == NULL);
    CuAssertTrue(testCase, group_getLink(group2) == NULL);
    cactusLinkTestTeardown(testCase);
}

void testLink_serialisation(CuTest* testCase) {
    cactusLinkTestSetup(testCase);
    int64_t i;
    void
            *vA =
                    binaryRepresentation_makeBinaryRepresentation(link2,
                            (void(*)(void *, void(*)(const void *, size_t,
                                    size_t))) link_writeBinaryRepresentation,
                            &i);
    CuAssertTrue(testCase, i > 0);
    link_destruct(link2);
    void *vA2 = vA;
    link2 = link_loadFromBinaryRepresentation(&vA2, chain);
    nestedTest = 1;
    testLink_getNextLink(testCase);
    testLink_getPreviousLink(testCase);
    testLink_getGroup(testCase);
    testLink_getLeft(testCase);
    testLink_getRight(testCase);
    testLink_getChain(testCase);
    //testLink_split(testCase); //can't do that test, because it disrupts stuff..
    nestedTest = 0;
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

void testLink_mergeIfTrivial(CuTest *testCase) {
    cactusLinkTestSetup(testCase);
    CuAssertTrue(testCase, !link_mergeIfTrivial(link1));
    CuAssertTrue(testCase, !link_mergeIfTrivial(link2));
    //CuAssertTrue(testCase, !link_mergeIfTrivial(link5));

    CuAssertTrue(testCase, chain_getLength(chain2) == 1);
    CuAssertTrue(testCase, flower_getBlockNumber(flower) == 4);
    CuAssertTrue(testCase, flower_getBlockEndNumber(flower) == 8);
    CuAssertTrue(testCase, flower_getStubEndNumber(flower) == 2);
    CuAssertTrue(testCase, flower_getGroupNumber(flower) == 5);
    CuAssertTrue(testCase, flower_getSegmentNumber(flower) == 3);
    CuAssertTrue(testCase, flower_getCapNumber(flower) == 6);

    CuAssertTrue(testCase, link_mergeIfTrivial(link4));

    CuAssertTrue(testCase, chain_getLength(chain2) == 0);
    CuAssertTrue(testCase, flower_getBlockNumber(flower) == 3);
    CuAssertTrue(testCase, flower_getBlockEndNumber(flower) == 6);
    CuAssertTrue(testCase, flower_getStubEndNumber(flower) == 2);
    CuAssertTrue(testCase, flower_getGroupNumber(flower) == 4);
    CuAssertTrue(testCase, flower_getSegmentNumber(flower) == 2);
    CuAssertTrue(testCase, flower_getCapNumber(flower) == 4);
    Segment *mergedSegment;
    Flower_SegmentIterator *it = flower_getSegmentIterator(flower);
    while((mergedSegment = flower_getNextSegment(it)) != NULL && (mergedSegment == segment3));
    flower_destructSegmentIterator(it);
    CuAssertTrue(testCase, mergedSegment != NULL);
    CuAssertTrue(testCase, mergedSegment != segment3);
    Block *mergedBlock = segment_getBlock(mergedSegment);
    CuAssertTrue(testCase, segment_getLength(mergedSegment) == 2);
    CuAssertTrue(testCase, block_getLength(mergedBlock) == 2);
    CuAssertTrue(testCase, block_isTrivialChain(mergedBlock));

    cactusLinkTestTeardown(testCase);
}

CuSuite* cactusLinkTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testLink_getNextLink);
    SUITE_ADD_TEST(suite, testLink_getPreviousLink);
    SUITE_ADD_TEST(suite, testLink_getGroup);
    SUITE_ADD_TEST(suite, testLink_getLeft);
    SUITE_ADD_TEST(suite, testLink_getRight);
    SUITE_ADD_TEST(suite, testLink_getChain);
    SUITE_ADD_TEST(suite, testLink_serialisation);
    SUITE_ADD_TEST(suite, testLink_split);
    SUITE_ADD_TEST(suite, testLink_isTrivial);
    SUITE_ADD_TEST(suite, testLink_mergeIfTrivial);
    SUITE_ADD_TEST(suite, testLink_construct);
    return suite;
}
