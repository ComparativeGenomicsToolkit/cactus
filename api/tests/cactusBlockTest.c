/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusBlocksTestShared.h"

static bool nestedTest = 0;

static void cactusBlockTestTeardown(const char *testName) {
    if (!nestedTest) {
        cactusBlocksTestSharedTeardown(testName);
    }
}

static void cactusBlockTestSetup(const char *testName) {
    if (!nestedTest) {
        cactusBlocksTestSharedSetup(testName);
    }
}

void testBlock_construct(CuTest* testCase) {
    nestedTest = 0;
    cactusBlockTestSetup(testCase->name);
    CuAssertTrue(testCase, block != NULL);
    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getName(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    CuAssertTrue(testCase, block_getName(block) != NULL_NAME);
    CuAssertTrue(testCase, flower_getBlock(flower, block_getName(block)) == block);
    CuAssertTrue(testCase, block_getName(block) == block_getName(block_getReverse(block)));
    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getOrientation(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    CuAssertTrue(testCase, block_getOrientation(block));
    CuAssertTrue(testCase, !block_getOrientation(block_getReverse(block)));
    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getReverse(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    CuAssertTrue(testCase, block_getReverse(block) != NULL);
    CuAssertTrue(testCase, block_getReverse(block_getReverse(block)) == block);
    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getLength(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    CuAssertIntEquals(testCase, 3, block_getLength(block));
    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getFlower(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    CuAssertTrue(testCase, block_getFlower(block) == flower);
    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getLeftEnd(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    End *leftEnd = block_get5End(block);
    CuAssertTrue(testCase, leftEnd != NULL);
    CuAssertTrue(testCase, end_getBlock(leftEnd) == block);
    CuAssertTrue(testCase, block_getOrientation(block) == end_getOrientation(leftEnd));
    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getRightEnd(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    End *rightEnd = block_get3End(block);
    CuAssertTrue(testCase, rightEnd != NULL);
    CuAssertTrue(testCase, end_getBlock(rightEnd) == block);
    CuAssertTrue(testCase, block_getOrientation(block) == end_getOrientation(rightEnd));
    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getInstanceNumber(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    CuAssertTrue(testCase, block_getInstanceNumber(block) == 3);
    Block *block2 = block_construct(1, flower);
    CuAssertTrue(testCase, block_getInstanceNumber(block2) == 0);
    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getInstance(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    CuAssertTrue(testCase, block_getInstance(block, segment_getName(rootSegment)) == segment_getReverse(rootSegment));
    CuAssertTrue(testCase, block_getInstance(block, segment_getName(leaf1Segment)) == leaf1Segment);
    CuAssertTrue(testCase, block_getInstance(block, segment_getName(leaf2Segment)) == segment_getReverse(leaf2Segment));

    CuAssertTrue(testCase, block_getInstance(block_getReverse(block), segment_getName(rootSegment)) == rootSegment);
    CuAssertTrue(testCase, block_getInstance(block_getReverse(block), segment_getName(leaf1Segment)) == segment_getReverse(leaf1Segment));
    CuAssertTrue(testCase, block_getInstance(block_getReverse(block), segment_getName(leaf2Segment)) == leaf2Segment);

    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getFirst(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    CuAssertTrue(testCase, block_getFirst(block) == segment_getReverse(rootSegment));
    CuAssertTrue(testCase, block_getFirst(block_getReverse(block)) == rootSegment);
    Block *block2 = block_construct(1, flower);
    CuAssertTrue(testCase, block_getFirst(block2) == NULL);
    CuAssertTrue(testCase, block_getFirst(block_getReverse(block2)) == NULL);
    cactusBlockTestTeardown(testCase->name);
}

void testBlock_instanceIterator(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);
    Block_InstanceIterator *iterator;
    iterator = block_getInstanceIterator(block);

    CuAssertTrue(testCase, iterator != NULL);
    CuAssertTrue(testCase, block_getNext(iterator) == segment_getReverse(rootSegment));
    CuAssertTrue(testCase, block_getNext(iterator) == leaf1Segment);
    CuAssertTrue(testCase, block_getNext(iterator) == segment_getReverse(leaf2Segment));
    CuAssertTrue(testCase, block_getNext(iterator) == NULL);

    block_destructInstanceIterator(iterator);

    iterator = block_getInstanceIterator(block_getReverse(block));
    CuAssertTrue(testCase, block_getNext(iterator) == rootSegment);
    CuAssertTrue(testCase, block_getNext(iterator) == segment_getReverse(leaf1Segment));
    CuAssertTrue(testCase, block_getNext(iterator) == leaf2Segment);
    CuAssertTrue(testCase, block_getNext(iterator) == NULL);

    block_destructInstanceIterator(iterator);

    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getChain(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);

    Block *block2 = block_construct(2, flower);
    Flower *flower2 = flower_construct(cactusDisk);
    end_copyConstruct(block_get3End(block), flower2);
    end_copyConstruct(block_get5End(block2), flower2);
    Group *group = group_construct(flower, flower2);
    Chain *chain = chain_construct(flower);
    link_construct(block_get3End(block), block_get5End(block2), group, chain);

    CuAssertTrue(testCase, block_getChain(block_construct(2, flower)) == NULL);
    CuAssertTrue(testCase, block_getChain(block) == chain);
    CuAssertTrue(testCase, block_getChain(block2) == chain);

    cactusBlockTestTeardown(testCase->name);
}

void testBlock_getSegmentForEvent(CuTest* testCase) {
    cactusBlockTestSetup(testCase->name);

    CuAssertPtrEquals(testCase, block_getSegmentForEvent(block_getReverse(block), event_getName(rootEvent)), rootSegment);
    Segment *segment = block_getSegmentForEvent(block, event_getName(leafEvent));
    CuAssertTrue(testCase, segment == leaf1Segment || segment == segment_getReverse(leaf2Segment));
    CuAssertTrue(testCase, block_getSegmentForEvent(block, NULL_NAME) == NULL);

    cactusBlockTestTeardown(testCase->name);
}

void testBlock_isTrivialChain(CuTest *testCase) {
    cactusBlockTestSetup(testCase->name);
    Group *group = group_construct2(flower);
    end_setGroup(block_get5End(block), group);
    end_setGroup(block_get3End(block), group);
    Chain *chain = chain_construct(flower);
    Group *group2 = group_construct2(flower);
    CuAssertTrue(testCase, block_isTrivialChain(block));
    Block *block1 = block_construct(1, flower);
    Block *block2 = block_construct(1, flower);
    end_setGroup(block_get5End(block2), group2);
    end_setGroup(block_get3End(block1), group2);
    link_construct(block_get3End(block1), block_get5End(block2), group2, chain);
    end_setGroup(block_get5End(block1), group);
    end_setGroup(block_get3End(block2), group);
    CuAssertTrue(testCase, block_isTrivialChain(block));
    CuAssertTrue(testCase, !block_isTrivialChain(block1));
    CuAssertTrue(testCase, !block_isTrivialChain(block2));
    cactusBlockTestTeardown(testCase->name);
}

CuSuite* cactusBlockTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testBlock_getName);
    SUITE_ADD_TEST(suite, testBlock_getOrientation);
    SUITE_ADD_TEST(suite, testBlock_getReverse);
    SUITE_ADD_TEST(suite, testBlock_getLength);
    SUITE_ADD_TEST(suite, testBlock_getFlower);
    SUITE_ADD_TEST(suite, testBlock_getLeftEnd);
    SUITE_ADD_TEST(suite, testBlock_getRightEnd);
    SUITE_ADD_TEST(suite, testBlock_getInstanceNumber);
    SUITE_ADD_TEST(suite, testBlock_getInstance);
    SUITE_ADD_TEST(suite, testBlock_getFirst);
    SUITE_ADD_TEST(suite, testBlock_instanceIterator);
    SUITE_ADD_TEST(suite, testBlock_getChain);
    SUITE_ADD_TEST(suite, testBlock_getSegmentForEvent);
    SUITE_ADD_TEST(suite, testBlock_isTrivialChain);
    SUITE_ADD_TEST(suite, testBlock_construct);
    return suite;
}
