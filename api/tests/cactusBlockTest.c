/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusBlocksTestShared.h"

static bool nestedTest = 0;

static void cactusBlockTestTeardown() {
    if (!nestedTest) {
        cactusBlocksTestSharedTeardown();
    }
}

static void cactusBlockTestSetup() {
    if (!nestedTest) {
        cactusBlocksTestSharedSetup();
    }
}

void testBlock_construct(CuTest* testCase) {
    nestedTest = 0;
    cactusBlockTestSetup();
    CuAssertTrue(testCase, block != NULL);
    cactusBlockTestTeardown();
}

void testBlock_getName(CuTest* testCase) {
    cactusBlockTestSetup();
    CuAssertTrue(testCase, block_getName(block) != NULL_NAME);
    CuAssertTrue(testCase, flower_getBlock(flower, block_getName(block)) == block);
    CuAssertTrue(testCase, block_getName(block) == block_getName(block_getReverse(block)));
    cactusBlockTestTeardown();
}

void testBlock_getOrientation(CuTest* testCase) {
    cactusBlockTestSetup();
    CuAssertTrue(testCase, block_getOrientation(block));
    CuAssertTrue(testCase, !block_getOrientation(block_getReverse(block)));
    cactusBlockTestTeardown();
}

void testBlock_getReverse(CuTest* testCase) {
    cactusBlockTestSetup();
    CuAssertTrue(testCase, block_getReverse(block) != NULL);
    CuAssertTrue(testCase, block_getReverse(block_getReverse(block)) == block);
    cactusBlockTestTeardown();
}

void testBlock_getLength(CuTest* testCase) {
    cactusBlockTestSetup();
    CuAssertIntEquals(testCase, 3, block_getLength(block));
    cactusBlockTestTeardown();
}

void testBlock_getFlower(CuTest* testCase) {
    cactusBlockTestSetup();
    CuAssertTrue(testCase, block_getFlower(block) == flower);
    cactusBlockTestTeardown();
}

void testBlock_getLeftEnd(CuTest* testCase) {
    cactusBlockTestSetup();
    End *leftEnd = block_get5End(block);
    CuAssertTrue(testCase, leftEnd != NULL);
    CuAssertTrue(testCase, end_getBlock(leftEnd) == block);
    CuAssertTrue(testCase, block_getOrientation(block) == end_getOrientation(leftEnd));
    cactusBlockTestTeardown();
}

void testBlock_getRightEnd(CuTest* testCase) {
    cactusBlockTestSetup();
    End *rightEnd = block_get3End(block);
    CuAssertTrue(testCase, rightEnd != NULL);
    CuAssertTrue(testCase, end_getBlock(rightEnd) == block);
    CuAssertTrue(testCase, block_getOrientation(block) == end_getOrientation(rightEnd));
    cactusBlockTestTeardown();
}

void testBlock_getInstanceNumber(CuTest* testCase) {
    cactusBlockTestSetup();
    CuAssertTrue(testCase, block_getInstanceNumber(block) == 3);
    Block *block2 = block_construct(1, flower);
    CuAssertTrue(testCase, block_getInstanceNumber(block2) == 0);
    cactusBlockTestTeardown();
}

void testBlock_getInstance(CuTest* testCase) {
    cactusBlockTestSetup();
    CuAssertTrue(testCase, block_getInstance(block, segment_getName(rootSegment)) == segment_getReverse(rootSegment));
    CuAssertTrue(testCase, block_getInstance(block, segment_getName(leaf1Segment)) == leaf1Segment);
    CuAssertTrue(testCase, block_getInstance(block, segment_getName(leaf2Segment)) == segment_getReverse(leaf2Segment));

    CuAssertTrue(testCase, block_getInstance(block_getReverse(block), segment_getName(rootSegment)) == rootSegment);
    CuAssertTrue(testCase, block_getInstance(block_getReverse(block), segment_getName(leaf1Segment)) == segment_getReverse(leaf1Segment));
    CuAssertTrue(testCase, block_getInstance(block_getReverse(block), segment_getName(leaf2Segment)) == leaf2Segment);

    cactusBlockTestTeardown();
}

void testBlock_getFirst(CuTest* testCase) {
    cactusBlockTestSetup();
    CuAssertTrue(testCase, block_getFirst(block) == segment_getReverse(rootSegment));
    CuAssertTrue(testCase, block_getFirst(block_getReverse(block)) == rootSegment);
    Block *block2 = block_construct(1, flower);
    CuAssertTrue(testCase, block_getFirst(block2) == NULL);
    CuAssertTrue(testCase, block_getFirst(block_getReverse(block2)) == NULL);
    cactusBlockTestTeardown();
}

void testBlock_getSetRootInstance(CuTest *testCase) {
    cactusBlockTestSetup();
    Block *block2 = block_construct(1, flower);
    CuAssertTrue(testCase, block_getRootInstance(block2) == NULL);
    block_destruct(block2);
    //block_setRootInstance(block, segment_getReverse(rootSegment)); //set in the constructor code of the test.
    CuAssertTrue(testCase, block_getRootInstance(block) == segment_getReverse(rootSegment));
    CuAssertTrue(testCase, end_getRootInstance(block_get5End(block)) == segment_get5Cap(segment_getReverse(rootSegment)));
    CuAssertTrue(testCase, end_getRootInstance(block_get3End(block)) == segment_get3Cap(segment_getReverse(rootSegment)));
    cactusBlockTestTeardown();
}

void testBlock_instanceIterator(CuTest* testCase) {
    cactusBlockTestSetup();
    Block_InstanceIterator *iterator;
    iterator = block_getInstanceIterator(block);

    CuAssertTrue(testCase, iterator != NULL);
    CuAssertTrue(testCase, block_getNext(iterator) == segment_getReverse(rootSegment));
    CuAssertTrue(testCase, block_getNext(iterator) == leaf1Segment);

    Block_InstanceIterator *iterator2 = block_copyInstanceIterator(iterator);

    CuAssertTrue(testCase, block_getNext(iterator) == segment_getReverse(leaf2Segment));
    CuAssertTrue(testCase, block_getNext(iterator) == NULL);
    CuAssertTrue(testCase, block_getPrevious(iterator) == segment_getReverse(leaf2Segment));
    CuAssertTrue(testCase, block_getPrevious(iterator) == leaf1Segment);
    CuAssertTrue(testCase, block_getPrevious(iterator) == segment_getReverse(rootSegment));
    CuAssertTrue(testCase, block_getPrevious(iterator) == NULL);

    CuAssertTrue(testCase, block_getNext(iterator2) == segment_getReverse(leaf2Segment));
    CuAssertTrue(testCase, block_getNext(iterator2) == NULL);
    CuAssertTrue(testCase, block_getPrevious(iterator2) == segment_getReverse(leaf2Segment));
    CuAssertTrue(testCase, block_getPrevious(iterator2) == leaf1Segment);
    CuAssertTrue(testCase, block_getPrevious(iterator2) == segment_getReverse(rootSegment));
    CuAssertTrue(testCase, block_getPrevious(iterator2) == NULL);

    block_destructInstanceIterator(iterator);
    block_destructInstanceIterator(iterator2);

    iterator = block_getInstanceIterator(block_getReverse(block));
    CuAssertTrue(testCase, block_getNext(iterator) == rootSegment);
    CuAssertTrue(testCase, block_getNext(iterator) == segment_getReverse(leaf1Segment));
    CuAssertTrue(testCase, block_getNext(iterator) == leaf2Segment);
    CuAssertTrue(testCase, block_getNext(iterator) == NULL);
    CuAssertTrue(testCase, block_getPrevious(iterator) == leaf2Segment);
    CuAssertTrue(testCase, block_getPrevious(iterator) == segment_getReverse(leaf1Segment));
    CuAssertTrue(testCase, block_getPrevious(iterator) == rootSegment);
    CuAssertTrue(testCase, block_getPrevious(iterator) == NULL);

    block_destructInstanceIterator(iterator);

    cactusBlockTestTeardown();
}

void testBlock_getChain(CuTest* testCase) {
    cactusBlockTestSetup();

    Block *block2 = block_construct(2, flower);
    Flower *flower2 = flower_construct(cactusDisk);
    eventTree_copyConstruct(flower_getEventTree(flower), flower2, NULL);
    end_copyConstruct(block_get3End(block), flower2);
    end_copyConstruct(block_get5End(block2), flower2);
    Group *group = group_construct(flower, flower2);
    Chain *chain = chain_construct(flower);
    link_construct(block_get3End(block), block_get5End(block2), group, chain);

    CuAssertTrue(testCase, block_getChain(block_construct(2, flower)) == NULL);
    CuAssertTrue(testCase, block_getChain(block) == chain);
    CuAssertTrue(testCase, block_getChain(block2) == chain);

    cactusBlockTestTeardown();
}

void testBlock_splitBlock(CuTest* testCase) {
    cactusBlockTestSetup();

    Block *leftBlock, *rightBlock;
    block_split(block, 2, &leftBlock, &rightBlock);
    CuAssertIntEquals(testCase, 2, block_getLength(leftBlock));
    CuAssertIntEquals(testCase, 1, block_getLength(rightBlock));
    CuAssertIntEquals(testCase, 3, block_getInstanceNumber(leftBlock));
    CuAssertIntEquals(testCase, 3, block_getInstanceNumber(rightBlock));

    Block *leftLeftBlock, *leftRightBlock;
    block_split(leftBlock, 1, &leftLeftBlock, &leftRightBlock);
    CuAssertIntEquals(testCase, 1, block_getLength(leftLeftBlock));
    CuAssertIntEquals(testCase, 1, block_getLength(leftRightBlock));
    CuAssertIntEquals(testCase, 3, block_getInstanceNumber(leftLeftBlock));
    CuAssertIntEquals(testCase, 3, block_getInstanceNumber(leftRightBlock));

    //doesn't currently check the instances.
    cactusBlockTestTeardown();
}

void testBlock_makeNewickString(CuTest *testCase) {
    assert(testCase != NULL);
    cactusBlockTestSetup();
    char *cA1 = block_makeNewickString(block, 1, 0);
    char *cA2 = block_makeNewickString(block, 1, 1);
    char *cA3 = block_makeNewickString(block, 0, 0);
    char *cA4 = block_makeNewickString(block, 0, 1);
    st_logDebug("I got the block tree string 1 0: %s\n", cA1);
    st_logDebug("I got the block tree string 1 1: %s\n", cA2);
    st_logDebug("I got the block tree string 0 0: %s\n", cA3);
    st_logDebug("I got the block tree string 0 1: %s\n", cA4);
    //CuAssertStrEquals(testCase, "(B,E)8;", cA);
    free(cA1);
    free(cA2);
    free(cA3);
    free(cA4);
    cactusBlockTestTeardown();
}

void testBlock_isTrivialChain(CuTest *testCase) {
    cactusBlockTestSetup();
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
    cactusBlockTestTeardown();
}

void testBlock_serialisation(CuTest* testCase) {
    cactusBlockTestSetup();
    Name rootInstanceName = segment_getName(rootSegment);
    Name leaf1InstanceName = segment_getName(leaf1Segment);
    Name leaf2InstanceName = segment_getName(leaf2Segment);
    int64_t i;
    void
            *vA =
                    binaryRepresentation_makeBinaryRepresentation(block,
                            (void(*)(void *, void(*)(const void *, size_t,
                                    size_t))) block_writeBinaryRepresentation,
                            &i);
    CuAssertTrue(testCase, i > 0);
    block_destruct(block);
    void *vA2 = vA;
    block = block_loadFromBinaryRepresentation(&vA2, flower);
    rootSegment
            = segment_getReverse(block_getInstance(block, rootInstanceName));
    leaf1Segment = block_getInstance(block, leaf1InstanceName);
    leaf2Segment = segment_getReverse(block_getInstance(block,
            leaf2InstanceName));
    free(vA);
    nestedTest = 1;
    testBlock_getName(testCase);
    testBlock_getOrientation(testCase);
    testBlock_getReverse(testCase);
    testBlock_getLength(testCase);
    testBlock_getFlower(testCase);
    testBlock_getLeftEnd(testCase);
    testBlock_getRightEnd(testCase);
    testBlock_getSetRootInstance(testCase);
    testBlock_getInstanceNumber(testCase);
    testBlock_getInstance(testCase);
    testBlock_getFirst(testCase);
    testBlock_instanceIterator(testCase);
    testBlock_makeNewickString(testCase);
    testBlock_getChain(testCase);
    nestedTest = 0;
    cactusBlockTestTeardown();
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
    SUITE_ADD_TEST(suite, testBlock_getSetRootInstance);
    SUITE_ADD_TEST(suite, testBlock_instanceIterator);
    SUITE_ADD_TEST(suite, testBlock_getChain);
    SUITE_ADD_TEST(suite, testBlock_splitBlock);
    SUITE_ADD_TEST(suite, testBlock_serialisation);
    SUITE_ADD_TEST(suite, testBlock_makeNewickString);
    SUITE_ADD_TEST(suite, testBlock_isTrivialChain);
    SUITE_ADD_TEST(suite, testBlock_construct);
    return suite;
}
