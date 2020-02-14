/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusEndsTestShared.h"

static bool nestedTest = 0;

static void cactusEndTestSetup(CuTest* testCase) {
    if (!nestedTest) {
        cactusEndsTestSharedSetup(testCase->name);
    }
}

static void cactusEndTestTeardown(CuTest* testCase) {
    if (!nestedTest) {
        cactusEndsTestSharedTeardown(testCase->name);
    }
}

void testEnd_construct(CuTest* testCase) {
    nestedTest = 0;
    cactusEndTestSetup(testCase);
    CuAssertTrue(testCase, end != NULL);
    cactusEndTestTeardown(testCase);
}

int64_t testEnd_copyConstructP(Event *event) {
    assert(event != NULL);
    return 1;
}

void testEnd_copyConstruct(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    Flower *flower2 = flower_construct(cactusDisk);
    sequence_construct(metaSequence, flower2);

    End *end2 = end_copyConstruct(end, flower2);
    CuAssertTrue(testCase, end_getName(end2) != NULL_NAME);
    CuAssertTrue(testCase, end_getName(end2) == end_getName(end));
    CuAssertTrue(testCase, flower_getEnd(flower2, end_getName(end2)) == end2);
    CuAssertTrue(testCase, cap_getName(end_getInstance(end2, cap_getName(rootCap))) == cap_getName(rootCap));
    CuAssertTrue(testCase, cap_getName(end_getInstance(end2, cap_getName(leaf1Cap))) == cap_getName(leaf1Cap));
    CuAssertTrue(testCase, cap_getName(end_getInstance(end2, cap_getName(leaf2Cap))) == cap_getName(leaf2Cap));
    cactusEndTestTeardown(testCase);
}

void testEnd_getName(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    CuAssertTrue(testCase, end_getName(end) != NULL_NAME);
    CuAssertTrue(testCase, flower_getEnd(flower, end_getName(end)) == end);
    CuAssertTrue(testCase, flower_getEnd(flower, end_getName(end_getReverse(end))) == end);
    cactusEndTestTeardown(testCase);
}

void testEnd_getOrientation(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    CuAssertTrue(testCase, end_getOrientation(end));
    CuAssertTrue(testCase, !end_getOrientation(end_getReverse(end)));
    cactusEndTestTeardown(testCase);
}

void testEnd_getReverse(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    CuAssertTrue(testCase, end_getReverse(end) != NULL);
    CuAssertTrue(testCase, end_getReverse(end_getReverse(end)) == end);
    cactusEndTestTeardown(testCase);
}

void testEnd_getSide(CuTest *testCase) {
    cactusEndTestSetup(testCase);
    Block *block = block_construct(10, flower);
    End *_5End = block_get5End(block);
    End *_3End = block_get3End(block);

    CuAssertTrue(testCase, end_getSide(end));
    CuAssertTrue(testCase, !end_getSide(end_getReverse(end)));

    CuAssertTrue(testCase, end_getSide(_5End));
    CuAssertTrue(testCase, !end_getSide(_3End));
    CuAssertTrue(testCase, end_getSide(end_getReverse(_3End)));
    CuAssertTrue(testCase, !end_getSide(end_getReverse(_5End)));

    cactusEndTestTeardown(testCase);
}

void testEnd_getFlower(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    CuAssertTrue(testCase, end_getFlower(end) == flower);
    cactusEndTestTeardown(testCase);
}

void testEnd_getBlock(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    Block *block = block_construct(10, flower);
    End *leftEnd = block_get5End(block);
    End *rightEnd = block_get3End(block);

    CuAssertTrue(testCase, end_getBlock(end) == NULL);
    CuAssertTrue(testCase, end_getBlock(end_getReverse(end)) == NULL);

    CuAssertTrue(testCase, end_getBlock(leftEnd) == block);
    CuAssertTrue(testCase, end_getBlock(end_getReverse(leftEnd)) == block_getReverse(block));
    CuAssertTrue(testCase, block_getOrientation(block) == end_getOrientation(leftEnd));

    CuAssertTrue(testCase, end_getBlock(rightEnd) == block);
    CuAssertTrue(testCase, end_getBlock(end_getReverse(rightEnd)) == block_getReverse(block));
    CuAssertTrue(testCase, block_getOrientation(block) == end_getOrientation(rightEnd));

    cactusEndTestTeardown(testCase);
}

void testEnd_getOtherBlockEnd(CuTest *testCase) {
    cactusEndTestSetup(testCase);
    Block *block = block_construct(10, flower);
    End *leftEnd = block_get5End(block);
    End *rightEnd = block_get3End(block);

    CuAssertTrue(testCase, end_getOtherBlockEnd(end) == NULL);
    CuAssertTrue(testCase, end_getOtherBlockEnd(end_getReverse(end)) == NULL);

    CuAssertTrue(testCase, end_getOtherBlockEnd(leftEnd) == rightEnd);
    CuAssertTrue(testCase, end_getOtherBlockEnd(rightEnd) == leftEnd);

    CuAssertTrue(testCase, end_getOtherBlockEnd(end_getReverse(leftEnd)) == end_getReverse(rightEnd));
    CuAssertTrue(testCase, end_getOtherBlockEnd(end_getReverse(rightEnd)) == end_getReverse(leftEnd));

    cactusEndTestTeardown(testCase);
}

void testEnd_getGroup(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    Flower *flower2 = flower_construct(cactusDisk);
    sequence_construct(metaSequence, flower2);
    End *end2 = end_copyConstruct(end, flower2);
    CuAssertTrue(testCase, end_getGroup(end) == NULL);
    Group *group = group_construct(flower, flower2);
    CuAssertTrue(testCase, end_getGroup(end) == group);
    CuAssertTrue(testCase, end_getGroup(end2) == NULL);
    cactusEndTestTeardown(testCase);
}

void testEnd_setGroup(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    Flower *flower2 = flower_construct(cactusDisk);
    Group *group2 = group_construct2(flower2);
    End *end2 = end_construct(1, flower2);
    End *end3 = end_construct(1, flower2);
    CuAssertTrue(testCase, group_getEndNumber(group2) == 0);
    CuAssertTrue(testCase, end_getGroup(end2) == NULL);
    CuAssertTrue(testCase, end_getGroup(end3) == NULL);
    end_setGroup(end2, group2);
    CuAssertTrue(testCase, group_getEndNumber(group2) == 1);
    CuAssertTrue(testCase, end_getGroup(end2) == group2);
    CuAssertTrue(testCase, group_getEnd(group2, end_getName(end2)) == end2);
    CuAssertTrue(testCase, end_getGroup(end3) == NULL);
    end_setGroup(end3, group2);
    CuAssertTrue(testCase, group_getEndNumber(group2) == 2);
    CuAssertTrue(testCase, end_getGroup(end2) == group2);
    CuAssertTrue(testCase, group_getEnd(group2, end_getName(end2)) == end2);
    CuAssertTrue(testCase, end_getGroup(end3) == group2);
    CuAssertTrue(testCase, group_getEnd(group2, end_getName(end3)) == end3);
    end_setGroup(end3, NULL);
    end_setGroup(end2, group2);
    CuAssertTrue(testCase, group_getEndNumber(group2) == 1);
    CuAssertTrue(testCase, end_getGroup(end2) == group2);
    CuAssertTrue(testCase, group_getEnd(group2, end_getName(end2)) == end2);
    CuAssertTrue(testCase, end_getGroup(end3) == NULL);
    cactusEndTestTeardown(testCase);
}

void testEnd_getInstanceNumber(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    End *end2 = end_construct(0, flower);
    CuAssertTrue(testCase, end_getInstanceNumber(end2) == 0);
    CuAssertTrue(testCase, end_getInstanceNumber(end) == 4);
    cactusEndTestTeardown(testCase);
}

void testEnd_getInstance(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    CuAssertTrue(testCase, end_getInstance(end, cap_getName(rootCap)) == cap_getReverse(rootCap));
    CuAssertTrue(testCase, end_getInstance(end, cap_getName(leaf1Cap)) == cap_getReverse(leaf1Cap));
    CuAssertTrue(testCase, end_getInstance(end, cap_getName(leaf2Cap)) == leaf2Cap);
    cactusEndTestTeardown(testCase);
}

void testEnd_getFirst(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    CuAssertTrue(testCase, end_getFirst(end) == cap_getReverse(rootCap));
    cactusEndTestTeardown(testCase);
}

void testEnd_getSetRootInstance(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    CuAssertTrue(testCase, end_getRootInstance(end) == cap_getReverse(rootCap));
    CuAssertTrue(testCase, end_getRootInstance(end_getReverse(end)) == rootCap);

    End *end2 = end_construct(0, flower);
    CuAssertTrue(testCase, end_getRootInstance(end2) == NULL);
    CuAssertTrue(testCase, end_getRootInstance(end_getReverse(end2)) == NULL);
    cactusEndTestTeardown(testCase);
}

void testEnd_instanceIterator(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    End_InstanceIterator *iterator = end_getInstanceIterator(end);
    CuAssertTrue(testCase, iterator != NULL);
    CuAssertTrue(testCase, end_getNext(iterator) == cap_getReverse(rootCap));
    CuAssertTrue(testCase, end_getNext(iterator) == cap_getReverse(leaf1Cap));

    End_InstanceIterator *iterator2 = end_copyInstanceIterator(iterator);

    CuAssertTrue(testCase, end_getNext(iterator) == leaf2Cap);
    CuAssertTrue(testCase, end_getNext(iterator) == cap_getReverse(leaf3Cap));
    CuAssertTrue(testCase, end_getNext(iterator) == NULL);
    CuAssertTrue(testCase, end_getPrevious(iterator) == cap_getReverse(leaf3Cap));
    CuAssertTrue(testCase, end_getPrevious(iterator) == leaf2Cap);
    CuAssertTrue(testCase, end_getPrevious(iterator) == cap_getReverse(leaf1Cap));
    CuAssertTrue(testCase, end_getPrevious(iterator) == cap_getReverse(rootCap));
    CuAssertTrue(testCase, end_getPrevious(iterator) == NULL);

    CuAssertTrue(testCase, end_getNext(iterator2) == leaf2Cap);
    CuAssertTrue(testCase, end_getNext(iterator2) == cap_getReverse(leaf3Cap));
    CuAssertTrue(testCase, end_getNext(iterator2) == NULL);
    CuAssertTrue(testCase, end_getPrevious(iterator2) == cap_getReverse(leaf3Cap));
    CuAssertTrue(testCase, end_getPrevious(iterator2) == leaf2Cap);
    CuAssertTrue(testCase, end_getPrevious(iterator2) == cap_getReverse(leaf1Cap));
    CuAssertTrue(testCase, end_getPrevious(iterator2) == cap_getReverse(rootCap));
    CuAssertTrue(testCase, end_getPrevious(iterator2) == NULL);

    end_destructInstanceIterator(iterator);
    end_destructInstanceIterator(iterator2);

    iterator = end_getInstanceIterator(end_getReverse(end));
    CuAssertTrue(testCase, end_getNext(iterator) == rootCap);
    CuAssertTrue(testCase, end_getNext(iterator) == leaf1Cap);
    CuAssertTrue(testCase, end_getNext(iterator) == cap_getReverse(leaf2Cap));
    CuAssertTrue(testCase, end_getNext(iterator) == leaf3Cap);
    CuAssertTrue(testCase, end_getNext(iterator) == NULL);
    CuAssertTrue(testCase, end_getPrevious(iterator) == leaf3Cap);
    CuAssertTrue(testCase, end_getPrevious(iterator) == cap_getReverse(leaf2Cap));
    CuAssertTrue(testCase, end_getPrevious(iterator) == leaf1Cap);
    CuAssertTrue(testCase, end_getPrevious(iterator) == rootCap);
    CuAssertTrue(testCase, end_getPrevious(iterator) == NULL);

    end_destructInstanceIterator(iterator);

    cactusEndTestTeardown(testCase);
}

void testEnd_isBlockOrStubEnd(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    CuAssertTrue(testCase, end_isStubEnd(end));
    CuAssertTrue(testCase, !end_isBlockEnd(end));
    Block *block = block_construct(2, flower);
    End *leftEnd = block_get5End(block);
    End *rightEnd = block_get3End(block);
    CuAssertTrue(testCase, end_isBlockEnd(leftEnd));
    CuAssertTrue(testCase, end_isBlockEnd(rightEnd));
    CuAssertTrue(testCase, !end_isStubEnd(leftEnd));
    CuAssertTrue(testCase, !end_isStubEnd(rightEnd));
    cactusEndTestTeardown(testCase);
}

void testEnd_isAttachedOrFree(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    End *end2 = end_construct(1, flower);
    End *end3 = end_construct(0, flower);
    Block *block = block_construct(2, flower);
    End *end4 = block_get5End(block);
    End *end5 = block_get3End(block);
    CuAssertTrue(testCase, end_isAttached(end2));
    CuAssertTrue(testCase, !end_isAttached(end3));
    CuAssertTrue(testCase, !end_isFree(end2));
    CuAssertTrue(testCase, end_isFree(end3));

    CuAssertTrue(testCase, !end_isAttached(end4));
    CuAssertTrue(testCase, !end_isAttached(end5));
    CuAssertTrue(testCase, end_isFree(end4));
    CuAssertTrue(testCase, end_isFree(end5));
    cactusEndTestTeardown(testCase);
}

void testEnd_getCapForEvent(CuTest* testCase) {
    cactusEndTestSetup(testCase);

    CuAssertPtrEquals(testCase, end_getCapForEvent(end_getReverse(end), event_getName(rootEvent)), rootCap);
    Cap *cap = end_getCapForEvent(end, event_getName(leafEvent));
    CuAssertTrue(testCase, cap == cap_getReverse(leaf1Cap) || cap == leaf2Cap || cap == cap_getReverse(leaf3Cap));
    CuAssertTrue(testCase, end_getCapForEvent(end, NULL_NAME) == NULL);

    cactusEndTestTeardown(testCase);
}

void testEnd_serialisation(CuTest* testCase) {
    cactusEndTestSetup(testCase);
    Name rootInstanceName = cap_getName(rootCap);
    Name leaf1InstanceName = cap_getName(leaf1Cap);
    Name leaf2InstanceName = cap_getName(leaf2Cap);
    Name leaf3InstanceName = cap_getName(leaf3Cap);
    int64_t i;
    void *vA = binaryRepresentation_makeBinaryRepresentation(end,
            (void(*)(void *, void(*)(const void *, size_t, size_t))) end_writeBinaryRepresentation, &i);
    CuAssertTrue(testCase, i > 0);
    end_destruct(end);
    void *vA2 = vA;
    end = end_loadFromBinaryRepresentation(&vA2, flower);
    rootCap = cap_getReverse(end_getInstance(end, rootInstanceName));
    leaf1Cap = cap_getReverse(end_getInstance(end, leaf1InstanceName));
    leaf2Cap = end_getInstance(end, leaf2InstanceName);
    leaf3Cap = cap_getReverse(end_getInstance(end, leaf3InstanceName));
    CuAssertTrue(testCase, leaf3Cap != NULL);
    free(vA);
    nestedTest = 1;
    testEnd_copyConstruct(testCase);
    testEnd_getName(testCase);
    testEnd_getOrientation(testCase);
    testEnd_getReverse(testCase);
    testEnd_getSide(testCase);
    testEnd_getFlower(testCase);
    testEnd_getBlock(testCase);
    testEnd_getOtherBlockEnd(testCase);
    testEnd_getGroup(testCase);
    testEnd_setGroup(testCase);
    testEnd_getInstanceNumber(testCase);
    testEnd_getInstance(testCase);
    testEnd_getFirst(testCase);
    testEnd_getSetRootInstance(testCase);
    testEnd_instanceIterator(testCase);
    testEnd_isBlockOrStubEnd(testCase);
    testEnd_isAttachedOrFree(testCase);
    testEnd_getCapForEvent(testCase);
    nestedTest = 0;
    cactusEndTestTeardown(testCase);
}

CuSuite* cactusEndTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testEnd_copyConstruct);
    SUITE_ADD_TEST(suite, testEnd_getName);
    SUITE_ADD_TEST(suite, testEnd_getOrientation);
    SUITE_ADD_TEST(suite, testEnd_getReverse);
    SUITE_ADD_TEST(suite, testEnd_getSide);
    SUITE_ADD_TEST(suite, testEnd_getFlower);
    SUITE_ADD_TEST(suite, testEnd_getBlock);
    SUITE_ADD_TEST(suite, testEnd_getOtherBlockEnd);
    SUITE_ADD_TEST(suite, testEnd_getGroup);
    SUITE_ADD_TEST(suite, testEnd_setGroup);
    SUITE_ADD_TEST(suite, testEnd_getInstanceNumber);
    SUITE_ADD_TEST(suite, testEnd_getInstance);
    SUITE_ADD_TEST(suite, testEnd_getFirst);
    SUITE_ADD_TEST(suite, testEnd_getSetRootInstance);
    SUITE_ADD_TEST(suite, testEnd_instanceIterator);
    SUITE_ADD_TEST(suite, testEnd_isBlockOrStubEnd);
    SUITE_ADD_TEST(suite, testEnd_isAttachedOrFree);
    SUITE_ADD_TEST(suite, testEnd_serialisation);
    SUITE_ADD_TEST(suite, testEnd_getCapForEvent);
    SUITE_ADD_TEST(suite, testEnd_construct);
    return suite;
}
