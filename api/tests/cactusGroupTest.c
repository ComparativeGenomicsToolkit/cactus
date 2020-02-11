/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk;
static Flower *flower;
static Flower *nestedFlower;
static End *end1;
static End *end2;
static End *end3;
static End *end4;
static End *nestedEnd1;
static End *nestedEnd2;
static Group *group, *group2;

static void cactusGroupTestTeardown(CuTest* testCase) {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
        cactusDisk = NULL;
    }
}

static void cactusGroupTestSetup(CuTest* testCase) {
    cactusGroupTestTeardown(testCase);
    cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
    flower = flower_construct(cactusDisk);
    nestedFlower = flower_construct(cactusDisk);
    end1 = end_construct2(0, 0, flower);
    end2 = end_construct(0, flower);
    end3 = end_construct(0, flower);
    nestedEnd1 = end_copyConstruct(end1, nestedFlower);
    nestedEnd2 = end_copyConstruct(end2, nestedFlower);
    group = group_construct(flower, nestedFlower);
    group2 = group_construct2(flower);
    end4 = end_construct(0, flower);
}

void testGroup_construct(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group != NULL);
    cactusGroupTestTeardown(testCase);
}

void testGroup_updateContainedEnds(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    end_copyConstruct(end3, nestedFlower);
    CuAssertTrue(testCase, group_getEndNumber(group) == 2);
    group_updateContainedEnds(group);
    CuAssertTrue(testCase, group_getEndNumber(group) == 3);
    CuAssertTrue(testCase, group_getEnd(group, end_getName(end1)) == end1);
    CuAssertTrue(testCase, group_getEnd(group, end_getName(end2)) == end2);
    CuAssertTrue(testCase, group_getEnd(group, end_getName(end3)) == end3);
    cactusGroupTestTeardown(testCase);
}

void testGroup_makeNonLeaf(CuTest *testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_isLeaf(group2));
    end_setGroup(end4, group2);
    group_makeNestedFlower(group2);
    CuAssertTrue(testCase, !group_isLeaf(group2));
    Flower *nestedFlower = group_getNestedFlower(group2);
    CuAssertTrue(testCase, nestedFlower != NULL);
    CuAssertTrue(testCase, !flower_builtBlocks(flower));
    CuAssertTrue(testCase, !flower_builtTrees(flower));
    CuAssertTrue(testCase, !flower_builtFaces(flower));
    CuAssertTrue(testCase, flower_getName(nestedFlower) == group_getName(group2));
    CuAssertTrue(testCase, flower_getParentGroup(nestedFlower) == group2);
    CuAssertTrue(testCase, flower_getEndNumber(nestedFlower) == 1);
    End *nestedEnd = flower_getFirstEnd(nestedFlower);
    CuAssertTrue(testCase, end_getName(end4) == end_getName(nestedEnd));
    CuAssertTrue(testCase, end_getGroup(nestedEnd) != NULL);
    CuAssertTrue(testCase, flower_getGroupNumber(nestedFlower) == 1);
    CuAssertTrue(testCase, flower_isTerminal(nestedFlower));
    cactusGroupTestTeardown(testCase);
}

void testGroup_addEnd(CuTest *testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_getEndNumber(group2) == 0);
    end_setGroup(end4, group2);
    CuAssertTrue(testCase, group_getEndNumber(group2) == 1);
    CuAssertTrue(testCase, end_getGroup(end4) == group2);
    CuAssertTrue(testCase, group_getEnd(group2, end_getName(end4)) == end4);
    cactusGroupTestTeardown(testCase);
}

void testGroup_isLeaf(CuTest *testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, !group_isLeaf(group));
    CuAssertTrue(testCase, group_isLeaf(group2));
    cactusGroupTestTeardown(testCase);
}

static Chain *setupChain() {
    Chain *chain = chain_construct(flower);
    end_setGroup(end1, group2);
    end_setGroup(end2, group2);
    link_construct(end1, end2, group2, chain);
    return chain;
}

void testGroup_getLink(CuTest *testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_getLink(group2) == NULL);
    Chain *chain = setupChain();
    CuAssertTrue(testCase, group_getLink(group2) == chain_getFirst(chain));
    chain_destruct(chain);
    CuAssertTrue(testCase, group_getLink(group2) == NULL);
    cactusGroupTestTeardown(testCase);
}

void testGroup_isTangle(CuTest *testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_isTangle(group2));
    Chain *chain = setupChain();
    CuAssertTrue(testCase, !group_isTangle(group2));
    chain_destruct(chain);
    CuAssertTrue(testCase, group_isTangle(group2));
    cactusGroupTestTeardown(testCase);
}

void testGroup_isLink(CuTest *testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, !group_isLink(group2));
    Chain *chain = setupChain();
    CuAssertTrue(testCase, group_isLink(group2));
    chain_destruct(chain);
    CuAssertTrue(testCase, !group_isLink(group2));
    cactusGroupTestTeardown(testCase);
}

void testGroup_getFlower(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_getFlower(group) == flower);
    cactusGroupTestTeardown(testCase);
}

void testGroup_getNestedFlowerName(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_getName(group) == flower_getName(nestedFlower));
    cactusGroupTestTeardown(testCase);
}

void testGroup_getNestedFlower(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_getNestedFlower(group) == nestedFlower);
    cactusGroupTestTeardown(testCase);
}

void testGroup_getChain(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_getLink(group) == NULL);
    Chain *chain = chain_construct(flower);
    link_construct(end1, end2, group, chain);
    CuAssertTrue(testCase, group_getLink(group) == chain_getFirst(chain));
    cactusGroupTestTeardown(testCase);
}

void testGroup_getEnd(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_getEnd(group, end_getName(end1)) == end1);
    CuAssertTrue(testCase, group_getEnd(group, end_getName(end2)) == end2);
    cactusGroupTestTeardown(testCase);
}

void testGroup_getEndNumber(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_getEndNumber(group) == 2);
    cactusGroupTestTeardown(testCase);
}

void testGroup_getAttachedStubAndBlockEndNumber(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    CuAssertTrue(testCase, group_getAttachedStubEndNumber(group) == 0);
    CuAssertTrue(testCase, group_getBlockEndNumber(group) == 0);
    CuAssertIntEquals(testCase, 2, group_getStubEndNumber(group));
    CuAssertIntEquals(testCase, 2, group_getFreeStubEndNumber(group));
    end_setGroup(end_construct(1, flower), group);
    end_setGroup(end_construct(0, flower), group);
    end_setGroup(end_construct(1, flower), group);
    CuAssertTrue(testCase, group_getAttachedStubEndNumber(group) == 2);
    CuAssertTrue(testCase, group_getBlockEndNumber(group) == 0);
    CuAssertTrue(testCase, group_getStubEndNumber(group) == 5);
    CuAssertTrue(testCase, group_getFreeStubEndNumber(group) == 3);
    Block *block = block_construct(1, flower);
    end_setGroup(block_get5End(block), group);
    end_setGroup(block_get3End(block), group);
    CuAssertTrue(testCase, group_getAttachedStubEndNumber(group) == 2);
    CuAssertTrue(testCase, group_getBlockEndNumber(group) == 2);
    CuAssertTrue(testCase, group_getStubEndNumber(group) == 5);
    CuAssertTrue(testCase, group_getFreeStubEndNumber(group) == 3);
    cactusGroupTestTeardown(testCase);
}

void testGroup_endIterator(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    Group_EndIterator *iterator = group_getEndIterator(group);

    CuAssertTrue(testCase, group_getNextEnd(iterator) == end1);
    CuAssertTrue(testCase, group_getNextEnd(iterator) == end2);
    CuAssertTrue(testCase, group_getNextEnd(iterator) == NULL);

    Group_EndIterator *iterator2 = group_copyEndIterator(iterator);

    CuAssertTrue(testCase, group_getPreviousEnd(iterator) == end2);
    CuAssertTrue(testCase, group_getPreviousEnd(iterator) == end1);
    CuAssertTrue(testCase, group_getPreviousEnd(iterator) == NULL);

    group_destructEndIterator(iterator);

    CuAssertTrue(testCase, group_getPreviousEnd(iterator2) == end2);
    CuAssertTrue(testCase, group_getPreviousEnd(iterator2) == end1);
    CuAssertTrue(testCase, group_getPreviousEnd(iterator2) == NULL);

    group_destructEndIterator(iterator2);

    cactusGroupTestTeardown(testCase);
}

void testGroup_getTotalBaseLength(CuTest *testCase) {
    cactusGroupTestSetup(testCase);

    CuAssertTrue(testCase, group_getTotalBaseLength(group) == 0);
    //cap_construct

    cactusGroupTestTeardown(testCase);
}

void testGroup_constructChainForLink(CuTest *testCase) {
    cactusGroupTestSetup(testCase);
    //Create a link group and test function works!
    Group *group3 = group_construct2(flower);
    End *_5End = end_construct2(1, 1, flower);
    End *_3End = end_construct2(0, 1, flower);
    End *stubEnd = end_construct2(1, 0, flower);
    end_setGroup(_5End, group3);
    end_setGroup(_3End, group3);
    end_setGroup(stubEnd, group3);
    CuAssertTrue(testCase, group_isTangle(group3));
    CuAssertTrue(testCase, group_getEndNumber(group3) == 3);
    group_constructChainForLink(group3);
    CuAssertTrue(testCase, group_isLink(group3));
    Link *link = group_getLink(group3);
    CuAssertTrue(testCase, link_get5End(link) == _5End);
    CuAssertTrue(testCase, link_get3End(link) == _3End);
    cactusGroupTestTeardown(testCase);
}

#if 0
void testGroup_mergeGroups(CuTest *testCase) {
 cactusGroupTestSetup(testCase);

 CuAssertTrue(testCase, group_getTotalBaseLength(group) == 0);
 //cap_construct
 Flower *parentFlower = flower_construct(cactusDisk);
 End *end1 = end_construct(0, parentFlower);
 End *end2 = end_construct(0, parentFlower);
 End *end3 = end_construct(0, parentFlower);
 End *end4 = end_construct(0, parentFlower);

 Group *parentGroup1 = group_construct2(parentFlower);
 Group *parentGroup2 = group_construct2(parentFlower);
 Group *parentGroup3 = group_construct2(parentFlower);
 Group *parentGroup4 = group_construct2(parentFlower);

 end_setGroup(end1, parentGroup1);
 end_setGroup(end2, parentGroup1);
 end_setGroup(end3, parentGroup2);
 end_setGroup(end4, parentGroup3);

 group_makeNonLeaf(parentGroup3);
 group_makeNonLeaf(parentGroup4);

 Group *mergedGroup1 = group_mergeGroups(parentGroup1, parentGroup2); //merge two leaf flowers
 Group *mergedGroup2 = group_mergeGroups(mergedGroup1, parentGroup3); //merge a leaf and non-leaf flower
 Group *finalGroup = group_mergeGroups(parentGroup4, mergedGroup2); //merge two non-leaf flowers

 CuAssertTrue(testCase, flower_getGroupNumber(parentFlower) == 1);
 CuAssertTrue(testCase, flower_getFirstGroup(parentFlower) == finalGroup);
 CuAssertTrue(testCase, flower_getEndNumber(parentFlower) == 4);
 CuAssertTrue(testCase, end_getGroup(end1) == finalGroup);
 CuAssertTrue(testCase, end_getGroup(end2) == finalGroup);
 CuAssertTrue(testCase, end_getGroup(end3) == finalGroup);
 CuAssertTrue(testCase, end_getGroup(end4) == finalGroup);

 CuAssertTrue(testCase, !group_isLeaf(finalGroup));
 Flower *nestedFlower = group_getNestedFlower(finalGroup);
 CuAssertTrue(testCase, flower_getEndNumber(nestedFlower) == 4);
 CuAssertTrue(testCase, flower_getEnd(nestedFlower, end_getName(end1)) != NULL);
 CuAssertTrue(testCase, flower_getEnd(nestedFlower, end_getName(end2)) != NULL);
 CuAssertTrue(testCase, flower_getEnd(nestedFlower, end_getName(end3)) != NULL);
 CuAssertTrue(testCase, flower_getEnd(nestedFlower, end_getName(end4)) != NULL);

 cactusGroupTestTeardown(testCase);
 }
#endif

void testGroup_serialisation(CuTest* testCase) {
    cactusGroupTestSetup(testCase);
    int64_t i;
    Name name = group_getName(group);
    void
            *vA =
                    binaryRepresentation_makeBinaryRepresentation(group,
                            (void(*)(void *, void(*)(const void *, size_t,
                                    size_t))) group_writeBinaryRepresentation,
                            &i);
    CuAssertTrue(testCase, i > 0);
    group_destruct(group);
    void *vA2 = vA;
    group = group_loadFromBinaryRepresentation(&vA2, flower);
    free(vA);
    CuAssertTrue(testCase, group_getName(group) == name);
    CuAssertTrue(testCase, group_getFlower(group) == flower);
    CuAssertTrue(testCase, group_getNestedFlower(group) == nestedFlower);
    Group_EndIterator *iterator = group_getEndIterator(group);
    CuAssertTrue(testCase, group_getNextEnd(iterator) == end1);
    CuAssertTrue(testCase, group_getNextEnd(iterator) == end2);
    CuAssertTrue(testCase, group_getNextEnd(iterator) == NULL);
    group_destructEndIterator(iterator);

    cactusGroupTestTeardown(testCase);
}

CuSuite* cactusGroupTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testGroup_updateContainedEnds);
    SUITE_ADD_TEST(suite, testGroup_addEnd);
    SUITE_ADD_TEST(suite, testGroup_isLeaf);
    SUITE_ADD_TEST(suite, testGroup_makeNonLeaf);
    SUITE_ADD_TEST(suite, testGroup_getLink);
    SUITE_ADD_TEST(suite, testGroup_isLink);
    SUITE_ADD_TEST(suite, testGroup_isTangle);
    SUITE_ADD_TEST(suite, testGroup_getFlower);
    SUITE_ADD_TEST(suite, testGroup_getNestedFlowerName);
    SUITE_ADD_TEST(suite, testGroup_getNestedFlower);
    SUITE_ADD_TEST(suite, testGroup_getChain);
    SUITE_ADD_TEST(suite, testGroup_getEnd);
    SUITE_ADD_TEST(suite, testGroup_getEndNumber);
    SUITE_ADD_TEST(suite, testGroup_getAttachedStubAndBlockEndNumber);
    SUITE_ADD_TEST(suite, testGroup_endIterator);
    SUITE_ADD_TEST(suite, testGroup_getTotalBaseLength);
    SUITE_ADD_TEST(suite, testGroup_constructChainForLink);
    //SUITE_ADD_TEST(suite, testGroup_mergeGroups);
    SUITE_ADD_TEST(suite, testGroup_serialisation);
    SUITE_ADD_TEST(suite, testGroup_construct);
    return suite;
}
