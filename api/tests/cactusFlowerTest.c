/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

/*
 * Global variables for test.
 */

static CactusDisk *cactusDisk = NULL;
static Flower *flower;
static EventTree *eventTree;
static Sequence *sequence;
static Sequence *sequence2;
static Sequence *sequence;
static Sequence *sequence2;
static End *end;
static End *end2;
static Block *block;
static Block *block2;
static Group *group;
static Group *group2;
static Chain *chain;
static Chain *chain2;
static Segment *segment;
static Segment *segment2;
static Cap *cap;
static Cap *cap2;

/*
 * Setup/teardown functions.
 */

static void cactusFlowerTestTeardown(CuTest* testCase) {
    if (cactusDisk != NULL) {
        cactusDisk_destruct(cactusDisk);
        cactusDisk = NULL;
        flower = NULL;
        eventTree = NULL;
        sequence = NULL;
        sequence = NULL;
    }
}

static void cactusFlowerTestSetup(CuTest* testCase) {
    cactusFlowerTestTeardown(testCase);
    cactusDisk = cactusDisk_construct();
    flower = flower_construct(cactusDisk);
    eventTree = eventTree_construct2(cactusDisk);
}

static void sequenceSetup() {
    sequence = sequence_construct(0, 10, "ACTGACTGAC", ">one",
            eventTree_getRootEvent(eventTree), cactusDisk);
    flower_addSequence(flower, sequence);
    sequence2 = sequence_construct(0, 10, "ACTGACTGAC", ">two",
            eventTree_getRootEvent(eventTree), cactusDisk);
    flower_addSequence(flower, sequence2);
}

static void endsSetup() {
    end = end_construct(1, flower);
    end2 = end_construct(1, flower);
}

static void capsSetup() {
    endsSetup();
    cap = cap_construct(end, eventTree_getRootEvent(eventTree));
    cap2 = cap_construct(end2, eventTree_getRootEvent(eventTree));
}

static void blocksSetup() {
    block = block_construct(1, flower);
    block2 = block_construct(2, flower);
}

void segmentsSetup() {
    blocksSetup();
    segment = segment_construct(block, eventTree_getRootEvent(eventTree));
    segment2 = segment_construct(block2, eventTree_getRootEvent(eventTree));
}

static void chainsSetup() {
    chain = chain_construct(flower);
    chain2 = chain_construct(flower);
}

static void groupsSetup() {
    group = group_construct(flower, flower_construct(cactusDisk));
    group2 = group_construct(flower, flower_construct(cactusDisk));
}

/*
 * This tests all the retrieval functions for each type of object.
 */
static void testObjectRetrieval(CuTest* testCase,
        int64_t(*getObjectNumberFn)(Flower *flower), void *(*getFirstObjectFn)(
                Flower *flower), Name(*objectGetNameFn)(void *),
        void *(*getObjectFn)(Flower *flower, Name name),
        void *(*constructIterator)(Flower *flower), void(*destructIterator)(
                void *iterator), void *(*getNext)(void *iterator),
        void *(*getPrevious)(void *iterator), void *(*copyIterator)(
                void *iterator), void *object, void *object2) {
    /*
     * Test number function
     */
    CuAssertTrue(testCase, getObjectNumberFn(flower) == 2);

    /*
     * Test get first function.
     */
    CuAssertTrue(testCase, getFirstObjectFn(flower) == object2);

    /*
     * Test get function
     */
    if (objectGetNameFn != NULL) {
        CuAssertTrue(testCase, getObjectFn(flower, objectGetNameFn(object)) == object);
        CuAssertTrue(testCase, getObjectFn(flower, objectGetNameFn(object2)) == object2);
    } else {
        assert(getObjectFn == NULL);
    }

    /*
     * Test iterator.
     */
    void *iterator = constructIterator(flower);
    CuAssertTrue(testCase, getNext(iterator) == object);
    CuAssertTrue(testCase, getNext(iterator) == object2);
    CuAssertTrue(testCase, getNext(iterator) == NULL);
    void *iterator2 = copyIterator(iterator);
    CuAssertTrue(testCase, getPrevious(iterator) == object2);
    CuAssertTrue(testCase, getPrevious(iterator) == object);
    CuAssertTrue(testCase, getPrevious(iterator) == NULL);
    destructIterator(iterator);
    CuAssertTrue(testCase, getPrevious(iterator2) == object2);
    CuAssertTrue(testCase, getPrevious(iterator2) == object);
    CuAssertTrue(testCase, getPrevious(iterator2) == NULL);
    destructIterator(iterator2);
}

/*
 * Now all the actual tests.
 */

void testFlower_constructAndDestruct(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    CuAssertTrue(testCase, flower != NULL);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_getName(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    CuAssertTrue(testCase, flower_getName(flower) != NULL_NAME);
    CuAssertTrue(testCase, cactusDisk_getFlower(cactusDisk, flower_getName(flower)) == flower);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_getCactusDisk(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    CuAssertTrue(testCase, flower_getCactusDisk(flower) == cactusDisk);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_getEventTree(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    CuAssertTrue(testCase, flower_getEventTree(flower) == eventTree);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_sequence(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    sequenceSetup();
    testObjectRetrieval(testCase,
            (int64_t(*)(Flower *flower)) flower_getSequenceNumber, (void *(*)(
                    Flower *)) flower_getFirstSequence,
            (Name(*)(void *)) sequence_getName,
            (void *(*)(Flower *, Name name)) flower_getSequence, (void *(*)(
                    Flower *)) flower_getSequenceIterator,
            (void(*)(void *)) flower_destructSequenceIterator, (void *(*)(
                    void *)) flower_getNextSequence,
            (void *(*)(void *)) flower_getPreviousSequence,
            (void *(*)(void *)) flower_copySequenceIterator,
            sequence, sequence2);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_cap(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    capsSetup();
    testObjectRetrieval(testCase,
            (int64_t(*)(Flower *flower)) flower_getCapNumber, (void *(*)(
                    Flower *)) flower_getFirstCap,
            (Name(*)(void *)) cap_getName,
            (void *(*)(Flower *, Name name)) flower_getCap,
            (void *(*)(Flower *)) flower_getCapIterator,
            (void(*)(void *)) flower_destructCapIterator,
            (void *(*)(void *)) flower_getNextCap,
            (void *(*)(void *)) flower_getPreviousCap,
            (void *(*)(void *)) flower_copyCapIterator, cap,
            cap2);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_end(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    endsSetup();
    testObjectRetrieval(testCase,
            (int64_t(*)(Flower *flower)) flower_getEndNumber, (void *(*)(
                    Flower *)) flower_getFirstEnd,
            (Name(*)(void *)) end_getName,
            (void *(*)(Flower *, Name name)) flower_getEnd,
            (void *(*)(Flower *)) flower_getEndIterator,
            (void(*)(void *)) flower_destructEndIterator,
            (void *(*)(void *)) flower_getNextEnd,
            (void *(*)(void *)) flower_getPreviousEnd,
            (void *(*)(void *)) flower_copyEndIterator, end,
            end2);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_chain(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    chainsSetup();
    testObjectRetrieval(testCase,
            (int64_t(*)(Flower *flower)) flower_getChainNumber, (void *(*)(
                    Flower *)) flower_getFirstChain,
            (Name(*)(void *)) chain_getName,
            (void *(*)(Flower *, Name name)) flower_getChain, (void *(*)(
                    Flower *)) flower_getChainIterator,
            (void(*)(void *)) flower_destructChainIterator,
            (void *(*)(void *)) flower_getNextChain,
            (void *(*)(void *)) flower_getPreviousChain,
            (void *(*)(void *)) flower_copyChainIterator, chain,
            chain2);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_getTrivialChainNumber(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    CuAssertIntEquals(testCase, 0, flower_getTrivialChainNumber(flower));
    chain = chain_construct(flower);
    chain2 = chain_construct(flower);
    CuAssertIntEquals(testCase, 0, flower_getTrivialChainNumber(flower));
    group = group_construct2(flower);
    group2 = group_construct2(flower);
    block = block_construct(1, flower);
    end_setGroup(block_get5End(block), group);
    end_setGroup(block_get3End(block), group2);
    CuAssertIntEquals(testCase, 1, flower_getTrivialChainNumber(flower));
    block2 = block_construct(1, flower);
    end_setGroup(block_get5End(block2), group2);
    end_setGroup(block_get3End(block2), group);
    CuAssertIntEquals(testCase, 2, flower_getTrivialChainNumber(flower));
    link_construct(block_get3End(block), block_get5End(block2), group2, chain);
    CuAssertIntEquals(testCase, 0, flower_getTrivialChainNumber(flower));
    link_construct(block_get3End(block2), block_get5End(block), group, chain);
    CuAssertIntEquals(testCase, 0, flower_getTrivialChainNumber(flower));
    cactusFlowerTestTeardown(testCase);
}

void testFlower_group(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    groupsSetup();
    testObjectRetrieval(testCase,
            (int64_t(*)(Flower *flower)) flower_getGroupNumber, (void *(*)(
                    Flower *)) flower_getFirstGroup,
            (Name(*)(void *)) group_getName,
            (void *(*)(Flower *, Name name)) flower_getGroup, (void *(*)(
                    Flower *)) flower_getGroupIterator,
            (void(*)(void *)) flower_destructGroupIterator,
            (void *(*)(void *)) flower_getNextGroup,
            (void *(*)(void *)) flower_getPreviousGroup,
            (void *(*)(void *)) flower_copyGroupIterator, group,
            group2);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_getEndNumber(CuTest *testCase) {
    /*
     * Tests the different end number functions.
     */
    cactusFlowerTestSetup(testCase);
    CuAssertTrue(testCase, flower_getEndNumber(flower) == 0);
    CuAssertTrue(testCase, flower_getBlockEndNumber(flower) == 0);
    CuAssertTrue(testCase, flower_getStubEndNumber(flower) == 0);
    CuAssertTrue(testCase, flower_getFreeStubEndNumber(flower) == 0);
    CuAssertTrue(testCase, flower_getAttachedStubEndNumber(flower) == 0);
    int64_t blockNumber = 10;
    int64_t freeStubEndNumber = 5;
    int64_t attachedStubEndNumber = 3;
    int64_t i;
    for (i = 0; i < blockNumber; i++) {
        block_construct(1, flower);
    }
    for (i = 0; i < freeStubEndNumber; i++) {
        end_construct(0, flower);
    }
    for (i = 0; i < attachedStubEndNumber; i++) {
        end_construct(1, flower);
    }

    CuAssertTrue(testCase, flower_getEndNumber(flower) == blockNumber*2 + freeStubEndNumber + attachedStubEndNumber);
    CuAssertTrue(testCase, flower_getBlockEndNumber(flower) == blockNumber*2);
    CuAssertTrue(testCase, flower_getStubEndNumber(flower) == freeStubEndNumber + attachedStubEndNumber);
    CuAssertTrue(testCase, flower_getFreeStubEndNumber(flower) == freeStubEndNumber);
    CuAssertTrue(testCase, flower_getAttachedStubEndNumber(flower) == attachedStubEndNumber);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_builtBlocks(CuTest *testCase) {
    cactusFlowerTestSetup(testCase);

    CuAssertTrue(testCase, !flower_builtBlocks(flower));
    flower_setBuiltBlocks(flower, 0);
    CuAssertTrue(testCase, !flower_builtBlocks(flower));
    flower_setBuiltBlocks(flower, 1);
    CuAssertTrue(testCase, flower_builtBlocks(flower));
    flower_setBuiltBlocks(flower, 0);
    CuAssertTrue(testCase, !flower_builtBlocks(flower));

    cactusFlowerTestTeardown(testCase);
}

void testFlower_isLeaf(CuTest *testCase) {
    cactusFlowerTestSetup(testCase);
    CuAssertTrue(testCase, flower_isLeaf(flower));
    Group *group = group_construct2(flower);
    CuAssertTrue(testCase, flower_isLeaf(flower));
    group_makeNestedFlower(group);
    CuAssertTrue(testCase, !flower_isLeaf(flower));
    cactusFlowerTestTeardown(testCase);
}

void testFlower_isTerminal(CuTest *testCase) {
    cactusFlowerTestSetup(testCase);
    CuAssertTrue(testCase, flower_isTerminal(flower));
    group_construct2(flower);
    CuAssertTrue(testCase, flower_isTerminal(flower));
    end_construct(0, flower);
    CuAssertTrue(testCase, flower_isTerminal(flower));
    block_construct(1, flower);
    CuAssertTrue(testCase, !flower_isTerminal(flower));
    cactusFlowerTestTeardown(testCase);
}

CuSuite* cactusFlowerTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testFlower_getName);
    SUITE_ADD_TEST(suite, testFlower_getCactusDisk);
    SUITE_ADD_TEST(suite, testFlower_getEventTree);
    SUITE_ADD_TEST(suite, testFlower_sequence);
    SUITE_ADD_TEST(suite, testFlower_cap);
    SUITE_ADD_TEST(suite, testFlower_end);
    SUITE_ADD_TEST(suite, testFlower_getEndNumber);
    SUITE_ADD_TEST(suite, testFlower_group);
    SUITE_ADD_TEST(suite, testFlower_chain);
    SUITE_ADD_TEST(suite, testFlower_getTrivialChainNumber);
    SUITE_ADD_TEST(suite, testFlower_builtBlocks);
    SUITE_ADD_TEST(suite, testFlower_isLeaf);
    SUITE_ADD_TEST(suite, testFlower_isTerminal);
    SUITE_ADD_TEST(suite, testFlower_constructAndDestruct);
    return suite;
}
