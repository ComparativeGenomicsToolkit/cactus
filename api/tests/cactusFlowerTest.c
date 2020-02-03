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
static MetaSequence *metaSequence;
static MetaSequence *metaSequence2;
static Sequence *sequence;
static Sequence *sequence2;
static End *end;
static End *end2;
static Block *block;
static Block *block2;
static Group *group;
static Group *group2;
static Face *face;
static Face *face2;
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
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
        cactusDisk = NULL;
        flower = NULL;
        eventTree = NULL;
        metaSequence = NULL;
        sequence = NULL;
    }
}

static void cactusFlowerTestSetup(CuTest* testCase) {
    cactusFlowerTestTeardown(testCase);
    cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
    flower = flower_construct(cactusDisk);
    eventTree = eventTree_construct2(cactusDisk);
}

static void sequenceSetup() {
    metaSequence = metaSequence_construct(0, 10, "ACTGACTGAC", ">one",
            event_getName(eventTree_getRootEvent(eventTree)), cactusDisk);
    sequence = sequence_construct(metaSequence, flower);
    metaSequence2 = metaSequence_construct(0, 10, "ACTGACTGAC", ">two",
            event_getName(eventTree_getRootEvent(eventTree)), cactusDisk);
    sequence2 = sequence_construct(metaSequence2, flower);
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

static void segmentsSetup() {
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

static void facesSetup() {
    face = face_construct(flower);
    face2 = face_construct(flower);
    if (flower_getFirstFace(flower) == face2) { //switch round, because we have no guarantted order for faces..
        Face *face3 = face2;
        face2 = face;
        face = face3;
    }
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
    CuAssertTrue(testCase, getFirstObjectFn(flower) == object);

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

void testFlower_segment(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    segmentsSetup();
    testObjectRetrieval(testCase,
            (int64_t(*)(Flower *flower)) flower_getSegmentNumber, (void *(*)(
                    Flower *)) flower_getFirstSegment,
            (Name(*)(void *)) segment_getName,
            (void *(*)(Flower *, Name name)) flower_getSegment, (void *(*)(
                    Flower *)) flower_getSegmentIterator,
            (void(*)(void *)) flower_destructSegmentIterator,
            (void *(*)(void *)) flower_getNextSegment,
            (void *(*)(void *)) flower_getPreviousSegment,
            (void *(*)(void *)) flower_copySegmentIterator,
            segment, segment2);
    cactusFlowerTestTeardown(testCase);
}

void testFlower_block(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    blocksSetup();
    testObjectRetrieval(testCase,
            (int64_t(*)(Flower *flower)) flower_getBlockNumber, (void *(*)(
                    Flower *)) flower_getFirstBlock,
            (Name(*)(void *)) block_getName,
            (void *(*)(Flower *, Name name)) flower_getBlock, (void *(*)(
                    Flower *)) flower_getBlockIterator,
            (void(*)(void *)) flower_destructBlockIterator,
            (void *(*)(void *)) flower_getNextBlock,
            (void *(*)(void *)) flower_getPreviousBlock,
            (void *(*)(void *)) flower_copyBlockIterator, block,
            block2);
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

void testFlower_face(CuTest* testCase) {
    cactusFlowerTestSetup(testCase);
    facesSetup();
    testObjectRetrieval(testCase,
            (int64_t(*)(Flower *flower)) flower_getFaceNumber, (void *(*)(
                    Flower *)) flower_getFirstFace, NULL, NULL, (void *(*)(
                    Flower *)) flower_getFaceIterator,
            (void(*)(void *)) flower_destructFaceIterator,
            (void *(*)(void *)) flower_getNextFace,
            (void *(*)(void *)) flower_getPreviousFace,
            (void *(*)(void *)) flower_copyFaceIterator, face, face2);
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

void testFlower_builtTrees(CuTest *testCase) {
    cactusFlowerTestSetup(testCase);

    CuAssertTrue(testCase, !flower_builtTrees(flower));
    flower_setBuiltTrees(flower, 0);
    CuAssertTrue(testCase, !flower_builtTrees(flower));
    flower_setBuiltTrees(flower, 1);
    CuAssertTrue(testCase, flower_builtTrees(flower));
    flower_setBuiltTrees(flower, 0);
    CuAssertTrue(testCase, !flower_builtTrees(flower));

    cactusFlowerTestTeardown(testCase);
}

void testFlower_builtFaces(CuTest *testCase) {
    cactusFlowerTestSetup(testCase);

    CuAssertTrue(testCase, !flower_builtFaces(flower));
    flower_setBuildFaces(flower, 0);
    CuAssertTrue(testCase, !flower_builtFaces(flower));
    flower_setBuildFaces(flower, 1);
    CuAssertTrue(testCase, flower_builtFaces(flower));
    flower_setBuildFaces(flower, 0);
    CuAssertTrue(testCase, !flower_builtFaces(flower));

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

void testFlower_removeIfRedundant(CuTest *testCase) {
    /*
     * Do a simple test to see if function can remove a redundant flower.
     */
    cactusFlowerTestSetup(testCase);
    endsSetup();

    //First construct a redundant flower from the root.
    Flower *flower2 = flower_construct(cactusDisk);
    Group *group = group_construct(flower, flower2);
    end_setGroup(end, group);
    end_setGroup(end2, group);

    //Now hang another couple of flowers of that.
    Flower *flower3 = flower_construct(cactusDisk);
    group_construct(flower2, flower3);

    //Now hang another flower of that.
    Group *group3b = group_construct2(flower2);

    //Finally hang one more flower on the end..
    Flower *flower4 = flower_construct(cactusDisk);
    group_construct(flower3, flower4);

    //Copy the ends into the flowers.
    end_copyConstruct(end, flower2);
    end_copyConstruct(end2, flower2);
    end_copyConstruct(end, flower3);
    end_setGroup(flower_getEnd(flower2, end_getName(end2)), group3b);
    end_copyConstruct(end, flower4);

    //st_uglyf("I got %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 "\n", flower_getName(flower), flower_getName(flower2), flower_getName(flower3), flower_getName(flower4));

    //Write the mess to disk.
    cactusDisk_write(cactusDisk);

    //Now test the removal function (check we get a negative on this leaf).
    CuAssertTrue(testCase, !flower_removeIfRedundant(flower4));
    //Check we can't remove the root..
    CuAssertTrue(testCase, !flower_removeIfRedundant(flower));

    //We will remove flower2

    //Before
    CuAssertTrue(testCase, flower_getGroupNumber(flower) == 1);
    CuAssertTrue(testCase, group_getFlower(flower_getParentGroup(flower2)) == flower);

    CuAssertTrue(testCase, flower_removeIfRedundant(flower2));

    //After, check the flower/group connections
    CuAssertTrue(testCase, flower_getGroupNumber(flower) == 2);
    CuAssertTrue(testCase, !flower_isLeaf(flower));
    CuAssertTrue(testCase, group_getFlower(flower_getParentGroup(flower3)) == flower);
    group3b = end_getGroup(end2);
    CuAssertTrue(testCase, group_getFlower(group3b) == flower);
    CuAssertTrue(testCase, group_isLeaf(group3b));
    CuAssertTrue(testCase, flower_getGroup(flower, flower_getName(flower3)) == flower_getParentGroup(flower3));
    //Check the ends..
    CuAssertTrue(testCase, flower_getEndNumber(flower) == 2);
    CuAssertTrue(testCase, flower_getEndNumber(flower3) == 1);
    CuAssertTrue(testCase, group_getEndNumber(group3b) == 1);
    CuAssertTrue(testCase, end_getGroup(end) == flower_getParentGroup(flower3));
    CuAssertTrue(testCase, end_getGroup(end2) == group3b);
    CuAssertTrue(testCase, flower_getEnd(flower3, end_getName(end)) != NULL);
    //Check the child of 3 is still okay..
    CuAssertTrue(testCase, group_getFlower(flower_getParentGroup(flower4)) == flower3);

    //Now do removal of flower3
    CuAssertTrue(testCase, !flower_removeIfRedundant(flower));
    CuAssertTrue(testCase, !flower_removeIfRedundant(flower4));
    CuAssertTrue(testCase, flower_removeIfRedundant(flower3));
    //Check groups again
    CuAssertTrue(testCase, flower_getGroupNumber(flower) == 2);
    CuAssertTrue(testCase, !flower_isLeaf(flower));
    CuAssertTrue(testCase, group_getFlower(flower_getParentGroup(flower4)) == flower);
    CuAssertTrue(testCase, group_getFlower(group3b) == flower);
    CuAssertTrue(testCase, flower_getGroup(flower, flower_getName(flower4)) == flower_getParentGroup(flower4));
    //Check the ends again..
    CuAssertTrue(testCase, flower_getEndNumber(flower) == 2);
    CuAssertTrue(testCase, flower_getEndNumber(flower4) == 1);
    CuAssertTrue(testCase, group_getEndNumber(group3b) == 1);
    CuAssertTrue(testCase, end_getGroup(end) == flower_getParentGroup(flower4));
    CuAssertTrue(testCase, end_getGroup(end2) == group3b);
    CuAssertTrue(testCase, flower_getEnd(flower4, end_getName(end)) != NULL);

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
    SUITE_ADD_TEST(suite, testFlower_segment);
    SUITE_ADD_TEST(suite, testFlower_block);
    SUITE_ADD_TEST(suite, testFlower_group);
    SUITE_ADD_TEST(suite, testFlower_chain);
    SUITE_ADD_TEST(suite, testFlower_getTrivialChainNumber);
    SUITE_ADD_TEST(suite, testFlower_face);
    //SUITE_ADD_TEST(suite, testFlower_mergeFlowers);
    SUITE_ADD_TEST(suite, testFlower_builtBlocks);
    SUITE_ADD_TEST(suite, testFlower_builtTrees);
    SUITE_ADD_TEST(suite, testFlower_builtFaces);
    SUITE_ADD_TEST(suite, testFlower_isLeaf);
    SUITE_ADD_TEST(suite, testFlower_isTerminal);
    SUITE_ADD_TEST(suite, testFlower_removeIfRedundant);
    SUITE_ADD_TEST(suite, testFlower_constructAndDestruct);
    return suite;
}
