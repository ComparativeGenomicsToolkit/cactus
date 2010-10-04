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
static Reference *reference;

/*
 * Setup/teardown functions.
 */

static void cactusFlowerTestTeardown() {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(cactusDisk);
        cactusDisk = NULL;
        flower = NULL;
        eventTree = NULL;
        metaSequence = NULL;
        sequence = NULL;
        reference = NULL;
    }
}

static void cactusFlowerTestSetup() {
    cactusFlowerTestTeardown();
    cactusDisk = testCommon_getTemporaryCactusDisk();
    flower = flower_construct(cactusDisk);
    eventTree = eventTree_construct2(flower);
    assert(flower_getReference(flower) == NULL);
    reference = reference_construct(flower);
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
static void testObjectRetrieval(CuTest* testCase, void(*setupFn)(),
        int32_t(*getObjectNumberFn)(Flower *flower), void *(*getFirstObjectFn)(
                Flower *flower), Name(*objectGetNameFn)(void *),
        void *(*getObjectFn)(Flower *flower, Name name),
        void *(*constructIterator)(Flower *flower), void(*destructIterator)(
                void *iterator), void *(*getNext)(void *iterator),
        void *(*getPrevious)(void *iterator), void *(*copyIterator)(
                void *iterator), void *object, void *object2) {
    cactusFlowerTestSetup();
    /*
     * Test number function
     */
    CuAssertTrue(testCase, getObjectNumberFn(flower) == 0);
    setupFn();
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
    cactusFlowerTestTeardown();
}

/*
 * Now all the actual tests.
 */

void testFlower_constructAndDestruct(CuTest* testCase) {
    cactusFlowerTestSetup();
    CuAssertTrue(testCase, flower != NULL);
    cactusFlowerTestTeardown();
}

void testFlower_getName(CuTest* testCase) {
    cactusFlowerTestSetup();
    CuAssertTrue(testCase, flower_getName(flower) != NULL_NAME);
    CuAssertTrue(testCase, cactusDisk_getFlower(cactusDisk, flower_getName(flower)) == flower);
    cactusFlowerTestTeardown();
}

void testFlower_getCactusDisk(CuTest* testCase) {
    cactusFlowerTestSetup();
    CuAssertTrue(testCase, flower_getCactusDisk(flower) == cactusDisk);
    cactusFlowerTestTeardown();
}

void testFlower_getEventTree(CuTest* testCase) {
    cactusFlowerTestSetup();
    CuAssertTrue(testCase, flower_getEventTree(flower) == eventTree);
    cactusFlowerTestTeardown();
}

void testFlower_getReference(CuTest* testCase) {
    cactusFlowerTestSetup();
    CuAssertTrue(testCase, flower_getReference(flower) == reference);
    cactusFlowerTestTeardown();
}

void testFlower_sequence(CuTest* testCase) {
    testObjectRetrieval(testCase, sequenceSetup,
            (int32_t(*)(Flower *flower)) flower_getSequenceNumber, (void *(*)(
                    Flower *)) flower_getFirstSequence,
            (Name(*)(void *)) sequence_getName,
            (void *(*)(Flower *, Name name)) flower_getSequence, (void *(*)(
                    Flower *)) flower_getSequenceIterator,
            (void(*)(void *)) flower_destructSequenceIterator, (void *(*)(
                    void *)) flower_getNextSequence,
            (void *(*)(void *)) flower_getPreviousSequence,
            (void *(*)(void *)) flower_copySequenceIterator,
            (void *)sequence, (void *)sequence2);
}

void testFlower_cap(CuTest* testCase) {
    testObjectRetrieval(testCase, capsSetup,
            (int32_t(*)(Flower *flower)) flower_getCapNumber, (void *(*)(
                    Flower *)) flower_getFirstCap,
            (Name(*)(void *)) cap_getName,
            (void *(*)(Flower *, Name name)) flower_getCap,
            (void *(*)(Flower *)) flower_getCapIterator,
            (void(*)(void *)) flower_destructCapIterator,
            (void *(*)(void *)) flower_getNextCap,
            (void *(*)(void *)) flower_getPreviousCap,
            (void *(*)(void *)) flower_copyCapIterator, (void *)cap,
            (void *)cap2);
}

void testFlower_end(CuTest* testCase) {
    testObjectRetrieval(testCase, endsSetup,
            (int32_t(*)(Flower *flower)) flower_getEndNumber, (void *(*)(
                    Flower *)) flower_getFirstEnd,
            (Name(*)(void *)) end_getName,
            (void *(*)(Flower *, Name name)) flower_getEnd,
            (void *(*)(Flower *)) flower_getEndIterator,
            (void(*)(void *)) flower_destructEndIterator,
            (void *(*)(void *)) flower_getNextEnd,
            (void *(*)(void *)) flower_getPreviousEnd,
            (void *(*)(void *)) flower_copyEndIterator, (void *)end,
            (void *)end2);
}

void testFlower_segment(CuTest* testCase) {
    testObjectRetrieval(testCase, segmentsSetup,
            (int32_t(*)(Flower *flower)) flower_getSegmentNumber, (void *(*)(
                    Flower *)) flower_getFirstSegment,
            (Name(*)(void *)) segment_getName,
            (void *(*)(Flower *, Name name)) flower_getSegment, (void *(*)(
                    Flower *)) flower_getSegmentIterator,
            (void(*)(void *)) flower_destructSegmentIterator,
            (void *(*)(void *)) flower_getNextSegment,
            (void *(*)(void *)) flower_getPreviousSegment,
            (void *(*)(void *)) flower_copySegmentIterator,
            (void *)segment, (void *)segment2);
}

void testFlower_block(CuTest* testCase) {
    testObjectRetrieval(testCase, blocksSetup,
            (int32_t(*)(Flower *flower)) flower_getBlockNumber, (void *(*)(
                    Flower *)) flower_getFirstBlock,
            (Name(*)(void *)) block_getName,
            (void *(*)(Flower *, Name name)) flower_getBlock, (void *(*)(
                    Flower *)) flower_getBlockIterator,
            (void(*)(void *)) flower_destructBlockIterator,
            (void *(*)(void *)) flower_getNextBlock,
            (void *(*)(void *)) flower_getPreviousBlock,
            (void *(*)(void *)) flower_copyBlockIterator, (void *)block,
            (void *)block2);
}

void testFlower_chain(CuTest* testCase) {
    testObjectRetrieval(testCase, chainsSetup,
            (int32_t(*)(Flower *flower)) flower_getChainNumber, (void *(*)(
                    Flower *)) flower_getFirstChain,
            (Name(*)(void *)) chain_getName,
            (void *(*)(Flower *, Name name)) flower_getChain, (void *(*)(
                    Flower *)) flower_getChainIterator,
            (void(*)(void *)) flower_destructChainIterator,
            (void *(*)(void *)) flower_getNextChain,
            (void *(*)(void *)) flower_getPreviousChain,
            (void *(*)(void *)) flower_copyChainIterator, (void *)chain,
            (void *)chain2);
}

void testFlower_getTrivialChainNumber(CuTest* testCase) {
    cactusFlowerTestSetup();
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
    cactusFlowerTestTeardown();
}

void testFlower_group(CuTest* testCase) {
    testObjectRetrieval(testCase, groupsSetup,
            (int32_t(*)(Flower *flower)) flower_getGroupNumber, (void *(*)(
                    Flower *)) flower_getFirstGroup,
            (Name(*)(void *)) group_getName,
            (void *(*)(Flower *, Name name)) flower_getGroup, (void *(*)(
                    Flower *)) flower_getGroupIterator,
            (void(*)(void *)) flower_destructGroupIterator,
            (void *(*)(void *)) flower_getNextGroup,
            (void *(*)(void *)) flower_getPreviousGroup,
            (void *(*)(void *)) flower_copyGroupIterator, (void *)group,
            (void *)group2);
}

void testFlower_face(CuTest* testCase) {
    testObjectRetrieval(testCase, facesSetup,
            (int32_t(*)(Flower *flower)) flower_getFaceNumber, (void *(*)(
                    Flower *)) flower_getFirstFace, NULL, NULL, (void *(*)(
                    Flower *)) flower_getFaceIterator,
            (void(*)(void *)) flower_destructFaceIterator,
            (void *(*)(void *)) flower_getNextFace,
            (void *(*)(void *)) flower_getPreviousFace,
            (void *(*)(void *)) flower_copyFaceIterator, (void *)face,
            (void *)face2);
}

void testFlower_getEndNumber(CuTest *testCase) {
    /*
     * Tests the different end number functions.
     */
    cactusFlowerTestSetup();
    CuAssertTrue(testCase, flower_getEndNumber(flower) == 0);
    CuAssertTrue(testCase, flower_getBlockEndNumber(flower) == 0);
    CuAssertTrue(testCase, flower_getStubEndNumber(flower) == 0);
    CuAssertTrue(testCase, flower_getFreeStubEndNumber(flower) == 0);
    CuAssertTrue(testCase, flower_getAttachedStubEndNumber(flower) == 0);
    int32_t blockNumber = 10;
    int32_t freeStubEndNumber = 5;
    int32_t attachedStubEndNumber = 3;
    int32_t i;
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
    cactusFlowerTestTeardown();
}

/*void testFlower_mergeFlowers(CuTest *testCase) {
 cactusFlowerTestSetup();
 //construct the flowers to merge...
 Flower *flower1 = flower_construct(cactusDisk);
 Flower *flower2 = flower_construct(cactusDisk);

 //construct flowers that are children of the flowers to merge.
 Flower *flower3 = flower_construct(cactusDisk);
 Flower *flower4 = flower_construct(cactusDisk);
 Flower *flower5 = flower_construct(cactusDisk);

 //Make event tree
 MetaEvent *internalMetaEvent = metaEvent_construct("INTERNAL", cactusDisk);
 MetaEvent *leafMetaEvent1 = metaEvent_construct("LEAF1", cactusDisk);
 MetaEvent *leafMetaEvent2 = metaEvent_construct("LEAF2", cactusDisk);
 Event *internalEvent = event_construct(internalMetaEvent, 0.5, eventTree_getRootEvent(eventTree), eventTree);
 Event *leafEvent1 = event_construct(leafMetaEvent1, 0.2, internalEvent, eventTree);
 event_construct(leafMetaEvent2, 1.3, internalEvent, eventTree);

 //Copy the event tree into the children (this is tested by the merge event function)..
 EventTree *eventTree1 = eventTree_copyConstruct(eventTree, flower1, NULL);
 eventTree_copyConstruct(eventTree, flower2, NULL);
 MetaEvent *unaryInternalMetaEvent1 = metaEvent_construct("UNARY1", cactusDisk);
 MetaEvent *unaryInternalMetaEvent2 = metaEvent_construct("UNARY2", cactusDisk);
 MetaEvent *unaryInternalMetaEvent3 = metaEvent_construct("UNARY3", cactusDisk);

 Event *internalEventChild = eventTree_getEvent(eventTree1, event_getName(internalEvent));
 Event *unaryEvent1 = event_construct2(unaryInternalMetaEvent1, 0.1,
 internalEventChild, eventTree_getEvent(eventTree1, event_getName(leafEvent1)), eventTree1);
 Event *unaryEvent2 = event_construct2(unaryInternalMetaEvent2, 0.1,
 internalEventChild, unaryEvent1, eventTree1);
 event_construct2(unaryInternalMetaEvent3, 0.1,
 internalEventChild, unaryEvent2, eventTree1);
 CuAssertTrue(testCase, eventTree_getEventNumber(eventTree1) == 7);

 //Make some sequences
 MetaSequence *metaSequence1 = metaSequence_construct(0, 5, "ACTGG", "one", event_getName(unaryEvent1), cactusDisk);
 MetaSequence *metaSequence2 = metaSequence_construct(0, 5, "CCCCC", "two", event_getName(unaryEvent2), cactusDisk);
 MetaSequence *metaSequence3 = metaSequence_construct(0, 5, "TTTTT", "three", event_getName(leafEvent1), cactusDisk);
 Sequence *sequence1 = sequence_construct(metaSequence1, flower1);
 Sequence *sequence2 = sequence_construct(metaSequence2, flower1);
 Sequence *sequence3 = sequence_construct(metaSequence3, flower2);

 //Make children and parent flowers and some groups..
 Group *group3 = group_construct(flower1, flower3);
 Group *group4 = group_construct(flower2, flower4);
 Group *group5 = group_construct(flower2, flower5);
 Group *group6 = group_construct2(flower1);
 Group *group7 = group_construct2(flower2);

 //Make some stubs to put in the ends..
 End *end1 = end_construct(0, flower1);
 End *end2 = end_construct(0, flower1);
 End *end3 = end_construct(0, flower2);
 end_setGroup(end1, group6);
 end_setGroup(end2, group6);
 end_setGroup(end3, group7);

 //Make a few caps
 Cap *cap1 = cap_construct2(end1, 0, 1, 1, sequence1);
 Cap *cap2 = cap_construct2(end1, 0, 1, 1, sequence2);
 cap_makeParentAndChild(cap2, cap1);
 Cap *cap3 = cap_construct2(end3, 0, 1, 1, sequence3);

 //Make some blocks..
 Block *block1 = block_construct(0, flower1);
 Block *block2 = block_construct(0, flower2);
 Block *block3 = block_construct(0, flower2);
 end_setGroup(block_get5End(block3), group7);

 //Make a segment
 Segment *segment1 = segment_construct(block1, eventTree_getEvent(eventTree1, metaEvent_getName(leafMetaEvent1)));

 //Make some chains...
 Chain *chain1 = chain_construct(flower1);
 Chain *chain2 = chain_construct(flower2);

 //Make a couple of links in the chain
 Link *link1 = link_construct(end1, end2, group6, chain1);
 Link *link2 = link_construct(end3, block_get5End(block3), group7, chain2);

 Name flowerName1 = flower_getName(flower1);
 Flower *flower6 = flower_mergeFlowers(flower1, flower2);

 //Check the events
 EventTree *eventTree3 = flower_getEventTree(flower6);
 CuAssertTrue(testCase, eventTree_getEventNumber(eventTree3) == 7);
 CuAssertTrue(testCase, eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent1)) != NULL);
 CuAssertTrue(testCase, eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent2)) != NULL);
 CuAssertTrue(testCase, eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent3)) != NULL);

 CuAssertTrue(testCase, event_getName(event_getParent(eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent1)))) == metaEvent_getName(unaryInternalMetaEvent2));
 CuAssertTrue(testCase, event_getName(event_getParent(eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent2)))) == metaEvent_getName(unaryInternalMetaEvent3));
 CuAssertTrue(testCase, event_getName(event_getParent(eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent3)))) == metaEvent_getName(internalMetaEvent));
 CuAssertTrue(testCase, event_getName(event_getParent(eventTree_getEvent(eventTree3, metaEvent_getName(leafMetaEvent1)))) == metaEvent_getName(unaryInternalMetaEvent1));

 //Check the sequences
 CuAssertTrue(testCase, flower_getSequence(flower6, metaSequence_getName(metaSequence2)) != NULL);
 CuAssertTrue(testCase, flower_getSequence(flower6, metaSequence_getName(metaSequence3)) != NULL);
 CuAssertTrue(testCase, flower_getSequence(flower6, metaSequence_getName(metaSequence1)) != NULL);

 CuAssertTrue(testCase, sequence_getEvent(flower_getSequence(flower6, metaSequence_getName(metaSequence1))) == eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent1)));
 CuAssertTrue(testCase, sequence_getEvent(flower_getSequence(flower6, metaSequence_getName(metaSequence2))) == eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent2)));
 CuAssertTrue(testCase, sequence_getEvent(flower_getSequence(flower6, metaSequence_getName(metaSequence3))) == eventTree_getEvent(eventTree3, metaEvent_getName(leafMetaEvent1)));

 //Check the groups..
 CuAssertTrue(testCase, flower_getGroupNumber(flower6) == 5);
 CuAssertTrue(testCase, flower_getGroup(flower6, group_getName(group3)) != NULL);
 CuAssertTrue(testCase, flower_getGroup(flower6, group_getName(group4)) != NULL);
 CuAssertTrue(testCase, flower_getGroup(flower6, group_getName(group5)) != NULL);
 CuAssertTrue(testCase, flower_getGroup(flower6, group_getName(group6)) != NULL);
 CuAssertTrue(testCase, flower_getGroup(flower6, group_getName(group7)) != NULL);
 CuAssertTrue(testCase, flower_getGroup(flower6, group_getName(group3)) == flower_getParentGroup(flower3));
 CuAssertTrue(testCase, flower_getGroup(flower6, group_getName(group4)) == flower_getParentGroup(flower4));
 CuAssertTrue(testCase, flower_getGroup(flower6, group_getName(group5)) == flower_getParentGroup(flower5));

 //Check the ends
 CuAssertTrue(testCase, flower_getEndNumber(flower6) == 9);
 CuAssertTrue(testCase, flower_getEnd(flower6, end_getName(end1)) == end1);
 CuAssertTrue(testCase, flower_getEnd(flower6, end_getName(end2)) == end2);
 CuAssertTrue(testCase, flower_getEnd(flower6, end_getName(end3)) == end3);
 CuAssertTrue(testCase, flower_getEnd(flower6, end_getName(block_get5End(block1))) == block_get5End(block1));
 CuAssertTrue(testCase, flower_getEnd(flower6, end_getName(block_get5End(block2))) == block_get5End(block2));
 CuAssertTrue(testCase, flower_getEnd(flower6, end_getName(block_get5End(block3))) == block_get5End(block3));
 CuAssertTrue(testCase, flower_getEnd(flower6, end_getName(block_get3End(block1))) == block_get3End(block1));
 CuAssertTrue(testCase, flower_getEnd(flower6, end_getName(block_get3End(block2))) == block_get3End(block2));
 CuAssertTrue(testCase, flower_getEnd(flower6, end_getName(block_get3End(block3))) == block_get3End(block3));

 //Check we have the right index of caps.
 CuAssertTrue(testCase, flower_getCapNumber(flower6) == 5);
 CuAssertTrue(testCase, flower_getCap(flower6, cap_getName(cap1)) == cap1);
 CuAssertTrue(testCase, flower_getCap(flower6, cap_getName(cap2)) == cap2);
 CuAssertTrue(testCase, flower_getCap(flower6, cap_getName(cap3)) == cap3);

 //Check the the events have been reassigned for the caps
 CuAssertTrue(testCase, cap_getEvent(cap1) == eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent1)));
 CuAssertTrue(testCase, cap_getEvent(cap2) == eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent2)));
 CuAssertTrue(testCase, cap_getEvent(cap3) == eventTree_getEvent(eventTree3, metaEvent_getName(leafMetaEvent1)));

 //Check the caps have the right sequences.
 CuAssertTrue(testCase, cap_getSequence(cap1) == flower_getSequence(flower6, metaSequence_getName(metaSequence1)));
 CuAssertTrue(testCase, cap_getSequence(cap2) == flower_getSequence(flower6, metaSequence_getName(metaSequence2)));
 CuAssertTrue(testCase, cap_getSequence(cap3) == flower_getSequence(flower6, metaSequence_getName(metaSequence3)));

 //Check the blocks
 CuAssertTrue(testCase, flower_getBlockNumber(flower6) == 3);
 CuAssertTrue(testCase, flower_getBlock(flower6, block_getName(block1)) == block1);
 CuAssertTrue(testCase, flower_getBlock(flower6, block_getName(block2)) == block2);
 CuAssertTrue(testCase, flower_getBlock(flower6, block_getName(block3)) == block3);

 //Check the segments
 CuAssertTrue(testCase, flower_getSegmentNumber(flower6) == 1);
 CuAssertTrue(testCase, flower_getSegment(flower6, segment_getName(segment1)) == segment1);
 CuAssertTrue(testCase, segment_getSequence(segment1) == NULL);
 CuAssertTrue(testCase, segment_getEvent(segment1) == eventTree_getEvent(eventTree3, event_getName(leafEvent1)));

 //Check the chains
 CuAssertTrue(testCase, flower_getChainNumber(flower6) == 2);
 CuAssertTrue(testCase, flower_getChain(flower6, chain_getName(chain1)) == chain1);
 CuAssertTrue(testCase, flower_getChain(flower6, chain_getName(chain2)) == chain2);
 CuAssertTrue(testCase, link_getChain(link1) == chain1);
 CuAssertTrue(testCase, link_getChain(link2) == chain2);
 CuAssertTrue(testCase, link_get5End(link1) == end1);
 CuAssertTrue(testCase, link_get3End(link1) == end2);
 CuAssertTrue(testCase, link_get5End(link2) == end3);
 CuAssertTrue(testCase, link_get3End(link2) == block_get5End(block3));

 //Check flower1 is no longer in the database anywhere...
 CuAssertTrue(testCase, cactusDisk_getFlower(cactusDisk, flowerName1) == NULL);

 //Check the merged flower is in the database.
 CuAssertTrue(testCase, cactusDisk_getFlower(cactusDisk, flower_getName(flower6)) == flower6);

 cactusFlowerTestTeardown();
 }*/

void testFlower_builtBlocks(CuTest *testCase) {
    cactusFlowerTestSetup();

    CuAssertTrue(testCase, !flower_builtBlocks(flower));
    flower_setBuiltBlocks(flower, 0);
    CuAssertTrue(testCase, !flower_builtBlocks(flower));
    flower_setBuiltBlocks(flower, 1);
    CuAssertTrue(testCase, flower_builtBlocks(flower));
    flower_setBuiltBlocks(flower, 0);
    CuAssertTrue(testCase, !flower_builtBlocks(flower));

    cactusFlowerTestTeardown();
}

void testFlower_builtTrees(CuTest *testCase) {
    cactusFlowerTestSetup();

    CuAssertTrue(testCase, !flower_builtTrees(flower));
    flower_setBuiltTrees(flower, 0);
    CuAssertTrue(testCase, !flower_builtTrees(flower));
    flower_setBuiltTrees(flower, 1);
    CuAssertTrue(testCase, flower_builtTrees(flower));
    flower_setBuiltTrees(flower, 0);
    CuAssertTrue(testCase, !flower_builtTrees(flower));

    cactusFlowerTestTeardown();
}

void testFlower_builtFaces(CuTest *testCase) {
    cactusFlowerTestSetup();

    CuAssertTrue(testCase, !flower_builtFaces(flower));
    flower_setBuildFaces(flower, 0);
    CuAssertTrue(testCase, !flower_builtFaces(flower));
    flower_setBuildFaces(flower, 1);
    CuAssertTrue(testCase, flower_builtFaces(flower));
    flower_setBuildFaces(flower, 0);
    CuAssertTrue(testCase, !flower_builtFaces(flower));

    cactusFlowerTestTeardown();
}

void testFlower_isLeaf(CuTest *testCase) {
    cactusFlowerTestSetup();
    CuAssertTrue(testCase, flower_isLeaf(flower));
    Group *group = group_construct2(flower);
    CuAssertTrue(testCase, flower_isLeaf(flower));
    group_makeNestedFlower(group);
    CuAssertTrue(testCase, !flower_isLeaf(flower));
    cactusFlowerTestTeardown();
}

void testFlower_isTerminal(CuTest *testCase) {
    cactusFlowerTestSetup();
    CuAssertTrue(testCase, flower_isTerminal(flower));
    group_construct2(flower);
    CuAssertTrue(testCase, flower_isTerminal(flower));
    end_construct(0, flower);
    CuAssertTrue(testCase, flower_isTerminal(flower));
    block_construct(1, flower);
    CuAssertTrue(testCase, !flower_isTerminal(flower));
    cactusFlowerTestTeardown();
}

void testFlower_removeIfRedundant(CuTest *testCase) {
    /*
     * Do a simple test to see if function can remove a redundant flower.
     */
    cactusFlowerTestSetup();
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

    //st_uglyf("I got %i %i %i %i\n", flower_getName(flower), flower_getName(flower2), flower_getName(flower3), flower_getName(flower4));

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

    cactusFlowerTestTeardown();
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
    SUITE_ADD_TEST(suite, testFlower_getReference);
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
