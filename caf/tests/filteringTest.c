#include "CuTest.h"
#include "sonLib.h"
#include "stCaf.h"
#include "stPinchGraphs.h"

static CactusDisk *cactusDisk;
static Flower *flower;
static EventTree *eventTree;
static Event *rootEvent;
static Event *ingroup1;
static Event *ingroup2;
static Event *ancestor;
static Event *outgroup1;
static Event *outgroup2;

// TODO: refactor and merge with the version in
// recoverableChainsTest.c

// Adds a thread with random nucleotides to the flower, and return its corresponding name in the pinch graph.
static Name addThreadToFlower(Flower *flower, Event *event, int64_t length) {
    char *dna = stRandom_getRandomDNAString(length, true, true, true);
    MetaSequence *metaSequence = metaSequence_construct(2, length, dna, "", event_getName(event), flower_getCactusDisk(flower));
    Sequence *sequence = sequence_construct(metaSequence, flower);

    End *end1 = end_construct2(0, 0, flower);
    End *end2 = end_construct2(1, 0, flower);
    Cap *cap1 = cap_construct2(end1, 1, 1, sequence);
    Cap *cap2 = cap_construct2(end2, length + 2, 1, sequence);
    cap_makeAdjacent(cap1, cap2);

    free(dna);
    return cap_getName(cap1);
}

static stCactusEdgeEnd *getChainEndFromBlock(stCactusGraph *cactusGraph,
                                             stPinchBlock *block) {
    stCactusEdgeEnd *result = NULL;
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(cactusGraph);
    stCactusNode *node;
    while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stCactusNodeEdgeEndIt endIt = stCactusNode_getEdgeEndIt(node);
        stCactusEdgeEnd *edgeEnd;
        while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&endIt)) != NULL) {
            stPinchEnd *pinchEnd = stCactusEdgeEnd_getObject(edgeEnd);
            if (stPinchEnd_getBlock(pinchEnd) == block) {
                result = edgeEnd;
                break;
            }
        }
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
    return result;
}

static void teardown() {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(cactusDisk);
        cactusDisk = NULL;
    }
}

// Set up a basic flower with a tree looking like:
// ((ingroup1, ingroup2)ancestor, outgroup1, outgroup2)root;
static void setup() {
    teardown();
    cactusDisk = testCommon_getTemporaryCactusDisk();
    flower = flower_construct(cactusDisk);
    // A group must be constructed because stCaf_setup expects a leaf group.
    group_construct2(flower);
    flower_check(flower);

    eventTree = eventTree_construct2(flower);
    rootEvent = eventTree_getRootEvent(eventTree);
    ancestor = event_construct3("ancestor", 0.2, rootEvent, eventTree);
    outgroup1 = event_construct3("outgroup1", 0.2, rootEvent, eventTree);
    event_setOutgroupStatus(outgroup1, true);
    outgroup2 = event_construct3("outgroup2", 0.2, rootEvent, eventTree);
    event_setOutgroupStatus(outgroup2, true);
    ingroup1 = event_construct3("ingroup1", 0.2, ancestor, eventTree);
    ingroup2 = event_construct3("ingroup2", 0.2, ancestor, eventTree);
}

static void testChainHasUnequalNumberOfIngroupCopies(CuTest *testCase) {
    setup();
    Name ingroup1Seq1 = addThreadToFlower(flower, ingroup1, 100);
    Name ingroup1Seq2 = addThreadToFlower(flower, ingroup1, 100);
    Name ingroup2Seq1 = addThreadToFlower(flower, ingroup2, 100);
    Name ingroup2Seq2 = addThreadToFlower(flower, ingroup2, 100);
    Name outgroup1Seq1 = addThreadToFlower(flower, outgroup1, 100);

    stPinchThreadSet *threadSet = stCaf_setup(flower);

    stPinchThread *ingroup1Thread1 = stPinchThreadSet_getThread(threadSet, ingroup1Seq1);
    stPinchThread *ingroup1Thread2 = stPinchThreadSet_getThread(threadSet, ingroup1Seq2);
    stPinchThread *ingroup2Thread1 = stPinchThreadSet_getThread(threadSet, ingroup2Seq1);
    stPinchThread *ingroup2Thread2 = stPinchThreadSet_getThread(threadSet, ingroup2Seq2);
    stPinchThread *outgroup1Thread1 = stPinchThreadSet_getThread(threadSet, outgroup1Seq1);

    // Block A: A block with 2 segments from each ingroup and no
    // outgroups--this should fail
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread1, 10, 10, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup1Thread2, 10, 10, 10, true);
    stPinchThread_pinch(ingroup1Thread2, ingroup2Thread2, 10, 10, 10, true);

    // Block B: A block with 2 segments from each ingroup (but ingroup
    // 1's two segments are on the same thread) and an outgroup
    // segment--this should also fail
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread1, 30, 30, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup1Thread1, 20, 30, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread2, 30, 30, 10, true);
    stPinchThread_pinch(ingroup1Thread1, outgroup1Thread1, 30, 30, 10, true);

    // Block C: A block with 2 segments from ingroup 1 and 1 segment
    // from ingroup 2--this should pass
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread1, 40, 40, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup1Thread2, 40, 40, 10, true);

    // Now create the cactus graph and test that the filter classifies
    // these blocks (which are guaranteed to end up in separate
    // chains) correctly.
    stCactusNode *startCactusNode;
    stList *deadEndComponent;
    stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, 0,
                                                                  0.0, true, 0);

    // Test block A, which should fail.
    stPinchBlock *blockA = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 10));
    stCactusEdgeEnd *chainA = getChainEndFromBlock(cactusGraph, blockA);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopies(chainA, flower) == false);

    // Test block B, which should fail.
    stPinchBlock *blockB = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 30));
    stCactusEdgeEnd *chainB = getChainEndFromBlock(cactusGraph, blockB);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopies(chainB, flower) == false);

    // Test block C, which should pass.
    stPinchBlock *blockC = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 40));
    stCactusEdgeEnd *chainC = getChainEndFromBlock(cactusGraph, blockC);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopies(chainC, flower) == true);

    teardown();
}

CuSuite* filteringTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testChainHasUnequalNumberOfIngroupCopies);
//    SUITE_ADD_TEST(suite, testChainHasUnequalNumberOfIngroupCopiesOrNoOutgroup);
    return suite;
}
