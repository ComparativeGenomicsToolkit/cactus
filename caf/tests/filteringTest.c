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
    Sequence *sequence = sequence_construct(2, length, dna, "", event, flower_getCactusDisk(flower));
    flower_addSequence(flower, sequence);

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

static void teardown(CuTest* testCase) {
    if (cactusDisk != NULL) {
        cactusDisk_destruct(cactusDisk);
        cactusDisk = NULL;
    }
}

// Set up a basic flower with a tree looking like:
// ((ingroup1, ingroup2)ancestor, outgroup1, outgroup2)root;
// If the includeOutgroups parameter is false outgroup1 and outgroup2
// aren't added.
static void setup(CuTest* testCase, bool includeOutgroups) {
    teardown(testCase);
    cactusDisk = cactusDisk_construct();
    eventTree_construct2(cactusDisk);
    flower = flower_construct(cactusDisk);
    // A group must be constructed because stCaf_setup expects a leaf group.
    group_construct2(flower);
    flower_check(flower);

    eventTree = eventTree_construct2(cactusDisk);
    rootEvent = eventTree_getRootEvent(eventTree);
    ancestor = event_construct3("ancestor", 0.2, rootEvent, eventTree);
    if (includeOutgroups) {
        outgroup1 = event_construct3("outgroup1", 0.2, rootEvent, eventTree);
        event_setOutgroupStatus(outgroup1, true);
        outgroup2 = event_construct3("outgroup2", 0.2, rootEvent, eventTree);
        event_setOutgroupStatus(outgroup2, true);
    }
    ingroup1 = event_construct3("ingroup1", 0.2, ancestor, eventTree);
    ingroup2 = event_construct3("ingroup2", 0.2, ancestor, eventTree);
}

static void testChainHasUnequalNumberOfIngroupCopies(CuTest *testCase) {
    setup(testCase, true);
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

    // Block D: A block with one segment from ingroup 1 and 1 from the
    // outgroup--this should pass
    stPinchThread_pinch(ingroup1Thread1, outgroup1Thread1, 50, 50, 10, true);

    // Now create the cactus graph and test that the filter classifies
    // these blocks (which are guaranteed to end up in separate
    // chains) correctly.
    stCactusNode *startCactusNode;
    stPinchComponent *deadEndComponent;
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

    // Test block D, which should pass.
    stPinchBlock *blockD = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 50));
    stCactusEdgeEnd *chainD = getChainEndFromBlock(cactusGraph, blockD);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopies(chainD, flower) == true);

    teardown(testCase);
}

static void testChainHasUnequalNumberOfIngroupCopiesOrNoOutgroup(CuTest *testCase) {
    setup(testCase, true);

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
    // outgroups--this should pass
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread1, 10, 10, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup1Thread2, 10, 10, 10, true);
    stPinchThread_pinch(ingroup1Thread2, ingroup2Thread2, 10, 10, 10, true);

    // Block B: A block with 2 segments from each ingroup (but ingroup
    // 1's two segments are on the same thread) and an outgroup
    // segment--this should fail
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread1, 30, 30, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup1Thread1, 20, 30, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread2, 30, 30, 10, true);
    stPinchThread_pinch(ingroup1Thread1, outgroup1Thread1, 30, 30, 10, true);

    // Block C: A block with 2 segments from ingroup 1 and 1 segment
    // from ingroup 2--this should pass
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread1, 40, 40, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup1Thread2, 40, 40, 10, true);

    // Block D: A block with one segment from ingroup 1 and 1 from the
    // outgroup--this should pass
    stPinchThread_pinch(ingroup1Thread1, outgroup1Thread1, 50, 50, 10, true);

    // Block E: Same as block C but with an outgroup -- this should pass
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread1, 60, 60, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup1Thread2, 60, 60, 10, true);
    stPinchThread_pinch(ingroup1Thread1, outgroup1Thread1, 60, 60, 10, true);

    // Now create the cactus graph and test that the filter classifies
    // these blocks (which are guaranteed to end up in separate
    // chains) correctly.
    stCactusNode *startCactusNode;
    stPinchComponent *deadEndComponent;
    stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, 0,
                                                                  0.0, true, 0);

    // Test block A, which should pass.
    stPinchBlock *blockA = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 10));
    stCactusEdgeEnd *chainA = getChainEndFromBlock(cactusGraph, blockA);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup(chainA, flower) == true);

    // Test block B, which should fail.
    stPinchBlock *blockB = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 30));
    stCactusEdgeEnd *chainB = getChainEndFromBlock(cactusGraph, blockB);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup(chainB, flower) == false);

    // Test block C, which should pass.
    stPinchBlock *blockC = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 40));
    stCactusEdgeEnd *chainC = getChainEndFromBlock(cactusGraph, blockC);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup(chainC, flower) == true);

    // Test block D, which should also pass.
    stPinchBlock *blockD = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 50));
    stCactusEdgeEnd *chainD = getChainEndFromBlock(cactusGraph, blockD);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup(chainD, flower) == true);

    // Test block E, which should also pass.
    stPinchBlock *blockE = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 60));
    stCactusEdgeEnd *chainE = getChainEndFromBlock(cactusGraph, blockE);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup(chainE, flower) == true);

    teardown(testCase);
}

// Test chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup under the
// case where there are no outgroups in the tree.
static void testChainHasUnequalNumberOfIngroupCopiesOrNoOutgroup_noOutgroups(CuTest *testCase) {
    setup(testCase, false);

    Name ingroup1Seq1 = addThreadToFlower(flower, ingroup1, 100);
    Name ingroup1Seq2 = addThreadToFlower(flower, ingroup1, 100);
    Name ingroup2Seq1 = addThreadToFlower(flower, ingroup2, 100);
    Name ingroup2Seq2 = addThreadToFlower(flower, ingroup2, 100);

    stPinchThreadSet *threadSet = stCaf_setup(flower);

    stPinchThread *ingroup1Thread1 = stPinchThreadSet_getThread(threadSet, ingroup1Seq1);
    stPinchThread *ingroup1Thread2 = stPinchThreadSet_getThread(threadSet, ingroup1Seq2);
    stPinchThread *ingroup2Thread1 = stPinchThreadSet_getThread(threadSet, ingroup2Seq1);
    stPinchThread *ingroup2Thread2 = stPinchThreadSet_getThread(threadSet, ingroup2Seq2);

    // Block A: A block with 2 segments from each ingroup and no
    // outgroups--this should fail
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread1, 10, 10, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup1Thread2, 10, 10, 10, true);
    stPinchThread_pinch(ingroup1Thread2, ingroup2Thread2, 10, 10, 10, true);

    // Block C: A block with 2 segments from ingroup 1 and 1 segment
    // from ingroup 2--this should pass
    stPinchThread_pinch(ingroup1Thread1, ingroup2Thread1, 40, 40, 10, true);
    stPinchThread_pinch(ingroup1Thread1, ingroup1Thread2, 40, 40, 10, true);

    // Now create the cactus graph and test that the filter classifies
    // these blocks (which are guaranteed to end up in separate
    // chains) correctly.
    stCactusNode *startCactusNode;
    stPinchComponent *deadEndComponent;
    stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 0, 0,
                                                                  0.0, true, 0);

    // Test block A, which should fail.
    stPinchBlock *blockA = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 10));
    stCactusEdgeEnd *chainA = getChainEndFromBlock(cactusGraph, blockA);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup(chainA, flower) == false);

    // Test block C, which should pass.
    stPinchBlock *blockC = stPinchSegment_getBlock(stPinchThread_getSegment(ingroup1Thread1, 40));
    stCactusEdgeEnd *chainC = getChainEndFromBlock(cactusGraph, blockC);
    CuAssertTrue(testCase, stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup(chainC, flower) == true);

    stPinchThreadSet_destruct(threadSet);

    teardown(testCase);
}

static void checkCycleFree(CuTest *testCase, stPinchThread *thread) {
    stPinchSegment *segment = stPinchThread_getFirst(thread);
    while (segment != NULL) {
        stPinchBlock *block = stPinchSegment_getBlock(segment);
        if (block != NULL) {
            stPinchBlockIt it = stPinchBlock_getSegmentIterator(block);
            stPinchSegment *segment2;
            while ((segment2 = stPinchBlockIt_getNext(&it)) != NULL) {
                if (segment2 != segment) {
                    CuAssertTrue(testCase, stPinchSegment_getThread(segment2) != thread);
                }
            }
        }
        segment = stPinchSegment_get3Prime(segment);
    }
}

static void testHGVMFiltering(CuTest *testCase) {
    for (int64_t testNum = 0; testNum < 50; testNum++) {
        int64_t numRandomPinches = 500;
        setup(testCase, false);
        Name ingroup1Seq1 = addThreadToFlower(flower, ingroup1, 100);
        Name ingroup1Seq2 = addThreadToFlower(flower, ingroup1, 100);
        addThreadToFlower(flower, ingroup1, 100);
        addThreadToFlower(flower, ingroup1, 100);
        addThreadToFlower(flower, ingroup2, 100);

        stPinchThreadSet *threadSet = stCaf_setup(flower);

        stSet *isolatedThreads = stSet_construct();
        stSet_insert(isolatedThreads, stPinchThreadSet_getThread(threadSet, ingroup1Seq1));
        stSet_insert(isolatedThreads, stPinchThreadSet_getThread(threadSet, ingroup1Seq2));
        stCaf_setThreadsToBeCycleFreeIsolatedComponents(threadSet, isolatedThreads);

        for (int64_t i = 0; i < numRandomPinches; i++) {
            stPinch pinch = stPinchThreadSet_getRandomPinch(threadSet);
            stPinchThread_filterPinch(stPinchThreadSet_getThread(threadSet, pinch.name1),
                                      stPinchThreadSet_getThread(threadSet, pinch.name2),
                                      pinch.start1,
                                      pinch.start2,
                                      pinch.length,
                                      pinch.strand,
                                      (bool(*)(stPinchSegment *, stPinchSegment *, void *))stCaf_filterToEnsureCycleFreeIsolatedComponents, flower);
            stSortedSet *threadComponents = stPinchThreadSet_getThreadComponents(threadSet);
            CuAssertTrue(testCase, stSortedSet_size(threadComponents) >= 2);
        }

        checkCycleFree(testCase, stPinchThreadSet_getThread(threadSet, ingroup1Seq1));
        checkCycleFree(testCase, stPinchThreadSet_getThread(threadSet, ingroup1Seq2));

        stSortedSet *threadComponents = stPinchThreadSet_getThreadComponents(threadSet);
        CuAssertTrue(testCase, stSortedSet_size(threadComponents) >= 2);

        stSet_destruct(isolatedThreads);
        stSortedSet_destruct(threadComponents);
        stPinchThreadSet_destruct(threadSet);

        teardown(testCase);
    }
}

// Random pinches through stCaf_filterByRepeatSpecies (<caf secondaryAlignmentFilter="sharedSpecies">) must never put
// two segments of one species in a block, while still building blocks of several species.
static void testRepeatSpeciesFiltering(CuTest *testCase) {
    for (int64_t testNum = 0; testNum < 50; testNum++) {
        int64_t numRandomPinches = 500;
        setup(testCase, true);
        addThreadToFlower(flower, ingroup1, 100);
        addThreadToFlower(flower, ingroup1, 100);
        addThreadToFlower(flower, ingroup2, 100);
        addThreadToFlower(flower, ingroup2, 100);
        addThreadToFlower(flower, outgroup1, 100);
        addThreadToFlower(flower, outgroup2, 100);

        stPinchThreadSet *threadSet = stCaf_setup(flower);

        for (int64_t i = 0; i < numRandomPinches; i++) {
            stPinch pinch = stPinchThreadSet_getRandomPinch(threadSet);
            stPinchThread_filterPinch(stPinchThreadSet_getThread(threadSet, pinch.name1),
                                      stPinchThreadSet_getThread(threadSet, pinch.name2),
                                      pinch.start1,
                                      pinch.start2,
                                      pinch.length,
                                      pinch.strand,
                                      (bool(*)(stPinchSegment *, stPinchSegment *, void *))stCaf_filterByRepeatSpecies, flower);
        }

        int64_t maxDegree = 0;
        stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
        stPinchBlock *block;
        while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
            stSet *events = stSet_construct();
            stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(block);
            stPinchSegment *segment;
            while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
                Event *event = stCaf_getEvent(segment, flower);
                CuAssertTrue(testCase, stSet_search(events, event) == NULL);
                stSet_insert(events, event);
            }
            stSet_destruct(events);
            if (stPinchBlock_getDegree(block) > maxDegree) {
                maxDegree = stPinchBlock_getDegree(block);
            }
        }
        CuAssertTrue(testCase, maxDegree > 2);

        stPinchThreadSet_destruct(threadSet);
        teardown(testCase);
    }
}

static stSet *getEventSet(stPinchSegment *segment) {
    stSet *events = stSet_construct();
    stPinchBlock *block = stPinchSegment_getBlock(segment);
    if (block == NULL) {
        stSet_insert(events, stCaf_getEvent(segment, flower));
        return events;
    }
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(block);
    while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        stSet_insert(events, stCaf_getEvent(segment, flower));
    }
    return events;
}

static CuTest *repeatSpeciesTestCase;

// stCaf_filterByRepeatSpecies, checked against the intersection of the two sides' sets of events
static bool checkedFilterByRepeatSpecies(stPinchSegment *segment1, stPinchSegment *segment2, Flower *flower) {
    stSet *events1 = getEventSet(segment1), *events2 = getEventSet(segment2);
    stSet *shared = stSet_getIntersection(events1, events2);
    bool expected = stSet_size(shared) > 0;
    stSet_destruct(events1);
    stSet_destruct(events2);
    stSet_destruct(shared);
    bool b = stCaf_filterByRepeatSpecies(segment1, segment2, flower);
    CuAssertTrue(repeatSpeciesTestCase, b == expected);
    CuAssertTrue(repeatSpeciesTestCase, stCaf_filterByRepeatSpecies(segment2, segment1, flower) == expected);
    return b;
}

// stCaf_filterByMultipleSpecies, checked against the sets of events: it caches its answer in the blocks' flags,
// which stCaf_filterByRepeatSpecies also sets
static bool checkedFilterByMultipleSpecies(stPinchSegment *segment1, stPinchSegment *segment2, Flower *flower) {
    bool expected = false;
    if (stPinchSegment_getBlock(segment1) != NULL && stPinchSegment_getBlock(segment2) != NULL) {
        stSet *events1 = getEventSet(segment1), *events2 = getEventSet(segment2);
        if (stPinchSegment_getBlock(segment1) == stPinchSegment_getBlock(segment2)) {
            expected = stPinchBlock_getLength(stPinchSegment_getBlock(segment1)) > 1 && stSet_size(events1) > 1;
        } else {
            expected = stSet_size(events1) > 1 && stSet_size(events2) > 1;
        }
        stSet_destruct(events1);
        stSet_destruct(events2);
    }
    bool b = stCaf_filterByMultipleSpecies(segment1, segment2, flower);
    CuAssertTrue(repeatSpeciesTestCase, b == expected);
    return b;
}

// stCaf_filterByRepeatSpecies has to have given each of two different blocks of degree 16
// (REPEAT_SPECIES_MIN_CACHED_DEGREE) or more a mask, which clears its modified flag, rather than scanned it
static void checkRepeatSpeciesMasks(stPinchBlock *block1, stPinchBlock *block2) {
    if (block1 != block2) {
        CuAssertTrue(repeatSpeciesTestCase, block1 == NULL || stPinchBlock_getDegree(block1) < 16 ||
                                            !stPinchBlock_getModifiedFlag(block1));
        CuAssertTrue(repeatSpeciesTestCase, block2 == NULL || stPinchBlock_getDegree(block2) < 16 ||
                                            !stPinchBlock_getModifiedFlag(block2));
    }
}

// stCaf_filterByRepeatSpecies keeps the events of each block of high degree between calls. It must still answer as
// the sets of events do while big blocks of repeat copies are pinched together, split and destroyed (their addresses
// then going to new blocks) in between, and while stCaf_filterByMultipleSpecies shares the blocks' flags. Half the
// tests have more events than fit in one word of a mask.
static void testRepeatSpeciesFilteringOfBigBlocks(CuTest *testCase) {
    repeatSpeciesTestCase = testCase;
    for (int64_t testNum = 0; testNum < 8; testNum++) {
        setup(testCase, true);
        int64_t speciesNumber = testNum % 4 < 2 ? 10 : 70, copies = 20, period = 10;
        Event *events[70] = { ingroup1, ingroup2, outgroup1, outgroup2 };
        for (int64_t i = 4; i < speciesNumber; i++) {
            char *name = stString_print("ingroup%" PRIi64, i);
            events[i] = event_construct3(name, 0.2, ancestor, eventTree);
            free(name);
        }
        for (int64_t i = 0; i < 2 * speciesNumber; i++) {
            addThreadToFlower(flower, events[i % speciesNumber], copies * period);
        }
        stPinchThreadSet *threadSet = stCaf_setup(flower);
        stList *threads = stList_construct(); // thread i is a copy of species i % speciesNumber
        for (int64_t i = 0; i < 2 * speciesNumber; i++) {
            stList_append(threads, NULL);
        }
        stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
        stPinchThread *thread;
        while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
            Event *event = stCaf_getEvent(stPinchThread_getFirst(thread), flower);
            int64_t i = 0;
            while (events[i] != event) {
                i++;
            }
            stList_set(threads, stList_get(threads, i) == NULL ? i : i + speciesNumber, thread);
        }
        int64_t maxDegree = 0;
        int64_t pinchNumber = 400 * speciesNumber; // so that each species gets big blocks
        for (int64_t i = 0; i < pinchNumber; i++) {
            // the same interval of two random copies of the repeat (or its reverse). The unfiltered pinches are
            // within a species, as a repeat family's primary alignments would be, and come first; then come the
            // secondaries between species, with some more unfiltered pinches, splits and block deaths among them
            double r = i < pinchNumber / 4 ? 0.45 * st_random() :
                       (i < 3 * pinchNumber / 8 ? 0.85 + 0.15 * st_random() : 0.15 + 0.85 * st_random());
            int64_t thread1Index = st_randomInt(0, 2 * speciesNumber);
            int64_t thread2Index = r < 0.45 ? (thread1Index + speciesNumber) % (2 * speciesNumber) :
                                   st_randomInt(0, 2 * speciesNumber);
            stPinchThread *thread1 = stList_get(threads, thread1Index), *thread2 = stList_get(threads, thread2Index);
            int64_t offset = st_randomInt(0, period), length = st_randomInt(1, period - offset + 1);
            int64_t start1 = stPinchThread_getStart(thread1) + 1 + st_randomInt(0, copies) * period + offset;
            bool strand = st_random() < 0.7;
            int64_t start2 = stPinchThread_getStart(thread2) + 1 + st_randomInt(0, copies) * period +
                             (strand ? offset : period - offset - length);
            if (r < 0.45) {
                stPinchThread_pinch(thread1, thread2, start1, start2, length, strand);
            } else if (r < 0.75) {
                stPinchThread_filterPinch(thread1, thread2, start1, start2, length, strand,
                                          (bool(*)(stPinchSegment *, stPinchSegment *, void *))checkedFilterByRepeatSpecies, flower);
            } else if (r < 0.85) {
                stPinchThread_filterPinch(thread1, thread2, start1, start2, length, strand,
                                          (bool(*)(stPinchSegment *, stPinchSegment *, void *))checkedFilterByMultipleSpecies, flower);
            } else if (r < 1.0 - 0.2 / speciesNumber) { // about 60 passes over all the blocks, whatever the species
                if (start1 < stPinchThread_getStart(thread1) + stPinchThread_getLength(thread1) - 1) {
                    stPinchThread_split(thread1, start1); // the block is split, whichever way its segment is oriented
                }
                checkedFilterByRepeatSpecies(stPinchThread_getSegment(thread1, start1),
                                             stPinchThread_getSegment(thread2, start2), flower);
            } else {
                stList *blocks = stList_construct();
                stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
                stPinchBlock *block;
                while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
                    if ((int64_t)stPinchBlock_getDegree(block) > maxDegree) {
                        maxDegree = stPinchBlock_getDegree(block);
                    }
                    stList_append(blocks, block);
                }
                // every block against another, each filter both before and after the other, as either may clear
                // the modified flag of a block that has changed
                for (int64_t k = 0; k < 3; k++) {
                    for (int64_t j = 0; j < stList_length(blocks); j++) {
                        stPinchSegment *segment1 = stPinchBlock_getFirst(stList_get(blocks, j));
                        stPinchSegment *segment2 = k == 1 ? segment1 :
                            stPinchBlock_getFirst(stList_get(blocks, st_randomInt(0, stList_length(blocks))));
                        if (k == 2) {
                            checkedFilterByMultipleSpecies(segment1, segment2, flower);
                        }
                        if (k != 1) {
                            checkedFilterByRepeatSpecies(segment1, segment2, flower);
                            checkRepeatSpeciesMasks(stPinchSegment_getBlock(segment1), stPinchSegment_getBlock(segment2));
                        } else {
                            checkedFilterByMultipleSpecies(segment1, segment1, flower);
                        }
                    }
                }
                for (int64_t j = 0; j < stList_length(blocks); j++) {
                    if (st_random() < 0.1) {
                        stPinchBlock_destruct(stList_get(blocks, j));
                    }
                }
                stList_destruct(blocks);
            }
        }
        CuAssertTrue(testCase, maxDegree >= 16);
        stList_destruct(threads);
        stPinchThreadSet_destruct(threadSet);
        stCaf_clearRepeatSpeciesCache();
        teardown(testCase);
    }
}

// Events that join the flower after stCaf_filterByRepeatSpecies sized its masks can leave it short of bits. It must
// then scan the blocks instead, and use masks again, with a word more, once its cache is cleared.
static void testRepeatSpeciesFilteringOfLateEvents(CuTest *testCase) {
    repeatSpeciesTestCase = testCase;
    setup(testCase, true);
    int64_t copies = 80, period = 10, lateEventNumber = 70;
    Name repeatName = addThreadToFlower(flower, ingroup1, copies * period);
    Name loneName = addThreadToFlower(flower, ingroup2, period);
    Name outgroupName = addThreadToFlower(flower, outgroup1, period);
    stPinchThreadSet *threadSet = stCaf_setup(flower);
    // a block of all the copies of a repeat of ingroup1, and two segments in no block
    stPinchThread *repeat = stPinchThreadSet_getThread(threadSet, repeatName);
    int64_t start = stPinchThread_getStart(repeat) + 1;
    for (int64_t i = 1; i < copies; i++) {
        stPinchThread_pinch(repeat, repeat, start, start + i * period, period, true);
    }
    stPinchSegment *repeatSegment = stPinchThread_getSegment(repeat, start);
    stPinchThread *lone = stPinchThreadSet_getThread(threadSet, loneName);
    stPinchSegment *loneSegment = stPinchThread_getSegment(lone, stPinchThread_getStart(lone) + 1);
    stPinchThread *outgroup = stPinchThreadSet_getThread(threadSet, outgroupName);
    stPinchSegment *outgroupSegment = stPinchThread_getSegment(outgroup, stPinchThread_getStart(outgroup) + 1);
    CuAssertTrue(testCase, !checkedFilterByRepeatSpecies(loneSegment, repeatSegment, flower)); // sizes the masks
    // a smaller block of a copy each of 70 new species, whose threads are added as stCaf_setup would, but without
    // their caps as user data, so that the filter looks them up in the flower
    stPinchThread *late = NULL;
    for (int64_t i = 0; i < lateEventNumber; i++) {
        char *name = stString_print("late%" PRIi64, i);
        Cap *cap = flower_getCap(flower, addThreadToFlower(flower, event_construct3(name, 0.2, ancestor, eventTree),
                                                           period));
        free(name);
        stPinchThread *thread = stPinchThreadSet_addThread(threadSet, cap_getName(cap), cap_getCoordinate(cap),
                                                           period + 2);
        if (late == NULL) {
            late = thread;
        } else {
            stPinchThread_pinch(late, thread, stPinchThread_getStart(late) + 1, stPinchThread_getStart(thread) + 1,
                                period, true);
        }
    }
    stPinchSegment *lateSegment = stPinchThread_getSegment(late, stPinchThread_getStart(late) + 1);
    // the late block runs the filter out of bits while it holds the repeat's mask, then it scans
    CuAssertTrue(testCase, !checkedFilterByRepeatSpecies(lateSegment, repeatSegment, flower));
    CuAssertTrue(testCase, !checkedFilterByRepeatSpecies(loneSegment, lateSegment, flower));
    CuAssertTrue(testCase, !checkedFilterByRepeatSpecies(repeatSegment, outgroupSegment, flower));
    // the blocks are joined through both filters
    stPinchThread_filterPinch(lone, late, stPinchThread_getStart(lone) + 1, stPinchThread_getStart(late) + 1, period,
                              true, (bool(*)(stPinchSegment *, stPinchSegment *, void *))checkedFilterByRepeatSpecies,
                              flower);
    CuAssertTrue(testCase, !checkedFilterByMultipleSpecies(lateSegment, repeatSegment, flower));
    stPinchThread_filterPinch(repeat, late, start, stPinchThread_getStart(late) + 1, period, true,
                              (bool(*)(stPinchSegment *, stPinchSegment *, void *))checkedFilterByRepeatSpecies, flower);
    CuAssertTrue(testCase, stPinchSegment_getBlock(lateSegment) == stPinchSegment_getBlock(repeatSegment));
    CuAssertTrue(testCase, checkedFilterByMultipleSpecies(lateSegment, repeatSegment, flower));
    CuAssertTrue(testCase, !checkedFilterByRepeatSpecies(outgroupSegment, lateSegment, flower));
    // once cleared, the masks are sized for every event and the joined block gets one
    stCaf_clearRepeatSpeciesCache();
    CuAssertTrue(testCase, !checkedFilterByRepeatSpecies(outgroupSegment, lateSegment, flower));
    checkRepeatSpeciesMasks(NULL, stPinchSegment_getBlock(lateSegment));
    stPinchThreadSet_destruct(threadSet);
    stCaf_clearRepeatSpeciesCache();
    teardown(testCase);
}

CuSuite* filteringTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testChainHasUnequalNumberOfIngroupCopies);
    SUITE_ADD_TEST(suite, testChainHasUnequalNumberOfIngroupCopiesOrNoOutgroup);
    SUITE_ADD_TEST(suite, testChainHasUnequalNumberOfIngroupCopiesOrNoOutgroup_noOutgroups);
    SUITE_ADD_TEST(suite, testHGVMFiltering);
    SUITE_ADD_TEST(suite, testRepeatSpeciesFiltering);
    SUITE_ADD_TEST(suite, testRepeatSpeciesFilteringOfBigBlocks);
    SUITE_ADD_TEST(suite, testRepeatSpeciesFilteringOfLateEvents);
    return suite;
}
