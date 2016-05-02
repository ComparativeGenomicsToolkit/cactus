#include "CuTest.h"
#include "sonLib.h"
#include "stCaf.h"
#include "stPinchGraphs.h"

// Adds a thread with random nucleotides to the flower, and return its corresponding name in the pinch graph.
static Name addThreadToFlower(Flower *flower, int64_t length) {
    char *dna = stRandom_getRandomDNAString(length, true, true, true);
    MetaSequence *metaSequence = metaSequence_construct(2, length, dna, "", 0, flower_getCactusDisk(flower));
    Sequence *sequence = sequence_construct(metaSequence, flower);

    End *end1 = end_construct2(0, 0, flower);
    End *end2 = end_construct2(1, 0, flower);
    Cap *cap1 = cap_construct2(end1, 1, 1, sequence);
    Cap *cap2 = cap_construct2(end2, length + 2, 1, sequence);
    cap_makeAdjacent(cap1, cap2);

    free(dna);
    return cap_getName(cap1);
}

// Test a simple case where two threads only align in a single
// chain--this should obviously not be recoverable.
static void testDoesNotRemoveIsolatedChain(CuTest *testCase) {
    CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk();
    Flower *flower = flower_construct2(0, cactusDisk);
    group_construct2(flower);
    flower_check(flower);

    Name thread1Name = addThreadToFlower(flower, 100);
    Name thread2Name = addThreadToFlower(flower, 100);
    stPinchThreadSet *threadSet = stCaf_setup(flower);
    stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, thread1Name);
    stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, thread2Name);
    // Three 10-bp blocks separated by 20 bp
    stPinchThread_pinch(thread1, thread2, 10, 10, 10, true);
    stPinchThread_pinch(thread1, thread2, 40, 40, 10, true);
    stPinchThread_pinch(thread1, thread2, 70, 70, 10, true);
    // There should now be 7 blocks -- one for each cap, and the blocks we just added
    CuAssertIntEquals(testCase, 7, stPinchThreadSet_getTotalBlockNumber(threadSet));
    // Run the remove-recoverable-chains code
    stCaf_meltRecoverableChains(flower, threadSet, true, 1000, NULL);
    // It shouldn't've removed any blocks
    CuAssertIntEquals(testCase, 7, stPinchThreadSet_getTotalBlockNumber(threadSet));

    stPinchThreadSet_destruct(threadSet);
    testCommon_deleteTemporaryCactusDisk(cactusDisk);
}

// If there are three threads, two of which share an indel, that indel
// should be considered recoverable.
static void testRemovesIndel(CuTest *testCase) {
    CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk();
    Flower *flower = flower_construct2(0, cactusDisk);
    group_construct2(flower);
    flower_check(flower);

    Name thread1Name = addThreadToFlower(flower, 100);
    Name thread2Name = addThreadToFlower(flower, 100);
    Name thread3Name = addThreadToFlower(flower, 100);
    stPinchThreadSet *threadSet = stCaf_setup(flower);
    stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, thread1Name);
    stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, thread2Name);
    stPinchThread *thread3 = stPinchThreadSet_getThread(threadSet, thread3Name);
    // Three 10-bp blocks separated by 20 bp. The two outer blocks
    // involve all three threads, but the middle one only involves
    // threads 1 and 2.
    stPinchThread_pinch(thread1, thread2, 10, 10, 10, true);
    stPinchThread_pinch(thread1, thread3, 10, 10, 10, true);
    stPinchThread_pinch(thread1, thread2, 40, 40, 10, true);
    stPinchThread_pinch(thread1, thread2, 70, 70, 10, true);
    stPinchThread_pinch(thread1, thread3, 70, 70, 10, true);
    // There should now be 9 blocks -- one for each cap, and the blocks we just added
    CuAssertIntEquals(testCase, 9, stPinchThreadSet_getTotalBlockNumber(threadSet));
    // Run the remove-recoverable-chains code
    stCaf_meltRecoverableChains(flower, threadSet, true, 1000, NULL);
    // One block should be missing
    CuAssertIntEquals(testCase, 8, stPinchThreadSet_getTotalBlockNumber(threadSet));
    // It should be the middle one (40-50)
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 45)) == NULL);

    stPinchThreadSet_destruct(threadSet);
    testCommon_deleteTemporaryCactusDisk(cactusDisk);
}

// If the alignment looks like this, with = representing aligned columns:
//
// thread 4     =
// thread 3   =-=-=
// thread 2 =-=-=-=-=
// thread 1 =-=-=-=-=
//
// then only the middle block should be kept. The older method would
// keep all blocks as all are telomere-adjacent.
static void testRecoverableTelomereAdjacentChainsNotKept(CuTest *testCase) {
    CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk();
    Flower *flower = flower_construct2(0, cactusDisk);
    group_construct2(flower);
    flower_check(flower);

    Name thread1Name = addThreadToFlower(flower, 100);
    Name thread2Name = addThreadToFlower(flower, 100);
    Name thread3Name = addThreadToFlower(flower, 100);
    Name thread4Name = addThreadToFlower(flower, 100);
    stPinchThreadSet *threadSet = stCaf_setup(flower);
    stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, thread1Name);
    stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, thread2Name);
    stPinchThread *thread3 = stPinchThreadSet_getThread(threadSet, thread3Name);
    stPinchThread *thread4 = stPinchThreadSet_getThread(threadSet, thread4Name);

    stPinchThread_pinch(thread1, thread2, 10, 10, 10, true);
    stPinchThread_pinch(thread1, thread2, 20, 20, 10, true);
    stPinchThread_pinch(thread1, thread3, 20, 20, 10, true);
    stPinchThread_pinch(thread1, thread2, 30, 30, 10, true);
    stPinchThread_pinch(thread1, thread3, 30, 30, 10, true);
    stPinchThread_pinch(thread1, thread4, 30, 30, 10, true);
    stPinchThread_pinch(thread1, thread2, 40, 40, 10, true);
    stPinchThread_pinch(thread1, thread3, 40, 40, 10, true);
    stPinchThread_pinch(thread1, thread2, 50, 50, 10, true);

    // There should now be 13 blocks -- one for each cap, and the blocks we just added
    CuAssertIntEquals(testCase, 13, stPinchThreadSet_getTotalBlockNumber(threadSet));
    // Run the remove-recoverable-chains code
    stCaf_meltRecoverableChains(flower, threadSet, true, 1000, NULL);
    // Four blocks should have been removed
    CuAssertIntEquals(testCase, 9, stPinchThreadSet_getTotalBlockNumber(threadSet));
    // The middle block (30-40) should be the one that still exists
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 35)) != NULL);

    stPinchThreadSet_destruct(threadSet);
    testCommon_deleteTemporaryCactusDisk(cactusDisk);
}

CuSuite *recoverableChainsTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testDoesNotRemoveIsolatedChain);
    SUITE_ADD_TEST(suite, testRemovesIndel);
    SUITE_ADD_TEST(suite, testRecoverableTelomereAdjacentChainsNotKept);
    return suite;
}
