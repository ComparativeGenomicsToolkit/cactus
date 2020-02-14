#include "CuTest.h"
#include "sonLib.h"
#include "stCaf.h"
#include "stPinchGraphs.h"

// Test a simple case where two threads only align in a single
// chain--this should obviously not be recoverable.
static void testDoesNotRemoveIsolatedChain(CuTest *testCase) {
    CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
    eventTree_construct2(cactusDisk);
    Flower *flower = flower_construct2(0, cactusDisk);
    group_construct2(flower);
    flower_check(flower);

    Name thread1Name = testCommon_addThreadToFlower(flower, "one", 100);
    Name thread2Name = testCommon_addThreadToFlower(flower, "two", 100);
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
    stCaf_meltRecoverableChains(flower, threadSet, true, 1000, NULL, 1, INT64_MAX);
    // It shouldn't've removed any blocks
    CuAssertIntEquals(testCase, 7, stPinchThreadSet_getTotalBlockNumber(threadSet));

    stPinchThreadSet_destruct(threadSet);
    testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
}

// If there are three threads, two of which share an indel, that indel
// should be considered recoverable.
static void testRemovesIndel(CuTest *testCase) {
    CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
    eventTree_construct2(cactusDisk);
    Flower *flower = flower_construct2(0, cactusDisk);
    group_construct2(flower);
    flower_check(flower);

    Name thread1Name = testCommon_addThreadToFlower(flower, "one", 100);
    Name thread2Name = testCommon_addThreadToFlower(flower, "two", 100);
    Name thread3Name = testCommon_addThreadToFlower(flower, "three", 100);
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
    stCaf_meltRecoverableChains(flower, threadSet, true, 1000, NULL, 1, INT64_MAX);
    // One block should be missing
    CuAssertIntEquals(testCase, 8, stPinchThreadSet_getTotalBlockNumber(threadSet));
    // It should be the middle one (40-50)
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 45)) == NULL);

    stPinchThreadSet_destruct(threadSet);
    testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
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
    CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
    eventTree_construct2(cactusDisk);
    Flower *flower = flower_construct2(0, cactusDisk);
    group_construct2(flower);
    flower_check(flower);

    Name thread1Name = testCommon_addThreadToFlower(flower, "one", 100);
    Name thread2Name = testCommon_addThreadToFlower(flower, "two", 100);
    Name thread3Name = testCommon_addThreadToFlower(flower, "three", 100);
    Name thread4Name = testCommon_addThreadToFlower(flower, "four", 100);
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
    stCaf_meltRecoverableChains(flower, threadSet, true, 1000, NULL, 1, INT64_MAX);
    // Four blocks should have been removed
    CuAssertIntEquals(testCase, 9, stPinchThreadSet_getTotalBlockNumber(threadSet));
    // The middle block (30-40) should be the one that still exists
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 35)) != NULL);

    stPinchThreadSet_destruct(threadSet);
    testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
}

CuSuite *recoverableChainsTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testDoesNotRemoveIsolatedChain);
    SUITE_ADD_TEST(suite, testRemovesIndel);
    SUITE_ADD_TEST(suite, testRecoverableTelomereAdjacentChainsNotKept);
    return suite;
}
