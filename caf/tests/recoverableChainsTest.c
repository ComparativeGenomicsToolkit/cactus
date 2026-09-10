#define _DEFAULT_SOURCE //for setenv
#include <stdlib.h>
#include "CuTest.h"
#include "sonLib.h"
#include "stCaf.h"
#include "stPinchGraphs.h"

// Test a simple case where two threads only align in a single
// chain--this should obviously not be recoverable.
static void testDoesNotRemoveIsolatedChain(CuTest *testCase) {
    CactusDisk *cactusDisk = cactusDisk_construct();
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
    cactusDisk_destruct(cactusDisk);
}

// If there are three threads, two of which share an indel, that indel
// should be considered recoverable.
static void testRemovesIndel(CuTest *testCase) {
    CactusDisk *cactusDisk = cactusDisk_construct();
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
    cactusDisk_destruct(cactusDisk);
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
    CactusDisk *cactusDisk = cactusDisk_construct();
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
    cactusDisk_destruct(cactusDisk);
}

// stCaf_meltChains joins trivial boundaries only where its deletions reach, on the argument that a
// chain melt cannot create any elsewhere. With CACTUS_CAF_CHECK_JOIN set it runs the full join after
// its own and aborts if that changes anything, so random graphs melted that way check the argument,
// including the boundaries at thread ends that stCaf_ensureEndsAreDistinct splits after every join.
static void testMeltChainsLeavesNoTrivialBoundary(CuTest *testCase) {
    setenv("CACTUS_CAF_CHECK_JOIN", "1", 1);
    for (int64_t test = 0; test < 100; test++) {
        CactusDisk *cactusDisk = cactusDisk_construct();
        eventTree_construct2(cactusDisk);
        Flower *flower = flower_construct2(0, cactusDisk);
        group_construct2(flower);
        int64_t threadNumber = st_randomInt(2, 6);
        Name threadNames[6];
        for (int64_t i = 0; i < threadNumber; i++) {
            char *header = stString_print("thread%" PRIi64, i);
            threadNames[i] = testCommon_addThreadToFlower(flower, header, st_randomInt(20, 200));
            free(header);
        }
        stPinchThreadSet *threadSet = stCaf_setup(flower);
        int64_t pinchNumber = st_randomInt(1, 60);
        for (int64_t i = 0; i < pinchNumber; i++) {
            stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, threadNames[st_randomInt(0, threadNumber)]);
            stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, threadNames[st_randomInt(0, threadNumber)]);
            // Leave the first and last base of each thread alone, as real alignments do: those are the
            // stub ends that stand for the caps
            int64_t start1 = st_randomInt(stPinchThread_getStart(thread1) + 1, stPinchThread_getStart(thread1) + stPinchThread_getLength(thread1) - 1);
            int64_t start2 = st_randomInt(stPinchThread_getStart(thread2) + 1, stPinchThread_getStart(thread2) + stPinchThread_getLength(thread2) - 1);
            int64_t maxLength = stPinchThread_getStart(thread1) + stPinchThread_getLength(thread1) - 1 - start1;
            int64_t maxLength2 = stPinchThread_getStart(thread2) + stPinchThread_getLength(thread2) - 1 - start2;
            if (maxLength2 < maxLength) {
                maxLength = maxLength2;
            }
            stPinchThread_pinch(thread1, thread2, start1, start2, st_randomInt(0, maxLength + 1), st_random() > 0.5);
        }
        stCaf_joinTrivialBoundaries(threadSet); //as after annealing
        int64_t blockNumber = stPinchThreadSet_getTotalBlockNumber(threadSet);
        stCaf_meltChains(flower, threadSet, 2, 0, INT64_MAX);
        int64_t destroyed = stCaf_meltChains(flower, threadSet, st_randomInt(3, 40), 1, st_randomInt(1, 100));
        CuAssertTrue(testCase, stPinchThreadSet_getTotalBlockNumber(threadSet) <= blockNumber - destroyed);
        stPinchThreadSet_destruct(threadSet);
        cactusDisk_destruct(cactusDisk);
    }
    unsetenv("CACTUS_CAF_CHECK_JOIN");
}

CuSuite *recoverableChainsTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testDoesNotRemoveIsolatedChain);
    SUITE_ADD_TEST(suite, testRemovesIndel);
    SUITE_ADD_TEST(suite, testRecoverableTelomereAdjacentChainsNotKept);
    SUITE_ADD_TEST(suite, testMeltChainsLeavesNoTrivialBoundary);
    return suite;
}
