#include "CuTest.h"
#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "rescue.h"

// Check that reading bed files works as expected.
static void test_getIngroupCoverage(CuTest *testCase) {

}

// Just check that running a rescue on pinch threads works correctly.
static void test_rescueRandomSequences(CuTest *testCase) {
    for (int64_t testNum = 0; testNum < 1000; testNum++) {
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomGraph();

        // First, for each thread in the set, generate a random coverage
        // array (independent of the coverage already in the random
        // graph).
        stHash *coveragesToRescue = stHash_construct2(NULL, free);
        stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
        stPinchThread *thread;
        while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
            int64_t threadStart = stPinchThread_getStart(thread);
            int64_t threadLen = stPinchThread_getLength(thread);
            bool *coverageArray = st_calloc(threadStart + threadLen, sizeof(bool));
            for (int64_t i = threadStart; i < threadLen; i++) {
                coverageArray[i] = st_random() < 0.3;
            }
            stHash_insert(coveragesToRescue, thread, coverageArray);
        }

        // Next, go through the graph and mark down those regions that
        // are already covered before rescue.
        stHash *regionsAlreadyCovered = stHash_construct2(NULL, free);
        threadIt = stPinchThreadSet_getIt(threadSet);
        while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
            int64_t threadStart = stPinchThread_getStart(thread);
            int64_t threadLen = stPinchThread_getLength(thread);
            bool *coverageArray = st_calloc(threadStart + threadLen, sizeof(bool));
            stPinchSegment *segment = stPinchThread_getFirst(thread);
            while (segment != NULL) {
                if (stPinchSegment_getBlock(segment) != NULL) {
                    int64_t start = stPinchSegment_getStart(segment);
                    int64_t end = start + stPinchSegment_getLength(segment);
                    assert(end - start <= threadLen);
                    for (int64_t i = start; i < end; i++) {
                        coverageArray[i] = 1;
                    }
                }
                segment = stPinchSegment_get3Prime(segment);
            }
            stHash_insert(regionsAlreadyCovered, thread, coverageArray);
        }

        // Run the rescue and make sure it worked.
        threadIt = stPinchThreadSet_getIt(threadSet);
        while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
            int64_t threadStart = stPinchThread_getStart(thread);
            int64_t threadLen = stPinchThread_getLength(thread);
            bool *coverageArray = stHash_search(coveragesToRescue, thread);
            assert(coverageArray != NULL);
            bool *alreadyCovered = stHash_search(regionsAlreadyCovered, thread);
            assert(alreadyCovered != NULL);
            rescueCoveredRegions(thread, coverageArray);
            stPinchSegment *segment = stPinchThread_getFirst(thread);
            while (segment != NULL) {
                int64_t start = stPinchSegment_getStart(segment);
                int64_t end = start + stPinchSegment_getLength(segment);
                assert(end <= threadStart + threadLen);
                bool covered = false;
                if (stPinchSegment_getBlock(segment) != NULL) {
                    covered = true;
                }
                for (int64_t i = start; i < end; i++) {
                    if (covered) {
                        // Make sure we haven't screwed up and made
                        // blocks where we shouldn't have.
                        CuAssertTrue(testCase, alreadyCovered[i] == 1 || coverageArray[i] == 1);
                    } else {
                        // Make sure we haven't taken away from
                        // the coverage that already existed.
                        CuAssertTrue(testCase, alreadyCovered[i] == 0 && coverageArray[i] == 0);
                    }
                }
                segment = stPinchSegment_get3Prime(segment);
            }
        }

        stHash_destruct(coveragesToRescue);
        stHash_destruct(regionsAlreadyCovered);
        stPinchThreadSet_destruct(threadSet);
    }
}

CuSuite *rescueTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_getIngroupCoverage);
    SUITE_ADD_TEST(suite, test_rescueRandomSequences);
    return suite;
}
