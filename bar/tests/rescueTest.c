#include "CuTest.h"
#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "rescue.h"

static void writeCoverageArrayToBed(int64_t name, bool *coverageArray,
                                    int64_t length, FILE *bedFile) {
    bool inCoveredRegion = false;
    for (int64_t i = 0; i < length; i++) {
        if (coverageArray[i] && !inCoveredRegion) {
            fprintf(bedFile, "%" PRIi64 "\t%" PRIi64 "\t", name, i);
            inCoveredRegion = true;
        } else if (!coverageArray[i] && inCoveredRegion) {
            fprintf(bedFile, "%" PRIi64 "\n", i);
            inCoveredRegion = false;
        }
    }
    if (inCoveredRegion) {
        fprintf(bedFile, "%" PRIi64 "\n", length);
    }
}

// Just check that running a rescue on pinch threads works correctly.
static void test_rescueRandomSequences(CuTest *testCase) {
    for (int64_t testNum = 0; testNum < 1000; testNum++) {
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomGraph();

        // First, for each thread in the set, generate a random coverage
        // array (independent of the coverage already in the random
        // graph).
        stHash *coveragesToRescue = stHash_construct2(NULL, free);
        char *bedFilePath = getTempFile();
        FILE *bedFile = fopen(bedFilePath, "w");
        stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
        stPinchThread *thread;
        while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
            printf("outgroup coverage for %" PRIi64 " (start %" PRIi64 ")\n", stPinchThread_getName(thread), stPinchThread_getStart(thread));
            int64_t threadStart = stPinchThread_getStart(thread);
            int64_t threadLen = stPinchThread_getLength(thread);
            bool *coverageArray = st_calloc(threadStart + threadLen, sizeof(bool));
            for (int64_t i = threadStart; i < threadStart + threadLen; i++) {
                coverageArray[i] = st_random() < 0.3;
                if(coverageArray[i]) {
                    printf("1");
                } else {
                    printf("0");
                }
            }
            printf("\n");
            stHash_insert(coveragesToRescue, thread, coverageArray);
            writeCoverageArrayToBed(stPinchThread_getName(thread), coverageArray, threadStart + threadLen, bedFile);
        }

        // Next, go through the graph and mark down those regions that
        // are already covered before rescue.
        stHash *regionsAlreadyCovered = stHash_construct2(NULL, free);
        threadIt = stPinchThreadSet_getIt(threadSet);
        while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
            printf("existing for %" PRIi64 "\n", stPinchThread_getName(thread));
            int64_t threadStart = stPinchThread_getStart(thread);
            int64_t threadLen = stPinchThread_getLength(thread);
            bool *coverageArray = st_calloc(threadStart + threadLen, sizeof(bool));
            stPinchSegment *segment = stPinchThread_getFirst(thread);
            while (segment != NULL) {
                int64_t start = stPinchSegment_getStart(segment);
                int64_t end = start + stPinchSegment_getLength(segment);
                assert(end - start <= threadLen);
                if (stPinchSegment_getBlock(segment) != NULL) {
                    for (int64_t i = start; i < end; i++) {
                        coverageArray[i] = 1;
                        printf("1");
                    }
                } else {
                    for (int64_t i = start; i < end; i++) {
                        printf("0");
                    }
                }
                segment = stPinchSegment_get3Prime(segment);
            }
            stHash_insert(regionsAlreadyCovered, thread, coverageArray);
            printf("\n");
        }

        // Sort the coverage BED.
        fclose(bedFile);
        char *sortedBedFilePath = getTempFile();
        printf("bed path: %s sorted bed path: %s\n",
               bedFilePath,
               sortedBedFilePath);
        char *cmd = stString_print("sort -k1n -k2n -k3n %s -o %s", bedFilePath,
                                   sortedBedFilePath);
        st_system(cmd);
        free(cmd);
        bedFile = fopen(sortedBedFilePath, "r");

        // Run the rescue and make sure it worked.
        threadIt = stPinchThreadSet_getIt(threadSet);
        bedRegion *curBedLine = readNextBedLine(bedFile);
        while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
            printf("rescued %" PRIi64 "\n", stPinchThread_getName(thread));
            int64_t threadStart = stPinchThread_getStart(thread);
            int64_t threadLen = stPinchThread_getLength(thread);
            bool *coverageArray = stHash_search(coveragesToRescue, thread);
            assert(coverageArray != NULL);
            bool *alreadyCovered = stHash_search(regionsAlreadyCovered, thread);
            assert(alreadyCovered != NULL);
            curBedLine = rescueCoveredRegions(thread, bedFile, curBedLine);
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
                        printf("1");
                        // Make sure we haven't screwed up and made
                        // blocks where we shouldn't have.
                        CuAssertTrue(testCase, alreadyCovered[i] == 1 || coverageArray[i] == 1);
                    } else {
                        printf("0");
                        // Make sure we haven't taken away from the
                        // coverage that already existed or missed
                        // regions we were supposed to rescue.
                        CuAssertTrue(testCase, alreadyCovered[i] == 0 && coverageArray[i] == 0);
                    }
                }
                segment = stPinchSegment_get3Prime(segment);
            }
            printf("\n");
        }

        stHash_destruct(coveragesToRescue);
        stHash_destruct(regionsAlreadyCovered);
        stPinchThreadSet_destruct(threadSet);
        stFile_rmrf(bedFilePath);
        stFile_rmrf(sortedBedFilePath);
        free(bedFilePath);
        free(sortedBedFilePath);
    }
}

CuSuite *rescueTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_rescueRandomSequences);
    return suite;
}
