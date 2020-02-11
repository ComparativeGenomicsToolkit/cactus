#include "CuTest.h"
#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "rescue.h"

// Get a little-endian bed region array as an mmaped coverage array
// would appear.
static bedRegion *getBedRegionArray(int64_t name, bool *coverageArray,
                                    int64_t length, bedRegion *array,
                                    size_t *numBeds, size_t *arraySize) {
    if (array == NULL) {
        *arraySize = 10;
        array = st_malloc(*arraySize * sizeof(bedRegion));
    }
    bool inCoveredRegion = false;
    bedRegion *curRegion = array + *numBeds;
    for (int64_t i = 0; i < length; i++) {
        if (coverageArray[i] && !inCoveredRegion) {
            curRegion->name = st_nativeInt64ToLittleEndian(name);
            curRegion->start = st_nativeInt64ToLittleEndian(i);
            inCoveredRegion = true;
        } else if (!coverageArray[i] && inCoveredRegion) {
            curRegion->stop = st_nativeInt64ToLittleEndian(i);
            inCoveredRegion = false;
            (*numBeds)++;
            if (*numBeds >= *arraySize) {
                *arraySize = *arraySize * 2 + 1;
                array = st_realloc(array, *arraySize * sizeof(bedRegion));
            }
            curRegion = array + *numBeds;
        }
    }
    if (inCoveredRegion) {
        curRegion->stop = st_nativeInt64ToLittleEndian(length);
        (*numBeds)++;
        if (*numBeds >= *arraySize) {
            *arraySize = *arraySize * 2 + 1;
            array = st_realloc(array, *arraySize * sizeof(bedRegion));
        }
        curRegion = array + *numBeds;
    }
    return array;
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
        bedRegion *bedRegionArray = NULL;
        size_t numBeds = 0, bedRegionArraySize = 0;
        while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
            st_logDebug("outgroup coverage for %" PRIi64 " (start %" PRIi64 ")\n", stPinchThread_getName(thread), stPinchThread_getStart(thread));
            int64_t threadStart = stPinchThread_getStart(thread);
            int64_t threadLen = stPinchThread_getLength(thread);
            bool *coverageArray = st_calloc(threadStart + threadLen, sizeof(bool));
            for (int64_t i = threadStart; i < threadStart + threadLen; i++) {
                coverageArray[i] = st_random() < 0.3;
                if(coverageArray[i]) {
                    st_logDebug("1");
                } else {
                    st_logDebug("0");
                }
            }
            st_logDebug("\n");
            stHash_insert(coveragesToRescue, thread, coverageArray);
            bedRegionArray = getBedRegionArray(stPinchThread_getName(thread),
                                               coverageArray,
                                               threadStart + threadLen,
                                               bedRegionArray, &numBeds,
                                               &bedRegionArraySize);
        }

        // Next, go through the graph and mark down those regions that
        // are already covered before rescue.
        stHash *regionsAlreadyCovered = stHash_construct2(NULL, free);
        threadIt = stPinchThreadSet_getIt(threadSet);
        while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
            st_logDebug("existing for %" PRIi64 "\n", stPinchThread_getName(thread));
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
                        st_logDebug("1");
                    }
                } else {
                    for (int64_t i = start; i < end; i++) {
                        st_logDebug("0");
                    }
                }
                segment = stPinchSegment_get3Prime(segment);
            }
            stHash_insert(regionsAlreadyCovered, thread, coverageArray);
            st_logDebug("\n");
        }

        // Sort the bedRegion array.
        qsort(bedRegionArray, numBeds, sizeof(bedRegion), (int (*)(const void *, const void *)) bedRegion_cmp);

        // Run the rescue and make sure it worked.
        threadIt = stPinchThreadSet_getIt(threadSet);
        while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
            st_logDebug("rescued %" PRIi64 "\n", stPinchThread_getName(thread));
            //int64_t threadStart = stPinchThread_getStart(thread);
            //int64_t threadLen = stPinchThread_getLength(thread);
            bool *coverageArray = stHash_search(coveragesToRescue, thread);
            CuAssertPtrNotNull(testCase, coverageArray);
            bool *alreadyCovered = stHash_search(regionsAlreadyCovered, thread);
            CuAssertPtrNotNull(testCase, alreadyCovered);
            rescueCoveredRegions(thread, bedRegionArray, numBeds,
                                 stPinchThread_getName(thread), 1, 0);
            stPinchSegment *segment = stPinchThread_getFirst(thread);
            while (segment != NULL) {
                int64_t start = stPinchSegment_getStart(segment);
                int64_t end = start + stPinchSegment_getLength(segment);
                bool covered = false;
                if (stPinchSegment_getBlock(segment) != NULL) {
                    covered = true;
                }
                for (int64_t i = start; i < end; i++) {
                    if (covered) {
                        st_logDebug("1");
                        // Make sure we haven't screwed up and made
                        // blocks where we shouldn't have.
                        // FIXME: skipped.
                        /* CuAssertTrue(testCase, alreadyCovered[i] == 1 || coverageArray[i] == 1); */
                    } else {
                        st_logDebug("0");
                        // Make sure we haven't taken away from the
                        // coverage that already existed or missed
                        // regions we were supposed to rescue.
                        // FIXME: skipped
                        /* CuAssertTrue(testCase, alreadyCovered[i] == 0 && coverageArray[i] == 0); */
                    }
                }
                segment = stPinchSegment_get3Prime(segment);
            }
            st_logDebug("\n");
        }

        stHash_destruct(coveragesToRescue);
        stHash_destruct(regionsAlreadyCovered);
        stPinchThreadSet_destruct(threadSet);
        free(bedRegionArray);
    }
}

CuSuite *rescueTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_rescueRandomSequences);
    return suite;
}
