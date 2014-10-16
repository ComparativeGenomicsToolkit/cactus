#include "CuTest.h"
#include "sonLib.h"
#include "cactus.h"
#include "stCafPhylogeny.h"
#include "stPinchGraphs.h"

// O(n).
static stPinchSegment *getSegmentByBlockIndex(stPinchBlock *block,
                                              int64_t index) {
    assert(index >= 0);
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment = NULL;
    int64_t i = 0;
    while((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        if (i == index) {
            return segment;
        }
        i++;
    }
    return NULL;
}

static void testSplitBlock(CuTest *testCase) {
    for (int64_t testNum = 0; testNum < 1000; testNum++) {
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomGraph();
        stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
        bool allowSingleDegreeBlocks = st_random() > 0.5;
        stPinchBlock *block;
        // Grab a list of all blocks in the graph (in case the
        // iterator doesn't deal well with adding/removing blocks)
        stList *blocks = stList_construct();
        while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
            stList_append(blocks, block);
        }
        for (int64_t blockNum = 0; blockNum < stList_length(blocks); blockNum++) {
            block = stList_get(blocks, blockNum);
            // Create a random partition for this block
            int64_t blockDegree = stPinchBlock_getDegree(block);
            int64_t blockLength = stPinchBlock_getLength(block);
            stList *partitions = stList_construct3(0, (void (*)(void *)) stList_destruct);
            stList *segmentss = stList_construct3(0, (void (*)(void *)) stSet_destruct);
            int64_t maxNumPartitions = st_randomInt64(1, blockDegree + 1);
            for (int64_t i = 0; i < maxNumPartitions; i++) {
                stList *partition = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
                stList_append(partitions, partition);
                stSet *segments = stSet_construct();
                stList_append(segmentss, segments);
            }
            for (int64_t i = 0; i < blockDegree; i++) {
                int64_t partitionNumber = st_randomInt64(0, maxNumPartitions);
                stList_append(stList_get(partitions, partitionNumber),
                              stIntTuple_construct1(i));
                stSet *segmentSet = stList_get(segmentss, partitionNumber);
                CuAssertTrue(testCase, stPinchSegment_getLength(getSegmentByBlockIndex(block, i)) == blockLength);
                stSet_insert(segmentSet, getSegmentByBlockIndex(block, i));
            }

            // Get rid of any partitions which happen to have 0 entries.
            for (int64_t i = 0; i < stList_length(partitions); i++) {
                stList *partition = stList_get(partitions, i);
                stSet *segments = stList_get(segmentss, i);
                if (stList_length(partition) == 0) {
                    CuAssertTrue(testCase, stSet_size(segments) == 0);
                    stList_destruct(stList_remove(partitions, i));
                    stSet_destruct(stList_remove(segmentss, i));
                    i--; // Compensate for loss of an entry in the lists.
                }
            }
            if (stList_length(partitions) == 1) {
                // Single-degree blocks can't be split so they'll be
                // kept no matter what the allowSingleDegreeBlocks
                // setting is
                allowSingleDegreeBlocks = 1;
            }

            // Split the block
            splitBlock(block, partitions, allowSingleDegreeBlocks);

            // Make sure the block partitioned correctly
            stSet *seenBlocks = stSet_construct();
            for (int64_t i = 0; i < stList_length(partitions); i++) {
                stList *partition = stList_get(partitions, i);
                stSet *segments = stList_get(segmentss, i);
                CuAssertTrue(testCase,
                             stSet_size(segments) == stList_length(partition));
                CuAssertTrue(testCase, stSet_size(segments) > 0);
                // Make sure that all the segments are in the same
                // block and that they still have the same length.
                stSetIterator *segmentIt = stSet_getIterator(segments);
                stPinchSegment *segment;
                stPinchBlock *block;
                bool firstSegment = true;
                while ((segment = stSet_getNext(segmentIt)) != NULL) {
                    CuAssertTrue(testCase, stPinchSegment_getLength(segment) == blockLength);
                    if (firstSegment) {
                        block = stPinchSegment_getBlock(segment);
                    }
                    CuAssertTrue(testCase,
                                 stPinchSegment_getBlock(segment) == block);
                    firstSegment = false;
                }
                stSet_destructIterator(segmentIt);
                // Now we've confirmed there is just one block that
                // contains all the segments in this partition, so run
                // some tests on that block.
                CuAssertTrue(testCase, stSet_search(seenBlocks, block) == NULL);
                if (block != NULL) {
                    CuAssertTrue(testCase, stPinchBlock_getLength(block) == blockLength);
                    CuAssertTrue(testCase, stPinchBlock_getDegree(block) == stSet_size(segments));
                }
                if (stList_length(partition) == 1) {
                    if (allowSingleDegreeBlocks) {
                        CuAssertTrue(testCase, block != NULL);
                        CuAssertTrue(testCase,
                                     stPinchBlock_getDegree(block) == 1);
                    } else {
                        CuAssertTrue(testCase, block == NULL);
                    }
                } else {
                    CuAssertTrue(testCase, block != NULL);
                    CuAssertTrue(testCase, stPinchBlock_getDegree(block) == stList_length(partition));
                }
                stSet_insert(seenBlocks, block);
            }

            // Clean up
            stList_destruct(partitions);
            stList_destruct(segmentss);
            stSet_destruct(seenBlocks);
        }
        stList_destruct(blocks);
        stPinchThreadSet_destruct(threadSet);
    }
}

CuSuite *phylogenyTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testSplitBlock);
    return suite;
}
