#include "CuTest.h"
#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCafPhylogeny.h"
#include "stCaf.h"
/*
static CactusDisk *cactusDisk;
static Flower *flower;
static EventTree *eventTree;
static stPinchThreadSet *threadSet;

static void setup(void) {
    cactusDisk = testCommon_getTemporaryCactusDisk();
    flower = flower_construct(cactusDisk);

    // Event tree: (((human,(mouse,rat)Anc3)Anc2,(cow,dog)Anc1)Anc0)ROOT;
    eventTree = eventTree_construct2(flower);
    event_construct3("Anc0", 0.0, eventTree_getRootEvent(eventTree), eventTree);
    event_construct3("Anc1", 1.0, eventTree_getEventByHeader(eventTree, "Anc0"), eventTree);
    event_construct3("dog", 1.0, eventTree_getEventByHeader(eventTree, "Anc1"), eventTree);
    event_construct3("cow", 1.0, eventTree_getEventByHeader(eventTree, "Anc1"), eventTree);
    event_construct3("Anc2", 1.0, eventTree_getEventByHeader(eventTree, "Anc0"), eventTree);
    event_construct3("Anc3", 1.0, eventTree_getEventByHeader(eventTree, "Anc2"), eventTree);
    event_construct3("rat", 1.0, eventTree_getEventByHeader(eventTree, "Anc3"), eventTree);
    event_construct3("mouse", 1.0, eventTree_getEventByHeader(eventTree, "Anc3"), eventTree);
    event_construct3("human", 1.0, eventTree_getEventByHeader(eventTree, "Anc2"), eventTree);

    threadSet = stCaf_setup(flower);
}

static void teardown(void) {
    testCommon_deleteTemporaryCactusDisk(cactusDisk);
}
*/

// Assume that the leaves of the gene tree are labeled according to
// their species names and produce a leafToSpecies hash.
static stHash *getLeafToSpeciesUsingLeafLabels(stTree *geneTree,
                                               stTree *speciesTree) {
    stHash *ret = stHash_construct();
    stList *bfQueue = stList_construct();
    stList_append(bfQueue, geneTree);
    while (stList_length(bfQueue) != 0) {
        stTree *gene = stList_pop(bfQueue);
        for (int64_t i = 0; i < stTree_getChildNumber(gene); i++) {
            stList_append(bfQueue, stTree_getChild(gene, i));
        }
        if (stTree_getChildNumber(gene) == 0) {
            stTree *species = stTree_findChild(speciesTree,
                                               stTree_getLabel(gene));
            assert(species != NULL);
            stHash_insert(ret, gene, species);
        }
    }
    return ret;
}

static void addArbitraryStIndexedTreeInfo_R(stTree *node, int64_t *curIndex,
                                            stHash *nodeToLabel) {
    if (stTree_getChildNumber(node) == 0) {
        stHash_insert(nodeToLabel, node, stString_copy(stTree_getLabel(node)));
        stTree_setLabel(node, stString_print_r("%" PRIi64, *curIndex));
        (*curIndex)++;
    }
    for (int64_t i = 0; i < stTree_getChildNumber(node); i++) {
        addArbitraryStIndexedTreeInfo_R(stTree_getChild(node, i), curIndex,
                                        nodeToLabel);
    }
}

// Add stIndexedTreeInfo, but without caring about which labels get
// assigned which matrix indices.
static void addArbitraryStIndexedTreeInfo(stTree *geneTree) {
    stHash *nodeToLabel = stHash_construct2(NULL, free);
    int64_t curIndex = 0;
    addArbitraryStIndexedTreeInfo_R(geneTree, &curIndex, nodeToLabel);
    stPhylogeny_addStIndexedTreeInfo(geneTree);
    stHashIterator *leafIt = stHash_getIterator(nodeToLabel);
    stTree *leaf = NULL;
    while ((leaf = stHash_getNext(leafIt)) != NULL) {
        stTree_setLabel(leaf, stHash_search(nodeToLabel, leaf));
    }
    stHash_destruct(nodeToLabel);
}

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

static void test_stCaf_splitBlock(CuTest *testCase) {
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
            stList *partitionedBlocks = stCaf_splitBlock(block, partitions, allowSingleDegreeBlocks);
            CuAssertIntEquals(testCase, stList_length(partitions), stList_length(partitionedBlocks));
            stList_destruct(partitionedBlocks);

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
                stPinchBlock *block = NULL;
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

static stList *setupTestChain(stPinchThreadSet *threadSet) {
    // Create a chain that looks like this:
    // thread 1: =1=>--=2=>=3=>--=2=>
    // thread 2: <2==--<1==<2==--<3==
    // thread 3: =3=>--=3=>=1=>--=1=>
    stPinchThread *thread1 = stPinchThreadSet_addThread(threadSet, 1, 0, 200);
    stPinchThread *thread2 = stPinchThreadSet_addThread(threadSet, 2, 0, 200);
    stPinchThread *thread3 = stPinchThreadSet_addThread(threadSet, 3, 0, 200);

    stPinchThread_pinch(thread1, thread2, 0, 90, 10, false);
    stPinchThread_pinch(thread1, thread3, 0, 0, 10, true);

    stPinchThread_pinch(thread1, thread3, 20, 20, 10, true);
    stPinchThread_pinch(thread2, thread1, 70, 20, 10, false);

    stPinchThread_pinch(thread3, thread2, 30, 60, 10, false);
    stPinchThread_pinch(thread3, thread1, 30, 30, 10, true);

    stPinchThread_pinch(thread3, thread2, 50, 40, 10, false);
    stPinchThread_pinch(thread3, thread1, 50, 50, 10, true);

    stList *chain = stList_construct();
    stList_append(chain, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 0)));
    stList_append(chain, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 20)));
    stList_append(chain, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 30)));
    stList_append(chain, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 50)));

    return chain;
}

static stList *setupTestChainWithChildChains(stPinchThreadSet *threadSet) {
    // Add some blocks in the adjacencies of the chain, so any
    // accidental assumptions about the blocks being directly adjacent
    // to each other are caught out.
    stList *chain = setupTestChain(threadSet);
    stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, 1);
    stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, 2);
    stPinchThread *thread3 = stPinchThreadSet_getThread(threadSet, 3);

    stPinchThread_pinch(thread1, thread2, 14, 84, 2, true);
    stPinchThread_pinch(thread2, thread3, 54, 44, 2, false);
    return chain;
}

/* static stList *setupTestChainWithTandemDup(stPinchThreadSet *threadSet) { */
/*     stList *chain = setupTestChain(threadSet); */
/*     stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, 1); */

/*     stPinchThread_pinch(thread1, thread1, 0, 150, 10, false); */
/*     stPinchThread_pinch(thread1, thread1, 20, 130, 10, false); */
/*     stPinchThread_pinch(thread1, thread1, 30, 120, 10, false); */
/*     stPinchThread_pinch(thread1, thread1, 50, 100, 10, false); */

/*     return chain; */
/* } */

static void test_stCaf_splitChainP(CuTest *testCase, stList *(*setup)(stPinchThreadSet *)) {
    stPinchThreadSet *threadSet = stPinchThreadSet_construct();
    stList *chain = setup(threadSet);
    int64_t oldBlockNumber = stPinchThreadSet_getTotalBlockNumber(threadSet);

    stPinchThread *thread1 = stPinchThreadSet_getThread(threadSet, 1);
    stPinchThread *thread2 = stPinchThreadSet_getThread(threadSet, 2);
    stPinchThread *thread3 = stPinchThreadSet_getThread(threadSet, 3);

    stList *partitions = stList_construct3(0, (void (*)(void *)) stList_destruct);
    stList *partition1 = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    stList_append(partition1, stIntTuple_construct1(0));
    stList_append(partition1, stIntTuple_construct1(2));
    stList_append(partitions, partition1);
    stList *partition2 = stList_construct();
    stList_append(partition2, stIntTuple_construct1(1));
    stList_append(partitions, partition2);

    stList *partitionedChains = stCaf_splitChain(chain, partitions, true);
    CuAssertIntEquals(testCase, stList_length(partitionedChains), stList_length(partitions));
    CuAssertIntEquals(testCase, stList_length(stList_get(partitionedChains, 0)), stList_length(chain));
    CuAssertIntEquals(testCase, stList_length(stList_get(partitionedChains, 1)), stList_length(chain));

    CuAssertIntEquals(testCase, oldBlockNumber + stList_length(chain), stPinchThreadSet_getTotalBlockNumber(threadSet));
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 0)) !=
                           stPinchSegment_getBlock(stPinchThread_getSegment(thread2, 90)));
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 0)) ==
                           stPinchSegment_getBlock(stPinchThread_getSegment(thread3, 0)));
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 20)) !=
                           stPinchSegment_getBlock(stPinchThread_getSegment(thread2, 70)));
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 20)) ==
                           stPinchSegment_getBlock(stPinchThread_getSegment(thread3, 20)));
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 30)) !=
                           stPinchSegment_getBlock(stPinchThread_getSegment(thread2, 60)));
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 30)) ==
                           stPinchSegment_getBlock(stPinchThread_getSegment(thread3, 30)));
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 50)) !=
                           stPinchSegment_getBlock(stPinchThread_getSegment(thread2, 40)));
    CuAssertTrue(testCase, stPinchSegment_getBlock(stPinchThread_getSegment(thread1, 50)) ==
                           stPinchSegment_getBlock(stPinchThread_getSegment(thread3, 50)));

    stPinchThreadSet_destruct(threadSet);
    stList_destruct(partitions);
    stList_destruct(chain);
}

static void test_stCaf_splitChain(CuTest *testCase) {
    // Test stCaf_splitChain on the simple test chain
    test_stCaf_splitChainP(testCase, setupTestChain);
    // Test on the same chain but with children
    test_stCaf_splitChainP(testCase, setupTestChainWithChildChains);
}

static void test_stCaf_findAndRemoveSplitBranches(CuTest *testCase) {
    stTree *speciesTree = stTree_parseNewickString("((human,(mouse,rat)Anc3)Anc2,(cow,dog)Anc1)Anc0;");
    // Here the reference event is Anc3.
    stSet *speciesToSplitOn = stSet_construct();
    stSet_insert(speciesToSplitOn, stTree_findChild(speciesTree, "Anc3"));
    stSet_insert(speciesToSplitOn, stTree_findChild(speciesTree, "Anc2"));
    stSet_insert(speciesToSplitOn, stTree_findChild(speciesTree, "Anc0"));
    stSortedSet *splitBranches = stSortedSet_construct3((int (*)(const void *, const void *)) stCaf_SplitBranch_cmp, free);

    // Create a tree that has a duplication on Anc3.
    stTree *geneTree = stTree_parseNewickString("((human,((mouse,rat)Anc3.0,(mouse,rat)Anc3.1)Anc3)Anc2,(cow,dog)Anc1)Anc0;");
    addArbitraryStIndexedTreeInfo(geneTree);
    stHash *leafToSpecies = getLeafToSpeciesUsingLeafLabels(geneTree,
                                                            speciesTree);
    stPhylogeny_reconcileAtMostBinary(geneTree, leafToSpecies, true);
    stCaf_findSplitBranches(NULL, geneTree, splitBranches, speciesToSplitOn);
    CuAssertTrue(testCase, stSortedSet_size(splitBranches) == 2);
    stCaf_removeSplitBranches(NULL, geneTree, speciesToSplitOn, splitBranches);
    CuAssertTrue(testCase, stSortedSet_size(splitBranches) == 0);

    stHash_destruct(leafToSpecies);
    stPhylogenyInfo_destructOnTree(geneTree);
    stTree_destruct(geneTree);
    stSet_destruct(speciesToSplitOn);
    stTree_destruct(speciesTree);
}

CuSuite *phylogenyTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stCaf_splitBlock);
    SUITE_ADD_TEST(suite, test_stCaf_splitChain);
    SUITE_ADD_TEST(suite, test_stCaf_findAndRemoveSplitBranches);
    return suite;
}
