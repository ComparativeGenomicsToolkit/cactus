/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stGiantComponent.h"
#include <math.h>

static stList *nodes = NULL;
static stList *edges;
static int64_t maxComponentSize;

static void teardown() {
    if (nodes != NULL) {
        stList_destruct(nodes);
        stList_destruct(edges);
        nodes = NULL;
    }
}

static void setup() {
    teardown();

    //Make nodes
    nodes = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    int64_t nodeNumber = st_randomInt(0, 1000);
    for (int64_t i = 0; i < nodeNumber; i++) {
        stList_append(nodes, stIntTuple_construct1( i));
    }

    //Make edges
    edges = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    float edgeProb = st_random();
    for (int64_t i = 0; i < nodeNumber; i++) {
        for (int64_t j = i; j < nodeNumber; j++) {
            if (st_random() <= edgeProb) {
                stList_append(edges, stIntTuple_construct3( st_randomInt(1, 100), i, j));
            }
        }
    }

    //Max component size
    maxComponentSize = 1 + log(nodeNumber) * 10; //(st_randomInt(0, nodeNumber+1);
}

static stHash *getComponents(stList *filteredEdges) {
    /*
     * A kind of stupid reimplementation of the greedy function, done just to trap typos.
     */
    stHash *nodesToComponents = stHash_construct3((uint64_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn, NULL, NULL);
    for (int64_t i = 0; i < stList_length(nodes); i++) {
        stIntTuple *node = stList_get(nodes, i);
        stSortedSet *component = stSortedSet_construct();
        stSortedSet_insert(component, node);
        stHash_insert(nodesToComponents, node, component);
    }
    for (int64_t i = 0; i < stList_length(filteredEdges); i++) {
        stIntTuple *edge = stList_get(filteredEdges, i);
        stIntTuple *node1 = stIntTuple_construct1( stIntTuple_get(edge, 1));
        stIntTuple *node2 = stIntTuple_construct1( stIntTuple_get(edge, 2));
        stSortedSet *component1 = stHash_search(nodesToComponents, node1);
        stSortedSet *component2 = stHash_search(nodesToComponents, node2);
        assert(component1 != NULL && component2 != NULL);
        if (component1 != component2) {
            stSortedSet *component3 = stSortedSet_getUnion(component1, component2);
            stSortedSetIterator *setIt = stSortedSet_getIterator(component3);
            stIntTuple *node3;
            while ((node3 = stSortedSet_getNext(setIt)) != NULL) {
                stHash_insert(nodesToComponents, node3, component3);
            }
            stSortedSet_destructIterator(setIt);
            stSortedSet_destruct(component1);
            stSortedSet_destruct(component2);
        }
        stIntTuple_destruct(node1);
        stIntTuple_destruct(node2);
    }
    return nodesToComponents;
}

static void checkComponents(CuTest *testCase, stList *filteredEdges) {
    stHash *nodesToComponents = getComponents(filteredEdges);
    //Check all components are smaller than threshold
    stList *components = stHash_getValues(nodesToComponents);
    for (int64_t i = 0; i < stList_length(components); i++) {
        stSortedSet *component = stList_get(components, i);
        CuAssertTrue(testCase, stSortedSet_size(component) <= maxComponentSize);
        CuAssertTrue(testCase, stSortedSet_size(component) >= 1);
    }
    //Check no edges can be added from those filtered.
    stSortedSet *filteredEdgesSet = stList_getSortedSet(filteredEdges, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    for (int64_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        if (stSortedSet_search(filteredEdgesSet, edge) == NULL) {
            stIntTuple *node1 = stIntTuple_construct1( stIntTuple_get(edge, 1));
            stIntTuple *node2 = stIntTuple_construct1( stIntTuple_get(edge, 2));
            stSortedSet *component1 = stHash_search(nodesToComponents, node1);
            stSortedSet *component2 = stHash_search(nodesToComponents, node2);
            CuAssertTrue(testCase, component1 != NULL && component2 != NULL);
            CuAssertTrue(testCase, component1 != component2);
            CuAssertTrue(testCase, stSortedSet_size(component1) + stSortedSet_size(component2) > maxComponentSize);
            stIntTuple_destruct(node1);
            stIntTuple_destruct(node2);
        }
    }
    stSortedSet_destruct(filteredEdgesSet);
    //Cleanup the components
    stSortedSet *componentsSet = stList_getSortedSet(components, NULL);
    stList_destruct(components);
    stSortedSet_setDestructor(componentsSet, (void(*)(void *)) stSortedSet_destruct);
    stSortedSet_destruct(componentsSet);
    stHash_destruct(nodesToComponents);
}

static void testBreakUpComponentGreedily(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        st_logInfo("Starting break up giant components random test %" PRIi64 "\n", test);
        setup();
        stList *edgesToDelete = stCaf_breakupComponentGreedily(nodes, edges, maxComponentSize);
        stSortedSet *edgesSet = stList_getSortedSet(edges, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
        stSortedSet *edgesToDeleteSet = stList_getSortedSet(edgesToDelete, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
        stSortedSet *filteredEdgesSet = stSortedSet_getDifference(edgesSet, edgesToDeleteSet);
        stList *filteredEdges = stSortedSet_getList(filteredEdgesSet);
        assert(stSortedSet_size(edgesToDeleteSet) + stSortedSet_size(filteredEdgesSet) == stSortedSet_size(edgesSet));
        checkComponents(testCase, filteredEdges);
        stSortedSet_destruct(edgesSet);
        stSortedSet_destruct(edgesToDeleteSet);
        stSortedSet_destruct(filteredEdgesSet);
        stList_destruct(filteredEdges);
        stList_destruct(edgesToDelete);
        teardown();
    }
}

static int64_t getSizeOfLargestAdjacencyComponent(stList *adjacencyComponents) {
    int64_t largestAdjacencyComponentSizeInGraph = 0;
    for (int64_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stList *adjacencyComponent = stList_get(adjacencyComponents, i);
        if (stList_length(adjacencyComponent) > largestAdjacencyComponentSizeInGraph) {
            largestAdjacencyComponentSizeInGraph = stList_length(adjacencyComponent);
        }
    }
    return largestAdjacencyComponentSizeInGraph;
}

static void testBreakUpPinchGraphAdjacencyComponentsGreedily(CuTest *testCase) {
    //return;
    for (int64_t test = 0; test < 10000; test++) {
        st_logInfo("Starting break up giant pinch graph components random test %" PRIi64 "\n", test);
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomGraph();
        int64_t totalNodes = 2 * stPinchThreadSet_getTotalBlockNumber(threadSet);
        float maximumAdjacencyComponentSizeRatio = st_random() * 10;
        int64_t maximumAdjacencyComponentSize = log(maximumAdjacencyComponentSizeRatio) * totalNodes;
        if (maximumAdjacencyComponentSize < 2) {
            maximumAdjacencyComponentSize = 2;
        }
        stPinchThreadSet_attachEnds(threadSet);
        stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents(threadSet);
        int64_t largestAdjacencyComponentSizeInGraph = getSizeOfLargestAdjacencyComponent(adjacencyComponents);
        st_logInfo(
                "We have a random pinch graph with %" PRIi64 " nodes and %" PRIi64 " adjacency components, the largest adjacency component has %" PRIi64 " nodes, with a ratio of %f we will break up adjacency components larger than %" PRIi64 " in size, this will result in a breakup: %" PRIi64 "\n",
                totalNodes, stList_length(adjacencyComponents), largestAdjacencyComponentSizeInGraph, maximumAdjacencyComponentSizeRatio,
                maximumAdjacencyComponentSize, largestAdjacencyComponentSizeInGraph > maximumAdjacencyComponentSize);
        stList_destruct(adjacencyComponents);
        stPinchThreadSet_detachEnds(threadSet);
        //Now do the actual breaking up
        stCaf_breakupComponentsGreedily(threadSet, maximumAdjacencyComponentSizeRatio);
        stPinchThreadSet_attachEnds(threadSet);
        adjacencyComponents = stPinchThreadSet_getAdjacencyComponents(threadSet);
        int64_t largestAdjacencyComponentSizeInGraphAfterBreakup = getSizeOfLargestAdjacencyComponent(adjacencyComponents);
        totalNodes = 2 * stPinchThreadSet_getTotalBlockNumber(threadSet);
        st_logInfo(
                "After splitting we have a pinch graph with %" PRIi64 " nodes and %" PRIi64 " adjacency components, the largest adjacency component has %" PRIi64 " nodes, with a ratio of %f that broke up adjacency components larger than %" PRIi64 " in size\n",
                totalNodes, stList_length(adjacencyComponents), largestAdjacencyComponentSizeInGraphAfterBreakup,
                maximumAdjacencyComponentSizeRatio, maximumAdjacencyComponentSize);
        //Cleanup
        stList_destruct(adjacencyComponents);
        stPinchThreadSet_detachEnds(threadSet);
        stPinchThreadSet_destruct(threadSet);
    }
}

typedef struct _pinch {
    int64_t thread1, thread2, start1, start2, length;
    bool strand;
} Pinch;

static void shuffle(int64_t *a, int64_t n) {
    for (int64_t i = n - 1; i > 0; i--) {
        int64_t j = st_randomInt(0, i + 1);
        int64_t t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
}

static stPinchThreadSet *buildGraph(Pinch *pinches, int64_t *pinchOrder, int64_t pinchNumber, int64_t *threadOrder,
                                    int64_t threadNumber, int64_t threadLength) {
    stPinchThreadSet *threadSet = stPinchThreadSet_construct();
    for (int64_t i = 0; i < threadNumber; i++) {
        stPinchThreadSet_addThread(threadSet, threadOrder[i], 0, threadLength);
    }
    for (int64_t i = 0; i < pinchNumber; i++) {
        Pinch *p = &pinches[pinchOrder[i]];
        stPinchThread_pinch(stPinchThreadSet_getThread(threadSet, p->thread1), stPinchThreadSet_getThread(threadSet, p->thread2),
                            p->start1, p->start2, p->length, p->strand);
    }
    // as annealing leaves it: where the pinching left a boundary between segments depends on the order of the pinches,
    // but the joined graph does not
    stPinchThreadSet_joinTrivialBoundaries(threadSet);
    return threadSet;
}

/*
 * The graph as a list of (thread, start, length, key thread, key start, orientation) for every segment, in thread and
 * then sequence order, where the key is the lowest (thread, start) of the segment's block (-1 if it has none) and the
 * orientation is relative to that segment: what the graph is, leaving out the order its blocks were made in
 */
static stList *getCanonicalGraph(stPinchThreadSet *threadSet, int64_t threadNumber) {
    stList *segments = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for (int64_t t = 0; t < threadNumber; t++) {
        stPinchThread *thread = stPinchThreadSet_getThread(threadSet, t);
        for (stPinchSegment *segment = stPinchThread_getFirst(thread); segment != NULL; segment = stPinchSegment_get3Prime(segment)) {
            int64_t keyThread = -1, keyStart = -1, orientation = 0;
            stPinchBlock *block = stPinchSegment_getBlock(segment);
            if (block != NULL) {
                stPinchSegment *keySegment = NULL;
                stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
                stPinchSegment *segment2;
                while ((segment2 = stPinchBlockIt_getNext(&blockIt)) != NULL) {
                    int64_t name = stPinchThread_getName(stPinchSegment_getThread(segment2));
                    if (keySegment == NULL || name < keyThread || (name == keyThread && stPinchSegment_getStart(segment2) < keyStart)) {
                        keySegment = segment2;
                        keyThread = name;
                        keyStart = stPinchSegment_getStart(segment2);
                    }
                }
                orientation = stPinchSegment_getBlockOrientation(segment) == stPinchSegment_getBlockOrientation(keySegment);
            }
            int64_t values[6] = { t, stPinchSegment_getStart(segment), stPinchSegment_getLength(segment), keyThread, keyStart, orientation };
            stList_append(segments, stIntTuple_constructN(6, values));
        }
    }
    return segments;
}

static bool sameGraph(stList *segments1, stList *segments2) {
    if (stList_length(segments1) != stList_length(segments2)) {
        return 0;
    }
    for (int64_t i = 0; i < stList_length(segments1); i++) {
        if (stIntTuple_cmpFn(stList_get(segments1, i), stList_get(segments2, i)) != 0) {
            return 0;
        }
    }
    return 1;
}

static void testBreakUpIsIndependentOfPinchOrder(CuTest *testCase) {
    /*
     * The same alignments pinched in a different order make the same graph, with its blocks made in a different order.
     * The breakup has to break it up the same way: cactus is deterministic only as long as the order the alignments
     * come in cannot change the result.
     *
     * Random pinches between a few short threads often align a base to another both ways round, through different
     * pinches, and then which way the graph keeps is the first pinch's: such graphs differ before any breakup and are
     * skipped.  The rest (about one in six) are compared.
     */
    int64_t tested = 0, brokenUp = 0;
    for (int64_t test = 0; test < 1000; test++) {
        int64_t threadNumber = st_randomInt(2, 6);
        int64_t threadLength = st_randomInt(50, 500);
        int64_t pinchNumber = st_randomInt(0, 300);
        Pinch *pinches = st_malloc(sizeof(Pinch) * (pinchNumber > 0 ? pinchNumber : 1));
        int64_t *pinchOrder1 = st_malloc(sizeof(int64_t) * (pinchNumber > 0 ? pinchNumber : 1));
        int64_t *pinchOrder2 = st_malloc(sizeof(int64_t) * (pinchNumber > 0 ? pinchNumber : 1));
        int64_t threadOrder1[6], threadOrder2[6];
        for (int64_t i = 0; i < pinchNumber; i++) {
            Pinch *p = &pinches[i];
            p->length = st_randomInt(1, 20);
            p->thread1 = st_randomInt(0, threadNumber);
            do {
                p->thread2 = st_randomInt(0, threadNumber);
            } while (p->thread2 == p->thread1);
            p->start1 = st_randomInt(0, threadLength - p->length + 1);
            p->start2 = st_randomInt(0, threadLength - p->length + 1);
            p->strand = st_random() > 0.5;
            pinchOrder1[i] = i;
            pinchOrder2[i] = i;
        }
        for (int64_t i = 0; i < threadNumber; i++) {
            threadOrder1[i] = i;
            threadOrder2[i] = i;
        }
        shuffle(pinchOrder2, pinchNumber);
        shuffle(threadOrder2, threadNumber);
        stPinchThreadSet *threadSet1 = buildGraph(pinches, pinchOrder1, pinchNumber, threadOrder1, threadNumber, threadLength);
        stPinchThreadSet *threadSet2 = buildGraph(pinches, pinchOrder2, pinchNumber, threadOrder2, threadNumber, threadLength);

        stList *before1 = getCanonicalGraph(threadSet1, threadNumber);
        stList *before2 = getCanonicalGraph(threadSet2, threadNumber);
        if (sameGraph(before1, before2)) {
            // a small ratio, so that the floor of 10 sets the maximum size and many graphs here get broken up
            float maximumAdjacencyComponentSizeRatio = st_random();
            stCaf_breakupComponentsGreedily(threadSet1, maximumAdjacencyComponentSizeRatio);
            stCaf_breakupComponentsGreedily(threadSet2, maximumAdjacencyComponentSizeRatio);
            stList *after1 = getCanonicalGraph(threadSet1, threadNumber);
            stList *after2 = getCanonicalGraph(threadSet2, threadNumber);
            CuAssertTrue(testCase, sameGraph(after1, after2));
            tested++;
            brokenUp += !sameGraph(before1, after1);
            stList_destruct(after1);
            stList_destruct(after2);
        }

        stList_destruct(before1);
        stList_destruct(before2);
        stPinchThreadSet_destruct(threadSet1);
        stPinchThreadSet_destruct(threadSet2);
        free(pinches);
        free(pinchOrder1);
        free(pinchOrder2);
    }
    st_logInfo("%" PRIi64 " of 1000 random graphs compared, %" PRIi64 " of them broken up\n", tested, brokenUp);
    // or the test tested nothing.  For scale: the old breakup, which numbered the ends in the order their blocks were
    // made, broke up 69 of the 166 graphs one run compared differently in the two orders
    CuAssertTrue(testCase, tested >= 50);
    CuAssertTrue(testCase, brokenUp >= 20);
}

CuSuite* giantComponentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testBreakUpComponentGreedily);
    SUITE_ADD_TEST(suite, testBreakUpPinchGraphAdjacencyComponentsGreedily);
    SUITE_ADD_TEST(suite, testBreakUpIsIndependentOfPinchOrder);
    return suite;
}
