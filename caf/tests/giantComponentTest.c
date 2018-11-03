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
        stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents(threadSet);
        int64_t largestAdjacencyComponentSizeInGraph = getSizeOfLargestAdjacencyComponent(adjacencyComponents);
        st_logInfo(
                "We have a random pinch graph with %" PRIi64 " nodes and %" PRIi64 " adjacency components, the largest adjacency component has %" PRIi64 " nodes, with a ratio of %f we will break up adjacency components larger than %" PRIi64 " in size, this will result in a breakup: %" PRIi64 "\n",
                totalNodes, stList_length(adjacencyComponents), largestAdjacencyComponentSizeInGraph, maximumAdjacencyComponentSizeRatio,
                maximumAdjacencyComponentSize, largestAdjacencyComponentSizeInGraph > maximumAdjacencyComponentSize);
        stList_destruct(adjacencyComponents);
        //Now do the actual breaking up
        stCaf_breakupComponentsGreedily(threadSet, maximumAdjacencyComponentSizeRatio);
        adjacencyComponents = stPinchThreadSet_getAdjacencyComponents(threadSet);
        int64_t largestAdjacencyComponentSizeInGraphAfterBreakup = getSizeOfLargestAdjacencyComponent(adjacencyComponents);
        totalNodes = 2 * stPinchThreadSet_getTotalBlockNumber(threadSet);
        st_logInfo(
                "After splitting we have a pinch graph with %" PRIi64 " nodes and %" PRIi64 " adjacency components, the largest adjacency component has %" PRIi64 " nodes, with a ratio of %f that broke up adjacency components larger than %" PRIi64 " in size\n",
                totalNodes, stList_length(adjacencyComponents), largestAdjacencyComponentSizeInGraphAfterBreakup,
                maximumAdjacencyComponentSizeRatio, maximumAdjacencyComponentSize);
        //Cleanup
        stList_destruct(adjacencyComponents);
        stPinchThreadSet_destruct(threadSet);
    }
}

CuSuite* giantComponentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testBreakUpComponentGreedily);
    SUITE_ADD_TEST(suite, testBreakUpPinchGraphAdjacencyComponentsGreedily);
    return suite;
}
