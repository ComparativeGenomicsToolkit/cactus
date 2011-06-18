/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "matchingAlgorithms.h"
#include "cycleConstrainedMatchingAlgorithms.h"

/*
 * Declarations of support functions defined elsewhere reference/src
 */

stList *getComponents2(stList *adjacencyEdges, stList *stubEdges,
        stList *chainEdges);

stList *filterToInclude(stList *list, stSortedSet *set);

stList *getComponents(stList *edges);

stList *mergeCycles(stList *components, stList *adjacencyEdges);

int32_t matchingWeight(stList *matching);

stList *splitMultipleStubCycles(stList *chosenEdges, stList *adjacencyEdges, stList *stubEdges,
        stList *chainEdges);

/*
 * The edges and nodes for the test.
 */

static int32_t nodeNumber = 0;
static stList *adjacencyEdges = NULL;
static stList *stubEdges = NULL;
static stList *chainEdges = NULL;

static void teardown() {
    //Gets rid of the random flower.
    if(adjacencyEdges != NULL) {
        stList_destruct(adjacencyEdges);
        stList_destruct(stubEdges);
        stList_destruct(chainEdges);
        adjacencyEdges = NULL;
        stubEdges = NULL;
        chainEdges = NULL;
        nodeNumber = 0;
    }
}

/*
 * Functions for fiddling with the random graphs generated.
 */

static void addWeightedEdge(int32_t node1, int32_t node2, int32_t weight, stList *edges) {
    stList_append(edges, stIntTuple_construct(3, node1, node2, weight));
}

static void addEdge(int32_t node1, int32_t node2, stList *edges) {
    stList_append(edges, stIntTuple_construct(2, node1, node2));
}

static bool nodeInSet(stSortedSet *nodes, int32_t node) {
    stIntTuple *i = stIntTuple_construct(1, node);
    bool b = stSortedSet_search(nodes, i) != NULL;
    stIntTuple_destruct(i);
    return b;
}

static void addNodeToSet(stSortedSet *nodes, int32_t node) {
    assert(!nodeInSet(nodes, node));
    stIntTuple *i = stIntTuple_construct(1, node);
    stSortedSet_insert(nodes, i);
}

static stSortedSet *getEmptyNodeSet() {
    return stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, (void (*)(void *))stIntTuple_destruct);
}

static int32_t getRandomNodeWithReplacement(int32_t nodeNumber) {
    return st_randomInt(0, nodeNumber);
}

static int32_t getRandomNodeWithoutReplacement(stSortedSet *nodes, int32_t nodeNumber) {
    if(stSortedSet_size(nodes) == nodeNumber) {
        assert(0);
    }
    int32_t node = getRandomNodeWithReplacement(nodeNumber);
    while(nodeInSet(nodes, node)) {
        node = getRandomNodeWithReplacement(nodeNumber);
    }
    addNodeToSet(nodes, node);
    return node;
}

static void logEdges(stList *edges, const char *edgesName) {
    for(int32_t i=0; i<stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        st_logDebug("Edge, type: %s, node1: %i, node2: %i\n", edgesName, stIntTuple_getPosition(edge, 0), stIntTuple_getPosition(edge, 1), (stIntTuple_length(edge) == 3 ? stIntTuple_getPosition(edge, 2) : INT32_MAX));
    }
}


static void setup() {
    /*
     * Creates a random graph to play with.
     */
    teardown();
    do {
        nodeNumber = st_randomInt(0, 100) + 2;
    } while(nodeNumber % 2 != 0);
    adjacencyEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    stubEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    chainEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    //Add first stub edge
    stSortedSet *nodeSet = getEmptyNodeSet();
    addEdge(getRandomNodeWithoutReplacement(nodeSet, nodeNumber), getRandomNodeWithoutReplacement(nodeSet, nodeNumber), stubEdges);
    //While there are nodes not paired in the set..
    while(stSortedSet_size(nodeSet) < nodeNumber) {
        addEdge(getRandomNodeWithoutReplacement(nodeSet, nodeNumber), getRandomNodeWithoutReplacement(nodeSet, nodeNumber), st_random() > 0.5 ? stubEdges : chainEdges);
    }
    stSortedSet_destruct(nodeSet);
    //Now randomly make a bunch of weighted adjacency edges..
    int32_t edgeNumber = st_randomInt(0, nodeNumber * 10);
    for(int32_t i=0; i<edgeNumber; i++) {
        int32_t from = getRandomNodeWithReplacement(nodeNumber);
        int32_t to = getRandomNodeWithReplacement(nodeNumber);
        if(from != to) {
            addWeightedEdge(from, to, st_randomInt(0, 100), adjacencyEdges);
        }
    }
    st_logInfo("We've created a random graph with %i nodes, %i stub edges, %i chain edges, %i adjacency edges", nodeNumber, stList_length(stubEdges), stList_length(chainEdges), stList_length(adjacencyEdges));
    logEdges(adjacencyEdges, "adjacencyEdges");
    logEdges(chainEdges, "chainEdges");
    logEdges(stubEdges, "stubEdges");
}


static void checkMatching(CuTest *testCase, stList *chosenEdges, bool makeStubsDisjoint) {
    /*
     * Check every node has one adjacency.
     */
    stSortedSet *nodeSet = getEmptyNodeSet();
    for(int32_t i=0; i<stList_length(chosenEdges); i++) {
        stIntTuple *edge = stList_get(chosenEdges, i);
        CuAssertTrue(testCase,!nodeInSet(nodeSet, stIntTuple_getPosition(edge, 0)));
        addNodeToSet(nodeSet, stIntTuple_getPosition(edge, 0));
        CuAssertTrue(testCase,!nodeInSet(nodeSet, stIntTuple_getPosition(edge, 1)));
        addNodeToSet(nodeSet, stIntTuple_getPosition(edge, 1));
    }
    CuAssertTrue(testCase,stSortedSet_size(nodeSet) == nodeNumber);
    stSortedSet_destruct(nodeSet);

    /*
     * Now check the cycles all contain at least one stub, and if make stubs disjoint is true,
     * contain only one stub.
     */
    stList *components = getComponents2(adjacencyEdges, stubEdges, chainEdges);
    stSortedSet *stubEdgeSet = stList_getSortedSet(stubEdges, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    for(int32_t i=0; stList_length(components); i++) {
        stList *component = stList_get(components, i);
        stList *stubComponent = filterToInclude(component, stubEdgeSet);
        CuAssertTrue(testCase, stList_length(stubComponent) > 0);
        if(makeStubsDisjoint) {
            CuAssertTrue(testCase, stList_length(stubComponent) == 1);
        }
        stList_destruct(stubComponent);
    }
    stSortedSet_destruct(stubEdgeSet);
    stList_destruct(components);
}

static void testGetComponentsSimple(CuTest *testCase) {
    /*
     * Gets a random cycle (the stub and chain edges), checks there is just one component.
     */
}

static void testGetComponents(CuTest *testCase) {
    /*
     * Gets a random graph (the adjacency edges), gets the components,
     * checks (1) each one is disjoint, (2) each one is a single component.
     */
}

static void testMergeCycles(CuTest *testCase) {
    /*
     * Gets a random set of cycles graph (the adjacency edges) and its components, and checks that the output is just one component.
     */
}

static void testSplitMultipleStubCycles(CuTest *testCase) {
    /*
     * Gets a random graph (the adjacency edges) and its components, and then
     */
}

static void testGetMatchingWithCyclicConstraints(CuTest *testCase,
        bool makeStubsDisjoint, stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber),
        bool isMaxCardinality) {
    /*
     * Creates random graphs, constructs matchings and checks they are valid.
     */
    for(int32_t i=0; i<100; i++) {
        setup();
        stList *cyclicMatching = getMatchingWithCyclicConstraints(nodeNumber, adjacencyEdges, stubEdges, chainEdges, makeStubsDisjoint, matchingAlgorithm);
        stList *matching = matchingAlgorithm(adjacencyEdges, nodeNumber);
        stList *greedyMatching = chooseMatching_greedy(adjacencyEdges, nodeNumber);
        checkMatching(testCase, cyclicMatching, makeStubsDisjoint);

        int32_t totalCyclicMatchingWeight = matchingWeight(cyclicMatching);
        int32_t totalMatchingWeight = matchingWeight(matching);
        int32_t totalGreedyWeight = matchingWeight(greedyMatching);

        st_logInfo("The total weight of the cyclic matching is %i, the total weight of the standard matching is %i, the total greedy weight is %i\n",
                totalCyclicMatchingWeight, totalMatchingWeight, totalGreedyWeight);

        if(isMaxCardinality) {
            CuAssertTrue(testCase,totalCyclicMatchingWeight + stList_length(stubEdges) + stList_length(chainEdges) -1 >= totalMatchingWeight);
        }

        stList_destruct(cyclicMatching);
        stList_destruct(matching);
        stList_destruct(greedyMatching);
        teardown();
    }
}

static void testGetMatchingWithCyclicConstraints_MaximumWeight_DontSweatJoinedStubs(CuTest *testCase) {
    st_logInfo("Running matching test without disjoint stubs and maximum weight matching");
    testGetMatchingWithCyclicConstraints(testCase, 1, chooseMatching_maximumWeightMatching, 0);
}

static void testGetMatchingWithCyclicConstraints_MaximumWeight_MakeStubsDisjoint(CuTest *testCase) {
    st_logInfo("Running matching test with disjoint stubs and maximum weight matching");
    testGetMatchingWithCyclicConstraints(testCase, 1, chooseMatching_maximumWeightMatching, 0);
}

static void testGetMatchingWithCyclicConstraints_MaximumCardinality_DontSweatJoinedStubs(CuTest *testCase) {
    st_logInfo("Running matching test without disjoint stubs and maximum weight matching");
    testGetMatchingWithCyclicConstraints(testCase, 1, chooseMatching_maximumCardinalityMatching, 1);
}

static void testGetMatchingWithCyclicConstraints_MaximumCardinality_MakeStubsDisjoint(CuTest *testCase) {
    st_logInfo("Running matching test with disjoint stubs and maximum weight matching");
    testGetMatchingWithCyclicConstraints(testCase, 1, chooseMatching_maximumCardinalityMatching, 1);
}


CuSuite* cyclesConstrainedMatchingAlgorithmsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testGetComponents);
    SUITE_ADD_TEST(suite, testMergeCycles);
    SUITE_ADD_TEST(suite, testSplitMultipleStubCycles);
    SUITE_ADD_TEST(suite, testGetMatchingWithCyclicConstraints_MaximumWeight_DontSweatJoinedStubs);
    SUITE_ADD_TEST(suite, testGetMatchingWithCyclicConstraints_MaximumWeight_MakeStubsDisjoint);
    SUITE_ADD_TEST(suite, testGetMatchingWithCyclicConstraints_MaximumCardinality_DontSweatJoinedStubs);
    SUITE_ADD_TEST(suite, testGetMatchingWithCyclicConstraints_MaximumCardinality_MakeStubsDisjoint);
    return suite;
}
