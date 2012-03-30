/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stMatchingAlgorithms.h"
#include "stCycleConstrainedMatchingAlgorithms.h"
#include "shared.h"

/*
 * Declarations of support functions defined elsewhere reference/src
 */

stList *getComponents2(stList *adjacencyEdges, stList *stubEdges,
        stList *chainEdges);

stList *mergeSimpleCycles(stList *cycles, stList *nonZeroWeightAdjacencyEdges, stSortedSet *allAdjacencyEdges);

void checkInputs(stSortedSet *nodes, stList *adjacencyEdges,
        stList *stubEdges, stList *chainEdges);

stList *getStubAndChainEdgeFreeComponents(stList *listOfLists, stList *list1,
        stList *list2);

/*
 * The edges and nodes for the test.
 */

static int32_t nodeNumber = 0;
static stList *adjacencyEdges = NULL;
static stList *stubEdges = NULL;
static stList *chainEdges = NULL;
static stSortedSet *nodes = NULL;

static void teardown() {
    //Gets rid of the random flower.
    if (adjacencyEdges != NULL) {
        stList_destruct(adjacencyEdges);
        stList_destruct(stubEdges);
        stList_destruct(chainEdges);
        stSortedSet_destruct(nodes);
        adjacencyEdges = NULL;
        stubEdges = NULL;
        chainEdges = NULL;
        nodeNumber = 0;
        nodes = NULL;
    }
}

/*
 * Functions for fiddling with the random graphs generated.
 */

static int32_t getRandomNodeWithReplacement(int32_t nodeNumber) {
    return st_randomInt(0, nodeNumber);
}

static int32_t getRandomNodeWithoutReplacement(stSortedSet *nodes,
        int32_t nodeNumber) {
    if (stSortedSet_size(nodes) == nodeNumber) {
        assert(0);
    }
    int32_t node = getRandomNodeWithReplacement(nodeNumber);
    while (nodeInSet(nodes, node)) {
        node = getRandomNodeWithReplacement(nodeNumber);
    }
    addNodeToSet(nodes, node);
    return node;
}

static stSortedSet *getEdgeSet(stList *edges) {
    return stList_getSortedSet(edges,
            (int(*)(const void *, const void *)) stIntTuple_cmpFn);
}

static void setup(int32_t maxNodeNumber) {
    /*
     * Creates a random graph to play with.
     */
    teardown();
    do {
        nodeNumber = st_randomInt(0, maxNodeNumber) + 2;
    } while (nodeNumber % 2 != 0);
    nodes = getEmptyNodeOrEdgeSetWithCleanup();
    for(int32_t i=0; i<nodeNumber; i++) {
        addNodeToSet(nodes, i);
    }
    stubEdges = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    chainEdges = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    //Add first stub edge
    stSortedSet *nodeSet = getEmptyNodeOrEdgeSetWithCleanup();
    addEdgeToList(getRandomNodeWithoutReplacement(nodeSet, nodeNumber),
            getRandomNodeWithoutReplacement(nodeSet, nodeNumber), stubEdges);
    //While there are nodes not paired in the set..
    while (stSortedSet_size(nodeSet) < nodeNumber) {
        addEdgeToList(getRandomNodeWithoutReplacement(nodeSet, nodeNumber),
                getRandomNodeWithoutReplacement(nodeSet, nodeNumber),
                st_random() > 0.5 ? stubEdges : chainEdges);
    }
    stSortedSet_destruct(nodeSet);
    //Now make a clique of bunch of adjacency edges..
    adjacencyEdges
            = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < nodeNumber; i++) {
        for (int32_t j = i + 1; j < nodeNumber; j++) {
            addWeightedEdgeToList(i, j, st_random() > 0.8 ? st_randomInt(0, 100) : 0,
                    adjacencyEdges);
        }
    }
    //Check the inputs.
    checkInputs(nodes, adjacencyEdges, stubEdges, chainEdges);
    st_logInfo(
            "We've created a random graph with %i nodes, %i stub edges, %i chain edges, %i adjacency edges",
            nodeNumber, stList_length(stubEdges), stList_length(chainEdges),
            stList_length(adjacencyEdges));
    logEdges(adjacencyEdges, "adjacencyEdges");
    logEdges(chainEdges, "chainEdges");
    logEdges(stubEdges, "stubEdges");
}

static void testGetComponentsSimple(CuTest *testCase) {
    /*
     * Simple test, one edge incident with each node, checks there are n/2 comopnents.
     */
    for (int32_t i = 0; i < 100; i++) {
        setup(50);
        stList *stubAndChainEdges = stList_construct();
        stList_appendAll(stubAndChainEdges, stubEdges);
        stList_appendAll(stubAndChainEdges, chainEdges);
        stSortedSet *stubAndChainEdgesSet = getEdgeSet(stubAndChainEdges);
        stList *components = getComponents(stubAndChainEdges);
        CuAssertTrue(testCase, stList_length(components) == nodeNumber / 2);
        stSortedSet *componentSet = getEdgeSet(stList_get(components, 0));
        stSortedSet_equals(stubAndChainEdgesSet, componentSet);
        //Cleanup
        stSortedSet_destruct(componentSet);
        stSortedSet_destruct(stubAndChainEdgesSet);
        stList_destruct(stubAndChainEdges);
        stList_destruct(components);
        teardown();
    }
}

static void testGetComponentsP(stList *component, stSortedSet *seenEdges,
        int32_t edgeIndex) {
    /*
     * Search the edges extending out the component, to show test if it is disjoint.
     */
    stList *stack = stList_construct();
    stList_append(stack, stList_get(component, edgeIndex));
    while (stList_length(stack) > 0) {
        stIntTuple *edge = stList_pop(stack);
        if (!stSortedSet_search(seenEdges, edge)) {
            stSortedSet_insert(seenEdges, edge);
            for (int32_t i = 0; i < stList_length(component); i++) {
                stIntTuple *edge2 = stList_get(component, i);
                if (edge2 != edge) {
                    if (stIntTuple_getPosition(edge, 0)
                            == stIntTuple_getPosition(edge2, 0)
                            || stIntTuple_getPosition(edge, 0)
                                    == stIntTuple_getPosition(edge2, 1)
                            || stIntTuple_getPosition(edge, 1)
                                    == stIntTuple_getPosition(edge2, 0)
                            || stIntTuple_getPosition(edge, 1)
                                    == stIntTuple_getPosition(edge2, 1)) {
                        stList_append(stack, stList_get(component, i));
                    }
                }
            }
        }
    }
    stList_destruct(stack);
}

static void testComponentIsNotDisjoint(CuTest *testCase, stList *component) {
    /*
     * Check that the edges form one connected component.
     */
    stSortedSet *edgeSet = getEmptyNodeOrEdgeSetWithoutCleanup();
    testGetComponentsP(component, edgeSet, 0);
    CuAssertTrue(testCase,
            stSortedSet_size(edgeSet) == stList_length(component));
    stSortedSet_destruct(edgeSet);
}

static bool itemInLists(stList *listOfLists, void *o) {
    /*
     * Returns non-zero iff o is in one if the lists.
     */
    for (int32_t j = 0; j < stList_length(listOfLists); j++) {
        stList *list = stList_get(listOfLists, j);
        if (stList_contains(list, o)) {
            return 1;
        }
    }
    return 0;
}

stList *getSubsetOfEdges(stList *originalEdges) {
    stList *edges = stList_construct();
    for (int32_t i = 0; i < stList_length(originalEdges); i++) {
        if (st_random() > 0.5) {
            stList_append(edges, stList_get(originalEdges, i));
        }
    }
    return edges;
}

static void testGetComponents(CuTest *testCase) {
    /*
     * Gets a random graph (the adjacency edges), gets the components,
     * checks (0) every edge is in one component,
     * (1) each one is disjoint, (2) each one is a single component.
     */
    for (int32_t i = 0; i < 100; i++) {
        setup(50);
        stList *allEdges = getSubsetOfEdges(adjacencyEdges);
        stList_appendAll(allEdges, chainEdges);
        stList_appendAll(allEdges, stubEdges);

        stList *components = getComponents(allEdges);

        //Check every edge is in one component
        for (int32_t i = 0; i < stList_length(allEdges); i++) {
            CuAssertTrue(testCase,
                    itemInLists(components, stList_get(allEdges, i)));
        }
        st_logDebug("We got %i components\n", stList_length(components));

        //Check the components are disjoint (share no nodes).
        stSortedSet *allNodes = getEmptyNodeOrEdgeSetWithoutCleanup();
        stList *nodeSets = stList_construct3(0,
                (void(*)(void *)) stSortedSet_destruct);
        for (int32_t j = 0; j < stList_length(components); j++) {
            stSortedSet *nodesSet =
                    getNodeSetOfEdges(stList_get(components, j));
            stSortedSet *intersection = stSortedSet_getIntersection(nodesSet,
                    allNodes);
            CuAssertTrue(testCase, stSortedSet_size(intersection) == 0);
            stSortedSet_destruct(intersection);
            //Update all nodes.
            stSortedSet *allNodes2 = stSortedSet_getUnion(nodesSet, allNodes);
            stSortedSet_destruct(allNodes);
            stList_append(nodeSets, nodesSet); //clean up at the end..
            allNodes = allNodes2;
        }
        CuAssertIntEquals(testCase, nodeNumber, stSortedSet_size(allNodes));
        stSortedSet_destruct(allNodes);
        stList_destruct(nodeSets);

        //Check the components are not internally disjoint.
        for (int32_t i = 0; i < stList_length(components); i++) {
            stList *component = stList_get(components, i);
            testComponentIsNotDisjoint(testCase, component);
        }

        stList_destruct(components);
        stList_destruct(allEdges);
        teardown();
    }
}

static void testMergeSimpleCycles(CuTest *testCase) {
    /*
     * Gets a random set of cycles and checks that the output
     * is just one component, containing all the nodes.
     */
    for (int32_t i = 0; i < 100; i++) {
        setup(50);
        //Build a set of simple cycles.
        stList *allEdges = stList_construct();
        stList_appendAll(allEdges, stubEdges);
        stList_appendAll(allEdges, chainEdges);
        stList *greedyMatching = chooseMatching_greedy(adjacencyEdges,
                nodeNumber);
        assert(stList_length(greedyMatching) == nodeNumber / 2);
        stList_appendAll(allEdges, greedyMatching);
        stList_destruct(greedyMatching);

        stList *simpleCycles = getComponents(allEdges);
        stList *adjacencyOnlySimpleCycles = getStubAndChainEdgeFreeComponents(
                simpleCycles, stubEdges, chainEdges);
        for (int32_t i = 0; i < stList_length(adjacencyOnlySimpleCycles); i++) {
            assert(stList_length(stList_get(simpleCycles, i)) > 0);
            assert(stList_length(stList_get(adjacencyOnlySimpleCycles, i)) > 0);
        }

        //Get merged simple cycles
        stList *nonZeroWeightAdjacencyEdges = getEdgesWithGreaterThanZeroWeight(adjacencyEdges);
        stSortedSet *allAdjacencyEdges = stList_getSortedSet(adjacencyEdges, (int (*)(const void *, const void *))stIntTuple_cmpFn);
        stList *simpleCycle = mergeSimpleCycles(adjacencyOnlySimpleCycles,
                nonZeroWeightAdjacencyEdges, allAdjacencyEdges);
        stList_appendAll(simpleCycle, stubEdges);
        stList_appendAll(simpleCycle, chainEdges);
        testComponentIsNotDisjoint(testCase, simpleCycle);
        stSortedSet *nodeSet = getNodeSetOfEdges(simpleCycle);
        CuAssertTrue(testCase, stSortedSet_size(nodeSet) == nodeNumber);

        //Cleanup
        stList_destruct(adjacencyOnlySimpleCycles);
        stList_destruct(simpleCycles);
        stList_destruct(simpleCycle);
        stSortedSet_destruct(nodeSet);
        stList_destruct(allEdges);
        stList_destruct(nonZeroWeightAdjacencyEdges);
        stSortedSet_destruct(allAdjacencyEdges);

        teardown();
    }
}

static void checkMatching(CuTest *testCase, stList *chosenEdges,
        bool makeStubsDisjoint) {
    /*
     * Check every node has one adjacency.
     */
    stSortedSet *nodeSet = getEmptyNodeOrEdgeSetWithCleanup();
    for (int32_t i = 0; i < stList_length(chosenEdges); i++) {
        stIntTuple *edge = stList_get(chosenEdges, i);
        CuAssertTrue(testCase,
                !nodeInSet(nodeSet, stIntTuple_getPosition(edge, 0)));
        addNodeToSet(nodeSet, stIntTuple_getPosition(edge, 0));
        CuAssertTrue(testCase,
                !nodeInSet(nodeSet, stIntTuple_getPosition(edge, 1)));
        addNodeToSet(nodeSet, stIntTuple_getPosition(edge, 1));
    }
    CuAssertIntEquals(testCase, nodeNumber, stSortedSet_size(nodeSet));
    stSortedSet_destruct(nodeSet);

    /*
     * Now check the cycles all contain at least one stub, and if make stubs disjoint is true,
     * contain only one stub.
     */
    stList *components = getComponents2(chosenEdges, stubEdges, chainEdges);
    stSortedSet *stubEdgeSet = getEdgeSet(stubEdges);
    for (int32_t i = 0; i < stList_length(components); i++) {
        stList *component = stList_get(components, i);
        stList *stubComponent = stList_filterToInclude(component, stubEdgeSet);
        CuAssertTrue(testCase, stList_length(stubComponent) > 0);
        if (makeStubsDisjoint) {
            CuAssertTrue(testCase, stList_length(stubComponent) == 1);
        }
        stList_destruct(stubComponent);
    }
    stSortedSet_destruct(stubEdgeSet);
    stList_destruct(components);
}

static void testEmptyCase(CuTest *testCase) {
    /*
     * First test the empty case..
     */
    stList *list = stList_construct();
    stSortedSet *set = stSortedSet_construct();
    getMatchingWithCyclicConstraints(set,
                    list, list, list, 1,
                    chooseMatching_greedy);
    stList_destruct(list);
}

static void testGetMatchingWithCyclicConstraints(CuTest *testCase,
        bool makeStubsDisjoint,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber),
        bool isMaxCardinality, int32_t maxNodeNumber, int32_t testNumber) {
    /*
     * Creates random graphs, constructs matchings and checks they are valid.
     */
    for (int32_t i = 0; i < 100; i++) {
        setup(maxNodeNumber);
        stList *matching = matchingAlgorithm(adjacencyEdges, nodeNumber);
        stList *greedyMatching = chooseMatching_greedy(adjacencyEdges,
                nodeNumber);
        stList *cyclicMatching = getMatchingWithCyclicConstraints(nodes,
                adjacencyEdges, stubEdges, chainEdges, makeStubsDisjoint,
                matchingAlgorithm); //Do this last, as may add to adjacencies.
        checkMatching(testCase, cyclicMatching, makeStubsDisjoint);

        int32_t totalCyclicMatchingWeight = matchingWeight(cyclicMatching);
        int32_t totalMatchingWeight = matchingWeight(matching);
        int32_t totalGreedyWeight = matchingWeight(greedyMatching);
        int32_t totalCyclicMatchingCardinality = matchingCardinality(
                cyclicMatching);
        int32_t totalMatchingCardinality = matchingCardinality(matching);
        int32_t totalGreedyCardinality = matchingCardinality(greedyMatching);

        st_logInfo(
                "The total weight of the cyclic matching is %i, the total weight of the standard matching is %i, the total greedy weight is %i\n",
                totalCyclicMatchingWeight, totalMatchingWeight,
                totalGreedyWeight);
        st_logInfo(
                "The total cardinality of the cyclic matching is %i, the total Cardinality of the standard matching is %i, the total greedy Cardinality is %i\n",
                totalCyclicMatchingCardinality, totalMatchingCardinality,
                totalGreedyCardinality);

        if (isMaxCardinality) {
            CuAssertTrue(
                    testCase,
                    totalCyclicMatchingCardinality + 2 * (stList_length(
                            stubEdges) - 1) + 2 * stList_length(chainEdges)
                            >= totalMatchingCardinality);
        }

        stList_destruct(cyclicMatching);
        stList_destruct(matching);
        stList_destruct(greedyMatching);
        teardown();
    }
}

static void testGetMatchingWithCyclicConstraints_MaximumWeight_DontSweatJoinedStubs(
        CuTest *testCase) {
    st_logInfo(
            "Running matching test without disjoint stubs and maximum weight matching");
    testGetMatchingWithCyclicConstraints(testCase, 0, chooseMatching_blossom5,
            0, 500, 10);
}

static void testGetMatchingWithCyclicConstraints_MaximumWeight_MakeStubsDisjoint(
        CuTest *testCase) {
    st_logInfo(
            "Running matching test with disjoint stubs and maximum weight matching");
    testGetMatchingWithCyclicConstraints(testCase, 1, chooseMatching_blossom5,
            0, 500, 10);
}

static void testGetMatchingWithCyclicConstraints_MaximumCardinality_DontSweatJoinedStubs(
        CuTest *testCase) {
    st_logInfo(
            "Running matching test without disjoint stubs and maximum weight matching");
    testGetMatchingWithCyclicConstraints(testCase, 0,
            chooseMatching_maximumCardinalityMatching, 1, 50, 50);
}

static void testGetMatchingWithCyclicConstraints_MaximumCardinality_MakeStubsDisjoint(
        CuTest *testCase) {
    st_logInfo(
            "Running matching test with disjoint stubs and maximum weight matching");
    testGetMatchingWithCyclicConstraints(testCase, 1,
            chooseMatching_maximumCardinalityMatching, 1, 50, 50);
}

CuSuite* cyclesConstrainedMatchingAlgorithmsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testEmptyCase);
    SUITE_ADD_TEST(suite, testGetComponentsSimple);
    SUITE_ADD_TEST(suite, testGetComponents);
    SUITE_ADD_TEST(suite, testMergeSimpleCycles);
    SUITE_ADD_TEST(suite,
            testGetMatchingWithCyclicConstraints_MaximumWeight_DontSweatJoinedStubs);
    SUITE_ADD_TEST(suite,
            testGetMatchingWithCyclicConstraints_MaximumWeight_MakeStubsDisjoint);
    SUITE_ADD_TEST(suite,
            testGetMatchingWithCyclicConstraints_MaximumCardinality_DontSweatJoinedStubs);
    SUITE_ADD_TEST(suite,
            testGetMatchingWithCyclicConstraints_MaximumCardinality_MakeStubsDisjoint);
    return suite;
}
