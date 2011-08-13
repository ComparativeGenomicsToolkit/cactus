/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "cactusReference.h"
#include "adjacencyProblem.h"
#include "cactusMatchingAlgorithms.h"
#include "checkEdges.h"

static stList *stubs;
static stList *chains;
static double *zMatrix;
static int32_t nodeNumber;

static void teardown() {
    if (nodeNumber != -1) {
        stList_destruct(stubs);
        stList_destruct(chains);
        free(zMatrix);
        nodeNumber = -1;
    }
}

static void setup() {
    teardown();
    assert(nodeNumber == -1);
    while(nodeNumber % 2 != 0) {
        nodeNumber = st_randomInt(0, 100);
    }
    assert(nodeNumber >= 0);
    assert(nodeNumber % 2 == 0);
    stubs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    chains = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int32_t i=0; i<nodeNumber/2; i++) {
        assert(nodeNumber/2 > 0);
        stIntTuple *edge = stIntTuple_construct(2, i, nodeNumber/2 + i);
        if(stList_length(stubs) == 0 || st_random(stubs) > 0.9) {
            stList_append(stubs, edge);
        }
        else {
            stList_append(chains, edge);
        }
    }
    zMatrix = st_calloc(nodeNumber*nodeNumber, sizeof(double));
    for(int32_t i=0; i<nodeNumber; i++) {
        for(int32_t j=i+1; j<nodeNumber; j++) {
            double score = st_random();
            zMatrix[i * nodeNumber + j] = score;
            zMatrix[j * nodeNumber + i] = score;
        }
    }
    st_logDebug("To test the adjacency problem we've created a problem with %i nodes %i stubs and %i chains\n", nodeNumber, stList_length(stubs), stList_length(chains));
}

static void checkIsValidReference(CuTest *testCase, stList *reference,
        double totalScore) {
    stList *chosenEdges = convertReferenceToAdjacencyEdges(reference);
    //Check that everyone has a partner.
    CuAssertIntEquals(testCase, nodeNumber, stList_length(chosenEdges) * 2);
    stSortedSet *nodes = stSortedSet_construct3((int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < nodeNumber; i++) {
        stSortedSet_insert(nodes, stIntTuple_construct(1, i));
    }
    checkEdges(chosenEdges, nodes, 1, 0);
    //Check that the score is correct
    double totalScore2 = calculateZScoreOfReference(reference, nodeNumber, zMatrix);
    CuAssertDblEquals(testCase, totalScore2, totalScore, 0.000001);
    //Check that the stubs are properly connected.
    stList *components = getComponents(chosenEdges);
    CuAssertIntEquals(testCase, stList_length(stubs), stList_length(components));
    CuAssertIntEquals(testCase, stList_length(stubs), stList_length(reference));
    for (int32_t i = 0; i < stList_length(stubs); i++) {
        stIntTuple *stub = stList_get(stubs, i);
        int32_t node1 = stIntTuple_getPosition(stub, 0);
        int32_t node2 = stIntTuple_getPosition(stub, 1);
        bool b = 0;
        for (int32_t j = 0; j < stList_length(components); j++) {
            stList *component = stList_get(components, j);
            stSortedSet *componentSet = stList_getSortedSet(component, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
            stIntTuple *node1A = stIntTuple_construct(1, node1);
            stIntTuple *node2A = stIntTuple_construct(1, node2);
            if(stSortedSet_search(componentSet, node1A) != NULL) {
                CuAssertTrue(testCase, stSortedSet_search(componentSet, node2A) != NULL);
                CuAssertTrue(testCase, !b);
                b = 1;
                stList_removeItem(components, component);
                stSortedSet_destruct(componentSet);
                stList_destruct(component);
            }
            else {
                CuAssertTrue(testCase, stSortedSet_search(componentSet, node2A) == NULL);
                stSortedSet_destruct(componentSet);
            }
            stIntTuple_destruct(node1A);
            stIntTuple_destruct(node2A);
        }
        CuAssertTrue(testCase, b);
    }
    CuAssertTrue(testCase, stList_length(components) == 0);
    //Cleanup
    stList_destruct(components);
    stSortedSet_destruct(nodes);
    stList_destruct(chosenEdges);
}

static void testMakeReferenceGreedily(CuTest *testCase) {
    for (int32_t i = 0; i < 100; i++) {
        setup();
        double totalScore;
        stList *reference = makeReferenceGreedily(stubs, chains, zMatrix,
                nodeNumber, &totalScore);
        checkIsValidReference(testCase, reference, totalScore);
        logReference(reference, nodeNumber, zMatrix, totalScore, "just greedy");
        teardown();
    }
}

static void testGibbsSamplingWithSimulatedAnnealing(CuTest *testCase,
        double(*temperatureFn)(double), bool pureGreedy) {
    for (int32_t i = 0; i < 100; i++) {
        setup();
        double totalScore;
        stList *reference = makeReferenceGreedily(stubs, chains, zMatrix,
                nodeNumber, &totalScore);
        checkIsValidReference(testCase, reference, totalScore);
        logReference(reference, nodeNumber, zMatrix, totalScore,
                "pre-annealing");
        int32_t permutations = st_randomInt(0, 100);
        gibbsSamplingWithSimulatedAnnealing(reference, chains, zMatrix,
                permutations, temperatureFn, pureGreedy);
        totalScore = calculateZScoreOfReference(reference, nodeNumber, zMatrix);
        checkIsValidReference(testCase, reference, totalScore);
        logReference(reference, nodeNumber, zMatrix, totalScore,
                "post-annealing");
        teardown();
    }
}

static void testGibbsSamplingWithSimulatedAnnealing_NoExponentiation_Greedy(
        CuTest *testCase) {
    st_logDebug("Running adjacency problem tests using gibbs sampling, but greedy sampling\n");
    testGibbsSamplingWithSimulatedAnnealing(testCase,
            NULL, 1);
}

static void testGibbsSamplingWithSimulatedAnnealing_NoExponentiation(
        CuTest *testCase) {
    st_logDebug("Running adjacency problem tests using gibbs sampling, but no exponentiation\n");
    testGibbsSamplingWithSimulatedAnnealing(testCase,
            NULL, 0);
}

static void testGibbsSamplingWithSimulatedAnnealing_ConstantTemperature(
        CuTest *testCase) {
    st_logDebug("Running adjacency problem tests using gibbs sampling, but with constant temperature\n");
    testGibbsSamplingWithSimulatedAnnealing(testCase, constantTemperatureFn, 0);
}

static void testGibbsSamplingWithSimulatedAnnealing_WithCooling(
        CuTest *testCase) {
    st_logDebug("Running adjacency problem tests using gibbs sampling, with exponentially decreasing temperature function\n");
    testGibbsSamplingWithSimulatedAnnealing(testCase,
            exponentiallyDecreasingTemperatureFn, 0);
}

static double calculateZScoreSlow(int32_t n, int32_t m, int32_t k, double theta) {
    double score = 0.0;
    for(int32_t i=0; i<n; i++) {
        for(int32_t j=0; j<m; j++) {
            score += pow(1.0 - theta, k + i + j);
        }
    }
    return score;
}

static void testCalculateZScore(CuTest *testCase) {
    for(int32_t i=0; i<100; i++) {
        int32_t n = st_randomInt(0, 100);
        int32_t m = st_randomInt(0, 100);
        int32_t k = st_randomInt(1, 10);
        double theta = st_random() > 0.1 ? st_random() : 0.0;
        double zScore = calculateZScore(n, m, k, theta);
        double zScoreSlow = calculateZScoreSlow(n, m, k, theta);
        st_logDebug("The slow computed score: %f the fast computed score: %f, n: %i m: %i k: %i, theta: %lf\n", zScoreSlow, zScore, n, m, k, theta);
        CuAssertDblEquals(testCase, zScoreSlow, zScore, 0.000001);
    }
}

CuSuite* adjacencyProblemTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testMakeReferenceGreedily);
    SUITE_ADD_TEST(suite,
                    testGibbsSamplingWithSimulatedAnnealing_NoExponentiation_Greedy);
    SUITE_ADD_TEST(suite,
                testGibbsSamplingWithSimulatedAnnealing_NoExponentiation);
    SUITE_ADD_TEST(suite,
            testGibbsSamplingWithSimulatedAnnealing_ConstantTemperature);
    SUITE_ADD_TEST(suite, testGibbsSamplingWithSimulatedAnnealing_WithCooling);
    SUITE_ADD_TEST(suite, testCalculateZScore);
    return suite;
}
