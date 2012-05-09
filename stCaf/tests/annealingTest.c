/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stCaf.h"
#include "stPinchGraphs.h"

void stCaf_anneal2(stPinchThreadSet *threadSet, stPinch *(*pinchIterator)(void *), void *extraArg);

void stCaf_annealBetweenAdjacencyComponents2(stPinchThreadSet *threadSet, stPinch *(*pinchIterator)(void *), void *extraArg);

static stPinch *randomPinch(void *extraArg) {
    if(st_random() < 0.01) {
        return NULL;
    }
    static stPinch pinch;
    pinch = stPinchThreadSet_getRandomPinch(extraArg);
    return &pinch;
}

static void testAnnealing(CuTest *testCase) {
    //return;
    for (int32_t test = 0; test < 100; test++) {
        st_logInfo("Starting annealing random test %i\n", test);
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomEmptyGraph();
        stCaf_anneal2(threadSet, randomPinch, threadSet);
    }
}

static void testAnnealingBetweenAdjacencyComponents(CuTest *testCase) {
    //return;
    for (int32_t test = 0; test < 100; test++) {
        st_logInfo("Starting annealing between adjacency components random test %i\n", test);
        stPinchThreadSet *threadSet = stPinchThreadSet_getRandomGraph();
        stCaf_annealBetweenAdjacencyComponents2(threadSet, randomPinch, threadSet);
    }
}

CuSuite* annealingTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testAnnealing);
    SUITE_ADD_TEST(suite, testAnnealingBetweenAdjacencyComponents);
    return suite;
}
