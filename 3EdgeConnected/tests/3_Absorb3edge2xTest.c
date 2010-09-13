/*
 * sonLibRandomTest.c
 *
 *  Created on: 22-Jun-2010
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "CuTest.h"
#include "3_Absorb3edge2x.h"

static void teardown() {

}

static void setup() {

}

static void test_3EdgeFunction(CuTest *testCase) {
    /*
     * Exercises the 3-edge function
     */
    setup();
    //stList *vertices = getRandomGraph();

    /*
     * Compute the 3-edge connected components.
     */
    //stList *threeEdgeConnectedComponents = computeThreeEdgeConnectedComponents(vertices);

    /*
     *
     */

    //stList *computeThreeEdgeConnectedComponents(stList *vertices);
    teardown();
}

CuSuite* threeEdgeTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_3EdgeFunction);
    return suite;
}
