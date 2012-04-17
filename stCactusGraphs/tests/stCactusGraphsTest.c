/*
 * stPinchGraphsTest.c
 *
 *  Created on: 11 Apr 2012
 *      Author: benedictpaten
 *
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stCactusGraphs.h"

static stCactusGraph *g;
static int32_t nO1=1, nO2=2, nO3=3, nO4=4;
static stCactusNode *n1, *n2, *n3;
static stCactusEdgeEnd *e12, *e21, *e23, *e32, *e13, *e31, *e24, *e42;

static void teardown() {
    if(g != NULL) {
        stCactusGraph_destruct(g);
        g = NULL;
    }
}

static void *mergeNodeObjects(void *a, void *b) {
    return a;
}

static void setup() {
    g = stCactusGraph_construct();
    n1 = stCactusNode_construct(g, &nO1);
    n2 = stCactusNode_construct(g, &nO2);
    n3 = stCactusNode_construct(g, &nO3);
    stCactusNode *n4 = stCactusNode_construct(g, &nO4);
    e12 = stCactusEdgeEnd_construct(g, n1, n2, &nO1, &nO2);
    e21 = stCactusEdgeEnd_getOtherEdgeEnd(e12);
    e23 = stCactusEdgeEnd_construct(g, n2, n3, &nO2, &nO3);
    e32 = stCactusEdgeEnd_getOtherEdgeEnd(e23);
    e13 = stCactusEdgeEnd_construct(g, n1, n3, &nO1, &nO3);
    e31 = stCactusEdgeEnd_getOtherEdgeEnd(e13);
    e24 = stCactusEdgeEnd_construct(g, n2, n4, &nO1, &nO4);
    e42 = stCactusEdgeEnd_getOtherEdgeEnd(e24);
    stCactusGraph_collapseToCactus(g, mergeNodeObjects, n1);
}

static void testStCactusNode(CuTest *testCase) {
    setup();
    CuAssertPtrEquals(testCase, &nO1, stCactusNode_getObject(n1));
    CuAssertPtrEquals(testCase, &nO2, stCactusNode_getObject(n2));
    CuAssertPtrEquals(testCase, &nO3, stCactusNode_getObject(n3));
    teardown();
}

CuSuite* stCactusGraphsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testStCactusNode);
    return suite;
}
