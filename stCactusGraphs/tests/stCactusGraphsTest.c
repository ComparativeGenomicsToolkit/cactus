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
static int32_t nO1 = 1, nO2 = 2, nO3 = 3, nO4 = 4, nO5 = 5, nO6 = 6;
static stCactusNode *n1, *n2, *n3, *n4, *n5, *n6;
static stCactusEdgeEnd *e12, *e21, *e23, *e32, *e13, *e31, *e24, *e42, *e11, *e11r, *e45, *e54, *e56, *e65, *e56r,
        *e65r;

static void teardown() {
    if (g != NULL) {
        stCactusGraph_destruct(g);
        g = NULL;
    }
}

static void *mergeNodeObjects(void *a, void *b) {
    return stIntTuple_cmpFn(a, b) == 1 ? b : a;
}

static void setup() {
    g = stCactusGraph_construct();
    n1 = stCactusNode_construct(g, &nO1);
    n2 = stCactusNode_construct(g, &nO2);
    n3 = stCactusNode_construct(g, &nO3);
    n4 = stCactusNode_construct(g, &nO4);
    n5 = stCactusNode_construct(g, &nO5);
    n6 = stCactusNode_construct(g, &nO6);
    e12 = stCactusEdgeEnd_construct(g, n1, n2, &nO1, &nO2);
    e21 = stCactusEdgeEnd_getOtherEdgeEnd(e12);
    e23 = stCactusEdgeEnd_construct(g, n2, n3, &nO2, &nO3);
    e32 = stCactusEdgeEnd_getOtherEdgeEnd(e23);
    e13 = stCactusEdgeEnd_construct(g, n1, n3, &nO1, &nO3);
    e31 = stCactusEdgeEnd_getOtherEdgeEnd(e13);
    e24 = stCactusEdgeEnd_construct(g, n2, n4, &nO1, &nO4);
    e42 = stCactusEdgeEnd_getOtherEdgeEnd(e24);
    e11 = stCactusEdgeEnd_construct(g, n1, n1, &nO1, &nO1);
    e11r = stCactusEdgeEnd_getOtherEdgeEnd(e11);
    e45 = stCactusEdgeEnd_construct(g, n4, n5, &nO4, &nO5);
    e54 = stCactusEdgeEnd_getOtherEdgeEnd(e45);
    e56 = stCactusEdgeEnd_construct(g, n5, n6, &nO5, &nO6);
    e65 = stCactusEdgeEnd_getOtherEdgeEnd(e56);
    e56r = stCactusEdgeEnd_construct(g, n5, n6, &nO5, &nO6);
    e65r = stCactusEdgeEnd_getOtherEdgeEnd(e56r);

    stCactusGraph_collapseToCactus(g, mergeNodeObjects, n1);
}

static void testStCactusNode(CuTest *testCase) {
    setup();
    //Test the get object
    CuAssertPtrEquals(testCase, &nO1, stCactusNode_getObject(n1));
    CuAssertPtrEquals(testCase, &nO2, stCactusNode_getObject(n2));
    CuAssertPtrEquals(testCase, &nO3, stCactusNode_getObject(n3));
    CuAssertPtrEquals(testCase, &nO4, stCactusNode_getObject(n4));

    //Test iterator
    stCactusNode_edgeEndIt it = stCactusNode_getEdgeEndIt(n1);
    CuAssertPtrEquals(testCase, e12, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e13, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e11, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e11r, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNode_edgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n2);
    CuAssertPtrEquals(testCase, e21, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e23, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e24, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNode_edgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n3);
    CuAssertPtrEquals(testCase, e32, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e31, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNode_edgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n4);
    CuAssertPtrEquals(testCase, e42, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e45, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNode_edgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n5);
    CuAssertPtrEquals(testCase, e54, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e56, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e56r, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNode_edgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n6);
    CuAssertPtrEquals(testCase, e65, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e65r, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNode_edgeEndIt_getNext(&it));

    teardown();
}

static void testStCactusGraph(CuTest *testCase) {
    setup();

    //stCactusGraph_getNode
    CuAssertPtrEquals(testCase, n1, stCactusGraph_getNode(g, &nO1));
    CuAssertPtrEquals(testCase, n2, stCactusGraph_getNode(g, &nO2));
    CuAssertPtrEquals(testCase, n3, stCactusGraph_getNode(g, &nO3));
    CuAssertPtrEquals(testCase, n4, stCactusGraph_getNode(g, &nO4));
    CuAssertPtrEquals(testCase, NULL, stCactusGraph_getNode(g, testCase));

    //stCactusGraphNodeIterator_construct
    //stCactusGraphNodeIterator_getNext
    //stCactusGraphNodeIterator_destruct
    stCactusGraphNodeIterator *nodeIt = stCactusGraphNodeIterator_construct(g);
    stCactusNode *n;
    stSortedSet *nodes = stSortedSet_construct();
    while ((n = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        CuAssertTrue(testCase, stSortedSet_search(nodes, n) == NULL);
        CuAssertTrue(testCase, n == n1 || n == n2 || n == n3 || n == n4 || n == n5 || n == n6);
        stSortedSet_insert(nodes, n);
    }
    CuAssertIntEquals(testCase, 6, stSortedSet_size(nodes));
    stCactusGraphNodeIterator_destruct(nodeIt);

    teardown();
}

static void basicTestsOfEdgeEnds(CuTest *testCase) {
    //stCactusEdgeEnd_getObject
    CuAssertPtrEquals(testCase, &nO1, stCactusEdgeEnd_getObject(e12));
    CuAssertPtrEquals(testCase, &nO1, stCactusEdgeEnd_getObject(e13));
    CuAssertPtrEquals(testCase, &nO1, stCactusEdgeEnd_getObject(e11));
    CuAssertPtrEquals(testCase, &nO1, stCactusEdgeEnd_getObject(e11r));
    CuAssertPtrEquals(testCase, &nO2, stCactusEdgeEnd_getObject(e21));
    CuAssertPtrEquals(testCase, &nO2, stCactusEdgeEnd_getObject(e23));
    CuAssertPtrEquals(testCase, &nO2, stCactusEdgeEnd_getObject(e24));
    CuAssertPtrEquals(testCase, &nO3, stCactusEdgeEnd_getObject(e32));
    CuAssertPtrEquals(testCase, &nO3, stCactusEdgeEnd_getObject(e31));
    CuAssertPtrEquals(testCase, &nO4, stCactusEdgeEnd_getObject(e42));

    CuAssertPtrEquals(testCase, &nO4, stCactusEdgeEnd_getObject(e45));
    CuAssertPtrEquals(testCase, &nO5, stCactusEdgeEnd_getObject(e54));
    CuAssertPtrEquals(testCase, &nO5, stCactusEdgeEnd_getObject(e56));
    CuAssertPtrEquals(testCase, &nO5, stCactusEdgeEnd_getObject(e56r));
    CuAssertPtrEquals(testCase, &nO6, stCactusEdgeEnd_getObject(e65r));
    CuAssertPtrEquals(testCase, &nO6, stCactusEdgeEnd_getObject(e65));

    //stCactusEdgeEnd_getNode
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getNode(e12));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getNode(e13));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getNode(e11));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getNode(e11r));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e21));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e23));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e24));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getNode(e32));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getNode(e31));
    CuAssertPtrEquals(testCase, n4, stCactusEdgeEnd_getNode(e42));

    //stCactusEdgeEnd_getOtherNode
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getOtherNode(e12));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getOtherNode(e13));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getOtherNode(e11));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getOtherNode(e11r));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getOtherNode(e21));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getOtherNode(e23));
    CuAssertPtrEquals(testCase, n4, stCactusEdgeEnd_getOtherNode(e24));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getOtherNode(e32));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getOtherNode(e31));
    CuAssertPtrEquals(testCase, n4, stCactusEdgeEnd_getOtherNode(e42));

    //stCactusEdgeEnd_getOtherEdgeEnd
    CuAssertPtrEquals(testCase, e21, stCactusEdgeEnd_getOtherEdgeEnd(e12));
    CuAssertPtrEquals(testCase, e31, stCactusEdgeEnd_getOtherEdgeEnd(e13));
    CuAssertPtrEquals(testCase, e11r, stCactusEdgeEnd_getOtherEdgeEnd(e11));
    CuAssertPtrEquals(testCase, e11, stCactusEdgeEnd_getOtherEdgeEnd(e11r));
    CuAssertPtrEquals(testCase, e12, stCactusEdgeEnd_getOtherEdgeEnd(e21));
    CuAssertPtrEquals(testCase, e32, stCactusEdgeEnd_getOtherEdgeEnd(e23));
    CuAssertPtrEquals(testCase, e42, stCactusEdgeEnd_getOtherEdgeEnd(e24));
    CuAssertPtrEquals(testCase, e23, stCactusEdgeEnd_getOtherEdgeEnd(e32));
    CuAssertPtrEquals(testCase, e13, stCactusEdgeEnd_getOtherEdgeEnd(e31));
    CuAssertPtrEquals(testCase, e24, stCactusEdgeEnd_getOtherEdgeEnd(e42));
}

static void chainStructureTests(CuTest *testCase) {
    //Now test the chain structure
    //stCactusEdgeEnd_getLink
    CuAssertPtrEquals(testCase, e13, stCactusEdgeEnd_getLink(e12));
    CuAssertPtrEquals(testCase, e12, stCactusEdgeEnd_getLink(e13));
    CuAssertPtrEquals(testCase, e11r, stCactusEdgeEnd_getLink(e11));
    CuAssertPtrEquals(testCase, e11, stCactusEdgeEnd_getLink(e11r));
    CuAssertPtrEquals(testCase, e23, stCactusEdgeEnd_getLink(e21));
    CuAssertPtrEquals(testCase, e21, stCactusEdgeEnd_getLink(e23));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e24));
    CuAssertPtrEquals(testCase, e31, stCactusEdgeEnd_getLink(e32));
    CuAssertPtrEquals(testCase, e31, stCactusEdgeEnd_getLink(e31));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e42));

    //stCactusEdgeEnd_getLinkOrientation
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_getLinkOrientation(e12));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e13));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_getLinkOrientation(e11));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e11r));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_getLinkOrientation(e21));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e23));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e24));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e32));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_getLinkOrientation(e31));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e42));

    //stCactusEdgeEnd_isChainEnd
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e12));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e13));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e11));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e11r));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e21));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e23));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e24));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e32));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e31));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e42));
}

static void testStCactusEdgeEnd(CuTest *testCase) {
    setup();
    basicTestsOfEdgeEnds(testCase);
    chainStructureTests(testCase);
    teardown();
}

static void testStCactusGraph_unmarkAndMarkCycles(CuTest *testCase) {
    setup();
    stCactusGraph_unmarkCycles(g);
    basicTestsOfEdgeEnds(testCase);

    //Now test the chain structure
    //stCactusEdgeEnd_getLink
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e12));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e13));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e11));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e11r));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e21));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e23));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e24));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e32));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e31));
    CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(e42));

    //stCactusEdgeEnd_getLinkOrientation
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e12));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e13));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e11));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e11r));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e21));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e23));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e24));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e32));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e31));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e42));

    //stCactusEdgeEnd_isChainEnd
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e12));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e13));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e11));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e11r));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e21));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e23));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e24));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e32));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e31));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e42));

    stCactusGraph_markCycles(g, n1);
    basicTestsOfEdgeEnds(testCase);
    chainStructureTests(testCase);

    teardown();
}

static void testStCactusGraph_collapseBridges(CuTest *testCase) {
    setup();
    stCactusGraph_collapseBridges(g, n1, mergeNodeObjects);
    n2 = stCactusGraph_getNode(g, &nO2);

    //Check pointers of nodes okay
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e21));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e23));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e24));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e42));

    //Check chains
    CuAssertPtrEquals(testCase, e23, stCactusEdgeEnd_getLink(e21));
    CuAssertPtrEquals(testCase, e21, stCactusEdgeEnd_getLink(e23));
    CuAssertPtrEquals(testCase, e42, stCactusEdgeEnd_getLink(e24));
    CuAssertPtrEquals(testCase, e24, stCactusEdgeEnd_getLink(e42));

    //Chain orientation
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_getLinkOrientation(e21));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e23));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_getLinkOrientation(e24));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(e42));

    //Chain end status
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e21));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e23));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e24));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e42));

    //Check edge end it.
    stCactusNode_edgeEndIt it = stCactusNode_getEdgeEndIt(n2);
    CuAssertPtrEquals(testCase, e21, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e23, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e24, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e42, stCactusNode_edgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNode_edgeEndIt_getNext(&it));

    //Check graph node it.
    stCactusGraphNodeIterator *nodeIt = stCactusGraphNodeIterator_construct(g);
    stCactusNode *n;
    stSortedSet *nodes = stSortedSet_construct();
    while ((n = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        CuAssertTrue(testCase, stSortedSet_search(nodes, n) == NULL);
        CuAssertTrue(testCase, n == n1 || n == n2 || n == n3);
        stSortedSet_insert(nodes, n);
    }
    CuAssertIntEquals(testCase, 3, stSortedSet_size(nodes));
    stCactusGraphNodeIterator_destruct(nodeIt);

    teardown();
}

static void testStCactusGraph_randomTest(CuTest *testCase) {
    //Creates a problem instances, then checks graph is okay by checking every edge
    //is properly connected, with right number of nodes and that everyone is in a chain
}

CuSuite* stCactusGraphsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testStCactusNode);
    SUITE_ADD_TEST(suite, testStCactusEdgeEnd);
    SUITE_ADD_TEST(suite, testStCactusGraph);
    SUITE_ADD_TEST(suite, testStCactusGraph_unmarkAndMarkCycles);
    SUITE_ADD_TEST(suite, testStCactusGraph_collapseBridges);
    SUITE_ADD_TEST(suite, testStCactusGraph_randomTest);
    return suite;
}
