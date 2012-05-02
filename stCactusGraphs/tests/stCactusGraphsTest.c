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
#include <stdlib.h>

static stCactusGraph *g;
static int32_t nO1 = 1, nO2 = 2, nO3 = 3, nO4 = 4, nO5 = 5, nO6 = 6, nO7 = 7, nO8 = 8, nO9 = 9;
static stCactusNode *n1, *n2, *n3, *n4, *n5, *n6, *n7, *n8, *n9;
static stCactusEdgeEnd *e12, *e21, *e23, *e32, *e13, *e31, *e24, *e42, *e11, *e11r, *e45, *e54, *e56, *e65, *e56r,
        *e65r, *e58, *e85, *e49, *e94, *e37, *e73, *e37r, *e73r;

static void teardown() {
    if (g != NULL) {
        stCactusGraph_destruct(g);
        g = NULL;
    }
}

static void *mergeNodeObjects(void *a, void *b) {
    int32_t i = *((int32_t *) a);
    int32_t j = *((int32_t *) b);
    return i > j ? b : a;
}

static void setup() {
    g = stCactusGraph_construct();
    n1 = stCactusNode_construct(g, &nO1);
    n2 = stCactusNode_construct(g, &nO2);
    n3 = stCactusNode_construct(g, &nO3);
    n4 = stCactusNode_construct(g, &nO4);
    n5 = stCactusNode_construct(g, &nO5);
    n6 = stCactusNode_construct(g, &nO6);
    n7 = stCactusNode_construct(g, &nO7);
    n8 = stCactusNode_construct(g, &nO8);
    n9 = stCactusNode_construct(g, &nO9);
    e12 = stCactusEdgeEnd_construct(g, n1, n2, &nO1, &nO2);
    e21 = stCactusEdgeEnd_getOtherEdgeEnd(e12);
    e23 = stCactusEdgeEnd_construct(g, n2, n3, &nO2, &nO3);
    e32 = stCactusEdgeEnd_getOtherEdgeEnd(e23);
    e13 = stCactusEdgeEnd_construct(g, n1, n3, &nO1, &nO3);
    e31 = stCactusEdgeEnd_getOtherEdgeEnd(e13);
    e24 = stCactusEdgeEnd_construct(g, n2, n4, &nO2, &nO4);
    e42 = stCactusEdgeEnd_getOtherEdgeEnd(e24);
    e11 = stCactusEdgeEnd_construct(g, n1, n1, &nO1, &nO1);
    e11r = stCactusEdgeEnd_getOtherEdgeEnd(e11);
    e45 = stCactusEdgeEnd_construct(g, n4, n5, &nO4, &nO5);
    e54 = stCactusEdgeEnd_getOtherEdgeEnd(e45);
    e56 = stCactusEdgeEnd_construct(g, n5, n6, &nO5, &nO6);
    e65 = stCactusEdgeEnd_getOtherEdgeEnd(e56);
    e56r = stCactusEdgeEnd_construct(g, n5, n6, &nO5, &nO6);
    e65r = stCactusEdgeEnd_getOtherEdgeEnd(e56r);
    e58 = stCactusEdgeEnd_construct(g, n5, n8, &nO5, &nO8);
    e85 = stCactusEdgeEnd_getOtherEdgeEnd(e58);
    e49 = stCactusEdgeEnd_construct(g, n4, n9, &nO4, &nO9);
    e94 = stCactusEdgeEnd_getOtherEdgeEnd(e49);
    e37 = stCactusEdgeEnd_construct(g, n3, n7, &nO3, &nO7);
    e73 = stCactusEdgeEnd_getOtherEdgeEnd(e37);
    e37r = stCactusEdgeEnd_construct(g, n3, n7, &nO3, &nO7);
    e73r = stCactusEdgeEnd_getOtherEdgeEnd(e37r);

    stCactusGraph_collapseToCactus(g, mergeNodeObjects, n1);
}

static void invariantNodeTests(CuTest *testCase) {
    //Test the get object
    CuAssertPtrEquals(testCase, &nO1, stCactusNode_getObject(n1));
    CuAssertPtrEquals(testCase, &nO2, stCactusNode_getObject(n2));
    CuAssertPtrEquals(testCase, &nO3, stCactusNode_getObject(n3));
    CuAssertPtrEquals(testCase, &nO5, stCactusNode_getObject(n5));
    CuAssertPtrEquals(testCase, &nO6, stCactusNode_getObject(n6));
    CuAssertPtrEquals(testCase, &nO7, stCactusNode_getObject(n7));

    //Test iterator
    stCactusNodeEdgeEndIt it = stCactusNode_getEdgeEndIt(n1);
    CuAssertPtrEquals(testCase, e12, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e13, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e11, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e11r, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n3);
    CuAssertPtrEquals(testCase, e32, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e31, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e37, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e37r, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n5);
    CuAssertPtrEquals(testCase, e54, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e56, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e56r, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e58, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n6);
    CuAssertPtrEquals(testCase, e65, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e65r, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n7);
    CuAssertPtrEquals(testCase, e73, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e73r, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));
}

static void testStCactusNode(CuTest *testCase) {
    setup();
    invariantNodeTests(testCase);

    //Test the get object
    CuAssertPtrEquals(testCase, &nO4, stCactusNode_getObject(n4));
    CuAssertPtrEquals(testCase, &nO8, stCactusNode_getObject(n8));
    CuAssertPtrEquals(testCase, &nO9, stCactusNode_getObject(n9));

    //Test iterator
    stCactusNodeEdgeEndIt it = stCactusNode_getEdgeEndIt(n2);
    CuAssertPtrEquals(testCase, e21, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e23, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e24, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n4);
    CuAssertPtrEquals(testCase, e42, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e45, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e49, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n8);
    CuAssertPtrEquals(testCase, e85, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));

    it = stCactusNode_getEdgeEndIt(n9);
    CuAssertPtrEquals(testCase, e94, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));

    teardown();
}

static void testStCactusGraph(CuTest *testCase) {
    setup();

    //stCactusGraph_getNode
    CuAssertPtrEquals(testCase, n1, stCactusGraph_getNode(g, &nO1));
    CuAssertPtrEquals(testCase, n2, stCactusGraph_getNode(g, &nO2));
    CuAssertPtrEquals(testCase, n3, stCactusGraph_getNode(g, &nO3));
    CuAssertPtrEquals(testCase, n4, stCactusGraph_getNode(g, &nO4));
    CuAssertPtrEquals(testCase, n5, stCactusGraph_getNode(g, &nO5));
    CuAssertPtrEquals(testCase, n6, stCactusGraph_getNode(g, &nO6));
    CuAssertPtrEquals(testCase, n7, stCactusGraph_getNode(g, &nO7));
    CuAssertPtrEquals(testCase, n8, stCactusGraph_getNode(g, &nO8));
    CuAssertPtrEquals(testCase, n9, stCactusGraph_getNode(g, &nO9));
    CuAssertPtrEquals(testCase, NULL, stCactusGraph_getNode(g, testCase));

    //stCactusGraphNodeIterator_construct
    //stCactusGraphNodeIterator_getNext
    //stCactusGraphNodeIterator_destruct
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(g);
    stCactusNode *n;
    stSortedSet *nodes = stSortedSet_construct();
    while ((n = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        CuAssertTrue(testCase, stSortedSet_search(nodes, n) == NULL);
        CuAssertTrue(testCase,
                n == n1 || n == n2 || n == n3 || n == n4 || n == n5 || n == n6 || n == n7 || n == n8 || n == n9);
        stSortedSet_insert(nodes, n);
    }
    CuAssertIntEquals(testCase, 9, stSortedSet_size(nodes));
    stCactusGraphNodeIterator_destruct(nodeIt);
    stSortedSet_destruct(nodes);
    teardown();
}

static void invariantEdgeEndTests(CuTest *testCase) {
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
    CuAssertPtrEquals(testCase, &nO3, stCactusEdgeEnd_getObject(e37));
    CuAssertPtrEquals(testCase, &nO3, stCactusEdgeEnd_getObject(e37r));
    CuAssertPtrEquals(testCase, &nO4, stCactusEdgeEnd_getObject(e42));
    CuAssertPtrEquals(testCase, &nO4, stCactusEdgeEnd_getObject(e45));
    CuAssertPtrEquals(testCase, &nO4, stCactusEdgeEnd_getObject(e49));
    CuAssertPtrEquals(testCase, &nO5, stCactusEdgeEnd_getObject(e54));
    CuAssertPtrEquals(testCase, &nO5, stCactusEdgeEnd_getObject(e56));
    CuAssertPtrEquals(testCase, &nO5, stCactusEdgeEnd_getObject(e56r));
    CuAssertPtrEquals(testCase, &nO5, stCactusEdgeEnd_getObject(e58));
    CuAssertPtrEquals(testCase, &nO6, stCactusEdgeEnd_getObject(e65r));
    CuAssertPtrEquals(testCase, &nO6, stCactusEdgeEnd_getObject(e65));
    CuAssertPtrEquals(testCase, &nO7, stCactusEdgeEnd_getObject(e73));
    CuAssertPtrEquals(testCase, &nO7, stCactusEdgeEnd_getObject(e73r));
    CuAssertPtrEquals(testCase, &nO8, stCactusEdgeEnd_getObject(e85));
    CuAssertPtrEquals(testCase, &nO9, stCactusEdgeEnd_getObject(e94));

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
    CuAssertPtrEquals(testCase, e73, stCactusEdgeEnd_getOtherEdgeEnd(e37));
    CuAssertPtrEquals(testCase, e73r, stCactusEdgeEnd_getOtherEdgeEnd(e37r));
    CuAssertPtrEquals(testCase, e24, stCactusEdgeEnd_getOtherEdgeEnd(e42));
    CuAssertPtrEquals(testCase, e54, stCactusEdgeEnd_getOtherEdgeEnd(e45));
    CuAssertPtrEquals(testCase, e94, stCactusEdgeEnd_getOtherEdgeEnd(e49));
    CuAssertPtrEquals(testCase, e45, stCactusEdgeEnd_getOtherEdgeEnd(e54));
    CuAssertPtrEquals(testCase, e65, stCactusEdgeEnd_getOtherEdgeEnd(e56));
    CuAssertPtrEquals(testCase, e65r, stCactusEdgeEnd_getOtherEdgeEnd(e56r));
    CuAssertPtrEquals(testCase, e85, stCactusEdgeEnd_getOtherEdgeEnd(e58));
    CuAssertPtrEquals(testCase, e56r, stCactusEdgeEnd_getOtherEdgeEnd(e65r));
    CuAssertPtrEquals(testCase, e56, stCactusEdgeEnd_getOtherEdgeEnd(e65));
    CuAssertPtrEquals(testCase, e37, stCactusEdgeEnd_getOtherEdgeEnd(e73));
    CuAssertPtrEquals(testCase, e37r, stCactusEdgeEnd_getOtherEdgeEnd(e73r));
    CuAssertPtrEquals(testCase, e58, stCactusEdgeEnd_getOtherEdgeEnd(e85));
    CuAssertPtrEquals(testCase, e49, stCactusEdgeEnd_getOtherEdgeEnd(e94));

    //stCactusEdgeEnd_getNode (those which don't change)
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getNode(e12));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getNode(e13));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getNode(e11));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getNode(e11r));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e21));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e23));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e24));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getNode(e32));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getNode(e31));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getNode(e37));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getNode(e37r));
    CuAssertPtrEquals(testCase, n5, stCactusEdgeEnd_getNode(e54));
    CuAssertPtrEquals(testCase, n5, stCactusEdgeEnd_getNode(e56));
    CuAssertPtrEquals(testCase, n5, stCactusEdgeEnd_getNode(e56r));
    CuAssertPtrEquals(testCase, n5, stCactusEdgeEnd_getNode(e58));
    CuAssertPtrEquals(testCase, n6, stCactusEdgeEnd_getNode(e65r));
    CuAssertPtrEquals(testCase, n6, stCactusEdgeEnd_getNode(e65));
    CuAssertPtrEquals(testCase, n7, stCactusEdgeEnd_getNode(e73));
    CuAssertPtrEquals(testCase, n7, stCactusEdgeEnd_getNode(e73r));

    //stCactusEdgeEnd_getOtherNode (those which don't change)
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getOtherNode(e12));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getOtherNode(e13));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getOtherNode(e11));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getOtherNode(e11r));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getOtherNode(e21));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getOtherNode(e23));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getOtherNode(e32));
    CuAssertPtrEquals(testCase, n1, stCactusEdgeEnd_getOtherNode(e31));
    CuAssertPtrEquals(testCase, n7, stCactusEdgeEnd_getOtherNode(e37));
    CuAssertPtrEquals(testCase, n7, stCactusEdgeEnd_getOtherNode(e37r));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getOtherNode(e42));
    CuAssertPtrEquals(testCase, n5, stCactusEdgeEnd_getOtherNode(e45));
    CuAssertPtrEquals(testCase, n6, stCactusEdgeEnd_getOtherNode(e56));
    CuAssertPtrEquals(testCase, n6, stCactusEdgeEnd_getOtherNode(e56r));
    CuAssertPtrEquals(testCase, n5, stCactusEdgeEnd_getOtherNode(e65r));
    CuAssertPtrEquals(testCase, n5, stCactusEdgeEnd_getOtherNode(e65));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getOtherNode(e73));
    CuAssertPtrEquals(testCase, n3, stCactusEdgeEnd_getOtherNode(e73r));
    CuAssertPtrEquals(testCase, n5, stCactusEdgeEnd_getOtherNode(e85));
}

static void edgeEndTests(CuTest *testCase) {
    invariantEdgeEndTests(testCase);

    //stCactusEdgeEnd_getNode (those which change)
    CuAssertPtrEquals(testCase, n4, stCactusEdgeEnd_getNode(e42));
    CuAssertPtrEquals(testCase, n4, stCactusEdgeEnd_getNode(e45));
    CuAssertPtrEquals(testCase, n4, stCactusEdgeEnd_getNode(e49));
    CuAssertPtrEquals(testCase, n8, stCactusEdgeEnd_getNode(e85));
    CuAssertPtrEquals(testCase, n9, stCactusEdgeEnd_getNode(e94));

    //stCactusEdgeEnd_getOtherNode (those which change)
    CuAssertPtrEquals(testCase, n4, stCactusEdgeEnd_getOtherNode(e24));
    CuAssertPtrEquals(testCase, n4, stCactusEdgeEnd_getOtherNode(e94));
    CuAssertPtrEquals(testCase, n8, stCactusEdgeEnd_getOtherNode(e58));
    CuAssertPtrEquals(testCase, n9, stCactusEdgeEnd_getOtherNode(e49));
    CuAssertPtrEquals(testCase, n4, stCactusEdgeEnd_getOtherNode(e54));
}

static void invariantChainStructureTests(CuTest *testCase) {
    //Now test the chain structure
    //stCactusEdgeEnd_getLink
    CuAssertPtrEquals(testCase, e13, stCactusEdgeEnd_getLink(e12));
    CuAssertPtrEquals(testCase, e12, stCactusEdgeEnd_getLink(e13));
    CuAssertPtrEquals(testCase, e11r, stCactusEdgeEnd_getLink(e11));
    CuAssertPtrEquals(testCase, e11, stCactusEdgeEnd_getLink(e11r));
    CuAssertPtrEquals(testCase, e23, stCactusEdgeEnd_getLink(e21));
    CuAssertPtrEquals(testCase, e21, stCactusEdgeEnd_getLink(e23));
    CuAssertPtrEquals(testCase, e31, stCactusEdgeEnd_getLink(e32));
    CuAssertPtrEquals(testCase, e32, stCactusEdgeEnd_getLink(e31));
    CuAssertPtrEquals(testCase, e37r, stCactusEdgeEnd_getLink(e37));
    CuAssertPtrEquals(testCase, e37, stCactusEdgeEnd_getLink(e37r));
    CuAssertPtrEquals(testCase, e56r, stCactusEdgeEnd_getLink(e56));
    CuAssertPtrEquals(testCase, e56, stCactusEdgeEnd_getLink(e56r));
    CuAssertPtrEquals(testCase, e65, stCactusEdgeEnd_getLink(e65r));
    CuAssertPtrEquals(testCase, e65r, stCactusEdgeEnd_getLink(e65));
    CuAssertPtrEquals(testCase, e73r, stCactusEdgeEnd_getLink(e73));
    CuAssertPtrEquals(testCase, e73, stCactusEdgeEnd_getLink(e73r));

    //stCactusEdgeEnd_getLinkOrientation (first of link edges)
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e12) != stCactusEdgeEnd_getLinkOrientation(e13));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e11) != stCactusEdgeEnd_getLinkOrientation(e11r));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e23) != stCactusEdgeEnd_getLinkOrientation(e21));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e31) != stCactusEdgeEnd_getLinkOrientation(e32));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e37) != stCactusEdgeEnd_getLinkOrientation(e37r));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e56) != stCactusEdgeEnd_getLinkOrientation(e56r));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e65) != stCactusEdgeEnd_getLinkOrientation(e65r));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e73) != stCactusEdgeEnd_getLinkOrientation(e73r));

    //stCactusEdgeEnd_getLinkOrientation (now of standard edges)
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e12) != stCactusEdgeEnd_getLinkOrientation(e21));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e11) != stCactusEdgeEnd_getLinkOrientation(e11r));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e23) != stCactusEdgeEnd_getLinkOrientation(e32));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e31) != stCactusEdgeEnd_getLinkOrientation(e13));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e37) != stCactusEdgeEnd_getLinkOrientation(e73));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e56) != stCactusEdgeEnd_getLinkOrientation(e65));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e56r) != stCactusEdgeEnd_getLinkOrientation(e65r));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e37r) != stCactusEdgeEnd_getLinkOrientation(e73r));

    //stCactusEdgeEnd_isChainEnd
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e12));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e13));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e11));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e11r));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e21));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e23));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e32));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e31));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e37));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e37r));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e56));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e56r));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e58));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e65r));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e65));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e73));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e73r));
}

static void chainStructureTests(CuTest *testCase) {
    invariantChainStructureTests(testCase);
    //Check the chain structure of the bridge ends is uninitialised
    stCactusEdgeEnd *edgeEnds[] = { e24, e42, e45, e49, e54, e58, e85, e94 };
    for (int32_t i = 0; i < 8; i++) {
        CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(edgeEnds[i]));
        CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(edgeEnds[i]));
        CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(edgeEnds[i]));
    }
}

static void testStCactusEdgeEnd(CuTest *testCase) {
    setup();
    edgeEndTests(testCase);
    chainStructureTests(testCase);
    teardown();
}

static void testStCactusGraph_unmarkAndMarkCycles(CuTest *testCase) {
    setup();
    stCactusGraph_unmarkCycles(g);
    edgeEndTests(testCase);
    invariantNodeTests(testCase);

    //Now test the chain structure
    //stCactusEdgeEnd_getLink
    stCactusEdgeEnd *edgeEnds[] = { e12, e13, e11, e11r, e21, e23, e24, e32, e31, e37, e37r, e42, e45, e49, e54, e56,
            e56r, e58, e65, e65r, e73, e73r, e85, e94 };

    for (int32_t i = 0; i < 24; i++) {
        CuAssertPtrEquals(testCase, NULL, stCactusEdgeEnd_getLink(edgeEnds[i]));
        CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_getLinkOrientation(edgeEnds[i]));
        CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(edgeEnds[i]));
    }

    stCactusGraph_markCycles(g, n1);
    edgeEndTests(testCase);
    chainStructureTests(testCase);
    invariantNodeTests(testCase);

    teardown();
}

static void testStCactusGraph_collapseBridges(CuTest *testCase) {
    setup();
    stCactusGraph_collapseBridges(g, n1, mergeNodeObjects);
    n2 = stCactusGraph_getNode(g, &nO2);

    invariantEdgeEndTests(testCase);
    invariantChainStructureTests(testCase);
    invariantNodeTests(testCase);

    //stCactusEdgeEnd_getNode (those which change)
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e42));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e45));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e49));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e85));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getNode(e94));

    //stCactusEdgeEnd_getOtherNode (those which change)
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getOtherNode(e24));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getOtherNode(e94));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getOtherNode(e58));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getOtherNode(e49));
    CuAssertPtrEquals(testCase, n2, stCactusEdgeEnd_getOtherNode(e54));

    //Check chains
    //stCactusEdgeEnd_getLink
    CuAssertPtrEquals(testCase, e42, stCactusEdgeEnd_getLink(e24));
    CuAssertPtrEquals(testCase, e24, stCactusEdgeEnd_getLink(e42));
    CuAssertPtrEquals(testCase, e85, stCactusEdgeEnd_getLink(e45));
    CuAssertPtrEquals(testCase, e94, stCactusEdgeEnd_getLink(e49));
    CuAssertPtrEquals(testCase, e58, stCactusEdgeEnd_getLink(e54));
    CuAssertPtrEquals(testCase, e54, stCactusEdgeEnd_getLink(e58));
    CuAssertPtrEquals(testCase, e45, stCactusEdgeEnd_getLink(e85));
    CuAssertPtrEquals(testCase, e49, stCactusEdgeEnd_getLink(e94));

    //stCactusEdgeEnd_getLinkOrientation
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e42) != stCactusEdgeEnd_getLinkOrientation(e24));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e85) != stCactusEdgeEnd_getLinkOrientation(e45));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e94) != stCactusEdgeEnd_getLinkOrientation(e49));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e58) != stCactusEdgeEnd_getLinkOrientation(e54));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e54) != stCactusEdgeEnd_getLinkOrientation(e45));
    CuAssertTrue(testCase, stCactusEdgeEnd_getLinkOrientation(e58) != stCactusEdgeEnd_getLinkOrientation(e85));

    //stCactusEdgeEnd_isChainEnd
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e24));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e42));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e45));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e49));
    CuAssertIntEquals(testCase, 0, stCactusEdgeEnd_isChainEnd(e54));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e85));
    CuAssertIntEquals(testCase, 1, stCactusEdgeEnd_isChainEnd(e94));

    //Check edge end it.
    stCactusNodeEdgeEndIt it = stCactusNode_getEdgeEndIt(n2);
    CuAssertPtrEquals(testCase, e21, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e23, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e24, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e42, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e45, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e49, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e94, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, e85, stCactusNodeEdgeEndIt_getNext(&it));
    CuAssertPtrEquals(testCase, NULL, stCactusNodeEdgeEndIt_getNext(&it));

    //Check graph node it.
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(g);
    stCactusNode *n;
    stSortedSet *nodes = stSortedSet_construct();
    while ((n = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        CuAssertTrue(testCase, stSortedSet_search(nodes, n) == NULL);
        CuAssertTrue(testCase, n == n1 || n == n2 || n == n3 || n == n5 || n == n6 || n == n7);
        stSortedSet_insert(nodes, n);
    }
    CuAssertIntEquals(testCase, 6, stSortedSet_size(nodes));
    stCactusGraphNodeIterator_destruct(nodeIt);
    stSortedSet_destruct(nodes);
    teardown();
}

static void testStCactusGraph_randomTest(CuTest *testCase) {
    //return;
    //Creates a problem instances, then checks graph is okay by checking every edge
    //is properly connected, with right number of nodes and that everyone is in a chain
    for (int32_t test = 0; test < 100; test++) {
        int32_t nodeNumber = st_randomInt(0, 1000);
        int32_t edgeNumber = nodeNumber > 0 ? st_randomInt(0, 1000) : 0;
        st_logInfo("We have %i edges and %i nodes in random test\n", edgeNumber, nodeNumber);
        stCactusGraph *g2 = stCactusGraph_construct();
        stList *nodeObjects = stList_construct3(0, free);
        for (int32_t i = 0; i < nodeNumber; i++) {
            int32_t *j = st_malloc(sizeof(int32_t));
            j[0] = i;
            stCactusNode_construct(g2, j);
            stList_append(nodeObjects, j);
        }
        stSortedSet *edgeEnds = stSortedSet_construct();
        if (nodeNumber > 0) {
            stList *includedNodeObjects = stList_construct(); //Edge construction ensures there is just one component containing edges.
            stList_append(includedNodeObjects, st_randomChoice(nodeObjects));
            for (int32_t i = 0; i < edgeNumber; i++) {
                void *nodeObject1 = st_randomChoice(includedNodeObjects);
                void *nodeObject2 = st_randomChoice(nodeObjects);
                if (!stList_contains(includedNodeObjects, nodeObject2)) {
                    stList_append(includedNodeObjects, nodeObject2);
                }
                stCactusNode *node1 = stCactusGraph_getNode(g2, nodeObject1);
                stCactusNode *node2 = stCactusGraph_getNode(g2, nodeObject2);
                stCactusEdgeEnd *edgeEnd = stCactusEdgeEnd_construct(g2, node1, node2, nodeObject1, nodeObject2);
                stSortedSet_insert(edgeEnds, edgeEnd);
                stSortedSet_insert(edgeEnds, stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd));
            }
            stCactusNode *startNode = stCactusGraph_getNode(g2, st_randomChoice(includedNodeObjects));
            stList_destruct(includedNodeObjects);
            CuAssertTrue(testCase, startNode != NULL);
            stCactusGraph_collapseToCactus(g2, mergeNodeObjects, startNode);
            stCactusGraph_collapseBridges(g2, startNode, mergeNodeObjects);
        }
        //Now iterate through nodes and check chains
        stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(g2);
        stCactusNode *node;
        stSortedSet *nodesInFinalGraph = stSortedSet_construct();
        while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
            //Check node has valid object
            CuAssertTrue(testCase, stList_contains(nodeObjects, stCactusNode_getObject(node)));
            CuAssertTrue(testCase, stSortedSet_search(nodesInFinalGraph, node) == NULL);
            stSortedSet_insert(nodesInFinalGraph, node);
        }
        stCactusGraphNodeIterator_destruct(nodeIt);
        nodeIt = stCactusGraphNodeIterator_construct(g2);
        int32_t edgeEndNumber = 0;
        while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
            CuAssertPtrEquals(testCase, node, stSortedSet_search(nodesInFinalGraph, node));
            stCactusNodeEdgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(node);
            stCactusEdgeEnd *edgeEnd;
            while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeEndIt)) != NULL) {
                edgeEndNumber++;
                //Check node connectivity
                CuAssertPtrEquals(testCase, edgeEnd, stSortedSet_search(edgeEnds, edgeEnd));
                CuAssertPtrEquals(testCase, node, stCactusEdgeEnd_getNode(edgeEnd));
                //Check other edge end connectivity
                stCactusEdgeEnd *otherEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd);
                CuAssertPtrEquals(testCase, otherEdgeEnd, stSortedSet_search(edgeEnds, otherEdgeEnd));
                CuAssertPtrEquals(testCase, edgeEnd, stCactusEdgeEnd_getOtherEdgeEnd(otherEdgeEnd));
                CuAssertTrue(testCase,
                        stSortedSet_search(nodesInFinalGraph, stCactusEdgeEnd_getNode(otherEdgeEnd)) != NULL);
                CuAssertTrue(testCase,
                        stCactusEdgeEnd_getLinkOrientation(edgeEnd) != stCactusEdgeEnd_getLinkOrientation(otherEdgeEnd));
                CuAssertPtrEquals(testCase, stCactusEdgeEnd_getOtherNode(edgeEnd),
                        stCactusEdgeEnd_getNode(otherEdgeEnd));
                //Check chain connectivity
                stCactusEdgeEnd *linkedEdgeEnd = stCactusEdgeEnd_getLink(edgeEnd);
                CuAssertPtrEquals(testCase, linkedEdgeEnd, stSortedSet_search(edgeEnds, linkedEdgeEnd));
                CuAssertPtrEquals(testCase, node, stCactusEdgeEnd_getNode(linkedEdgeEnd));
                CuAssertPtrEquals(testCase, edgeEnd, stCactusEdgeEnd_getLink(linkedEdgeEnd));
                CuAssertTrue(testCase, stCactusEdgeEnd_isChainEnd(edgeEnd) == stCactusEdgeEnd_isChainEnd(linkedEdgeEnd));
                CuAssertTrue(
                        testCase,
                        stCactusEdgeEnd_getLinkOrientation(edgeEnd)
                                != stCactusEdgeEnd_getLinkOrientation(linkedEdgeEnd));
                //Properties that must be true if chain is a cycle
                CuAssertTrue(testCase, linkedEdgeEnd != edgeEnd);
                CuAssertTrue(testCase, otherEdgeEnd != edgeEnd);
            }
        }
        stCactusGraphNodeIterator_destruct(nodeIt);
        CuAssertIntEquals(testCase, edgeNumber * 2, edgeEndNumber);
        //Check each chain is a simple cycle with one orientation
        nodeIt = stCactusGraphNodeIterator_construct(g2);
        while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
            stCactusNodeEdgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(node);
            stCactusEdgeEnd *edgeEnd;
            while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeEndIt)) != NULL) {
                stSortedSet *nodesOnCycle = stSortedSet_construct();
                stSortedSet_insert(nodesOnCycle, stCactusEdgeEnd_getNode(edgeEnd));
                stCactusEdgeEnd *chainEdgeEnd = edgeEnd;
                bool chainEnd = 0;
                while (1) {
                    chainEdgeEnd = stCactusEdgeEnd_getLink(stCactusEdgeEnd_getOtherEdgeEnd(chainEdgeEnd));
                    if (stCactusEdgeEnd_isChainEnd(chainEdgeEnd)) {
                        CuAssertTrue(testCase, chainEnd == 0);
                        chainEnd = 1;
                    }
                    if (chainEdgeEnd == edgeEnd) {
                        break;
                    }
                    CuAssertPtrEquals(testCase, NULL,
                            stSortedSet_search(nodesOnCycle, stCactusEdgeEnd_getNode(chainEdgeEnd)));
                    stSortedSet_insert(nodesOnCycle, stCactusEdgeEnd_getNode(chainEdgeEnd));
                }
                CuAssertTrue(testCase, chainEnd);
                stSortedSet_destruct(nodesOnCycle);
            }
        }
        stCactusGraphNodeIterator_destruct(nodeIt);

        stCactusGraph_destruct(g2);
        stSortedSet_destruct(edgeEnds);
        stList_destruct(nodeObjects);
        stSortedSet_destruct(nodesInFinalGraph);
    }
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
