/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "giantComponent.h"
#include <math.h>

static stList *nodes = NULL;
static stList *edges;
static int32_t maxComponentSize;

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
    nodes = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    int32_t nodeNumber = st_randomInt(0, 1000);
    for(int32_t i=0; i<nodeNumber; i++) {
        stList_append(nodes, stIntTuple_construct(1, i));
    }

    //Make edges
    edges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    float edgeProb = st_random();
    for(int32_t i=0; i<nodeNumber; i++) {
        for(int32_t j=i; j<nodeNumber; j++) {
            if(st_random() <= edgeProb) {
                stList_append(edges, stIntTuple_construct(3, st_randomInt(1, 100), i, j));
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
    stHash *nodesToComponents = stHash_construct3((uint32_t (*)(const void *))stIntTuple_hashKey, (int (*)(const void *, const void *))stIntTuple_equalsFn, NULL, NULL);
    for(int32_t i=0; i<stList_length(nodes); i++) {
        stIntTuple *node = stList_get(nodes, i);
        stSortedSet *component = stSortedSet_construct();
        stSortedSet_insert(component, node);
        stHash_insert(nodesToComponents, node, component);
    }
    for(int32_t i=0; i<stList_length(filteredEdges); i++) {
        stIntTuple *edge = stList_get(filteredEdges, i);
        stIntTuple *node1 = stIntTuple_construct(1, stIntTuple_getPosition(edge, 1));
        stIntTuple *node2 = stIntTuple_construct(1, stIntTuple_getPosition(edge, 2));
        stSortedSet *component1 = stHash_search(nodesToComponents, node1);
        stSortedSet *component2 = stHash_search(nodesToComponents, node2);
        assert(component1 != NULL && component2 != NULL);
        if(component1 != component2) {
            stSortedSet *component3 = stSortedSet_getUnion(component1, component2);
            stSortedSetIterator *setIt = stSortedSet_getIterator(component3);
            stIntTuple *node3;
            while((node3 = stSortedSet_getNext(setIt)) != NULL) {
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
    for(int32_t i=0; i<stList_length(components); i++) {
        stSortedSet *component = stList_get(components, i);
        CuAssertTrue(testCase, stSortedSet_size(component) <= maxComponentSize);
        CuAssertTrue(testCase, stSortedSet_size(component) >= 1);
    }
    //Check no edges can be added from those filtered.
    stSortedSet *filteredEdgesSet = stList_getSortedSet(filteredEdges, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    for(int32_t i=0; i<stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        if(stSortedSet_search(filteredEdgesSet, edge) == NULL) {
            stIntTuple *node1 = stIntTuple_construct(1, stIntTuple_getPosition(edge, 1));
            stIntTuple *node2 = stIntTuple_construct(1, stIntTuple_getPosition(edge, 2));
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
    stSortedSet_setDestructor(componentsSet, (void (*)(void *))stSortedSet_destruct);
    stSortedSet_destruct(componentsSet);
    stHash_destruct(nodesToComponents);
}

static void testBreakUpComponentGreedily(CuTest *testCase) {
    for (int32_t i = 0; i < 100; i++) {
        setup();
        stList *filteredEdges = breakUpComponentGreedily(nodes, edges, maxComponentSize);
        checkComponents(testCase, filteredEdges);
        stList_destruct(filteredEdges);
        teardown();
    }
}

CuSuite* giantComponentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testBreakUpComponentGreedily);
    return suite;
}
