/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * sonLibRandomTest.c
 *
 *  Created on: 22-Jun-2010
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "CuTest.h"
#include "3_Absorb3edge2x.h"

void addEdgeToList(stList *vertices, int32_t vertex1, int32_t vertex2) {
    stList *edges = stList_get(vertices, vertex1);
    stList_append(edges, stIntTuple_construct(1, vertex2));
}

stList *getRandomGraph(int32_t vertexNumber, int32_t edgeNumber) {
    stList *vertices = stList_construct3(0, (void (*)(void *))stList_destruct);
    //Add vertices.
    for(int32_t i=0; i<vertexNumber; i++) {
        stList_append(vertices, stList_construct3(0, (void (*)(void *))stIntTuple_destruct));
    }
    //Add edges.
    while(edgeNumber-- > 0) {
        int32_t vertex1 = st_randomInt(0, vertexNumber);
        int32_t vertex2 = st_randomInt(0, vertexNumber);
        addEdgeToList(vertices, vertex1, vertex2);
        addEdgeToList(vertices, vertex2, vertex1);
    }
    return vertices;
}

static void test_3EdgeFunction(CuTest *testCase) {
    /*
     * Exercises the 3-edge function
     */
    for(int32_t test=0; test<100; test++) {

        /*
         * Get a random graph to test.
         */
        int32_t vertexNumber = st_randomInt(0, 100);
        int32_t edgeNumber = vertexNumber > 0 ? st_randomInt(0, vertexNumber * vertexNumber) : 0;
        stList *vertices = getRandomGraph(vertexNumber, edgeNumber);

        /*
         * Compute the 3-edge connected components.
         */
        stList *threeEdgeConnectedComponents = computeThreeEdgeConnectedComponents(vertices);

        /*
         * Check the output.
         */
        stSortedSet *seen = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, NULL);
        stListIterator *it = stList_getIterator(threeEdgeConnectedComponents);
        stList *threeEdgeConnectedComponent;
        while((threeEdgeConnectedComponent = stList_getNext(it)) != NULL) {
            stListIterator *it2 = stList_getIterator(threeEdgeConnectedComponent);
            stIntTuple *vertex;
            while((vertex = stList_getNext(it2)) != NULL) {
                CuAssertTrue(testCase, stSortedSet_search(seen, vertex) == NULL);
                stSortedSet_insert(seen, vertex);

                /*
                 * Now check other vertices in component are 3-edge connected to this vertex.
                 */

                /*
                 * And check other vertices are not three edge connected.
                 */
            }
            stList_destructIterator(it2);
        }
        stList_destructIterator(it);
        CuAssertTrue(testCase, stList_length(vertices) == stSortedSet_size(seen));

        /*
         * Cleanup
         */
        stList_destruct(vertices);
        stSortedSet_destruct(seen);
    }
}

CuSuite* threeEdgeTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_3EdgeFunction);
    return suite;
}
