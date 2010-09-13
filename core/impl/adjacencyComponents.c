/*
 * adjacencyComponents.c
 *
 *  Created on: 10 Sep 2010
 *      Author: benedictpaten
 */

#include "adjacencyComponents.h"

static void getAdjacencyComponentsP2(struct PinchVertex *vertex,
        bool (*passThroughEdgeFn)(struct PinchEdge *), stSortedSet *seen,
        stSortedSet *adjacencyComponent, stList *stack) {
    //add it to the hash
    assert(stSortedSet_search(seen, vertex) == NULL);
    stSortedSet_insert(seen, vertex);

    //add it to the component
    stSortedSet_insert(adjacencyComponent, vertex);

    //search across blackedges.
    if (lengthBlackEdges(vertex) > 0) {
        struct PinchEdge *edge = getFirstBlackEdge(vertex);
        struct PinchVertex *vertex2 = edge->to;
        if (passThroughEdgeFn(edge)) {
            if (stSortedSet_search(seen, vertex2) == NULL) {
                stList_append(stack, vertex2);
            }
            else {
                assert(stSortedSet_search(adjacencyComponent, vertex2) != NULL);
            }
        }
    }

    //search across the greyedges.
    void *greyEdgeIterator = getGreyEdgeIterator(vertex);
    struct PinchVertex *vertex2;
    while ((vertex2 = getNextGreyEdge(vertex, greyEdgeIterator)) != NULL) {
        if (stSortedSet_search(seen, vertex2) == NULL) {
            stList_append(stack, vertex2);
        }
        else {
            assert(stSortedSet_search(adjacencyComponent, vertex2) != NULL);
        }
    }
    destructGreyEdgeIterator(greyEdgeIterator);
}

stList *getAdjacencyComponents2(struct PinchGraph *pinchGraph,
        bool (*passThroughEdgeFn)(struct PinchEdge *)) {
    /*
     * find recursive components (each recursive component is represented as a series of vertices)
     * do as series of DFS on each connected component, not traversing long black edges.
     */
    //allocate stuff.
    stSortedSet *seen = stSortedSet_construct();
    stList *adjacencyComponents = stList_construct3(0, (void (*)(void *))stSortedSet_destruct);
    stList *stack = stList_construct();

    for (int32_t i = 0; i < pinchGraph->vertices->length; i++) {
        struct PinchVertex *vertex = pinchGraph->vertices->list[i];

        //if not seen
        if (stSortedSet_search(seen, vertex) == NULL) {
            //get component
            stSortedSet *adjacencyComponent = stSortedSet_construct();
            //add to final components
            stList_append(adjacencyComponents, adjacencyComponent);

            //now build the component via a (recursive) function -- now changed to have a stack which we edit to maintain the list.
            assert(stList_length(stack) == 0);
            stList_append(stack, vertex);
            while (stList_length(stack) > 0) {
                struct PinchVertex *vertex2 = stList_pop(stack);
                if (stSortedSet_search(seen, vertex2) == NULL) {
                    getAdjacencyComponentsP2(vertex2, passThroughEdgeFn, seen, adjacencyComponent, stack);
                }
            }
        }
    }
    //clean up
    stSortedSet_destruct(seen);
    stList_destruct(stack);

    return adjacencyComponents;
}

static bool getAdjacencyComponents_passThroughDegree1Edges(struct PinchEdge *edge) {
    assert(lengthBlackEdges(edge->from) > 0);
    return lengthBlackEdges(edge->from) == 1 && (!isAStub(edge));
}

stList *getAdjacencyComponents(struct PinchGraph *pinchGraph) {
    return getAdjacencyComponents2(pinchGraph,
        getAdjacencyComponents_passThroughDegree1Edges);
}

stHash *getVertexToAdjacencyComponentHash(struct PinchGraph *pinchGraph, stList *adjacencyComponents) {
    stHash *vertexToAdjacencyComponentHash = stHash_construct2(NULL, (void (*)(void *))stIntTuple_destruct);
    for(int32_t i=0; i<stList_length(adjacencyComponents); i++) {
        stSortedSet *adjacencyComponent = stList_get(adjacencyComponents, i);
        stSortedSetIterator *it = stSortedSet_getIterator(adjacencyComponent);
        struct PinchVertex *vertex;
        while((vertex = stSortedSet_getNext(it)) != NULL) {
            stHash_insert(vertexToAdjacencyComponentHash, vertex, stIntTuple_construct(1, i));
        }
        stSortedSet_destructIterator(it);
    }
    assert(pinchGraph->vertices->length == stHash_size(vertexToAdjacencyComponentHash));
    return vertexToAdjacencyComponentHash;
}

stList *getAdjacencyComponentGraph(struct PinchGraph *pinchGraph, stList *adjacencyComponents,
        stHash *vertexToAdjacencyComponentsHash) {
    stList *vertices = stList_construct3(0, (void(*)(void *)) stList_destruct);
    //Write number of nodes.
    int32_t vertexCounter = 0;
    for (int32_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stSortedSet *adjacencyComponent = stList_get(adjacencyComponents, i);
        stList *edges = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
        stList_append(vertices, edges);
        stSortedSetIterator *it = stSortedSet_getIterator(adjacencyComponent);
        struct PinchVertex *vertex;
        while((vertex = stSortedSet_getNext(it)) != NULL) {
            vertexCounter++;
            //The black edges
            if (lengthBlackEdges(vertex) > 0) {
                struct PinchEdge *edge = getFirstBlackEdge(vertex);
                assert(stIntTuple_getPosition(stHash_search(vertexToAdjacencyComponentsHash, edge->from), 0) == i);
                int32_t k = stIntTuple_getPosition(stHash_search(vertexToAdjacencyComponentsHash, edge->to), 0);
                stList_append(edges, stIntTuple_construct(1, k));
            }
            else {
                assert(vertex->vertexID == 0);
            }
        }
        stSortedSet_destructIterator(it);
    }
    assert(vertexCounter == pinchGraph->vertices->length);
    return vertices;
}

bool adjacencyComponentsAreWithinNEdges(int32_t adjacencyComponent1, int32_t adjacencyComponent2, stList *adjacencyComponentGraph,
        int32_t n) {
    assert(n >= 0);
    if(adjacencyComponent1 == adjacencyComponent2) {
        return 1;
    }
    if(n > 0) {
        stSortedSet *edges = stList_get(adjacencyComponentGraph, adjacencyComponent1);
        stSortedSetIterator *it = stSortedSet_getIterator(edges);
        stIntTuple *edge;
        while((edge = stSortedSet_getNext(it)) != NULL) {
            if(adjacencyComponentsAreWithinNEdges(stIntTuple_getPosition(edge, 0), adjacencyComponent2, adjacencyComponentGraph, n-1)) {
                stSortedSet_destructIterator(it);
                return 1;
            }
        }
        stSortedSet_destructIterator(it);
    }
    return 0;
}
