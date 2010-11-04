/*
 * adjacencyComponents.c
 *
 *  Created on: 10 Sep 2010
 *      Author: benedictpaten
 */

#include "adjacencyComponents.h"

static inline void getAdjacencyComponentsP2(struct PinchVertex *vertex,
        bool (*passThroughEdgeFn)(struct PinchEdge *), bool *seen,
        stSortedSet *adjacencyComponent, stList *stack) {
    //add it to the hash
    assert(!seen[vertex->vertexID]);
    seen[vertex->vertexID] = 1;

    //add it to the component
    stSortedSet_insert(adjacencyComponent, vertex);

    //search across blackedges.
    if (lengthBlackEdges(vertex) > 0) {
        struct PinchEdge *edge = getFirstBlackEdge(vertex);
        struct PinchVertex *vertex2 = edge->to;
        if (passThroughEdgeFn(edge)) {
            if (!seen[vertex2->vertexID]) {
                stList_append(stack, vertex2);
            }
#ifdef BEN_DEBUG
            else {
                assert(stSortedSet_search(adjacencyComponent, vertex2) != NULL);
            }
#endif
        }
    }

    //search across the greyedges.
    void *greyEdgeIterator = getGreyEdgeIterator(vertex);
    struct PinchVertex *vertex2;
    while ((vertex2 = getNextGreyEdge(vertex, greyEdgeIterator)) != NULL) {
        if (!seen[vertex2->vertexID]) {
            stList_append(stack, vertex2);
        }
#ifdef BEN_DEBUG
        else {
            assert(stSortedSet_search(adjacencyComponent, vertex2) != NULL);
        }
#endif
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
    bool *seen = st_calloc(pinchGraph->vertices->length, sizeof(bool));
    stList *adjacencyComponents = stList_construct3(0, (void (*)(void *))stSortedSet_destruct);
    stList *stack = stList_construct();

    for (int32_t i = 0; i < pinchGraph->vertices->length; i++) {
        struct PinchVertex *vertex = pinchGraph->vertices->list[i];

        //if not seen
        if (!seen[vertex->vertexID]) {
            //get component
            stSortedSet *adjacencyComponent = stSortedSet_construct();
            //add to final components
            stList_append(adjacencyComponents, adjacencyComponent);

            //now build the component via a (recursive) function -- now changed to have a stack which we edit to maintain the list.
            assert(stList_length(stack) == 0);
            while(1) {
                if (!seen[vertex->vertexID]) {
                    getAdjacencyComponentsP2(vertex, passThroughEdgeFn, seen, adjacencyComponent, stack);
                }
                if(stList_length(stack) == 0) {
                    break;
                }
                vertex = stList_pop(stack);
            }
        }
    }
    //clean up
    free(seen);
    stList_destruct(stack);

    return adjacencyComponents;
}

bool passThroughDegree1EdgesFn(struct PinchEdge *edge) {
    assert(lengthBlackEdges(edge->from) > 0);
    return lengthBlackEdges(edge->from) == 1 && (!isAStub(edge));
}

bool doNotPassThroughDegree1EdgesFn(struct PinchEdge *edge) {
    assert(lengthBlackEdges(edge->from) > 0);
    return 0;
}

stList *getAdjacencyComponents(struct PinchGraph *pinchGraph) {
    return getAdjacencyComponents2(pinchGraph,
        passThroughDegree1EdgesFn);
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

