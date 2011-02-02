/*
 * orientation.c
 *
 * Orients chain with respect to one another.
 *
 *  Created on: 19 Oct 2010
 *      Author: benedictpaten
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "pinchGraph.h"
#include "cactusGraph.h"
#include "cactus.h"
#include "cactusFlowerFunctions.h"
#include "avl.h"
#include "sonLib.h"
#include "adjacencyComponents.h"

static End *getStubEnd(struct PinchVertex *vertex, Flower *flower, struct hashtable *endNamesHash) {
    /*
     * Gets the inherited end associated with the vertex.
     */
    if (!vertex_isEnd(vertex)) {
        return NULL;
    }
    char *nameString = hashtable_search(endNamesHash, vertex);
#ifdef BEN_DEBUG
    assert(nameString != NULL);
#endif
    End *end = flower_getEnd(flower, cactusMisc_stringToName(nameString));
#ifdef BEN_DEBUG
    assert(end != NULL);
    assert(end_getOrientation(end));
#endif
    return end;
}

void getAdjacentVerticesP(struct PinchVertex *vertex, stSortedSet *includedVertices, struct PinchGraph *pinchGraph,
        Flower *flower, struct hashtable *endNamesHash, stHash *orientationHash, stSortedSet *adjacentVertices) {
    void *blackEdgeIterator = getBlackEdgeIterator(vertex);
    struct PinchEdge *edge;
    while ((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
#ifdef BEN_DEBUG
        assert(edge->rEdge->to == vertex);
#endif
        struct PinchEdge *edge2 = getNextEdge(pinchGraph, edge->rEdge, flower);
#ifdef BEN_DEBUG
        assert(edge2 != NULL);
#endif
        while (stSortedSet_search(includedVertices, edge2->from) == NULL && stHash_search(orientationHash, edge2->from)
                == NULL) {
#ifdef BEN_DEBUG
            assert(stHash_search(orientationHash, edge2->to) == NULL);
#endif
            if (vertex_isEnd(edge2->from)) {
#ifdef BEN_DEBUG
                assert(vertex_isDeadEnd(edge2->to));
                assert(stSortedSet_search(includedVertices, edge2->to) != NULL);
#endif
                break;
            }
#ifdef BEN_DEBUG
            assert(stSortedSet_search(includedVertices, edge2->to) == NULL);
#endif
            edge2 = getNextEdge(pinchGraph, edge2, flower);
            assert(edge2 != NULL);
        }
        if (vertex_isEnd(edge2->from)) {
            assert(vertex_isDeadEnd(edge2->to));
            End *end = getStubEnd(edge2->from, flower, endNamesHash);
#ifdef BEN_DEBUG
            assert(end != NULL);
            assert(!end_isBlockEnd(end));
#endif
            if (end_isAttached(end)) { //We choose to ignore free stubs because they will not form chains
                stSortedSet_insert(adjacentVertices, edge2->from);
            }
        } else {
            stSortedSet_insert(adjacentVertices, edge2->from);
        }
    }
    destructBlackEdgeIterator(blackEdgeIterator);
}

static stSortedSet *getAdjacentVertices(struct PinchVertex *vertex, stSortedSet *includedVertices,
        struct PinchGraph *pinchGraph, Flower *flower, struct hashtable *endNamesHash, stHash *orientationHash) {
    /*
     * Gets a list of the vertices that are 'adjacent' to this vertex. That is vertices that reachable by a sequence from the vertex, or from the
     * vertex and one of its directly adjacent vertices.
     */
    stSortedSet *adjacentVertices = stSortedSet_construct();
    if (vertex_isDeadEnd(vertex)) {
        return adjacentVertices;
    }

    getAdjacentVerticesP(vertex, includedVertices, pinchGraph, flower, endNamesHash, orientationHash, adjacentVertices);
    //No need to consider itself again.
    if (stSortedSet_search(adjacentVertices, vertex) != NULL) {
        stSortedSet_remove(adjacentVertices, vertex);
    }

    stList *list = stSortedSet_getList(adjacentVertices);
    while (stList_length(list) > 0) {
        getAdjacentVerticesP(stList_pop(list), includedVertices, pinchGraph, flower, endNamesHash, orientationHash,
                adjacentVertices);
    }
    stList_destruct(list);

    if (stSortedSet_search(adjacentVertices, vertex) != NULL) {
        stSortedSet_remove(adjacentVertices, vertex);
    }

    return adjacentVertices;
}

static stHash *getAdjacentVerticesHash(stSortedSet *endVertices, struct PinchGraph *pinchGraph, Flower *flower,
        struct hashtable *endNamesHash, stHash *orientationHash) {
    //Get the adjacent vertices to each vertex
    stHash *adjacentVertices = stHash_construct2(NULL, (void(*)(void *)) stSortedSet_destruct);
    stSortedSetIterator *it = stSortedSet_getIterator(endVertices);
    struct PinchVertex *vertex;
    while ((vertex = stSortedSet_getNext(it)) != NULL) {
        stHash_insert(adjacentVertices, vertex, getAdjacentVertices(vertex, endVertices, pinchGraph, flower,
                endNamesHash, orientationHash));
    }
    stSortedSet_destructIterator(it);
    return adjacentVertices;
}

static stList *getPetals(struct List *biConnectedComponents, int32_t *i) {
    /*
     * Gets the next several chains from the list with the same source vertex, returns NULL when the
     * index i is at the end of the list.
     */
    if (*i >= biConnectedComponents->length) {
        return NULL;
    }
    stList *petals = stList_construct();
    struct List *biConnectedComponent = biConnectedComponents->list[*i];
    stList_append(petals, biConnectedComponent);
#ifdef BEN_DEBUG
    assert(biConnectedComponent->length > 0);
#endif
    struct CactusVertex *vertex = ((struct CactusEdge *) biConnectedComponent->list[0])->from;
#ifdef BEN_DEBUG
    assert(vertex == ((struct CactusEdge *)biConnectedComponent->list[biConnectedComponent->length-1])->to);
#endif
    while (++(*i) < biConnectedComponents->length) {
        biConnectedComponent = biConnectedComponents->list[(*i)];
#ifdef BEN_DEBUG
        assert(biConnectedComponent->length > 0);
#endif
        struct CactusVertex *vertex2 = ((struct CactusEdge *) biConnectedComponent->list[0])->from;
#ifdef BEN_DEBUG
        assert(vertex2 == ((struct CactusEdge *)biConnectedComponent->list[biConnectedComponent->length-1])->to);
#endif
        if (vertex != vertex2) {
            break;
        }
        stList_append(petals, biConnectedComponent);
    }
    return petals;
}

static stList *getEndVertices(stList *chains, struct PinchGraph *pinchGraph) {
    /*
     * Gets the vertices in the pinch graph representing the ends of the given set of chains.
     */
#ifdef BEN_DEBUG
    assert(stList_length(chains) > 0);
#endif
    stList *endVertices = stList_construct();
    for (int32_t i = 0; i < stList_length(chains); i++) {
        struct List *biConnectedComponent = stList_get(chains, i);
#ifdef BEN_DEBUG
        assert(biConnectedComponent->length > 0);
#endif
        stList_append(endVertices, cactusEdgeToFirstPinchEdge(biConnectedComponent->list[0], pinchGraph)->from);
        stList_append(endVertices, cactusEdgeToFirstPinchEdge(biConnectedComponent->list[biConnectedComponent->length
                - 1], pinchGraph)->to);
    }
    return endVertices;
}

static void orientP(struct PinchVertex *vertex1, struct PinchVertex *vertex2, stHash *orientationHash, bool orientation) {
#ifdef BEN_DEBUG
    assert(stHash_search(orientationHash, vertex1) == NULL);
    assert(stHash_search(orientationHash, vertex2) == NULL);
#endif
    stHash_insert(orientationHash, vertex1, stIntTuple_construct(1, orientation));
    stHash_insert(orientationHash, vertex2, stIntTuple_construct(1, !orientation));
}

static bool orient(struct PinchVertex *vertex1, struct PinchVertex *vertex2, stHash *orientationHash,
        stHash *adjacentVertices, Flower *flower, struct hashtable *endNamesHash, bool deadLock) {
    /*
     * Orients the ends.
     */
#ifdef BEN_DEBUG
    assert(stHash_search(orientationHash, vertex1) == NULL);
    assert(stHash_search(orientationHash, vertex2) == NULL);
#endif
    stSortedSet *adjacentVertices1 = stHash_search(adjacentVertices, vertex1);
    stSortedSet *adjacentVertices2 = stHash_search(adjacentVertices, vertex2);
#ifdef BEN_DEBUG
    assert(adjacentVertices1 != NULL);
    assert(adjacentVertices2 != NULL);
#endif
    if (vertex_isDeadEnd(vertex1)) {
        struct PinchVertex *vertex3 = getFirstBlackEdge(vertex1)->to;
#ifdef BEN_DEBUG
        assert(vertex_isEnd(vertex3));
#endif
        End *end = getStubEnd(vertex3, flower, endNamesHash);
#ifdef BEN_DEBUG
        assert(end != NULL);
        assert(!end_isBlockEnd(end));
        if (vertex3 != vertex2) {
            assert(end_isAttached(end));
            assert(!vertex_isEnd(vertex2));
            if (vertex_isDeadEnd(vertex2)) {
                assert(lengthBlackEdges(vertex2) > 0);
                struct PinchVertex *vertex4 = getFirstBlackEdge(vertex2)->to;
                assert(vertex_isEnd(vertex4));
                End *end2 = getStubEnd(vertex4, flower, endNamesHash);
                assert(end_isAttached(end2));
                assert(end_getSide(end) == !end_getSide(end2));
            } else {
                if (stSortedSet_size(adjacentVertices2) == 1) { //Is in a pseudo chain
                    stIntTuple *i = stHash_search(orientationHash, stSortedSet_getFirst(
                                    adjacentVertices2));
                    if (i != NULL) {
                        assert(!end_getSide(end) == stIntTuple_getPosition(i, 0));
                    }
                }
            }
        }
#endif
        orientP(vertex1, vertex2, orientationHash, !end_getSide(end));
        return 1;
    }
    if (vertex_isEnd(vertex1)) {
#ifdef BEN_DEBUG
        assert(vertex_isDeadEnd(vertex2));
#endif
        End *end = getStubEnd(vertex1, flower, endNamesHash);
#ifdef BEN_DEBUG
        assert(end != NULL);
        assert(!end_isBlockEnd(end));
#endif
        orientP(vertex1, vertex2, orientationHash, end_getSide(end));
        return 1;
    }
#ifdef BEN_DEBUG
    assert(!vertex_isEnd(vertex2));
#endif
    if (stSortedSet_size(adjacentVertices1) == 1) { //Is in a pseudo chain
        stIntTuple *i = stHash_search(orientationHash, stSortedSet_getFirst(adjacentVertices1));
        if (i != NULL) {
#ifdef BEN_DEBUG
            if(vertex_isDeadEnd(vertex2)) {
                struct PinchVertex *vertex3 = getFirstBlackEdge(vertex2)->to;
                assert(vertex_isEnd(vertex3));
                End *end = getStubEnd(vertex3, flower, endNamesHash);
                assert(end != NULL);
                assert(end_getSide(end) == !stIntTuple_getPosition(i, 0));
            } else if(stSortedSet_size(adjacentVertices2) == 1) {
                stIntTuple *j = stHash_search(orientationHash, stSortedSet_getFirst(
                                adjacentVertices2));
                if(j != NULL) {
                    assert(stIntTuple_getPosition(i, 0) == !stIntTuple_getPosition(j, 0));
                }
            }
#endif
            orientP(vertex1, vertex2, orientationHash, !stIntTuple_getPosition(i, 0));
            return 1;
        }
    } else if ((deadLock || stSortedSet_size(adjacentVertices2) != 1) && !vertex_isDeadEnd(vertex2)) { //Does not have an end or a pseudo chain at either end, so we are free to choose
        orientP(vertex1, vertex2, orientationHash, 1);
        return 1;
    }
    return 0;
}

void propagateOrientation(struct List *chain, stHash *orientationHash, struct PinchGraph *pinchGraph) {
    assert(chain->length > 0);
    stIntTuple *i = stHash_search(orientationHash, cactusEdgeToFirstPinchEdge(((struct CactusEdge *) chain->list[0]),
            pinchGraph)->from);
    assert(i != NULL);
    bool orientation = stIntTuple_getPosition(i, 0);
#ifdef BEN_DEBUG
    i = stHash_search(orientationHash,
            cactusEdgeToFirstPinchEdge(
                    ((struct CactusEdge *) chain->list[chain->length - 1]),
                    pinchGraph)->to);
    assert(i != NULL);
    assert(stIntTuple_getPosition(i, 0) == !orientation);
#endif
    for (int32_t j = 1; j < chain->length; j++) {
        struct PinchVertex *vertex =
                cactusEdgeToFirstPinchEdge(((struct CactusEdge *) chain->list[j]), pinchGraph)->from;
        assert(stHash_search(orientationHash, vertex) == NULL);
        stHash_insert(orientationHash, vertex, stIntTuple_construct(1, orientation));
    }
    for (int32_t j = 0; j < chain->length - 1; j++) {
        struct PinchVertex *vertex = cactusEdgeToFirstPinchEdge(((struct CactusEdge *) chain->list[j]), pinchGraph)->to;
        assert(stHash_search(orientationHash, vertex) == NULL);
        stHash_insert(orientationHash, vertex, stIntTuple_construct(1, !orientation));
    }
#ifdef BEN_DEBUG
    for(int32_t j=0; j<chain->length; j++) {
        struct PinchEdge *edge = cactusEdgeToFirstPinchEdge(
                ((struct CactusEdge *) chain->list[j]), pinchGraph);
        assert(stHash_search(orientationHash, edge->from) != NULL);
        assert(stHash_search(orientationHash, edge->to) != NULL);
    }
#endif
}

stHash *buildOrientationHash(struct List *biConnectedComponents, struct PinchGraph *pinchGraph, Flower *flower,
        struct hashtable *endNamesHash) {
    /*
     * Assigns an orientation to every pinch graph vertex that will be in the final cactus.
     */
    stHash *orientationHash = stHash_construct2(NULL, (void(*)(void *)) stIntTuple_destruct);
    int32_t i = 0;
    stList *petals = NULL;
    while ((petals = getPetals(biConnectedComponents, &i)) != NULL) {
        stList *endVertices = getEndVertices(petals, pinchGraph);
        stSortedSet *endVerticesSet = stList_getSortedSet(endVertices, NULL);
#ifdef BEN_DEBUG
        assert(stList_length(endVertices) > 1); //Must be at least 2 ends.
        assert(stList_length(endVertices) % 2 == 0); //Must be equal number.
        assert(stList_length(endVertices) == stSortedSet_size(endVerticesSet));
        int32_t endVerticesSize = stSortedSet_size(endVerticesSet);
#endif
        stHash *adjacentVerticesHash = getAdjacentVerticesHash(endVerticesSet, pinchGraph, flower, endNamesHash,
                orientationHash);
        //Now fill in the orientations
        bool deadLock = 0;
        while (stList_length(endVertices) > 0) {
            stList *endVertices2 = stList_construct();
            for (int32_t j = 0; j < stList_length(endVertices); j += 2) {
                struct PinchVertex *vertex1 = stList_get(endVertices, j);
                struct PinchVertex *vertex2 = stList_get(endVertices, j + 1);
                if (!orient(vertex1, vertex2, orientationHash, adjacentVerticesHash, flower, endNamesHash, deadLock) && !orient(
                        vertex2, vertex1, orientationHash, adjacentVerticesHash, flower, endNamesHash, deadLock)) {
                    stList_append(endVertices2, vertex1);
                    stList_append(endVertices2, vertex2);
                }
                else {
                    deadLock = 0;
                }
            }
            deadLock = stList_length(endVertices) == stList_length(endVertices2);
            stList_destruct(endVertices);
            endVertices = endVertices2;
        }
#ifdef BEN_DEBUG
        assert(endVerticesSize == stSortedSet_size(endVerticesSet));
        stSortedSetIterator *it = stSortedSet_getIterator(endVerticesSet);
        struct PinchVertex *vertex;
        while((vertex = stSortedSet_getNext(it)) != NULL) {
            stIntTuple *j = stHash_search(orientationHash, vertex);
            assert(j != NULL);
            bool isFreeStub = 0;
            if(vertex_isEnd(vertex)) {
                End *end = getStubEnd(vertex, flower, endNamesHash);
                assert(end != NULL);
                assert(end_getSide(end) == stIntTuple_getPosition(j, 0));
                assert(!end_isBlockEnd(end));
                isFreeStub = end_isFree(end);
            }
            stSortedSet *adjacentVertices = stHash_search(adjacentVerticesHash, vertex);
            assert(adjacentVertices != NULL);
            if(!vertex_isDeadEnd(vertex) && stSortedSet_size(adjacentVertices) == 1 && !isFreeStub) {
                stIntTuple *k = stHash_search(orientationHash, stSortedSet_getFirst(adjacentVertices));
                assert(k != NULL);
                assert(stIntTuple_getPosition(j, 0) == !stIntTuple_getPosition(k, 0));
            }
        }
        stSortedSet_destructIterator(it);
#endif
        //Now finally propagate the orientation information through the chains
        for (int32_t j = 0; j < stList_length(petals); j++) {
            struct List *chain = stList_get(petals, j);
            propagateOrientation(chain, orientationHash, pinchGraph);
        }
        stSortedSet_destruct(endVerticesSet);
        stList_destruct(endVertices);
        stHash_destruct(adjacentVerticesHash);
        stList_destruct(petals);
    }
#ifdef BEN_DEBUG
    int32_t l = 0;
    for(int32_t j=0; j<biConnectedComponents->length; j++) {
        struct List *biConnectedComponent = biConnectedComponents->list[j];
        for(int32_t k=0; k<biConnectedComponent->length; k++) {
            struct PinchEdge *edge = cactusEdgeToFirstPinchEdge(biConnectedComponent->list[k], pinchGraph);
            assert(stHash_search(orientationHash, edge->from) != NULL);
            assert(stHash_search(orientationHash, edge->to) != NULL);
            l+=2;
        }
    }
    assert(l == stHash_size(orientationHash));
#endif
    return orientationHash;
}
