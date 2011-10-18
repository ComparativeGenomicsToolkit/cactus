/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
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
#include "3_Absorb3edge2x.h"
#include "adjacencyComponents.h"
#include "pinchGraphManipulation.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Cactus graph data structures for representing the
//basic cactus graph.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct CactusVertex *constructCactusVertex() {
    struct CactusVertex *cactusVertex;

    cactusVertex = st_malloc(sizeof(struct CactusVertex));
    cactusVertex->edges = constructEmptyList(0, (void(*)(void *)) destructCactusEdge);

    return cactusVertex;
}

void destructCactusVertex(struct CactusVertex *vertex) {
    destructList(vertex->edges);
    free(vertex);
}

struct CactusEdge *constructCactusEdge(struct List *pieces) {
    struct CactusEdge *edge;
    struct CactusEdge *rEdge;
    struct Piece *piece;
    int32_t i;

    edge = st_malloc(sizeof(struct CactusEdge));
    rEdge = st_malloc(sizeof(struct CactusEdge));
    edge->rEdge = rEdge;
    rEdge->rEdge = edge;

    edge->pieces = constructEmptyList(pieces->length, NULL);
    edge->rEdge->pieces = constructEmptyList(pieces->length, NULL);

    for (i = 0; i < pieces->length; i++) {
        piece = pieces->list[i];
        edge->pieces->list[i] = piece;
        edge->rEdge->pieces->list[i] = piece->rPiece;
    }

    return edge;
}

struct CactusEdge *constructCactusEdge2(struct List *pieces, struct CactusVertex *from, struct CactusVertex *to) {
    struct CactusEdge *cactusEdge = constructCactusEdge(pieces);

    listAppend(from->edges, cactusEdge);
    cactusEdge->from = from;
    cactusEdge->rEdge->to = from;

    listAppend(to->edges, cactusEdge->rEdge);
    cactusEdge->to = to;
    cactusEdge->rEdge->from = to;

    return cactusEdge;
}

void destructCactusEdge(struct CactusEdge *edge) {
    //the rEdge is dealt with by its respective from vertex.
    destructList(edge->pieces);
    free(edge);
}

struct CactusGraph *constructCactusGraph(struct PinchGraph *pinchGraph, struct List *threeEdgeConnectedComponents,
        bool(*passThroughEdgeFn)(struct PinchEdge *)) {
    struct CactusGraph *cactusGraph;
    struct CactusEdge *cactusEdge;
    struct CactusVertex *cactusVertex;
    struct CactusVertex *cactusVertex2;
    int32_t i, j;
    struct List *list;
    struct PinchVertex *pinchVertex;
    struct PinchVertex *pinchVertex2;
    struct PinchEdge *pinchEdge;
    struct List *pinchVertexToCactusVertex;
    struct List *emptyList;
    struct List *list2;

    cactusGraph = st_malloc(sizeof(struct CactusGraph));
    cactusGraph->vertices = constructEmptyList(0, (void(*)(void *)) destructCactusVertex);

#ifdef BEN_DEBUG
    assert(threeEdgeConnectedComponents->length > 0);
    list = threeEdgeConnectedComponents->list[0];
    j = FALSE;
    for (i = 0; i < list->length; i++) {
        pinchVertex = list->list[i];
        if (pinchVertex->vertexID == 0) {
            j = TRUE;
        }
    }
    assert(j == TRUE);
#endif

    pinchVertexToCactusVertex = constructEmptyList(pinchGraph->vertices->length, NULL);
    assert(pinchVertexToCactusVertex->length == pinchGraph->vertices->length);
    for (i = 0; i < pinchVertexToCactusVertex->length; i++) {
        pinchVertexToCactusVertex->list[i] = NULL;
    }

    for (i = 0; i < threeEdgeConnectedComponents->length; i++) {
        cactusVertex = constructCactusVertex();
        cactusVertex->vertexID = i;
        listAppend(cactusGraph->vertices, cactusVertex);

        list = threeEdgeConnectedComponents->list[i];
#ifdef BEN_DEBUG
        //checks all components are non empty.
        assert(list->length> 0);
#endif

        for (j = 0; j < list->length; j++) {
            pinchVertex = list->list[j];
#ifdef BEN_DEBUG
            //checks that there is only one reference to the pinch vertex in the list of vertices.
            assert(pinchVertexToCactusVertex->list[pinchVertex->vertexID] == NULL);
#endif
            pinchVertexToCactusVertex->list[pinchVertex->vertexID] = cactusVertex;
        }
    }

#ifdef BEN_DEBUG
    //checks all vertices in pinch graph are represented in cactus vertex.
    for (i = 0; i < pinchVertexToCactusVertex->length; i++) {
        assert(pinchVertexToCactusVertex->list[i] != NULL);
    }
#endif

    emptyList = constructEmptyList(0, NULL);
    for (i = 0; i < threeEdgeConnectedComponents->length; i++) {
        list = threeEdgeConnectedComponents->list[i];
        cactusVertex = cactusGraph->vertices->list[i];

        for (j = 0; j < list->length; j++) {
            pinchVertex = list->list[j];

            //black edges
            if (lengthBlackEdges(pinchVertex) > 0) {
                pinchEdge = getFirstBlackEdge(pinchVertex);
                if (!passThroughEdgeFn(pinchEdge)) {
                    pinchVertex2 = pinchEdge->to;
                    cactusVertex2 = pinchVertexToCactusVertex->list[pinchVertex2->vertexID];

                    if (pinchVertex->vertexID < pinchVertex2->vertexID) {
                        //getting the ordered pieces.
                        list2 = constructEmptyList(0, NULL);
                        void *blackEdgeIterator = getBlackEdgeIterator(pinchVertex);
                        while ((pinchEdge = getNextBlackEdge(pinchVertex, blackEdgeIterator)) != NULL) {
                            listAppend(list2, pinchEdge->piece);
                        }
                        destructBlackEdgeIterator(blackEdgeIterator);

                        cactusEdge = constructCactusEdge2(list2, cactusVertex, cactusVertex2);
                        destructList(list2);
                    }
                }
            }
#ifdef BEN_DEBUG
            else {
                assert(pinchVertex->vertexID == 0);
            }
#endif

#ifdef BEN_DEBUG
            //grey edges
            void *greyEdgeIterator = getGreyEdgeIterator(pinchVertex);
            while ((pinchVertex2 = getNextGreyEdge(pinchVertex,
                                    greyEdgeIterator)) != NULL) {
                cactusVertex2
                = pinchVertexToCactusVertex->list[pinchVertex2->vertexID];

                if (cactusVertex != cactusVertex2 && cactusVertex
                        < cactusVertex2) {
                    assert(0); //We always merge adjacency components, so this should never happen!
                    //cactusEdge = constructCactusEdge2(emptyList, cactusVertex,
                    //        cactusVertex2);
                }
            }
            destructGreyEdgeIterator(greyEdgeIterator);
#endif
        }
    }

    //cleanup
    destructList(emptyList);
    destructList(pinchVertexToCactusVertex);

    return cactusGraph;
}

void destructCactusGraph(struct CactusGraph *cactusGraph) {
    destructList(cactusGraph->vertices);
    free(cactusGraph);
}

int32_t cactusGraph_getEdgeNumber(struct CactusGraph *cactusGraph) {
    stSortedSet *sortedSet = stSortedSet_construct();
    for (int32_t i = 0; i < cactusGraph->vertices->length; i++) {
        struct CactusVertex *cactusVertex = cactusGraph->vertices->list[i];
        for (int32_t j = 0; j < cactusVertex->edges->length; j++) {
            struct CactusEdge *cactusEdge = cactusVertex->edges->list[j];
            stSortedSet_insert(sortedSet, cactusEdge);
            stSortedSet_insert(sortedSet, cactusEdge->rEdge);
        }
    }
    int32_t k = stSortedSet_size(sortedSet);
    stSortedSet_destruct(sortedSet);
    assert(k % 2 == 0);
    return k / 2;
}

void checkCactusGraph(struct PinchGraph *pinchGraph, struct List *threeEdgeConnectedComponents,
        struct CactusGraph *cactusGraph) {
    assert(pinchGraph != NULL);
    assert(threeEdgeConnectedComponents != NULL);
    assert(cactusGraph != NULL);
#ifdef BEN_ULTRA_DEBUG
    struct hashtable *hashTable;
    int32_t i, j, k, l, m;
    struct CactusVertex *cactusVertex;
    struct CactusVertex *cactusVertex2;
    struct CactusEdge *cactusEdge;
    struct CactusEdge *cactusEdge2;
    struct PinchVertex *pinchVertex;
    struct PinchVertex *pinchVertex2;
    struct PinchEdge *pinchEdge;
    struct Piece *piece;
    struct List *list;
    struct List *biConnectedComponents;
    /*
     * Method to check that pinch graph is okay.
     */

    hashTable = create_hashtable(cactusGraph->vertices->length*10,
            hashtable_key, hashtable_equalKey,
            NULL, NULL);

    //does general consistency checks on graph.
    assert(threeEdgeConnectedComponents->length == cactusGraph->vertices->length);
    for(i=0; i<cactusGraph->vertices->length; i++) {
        cactusVertex = cactusGraph->vertices->list[i];
        assert(cactusVertex->vertexID == i);

        list = threeEdgeConnectedComponents->list[i];
        assert(list->length> 0);

        for(j=0; j<list->length; j++) {
            hashtable_insert(hashTable, list->list[j], cactusVertex);
        }

        for(k=0; k<cactusVertex->edges->length; k++) {
            cactusEdge = cactusVertex->edges->list[k];
            assert(cactusEdge->from == cactusVertex);
            assert(cactusEdge->rEdge->to == cactusVertex);
        }
    }

    //check 0 vertex is special in cactus graph..
    pinchVertex = pinchGraph->vertices->list[0];
    //check vertex is in right component.
    cactusVertex = hashtable_search(hashTable, pinchVertex);
    assert(cactusVertex->vertexID == 0);

    //checks every edge in pinch graph is represented, and connected to its appropriate
    //3 edge connected component vertex.
    for(i=0; i<pinchGraph->vertices->length; i++) {
        pinchVertex = pinchGraph->vertices->list[i];
        //check vertex is in right component.
        cactusVertex = hashtable_search(hashTable, pinchVertex);
        assert(cactusVertex != NULL);
        list = threeEdgeConnectedComponents->list[cactusVertex->vertexID];
        k = FALSE;
        for(j=0; j<list->length; j++) {
            if(list->list[j] == pinchVertex) {
                k = TRUE;
            }
        }
        assert(k == TRUE);

        //black edges.
        void *blackEdgeIterator = getBlackEdgeIterator(pinchVertex);
        while((pinchEdge = getNextBlackEdge(pinchVertex, blackEdgeIterator)) != NULL) {
            cactusVertex2 = hashtable_search(hashTable, pinchEdge->to);
            assert(cactusVertex2 != NULL);

            m = FALSE;
            for(k=0; k<cactusVertex->edges->length; k++) {
                cactusEdge = cactusVertex->edges->list[k];
                assert(cactusEdge->from == cactusVertex);
                if(cactusEdge->to == cactusVertex2) {
                    for(l=0; l<cactusEdge->pieces->length; l++) {
                        piece = cactusEdge->pieces->list[l];
                        if(pieceComparator(pinchEdge->piece, piece) == 0) {
                            assert(m == FALSE);
                            m = TRUE;
                        }
                    }
                }
            }
            assert(m == TRUE);
        }
        destructBlackEdgeIterator(blackEdgeIterator);

        //grey edges.
        void *greyEdgeIterator = getGreyEdgeIterator(pinchVertex);
        while((pinchVertex2 = getNextGreyEdge(pinchVertex, greyEdgeIterator)) != NULL) {
            cactusVertex2 = hashtable_search(hashTable, pinchVertex2);
            assert(cactusVertex2 != NULL);

            m = 0;
            for(k=0; k<cactusVertex->edges->length; k++) {
                cactusEdge = cactusVertex->edges->list[k];
                assert(cactusEdge->from == cactusVertex);
                if(cactusEdge->to == cactusVertex2) {
                    if(cactusEdge->pieces->length == 0) {
                        m += 1;
                    }
                }
            }
            if(cactusVertex != cactusVertex2) {
                assert(m> 0);
            }
            else {
                assert(m == 0);
            }
        }
        destructGreyEdgeIterator(greyEdgeIterator);

    }
    hashtable_destroy(hashTable, FALSE, FALSE);

    //check bi connected ness of the graph!
    biConnectedComponents = computeBiConnectedComponents(cactusGraph);

    hashTable = create_hashtable(cactusGraph->vertices->length*10,
            hashtable_key, hashtable_equalKey,
            NULL, NULL);

    //check each edge is in exactly one biConnectedComponent
    for(i=0; i<biConnectedComponents->length; i++) {
        list = biConnectedComponents->list[i];
        for(j=0; j<list->length; j++) {
            cactusEdge = list->list[j];
            assert(hashtable_search(hashTable, cactusEdge) == NULL);
            assert(hashtable_search(hashTable, cactusEdge->rEdge) == NULL);
            hashtable_insert(hashTable, cactusEdge, cactusEdge);
        }
    }
    for(i=0; i<cactusGraph->vertices->length; i++) {
        cactusVertex = cactusGraph->vertices->list[i];
        for(j=0; j<cactusVertex->edges->length; j++) {
            cactusEdge = cactusVertex->edges->list[j];
            assert(hashtable_search(hashTable, cactusEdge) != NULL || hashtable_search(hashTable, cactusEdge->rEdge) != NULL);
        }
    }
    hashtable_destroy(hashTable, FALSE, FALSE);

    hashTable = create_hashtable(cactusGraph->vertices->length*10,
            hashtable_key, hashtable_equalKey,
            NULL, NULL);

    //check every biconnected component is either a single edge or a cycle.
    for(i=0; i<biConnectedComponents->length; i++) {
        list = biConnectedComponents->list[i];

        if(list->length> 1) {
            //now traverse list;
            cactusEdge = list->list[0];

            for(j=1; j<list->length; j++) {
                assert(hashtable_search(hashTable, cactusEdge) == NULL);
                hashtable_insert(hashTable, cactusEdge, cactusEdge);

                for(k=0; k<list->length; k++) {
                    cactusEdge2 = list->list[k];
                    if(cactusEdge2->from == cactusEdge->to) {
                        break;
                    }
                }
                assert(cactusEdge2->from == cactusEdge->to);
                cactusEdge = cactusEdge2;
            }
            assert(hashtable_search(hashTable, cactusEdge) == NULL);
            hashtable_insert(hashTable, cactusEdge, cactusEdge);

            cactusEdge2 = list->list[0];
            assert(cactusEdge2->from == cactusEdge->to);
        }
    }

    //Clean up all remainders.
    destructList(biConnectedComponents);
    hashtable_destroy(hashTable, FALSE, FALSE);

#endif
}

void checkCactusContainsOnly2EdgeConnectedComponents(struct CactusGraph *cactusGraph) {
    assert(cactusGraph != NULL);
#ifdef BEN_ULTRA_DEBUG
    struct CactusEdge *edge;
    struct CactusEdge *edge2;
    struct List *list;
    struct List *biConnectedComponents;
    struct hashtable *hashTable;

    int32_t i, j, k;

    ////////////////////////////////////////////////
    //(1) Get bi-connected components.
    ////////////////////////////////////////////////
    biConnectedComponents = computeBiConnectedComponents(cactusGraph);
    st_logDebug("Constructed the biconnected components for making the flower\n");

    ////////////////////////////////////////////////
    //(2) Check every component is a cycle.
    ////////////////////////////////////////////////

    hashTable = create_hashtable(cactusGraph->vertices->length*10,
            hashtable_key, hashtable_equalKey,
            NULL, NULL);

    for(i=0; i<biConnectedComponents->length; i++) {
        list = biConnectedComponents->list[i];

        assert(list->length> 0);

        edge = list->list[0];
        for(j=1; j<list->length; j++) {
            assert(hashtable_search(hashTable, edge) == NULL);
            hashtable_insert(hashTable, edge, edge);

            for(k=0; k<list->length; k++) {
                edge2 = list->list[k];
                if(edge2->from == edge->to) {
                    break;
                }
            }
            assert(edge2->from == edge->to);
            edge = edge2;
        }
        assert(hashtable_search(hashTable, edge) == NULL);
        hashtable_insert(hashTable, edge, edge);

        //check the cycle matches up.
        edge2 = list->list[0];
        assert(edge2->from == edge->to);
    }

    hashtable_destroy(hashTable, FALSE, FALSE);
    destructList(biConnectedComponents);
    st_logDebug("Checked that all edges in the cactus graph are in a single simple cycle.\n");
#endif
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method to calculate 2-edge connected components
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void computeBiConnectedComponents_P(struct hashtable *flag, int32_t count, int32_t *dFN, int32_t *low,
        struct List *stack, int32_t *father, struct List *biConnnectedComponents, struct CactusVertex *v) {
    int32_t i;
    struct CactusVertex *w;
    struct CactusEdge *edge;
    struct List *list;

    hashtable_insert(flag, v, v);
    count++;
    dFN[v->vertexID] = count;
    low[v->vertexID] = count;

    for (i = 0; i < v->edges->length; i++) {
        edge = v->edges->list[i];
        w = edge->to;

        assert(hashtable_search(flag, edge) == NULL);

        //self edges are special case
        if (w == v) {
            if (hashtable_search(flag, edge->rEdge) == NULL) {
                hashtable_insert(flag, edge, edge);
                list = constructEmptyList(0, NULL);
                listAppend(list, edge);
                listAppend(biConnnectedComponents, list);
            }
            continue;
        }

        if (hashtable_search(flag, edge->rEdge) == NULL) {
            hashtable_insert(flag, edge, edge);
            listAppend(stack, edge);
        }

        if (hashtable_search(flag, w) == NULL) {
            father[w->vertexID] = v->vertexID;
            computeBiConnectedComponents_P(flag, count, dFN, low, stack, father, biConnnectedComponents, w);
            if (low[w->vertexID] >= dFN[v->vertexID]) { //is articulation point
                list = constructEmptyList(0, NULL);
                while (stack->length > 0) {
                    listAppend(list, stack->list[--stack->length]);
                    if (list->list[list->length - 1] == edge) {
                        break;
                    }
                }
                listAppend(biConnnectedComponents, list);
            }
            low[v->vertexID] = low[v->vertexID] < low[w->vertexID] ? low[v->vertexID] : low[w->vertexID];
        } else {
            if (w->vertexID != father[v->vertexID]) {
                low[v->vertexID] = low[v->vertexID] < dFN[w->vertexID] ? low[v->vertexID] : dFN[w->vertexID];
            }
        }
    }
}

struct List *computeBiConnectedComponents(struct CactusGraph *cactusGraph) {
    /*
     * Computes the set of bi-connected components, as lists of edges, for the cactus graph.
     * Bridge edges are components containing single edges.
     * See Dan Hirschberg's explanation..
     * See http://www.ics.uci.edu/~dan/class/161/notes/8/Bicomps.html for an explantion
     */
    struct hashtable *flag;
    int32_t count;
    int32_t *dFN;
    int32_t *low;
    struct List *stack;
    int32_t *father;
    struct List *biConnectedComponents;
    int32_t i;

    flag = create_hashtable(cactusGraph->vertices->length * 10, hashtable_key, hashtable_equalKey, NULL, NULL);
    count = 0;
    dFN = st_malloc(sizeof(int32_t) * cactusGraph->vertices->length);
    low = st_malloc(sizeof(int32_t) * cactusGraph->vertices->length);
    father = st_malloc(sizeof(int32_t) * cactusGraph->vertices->length);
    for (i = 0; i < cactusGraph->vertices->length; i++) {
        dFN[i] = -1;
        low[i] = -1;
        father[i] = -1;
    }
    stack = constructEmptyList(0, NULL);
    biConnectedComponents = constructEmptyList(0, (void(*)(void *)) destructList);

    computeBiConnectedComponents_P(flag, count, dFN, low, stack, father, biConnectedComponents,
            cactusGraph->vertices->list[0]);

    hashtable_destroy(flag, FALSE, FALSE);
    free(dFN);
    free(low);
    free(father);
    destructList(stack);
    return biConnectedComponents;
}

int32_t *sortBiConnectedComponents_vertexOrdering;

int sortBiConnectedComponentsP(struct CactusEdge **edge, struct CactusEdge **edge2) {
    return sortBiConnectedComponents_vertexOrdering[(*edge)->from->vertexID]
            - sortBiConnectedComponents_vertexOrdering[(*edge2)->from->vertexID];
}

struct List *computeSortedBiConnectedComponents(struct CactusGraph *cactusGraph) {
    struct List *biConnectedComponents;
    struct List *biConnectedComponent;
    int32_t i;

#ifdef BEN_ULTRA_DEBUG
    struct CactusEdge *edge;
    struct CactusEdge *edge2;
    struct CactusVertex *vertex;
    int32_t j;
#endif

    ////////////////////////////////////////////////
    //(1) Get bi-connected components.
    ////////////////////////////////////////////////
    biConnectedComponents = computeBiConnectedComponents(cactusGraph);
    st_logDebug("Constructed the biconnected components for making the flower\n");

    ////////////////////////////////////////////////
    //(1) Get DFS ordering on nodes
    ////////////////////////////////////////////////
    sortBiConnectedComponents_vertexOrdering = getDFSDiscoveryTimes(cactusGraph);
    st_logDebug("Got the DFS discovery times\n");

    ////////////////////////////////////////////////
    //(2) Now sort each bi connected component
    ////////////////////////////////////////////////
    for (i = 0; i < biConnectedComponents->length; i++) {
        biConnectedComponent = biConnectedComponents->list[i];
        qsort(biConnectedComponent->list, biConnectedComponent->length, sizeof(void *),
                (int(*)(const void *v, const void *)) sortBiConnectedComponentsP);
#ifdef BEN_ULTRA_DEBUG
        edge = biConnectedComponent->list[0];
        vertex = edge->from;
        for(j=0; j+1<biConnectedComponent->length; j++) {
            edge = biConnectedComponent->list[j];
            edge2 = biConnectedComponent->list[j+1];
            assert(edge2->from->vertexID == edge->to->vertexID);
            assert(edge->to != vertex);
            assert(edge2->from != vertex);
            assert(sortBiConnectedComponents_vertexOrdering[edge->to->vertexID]> sortBiConnectedComponents_vertexOrdering[edge->from->vertexID]);
        }
        edge = biConnectedComponent->list[biConnectedComponent->length-1];
        assert(edge->to->vertexID == vertex->vertexID);
#endif
    }
    st_logDebug("Sorted each bi-connected component\n");

    /*
     * Cleanup
     */
    free(sortBiConnectedComponents_vertexOrdering);

    return biConnectedComponents;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method to compute DFS discovery time for vertices.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t getDFSDiscoveryTimesP(struct CactusGraph *cactusGraph, struct CactusVertex *vertex, int32_t counter,
        int32_t *vertexOrdering) {
    int32_t i;
    struct CactusEdge *edge;

    assert(vertexOrdering[vertex->vertexID] == -1);
    vertexOrdering[vertex->vertexID] = counter++;
    for (i = 0; i < vertex->edges->length; i++) {
        edge = vertex->edges->list[i];
        if (vertexOrdering[edge->to->vertexID] == -1) {
            counter = getDFSDiscoveryTimesP(cactusGraph, edge->to, counter, vertexOrdering);
        }
    }
    return counter;
}

int32_t *getDFSDiscoveryTimes(struct CactusGraph *cactusGraph) {
    int32_t i;
    int32_t *vertexOrdering;

    vertexOrdering = st_malloc(sizeof(int32_t) * cactusGraph->vertices->length);
    for (i = 0; i < cactusGraph->vertices->length; i++) {
        vertexOrdering[i] = -1;
    }

    getDFSDiscoveryTimesP(cactusGraph, cactusGraph->vertices->list[0], 0, vertexOrdering);

    for (i = 0; i < cactusGraph->vertices->length; i++) {
        assert(vertexOrdering[i] >= 0);
        assert(vertexOrdering[i] < cactusGraph->vertices->length);
    }

    return vertexOrdering;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods to interact with the 3-edge connected component
//code.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct List *readThreeEdgeComponents(struct PinchGraph *pinchGraph, stList *adjacencyComponents,
        stList *threeEdgeConnectedAdjacencyComponents) {
    /*
     * Reads in the three edge connected components written out by the three
     * edge script.
     */
    int32_t i, j;
    struct PinchVertex *vertex;

#ifdef BEN_DEBUG
    int32_t l = 0;
#endif
    struct List *threeEdgeConnectedComponents = constructEmptyList(0, (void(*)(void *)) destructList);
    for (i = 0; i < stList_length(threeEdgeConnectedAdjacencyComponents); i++) {
        stList *threeEdgeConnectedAdjacencyComponent = stList_get(threeEdgeConnectedAdjacencyComponents, i);
        struct List *threeEdgeConnectedComponent = constructEmptyList(0, NULL);
        listAppend(threeEdgeConnectedComponents, threeEdgeConnectedComponent);
        for (j = 0; j < stList_length(threeEdgeConnectedAdjacencyComponent); j++) {
            stSortedSet *adjacencyComponent = stList_get(adjacencyComponents,
                    stIntTuple_getPosition(stList_get(threeEdgeConnectedAdjacencyComponent, j), 0));
            assert(adjacencyComponent != NULL);
            stSortedSetIterator *it = stSortedSet_getIterator(adjacencyComponent);
            while ((vertex = stSortedSet_getNext(it)) != NULL) {
                listAppend(threeEdgeConnectedComponent, vertex);
#ifdef BEN_DEBUG
                l++;
#endif
            }
            stSortedSet_destructIterator(it);
        }
    }
#ifdef BEN_DEBUG
    assert(l == pinchGraph->vertices->length);
    int32_t k = FALSE;
#endif
    for (i = 0; i < threeEdgeConnectedComponents->length; i++) {
        struct List *threeEdgeConnectedComponent = threeEdgeConnectedComponents->list[i];
        for (j = 0; j < threeEdgeConnectedComponent->length; j++) {
            vertex = threeEdgeConnectedComponent->list[j];
            if (vertex->vertexID == 0) {
#ifdef BEN_DEBUG
                assert(k == FALSE);
                k = TRUE;
#endif
                struct List *list3 = threeEdgeConnectedComponents->list[0];
                threeEdgeConnectedComponents->list[0] = threeEdgeConnectedComponent;
                threeEdgeConnectedComponents->list[i] = list3;
            }
        }
    }
#ifdef BEN_DEBUG
    assert(k == TRUE);
#endif

    return threeEdgeConnectedComponents;
}

void writeOutCactusGraph(struct CactusGraph *cactusGraph, struct PinchGraph *pinchGraph, FILE *fileHandle) {
    assert(pinchGraph != NULL);
    /*
     * Writes out a graph in 'dot' format, compatible with graphviz.
     *
     * The associated function 'writeOutEdgePieces' gives a way of decoding
     * the pieces associated with the black (piece containing) edges
     * of the graph.
     */
    int32_t i, j;
    struct CactusVertex *vertex;
    struct CactusEdge *edge;

    //Write the preliminaries.
    fprintf(fileHandle, "graph G {\n");
    fprintf(fileHandle, "overlap=false\n");
    //Write the vertices.
    for (i = 0; i < cactusGraph->vertices->length; i++) {
        vertex = cactusGraph->vertices->list[i];
#ifdef BEN_DEBUG
        assert(vertex->vertexID == i);
#endif
        fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n", vertex->vertexID, vertex->vertexID);
    }

    fprintf(fileHandle, "edge[color=black,len=2.5,weight=100,dir=forward];\n");
    for (i = 0; i < cactusGraph->vertices->length; i++) {
        vertex = cactusGraph->vertices->list[i];
        for (j = 0; j < vertex->edges->length; j++) {
            edge = vertex->edges->list[j];
#ifdef BEN_DEBUG
            assert(edge != edge->rEdge);
#endif
            if (edge > edge->rEdge) {
                fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n", edge->from->vertexID, edge->to->vertexID);
            }
        }
    }
    fprintf(fileHandle, "}\n");
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Script to construct the cactus graph
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


struct CactusGraph *computeCactusGraph(struct PinchGraph *pinchGraph, bool(*passThroughEdgeFn)(struct PinchEdge *)) {
    ///////////////////////////////////////////////////////////////////////////
    // Run the three-edge connected component algorithm to identify
    // three edge connected components.
    ///////////////////////////////////////////////////////////////////////////

    stList *adjacencyComponents = getAdjacencyComponents2(pinchGraph, passThroughEdgeFn);
    stHash *vertexToAdjacencyComponentHash = getVertexToAdjacencyComponentHash(pinchGraph, adjacencyComponents);
    stList *adjacencyComponentGraph = getAdjacencyComponentGraph(pinchGraph, adjacencyComponents,
            vertexToAdjacencyComponentHash);

    stList *threeEdgeConnectedAdjacencyComponents = computeThreeEdgeConnectedComponents(adjacencyComponentGraph);
    st_logInfo("Seems to have successfully run the three edge command: %i\n",
            stList_length(threeEdgeConnectedAdjacencyComponents));

    //Parse results (the three edge connected components).
    struct List *threeEdgeConnectedComponents = readThreeEdgeComponents(pinchGraph, adjacencyComponents,
            threeEdgeConnectedAdjacencyComponents);
    st_logInfo("Read in the three edge components\n");

    //Cleanup
    stList_destruct(threeEdgeConnectedAdjacencyComponents);
    stList_destruct(adjacencyComponents);
    stHash_destruct(vertexToAdjacencyComponentHash);
    stList_destruct(adjacencyComponentGraph);

    ///////////////////////////////////////////////////////////////////////////
    // Merge vertices in each three edge connected component to convert the main graph
    // into a 'cactus' graph.
    ///////////////////////////////////////////////////////////////////////////

    //Collapse the graph to cactus tree
    struct CactusGraph *cactusGraph = constructCactusGraph(pinchGraph, threeEdgeConnectedComponents, passThroughEdgeFn);
    checkCactusGraph(pinchGraph, threeEdgeConnectedComponents, cactusGraph);

    //Cleanup
    destructList(threeEdgeConnectedComponents);
    return cactusGraph;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods to manipulate the cactus graph
//
//First methods are for circularising stems.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void cactusVertex_merge(struct CactusGraph *cactusGraph, struct CactusVertex *one, struct CactusVertex *two) {
    /*
     * Merges two into one, destroying two and adding to one.
     */
    assert(one != two);
    assert(one->vertexID != two->vertexID);
    //redirect all edges to one.
    while (two->edges->length > 0) { //handles self loops on second because it contains both orientations
        struct CactusEdge *edge = two->edges->list[--two->edges->length];
#if BEN_DEBUG
        assert(edge->from == two);
        assert(edge->rEdge->to == two);
        assert(!listContains(one->edges, edge));
#endif
        edge->from = one;
        edge->rEdge->to = one;
        listAppend(one->edges, edge);
    }
    //switch the vertex id's of the highest cactus vertex and two.
    struct CactusVertex *three = cactusGraph->vertices->list[cactusGraph->vertices->length - 1];
    assert(three->vertexID == cactusGraph->vertices->length - 1);
    three->vertexID = two->vertexID;
    assert(cactusGraph->vertices->list[two->vertexID] == two);
    cactusGraph->vertices->list[two->vertexID] = three;
    cactusGraph->vertices->length--;

    //now destruct two, but not its list of edges.
    destructCactusVertex(two);
}

static void circulariseStemsP(struct CactusGraph *cactusGraph, struct CactusEdge *edge,
        struct CactusVertex *sourceVertex, struct hashtable *stemHash, struct hashtable *seen,
        struct List *verticesToMerge, struct PinchGraph *pinchGraph, Flower *flower) {
    int32_t i, j;
    struct CactusEdge *edge2;

    if (hashtable_search(seen, edge->to) == NULL) {
        hashtable_insert(seen, edge->to, edge->to);
        assert(hashtable_search(seen, edge->from) != NULL);
        assert(edge->from != edge->to);
        if (isAFreeStubCactusEdge(edge, pinchGraph, flower)) { //is a free stub end, so we circularise it into its own component.
            assert(hashtable_search(stemHash, edge) == NULL);
            assert(hashtable_search(stemHash, edge->rEdge) == NULL);
            listAppend(verticesToMerge, edge->from);
            listAppend(verticesToMerge, edge->to);
        } else if (hashtable_search(stemHash, edge) != NULL) { //is a stem
            j = 0;
            for (i = 0; i < edge->to->edges->length; i++) {
                edge2 = edge->to->edges->list[i];
                assert(edge2->from == edge->to);
                assert(edge2 != edge);
                if (hashtable_search(stemHash, edge2) != NULL) {//is a stem
                    assert(edge2->from != edge2->to);
                    j++;
                }
            }
            assert(j > 0);
            if (j != 2) { //if is either branch or leaf.
                assert(listContains(edge->to->edges, edge->rEdge));
                listAppend(verticesToMerge, sourceVertex);
                listAppend(verticesToMerge, edge->to);
            }
        } else {
            sourceVertex = edge->to;
        }
        for (i = 0; i < edge->to->edges->length; i++) { //call recursively.
            edge2 = edge->to->edges->list[i];
            circulariseStemsP(cactusGraph, edge2, sourceVertex, stemHash, seen, verticesToMerge, pinchGraph, flower);
        }
    }
}

void circulariseStems(struct CactusGraph *cactusGraph, struct PinchGraph *pinchGraph, Flower *flower) {
    struct List *biConnectedComponent;
    struct List *biConnectedComponents;
    struct hashtable *stemHash;
    struct hashtable *seen;
    struct CactusEdge *edge;
    struct CactusVertex *vertex;
    struct List *cactusGraphMerges;
    int32_t i;

    st_logDebug("Circularising the stems\n");

    ////////////////////////////////////////////////
    //(1) Get bi-connected components.
    ////////////////////////////////////////////////
    biConnectedComponents = computeBiConnectedComponents(cactusGraph);
    st_logDebug("Constructed the biconnected components for making the flower\n");

    ////////////////////////////////////////////////
    //(2) Put stems in a hash
    ////////////////////////////////////////////////

    stemHash = create_hashtable(cactusGraph->vertices->length * 2, hashtable_key, hashtable_equalKey, NULL, NULL);
    for (i = 0; i < biConnectedComponents->length; i++) {
        biConnectedComponent = biConnectedComponents->list[i];
        if (biConnectedComponent->length == 1) {
            edge = biConnectedComponent->list[0];
            if (edge->from != edge->to) { //must be a stem as not a self loop
                if (!isAFreeStubCactusEdge(edge, pinchGraph, flower)) { //Is not a free stub end
                    hashtable_insert(stemHash, edge, edge);
                    hashtable_insert(stemHash, edge->rEdge, edge->rEdge);
                }
            }
        }
    }
    st_logDebug("Put the stems in a hash\n");

    ////////////////////////////////////////////////
    //(3) Do DFS on graph to find vertex merges that we wish to make
    ////////////////////////////////////////////////

    seen = create_hashtable(cactusGraph->vertices->length * 2, hashtable_key, hashtable_equalKey, NULL, NULL);
    vertex = cactusGraph->vertices->list[0];
    hashtable_insert(seen, vertex, vertex);
    cactusGraphMerges = constructEmptyList(0, NULL);
    for (i = 0; i < vertex->edges->length; i++) {
        circulariseStemsP(cactusGraph, vertex->edges->list[i], vertex, stemHash, seen, cactusGraphMerges, pinchGraph,
                flower);
    }
    st_logDebug("Done the DFS\n");

    ////////////////////////////////////////////////
    //(4) Do vertex merges.
    ////////////////////////////////////////////////

    assert(hashtable_count(seen) == cactusGraph->vertices->length);
    assert((cactusGraphMerges->length % 2) == 0);
    while (cactusGraphMerges->length > 0) {
        struct CactusVertex *stemVertex = cactusGraphMerges->list[--cactusGraphMerges->length];
        assert(cactusGraphMerges->length > 0);
        struct CactusVertex *sourceVertex = cactusGraphMerges->list[--cactusGraphMerges->length];
#ifdef BEN_DEBUG
        assert(cactusGraph->vertices->list[sourceVertex->vertexID] == sourceVertex);
        assert(cactusGraph->vertices->list[stemVertex->vertexID] == stemVertex);
        assert(!listContains(cactusGraphMerges, stemVertex));
#endif
        cactusVertex_merge(cactusGraph, sourceVertex, stemVertex);
    }
    assert(cactusGraphMerges->length == 0);

    ////////////////////////////////////////////////
    //(4) Do loop back merges of free stub ends
    ////////////////////////////////////////////////

    for (i = 0; i < cactusGraph->vertices->length; i++) {
        vertex = cactusGraph->vertices->list[i];
        for (int32_t j = 0; j < vertex->edges->length; j++) {
            edge = vertex->edges->list[j];
            if (edge->from != edge->to && isAFreeStubCactusEdge(edge, pinchGraph, flower)) {
                cactusVertex_merge(cactusGraph, edge->from, edge->to);
            }
        }
    }

    ////////////////////////////////////////////////
    //(5) Cleanup
    ////////////////////////////////////////////////

    hashtable_destroy(stemHash, FALSE, FALSE);
    hashtable_destroy(seen, FALSE, FALSE);
    destructList(biConnectedComponents);
    destructList(cactusGraphMerges);

#ifdef BEN_DEBUG
    biConnectedComponents = computeBiConnectedComponents(cactusGraph);
    for (i = 0; i < biConnectedComponents->length; i++) {
        biConnectedComponent = biConnectedComponents->list[i];
        assert(biConnectedComponent->length > 0);
        if (biConnectedComponent->length == 1) {
            edge = biConnectedComponent->list[0];
            assert(edge->from == edge->to);
        }
        for (int32_t j = 0; j < biConnectedComponent->length; j++) {
            edge = biConnectedComponent->list[j];
            if (isAFreeStubCactusEdge(edge, pinchGraph, flower)) {
                assert(biConnectedComponent->length == 1);
            }
        }
    }
    destructList(biConnectedComponents);
#endif
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Choosing which blocks in the chains to keep.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


float treeCoverage2(struct CactusEdge *cactusEdge, Flower *flower, struct PinchGraph *pinchGraph) {
    /*
     * Returns the proportion of the tree covered by the block.
     */
#ifdef BEN_DEBUG
    assert(!isAStubCactusEdge(cactusEdge, pinchGraph));
#endif
    return treeCoverage(cactusEdgeToFirstPinchEdge(cactusEdge, pinchGraph)->from, flower);
}

stSortedSet *getEventStrings(struct PinchVertex *vertex, Flower *flower) {
    void *blackEdgeIterator = getBlackEdgeIterator(vertex);
    struct PinchEdge *edge;
    stSortedSet *eventStrings = stSortedSet_construct3((int(*)(const void *, const void *)) strcmp, free);
    while ((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
        struct Piece *piece = edge->piece;
        Sequence *sequence = flower_getSequence(flower, piece->contig);
        assert(sequence != NULL);
        Event *event = sequence_getEvent(sequence);
        assert(event != NULL);
        if (stSortedSet_search(eventStrings, (void *) event_getHeader(event)) == NULL) {
            stSortedSet_insert(eventStrings, stString_copy(event_getHeader(event)));
        }
    }
    destructBlackEdgeIterator(blackEdgeIterator);
    return eventStrings;
}

static bool containsRequiredSpecies(struct CactusEdge *cactusEdge, Flower *flower, struct PinchGraph *pinchGraph,
        stList *listOfSetsOfRequiredSpecies, /* A block's segments must have an event with at least coverage (given for each set individually) number of elements in each of these set (unless it is NULL) */
        stList *listOfRequiredSpeciesCoverages) {
    /*
     * Returns the proportion of the tree covered by the block.
     */
#ifdef BEN_DEBUG
    assert(!isAStubCactusEdge(cactusEdge, pinchGraph));
#endif
    stSortedSet *eventStrings = getEventStrings(cactusEdgeToFirstPinchEdge(cactusEdge, pinchGraph)->from, flower);
    for(int32_t i=0; i<stList_length(listOfRequiredSpeciesCoverages); i++) {
        int32_t coverage = stIntTuple_getPosition(stList_get(listOfRequiredSpeciesCoverages, i), 0);
        assert(coverage > 0);
        stSortedSet *requiredSpecies = stList_get(listOfSetsOfRequiredSpecies, i);
        stSortedSet *commonEvents = stSortedSet_getIntersection(requiredSpecies, eventStrings);
        if(stSortedSet_size(commonEvents) < coverage) {
            stSortedSet_destruct(commonEvents);
            stSortedSet_destruct(eventStrings);
            return 0;
        }
        stSortedSet_destruct(commonEvents);
    }
    stSortedSet_destruct(eventStrings);
    return 1;
}

static bool containsMultipleCopiesOfSpecies(struct CactusEdge *cactusEdge, Flower *flower,
        struct PinchGraph *pinchGraph, stSortedSet *singleCopySpecies) {
    /*
     * Returns 1 iff the edge contains more than copy of any of the given sequences.
     */
    struct PinchVertex *vertex = cactusEdgeToFirstPinchEdge(cactusEdge, pinchGraph)->from;
    void *blackEdgeIterator = getBlackEdgeIterator(vertex);
    struct PinchEdge *edge;
    stSortedSet *eventStrings = stSortedSet_construct3((int(*)(const void *, const void *)) strcmp, free);
    bool b = 0;
    while ((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
        struct Piece *piece = edge->piece;
        Sequence *sequence = flower_getSequence(flower, piece->contig);
        assert(sequence != NULL);
        Event *event = sequence_getEvent(sequence);
        assert(event != NULL);
        if (stSortedSet_search(singleCopySpecies, (void *) event_getHeader(event)) != NULL) {
            if (stSortedSet_search(eventStrings, (void *) event_getHeader(event)) == NULL) {
                stSortedSet_insert(eventStrings, stString_copy(event_getHeader(event)));
            } else {
               b = 1;
               break;
            }
        }
    }
    destructBlackEdgeIterator(blackEdgeIterator);
    stSortedSet_destruct(eventStrings);
    return b;
}

int32_t chainLength(struct List *biConnectedComponent, int32_t includeStubs, struct PinchGraph *pinchGraph) {
    /*
     * Get the number of links in the chain.
     */
    int32_t i, j;
    struct CactusEdge *cactusEdge;
    i = 0;
    for (j = 0; j < biConnectedComponent->length; j++) {
        cactusEdge = biConnectedComponent->list[j];
        if (includeStubs || !isAStubCactusEdge(cactusEdge, pinchGraph)) {
            i++;
        }
    }
    return i;
}

int32_t maxChainDegreeOfNonStubBlocks(struct List *biConnectedComponent, struct PinchGraph *pinchGraph) {
    int32_t i = 0;
    for (int32_t j = 0; j < biConnectedComponent->length; j++) {
        struct CactusEdge *cactusEdge = biConnectedComponent->list[j];
        if (!isAStubCactusEdge(cactusEdge, pinchGraph) && cactusEdge->pieces->length > i) {
            i = cactusEdge->pieces->length;
        }
    }
    return i;
}

int32_t chainBaseLength(struct List *biConnectedComponent, struct PinchGraph *pinchGraph) {
    int64_t i, j;
    struct CactusEdge *cactusEdge;
    struct Piece *piece;

    i = 0;
    for (j = 0; j < biConnectedComponent->length; j++) {
        cactusEdge = biConnectedComponent->list[j];
        if (!isAStubCactusEdge(cactusEdge, pinchGraph)) {
            piece = cactusEdge->pieces->list[0];
            i += ((int64_t)piece->end) - ((int64_t)piece->start + 1);
        }
    }
    assert(i <= INT32_MAX);
    return (int32_t)i;
}

stSortedSet *getPinchVerticesSet(stSortedSet *cactusEdges, struct PinchGraph *pinchGraph) {
    stSortedSet *pinchVerticesSet = stSortedSet_construct();
    stSortedSetIterator *it = stSortedSet_getIterator(cactusEdges);
    struct CactusEdge *edge;
    while ((edge = stSortedSet_getNext(it)) != NULL) {
        struct PinchEdge *pinchEdge = cactusEdgeToFirstPinchEdge(edge, pinchGraph);
#ifdef BEN_DEBUG
        assert(stSortedSet_search(pinchVerticesSet, pinchEdge->from) == NULL);
        assert(stSortedSet_search(pinchVerticesSet, pinchEdge->to) == NULL);
#endif
        stSortedSet_insert(pinchVerticesSet, pinchEdge->to);
        stSortedSet_insert(pinchVerticesSet, pinchEdge->from);
    }
    stSortedSet_destructIterator(it);
    return pinchVerticesSet;
}

stSortedSet *filterBlocksByTreeCoverageAndLength(struct List *biConnectedComponents, Flower *flower,
        float minimumTreeCoverage, /*Minimum tree coverage to be included (>=) */
        int32_t minimumBlockDegree, /*The minimum number of segments in a block to be included (>=)*/
        int32_t minimumBlockLength, /*The minimum length of an block to be included (>=)*/
        int32_t minimumChainLength, /* Minimum chain length to be included (>=)*/
        stList *listOfSetsOfRequiredSpecies, /* A block's segments must have an event with at least coverage (given for each set individually) number of elements in each of these set (unless it is NULL) */
        stList *listOfRequiredSpeciesCoverages,
        stSortedSet *singleCopySpecies, /* A block must not have multiple (more than one) segments from any of the following species (unless it is NULL)*/
        struct PinchGraph *pinchGraph) {
    /*
     * Filters blocks in chains by base length and tree coverage.
     *
     * Returns a list of all accepted blocks, excluding stubs.
     */
    stSortedSet *chosenBlocks = stSortedSet_construct();
    for (int32_t i = 0; i < biConnectedComponents->length; i++) {
        struct List *biConnectedComponent = biConnectedComponents->list[i];
        if (minimumChainLength <= 0 || chainBaseLength(biConnectedComponent, pinchGraph) >= minimumChainLength) {
            for (int32_t j = 0; j < biConnectedComponent->length; j++) {
                struct CactusEdge *cactusEdge = biConnectedComponent->list[j];
                if (!isAStubCactusEdge(cactusEdge, pinchGraph)) {
                    assert(cactusEdge->pieces->length > 0);
                    if (minimumBlockDegree <= 0 || cactusEdge->pieces->length >= minimumBlockDegree) {
                        if (minimumTreeCoverage <= 0.0 || treeCoverage2(cactusEdge, flower, pinchGraph)
                                >= minimumTreeCoverage) {
                            struct Piece *piece = cactusEdge->pieces->list[0];
                            if (minimumBlockLength <= 0 || piece->end - piece->start + 1 >= minimumBlockLength) {
                                if (listOfSetsOfRequiredSpecies == NULL || containsRequiredSpecies(cactusEdge, flower, pinchGraph,
                                        listOfSetsOfRequiredSpecies, listOfRequiredSpeciesCoverages)) {
                                    if (singleCopySpecies == NULL || !containsMultipleCopiesOfSpecies(cactusEdge,
                                            flower, pinchGraph, singleCopySpecies)) {
                                        stSortedSet_insert(chosenBlocks, cactusEdge);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return chosenBlocks;
}

void logTheChosenBlockSubset(struct List *biConnectedComponents, stSortedSet *chosenBlocks,
        struct PinchGraph *pinchGraph, Flower *flower) {
    /*
     * Produces logging information about the chosen blocks.
     */
    int32_t i, j;
    struct List *biConnectedComponent;
    struct CactusEdge *cactusEdge;
    struct Piece *piece;
    float totalBlockScore = 0.0;
    float totalBlockLength = 0.0;
    float totalBaseLengthOfAllBlocks = 0.0;
    float totalNumberOfAllBlocks = 0.0;
    float totalNumberOfStubBlocks = 0.0;
    float averagePieceNumber = 0.0;
    float averagePieceNumberOfAllBlocks = 0.0;
    for (i = 0; i < biConnectedComponents->length; i++) {
        biConnectedComponent = biConnectedComponents->list[i];
        totalBaseLengthOfAllBlocks += chainBaseLength(biConnectedComponent, pinchGraph);
        totalNumberOfAllBlocks += chainLength(biConnectedComponent, FALSE, pinchGraph);
        totalNumberOfStubBlocks += chainLength(biConnectedComponent, TRUE, pinchGraph) - chainLength(
                biConnectedComponent, FALSE, pinchGraph);
        for (j = 0; j < biConnectedComponent->length; j++) {
            cactusEdge = biConnectedComponent->list[j];
            if (!isAStubCactusEdge(cactusEdge, pinchGraph)) {
                averagePieceNumberOfAllBlocks += cactusEdge->pieces->length;
            }
        }
    }
    j = 0;
    stSortedSetIterator *it = stSortedSet_getIterator(chosenBlocks);
    while ((cactusEdge = stSortedSet_getNext(it)) != NULL) {
        if (!isAStubCactusEdge(cactusEdge, pinchGraph)) {
            totalBlockScore += treeCoverage2(cactusEdge, flower, pinchGraph);
            piece = cactusEdge->pieces->list[0];
            totalBlockLength += piece->end - piece->start + 1;
            averagePieceNumber += cactusEdge->pieces->length;
            j++;
        }
    }
    stSortedSet_destructIterator(it);
    st_logInfo(
            "Chosen block subset composed of %i blocks, of average length %f and average tree coverage %f, total base length of all blocks: %f, total number of all blocks %f, average length of all blocks: %f, total number of stub blocks: %f, average piece number of chosen blocks: %f, average piece number of all blocks: %f\n",
            stSortedSet_size(chosenBlocks), totalBlockLength / j, totalBlockScore / j, totalBaseLengthOfAllBlocks,
            totalNumberOfAllBlocks, totalBaseLengthOfAllBlocks / totalNumberOfAllBlocks, totalNumberOfStubBlocks,
            averagePieceNumber / j, averagePieceNumberOfAllBlocks / totalNumberOfAllBlocks);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Cactus graph misc functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchEdge *cactusEdgeToFirstPinchEdge(struct CactusEdge *edge, struct PinchGraph *pinchGraph) {
    struct Piece *piece;
#ifdef BEN_DEBUG
    assert(edge->pieces->length> 0);
#endif
    piece = edge->pieces->list[0];
    return getContainingBlackEdge(pinchGraph, piece->contig, piece->start);
}

int32_t isAStubCactusEdge(struct CactusEdge *edge, struct PinchGraph *pinchGraph) {
    return isAStub(cactusEdgeToFirstPinchEdge(edge, pinchGraph));
}

int32_t isAFreeStubCactusEdge(struct CactusEdge *edge, struct PinchGraph *pinchGraph, Flower *flower) {
    if (isAStubCactusEdge(edge, pinchGraph)) {
        struct PinchEdge *pinchEdge;
        assert(edge->pieces->length >= 1);
        pinchEdge = cactusEdgeToFirstPinchEdge(edge, pinchGraph);
        assert(isAStub(pinchEdge));
        assert(vertex_isDeadEnd(pinchEdge->from) || vertex_isDeadEnd(pinchEdge->to));
        assert(!(vertex_isDeadEnd(pinchEdge->from) && vertex_isDeadEnd(pinchEdge->to)));
        Cap *cap = flower_getCap(flower, pinchEdge->piece->contig);
        assert(cap != NULL);
        End *end = cap_getEnd(cap);
        return end_isStubEnd(end) && end_isFree(end);
    }
    return 0;
}

struct hashtable *createHashColouringPinchEdgesByChains(struct PinchGraph *pinchGraph,
        struct List *biConnectComponentsList) {
    struct List *biConnectedComponent;
    struct CactusEdge *cactusEdge;
    struct Piece *piece;
    int32_t i, j, k;

    //Put the chain pieces in a hash to colour the black edges of the pinch graph.
    struct hashtable *hash = create_hashtable(pinchGraph->vertices->length * 10, hashtable_key, hashtable_equalKey,
            NULL, (void(*)(void *)) destructInt);

    if (biConnectComponentsList != NULL) {
        for (i = 0; i < biConnectComponentsList->length; i++) {
            biConnectedComponent = biConnectComponentsList->list[i];
            for (k = 0; k < biConnectedComponent->length; k++) {
                cactusEdge = biConnectedComponent->list[k];
                for (j = 0; j < cactusEdge->pieces->length; j++) {
                    piece = cactusEdge->pieces->list[j];
                    hashtable_insert(hash, piece, constructInt(i));
                    hashtable_insert(hash, piece->rPiece, constructInt(i));
                }
            }
        }
    }

    return hash;
}

