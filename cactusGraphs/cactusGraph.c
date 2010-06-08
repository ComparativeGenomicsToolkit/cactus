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
    cactusVertex->edges = constructEmptyList(0,
            (void(*)(void *)) destructCactusEdge);

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

struct CactusEdge *constructCactusEdge2(struct List *pieces,
        struct CactusVertex *from, struct CactusVertex *to) {
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

struct CactusGraph *constructCactusGraph(struct PinchGraph *pinchGraph,
        struct List *threeEdgeConnectedComponents) {
    struct CactusGraph *cactusGraph;
    struct CactusEdge *cactusEdge;
    struct CactusVertex *cactusVertex;
    struct CactusVertex *cactusVertex2;
    int32_t i, j;
    void *greyEdgeIterator;
    struct List *list;
    struct PinchVertex *pinchVertex;
    struct PinchVertex *pinchVertex2;
    struct PinchEdge *pinchEdge;
    struct List *pinchVertexToCactusVertex;
    struct List *emptyList;
    struct List *list2;

    cactusGraph = st_malloc(sizeof(struct CactusGraph));
    cactusGraph->vertices = constructEmptyList(0,
            (void(*)(void *)) destructCactusVertex);

#ifdef BEN_DEBUG
    list = threeEdgeConnectedComponents->list[0];
    j = FALSE;
    for (i = 0; i < list->length; i++) {
        pinchVertex = list->list[0];
        if (pinchVertex->vertexID == 0) {
            j = TRUE;
        }
    }
    assert(j == TRUE);
#endif

    pinchVertexToCactusVertex = constructEmptyList(
            pinchGraph->vertices->length, NULL);
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
            pinchVertexToCactusVertex->list[pinchVertex->vertexID]
                    = cactusVertex;
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
                pinchVertex2 = pinchEdge->to;
                cactusVertex2
                        = pinchVertexToCactusVertex->list[pinchVertex2->vertexID];

                if (pinchVertex->vertexID < pinchVertex2->vertexID) {
                    //getting the ordered pieces.
                    list2 = constructEmptyList(0, NULL);
                    void *blackEdgeIterator = getBlackEdgeIterator(pinchVertex);
                    while ((pinchEdge = getNextBlackEdge(pinchVertex,
                            blackEdgeIterator)) != NULL) {
                        listAppend(list2, pinchEdge->piece);
                    }
                    destructBlackEdgeIterator(blackEdgeIterator);

                    cactusEdge = constructCactusEdge2(list2, cactusVertex,
                            cactusVertex2);
                    destructList(list2);
                }
            }
#ifdef BEN_DEBUG
            else {
                assert(pinchVertex->vertexID == 0);
            }
#endif

            //grey edges
            greyEdgeIterator = getGreyEdgeIterator(pinchVertex);
            while ((pinchVertex2 = getNextGreyEdge(pinchVertex,
                    greyEdgeIterator)) != NULL) {
                cactusVertex2
                        = pinchVertexToCactusVertex->list[pinchVertex2->vertexID];

                if (cactusVertex != cactusVertex2 && cactusVertex
                        < cactusVertex2) {
                    cactusEdge = constructCactusEdge2(emptyList, cactusVertex,
                            cactusVertex2);
                }
            }
            destructGreyEdgeIterator(greyEdgeIterator);
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
    for(int32_t i=0; i<cactusGraph->vertices->length; i++) {
        struct CactusVertex *cactusVertex = cactusGraph->vertices->list[i];
        for(int32_t j=0; j<cactusVertex->edges->length; j++) {
            struct CactusEdge *cactusEdge = cactusVertex->edges->list[j];
            stSortedSet_insert(sortedSet, cactusEdge);
            stSortedSet_insert(sortedSet, cactusEdge->rEdge);
        }
    }
    int32_t k = stSortedSet_size(sortedSet);
    stSortedSet_destruct(sortedSet);
    assert(k % 2 == 0);
    return k/2;
}

void checkCactusGraph(struct PinchGraph *pinchGraph,
        struct List *threeEdgeConnectedComponents,
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

void checkCactusContainsOnly2EdgeConnectedComponents(
        struct CactusGraph *cactusGraph) {
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
    st_logDebug("Constructed the biconnected components for making the net\n");

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

void computeBiConnectedComponents_P(struct hashtable *flag, int32_t count,
        int32_t *dFN, int32_t *low, struct List *stack, int32_t *father,
        struct List *biConnnectedComponents, struct CactusVertex *v) {
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
            computeBiConnectedComponents_P(flag, count, dFN, low, stack,
                    father, biConnnectedComponents, w);
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
            low[v->vertexID]
                    = low[v->vertexID] < low[w->vertexID] ? low[v->vertexID]
                            : low[w->vertexID];
        } else {
            if (w->vertexID != father[v->vertexID]) {
                low[v->vertexID]
                        = low[v->vertexID] < dFN[w->vertexID] ? low[v->vertexID]
                                : dFN[w->vertexID];
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

    flag = create_hashtable(cactusGraph->vertices->length * 10, hashtable_key,
            hashtable_equalKey, NULL, NULL);
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
    biConnectedComponents = constructEmptyList(0,
            (void(*)(void *)) destructList);

    computeBiConnectedComponents_P(flag, count, dFN, low, stack, father,
            biConnectedComponents, cactusGraph->vertices->list[0]);

    hashtable_destroy(flag, FALSE, FALSE);
    free(dFN);
    free(low);
    free(father);
    destructList(stack);

    return biConnectedComponents;
}

int32_t *sortBiConnectedComponents_vertexOrdering;

int sortBiConnectedComponentsP(struct CactusEdge **edge,
        struct CactusEdge **edge2) {
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
    st_logDebug("Constructed the biconnected components for making the net\n");

    ////////////////////////////////////////////////
    //(1) Get DFS ordering on nodes
    ////////////////////////////////////////////////
    sortBiConnectedComponents_vertexOrdering
            = getDFSDiscoveryTimes(cactusGraph);
    st_logDebug("Got the DFS discovery times\n");

    ////////////////////////////////////////////////
    //(2) Now sort each bi connected component
    ////////////////////////////////////////////////
    for (i = 0; i < biConnectedComponents->length; i++) {
        biConnectedComponent = biConnectedComponents->list[i];
        qsort(
                biConnectedComponent->list,
                biConnectedComponent->length,
                sizeof(void *),
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
    return biConnectedComponents;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method to compute DFS discovery time for vertices.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t getDFSDiscoveryTimesP(struct CactusGraph *cactusGraph,
        struct CactusVertex *vertex, int32_t counter, int32_t *vertexOrdering) {
    int32_t i;
    struct CactusEdge *edge;

    assert(vertexOrdering[vertex->vertexID] == -1);
    vertexOrdering[vertex->vertexID] = counter++;
    for (i = 0; i < vertex->edges->length; i++) {
        edge = vertex->edges->list[i];
        if (vertexOrdering[edge->to->vertexID] == -1) {
            counter = getDFSDiscoveryTimesP(cactusGraph, edge->to, counter,
                    vertexOrdering);
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

    getDFSDiscoveryTimesP(cactusGraph, cactusGraph->vertices->list[0], 0,
            vertexOrdering);

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

stList *writeOut3EdgeGraph(struct PinchGraph *pinchGraph,
        struct List *greyEdgeComponents) {
    /*
     * Writes a format compatible with the 3-edge connected component algorithm.
     */
    struct PinchVertex *vertex;
    struct PinchEdge *edge;
    int32_t i, j, k;
    struct List *component;
    struct hashtable *vertexHash;
    stList *vertices = stList_construct3(0,
            (void(*)(void *)) destructIntList);

    //setup vertex to grey edge component hash
    vertexHash = create_hashtable(pinchGraph->vertices->length * 2,
            hashtable_key, hashtable_equalKey, NULL,
            (void(*)(void *)) destructInt);

    for (i = 0; i < greyEdgeComponents->length; i++) {
        component = greyEdgeComponents->list[i];
        for (j = 0; j < component->length; j++) {
            vertex = component->list[j];
            hashtable_insert(vertexHash, vertex, constructInt(i));
        }
    }
#ifdef BEN_ULTRA_DEBUG
    assert((int32_t)hashtable_count(vertexHash) == pinchGraph->vertices->length);

    for(i=0; i<pinchGraph->vertices->length; i++) {
        vertex = pinchGraph->vertices->list[i];
        j = *((int32_t *)hashtable_search(vertexHash, vertex));
        void *greyEdgeIterator = getGreyEdgeIterator(vertex);
        while((vertex2 = getNextGreyEdge(vertex, greyEdgeIterator)) != NULL) {
            l = *((int32_t *)hashtable_search(vertexHash, vertex2));
            assert(j == l);
        }
        destructGreyEdgeIterator(greyEdgeIterator);
    }
#endif

    //Write number of nodes.
    for (i = 0; i < greyEdgeComponents->length; i++) {
        component = greyEdgeComponents->list[i];
        struct IntList *edges = constructEmptyIntList(0);
        stList_append(vertices, edges);
        for (j = 0; j < component->length; j++) {
            vertex = component->list[j];

            //The black edges
            if (lengthBlackEdges(vertex) > 0) {
                edge = getFirstBlackEdge(vertex);
                k = *((int32_t *) hashtable_search(vertexHash, edge->to));
                intListAppend(edges, k + 1);
            }
#ifdef BEN_DEBUG
            else {
                assert(vertex->vertexID == 0);
            }
#endif
        }
    }
    hashtable_destroy(vertexHash, TRUE, FALSE);

    return vertices;
}

struct List *readThreeEdgeComponents(struct PinchGraph *pinchGraph,
        struct List *greyEdgeComponents, stList *threeEdgeComponents) {
    /*
     * Reads in the three edge connected components written out by the three
     * edge script.
     */
    int32_t i, j, k;
    struct List *list;
    struct List *list2;
    stList *list3;
    struct List *component;
    struct PinchVertex *vertex;

#ifdef BEN_DEBUG
    int32_t l = 0;
#endif
    list = constructEmptyList(0, (void(*)(void *)) destructList);
    for (i = 0; i < stList_length(threeEdgeComponents); i++) {
        list3 = stList_get(threeEdgeComponents, i);
        list2 = constructEmptyList(0, NULL);
        listAppend(list, list2);
        for (j = 0; j < stList_length(list3); j++) {
            component = greyEdgeComponents->list[stIntTuple_getPosition(stList_get(list3, j), 0) - 1];
            assert(component != NULL);
            for (k = 0; k < component->length; k++) {
                vertex = component->list[k];
                listAppend(list2, vertex);
#ifdef BEN_DEBUG
                l++;
#endif
            }
        }
    }
#ifdef BEN_DEBUG
    assert(l == pinchGraph->vertices->length);
    k = FALSE;
#endif
    for (i = 0; i < list->length; i++) {
        list2 = list->list[i];
        for (j = 0; j < list2->length; j++) {
            vertex = list2->list[j];
            if (vertex->vertexID == 0) {
#ifdef BEN_DEBUG
                assert(k == FALSE);
                k = TRUE;
#endif
                list3 = list->list[0];
                list->list[0] = list2;
                list->list[i] = list3;
            }
        }
    }
#ifdef BEN_DEBUG
    assert(k == TRUE);
#endif

    //destroy stuff
    return list;
}

void writeOutCactusGraph(struct CactusGraph *cactusGraph,
        struct PinchGraph *pinchGraph, FILE *fileHandle) {
    assert(pinchGraph != NULL);
    /*
     * Writes out a graph in 'dot' format, compatible with graphviz.
     *
     * The associated function 'writeOutEdgePieces' gives a way of decoding
     * the pieces associated with the black (piece containing) edges
     * of the graph.
     */
    int32_t i, j; //, k;
    struct CactusVertex *vertex;
    struct CactusEdge *edge;
    //struct PinchEdge *pinchEdge;
    //struct Piece *piece;

    //Write the preliminaries.
    fprintf(fileHandle, "graph G {\n");
    fprintf(fileHandle, "overlap=false\n");
    //Write the vertices.
    for (i = 0; i < cactusGraph->vertices->length; i++) {
        vertex = cactusGraph->vertices->list[i];
#ifdef BEN_DEBUG
        assert(vertex->vertexID == i);
#endif
        fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n",
                vertex->vertexID, vertex->vertexID);
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
                /*if(edge->pieces->length > 0) {
                 for(k=0; k<edge->pieces->length; k++) {
                 piece = edge->pieces->list[k];
                 pinchEdge = getContainingBlackEdge(pinchGraph, piece->contig, piece->start);
                 fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n [label=\"" INT_STRING ":" INT_STRING ":%s\"];\n",
                 edge->from->vertexID, edge->to->vertexID, piece->start, piece->end, netMisc_nameToStringStatic(piece->contig));
                 }
                 }
                 else {*/
                fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n",
                        edge->from->vertexID, edge->to->vertexID);
                //}
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

void computeCactusGraph(struct PinchGraph *pinchGraph,
        struct CactusGraph **cactusGraph,
        struct List **threeEdgeConnectedComponents) {
    struct PinchVertex *vertex;
    stList *list;
    int32_t i, j;
    struct List *greyEdgeComponents, *list2;
    stList *vertices;

    ///////////////////////////////////////////////////////////////////////////
    // Run the three-edge connected component algorithm to identify
    // three edge connected components.
    ///////////////////////////////////////////////////////////////////////////

    greyEdgeComponents = getRecursiveComponents(pinchGraph, NULL);

    vertices = writeOut3EdgeGraph(pinchGraph, greyEdgeComponents);

    list = computeThreeEdgeConnectedComponents(vertices);
    st_logInfo("Seems to have successfully run the three edge command: %i\n",
            stList_length(list));

    //Parse results (the three edge connected components).
    *threeEdgeConnectedComponents = readThreeEdgeComponents(pinchGraph,
            greyEdgeComponents, list);
    st_logInfo("Read in the three edge components\n");
    stList_destruct(list);

    for (i = 0; i < (*threeEdgeConnectedComponents)->length; i++) {
        list2 = (*threeEdgeConnectedComponents)->list[i];
        st_logDebug("3 edge component : " INT_STRING " ", i);
        for (j = 0; j < list2->length; j++) {
            vertex = list2->list[j];
            st_logDebug(" vertex, " INT_STRING " ", vertex->vertexID);
        }
        st_logDebug("\n");
    }
    //Cleanup the three edge input/output fil
    destructList(greyEdgeComponents);
    stList_destruct(vertices);

    ///////////////////////////////////////////////////////////////////////////
    // Merge vertices in each three edge connected component to convert the main graph
    // into a 'cactus' graph.
    ///////////////////////////////////////////////////////////////////////////

    //Collapse the graph to cactus tree
    *cactusGraph = constructCactusGraph(pinchGraph,
            *threeEdgeConnectedComponents);
    checkCactusGraph(pinchGraph, *threeEdgeConnectedComponents, *cactusGraph);
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

void cactusVertex_merge(struct CactusGraph *cactusGraph,
        struct CactusVertex *one, struct CactusVertex *two) {
    /*
     * Merges two into one, destroying two and adding to one.
     */
    assert(one != two);
    assert(one->vertexID != two->vertexID);
    //redirect all edges to one.
    int32_t i;
    for (i = 0; i < two->edges->length; i++) {
        struct CactusEdge *edge = two->edges->list[i];
#if BEN_ULTRA_DEBUG
        assert(edge->from == two);
        assert(edge->rEdge->to == two);
        assert(!listContains(one->edges, edge));
#endif
        edge->from = one;
        edge->rEdge->to = one;
        listAppend(one->edges, edge);
    }
    //switch the vertex id's of the highest cactus vertex and two.
    struct CactusVertex *three =
            cactusGraph->vertices->list[cactusGraph->vertices->length - 1];
    assert(three->vertexID == cactusGraph->vertices->length-1);
    three->vertexID = two->vertexID;
    assert(cactusGraph->vertices->list[two->vertexID] == two);
    cactusGraph->vertices->list[two->vertexID] = three;
    cactusGraph->vertices->length--;

    //now destruct two, but not its list of edges.
    two->edges->destructElement = NULL;
    destructCactusVertex(two);
}

void circulariseStemsP(struct CactusGraph *cactusGraph,
        struct CactusEdge *edge, struct CactusVertex *sourceVertex,
        struct hashtable *stemHash, struct hashtable *seen,
        struct List *edgesToMerge) {
    int32_t i, j;
    struct CactusEdge *edge2;

    if (hashtable_search(seen, edge->to) == NULL) {
        hashtable_insert(seen, edge->to, edge->to);
        if (hashtable_search(stemHash, edge) != NULL) { //is a stem
            j = 0;
            for (i = 0; i < edge->to->edges->length; i++) {
                edge2 = edge->to->edges->list[i];
                if (hashtable_search(stemHash, edge) != NULL) {//is a stem
                    j++;
                }
            }
            assert(j> 0);
            if (j != 2) { //if is either branch or leaf.
                listAppend(edgesToMerge, sourceVertex);
                listAppend(edgesToMerge, edge->to);
                sourceVertex = edge->to;
            }
        } else {
            sourceVertex = edge->to;
        }
        for (i = 0; i < edge->to->edges->length; i++) { //call recursively.
            edge2 = edge->to->edges->list[i];
            circulariseStemsP(cactusGraph, edge2, sourceVertex, stemHash, seen,
                    edgesToMerge);
        }
    }
}

void circulariseStems(struct CactusGraph *cactusGraph) {
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
    st_logDebug("Constructed the biconnected components for making the net\n");

    ////////////////////////////////////////////////
    //(2) Put stems in a hash
    ////////////////////////////////////////////////

    stemHash = create_hashtable(cactusGraph->vertices->length * 2,
            hashtable_key, hashtable_equalKey, NULL, NULL);
    for (i = 0; i < biConnectedComponents->length; i++) {
        biConnectedComponent = biConnectedComponents->list[i];
        if (biConnectedComponent->length == 1) {
            edge = biConnectedComponent->list[0];
            if (edge->from != edge->to) { //must be a stem as not a self loop
                hashtable_insert(stemHash, edge, edge);
                hashtable_insert(stemHash, edge->rEdge, edge->rEdge);
            }
        }
    }

    st_logDebug("Put the stems in a hash\n");

    ////////////////////////////////////////////////
    //(3) Do DFS on graph to find vertex merges that we wish to make
    ////////////////////////////////////////////////

    seen = create_hashtable(cactusGraph->vertices->length * 2, hashtable_key,
            hashtable_equalKey, NULL, NULL);
    vertex = cactusGraph->vertices->list[0];
    hashtable_insert(seen, vertex, vertex);
    cactusGraphMerges = constructEmptyList(0, NULL);
    for (i = 0; i < vertex->edges->length; i++) {
        circulariseStemsP(cactusGraph, vertex->edges->list[i], vertex,
                stemHash, seen, cactusGraphMerges);
    }
    st_logDebug("Done the DFS\n");

    ////////////////////////////////////////////////
    //(4) Do vertex merges.
    ////////////////////////////////////////////////

    assert((cactusGraphMerges->length % 2) == 0);
    for (i = cactusGraphMerges->length - 1; i >= 0; i -= 2) {
        assert(i> 0);
        cactusVertex_merge(cactusGraph, cactusGraphMerges->list[i - 1],
                cactusGraphMerges->list[i]);
    }

    ////////////////////////////////////////////////
    //(5) Cleanup
    ////////////////////////////////////////////////

    hashtable_destroy(stemHash, FALSE, FALSE);
    hashtable_destroy(seen, FALSE, FALSE);
    destructList(biConnectedComponents);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Choosing which blocks in the chains to keep.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


float treeCoverage2(struct CactusEdge *cactusEdge, Net *net,
        struct PinchGraph *pinchGraph) {
    /*
     * Returns the proportion of the tree covered by the block.
     */
#ifdef BEN_DEBUG
    assert(!isAStubCactusEdge(cactusEdge, pinchGraph));
#endif
    return treeCoverage(
            cactusEdgeToFirstPinchEdge(cactusEdge, pinchGraph)->from, net);
}

int32_t chainLength(struct List *biConnectedComponent, int32_t includeStubs,
        struct PinchGraph *pinchGraph) {
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

int32_t chainBaseLength(struct List *biConnectedComponent,
        struct PinchGraph *pinchGraph) {
    /*
     * Get the number of links in the chain.
     */
    int32_t i, j;
    struct CactusEdge *cactusEdge;
    struct Piece *piece;

    i = 0;
    for (j = 0; j < biConnectedComponent->length; j++) {
        cactusEdge = biConnectedComponent->list[j];
        if (!isAStubCactusEdge(cactusEdge, pinchGraph)) {
            piece = cactusEdge->pieces->list[0];
            i += piece->end - piece->start + 1;
        }
    }
    st_uglyf("I have a chain with %i blocks and %i lengths\n",
            biConnectedComponent->length, i);
    return i;
}

stSortedSet *filterBlocksByTreeCoverageAndLength(
        struct List *biConnectedComponents, Net *net,
        float minimumTreeCoverage, /*Minimum tree coverage to be included (>=) */
        int32_t minimumBlockDegree, /*The minimum number of segments in a block to be included (>=)*/
        int32_t minimumBlockLength, /*The minimum length of an block to be included (>=)*/
        int32_t minimumChainLength, /* Minimum chain length to be included (>=)*/
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
                    assert(cactusEdge->pieces->length> 0);
                    if(minimumBlockDegree <= 0 || cactusEdge->pieces->length >= minimumBlockDegree) {
                        if(minimumTreeCoverage <= 0.0 || treeCoverage2(cactusEdge, net, pinchGraph) >= minimumTreeCoverage) {
                            struct Piece *piece = cactusEdge->pieces->list[0];
                            if (minimumBlockLength <= 0 || piece->end - piece->start + 1 >= minimumBlockLength) {
                                stSortedSet_insert(chosenBlocks, cactusEdge);
                            }
                        }
                    }
                }
            }
        }
    }
    return chosenBlocks;
}

void logTheChosenBlockSubset(struct List *biConnectedComponents,
        stSortedSet *chosenBlocks, struct PinchGraph *pinchGraph, Net *net) {
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
        totalBaseLengthOfAllBlocks += chainBaseLength(biConnectedComponent,
                pinchGraph);
        totalNumberOfAllBlocks += chainLength(biConnectedComponent,
                FALSE, pinchGraph);
        totalNumberOfStubBlocks += chainLength(biConnectedComponent,
                TRUE, pinchGraph) - chainLength(biConnectedComponent,
                FALSE, pinchGraph);
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
            totalBlockScore += treeCoverage2(cactusEdge, net, pinchGraph);
            piece = cactusEdge->pieces->list[0];
            totalBlockLength += piece->end - piece->start + 1;
            averagePieceNumber += cactusEdge->pieces->length;
            j++;
        }
    }
    stSortedSet_destructIterator(it);
    st_logInfo(
            "Chosen block subset composed of %i blocks, of average length %f and average tree coverage %f, total base length of all blocks: %f, total number of all blocks %f, average length of all blocks: %f, total number of stub blocks: %f, average piece number of chosen blocks: %f, average piece number of all blocks: %f\n",
            stSortedSet_size(chosenBlocks), totalBlockLength / j, totalBlockScore / j,
            totalBaseLengthOfAllBlocks, totalNumberOfAllBlocks,
            totalBaseLengthOfAllBlocks / totalNumberOfAllBlocks,
            totalNumberOfStubBlocks, averagePieceNumber / j,
            averagePieceNumberOfAllBlocks / totalNumberOfAllBlocks);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Cactus graph misc functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchEdge *cactusEdgeToFirstPinchEdge(struct CactusEdge *edge,
        struct PinchGraph *pinchGraph) {
    struct Piece *piece;
#ifdef BEN_DEBUG
    assert(edge->pieces->length> 0);
#endif
    piece = edge->pieces->list[0];
    return getContainingBlackEdge(pinchGraph, piece->contig, piece->start);
}

int32_t isAStubCactusEdge(struct CactusEdge *edge,
        struct PinchGraph *pinchGraph) {
    return isAStub(cactusEdgeToFirstPinchEdge(edge, pinchGraph));
}

struct hashtable *createHashColouringPinchEdgesByChains(
        struct PinchGraph *pinchGraph, struct List *biConnectComponentsList) {
    struct List *biConnectedComponent;
    struct CactusEdge *cactusEdge;
    struct Piece *piece;
    int32_t i, j, k;

    //Put the chain pieces in a hash to colour the black edges of the pinch graph.
    struct hashtable *hash = create_hashtable(
            pinchGraph->vertices->length * 10, hashtable_key,
            hashtable_equalKey, NULL, (void(*)(void *)) destructInt);

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

