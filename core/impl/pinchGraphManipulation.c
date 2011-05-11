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
#include <float.h>

#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "ctype.h"
#include "adjacencyComponents.h"
#include "pinchGraph.h"
#include "pinchGraphManipulation.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Procedures for removing homology between over aligned edges.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

//for each vertex with black edge degree greater than X add to list + hash.
//for each member of this list:

//remove each edge from vertex.
//if edge is followed by an edge starting from a another high degree vertex or a vertex with only one black edge, then
//	merge the edge.
//else:
//create a new vertex and rejoin the edge.

void removeTrivialGreyEdge(struct PinchGraph *graph,
        struct PinchVertex *vertex1, struct PinchVertex *vertex2,
        Flower *flower) {
#ifdef BEN_DEBUG
    assert(lengthBlackEdges(vertex1) == lengthBlackEdges(vertex2));
    assert(lengthGreyEdges(vertex1) == 1);
    assert(lengthGreyEdges(vertex2) == 1);
    assert(getFirstGreyEdge(vertex1) == vertex2);
    assert(getFirstGreyEdge(vertex2) == vertex1);
#endif

    //if(lengthBlackEdges(vertex1) > 1) {
    //	assert(FALSE);
    //}

    //For each black edge to vertex1 find consecutive edge from vertex2, then join.
    while (lengthBlackEdges(vertex1) > 0) {
        assert(lengthBlackEdges(vertex1) == lengthBlackEdges(vertex2));

        struct PinchEdge *edge1 = getFirstBlackEdge(vertex1);
#ifdef BEN_DEBUG
        assert(edge1 != NULL);
        assert(!isAStub(edge1));
#endif
        edge1 = edge1->rEdge;
#ifdef BEN_DEBUG
        assert(edge1->to == vertex1);
#endif
        //first find the grey edge to attach to the new vertex we're about to create

        struct PinchEdge *edge2 = getNextEdge(graph, edge1, flower);
#ifdef BEN_DEBUG
        assert(edge2 != NULL);
        assert(!isAStub(edge2));
        assert(edge2->from == vertex2);
#endif

        struct PinchEdge *edge3 = constructPinchEdge(constructPiece(
                edge1->piece->contig, edge1->piece->start, edge2->piece->end));
        connectPinchEdge(edge3, edge1->from, edge2->to);

        //Remove the old edges
        removePinchEdgeFromGraphAndDestruct(graph, edge1);
        removePinchEdgeFromGraphAndDestruct(graph, edge2);

        //Add the new pinch edge to the graph after removing the old edges from the graph.
        addPinchEdgeToGraph(graph, edge3);
    }

    //Destruct the old vertices.
#ifdef BEN_DEBUG
    assert(lengthBlackEdges(vertex1) == 0);
    assert(lengthBlackEdges(vertex2) == 0);
#endif
    removeVertexFromGraphAndDestruct(graph, vertex1);
    removeVertexFromGraphAndDestruct(graph, vertex2);
}

void removeTrivialGreyEdgeComponents(struct PinchGraph *graph,
        struct List *listOfVertices, Flower *flower) {
    /*
     * Finds cases where two vertices are linked by adjacency, and have no other adjacencies,
     * to remove them from the graph.
     */
    struct List *list;
    int32_t i;
    struct PinchVertex *vertex1;
    struct PinchVertex *vertex2;
    struct PinchEdge *edge1;
    struct PinchEdge *edge2;

    //Build the list of trivial components.
    list = constructEmptyList(0, NULL);
    for (i = 0; i < listOfVertices->length; i++) {
        vertex1 = listOfVertices->list[i];
        if (lengthGreyEdges(vertex1) == 1 && lengthBlackEdges(vertex1) > 0) {
            edge1 = getFirstBlackEdge(vertex1);
            vertex2 = getFirstGreyEdge(vertex1);
            if (lengthGreyEdges(vertex2) == 1 && lengthBlackEdges(vertex2) > 0) {
                edge2 = getFirstBlackEdge(vertex2);
                if (!isAStub(edge1) && !isAStub(edge2)) {
                    if (vertex1->vertexID < vertex2->vertexID) { //Avoid treating self loops (equal) and dealing with trivial grey components twice.
                        listAppend(list, vertex1);
                    }
                }
            }
        }
    }

    //Remove the trivial components.
    for (i = 0; i < list->length; i++) {
        vertex1 = list->list[i];
        vertex2 = getFirstGreyEdge(vertex1);
        removeTrivialGreyEdge(graph, vertex1, vertex2, flower);
    }

    //cleanup
    destructList(list);
}

void splitMultipleBlackEdgesFromVertex(struct PinchGraph *pinchGraph,
        struct PinchVertex *vertex, struct List *newVerticesList,
        Flower *flower) {
    /*
     * Splits multiple black edges from the vertex, so that vertex is incidental with only one black and grey
     * edge.
     */
    int32_t j;
    struct PinchEdge *edge;
    struct PinchVertex *vertex2;
    struct PinchVertex *vertex3;
    struct List *list;

#ifdef BEN_DEBUG
    assert(vertex == pinchGraph->vertices->list[vertex->vertexID]);
    assert(lengthBlackEdges(vertex) > 0);
    assert(!vertex_isDeadEnd(vertex));
    assert(!vertex_isEnd(vertex));
#endif

    list = constructEmptyList(0, NULL);
    while (lengthBlackEdges(vertex) > 0) {
        edge = getFirstBlackEdge(vertex);
        //first find the grey edge to attach to the new vertex we're about to create
        vertex3 = getNextEdge(pinchGraph, edge->rEdge, flower)->from;
        listAppend(list, vertex3); //can't detach the old vertices yet

        assert(popBlackEdge(vertex) == edge); //detaches edge from vertex.
#ifdef BEN_DEBUG
        assert(!isAStub(edge));
#endif
        //make a new vertex
        vertex2 = constructPinchVertex(pinchGraph, -1, 0, 0);
        listAppend(newVerticesList, vertex2);

        //attach the new vertex to the black edges.
        edge->from = vertex2;
        edge->rEdge->to = vertex2;
        insertBlackEdge(vertex2, edge);

        //finally connect the two new vertices.
        connectVertices(vertex2, vertex3);
    }
    for (j = 0; j < list->length; j++) {
        vertex3 = list->list[j];
        if (containsGreyEdge(vertex3, vertex)) { //it may have already been detached.
            removeGreyEdge(vertex3, vertex);
        }
    }
    //now remove the old vertex
    removeVertexFromGraphAndDestruct(pinchGraph, vertex);
    destructList(list);
}

void removeOverAlignedEdges_P(struct PinchVertex *vertex,
        int32_t extensionSteps, struct List *list, struct hashtable *hash) {
    void *greyEdgeIterator = getGreyEdgeIterator(vertex);
    struct PinchVertex *vertex2;
    struct PinchVertex *vertex3;

    int32_t distance = *(int32_t *) hashtable_search(hash, vertex);
    if (distance < extensionSteps) {
        while ((vertex2 = getNextGreyEdge(vertex, greyEdgeIterator)) != NULL) {
            if (lengthBlackEdges(vertex2) > 0) {
                struct PinchEdge *edge = getFirstBlackEdge(vertex2);
                int32_t length = edge->piece->end - edge->piece->start + 1;
                if (!isAStub(edge)) {
                    int32_t *i = hashtable_search(hash, vertex2);
                    vertex3 = edge->to;
                    int32_t *j = hashtable_search(hash, vertex3);
                    if (i == NULL) {
                        assert(j == NULL);
                        listAppend(list,
                                vertex2->vertexID > vertex3->vertexID ? vertex3
                                        : vertex2);
                        hashtable_insert(hash, vertex2, constructInt(distance));
                        hashtable_insert(hash, vertex3, constructInt(distance
                                + length));
                    } else {
                        *i = *i > distance ? distance : *i;
                        assert(j != NULL);
                        *j = *j > distance + length ? distance + length : *j;
                    }
                }
            }
        }
    }
    destructGreyEdgeIterator(greyEdgeIterator);
}

void removeOverAlignedEdges(struct PinchGraph *pinchGraph,
        float minimumTreeCoverage, int32_t maxDegree,
        struct List *extraEdgesToUndo, int32_t extensionSteps, Flower *flower) {
    /*
     * Method splits black edges from the graph with degree higher than a given number of sequences.
     */
    int32_t i, j, k;
    struct List *list;
    struct List *list2;
    struct PinchEdge *edge;
    struct PinchVertex *vertex;
    struct PinchVertex *vertex2;
    struct hashtable *hash;

    list = constructEmptyList(0, NULL);

    /*
     * Finds the list of edges to undo (those which meet the criteria).
     */
    hash = create_hashtable(0, hashtable_key, hashtable_equalKey, NULL, free);
    for (i = 0; i < pinchGraph->vertices->length; i++) {
        vertex = pinchGraph->vertices->list[i];
        if (lengthBlackEdges(vertex) >= 1
                && !isAStub(getFirstBlackEdge(vertex))) {
            if (lengthBlackEdges(vertex) > maxDegree
                    || (minimumTreeCoverage > 0.0 && treeCoverage(vertex,
                            flower) < minimumTreeCoverage)) { //has a high degree and is not a stub/cap
                vertex2 = getFirstBlackEdge(vertex)->to;
                if (vertex->vertexID < vertex2->vertexID) {
                    hashtable_insert(hash, vertex, constructInt(0));
                    hashtable_insert(hash, vertex2, constructInt(0));
                    listAppend(list, vertex);
                }
            }
        }
    }

    /*
     * This adds a bunch of extra edges to the list which should be undone. It ignores stub edges
     * and duplicates.
     */
    if (extraEdgesToUndo != NULL) {
        for (i = 0; i < extraEdgesToUndo->length; i++) {
            edge = extraEdgesToUndo->list[i];
            if (!isAStub(edge)) {
                if (edge->from->vertexID > edge->to->vertexID) {
                    edge = edge->rEdge;
                }
                if (hashtable_search(hash, edge->from) == NULL) {
                    assert(hashtable_search(hash, edge->to) == NULL);
                    hashtable_insert(hash, edge->from, constructInt(0));
                    hashtable_insert(hash, edge->to, constructInt(0));
                    listAppend(list, edge->from);
                } else {
                    assert(hashtable_search(hash, edge->to) != NULL);
                }
            }
        }
    }

    st_logDebug(
            "Got the initial list of over-aligned black edges to undo, total: %i\n",
            list->length);

    if (extensionSteps > 0) {
        i = 0, k = 10;
        while (list->length != i || k-- > 0) { //k term to ensure distances have been propagated
            assert(list->length >= i);
            i = list->length;
            list2 = listCopy(list); //just use the vertices in the existing list
            for (j = 0; j < list2->length; j++) {
                vertex = list2->list[j];
                vertex2 = getFirstBlackEdge(vertex)->to;
                removeOverAlignedEdges_P(vertex, extensionSteps, list, hash);
                removeOverAlignedEdges_P(vertex2, extensionSteps, list, hash);
            }
            destructList(list2);
        }
    }

    //now remove all single black edge connected vertices
    list2 = constructEmptyList(0, NULL);
    for (i = 0; i < list->length; i++) {
        vertex = list->list[i];
        if (lengthBlackEdges(vertex) > 1) {
            listAppend(list2, vertex);
        }
    }
    destructList(list);
    list = list2;

    st_logDebug("Got the list of black edges to undo, total length: %i!\n",
            list->length);

    list2 = constructEmptyList(0, NULL);
    for (i = 0; i < list->length; i++) {
        vertex = list->list[i];
        vertex2 = getFirstBlackEdge(vertex)->to;
        list2->length = 0;
        splitMultipleBlackEdgesFromVertex(pinchGraph, vertex, list2, flower);
        splitMultipleBlackEdgesFromVertex(pinchGraph, vertex2, list2, flower);
        removeTrivialGreyEdgeComponents(pinchGraph, list2, flower); //now get rid of any trivial components
    }

    destructList(list);
    destructList(list2);
    hashtable_destroy(hash, 1, 0);
}

void trimEdges(struct PinchGraph *pinchGraph, int32_t trim, Flower *flower) {
    assert(trim >= 0);

    if (trim == 0) {
        return;
    }

    /*
     * Get the list of edges to be trimmed.
     */
    stList *listOfEdges = stList_construct();
    for (int32_t i = 0; i < pinchGraph->vertices->length; i++) {
        struct PinchVertex *vertex = pinchGraph->vertices->list[i];
        if (lengthBlackEdges(vertex) >= 1
                && !isAStub(getFirstBlackEdge(vertex))) {
            struct PinchEdge *edge = getFirstBlackEdge(vertex);
            assert(edge->piece->end - edge->piece->start >= 0);
            struct PinchVertex *vertex2 = edge->to;
            if (vertex->vertexID < vertex2->vertexID) {
                if (edge->piece->end - edge->piece->start + 1 > 2 * trim) { //When trimmed there will be a block left of length 1 or greater.
                    stList_append(listOfEdges, edge);
                }
            }
        }
    }

    /*
     * Now trim the edges.
     */
    stHash *vertexToAdjacencyComponent = stHash_construct();
    for (int32_t i = 0; i < stList_length(listOfEdges); i++) {
        struct PinchEdge *edge = stList_get(listOfEdges, i);
        int32_t position1 = edge->piece->start + trim;
        int32_t position2 = edge->piece->end - trim;
        assert(position1 <= position2); //This must be the case as there must be a block of length 1 or greater in between the two trimmed regions.
        Name contig = edge->piece->contig;
        splitEdge(pinchGraph, contig, position1, LEFT,
                vertexToAdjacencyComponent);
        splitEdge(pinchGraph, contig, position2, RIGHT,
                vertexToAdjacencyComponent);
        edge = getContainingBlackEdge(pinchGraph, contig, position1);
//#ifdef BEN_DEBUG
        assert(edge != NULL);
        assert(edge->piece->contig == contig);
        assert(edge->piece->start == position1);
        assert(edge->piece->end == position2);
//#endif
        stList_set(listOfEdges, i, edge); //We keep just the middle bit of the edge.
    }
    stHash_destruct(vertexToAdjacencyComponent);

    /*
     * Get the list of vertices incident with black edges we won't remove.
     */
    stSortedSet *vertices = stSortedSet_construct();
    for (int32_t i = 0; i < stList_length(listOfEdges); i++) {
        struct PinchEdge *edge = stList_get(listOfEdges, i);
//#ifdef BEN_DEBUG
        assert(stSortedSet_search(vertices, edge->from) == NULL);
        assert(stSortedSet_search(vertices, edge->to) == NULL);
        assert(edge->from != edge->to);
//#endif
        stSortedSet_insert(vertices, edge->from);
        stSortedSet_insert(vertices, edge->to);
    }

    /*
     * Convert the list of vertices that
     */
    struct List *listOfEdgesToRemove = constructEmptyList(0, NULL);
    for (int32_t i = 0; i < pinchGraph->vertices->length; i++) {
        struct PinchVertex *vertex = pinchGraph->vertices->list[i];
        if (lengthBlackEdges(vertex) >= 1
                && !isAStub(getFirstBlackEdge(vertex)) && stSortedSet_search(
                vertices, vertex) == NULL) {//is not a stub, is of degree greater than 1 and is not the middle of one of the trimmed edges (so it must be removed)
//#ifdef BEN_DEBUG
            struct PinchEdge *edge = getFirstBlackEdge(vertex);
            //st_uglyf("I got %i %i %i %i\n", edge->piece->end, edge->piece->start, edge->piece->end - edge->piece->start + 1, trim);
            assert(edge->piece->end - edge->piece->start + 1 <= 2 * trim);
//#endif
            listAppend(listOfEdgesToRemove, getFirstBlackEdge(vertex));
        }
    }

    /*
     * Finally remove all the trimmed edges.
     */
    st_logInfo(
            "Before trimming the graph, the graph has %i vertices and %i black edges\n",
            pinchGraph->vertices->length, avl_count(pinchGraph->edges));
    removeOverAlignedEdges(pinchGraph, 0.0, INT32_MAX, listOfEdgesToRemove, 0,
            flower);
    st_logInfo(
            "After trimming the graph, the graph has %i vertices and %i black edges\n",
            pinchGraph->vertices->length, avl_count(pinchGraph->edges));

    removeTrivialGreyEdgeComponents(pinchGraph, pinchGraph->vertices, flower); //remove any pointless adjacencies.
    st_logInfo(
            "After removing the trivial graph components the graph has %i vertices and %i black edges\n",
            pinchGraph->vertices->length, avl_count(pinchGraph->edges));
    checkPinchGraph(pinchGraph);

    /*
     * Cleanup remaining
     */
    stList_destruct(listOfEdges);
    destructList(listOfEdgesToRemove);
    stSortedSet_destruct(vertices);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method for linking the stub components to the
//sink component.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

bool linkStubComponentsToTheSinkComponent_passThroughFn(struct PinchEdge *edge) {
    assert(edge != NULL);
    return 1;
}

void linkStubComponentsToTheSinkComponent(struct PinchGraph *pinchGraph,
        Flower *flower, int32_t attachEnds) {
    struct PinchVertex *vertex;
    struct PinchVertex *sinkVertex;
    int32_t i, k, l;
    struct PinchEdge *edge;
    Sequence *sequence;
    Sequence *longestSequence;
    Cap *cap;

    //isolate the separate graph components using the components method
    stList *adjacencyComponents = getAdjacencyComponents2(pinchGraph,
            linkStubComponentsToTheSinkComponent_passThroughFn);

    sinkVertex = pinchGraph->vertices->list[0];
    //for each non-sink component select a random stub to link to the sink vertex.
    k = 0;
    l = 0;
    for (i = 0; i < stList_length(adjacencyComponents); i++) {
        stSortedSet *adjacencyComponent = stList_get(adjacencyComponents, i);
        assert(stSortedSet_size(adjacencyComponent) > 0);
        if (stSortedSet_search(adjacencyComponent, sinkVertex) == NULL) {
            //Get the longest sequence contained in the component and attach
            //its two ends to the source vertex.
            //Make the the two ends attached end_makeAttached(end) / end_makeUnattached(end)
            longestSequence = NULL;
            stSortedSetIterator *it = stSortedSet_getIterator(
                    adjacencyComponent);
            while ((vertex = stSortedSet_getNext(it)) != NULL) {
                if (vertex_isDeadEnd(vertex)) {
                    assert(lengthGreyEdges(vertex) == 0);
                    assert(lengthBlackEdges(vertex) == 1);
                    edge = getFirstBlackEdge(vertex);
                    cap = flower_getCap(flower, edge->piece->contig);
                    assert(cap != NULL);
                    sequence = cap_getSequence(cap);
                    assert(sequence != NULL);
                    if (longestSequence == NULL || sequence_getLength(sequence)
                            > sequence_getLength(longestSequence)) {
                        longestSequence = sequence;
                    }
                }
            }
            stSortedSet_destructIterator(it);
            //assert(0);
            assert(longestSequence != NULL);
            it = stSortedSet_getIterator(adjacencyComponent);
            while ((vertex = stSortedSet_getNext(it)) != NULL) {
                if (vertex_isDeadEnd(vertex)) {
                    assert(lengthGreyEdges(vertex) == 0);
                    assert(lengthBlackEdges(vertex) == 1);
                    edge = getFirstBlackEdge(vertex);
                    cap = flower_getCap(flower, edge->piece->contig);
                    assert(cap != NULL);
                    End *end = cap_getEnd(cap);
                    assert(end_isStubEnd(end));
                    assert(end_isFree(end));
                    sequence = cap_getSequence(cap);
                    assert(sequence != NULL);
                    if (sequence == longestSequence) {
                        if (!flower_hasParentGroup(flower) && attachEnds) {
                            end_makeAttached(end);
                        }
                        connectVertices(vertex, sinkVertex);
                        k++;
                    }
                }
            }
            stSortedSet_destructIterator(it);
        }
    }

#ifdef BEN_DEBUG
    assert(k == 2*(stList_length(adjacencyComponents)-1));
#endif

    //clean up
    stList_destruct(adjacencyComponents);
}

void unlinkStubComponentsFromTheSinkComponent(struct PinchGraph *pinchGraph,
        Flower *flower) {
    for (int32_t i = 0; i < pinchGraph->vertices->length; i++) {
        struct PinchVertex *vertex = pinchGraph->vertices->list[i];
        if (vertex_isDeadEnd(vertex)) {
            assert(lengthBlackEdges(vertex) >= 1);
            struct PinchEdge *pinchEdge = getFirstBlackEdge(vertex);
            Cap *cap = flower_getCap(flower, pinchEdge->piece->contig);
            assert(cap != NULL);
            End *end = cap_getEnd(cap);
            assert(end_isStubEnd(end));
            if (end_isStubEnd(end) && end_isFree(end)) {
                if (lengthGreyEdges(vertex) == 1) { //is attached to the origin node.
                    assert(getFirstGreyEdge(vertex)
                            == pinchGraph->vertices->list[0]);
                    disconnectVertices(vertex, pinchGraph->vertices->list[0]);
                } else {
                    assert(lengthGreyEdges(vertex) == 0);
                }
            }
        }
    }
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method for assessing how much of the event tree the
//given set of connected pinch edges covers.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


float treeCoverage(struct PinchVertex *vertex, Flower *flower) {
    /*
     * Returns the proportion of the tree covered by the block.
     */
    struct Piece *piece;
    EventTree *eventTree;
    Event *event;
    Event *commonAncestorEvent;
    struct hashtable *hash;
    float treeCoverage;
    Sequence *sequence;

#ifdef BEN_DEBUG
    assert(lengthBlackEdges(vertex) > 0);
    assert(!isAStub(getFirstBlackEdge(vertex)));
#endif

    eventTree = flower_getEventTree(flower);
    commonAncestorEvent = NULL;
    void *blackEdgeIterator = getBlackEdgeIterator(vertex);
    struct PinchEdge *edge;
    while ((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
        piece = edge->piece;
        sequence = flower_getSequence(flower, piece->contig);
        assert(sequence != NULL);
        event = sequence_getEvent(sequence);
        assert(event != NULL);
        commonAncestorEvent = commonAncestorEvent == NULL ? event
                : eventTree_getCommonAncestor(event, commonAncestorEvent);
    }
    destructBlackEdgeIterator(blackEdgeIterator);
    assert(commonAncestorEvent != NULL);
    treeCoverage = 0.0;
    hash = create_hashtable(eventTree_getEventNumber(eventTree) * 2,
            hashtable_key, hashtable_equalKey, NULL, NULL);

    blackEdgeIterator = getBlackEdgeIterator(vertex);
    while ((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
        piece = edge->piece;
        sequence = flower_getSequence(flower, piece->contig);
        assert(sequence != NULL);
        event = sequence_getEvent(sequence);
        assert(event != NULL);
        while (event != commonAncestorEvent && hashtable_search(hash, event)
                == NULL) {
            treeCoverage += event_getBranchLength(event);
            hashtable_insert(hash, event, event);
            event = event_getParent(event);
#ifdef BEN_DEBUG
            assert(event != NULL);
#endif
        }
    }
    destructBlackEdgeIterator(blackEdgeIterator);
    hashtable_destroy(hash, FALSE, FALSE);
    float wholeTreeCoverage = event_getSubTreeBranchLength(event_getChild(
            eventTree_getRootEvent(eventTree), 0));
    assert(wholeTreeCoverage >= 0.0);
    if (wholeTreeCoverage <= 0.0) { //deal with case all leaf branches are not empty.
        return 0.0;
    }
    treeCoverage /= wholeTreeCoverage;
    assert(treeCoverage >= -0.001);
    assert(treeCoverage <= 1.0001);
    return treeCoverage;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for pinching the graph.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct VertexChain {
    /*
     * A holder type for getting chains of vertices, see getChainOfVertices().
     */
    struct List *listOfVertices;
    struct IntList *coordinates;
    struct IntList *leftsOrRights;
};

struct VertexChain *constructVertexChain() {
    struct VertexChain *vertexChain;

    vertexChain = st_malloc(sizeof(struct VertexChain));
    vertexChain->listOfVertices = constructEmptyList(0, NULL);
    vertexChain->coordinates = constructEmptyIntList(0);
    vertexChain->leftsOrRights = constructEmptyIntList(0);

    return vertexChain;
}

void resetVertexChain(struct VertexChain *vertexChain) {
    vertexChain->listOfVertices->length = 0;
    vertexChain->coordinates->length = 0;
    vertexChain->leftsOrRights->length = 0;
}

void getChainOfVertices(struct VertexChain *vertexChain,
        struct PinchGraph *graph, struct Piece *piece,
        stHash *vertexToAdjacencyComponents) {
    struct PinchVertex *vertex;
    struct PinchEdge *edge;

    resetVertexChain(vertexChain);

    //do any adjustments off the bat
    splitEdge(graph, piece->contig, piece->start, LEFT,
            vertexToAdjacencyComponents);
    splitEdge(graph, piece->contig, piece->end, RIGHT,
            vertexToAdjacencyComponents);

    vertex = splitEdge(graph, piece->contig, piece->start, LEFT,
            vertexToAdjacencyComponents);
    intListAppend(vertexChain->coordinates, 0);
    intListAppend(vertexChain->leftsOrRights, LEFT);
    listAppend(vertexChain->listOfVertices, vertex);

    //follow chain to get remaining vertices
    edge = getContainingBlackEdge(graph, piece->contig, piece->start);
    while (edge->piece->end < piece->end) {
        intListAppend(vertexChain->coordinates, edge->piece->end - piece->start);
        intListAppend(vertexChain->leftsOrRights, RIGHT);
        listAppend(vertexChain->listOfVertices, edge->to);

        edge = getContainingBlackEdge(graph, piece->contig, edge->piece->end
                + 1);
        intListAppend(vertexChain->coordinates, edge->piece->start
                - piece->start);
        intListAppend(vertexChain->leftsOrRights, LEFT);
        listAppend(vertexChain->listOfVertices, edge->from);
    }

    //add the second vertex.
    vertex = splitEdge(graph, piece->contig, piece->end, RIGHT,
            vertexToAdjacencyComponents);
    intListAppend(vertexChain->coordinates, piece->end - piece->start);
    intListAppend(vertexChain->leftsOrRights, RIGHT);
    listAppend(vertexChain->listOfVertices, vertex);

#ifdef BEN_DEBUG
    //now return the list (the other lists are assigned implicitly)!
    assert(vertexChain->listOfVertices->length == vertexChain->coordinates->length);
    assert(vertexChain->leftsOrRights->length == vertexChain->listOfVertices->length);
#endif
}

int32_t pinchMergePiece_P(struct VertexChain *vertexChain1,
        struct VertexChain *vertexChain2) {
    int32_t i, j, k;

    if (vertexChain1->listOfVertices->length
            != vertexChain2->listOfVertices->length) {
        return FALSE;
    }
    for (i = 0; i < vertexChain1->coordinates->length; i++) {
        j = vertexChain1->coordinates->list[i];
        k = vertexChain2->coordinates->list[i];
        if (j != k) {
            return FALSE;
        }

        j = vertexChain1->leftsOrRights->list[i];
        k = vertexChain2->leftsOrRights->list[i];
#ifdef BEN_DEBUG
        assert(j == LEFT || j == RIGHT);
        assert(k == LEFT || k == RIGHT);
#endif
        if (j != k) {
            return FALSE;
        }
    }
    return TRUE;
}

void pinchMergePiece_getChainOfVertices(struct PinchGraph *graph,
        struct Piece *piece1, struct Piece *piece2,
        struct VertexChain *vertexChain1, struct VertexChain *vertexChain2,
        stHash *vertexToAdjacencyComponents) {
    int32_t i, j, k;
    getChainOfVertices(vertexChain1, graph, piece1, vertexToAdjacencyComponents);
    getChainOfVertices(vertexChain2, graph, piece2, vertexToAdjacencyComponents);

    while (pinchMergePiece_P(vertexChain1, vertexChain2) == FALSE) {
        /*
         * match up the set of vertices for each chain
         */
        for (i = 0; i < vertexChain1->coordinates->length; i++) {
            j = vertexChain1->coordinates->list[i];
            k = vertexChain1->leftsOrRights->list[i];
            //now search if there is an equivalent vertex.
            splitEdge(graph, piece2->contig, piece2->start + j, k,
                    vertexToAdjacencyComponents);
        }
        for (i = 0; i < vertexChain2->coordinates->length; i++) {
            j = vertexChain2->coordinates->list[i];
            k = vertexChain2->leftsOrRights->list[i];
            //now search if there is an equivalent vertex.
            splitEdge(graph, piece1->contig, piece1->start + j, k,
                    vertexToAdjacencyComponents);
        }

        getChainOfVertices(vertexChain1, graph, piece1,
                vertexToAdjacencyComponents);
        getChainOfVertices(vertexChain2, graph, piece2,
                vertexToAdjacencyComponents);
    }

}

struct VertexChain *pMS_vertexChain1 = NULL;
struct VertexChain *pMS_vertexChain2 = NULL;

void pinchMergePiece(struct PinchGraph *graph, struct Piece *piece1,
        struct Piece *piece2, stHash *vertexToAdjacencyComponents) {
    /*
     * Pinches the graph (with the minimum number of required pinches, to
     * represent the contiguous alignment of the two pieces.
     *
     * Pieces have to be of equal length.
     */
    int32_t i, j, k;
    Name contig;
    struct PinchVertex *vertex1;
    struct PinchVertex *vertex2;
    struct PinchVertex *vertex3;
    struct PinchVertex *vertex4;
    struct PinchVertex *vertex5;
    struct PinchEdge *edge;

    if (pMS_vertexChain1 == NULL) {
#ifdef BEN_DEBUG
        assert(pMS_vertexChain2 == NULL);
#endif
        pMS_vertexChain1 = constructVertexChain();
        pMS_vertexChain2 = constructVertexChain();
    }

    /*
     * Check pieces are of the same length (the current (temporary assumption))
     */
#ifdef BEN_DEBUG
    assert(piece1->end - piece1->start == piece2->end - piece2->start);
#endif

    /*
     * run through each chain finding the list of vertices.
     */
    splitEdge(graph, piece1->contig, piece1->start, LEFT,
            vertexToAdjacencyComponents);
    splitEdge(graph, piece1->contig, piece1->end, RIGHT,
            vertexToAdjacencyComponents);
    splitEdge(graph, piece2->contig, piece2->start, LEFT,
            vertexToAdjacencyComponents);
    splitEdge(graph, piece2->contig, piece2->end, RIGHT,
            vertexToAdjacencyComponents);

    pinchMergePiece_getChainOfVertices(graph, piece1, piece2, pMS_vertexChain1,
            pMS_vertexChain2, vertexToAdjacencyComponents);

#ifdef BEN_ULTRA_DEBUG
    //do some debug checks
    assert(pMS_vertexChain1->listOfVertices->length == pMS_vertexChain2->listOfVertices->length);
    assert((pMS_vertexChain1->listOfVertices->length % 2) == 0);
    for(i=0; i<pMS_vertexChain1->listOfVertices->length; i+=2) {
        vertex1 = pMS_vertexChain1->listOfVertices->list[i];
        edge = getFirstBlackEdge(vertex1);
        assert(edge->to == pMS_vertexChain1->listOfVertices->list[i+1]);
    }
    for(i=0; i<pMS_vertexChain2->listOfVertices->length; i+=2) {
        vertex2 = pMS_vertexChain2->listOfVertices->list[i];
        edge = getFirstBlackEdge(vertex2);
        assert(edge->to == pMS_vertexChain2->listOfVertices->list[i+1]);
    }
#endif

    /*
     * Determine if we should proceed with the merge by checking if all the
     * pieces are in the same, else quit.
     */
    for (i = 0; i < pMS_vertexChain1->listOfVertices->length; i++) {
        vertex1 = pMS_vertexChain1->listOfVertices->list[i];
        vertex2 = pMS_vertexChain2->listOfVertices->list[i];
        if (vertex1 == vertex2) {
            continue;
        }
        void *adjacencyComponent1 = stHash_search(vertexToAdjacencyComponents,
                vertex1);
        void *adjacencyComponent2 = stHash_search(vertexToAdjacencyComponents,
                vertex2);
        if (adjacencyComponent1 == NULL || adjacencyComponent2 == NULL
                || adjacencyComponent1 != adjacencyComponent2) {
            return;
        }
    }

    /*
     * Merge the lists of vertices to do the final merge.
     */
    for (i = 0; i < pMS_vertexChain1->listOfVertices->length;) {
        vertex1 = pMS_vertexChain1->listOfVertices->list[i];
        vertex2 = pMS_vertexChain2->listOfVertices->list[i];

#ifdef BEN_DEBUG
        assert(lengthBlackEdges(vertex1) > 0);
#endif
        //check if the two vertices are the ends of same piece.
        edge = getFirstBlackEdge(vertex1);
        if (edge->to == vertex2) {
            //if edge piece is of length greater than one.
            if (edge->piece->end - edge->piece->start > 0) {
                j = (edge->piece->end - edge->piece->start + 1) / 2
                        + edge->piece->start - 1;
                k = 1 + ((edge->piece->end - edge->piece->start + 1) % 2);

                contig = edge->piece->contig;
                vertex4 = splitEdge(graph, contig, j, RIGHT,
                        vertexToAdjacencyComponents);
                vertex5 = splitEdge(graph, contig, j + k, LEFT,
                        vertexToAdjacencyComponents);

#ifdef BEN_DEBUG
                //debug checks
                assert(lengthGreyEdges(vertex4) == 1);
                assert(lengthGreyEdges(vertex5) == 1);
                if (k == 1) {
                    assert(getFirstGreyEdge(vertex4) == vertex5);
                    assert(getFirstGreyEdge(vertex5) == vertex4);
                }
#endif
                /*
                 * The new vertices are not in the chain, so we re parse the vertex chain and start again.
                 */
                pinchMergePiece_getChainOfVertices(graph, piece1, piece2,
                        pMS_vertexChain1, pMS_vertexChain2,
                        vertexToAdjacencyComponents);

                i = 0;
                continue;
            }
            // Else we do nothing, as we can't have self black edges and move one.
        } else {
            if (vertex1 != vertex2) {
                assert(stHash_search(vertexToAdjacencyComponents, vertex1)
                        != NULL);
                assert(stHash_search(vertexToAdjacencyComponents, vertex2)
                        != NULL);
                /*
                 * We have randomly chosen one of the vertex adjacency components..
                 */
                void *adjacencyComponent = stHash_remove(
                        vertexToAdjacencyComponents, vertex1);
                assert(adjacencyComponent != NULL);
                assert(adjacencyComponent == stHash_remove(
                        vertexToAdjacencyComponents, vertex2));
                assert(stHash_search(vertexToAdjacencyComponents, vertex1)
                        == NULL);
                assert(stHash_search(vertexToAdjacencyComponents, vertex2)
                        == NULL);
                vertex3 = mergeVertices(graph, vertex1, vertex2);
                stHash_insert(vertexToAdjacencyComponents, vertex3,
                        adjacencyComponent);
                for (j = i + 1; j < pMS_vertexChain1->listOfVertices->length; j++) {
                    if (pMS_vertexChain1->listOfVertices->list[j] == vertex1
                            || pMS_vertexChain1->listOfVertices->list[j]
                                    == vertex2) {
                        pMS_vertexChain1->listOfVertices->list[j] = vertex3;
                    }
                    if (pMS_vertexChain2->listOfVertices->list[j] == vertex1
                            || pMS_vertexChain2->listOfVertices->list[j]
                                    == vertex2) {
                        pMS_vertexChain2->listOfVertices->list[j] = vertex3;
                    }
                }
            }
        }
        i++;
    }
    //Done the merging of the vertices

#ifdef BEN_ULTRA_DEBUG
    /*
     * Do debug checks that the merge went okay.
     */

    getChainOfVertices(pMS_vertexChain1, graph, piece1);
    getChainOfVertices(pMS_vertexChain2, graph, piece2);

    for(i=0; i<pMS_vertexChain1->listOfVertices->length; i++) {
        vertex1 = pMS_vertexChain1->listOfVertices->list[i];
        void *blackEdgeIterator = getBlackEdgeIterator(vertex1);
        edge = getNextBlackEdge(vertex1, blackEdgeIterator);
        k = edge->piece->end - edge->piece->start;
        vertex2 = edge->to;
        while((edge = getNextBlackEdge(vertex1, blackEdgeIterator)) != NULL) {
            assert(edge->piece->end - edge->piece->start == k);
            assert(edge->to == vertex2);
        }
        destructBlackEdgeIterator(blackEdgeIterator);
    }

    for(i=0; i<pMS_vertexChain2->listOfVertices->length; i++) {
        vertex1 = pMS_vertexChain2->listOfVertices->list[i];
        void *blackEdgeIterator = getBlackEdgeIterator(vertex1);
        edge = getNextBlackEdge(vertex1, blackEdgeIterator);
        k = edge->piece->end - edge->piece->start;
        vertex2 = edge->to;
        while((edge = getNextBlackEdge(vertex1, blackEdgeIterator)) != NULL) {
            assert(edge->piece->end - edge->piece->start == k);
            assert(edge->to == vertex2);
        }
        destructBlackEdgeIterator(blackEdgeIterator);
    }

    assert(pMS_vertexChain1->listOfVertices->length == pMS_vertexChain2->listOfVertices->length);
    for(i=0; i<pMS_vertexChain1->listOfVertices->length; i++) {
        vertex1 = pMS_vertexChain1->listOfVertices->list[i];
        vertex2 = pMS_vertexChain2->listOfVertices->list[i];
        if(vertex1 != vertex2) {
            edge = getFirstBlackEdge(vertex1);
            assert(edge->to == vertex2);
            edge = getFirstBlackEdge(vertex2);
            assert(edge->to == vertex1);
        }
        j = pMS_vertexChain1->coordinates->list[i];
        k = pMS_vertexChain2->coordinates->list[i];
        assert(j == k);

        j = pMS_vertexChain1->leftsOrRights->list[i];
        k = pMS_vertexChain2->leftsOrRights->list[i];
        assert(j == k);
        assert(j == LEFT || j == RIGHT);
    }
#endif
}

int32_t pinchMerge_getContig(char *contig, int32_t start,
        struct hashtable *contigStringToContigIndex) {
    int32_t i, j, k;
    int32_t *iA;
    struct List *list;

    list = hashtable_search(contigStringToContigIndex, contig);
#ifdef BEN_DEBUG
    assert(list != NULL);
    assert(list->length > 0);
#endif
    k = 0;
    j = INT_MAX;
    for (i = 0; i < list->length; i++) {
        iA = list->list[i];
        if (iA[0] <= start && start - iA[0] < j) {
            j = start - iA[0];
            k = iA[1];
        }
    }
#ifdef BEN_DEBUG
    assert(j != INT_MAX);
#endif
    return k;
}

void pinchMerge(struct PinchGraph *graph, struct PairwiseAlignment *pA,
        void(*addFunction)(struct PinchGraph *pinchGraph, struct Piece *,
                struct Piece *, stHash *, void *), void *extraParameter,
        stHash *vertexToAdjacencyComponents) {
    /*
     * Method to pinch together the graph using all the aligned matches in the
     * input alignment.
     */
    int32_t i, j, k;
    Name contig1, contig2;
    struct AlignmentOperation *op;
    static struct Piece piece1;
    static struct Piece rPiece1;
    static struct Piece piece2;
    static struct Piece rPiece2;

    //links the static pieces
    piece1.rPiece = &rPiece1;
    rPiece1.rPiece = &piece1;

    piece2.rPiece = &rPiece2;
    rPiece2.rPiece = &piece2;

    j = pA->start1;
    k = pA->start2;

    contig1 = cactusMisc_stringToName(pA->contig1);
    contig2 = cactusMisc_stringToName(pA->contig2);

    //logPairwiseAlignment(pA);

    for (i = 0; i < pA->operationList->length; i++) {
        op = pA->operationList->list[i];
        if (op->opType == PAIRWISE_MATCH) {
            if (op->length >= 1) { //deal with the possibility of a zero length match (strange, but not illegal)
                if (pA->strand1) {
                    piece_recycle(&piece1, contig1, j, j + op->length - 1);
                } else {
                    piece_recycle(&piece1, contig1, -(j - 1), -(j - op->length));
                }
                if (pA->strand2) {
                    piece_recycle(&piece2, contig2, k, k + op->length - 1);
                } else {
                    piece_recycle(&piece2, contig2, -(k - 1), -(k - op->length));
                }
                addFunction(graph, &piece1, &piece2,
                        vertexToAdjacencyComponents, extraParameter);
            }
        }
        if (op->opType != PAIRWISE_INDEL_Y) {
            j += pA->strand1 ? op->length : -op->length;
        }
        if (op->opType != PAIRWISE_INDEL_X) {
            k += pA->strand2 ? op->length : -op->length;
        }
    }

    assert(j == pA->end1);
    assert(k == pA->end2);
}
