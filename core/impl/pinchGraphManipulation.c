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
    assert(lengthBlackEdges(vertex1) == lengthBlackEdges(vertex2));
    assert(lengthGreyEdges(vertex1) == 1);
    assert(lengthGreyEdges(vertex2) == 1);
    assert(getFirstGreyEdge(vertex1) == vertex2);
    assert(getFirstGreyEdge(vertex2) == vertex1);

    //if(lengthBlackEdges(vertex1) > 1) {
    //	assert(FALSE);
    //}

    //For each black edge to vertex1 find consecutive edge from vertex2, then join.
    while (lengthBlackEdges(vertex1) > 0) {
        assert(lengthBlackEdges(vertex1) == lengthBlackEdges(vertex2));

        struct PinchEdge *edge1 = getFirstBlackEdge(vertex1);
        assert(edge1 != NULL);
        assert(!isAStub(edge1));
        edge1 = edge1->rEdge;
        assert(edge1->to == vertex1);
        //first find the grey edge to attach to the new vertex we're about to create

        struct PinchEdge *edge2 = getNextEdge(graph, edge1, flower);
        assert(edge2 != NULL);
        assert(!isAStub(edge2));
        assert(edge2->from == vertex2);

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
    assert(lengthBlackEdges(vertex1) == 0);
    assert(lengthBlackEdges(vertex2) == 0);
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

    hash = create_hashtable(0, hashtable_key, hashtable_equalKey, NULL, free);
    for (i = 0; i < pinchGraph->vertices->length; i++) {
        vertex = pinchGraph->vertices->list[i];
        if (lengthBlackEdges(vertex) >= 1
                && !isAStub(getFirstBlackEdge(vertex))) {
            if (lengthBlackEdges(vertex) > maxDegree || treeCoverage(vertex,
                    flower) < minimumTreeCoverage) { //has a high degree and is not a stub/cap
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
     * This adds a bunch of extra edges to the list which should be undone. It ignored stub edges
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
    stList *adjacencyComponents = getAdjacencyComponents2(pinchGraph, linkStubComponentsToTheSinkComponent_passThroughFn);

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
            stSortedSetIterator *it = stSortedSet_getIterator(adjacencyComponent);
            while((vertex = stSortedSet_getNext(it)) != NULL) {
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
            while((vertex = stSortedSet_getNext(it)) != NULL) {
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
                        if(attachEnds) {
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
                if(lengthGreyEdges(vertex) == 1) { //is attached to the origin node.
                    assert(getFirstGreyEdge(vertex) == pinchGraph->vertices->list[0]);
                    disconnectVertices(vertex, pinchGraph->vertices->list[0]);
                }
                else {
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
    if (treeCoverage <= -0.001 || treeCoverage >= 1.001) {
        st_uglyf("The tree coverage for this case is: %f, %f \n", treeCoverage,
                wholeTreeCoverage);
    }
    assert(treeCoverage >= -0.001);
    assert(treeCoverage <= 1.0001);
    return treeCoverage;
}
