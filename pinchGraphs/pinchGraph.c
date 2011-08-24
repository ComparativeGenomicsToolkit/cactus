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
#include "cactus.h"
#include "pairwiseAlignment.h"
//#include "cactusGraph.h"

/*
 * Basic stuff for building and manipulating pinch graphs
 */

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Pieces
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void piece_recycle(struct Piece *piece, Name contig, int32_t start, int32_t end) {
    piece->contig = contig;

    piece->start = start;
    piece->end = end;

    piece->rPiece->contig = piece->contig;
    piece->rPiece->start = -end;
    piece->rPiece->end = -start;
}

struct Piece *constructPiece(Name contig, int32_t start, int32_t end) {
    struct Piece *piece;
    struct Piece *rPiece;

    piece = (struct Piece *) st_malloc(sizeof(struct Piece));
    rPiece = (struct Piece *) st_malloc(sizeof(struct Piece));

    piece->rPiece = rPiece;
    rPiece->rPiece = piece;

    piece_recycle(piece, contig, start, end);
    return piece;
}

void destructPiece(struct Piece *piece) {
    free(piece->rPiece);
    free(piece);
}

void logPiece(struct Piece *piece) {
    st_logDebug("Contig : %s, start : " INT_STRING ", end : " INT_STRING "\n",
            cactusMisc_nameToStringStatic(piece->contig), piece->start,
            piece->end);
}

int pieceComparatorPointers(struct Piece **piece1, struct Piece **piece2) {
    return pieceComparator(*piece1, *piece2);
}

int pieceComparator(struct Piece *piece1, struct Piece *piece2) {
    /*
     * Compares two pieces to allow an ordering on any two pieces
     * to be constructued.
     */
    //Compare the contigs
    int32_t i = cactusMisc_nameCompare(piece1->contig, piece2->contig);
    if (i != 0) {
        return i;
    }
    //Check if overlap.
    if (piece1->start <= piece2->start) {
        if (piece1->end >= piece2->start) {
            return 0;
        }
    } else {
        if (piece2->end >= piece1->start) {
            return 0;
        }
    }
    //Back to business
    if (piece1->start < piece2->start) {
        return -1;
    }
    return 1;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Vertex methods
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchVertex *constructPinchVertex(struct PinchGraph *graph,
        int32_t vertexID, bool isEnd, bool isDeadEnd) {
    struct PinchVertex *pinchVertex;

    pinchVertex = st_malloc(sizeof(struct PinchVertex));

    pinchVertex->blackEdges = constructBlackEdges();
    pinchVertex->greyEdges = constructGreyEdges();
    if (vertexID < 0) {
        pinchVertex->vertexID = graph->vertices->length;
        listAppend(graph->vertices, pinchVertex);
    } else {
        pinchVertex->vertexID = vertexID;
    }
    assert((!isEnd && !isDeadEnd) || (isEnd && !isDeadEnd) || (!isEnd
            && isDeadEnd));
    pinchVertex->isEnd = isEnd;
    pinchVertex->isDeadEnd = isDeadEnd;

    return pinchVertex;
}

bool vertex_isEnd(struct PinchVertex *vertex) {
    return vertex->isEnd;
}

bool vertex_isDeadEnd(struct PinchVertex *vertex) {
    return vertex->isDeadEnd;
}

void destructPinchVertex(struct PinchVertex *pinchVertex) {
    //Does not destroy the connected edges, this is left up to the containing pinch-graph
    destructBlackEdges(pinchVertex->blackEdges);
    destructGreyEdges(pinchVertex->greyEdges);
    free(pinchVertex);
}

void removeVertexFromGraphAndDestruct(struct PinchGraph *graph,
        struct PinchVertex *vertex) {
    /*
     * Destructs pinch vertex and removes its reference from the parent graph structure.
     */
    struct PinchVertex *highestVertex;

    highestVertex = graph->vertices->list[graph->vertices->length - 1];
    assert(highestVertex->vertexID == graph->vertices->length - 1);

    highestVertex->vertexID = vertex->vertexID;
    graph->vertices->list[highestVertex->vertexID] = highestVertex;

    graph->vertices->length--;
    destructPinchVertex(vertex);
}

void mergeVerticesBlackEdges(struct PinchVertex *oldVertex,
        struct PinchVertex *newVertex) {
    /*
     * Merges black vertices.
     */
    struct PinchEdge *edge;
    void *blackEdgeIterator = getBlackEdgeIterator(oldVertex);
    while ((edge = getNextBlackEdge(oldVertex, blackEdgeIterator)) != NULL) {
#ifdef BEN_DEBUG
        assert(edge->from == oldVertex);
        assert(edge->to != newVertex);
#endif
        edge->from = newVertex;
        edge->rEdge->to = newVertex;
        insertBlackEdge(newVertex, edge);
    }
    destructBlackEdgeIterator(blackEdgeIterator);
}

void mergeVerticesGreyEdges(struct PinchVertex *oldVertex,
        struct PinchVertex *newVertex) {
    /*
     * Merges grey vertices.
     */
    struct PinchVertex *vertex;
    void *greyEdgeIterator = getGreyEdgeIterator(oldVertex);
    while ((vertex = getNextGreyEdge(oldVertex, greyEdgeIterator)) != NULL) {
        if (vertex != oldVertex) { //don't deal with self edge here, as will disrupt the iterator
            //Replace references in connected vertex.
            if (containsGreyEdge(vertex, oldVertex)) {
                removeGreyEdge(vertex, oldVertex);
            }
            connectVertices(vertex, newVertex);
        }
    }
    destructGreyEdgeIterator(greyEdgeIterator);
    if (containsGreyEdge(oldVertex, oldVertex)) { //now deal with self edge in vertex2.
        connectVertices(newVertex, newVertex);
    }
    assert(!containsGreyEdge(newVertex, oldVertex));
}

struct PinchVertex *mergeVertices(struct PinchGraph *graph,
        struct PinchVertex *vertex1, struct PinchVertex *vertex2) {
    /*
     * Method to merge two vertices together to create one vertex in the graph.
     */
    struct PinchVertex *vertex3;
    if (vertex1 != vertex2) {
        if (lengthBlackEdges(vertex1) + lengthGreyEdges(vertex1)
                < lengthBlackEdges(vertex2) + lengthGreyEdges(vertex2)) {
            //This switch makes a huge difference to the performance of the algorithm.
            //Because it is unlikely that two vertices with very high degree
            //which are homologous
            //will exist in the graph for long without being merged. Hence, we always choose to merge
            //the smaller node into the larger.
            vertex3 = vertex1;
            vertex1 = vertex2;
            vertex2 = vertex3;
        }
        //we will delete the highest vertex, exchange it for vertex2
        assert(graph->vertices->length > 0);
        vertex3 = graph->vertices->list[--graph->vertices->length];
        vertex3->vertexID = vertex2->vertexID;
        graph->vertices->list[vertex3->vertexID] = vertex3;

        //black edges
        mergeVerticesBlackEdges(vertex2, vertex1);
        //grey edges
        mergeVerticesGreyEdges(vertex2, vertex1);
        destructPinchVertex(vertex2);
    }
    return vertex1;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for accessing adjacency (grey) edges
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static int32_t oComparator(const void *o1, const void *o2, void *a) {
    /*
     * Compares the objects by there address.
     */
    assert(a == NULL);
    return o1 > o2 ? 1 : o1 < o2 ? -1 : 0;
}

void *constructGreyEdges() {
    return avl_create(oComparator, NULL, NULL);
}

void destructGreyEdges(void *greyEdges) {
    avl_destroy(greyEdges, NULL);
}

int32_t lengthGreyEdges(struct PinchVertex *vertex) {
    return avl_count(vertex->greyEdges);
}

struct PinchVertex *getFirstGreyEdge(struct PinchVertex *vertex) {
    assert(lengthGreyEdges(vertex) > 0);
    static struct avl_traverser iterator;
    avl_t_init(&iterator, vertex->greyEdges);
    return avl_t_first(&iterator, vertex->greyEdges);
}

int32_t containsGreyEdge(struct PinchVertex *vertex,
        struct PinchVertex *vertex2) {
    return avl_find(vertex->greyEdges, vertex2) != NULL;
}

void insertGreyEdge(struct PinchVertex *vertex, struct PinchVertex *vertex2) {
    if (!containsGreyEdge(vertex, vertex2)) {
        avl_insert(vertex->greyEdges, vertex2);
        assert(containsGreyEdge(vertex, vertex2));
    }
}

void removeGreyEdge(struct PinchVertex *vertex, struct PinchVertex *vertex2) {
    assert(containsGreyEdge(vertex, vertex2));
    avl_delete(vertex->greyEdges, vertex2);
    assert(!containsGreyEdge(vertex, vertex2));
}

struct PinchVertex *popGreyEdge(struct PinchVertex *vertex) {
    assert(lengthGreyEdges(vertex) > 0);
    struct PinchVertex *vertex2 = getFirstGreyEdge(vertex);
    removeGreyEdge(vertex, vertex2);
    return vertex2;
}

void *getGreyEdgeIterator(struct PinchVertex *vertex) {
    struct avl_traverser *iterator;
    iterator = st_malloc(sizeof(struct avl_traverser));
    avl_t_init(iterator, vertex->greyEdges);
    return iterator;
}

struct PinchVertex *getNextGreyEdge(struct PinchVertex *vertex, void *iterator) {
    assert(vertex != NULL);
    return avl_t_next(iterator);
}

void destructGreyEdgeIterator(void *iterator) {
    free(iterator);
}

void connectVertices(struct PinchVertex *vertex, struct PinchVertex *vertex2) {
    insertGreyEdge(vertex, vertex2);
    insertGreyEdge(vertex2, vertex);
}

void disconnectVertices(struct PinchVertex *vertex, struct PinchVertex *vertex2) {
    removeGreyEdge(vertex, vertex2);
    removeGreyEdge(vertex2, vertex);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for accessing pinch (black) edges
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void *constructBlackEdges() {
    return avl_create(oComparator, NULL, NULL);
    //return avl_create((int32_t (*)(const void *o1, const void *o2, void *a))edgeComparator, NULL, NULL);
}

void destructBlackEdges(void *blackEdges) {
    avl_destroy(blackEdges, NULL);
}

int32_t lengthBlackEdges(struct PinchVertex *vertex) {
    return avl_count(vertex->blackEdges);
}

struct PinchEdge *getFirstBlackEdge(struct PinchVertex *vertex) {
    assert(lengthBlackEdges(vertex) > 0);
    static struct avl_traverser iterator;
    avl_t_init(&iterator, vertex->blackEdges);
    return avl_t_first(&iterator, vertex->blackEdges);
}

int32_t containsBlackEdge(struct PinchVertex *vertex, struct PinchEdge *edge) {
    return avl_find(vertex->blackEdges, edge) != NULL;
}

void insertBlackEdge(struct PinchVertex *vertex, struct PinchEdge *edge) {
    if (!containsBlackEdge(vertex, edge)) {
        avl_insert(vertex->blackEdges, edge);
        assert(containsBlackEdge(vertex, edge));
    }
}

void removeBlackEdge(struct PinchVertex *vertex, struct PinchEdge *edge) {
    assert(containsBlackEdge(vertex, edge));
    avl_delete(vertex->blackEdges, edge);
    assert(!containsBlackEdge(vertex, edge));
}

struct PinchEdge *popBlackEdge(struct PinchVertex *vertex) {
    assert(lengthBlackEdges(vertex) > 0);
    struct PinchEdge *edge = getFirstBlackEdge(vertex);
    removeBlackEdge(vertex, edge);
    return edge;
}

void *getBlackEdgeIterator(struct PinchVertex *vertex) {
    struct avl_traverser *iterator;
    iterator = st_malloc(sizeof(struct avl_traverser));
    avl_t_init(iterator, vertex->blackEdges);
    return iterator;
}

struct PinchEdge *getNextBlackEdge(struct PinchVertex *vertex, void *iterator) {
    assert(vertex != NULL);
    return avl_t_next(iterator);
}

void destructBlackEdgeIterator(void *iterator) {
    free(iterator);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Edge methods
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


struct PinchEdge *constructPinchEdge(struct Piece *piece) {
    struct PinchEdge *pinchEdge;
    struct PinchEdge *rPinchEdge;

    pinchEdge = st_malloc(sizeof(struct PinchEdge));
    rPinchEdge = st_malloc(sizeof(struct PinchEdge));

    pinchEdge->rEdge = rPinchEdge;
    rPinchEdge->rEdge = pinchEdge;

    pinchEdge->piece = piece;
    rPinchEdge->piece = piece->rPiece;

    return pinchEdge;
}

void connectPinchEdge(struct PinchEdge *edge, struct PinchVertex *from,
        struct PinchVertex *to) {
    /*
     * Connected the pinch edge to its to and from vertices.
     */
    edge->from = from;
    edge->rEdge->to = from;
    edge->to = to;
    edge->rEdge->from = to;
    insertBlackEdge(edge->from, edge);
    insertBlackEdge(to, edge->rEdge);
}

void addPinchEdgeToGraph(struct PinchGraph *graph, struct PinchEdge *edge) {
    /*
     * Adds the edge to the pinch graph edge tree.
     */
    avl_insert(graph->edges, edge);
    avl_insert(graph->edges, edge->rEdge);
}

void destructPinchEdge(struct PinchEdge *pinchEdge) {
    //Destroys the contained piece but not the rEdge, so must be called on both.
    free(pinchEdge->rEdge);
    destructPiece(pinchEdge->piece);
    free(pinchEdge);
}

void removePinchEdgeFromGraphAndDestruct(struct PinchGraph *graph,
        struct PinchEdge *edge) {
    /*
     * Destroys the edge and removes it from the graph.
     */
    avl_delete(graph->edges, edge);
    avl_delete(graph->edges, edge->rEdge);
    if (containsBlackEdge(edge->from, edge)) {
        removeBlackEdge(edge->from, edge);
    }
    if (containsBlackEdge(edge->to, edge->rEdge)) {
        removeBlackEdge(edge->to, edge->rEdge);
    }
    destructPinchEdge(edge);
}

int32_t edgeComparator(struct PinchEdge *edge1, struct PinchEdge *edge2,
        void *o) {
    /*
     * Compares the pieces by their pieces.
     */
    assert(o == NULL);
    return pieceComparator(edge1->piece, edge2->piece);
}

struct PinchEdge *getContainingBlackEdge(struct PinchGraph *graph, Name contig,
        int32_t position) {
    /*
     * Gets edge containing the considered piece.
     */
    static struct PinchEdge edge;
    struct PinchEdge *edge2;
    static struct Piece piece;
    //Edge is a wrapper for the piece - this is all a bit silly.
    edge.piece = &piece;
    piece.contig = contig;
    piece.start = position;
    piece.end = position;
    //Now get the edge.
    edge2 = avl_find(graph->edges, &edge);
    //assert(edge2 != NULL);
    //if(edge2 == NULL) {
        //st_uglyf("The edge is null\n");
        //st_uglyf("The number of elements in the set %i\n", avl_count(graph->edges));
        //st_uglyf("The name %lli the position %i\n", contig, position);
    //}
    return edge2;
}

struct PinchEdge *getNextEdge(struct PinchGraph *graph, struct PinchEdge *edge,
        Flower *flower) {
    /*
     * Gets the next in the sequence from the given edge, used for traversing paths in the pinch graph.
     */
    struct PinchEdge *edge2;
    struct PinchVertex *vertex2;
    Cap *cap;

    Name contig = edge->piece->contig;

    assert(!vertex_isDeadEnd(edge->to));
    if (vertex_isDeadEnd(edge->from)) { //We're starting from a dead end
        assert(vertex_isEnd(edge->to));
        cap = flower_getCap(flower, edge->piece->contig);
        assert(cap != NULL);
        contig = sequence_getName(cap_getSequence(cap));
    }

    edge2 = getContainingBlackEdge(graph, contig, edge->piece->end
            + 1);
    if (edge2 != NULL) {
        return edge2;
    }
    void *iterator = getGreyEdgeIterator(edge->to);
    while ((vertex2 = getNextGreyEdge(edge->to, iterator)) != NULL) {
        if (vertex_isEnd(vertex2)) {
            void *iterator2 = getBlackEdgeIterator(vertex2);
            while ((edge2 = getNextBlackEdge(vertex2, iterator2)) != NULL) {
                if (edge2->piece->start == edge->piece->end + 1) {
                    cap = flower_getCap(flower, edge2->piece->contig);
                    assert(cap != NULL);
                    if (sequence_getName(cap_getSequence(cap))
                            == contig) {
                        destructGreyEdgeIterator(iterator);
                        destructBlackEdgeIterator(iterator2);
                        return edge2;
                    }
                }
            }
            destructBlackEdgeIterator(iterator2);
        }
    }
    destructGreyEdgeIterator(iterator);
    assert(0);
    return NULL;
}

void splitEdge_P(struct PinchGraph *graph, struct PinchEdge *edge,
        int32_t position, struct PinchVertex *vertex1,
        struct PinchVertex *vertex2) {
    struct PinchEdge *edge1;
    struct PinchEdge *edge2;

#ifdef BEN_DEBUG
    assert(edge->piece->start < position);
#endif

    //Otherwise create new pinch point.
    //Split piece
    edge1 = constructPinchEdge(constructPiece(edge->piece->contig,
            edge->piece->start, position - 1));
    connectPinchEdge(edge1, edge->from, vertex1);

    edge2 = constructPinchEdge(constructPiece(edge->piece->contig, position,
            edge->piece->end));
    connectPinchEdge(edge2, vertex2, edge->to);

#ifdef BEN_DEBUG
    assert(containsBlackEdge(edge->from, edge));
    assert(containsBlackEdge(edge->to, edge->rEdge));
#endif
    removePinchEdgeFromGraphAndDestruct(graph, edge);

    //add to graph after deleting old edge, or will have overlap.
    addPinchEdgeToGraph(graph, edge1);
    addPinchEdgeToGraph(graph, edge2);
}

struct List *sE_list = NULL;

struct PinchVertex *splitEdge(struct PinchGraph *graph, Name contig,
        int32_t position, int32_t leftOrRight, stHash *vertexToAdjacencyComponent) {
    /*
     * This function splits the edges (including all those in the same aligned edge).
     *
     * The position defines the adjacent position to the vertex being created/or sought.
     *
     * The left or right specifies the position of the vertex.
     *
     * For example, in the chain of positions X;Y;Z, the request splitEdge(g, Y, LEFT) would
     * retrieve/create the vertex between X and Y nearest to Y (on Y's left side). Conversely splitEdge(g, Y, RIGHT)
     * would get/create the vertex between Y and Z nearest Z (on Z's right side).
     */
    int32_t i, j;
    struct PinchVertex *vertex1;
    struct PinchVertex *vertex2;
    void *blackEdgeIterator;
    struct PinchEdge *edge;
    struct PinchEdge *edge2;

    edge = getContainingBlackEdge(graph, contig, position);

#ifdef BEN_DEBUG
    assert(edge != NULL);
    assert(leftOrRight == LEFT || leftOrRight == RIGHT);
#endif

    //If pinch point already at the start..
    if (edge->piece->start == position && leftOrRight == LEFT) {
        return edge->from;
    }

    //If pinch point already at the end..
    if (edge->piece->end == position && leftOrRight == RIGHT) {
        return edge->to;
    }

    if (sE_list == NULL) {
        sE_list = constructEmptyList(0, NULL);
    }
    sE_list->length = 0;

    void *adjacencyComponent1 = stHash_search(vertexToAdjacencyComponent, edge->from);
    void *adjacencyComponent2 = stHash_search(vertexToAdjacencyComponent, edge->to);
    void *adjacencyComponent3;
    if(adjacencyComponent1 == adjacencyComponent2) {
        adjacencyComponent3 = adjacencyComponent1;
    }
    else {
        adjacencyComponent3 = NULL;
    }

    //For each of the aligned edges, do the split.
    blackEdgeIterator = getBlackEdgeIterator(edge->from);
    while ((edge2 = getNextBlackEdge(edge->from, blackEdgeIterator)) != NULL) {
        listAppend(sE_list, edge2);
    }
    destructBlackEdgeIterator(blackEdgeIterator);

#ifdef BEN_ULTRA_DEBUG
    //check every edge has the same absolute length.
    blackEdgeIterator = getBlackEdgeIterator(edge->from);
    while((edge2 = getNextBlackEdge(edge->from, blackEdgeIterator)) != NULL) {
        assert(edge2->piece->end - edge2->piece->start == edge->piece->end - edge->piece->start);
    }
    destructBlackEdgeIterator(blackEdgeIterator);
#endif

    //Construct the vertices these things will be attached to.
    vertex1 = constructPinchVertex(graph, -1, 0, 0);
    vertex2 = constructPinchVertex(graph, -1, 0, 0);
    //Add grey edges to connect new vertices.
    connectVertices(vertex1, vertex2);

    j = position - edge->piece->start + (leftOrRight == LEFT ? 0 : 1);
    assert(j > 0);
    for (i = 0; i < sE_list->length; i++) {
        edge2 = sE_list->list[i];
#ifdef BEN_DEBUG
        assert(j + edge2->piece->start <= edge2->piece->end);
#endif
        splitEdge_P(graph, edge2, j + edge2->piece->start, vertex1, vertex2);
    }

    //Add to adjacency component hash
    stHash_insert(vertexToAdjacencyComponent, vertex1, adjacencyComponent3);
    stHash_insert(vertexToAdjacencyComponent, vertex2, adjacencyComponent3);

    //Return either left or right new vertex, dependent on which was made.
    return leftOrRight == LEFT ? vertex2 : vertex1;
}

int32_t isAStub(struct PinchEdge *edge) {
    return vertex_isEnd(edge->from) || vertex_isEnd(edge->to);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//The actual graph construction methods (the constructor has been moved
//to the reconstructionTree file, as it builds from an instance of
//a reconstruction problem.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchGraph *pinchGraph_construct() {
    struct PinchGraph *pinchGraph;

    pinchGraph = st_malloc(sizeof(struct PinchGraph));
    pinchGraph->edges = avl_create((int32_t(*)(const void *, const void *,
            void *a)) edgeComparator, NULL, NULL);
    pinchGraph->vertices = constructEmptyList(0,
            (void(*)(void *)) destructPinchVertex);
    constructPinchVertex(pinchGraph, -1, 0, 0);

    return pinchGraph;
}

void destructPinchGraph_1(struct PinchEdge *edge, void *o) {
    assert(o == NULL);
    if (edge->from->vertexID > edge->to->vertexID) {
        destructPiece(edge->piece);
    }
    free(edge);
}

void destructPinchGraph(struct PinchGraph *pinchGraph) {
    avl_destroy(pinchGraph->edges,
            (void(*)(void *avl_item, void *avl_param)) destructPinchGraph_1);
    //destroylist(pinchGraph->edges);
    destructList(pinchGraph->vertices);
    free(pinchGraph);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for writing out the basic pinch graph.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

char *getColour(struct hashtable *hash, void *thing) {
    static char cA[50];
    char *colour;
    int32_t j;

    if (hashtable_search(hash, thing) != NULL) {
        j = *((int32_t *) hashtable_search(hash, thing));
        j = j % 13;
    } else {
        j = 100;
    }

    switch (j) {
        case 0:
            colour = "red";
            break;
        case 1:
            colour = "blue";
            break;
        case 2:
            colour = "green";
            break;
        case 3:
            colour = "yellow";
            break;
        case 4:
            colour = "cyan";
            break;
        case 5:
            colour = "magenta";
            break;
        case 6:
            colour = "orange";
            break;
        case 7:
            colour = "purple";
            break;
        case 8:
            colour = "brown";
            break;
        case 9:
            colour = "palegreen";
            break;
        case 10:
            colour = "lightpink";
            break;
        case 11:
            colour = "lightcyan";
            break;
        case 12:
            colour = "black";
            break;
        default:
            colour = "grey80";
    }
    strcpy(cA, colour);
    return cA;
}

void writeOutPinchGraphWithChains(struct PinchGraph *pinchGraph,
        struct hashtable *edgeColours, struct List *groups, FILE *fileHandle) {
    /*
     * Writes out a graph in 'dot' format, compatible with graphviz.
     *
     * The associated function 'writeOutEdgePieces' gives a way of decoding
     * the pieces associated with the black (piece containing) edges
     * of the graph.
     */
    int32_t i, j, k;
    struct PinchVertex *vertex;
    struct PinchVertex *vertex2;
    struct PinchEdge *edge;
    struct hashtable *hash2;
    struct List *group;
    char *colour;

    st_logDebug("Writing the pinch graph\n");

    hash2 = create_hashtable(pinchGraph->vertices->length * 10, hashtable_key,
            hashtable_equalKey, NULL, (void(*)(void *)) destructInt);

    if (groups != NULL) {
        for (i = 0; i < groups->length; i++) {
            group = groups->list[i];
            for (j = 0; j < group->length; j++) {
                vertex = group->list[j];
                hashtable_insert(hash2, vertex, constructInt(i));
            }
        }
    }

    st_logDebug("Done the preliminaries\n");

    //Write the preliminaries.
    fprintf(fileHandle, "graph G {\n");
    fprintf(fileHandle, "overlap=false\n");
    //Write the vertices.
    for (i = 0; i < pinchGraph->vertices->length; i++) {
        vertex = pinchGraph->vertices->list[i];
#ifdef BEN_DEBUG
        assert(vertex->vertexID == i);
#endif
        colour = getColour(hash2, vertex);
        fprintf(
                fileHandle,
                "node[width=0.3,height=0.3,shape=circle,colour=%s,fontsize=14];\n",
                colour);
        fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n",
                vertex->vertexID, vertex->vertexID);
    }

    st_logDebug("Written the vertices\n");

    for (i = 0; i < pinchGraph->vertices->length; i++) {
        vertex = pinchGraph->vertices->list[i];
        if (lengthBlackEdges(vertex) > 0) {
            edge = getFirstBlackEdge(vertex);
#ifdef BEN_DEBUG
            assert(vertex->vertexID != edge->to->vertexID);
#endif
            if (vertex->vertexID < edge->to->vertexID) {
                colour = getColour(edgeColours, edge->piece);
                fprintf(fileHandle,
                        "edge[color=%s,len=2.5,weight=100,dir=forward];\n",
                        colour);

                void *blackEdgeIterator = getBlackEdgeIterator(vertex);
                while ((edge = getNextBlackEdge(vertex, blackEdgeIterator))
                        != NULL) {
                    fprintf(
                            fileHandle,
                            "n" INT_STRING "n -- n" INT_STRING "n [label=\"" INT_STRING ":" INT_STRING ":%s\"];\n",
                            edge->from->vertexID, edge->to->vertexID,
                            edge->piece->start, edge->piece->end,
                            cactusMisc_nameToStringStatic(edge->piece->contig));
                }
                destructBlackEdgeIterator(blackEdgeIterator);
            }
        }
    }

    st_logDebug("Written the black edges\n");

    //Write the grey edges.
    fprintf(fileHandle, "edge[color=grey20,len=2.5,weight=100,dir=none];\n");
    for (i = 0; i < pinchGraph->vertices->length; i++) {
        vertex = pinchGraph->vertices->list[i];
        k = 0;
        void *greyEdgeIterator = getGreyEdgeIterator(vertex);
        while ((vertex2 = getNextGreyEdge(vertex, greyEdgeIterator)) != NULL) {
            if (vertex->vertexID < vertex2->vertexID) {
                fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n",
                        vertex->vertexID, vertex2->vertexID);
            }
            if (vertex->vertexID == vertex2->vertexID) {
                //Get one for every edge.
                if (k == 0) {
                    fprintf(fileHandle,
                            "n" INT_STRING "n -- n" INT_STRING "n;\n",
                            vertex->vertexID, vertex2->vertexID);
                    k = 1;
                } else {
                    k = 0;
                }
            }
        }
        destructGreyEdgeIterator(greyEdgeIterator);
    }
    fprintf(fileHandle, "}\n");

    st_logDebug("Written the grey edges\n");

    hashtable_destroy(hash2, TRUE, FALSE);

    st_logDebug("Written the pinch graph\n");
}
