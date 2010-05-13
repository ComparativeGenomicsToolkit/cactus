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

	piece = (struct Piece *)mallocLocal(sizeof(struct Piece));
	rPiece = (struct Piece *)mallocLocal(sizeof(struct Piece));

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
	logDebug("Contig : " INT_STRING ", start : " INT_STRING ", end : " INT_STRING "\n",
			piece->contig, piece->start, piece->end);
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
	if(piece1->contig > piece2->contig) {
		return 1;
	}
	if(piece1->contig < piece2->contig) {
		return -1;
	}
	//Check if overlap.
	if(piece1->start <= piece2->start) {
		if(piece1->end >= piece2->start) {
			return 0;
		}
	}
	else {
		if(piece2->end >= piece1->start) {
			return 0;
		}
	}
	//Back to business
	if(piece1->start < piece2->start) {
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

struct PinchVertex *constructPinchVertex(struct PinchGraph *graph, int32_t vertexID, bool isEnd, bool isDeadEnd) {
	struct PinchVertex *pinchVertex;

	pinchVertex = mallocLocal(sizeof(struct PinchVertex));

	pinchVertex->blackEdges = constructBlackEdges();
	pinchVertex->greyEdges = constructGreyEdges();
	if(vertexID < 0) {
		pinchVertex->vertexID = graph->vertices->length;
		listAppend(graph->vertices, pinchVertex);
	}
	else {
		pinchVertex->vertexID = vertexID;
	}
	assert((!isEnd && !isDeadEnd) || (isEnd && !isDeadEnd) || (!isEnd && isDeadEnd));
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

void removeVertexFromGraphAndDestruct(struct PinchGraph *graph, struct PinchVertex *vertex) {
	/*
	 * Destructs pinch vertex and removes its reference from the parent graph structure.
	 */
	struct PinchVertex *highestVertex;

	highestVertex = graph->vertices->list[graph->vertices->length-1];
	assert(highestVertex->vertexID == graph->vertices->length-1);

	highestVertex->vertexID = vertex->vertexID;
	graph->vertices->list[highestVertex->vertexID] = highestVertex;

	graph->vertices->length--;
	destructPinchVertex(vertex);
}

void mergeVerticesBlackEdges(struct PinchVertex *oldVertex, struct PinchVertex *newVertex) {
	/*
	 * Merges black vertices.
	 */
	struct PinchEdge *edge;
	void *blackEdgeIterator = getBlackEdgeIterator(oldVertex);
	while((edge = getNextBlackEdge(oldVertex, blackEdgeIterator)) != NULL) {
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

void mergeVerticesGreyEdges(struct PinchVertex *oldVertex, struct PinchVertex *newVertex) {
	/*
	* Merges grey vertices.
	*/
	struct PinchVertex *vertex;
	void *greyEdgeIterator=getGreyEdgeIterator(oldVertex);
	while((vertex = getNextGreyEdge(oldVertex, greyEdgeIterator)) != NULL) {
		if(vertex != oldVertex) { //don't deal with self edge here, as will disrupt the iterator
			//Replace references in connected vertex.
			if(containsGreyEdge(vertex, oldVertex)) {
				removeGreyEdge(vertex, oldVertex);
			}
			connectVertices(vertex, newVertex);
		}
	}
	destructGreyEdgeIterator(greyEdgeIterator);
	if(containsGreyEdge(oldVertex, oldVertex)) { //now deal with self edge in vertex2.
		connectVertices(newVertex, newVertex);
	}
	assert(!containsGreyEdge(newVertex, oldVertex));
}

struct PinchVertex *mergeVertices(struct PinchGraph *graph, struct PinchVertex *vertex1, struct PinchVertex *vertex2) {
	/*
	 * Method to merge two vertices together to create one vertex in the graph.
	 */
	struct PinchVertex *vertex3;
	if (vertex1 != vertex2) {
		if(lengthBlackEdges(vertex1) + lengthGreyEdges(vertex1) < lengthBlackEdges(vertex2) + lengthGreyEdges(vertex2)) {
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

int32_t containsGreyEdge(struct PinchVertex *vertex, struct PinchVertex *vertex2) {
	return avl_find(vertex->greyEdges, vertex2) != NULL;
}

void insertGreyEdge(struct PinchVertex *vertex, struct PinchVertex *vertex2) {
	if(!containsGreyEdge(vertex, vertex2)) {
		avl_insert(vertex->greyEdges, vertex2);
		assert(containsGreyEdge(vertex, vertex2));
	}
}

void removeGreyEdge(struct PinchVertex *vertex, struct PinchVertex *vertex2)  {
	assert(containsGreyEdge(vertex, vertex2));
	avl_delete(vertex->greyEdges, vertex2);
	assert(!containsGreyEdge(vertex, vertex2));
}

struct PinchVertex *popGreyEdge(struct PinchVertex *vertex)  {
	assert(lengthGreyEdges(vertex) > 0);
	struct PinchVertex *vertex2 = getFirstGreyEdge(vertex);
	removeGreyEdge(vertex, vertex2);
	return vertex2;
}

void *getGreyEdgeIterator(struct PinchVertex *vertex) {
	struct avl_traverser *iterator;
	iterator = mallocLocal(sizeof(struct avl_traverser));
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
	if(!containsBlackEdge(vertex, edge)) {
		avl_insert(vertex->blackEdges, edge);
		assert(containsBlackEdge(vertex, edge));
	}
}

void removeBlackEdge(struct PinchVertex *vertex, struct PinchEdge *edge)  {
	assert(containsBlackEdge(vertex, edge));
	avl_delete(vertex->blackEdges, edge);
	assert(!containsBlackEdge(vertex, edge));
}

struct PinchEdge *popBlackEdge(struct PinchVertex *vertex)  {
	assert(lengthBlackEdges(vertex) > 0);
	struct PinchEdge *edge = getFirstBlackEdge(vertex);
	removeBlackEdge(vertex, edge);
	return edge;
}

void *getBlackEdgeIterator(struct PinchVertex *vertex) {
	struct avl_traverser *iterator;
	iterator = mallocLocal(sizeof(struct avl_traverser));
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

	pinchEdge = mallocLocal(sizeof(struct PinchEdge));
	rPinchEdge = mallocLocal(sizeof(struct PinchEdge));

	pinchEdge->rEdge = rPinchEdge;
	rPinchEdge->rEdge = pinchEdge;

	pinchEdge->piece = piece;
	rPinchEdge->piece = piece->rPiece;

	return pinchEdge;
}

void connectPinchEdge(struct PinchEdge *edge, struct PinchVertex *from, struct PinchVertex *to) {
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

void removePinchEdgeFromGraphAndDestruct(struct PinchGraph *graph, struct PinchEdge *edge) {
	/*
	 * Destroys the edge and removes it from the graph.
	 */
	avl_delete(graph->edges, edge);
	avl_delete(graph->edges, edge->rEdge);
	if(containsBlackEdge(edge->from, edge)) {
		removeBlackEdge(edge->from, edge);
	}
	if(containsBlackEdge(edge->to, edge->rEdge)) {
		removeBlackEdge(edge->to, edge->rEdge);
	}
	destructPinchEdge(edge);
}

int32_t edgeComparator(struct PinchEdge *edge1, struct PinchEdge *edge2, void *o) {
	/*
	 * Compares the pieces by their pieces.
	 */
	assert(o == NULL);
	return pieceComparator(edge1->piece, edge2->piece);
}

struct PinchEdge *getContainingBlackEdge(struct PinchGraph *graph, Name contig, int32_t position) {
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
	return edge2;
}

struct PinchEdge *getNextEdge(struct PinchGraph *graph, struct PinchEdge *edge, Net *net) {
	/*
	 * Gets the next in the sequence from the given edge, used for traversing paths in the pinch graph.
	 */
	struct PinchEdge *edge2;
	struct PinchVertex *vertex2;
	Cap *cap;

	edge2 = getContainingBlackEdge(graph, edge->piece->contig, edge->piece->end+1);
	if(edge2 != NULL) {
		return edge2;
	}
	void *iterator = getGreyEdgeIterator(edge->to);
	while((vertex2 = getNextGreyEdge(edge->to, iterator)) != NULL) {
		if(vertex_isEnd(vertex2)) {
			void *iterator2 = getBlackEdgeIterator(vertex2);
			while((edge2 = getNextBlackEdge(vertex2, iterator2)) != NULL) {
				if(edge2->piece->start == edge->piece->end+1) {
					cap = net_getCap(net, edge2->piece->contig);
					assert(cap != NULL);
					if(sequence_getName(cap_getSequence(cap)) == edge->piece->contig) {
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
	exitOnFailure(0, "Failed to find correct edge\n");
	return NULL;
}

void splitEdge_P(struct PinchGraph *graph,
				 struct PinchEdge *edge, int32_t position,
				 struct PinchVertex *vertex1, struct PinchVertex *vertex2) {
	struct PinchEdge *edge1;
	struct PinchEdge *edge2;

#ifdef BEN_DEBUG
	assert(edge->piece->start < position);
#endif

	//Otherwise create new pinch point.
	//Split piece
	edge1 = constructPinchEdge(constructPiece(edge->piece->contig,
								edge->piece->start, position-1));
	connectPinchEdge(edge1, edge->from, vertex1);

	edge2 = constructPinchEdge(constructPiece(edge->piece->contig,
								position, edge->piece->end));
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
							  int32_t position, int32_t leftOrRight) {
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
	if(edge->piece->start == position && leftOrRight == LEFT) {
		return edge->from;
	}

	//If pinch point already at the end..
	if(edge->piece->end == position && leftOrRight == RIGHT) {
		return edge->to;
	}

	if(sE_list == NULL) {
		sE_list = constructEmptyList(0, NULL);
	}
	sE_list->length = 0;

	//For each of the aligned edges, do the split.
	blackEdgeIterator = getBlackEdgeIterator(edge->from);
	while((edge2 = getNextBlackEdge(edge->from, blackEdgeIterator)) != NULL) {
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
	for(i=0; i<sE_list->length; i++) {
		edge2 = sE_list->list[i];
#ifdef BEN_DEBUG
		assert(j + edge2->piece->start <= edge2->piece->end);
#endif
		splitEdge_P(graph, edge2, j + edge2->piece->start, vertex1, vertex2);
	}

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

	pinchGraph = mallocLocal(sizeof(struct PinchGraph));
	pinchGraph->edges = avl_create((int32_t (*)(const void *, const void *, void *a))edgeComparator, NULL, NULL);
	pinchGraph->vertices = constructEmptyList(0, (void (*)(void *))destructPinchVertex);
	constructPinchVertex(pinchGraph, -1, 0, 0);

	return pinchGraph;
}

void destructPinchGraph_1(struct PinchEdge *edge, void *o) {
	assert(o == NULL);
	if(edge->from->vertexID > edge->to->vertexID) {
		destructPiece(edge->piece);
	}
	free(edge);
}

void destructPinchGraph(struct PinchGraph *pinchGraph) {
	avl_destroy(pinchGraph->edges, (void (*) (void *avl_item, void *avl_param))destructPinchGraph_1);
	//destroylist(pinchGraph->edges);
	destructList(pinchGraph->vertices);
	free(pinchGraph);
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

	vertexChain = mallocLocal(sizeof(struct VertexChain));
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
						struct PinchGraph *graph,
					    struct Piece *piece) {
	struct PinchVertex *vertex;
	struct PinchEdge *edge;

	resetVertexChain(vertexChain);

	//do any adjustments off the bat
	splitEdge(graph, piece->contig, piece->start, LEFT);
	splitEdge(graph, piece->contig, piece->end, RIGHT);

	vertex = splitEdge(graph, piece->contig, piece->start, LEFT);
	intListAppend(vertexChain->coordinates, 0);
	intListAppend(vertexChain->leftsOrRights, LEFT);
	listAppend(vertexChain->listOfVertices, vertex);

	//follow chain to get remaining vertices
	edge = getContainingBlackEdge(graph, piece->contig, piece->start);
	while(edge->piece->end < piece->end) {
		intListAppend(vertexChain->coordinates, edge->piece->end - piece->start);
		intListAppend(vertexChain->leftsOrRights, RIGHT);
		listAppend(vertexChain->listOfVertices, edge->to);

		edge = getContainingBlackEdge(graph, piece->contig, edge->piece->end+1);
		intListAppend(vertexChain->coordinates, edge->piece->start - piece->start);
		intListAppend(vertexChain->leftsOrRights, LEFT);
		listAppend(vertexChain->listOfVertices, edge->from);
	}

	//add the second vertex.
	vertex = splitEdge(graph, piece->contig, piece->end, RIGHT);
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

	if(vertexChain1->listOfVertices->length != vertexChain2->listOfVertices->length) {
		return FALSE;
	}
	for(i=0; i<vertexChain1->coordinates->length; i++) {
		j = vertexChain1->coordinates->list[i];
		k = vertexChain2->coordinates->list[i];
		if(j != k) {
			return FALSE;
		}

		j = vertexChain1->leftsOrRights->list[i];
		k = vertexChain2->leftsOrRights->list[i];
#ifdef BEN_DEBUG
		assert(j == LEFT || j == RIGHT);
		assert(k == LEFT || k == RIGHT);
#endif
		if(j != k) {
			return FALSE;
		}
	}
	return TRUE;
}

void updateVertexAdjacencyComponentLabels(struct hashtable *vertexAdjacencyComponents,
		struct PinchVertex *vertex,
		struct PinchGraph *pinchGraph) {
	/*
	 * Method establishes which adjacency component the vertex belongs in.
	 */
	int32_t *i, j;
	i = hashtable_search(vertexAdjacencyComponents, vertex);
	if(i == NULL) {
		struct List *list = constructEmptyList(0, NULL);
		while(1) {
			listAppend(list, vertex);
			assert(lengthGreyEdges(vertex) == 1);
			struct PinchVertex *vertex2 = getFirstGreyEdge(vertex);
			i = hashtable_search(vertexAdjacencyComponents, vertex2);
			if(i != NULL) {
				break;
			}
			listAppend(list, vertex2);
			assert(lengthBlackEdges(vertex2) > 0);
			vertex = getFirstBlackEdge(vertex2)->to;
			i = hashtable_search(vertexAdjacencyComponents, vertex);
			if(i != NULL) {
				break;
			}
		}
		for(j=0; j<list->length; j++) {
			assert(hashtable_search(vertexAdjacencyComponents, list->list[j]) == NULL);
			hashtable_insert(vertexAdjacencyComponents, list->list[j], constructInt(*i));
		}
		destructList(list);
	}
}

void updateVertexAdjacencyComponentLabelsForChain(struct hashtable *vertexAdjacencyComponents,
		struct VertexChain *vertexChain,
		struct PinchGraph *pinchGraph) {
	/*
	 * Method runs through the vertices in the vertex chain and ensures each vertex has a label.
	 */
	int32_t i;
	for(i=0; i<vertexChain->listOfVertices->length; i++) {
		updateVertexAdjacencyComponentLabels(vertexAdjacencyComponents, vertexChain->listOfVertices->list[i], pinchGraph);
	}
}

void pinchMergePiece_getChainOfVertices(struct PinchGraph *graph,
										  struct Piece *piece1,
										  struct Piece *piece2,
										  struct VertexChain *vertexChain1,
										  struct VertexChain *vertexChain2,
										  struct hashtable *vertexAdjacencyComponents) {
	int32_t i, j, k;
	getChainOfVertices(vertexChain1, graph, piece1);
	getChainOfVertices(vertexChain2, graph, piece2);

	while(pinchMergePiece_P(vertexChain1, vertexChain2) == FALSE) {
		/*
		 * match up the set of vertices for each chain
		 */
		for(i=0; i<vertexChain1->coordinates->length; i++) {
			j = vertexChain1->coordinates->list[i];
			k = vertexChain1->leftsOrRights->list[i];
			//now search if there is an equivalent vertex.
			splitEdge(graph, piece2->contig, piece2->start + j, k);
		}
		for(i=0; i<vertexChain2->coordinates->length; i++) {
			j = vertexChain2->coordinates->list[i];
			k = vertexChain2->leftsOrRights->list[i];
			//now search if there is an equivalent vertex.
			splitEdge(graph, piece1->contig, piece1->start + j, k);
		}

		getChainOfVertices(vertexChain1, graph, piece1);
		getChainOfVertices(vertexChain2, graph, piece2);
	}

	/*
	 * Label the new vertices in the chain with adjacency component labels.
	 */
	updateVertexAdjacencyComponentLabelsForChain(vertexAdjacencyComponents, vertexChain1, graph);
	updateVertexAdjacencyComponentLabelsForChain(vertexAdjacencyComponents, vertexChain2, graph);
}

struct VertexChain *pMS_vertexChain1 = NULL;
struct VertexChain *pMS_vertexChain2 = NULL;

void pinchMergePiece(struct PinchGraph *graph,
					   struct Piece *piece1,
					   struct Piece *piece2,
					   struct hashtable *vertexAdjacencyComponents) {
	/*
	 * Pinches the graph (with the minimum number of required pinches, to
	 * represent the contiguous alignment of the two pieces.
	 *
	 * Pieces have to be of equal length.
	 */
	int32_t i, j, k, contig;
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
	splitEdge(graph, piece1->contig, piece1->start, LEFT);
	splitEdge(graph, piece1->contig, piece1->end, RIGHT);
	splitEdge(graph, piece2->contig, piece2->start, LEFT);
	splitEdge(graph, piece2->contig, piece2->end, RIGHT);

	pinchMergePiece_getChainOfVertices(graph, piece1, piece2, pMS_vertexChain1, pMS_vertexChain2, vertexAdjacencyComponents);

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
	for(i=0; i<pMS_vertexChain1->listOfVertices->length; i++) {
		assert(hashtable_search(vertexAdjacencyComponents, pMS_vertexChain1->listOfVertices->list[i]) != NULL);
		assert(hashtable_search(vertexAdjacencyComponents, pMS_vertexChain2->listOfVertices->list[i]) != NULL);
		j = *(int32_t *)hashtable_search(vertexAdjacencyComponents, pMS_vertexChain1->listOfVertices->list[i]);
		k = *(int32_t *)hashtable_search(vertexAdjacencyComponents, pMS_vertexChain2->listOfVertices->list[i]);
		if(j != k) {
			return;
		}
	}

	/*
	 * Merge the lists of vertices to do the final merge.
	 */
	for(i=0; i<pMS_vertexChain1->listOfVertices->length;) {
		vertex1 = pMS_vertexChain1->listOfVertices->list[i];
		vertex2 = pMS_vertexChain2->listOfVertices->list[i];

#ifdef BEN_DEBUG
		assert(lengthBlackEdges(vertex1) > 0);
#endif
		//check if the two vertices are the ends of same piece.
		edge = getFirstBlackEdge(vertex1);
		if(edge->to == vertex2) {
			//if edge piece is of length greater than one.
			if(edge->piece->end - edge->piece->start > 0) {
				j = (edge->piece->end - edge->piece->start + 1) / 2 + edge->piece->start - 1;
				k = 1 + ((edge->piece->end - edge->piece->start + 1) % 2);

				contig = edge->piece->contig;
				vertex4 = splitEdge(graph, contig, j, RIGHT);
				vertex5 = splitEdge(graph, contig, j+k, LEFT);

#ifdef BEN_DEBUG
				//debug checks
				assert(lengthGreyEdges(vertex4) == 1);
				assert(lengthGreyEdges(vertex5) == 1);
				if(k == 1) {
					assert(getFirstGreyEdge(vertex4) == vertex5);
					assert(getFirstGreyEdge(vertex5) == vertex4);
				}
#endif
				/*
				 * The new vertices are not in the chain, so we re parse the vertex chain and start again.
				 */
				pinchMergePiece_getChainOfVertices(graph, piece1, piece2,
						pMS_vertexChain1, pMS_vertexChain2, vertexAdjacencyComponents);

				i = 0;
				continue;
			}
			else {
				/*
				 * In this case we do nothing, as we can't have self black edges and move one.
				 */
				i++;
				continue;
			}
		}
		else {
			assert(hashtable_search(vertexAdjacencyComponents, vertex1) != NULL);
			assert(hashtable_search(vertexAdjacencyComponents, vertex2) != NULL);
			k = *(int32_t *)hashtable_search(vertexAdjacencyComponents, vertex1);
			assert(k == *(int32_t *)hashtable_search(vertexAdjacencyComponents, vertex2));
			destructInt(hashtable_remove(vertexAdjacencyComponents, vertex1, 0));
			destructInt(hashtable_remove(vertexAdjacencyComponents, vertex2, 0));
			vertex3 = mergeVertices(graph, vertex1, vertex2);
			hashtable_insert(vertexAdjacencyComponents, vertex3, constructInt(k));
		}

		for(j=i+1; j<pMS_vertexChain1->listOfVertices->length; j++) {
			if(pMS_vertexChain1->listOfVertices->list[j] == vertex1 || pMS_vertexChain1->listOfVertices->list[j] == vertex2) {
				pMS_vertexChain1->listOfVertices->list[j] = vertex3;
			}
			if(pMS_vertexChain2->listOfVertices->list[j] == vertex1 || pMS_vertexChain2->listOfVertices->list[j] == vertex2) {
				pMS_vertexChain2->listOfVertices->list[j] = vertex3;
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

int32_t pinchMerge_getContig(char *contig, int32_t start, struct hashtable *contigStringToContigIndex) {
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
	for(i=0; i<list->length; i++) {
		iA = list->list[i];
		if(iA[0] <= start && start - iA[0] < j) {
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
		void (*addFunction)(struct PinchGraph *pinchGraph, struct Piece *, struct Piece *, struct hashtable *, void *),
		void *extraParameter,
		struct hashtable *vertexAdjacencyComponents) {
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

	contig1 = netMisc_stringToName(pA->contig1);
	contig2 = netMisc_stringToName(pA->contig2);

	logPairwiseAlignment(pA);

	for(i=0; i<pA->operationList->length; i++) {
		op = pA->operationList->list[i];
		if(op->opType == PAIRWISE_MATCH) {
			assert(op->length >= 1);
			if(pA->strand1) {
				piece_recycle(&piece1, contig1, j, j+op->length-1);
			}
			else {
				piece_recycle(&piece1, contig1, -(j-1), -(j-op->length));
			}
			if(pA->strand2) {
				piece_recycle(&piece2, contig2, k, k+op->length-1);
			}
			else {
				piece_recycle(&piece2, contig2, -(k-1), -(k-op->length));
			}
			addFunction(graph, &piece1, &piece2, vertexAdjacencyComponents, extraParameter);
		}
		if(op->opType != PAIRWISE_INDEL_Y) {
			j += pA->strand1 ? op->length : -op->length;
		}
		if(op->opType != PAIRWISE_INDEL_X) {
			k += pA->strand2 ? op->length : -op->length;
		}
	}

	assert(j == pA->end1);
	assert(k == pA->end2);
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

	if(hashtable_search(hash, thing) != NULL) {
		j = *((int32_t *)hashtable_search(hash, thing));
		j = j % 13;
	}
	else {
		j = 100;
	}

	switch(j) {
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

void writeOutPinchGraphWithChains(struct PinchGraph *pinchGraph, struct hashtable *edgeColours,
								  struct List *groups, FILE *fileHandle) {
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

	logDebug("Writing the pinch graph\n");

	hash2 = create_hashtable(pinchGraph->vertices->length*10,
							 hashtable_key, hashtable_equalKey,
							 NULL, (void (*)(void *))destructInt);

	if(groups != NULL) {
		for(i=0; i<groups->length; i++) {
			group = groups->list[i];
			for(j=0; j<group->length; j++) {
				vertex = group->list[j];
				hashtable_insert(hash2, vertex, constructInt(i));
			}
		}
	}

	logDebug("Done the preliminaries\n");

	//Write the preliminaries.
	fprintf(fileHandle, "graph G {\n");
	fprintf(fileHandle, "overlap=false\n");
	//Write the vertices.
	for(i=0; i<pinchGraph->vertices->length;i++) {
		vertex = pinchGraph->vertices->list[i];
#ifdef BEN_DEBUG
		assert(vertex->vertexID == i);
#endif
		colour = getColour(hash2, vertex);
		fprintf(fileHandle, "node[width=0.3,height=0.3,shape=circle,colour=%s,fontsize=14];\n", colour);
		fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n", vertex->vertexID, vertex->vertexID);
	}

	logDebug("Written the vertices\n");

	for(i=0; i<pinchGraph->vertices->length;i++) {
		vertex = pinchGraph->vertices->list[i];
		if(lengthBlackEdges(vertex) > 0) {
			edge = getFirstBlackEdge(vertex);
#ifdef BEN_DEBUG
			assert(vertex->vertexID != edge->to->vertexID);
#endif
			if(vertex->vertexID < edge->to->vertexID) {
				colour = getColour(edgeColours, edge->piece);
				fprintf(fileHandle, "edge[color=%s,len=2.5,weight=100,dir=forward];\n", colour);

				void *blackEdgeIterator = getBlackEdgeIterator(vertex);
				while((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
					fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n [label=\"" INT_STRING ":" INT_STRING ":%s\"];\n",
							edge->from->vertexID, edge->to->vertexID, edge->piece->start, edge->piece->end,
							netMisc_nameToStringStatic(edge->piece->contig));
				}
				destructBlackEdgeIterator(blackEdgeIterator);
			}
		}
	}

	logDebug("Written the black edges\n");

	//Write the grey edges.
	fprintf(fileHandle, "edge[color=grey20,len=2.5,weight=100,dir=none];\n");
	for(i=0; i<pinchGraph->vertices->length;i++) {
		vertex = pinchGraph->vertices->list[i];
		k = 0;
		void *greyEdgeIterator=getGreyEdgeIterator(vertex);
		while((vertex2 = getNextGreyEdge(vertex, greyEdgeIterator)) != NULL) {
			if(vertex->vertexID < vertex2->vertexID) {
				fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n", vertex->vertexID, vertex2->vertexID);
			}
			if(vertex->vertexID == vertex2->vertexID) {
				//Get one for every edge.
				if(k == 0) {
					fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n", vertex->vertexID, vertex2->vertexID);
					k = 1;
				}
				else {
					k = 0;
				}
			}
		}
		destructGreyEdgeIterator(greyEdgeIterator);
	}
	fprintf(fileHandle, "}\n");

	logDebug("Written the grey edges\n");

	hashtable_destroy(hash2, TRUE, FALSE);

	logDebug("Written the pinch graph\n");
}
