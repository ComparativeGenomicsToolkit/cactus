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
#include "net.h"
#include "pairwiseAlignment.h"
#include "cactusGraph.h"

/*
 * Basic stuff for building and manipulating pinch graphs
 */

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Segments
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void segment_recycle(struct Segment *segment, int32_t contig, int32_t start, int32_t end) {
	segment->contig = contig;

	segment->start = start;
	segment->end = end;

	segment->rSegment->contig = segment->contig;
	segment->rSegment->start = -end;
	segment->rSegment->end = -start;
}

struct Segment *constructSegment(int32_t contig, int32_t start, int32_t end) {
	struct Segment *segment;
	struct Segment *rSegment;

	segment = (struct Segment *)mallocLocal(sizeof(struct Segment));
	rSegment = (struct Segment *)mallocLocal(sizeof(struct Segment));

	segment->rSegment = rSegment;
	rSegment->rSegment = segment;

	segment_recycle(segment, contig, start, end);
	return segment;
}

void destructSegment(struct Segment *segment) {
	free(segment->rSegment);
	free(segment);
}

void logSegment(struct Segment *segment) {
	logDebug("Contig : " INT_STRING ", start : " INT_STRING ", end : " INT_STRING "\n",
			segment->contig, segment->start, segment->end);
}

int segmentComparatorPointers(struct Segment **segment1, struct Segment **segment2) {
	return segmentComparator(*segment1, *segment2);
}

int segmentComparator(struct Segment *segment1, struct Segment *segment2) {
	/*
	 * Compares two segments to allow an ordering on any two segments
	 * to be constructued.
	 */
	//Compare the contigs
	if(segment1->contig > segment2->contig) {
		return 1;
	}
	if(segment1->contig < segment2->contig) {
		return -1;
	}
	//Check if overlap.
	if(segment1->start <= segment2->start) {
		if(segment1->end >= segment2->start) {
			return 0;
		}
	}
	else {
		if(segment2->end >= segment1->start) {
			return 0;
		}
	}
	//Back to business
	if(segment1->start < segment2->start) {
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

struct PinchVertex *constructPinchVertex(struct PinchGraph *graph, int32_t vertexID) {
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
	return pinchVertex;
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

int32_t oComparator(const void *o1, const void *o2, void *a) {
	/*
	 * Compares the objects by there address.
	 */
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


struct PinchEdge *constructPinchEdge(struct Segment *segment) {
	struct PinchEdge *pinchEdge;
	struct PinchEdge *rPinchEdge;

	pinchEdge = mallocLocal(sizeof(struct PinchEdge));
	rPinchEdge = mallocLocal(sizeof(struct PinchEdge));

	pinchEdge->rEdge = rPinchEdge;
	rPinchEdge->rEdge = pinchEdge;

	pinchEdge->segment = segment;
	rPinchEdge->segment = segment->rSegment;

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
	//Destroys the contained segment but not the rEdge, so must be called on both.
	free(pinchEdge->rEdge);
	destructSegment(pinchEdge->segment);
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
	 * Compares the segments by their segments.
	 */
	return segmentComparator(edge1->segment, edge2->segment);
}

struct PinchEdge *getContainingBlackEdge(struct PinchGraph *graph, int32_t contig, int32_t position) {
	/*
	 * Gets edge containing the considered segment.
	 */
	static struct PinchEdge edge;
	struct PinchEdge *edge2;
	static struct Segment segment;
	//Edge is a wrapper for the segment - this is all a bit silly.
	edge.segment = &segment;
	segment.contig = contig;
	segment.start = position;
	segment.end = position;
	//Now get the edge.
	edge2 = avl_find(graph->edges, &edge);
	return edge2;
}

struct PinchEdge *getNextEdge(struct PinchGraph *graph, struct PinchEdge *edge) {
	/*
	 * Gets the next in the sequence from the given edge, used for traversing paths in the pinch graph.
	 */
	struct PinchEdge *edge2;

	edge2 = getContainingBlackEdge(graph, edge->segment->contig - (edge->segment->contig % 3) + 2, edge->segment->end+1);
	if(edge2 == NULL) {
		edge2 = getContainingBlackEdge(graph, edge->segment->contig - (edge->segment->contig % 3) + 1, edge->segment->end+1);
	}
	if(edge2 == NULL) {
		edge2 = getContainingBlackEdge(graph, edge->segment->contig - (edge->segment->contig % 3), edge->segment->end+1);
	}
#ifdef BEN_DEBUG
	assert(edge2 != NULL);
#endif
	return edge2;
}

void splitEdge_P(struct PinchGraph *graph,
				 struct PinchEdge *edge, int32_t position,
				 struct PinchVertex *vertex1, struct PinchVertex *vertex2) {
	struct PinchEdge *edge1;
	struct PinchEdge *edge2;

#ifdef BEN_DEBUG
	assert(edge->segment->start < position);
#endif

	//Otherwise create new pinch point.
	//Split segment
	edge1 = constructPinchEdge(constructSegment(edge->segment->contig,
								edge->segment->start, position-1));
	connectPinchEdge(edge1, edge->from, vertex1);

	edge2 = constructPinchEdge(constructSegment(edge->segment->contig,
								position, edge->segment->end));
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

struct PinchVertex *splitEdge(struct PinchGraph *graph, int32_t contig,
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
	if(edge->segment->start == position && leftOrRight == LEFT) {
		return edge->from;
	}

	//If pinch point already at the end..
	if(edge->segment->end == position && leftOrRight == RIGHT) {
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

#ifdef BEN_DEBUG
	//check every edge has the same absolute length.
	blackEdgeIterator = getBlackEdgeIterator(edge->from);
	while((edge2 = getNextBlackEdge(edge->from, blackEdgeIterator)) != NULL) {
		assert(edge2->segment->end - edge2->segment->start == edge->segment->end - edge->segment->start);
	}
	destructBlackEdgeIterator(blackEdgeIterator);
#endif

	//Construct the vertices these things will be attached to.
	vertex1 = constructPinchVertex(graph, -1);
	vertex2 = constructPinchVertex(graph, -1);
	//Add grey edges to connect new vertices.
	connectVertices(vertex1, vertex2);

	j = position - edge->segment->start + (leftOrRight == LEFT ? 0 : 1);
	assert(j > 0);
	for(i=0; i<sE_list->length; i++) {
		edge2 = sE_list->list[i];
#ifdef BEN_DEBUG
		assert(j + edge2->segment->start <= edge2->segment->end);
#endif
		splitEdge_P(graph, edge2, j + edge2->segment->start, vertex1, vertex2);
	}

	//Return either left or right new vertex, dependent on which was made.
	return leftOrRight == LEFT ? vertex2 : vertex1;
}

int32_t isAStubOrCap(struct PinchEdge *edge) {
	return edge->segment->contig % 3 != 1;
}

int32_t isADeadEnd(struct PinchVertex *vertex) {
	/*
	 * Is dead end if is sink vertex, is attached by a grey edge to the sink vertex
	 * or has no grey edges.
	 */
	return lengthBlackEdges(vertex) == 0 ||
	lengthGreyEdges(vertex) == 0 ||
	(lengthGreyEdges(vertex) == 1 && getFirstGreyEdge(vertex)->vertexID == 0);
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

	pinchGraph->edges = avl_create((int32_t (*)(const void *, const void *, void *a))edgeComparator, NULL, NULL);
	pinchGraph->vertices = constructEmptyList(0, (void (*)(void *))destructPinchVertex);
	constructPinchVertex(pinchGraph, -1);

	return pinchGraph;
}

void destructPinchGraph_1(struct PinchEdge *edge, void *o) {
	if(edge->from->vertexID > edge->to->vertexID) {
		destructSegment(edge->segment);
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
					    struct Segment *segment) {
	struct PinchVertex *vertex;
	struct PinchEdge *edge;

	resetVertexChain(vertexChain);

	//do any adjustments off the bat
	splitEdge(graph, segment->contig, segment->start, LEFT);
	splitEdge(graph, segment->contig, segment->end, RIGHT);

	vertex = splitEdge(graph, segment->contig, segment->start, LEFT);
	intListAppend(vertexChain->coordinates, 0);
	intListAppend(vertexChain->leftsOrRights, LEFT);
	listAppend(vertexChain->listOfVertices, vertex);

	//follow chain to get remaining vertices
	edge = getContainingBlackEdge(graph, segment->contig, segment->start);
	while(edge->segment->end < segment->end) {
		intListAppend(vertexChain->coordinates, edge->segment->end - segment->start);
		intListAppend(vertexChain->leftsOrRights, RIGHT);
		listAppend(vertexChain->listOfVertices, edge->to);

		edge = getContainingBlackEdge(graph, segment->contig, edge->segment->end+1);
		intListAppend(vertexChain->coordinates, edge->segment->start - segment->start);
		intListAppend(vertexChain->leftsOrRights, LEFT);
		listAppend(vertexChain->listOfVertices, edge->from);
	}

	//add the second vertex.
	vertex = splitEdge(graph, segment->contig, segment->end, RIGHT);
	intListAppend(vertexChain->coordinates, segment->end - segment->start);
	intListAppend(vertexChain->leftsOrRights, RIGHT);
	listAppend(vertexChain->listOfVertices, vertex);

#ifdef BEN_DEBUG
	//now return the list (the other lists are assigned implicitly)!
	assert(vertexChain->listOfVertices->length == vertexChain->coordinates->length);
	assert(vertexChain->leftsOrRights->length == vertexChain->listOfVertices->length);
#endif
}

int32_t pinchMergeSegment_P(struct VertexChain *vertexChain1,
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

void pinchMergeSegment_getChainOfVertices(struct PinchGraph *graph,
										  struct Segment *segment1,
										  struct Segment *segment2,
										  struct VertexChain *vertexChain1,
										  struct VertexChain *vertexChain2) {
	int32_t i, j, k;

	getChainOfVertices(vertexChain1, graph, segment1);
	getChainOfVertices(vertexChain2, graph, segment2);

	while(pinchMergeSegment_P(vertexChain1, vertexChain2) == FALSE) {
		/*
		 * match up the set of vertices for each chain
		 */
		for(i=0; i<vertexChain1->coordinates->length; i++) {
			j = vertexChain1->coordinates->list[i];
			k = vertexChain1->leftsOrRights->list[i];
			//now search if there is an equivalent vertex.
			splitEdge(graph, segment2->contig, segment2->start + j, k);
		}
		for(i=0; i<vertexChain2->coordinates->length; i++) {
			j = vertexChain2->coordinates->list[i];
			k = vertexChain2->leftsOrRights->list[i];
			//now search if there is an equivalent vertex.
			splitEdge(graph, segment1->contig, segment1->start + j, k);
		}

		getChainOfVertices(vertexChain1, graph, segment1);
		getChainOfVertices(vertexChain2, graph, segment2);
	}
}


struct VertexChain *pMS_vertexChain1 = NULL;
struct VertexChain *pMS_vertexChain2 = NULL;

void pinchMergeSegment(struct PinchGraph *graph,
					   struct Segment *segment1,
					   struct Segment *segment2) {
	/*
	 * Pinches the graph (with the minimum number of required pinches, to
	 * represent the contiguous alignment of the two segments.
	 *
	 * Segments have to be of equal length.
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
	 * Check segments are of the same length (the current (temporary assumption))
	 */
#ifdef BEN_DEBUG
	assert(segment1->end - segment1->start == segment2->end - segment2->start);
#endif

	/*
	 * run through each chain finding the list of vertices.
	 */
	splitEdge(graph, segment1->contig, segment1->start, LEFT);
	splitEdge(graph, segment1->contig, segment1->end, RIGHT);
	splitEdge(graph, segment2->contig, segment2->start, LEFT);
	splitEdge(graph, segment2->contig, segment2->end, RIGHT);

	pinchMergeSegment_getChainOfVertices(graph, segment1, segment2, pMS_vertexChain1, pMS_vertexChain2);

#ifdef BEN_DEBUG
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
	 * Merge the lists of vertices to do the final merge.
	 */
	for(i=0; i<pMS_vertexChain1->listOfVertices->length;) {
		vertex1 = pMS_vertexChain1->listOfVertices->list[i];
		vertex2 = pMS_vertexChain2->listOfVertices->list[i];

#ifdef BEN_DEBUG
		assert(lengthBlackEdges(vertex1) > 0);
#endif
		//check if the two vertices are the ends of same segment.
		edge = getFirstBlackEdge(vertex1);
		if(edge->to == vertex2) {
			//if edge segment is of length greater than one.
			if(edge->segment->end - edge->segment->start > 0) {
				j = (edge->segment->end - edge->segment->start + 1) / 2 + edge->segment->start - 1;
				k = 1 + ((edge->segment->end - edge->segment->start + 1) % 2);

				contig = edge->segment->contig;
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
				pinchMergeSegment_getChainOfVertices(graph, segment1, segment2, pMS_vertexChain1, pMS_vertexChain2);
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
			vertex3 = mergeVertices(graph, vertex1, vertex2);
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

#ifdef BEN_DEBUG
	/*
	 * Do debug checks that the merge went okay.
	 */

	getChainOfVertices(pMS_vertexChain1, graph, segment1);
	getChainOfVertices(pMS_vertexChain2, graph, segment2);

	for(i=0; i<pMS_vertexChain1->listOfVertices->length; i++) {
		vertex1 = pMS_vertexChain1->listOfVertices->list[i];
		void *blackEdgeIterator = getBlackEdgeIterator(vertex1);
		edge = getNextBlackEdge(vertex1, blackEdgeIterator);
		k = edge->segment->end - edge->segment->start;
		vertex2 = edge->to;
		while((edge = getNextBlackEdge(vertex1, blackEdgeIterator)) != NULL) {
			assert(edge->segment->end - edge->segment->start == k);
			assert(edge->to == vertex2);
		}
		destructBlackEdgeIterator(blackEdgeIterator);
	}

	for(i=0; i<pMS_vertexChain2->listOfVertices->length; i++) {
		vertex1 = pMS_vertexChain2->listOfVertices->list[i];
		void *blackEdgeIterator = getBlackEdgeIterator(vertex1);
		edge = getNextBlackEdge(vertex1, blackEdgeIterator);
		k = edge->segment->end - edge->segment->start;
		vertex2 = edge->to;
		while((edge = getNextBlackEdge(vertex1, blackEdgeIterator)) != NULL) {
			assert(edge->segment->end - edge->segment->start == k);
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

void pinchMerge(struct PinchGraph *graph, struct PairwiseAlignment *pA, struct hashtable *contigStringToContigIndex) {
	/*
	 * Method to pinch together the graph using all the aligned matches in the
	 * input alignment.
	 */
	int32_t i, j, k, contig1, contig2;
	struct AlignmentOperation *op;
	static struct Segment segment1;
	static struct Segment rSegment1;
	static struct Segment segment2;
	static struct Segment rSegment2;

	//links the static segments
	segment1.rSegment = &rSegment1;
	rSegment1.rSegment = &segment1;

	segment2.rSegment = &rSegment2;
	rSegment2.rSegment = &segment2;

	j = pA->start1+2;
	k = pA->start2+2;

	contig1 = pinchMerge_getContig(pA->contig1, j, contigStringToContigIndex);
	contig2 = pinchMerge_getContig(pA->contig2, k, contigStringToContigIndex);

	for(i=0; i<pA->operationList->length; i++) {
		op = pA->operationList->list[i];
		if(op->opType == PAIRWISE_MATCH) {
			if(pA->strand1 == '+') {
				segment_recycle(&segment1, contig1, j, j+op->length-1);
			}
			else {
				segment_recycle(&segment1, contig1, -(j+op->length-1), -j);
			}
			if(pA->strand2 == '+') {
				segment_recycle(&segment2, contig2, k, k+op->length-1);
			}
			else {
				segment_recycle(&segment2, contig2, -(k+op->length-1), -k);
			}
			pinchMergeSegment(graph, &segment1, &segment2);
		}
		if(op->opType != PAIRWISE_INDEL_Y) {
			j += op->length;
		}
		if(op->opType != PAIRWISE_INDEL_X) {
			k += op->length;
		}
	}
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for writing out the basic pinch graph.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

char *removeInstance(const char *name) {
	/*
	 * Returns the name of the element from a element instance string.
	 *
	 * All element instances are of the form element_name.instance_id, we therefore
	 * simply remove the .instance_id suffix
	 *
	 */
	char *cA;
	char *cA2;
	int32_t i;

	cA = (char *)name;
	cA2 = (char *)mallocLocal(sizeof(char)*(strlen(cA)+1));

	i=0;
	while(cA[i] != '.' && cA[i] != '\0') {
		cA2[i] = cA[i];
		i++;
	}

#ifdef BEN_DEBUG
	assert(i <= (int32_t)strlen(cA));
#endif

	cA2[i] = '\0';
	return cA2;
}

char *getInstance(const char *name) {
	/*
	 * Returns the instance name of the element instance from a element instance string.
	 *
	 * All element instances are of the form element_name.instance_id, we therefore
	 * simply remove the element_name. prefix
	 *
	 */
	char *cA;
	char *cA2;
	char *cA3;
#ifdef BEN_DEBUG
	assert(name != NULL);
#endif
	cA = stringCopy(name);
	cA2 = cA;
#ifdef BEN_DEBUG
	assert(strlen(cA2) > 0);
#endif
	while(cA2[0] != '.') {
		cA2++;
#ifdef BEN_DEBUG
		assert(strlen(cA2) > 0);
#endif
	}
#ifdef BEN_DEBUG
	assert(cA2[0] == '.');
#endif
	cA3 = stringCopy(cA2+1);
	free(cA);
	return cA3;
}

struct hashtable *getNames(struct PinchGraph *pinchGraph, struct List *contigIndexToContigStrings, const char *atomNamePrefix) {
	/*
	 * Constructs names for each element (vertex and edge).
	 *
	 * Every edge is given a instance name, which is either an atom instance name, or an existing cap instance name,
	 * depending on if the edge is part of an atom or cap.
	 *
	 * Every vertex except 'dead' vertices are given cap names.
	 */
	int32_t i, j;
	struct PinchVertex *vertex;
	struct PinchEdge *edge;
	struct hashtable *names;
	char *cA;

	names = create_hashtable(pinchGraph->vertices->length*10,
							 hashtable_key, hashtable_equalKey,
							 NULL, free);

	for(i=0; i<pinchGraph->vertices->length; i++) { //excludes the sink vertex.
		vertex = (struct PinchVertex *)pinchGraph->vertices->list[i];

		if(isADeadEnd(vertex) == FALSE) {
#ifdef BEN_DEBUG
			assert(lengthBlackEdges(vertex) > 0);
#endif
			//for each edge
			edge = getFirstBlackEdge(vertex);
			if(isAStubOrCap(edge)) { //if the vertex is part of a stub or cap return the stub/cap name
				void *blackEdgeIterator = getBlackEdgeIterator(vertex);
				while((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) { //for each edge create atom/cap instance name
#ifdef BEN_DEBUG
					assert(isAStubOrCap(edge) == TRUE);
					assert(hashtable_search(names, edge) == NULL);
					assert(hashtable_search(names, edge->rEdge) == NULL);
#endif
					cA = stringCopy((char *)contigIndexToContigStrings->list[edge->segment->contig]);
					hashtable_insert(names, edge, cA);
				}
				destructBlackEdgeIterator(blackEdgeIterator);
			}
			else {
				if(edge->from->vertexID < edge->to->vertexID) {
					j = 0;
					void *blackEdgeIterator = getBlackEdgeIterator(vertex);
					while((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) { //for each edge create atom/cap instance name
#ifdef BEN_DEBUG
						assert(edge->from->vertexID < edge->to->vertexID);
						assert(isAStubOrCap(edge) == FALSE);
						assert(hashtable_search(names, edge) == NULL);
						assert(hashtable_search(names, edge->rEdge) == NULL);
#endif
						//return a name using the given prefix plus the smaller vertex of the atom
						cA = (char *)mallocLocal(sizeof(char)*(50 + strlen(atomNamePrefix)));
						sprintf(cA, "%s%i.%i", atomNamePrefix, edge->from->vertexID, j);
						hashtable_insert(names, edge, cA);

						cA = (char *)mallocLocal(sizeof(char)*(50 + strlen(atomNamePrefix)));
						sprintf(cA, "-%s%i.%i", atomNamePrefix, edge->from->vertexID, j);
						hashtable_insert(names, edge->rEdge, cA);
						j++;
					}
					destructBlackEdgeIterator(blackEdgeIterator);
				}
			}

			//for each vertex get the cap name, unless is wrong end of stub/cap, else do not add.
			edge = getFirstBlackEdge(vertex);
			if(isAStubOrCap(edge)) {
				cA = removeInstance((char *)contigIndexToContigStrings->list[edge->segment->contig]);
			}
			else {
				cA = (char *)mallocLocal(sizeof(char)*(50 + strlen(atomNamePrefix)));
				if(edge->from->vertexID < edge->to->vertexID) {
					sprintf(cA, "%s%i]", atomNamePrefix, edge->from->vertexID);
				}
				else {
					sprintf(cA, "[%s%i", atomNamePrefix, edge->to->vertexID);
				}
			}
#ifdef BEN_DEBUG
			assert(hashtable_search(names, vertex) == NULL);
#endif
			hashtable_insert(names, vertex, cA);
		}
	}

	return names;
}

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

void writeOutPinchGraphWithChains(struct PinchGraph *pinchGraph,
								  struct List *biConnectComponentsList,
								  struct List *adjacencyComponents,
								  struct List *contigIndexToContigStrings,
								  const char *uniqueNamePrefix,
								  FILE *fileHandle) {
	/*
	 * Writes out a graph in 'dot' format, compatible with graphviz.
	 *
	 * The associated function 'writeOutEdgeSegments' gives a way of decoding
	 * the segments associated with the black (segment containing) edges
	 * of the graph.
	 */
	int32_t i, j, k;
	struct PinchVertex *vertex;
	struct PinchVertex *vertex2;
	struct PinchEdge *edge;
	struct hashtable *hash;
	struct List *chains;
	struct List *adjacencyComponent;
	struct List *biConnectedComponent;
	struct CactusEdge *cactusEdge;
	struct Segment *segment;
	char *colour;
	char *name;
	struct hashtable *names;

	logDebug("Writing the pinch graph\n");

	names = getNames(pinchGraph, contigIndexToContigStrings, uniqueNamePrefix);

	//Put the chain segments in a hash to colour the black edges of the pinch graph.
	hash = create_hashtable(pinchGraph->vertices->length*10,
								hashtable_key, hashtable_equalKey,
							   NULL, (void (*)(void *))destructInt);

	if(biConnectComponentsList != NULL) {
		for(i=0; i<biConnectComponentsList->length;i++) {
			biConnectedComponent = biConnectComponentsList->list[i];
			for(k=0; k<biConnectedComponent->length; k++) {
				cactusEdge = biConnectedComponent->list[k];
				for(j=0; j<cactusEdge->segments->length; j++) {
					segment = cactusEdge->segments->list[j];
					hashtable_insert(hash, segment, constructInt(i));
					hashtable_insert(hash, segment->rSegment, constructInt(i));
				}
			}
		}
	}

	if(adjacencyComponents != NULL) {
		for(i=0; i<adjacencyComponents->length; i++) {
			adjacencyComponent = adjacencyComponents->list[i];
			for(j=0; j<adjacencyComponent->length; j++) {
				vertex = adjacencyComponent->list[j];
				hashtable_insert(hash, vertex, constructInt(i));
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
		colour = getColour(hash, vertex);
		fprintf(fileHandle, "node[width=0.3,height=0.3,shape=circle,colour=%s,fontsize=14];\n", colour);
		name = (char *)hashtable_search(names, vertex);
		if(name == NULL) {
			fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n", vertex->vertexID, vertex->vertexID);
		}
		else {
			fprintf(fileHandle, "n" INT_STRING "n [label=\"%s\"];\n", vertex->vertexID, name);
		}
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
				colour = getColour(hash, edge->segment);
				fprintf(fileHandle, "edge[color=%s,len=2.5,weight=100,dir=forward];\n", colour);

				void *blackEdgeIterator = getBlackEdgeIterator(vertex);
				while((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
					name = (char *)hashtable_search(names, edge);
					if(name == NULL) {
						name = (char *)hashtable_search(names, edge->rEdge);
					}
#ifdef BEN_DEBUG
					assert(name != NULL);
#endif
					fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n [label=\"" INT_STRING "." INT_STRING ":" INT_STRING ":%s\"];\n",
							edge->from->vertexID, edge->to->vertexID, edge->segment->contig, edge->segment->start, edge->segment->end, name);
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

	destructList(chains);
	hashtable_destroy(hash, TRUE, FALSE);
	hashtable_destroy(names, TRUE, FALSE);

	logDebug("Written the pinch graph\n");
}
