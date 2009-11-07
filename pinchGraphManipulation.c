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

#include "pinchGraph.h"

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


struct PinchEdge *removeTrivialGreyEdge(struct PinchGraph *graph, struct PinchEdge *edge1, struct PinchEdge *edge2) {
	/*
	 * Removes a trivial grey edge component of two vertices from the graph.
	 */
	struct PinchEdge *edge3;
	struct PinchVertex *vertex1;
	struct PinchVertex *vertex2;

#ifdef BEN_DEBUG
	assert(lengthGreyEdges(edge1->to) == 1);
	assert(lengthBlackEdges(edge1->to) == 1);
	assert(lengthGreyEdges(edge2->from) == 1);
	assert(lengthBlackEdges(edge2->from) == 1);
	assert(edge1->segment->contig == edge2->segment->contig);
	assert(edge1->segment->end+1 == edge2->segment->start);
#endif

	edge3 = constructPinchEdge(constructSegment(edge1->segment->contig, edge1->segment->start, edge2->segment->end));
	connectPinchEdge(edge3, edge1->from, edge2->to);

	//Remove the old edges
	vertex1 = edge1->to;
	vertex2 = edge2->from;
	removePinchEdgeFromGraphAndDestruct(graph, edge1);
	removePinchEdgeFromGraphAndDestruct(graph, edge2);

	//Destruct the old vertices, after destructing the edges 1 and 2.
	removeVertexFromGraphAndDestruct(graph, vertex1);
	removeVertexFromGraphAndDestruct(graph, vertex2);

	//Add the new pinch edge to the graph after removing the old edges from the graph.
	addPinchEdgeToGraph(graph, edge3);

	return edge3;
}

void removeTrivialGreyEdgeComponents(struct PinchGraph *graph, struct List *listOfVertices) {
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
	for(i=0; i<listOfVertices->length; i++) {
		vertex1 = listOfVertices->list[i];
		if(lengthBlackEdges(vertex1) == 1 && lengthGreyEdges(vertex1) == 1) {
			edge1 = getFirstBlackEdge(vertex1);;
			vertex2 = getFirstGreyEdge(vertex1);
			if(lengthBlackEdges(vertex2) == 1 && lengthGreyEdges(vertex2) == 1) {
				edge2 = getFirstBlackEdge(vertex2);
				if(!isAStubOrCap(edge1) && !isAStubOrCap(edge2)) {
					if(vertex1->vertexID < vertex2->vertexID) { //Avoid treating self loops (equal) and dealing with trivial grey components twice.
						listAppend(list, vertex1);
					}
				}
			}
		}
	}

	//Remove the trivial components.
	for(i=0; i<list->length; i++) {
		vertex1 = list->list[i];
#ifdef BEN_DEBUG
		assert(lengthGreyEdges(vertex1) == 1);
		assert(lengthBlackEdges(vertex1) == 1);
#endif
		vertex2 = getFirstGreyEdge(vertex1);
#ifdef BEN_DEBUG
		assert(lengthGreyEdges(vertex2) == 1);
		assert(lengthBlackEdges(vertex2) == 1);
#endif
		removeTrivialGreyEdge(graph, getFirstBlackEdge(vertex1)->rEdge, getFirstBlackEdge(vertex2));
	}

	//cleanup
	destructList(list);
}

void splitMultipleBlackEdgesFromVertex(struct PinchGraph *pinchGraph, struct PinchVertex *vertex,
		struct List *newVerticesList, Net *net) {
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
	while(lengthBlackEdges(vertex) > 0) {
		edge = getFirstBlackEdge(vertex);
		//first find the grey edge to attach to the new vertex we're about to create
		vertex3 = getNextEdge(pinchGraph, edge->rEdge, net)->from;
		listAppend(list, vertex3); //can't detach the old vertices yet

		assert(popBlackEdge(vertex) == edge); //detaches edge from vertex.
#ifdef BEN_DEBUG
		assert(!isAStubOrCap(edge));
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
	for(j=0; j<list->length; j++) {
		vertex3 = list->list[j];
		if(containsGreyEdge(vertex3, vertex)) { //it may have already been detached.
			removeGreyEdge(vertex3, vertex);
		}
	}
	//now remove the old vertex
	removeVertexFromGraphAndDestruct(pinchGraph, vertex);
	destructList(list);
}

void removeOverAlignedEdges(struct PinchGraph *pinchGraph, int32_t degree, Net *net) {
	/*
	 * Method splits black edges from the graph with degree higher than a given number of sequences.
	 */
	int32_t i;
	struct List *list;
	struct List *list2;
	struct PinchVertex *vertex;
	struct PinchVertex *vertex2;

	list = constructEmptyList(0, NULL);
	list2 = constructEmptyList(0, NULL);

	for(i=0; i<pinchGraph->vertices->length; i++) {
		vertex = pinchGraph->vertices->list[i];
		if(lengthBlackEdges(vertex) > degree && isAStubOrCap(getFirstBlackEdge(vertex)) == FALSE) { //has a high degree and is not a stub/cap
			vertex2 = getFirstBlackEdge(vertex)->to;
			if(vertex->vertexID < vertex2->vertexID) {
				listAppend(list, vertex);
			}
		}
	}

	for(i=0; i<list->length; i++) {
		vertex = list->list[i];
		vertex2 = getFirstBlackEdge(vertex)->to;
		list2->length = 0;
		splitMultipleBlackEdgesFromVertex(pinchGraph, vertex, list2, net);
		splitMultipleBlackEdgesFromVertex(pinchGraph, vertex2, list2, net);
		removeTrivialGreyEdgeComponents(pinchGraph, list2); //now get rid of any trivial components
	}

	destructList(list);
	destructList(list2);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method for getting graph components
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void getRecursiveComponents_P(struct PinchVertex *vertex, int32_t (*excludedEdgesFn)(void *), struct hashtable *seen,
							  struct List *component, struct List *stack) {
	//now we are going to search from the vertex to find its associated component.
	struct PinchEdge *edge;
	struct PinchVertex *vertex2;
	//add it to the hash
	assert(hashtable_search(seen, vertex) == NULL);
	hashtable_insert(seen, vertex, vertex);

	//add it to the component
	listAppend(component, vertex);

	//search across blackedges.
	if(lengthBlackEdges(vertex) > 0) {
		edge = getFirstBlackEdge(vertex);
		vertex2 = edge->to;
		if(!excludedEdgesFn(edge)) {
			if(hashtable_search(seen, vertex2) == NULL) {
				listAppend(stack, vertex2);
				//getRecursiveComponents_P(vertex2, excludedEdges, seen, component);
			}
#ifdef BEN_ULTRA_DEBUG
			else {
				//check vertex2 is already in the component.
				int32_t k;
				int32_t l = FALSE;
				for(k=0; k<component->length; k++) {
					if(component->list[k] == vertex2) {
						l = TRUE;
					}
				}
				assert(l == TRUE);
			}
#endif
		}
	}

	//search across the greyedges.
	void *greyEdgeIterator=getGreyEdgeIterator(vertex);
	while((vertex2 = getNextGreyEdge(vertex, greyEdgeIterator)) != NULL) {
		if(hashtable_search(seen, vertex2) == NULL) {
			listAppend(stack, vertex2);
			//getRecursiveComponents_P(vertex2, excludedEdges, seen, component);
		}
#ifdef BEN_ULTRA_DEBUG
		else {
			//check vertex2 is already in the component.
			int32_t k;
			int32_t l = FALSE;
			for(k=0; k<component->length; k++) {
				if(component->list[k] == vertex2) {
					l = TRUE;
				}
			}
			assert(l == TRUE);
		}
#endif
	}
	destructGreyEdgeIterator(greyEdgeIterator);
}

struct List *getRecursiveComponents(struct PinchGraph *pinchGraph, int32_t (*excludedEdgesFn)(void *)) {
	/*
	 * find recursive components (each recursive component is represented as a series of vertices)
	 * do as series of DFS on each connected component, not traversing long black edges.
	 */
	struct hashtable *seen;
	int32_t i;
	struct List *components;
	struct List *component;
	struct PinchVertex *vertex;
	struct PinchVertex *vertex2;
	struct List *stack;

	//allocate stuff.
	seen = create_hashtable(0, hashtable_key, hashtable_equalKey, NULL, NULL);
	components = constructEmptyList(0, (void (*)(void *))destructList);
	stack = constructEmptyList(0, NULL);

	for(i=0; i<pinchGraph->vertices->length; i++) {
		vertex = pinchGraph->vertices->list[i];

		//if not seen
		if(hashtable_search(seen, vertex) == NULL) {
			//get component
			component = constructEmptyList(0, NULL);
			//add to final components
			listAppend(components, component);

			//now build the component via a (recursive) function -- now changed to have a stack which we edit to maintain the list.
			stack->length = 0;
			listAppend(stack, vertex);
			while(stack->length > 0) {
				vertex2 = stack->list[--stack->length];
				if(hashtable_search(seen, vertex2) == NULL) {
					getRecursiveComponents_P(vertex2, excludedEdgesFn, seen, component, stack);
				}
			}
		}
	}
	//clean up
	hashtable_destroy(seen, FALSE, FALSE);
	destructList(stack);

	return components;
}

struct hashtable *getRecursiveComponents2_excludedEdgesHash;
int32_t getRecursiveComponents2_excludedEdgesFn(void *o) {
	struct PinchEdge *edge;
	static int32_t iA[2];

	edge = o;
	iA[0] = edge->from->vertexID;
	iA[1] = edge->to->vertexID;
	return hashtable_search(getRecursiveComponents2_excludedEdgesHash, iA) != NULL;
}

struct List *getRecursiveComponents2(struct PinchGraph *pinchGraph, struct List *edgesToExclude) {
	/*
	 * Gets the adjacency components, given the list of edges to exlude
	 */
	struct PinchEdge *edge;
	struct List *adjacencyComponents;
	int32_t i;
	int32_t *iA;

	getRecursiveComponents2_excludedEdgesHash = create_hashtable(0, hashtable_intPairHashKey,
			hashtable_intPairEqualKey, (void (*)(void *))destructIntPair, NULL);
	for(i=0; i<edgesToExclude->length; i++) { //build the excluded edges hash.
		edge = edgesToExclude->list[i];
		iA = constructIntPair(edge->from->vertexID, edge->to->vertexID);
		if(hashtable_search(getRecursiveComponents2_excludedEdgesHash, iA) == NULL) {
			hashtable_insert(getRecursiveComponents2_excludedEdgesHash, iA, iA);
		}
		else {
			destructIntPair(iA);
		}
	}
	adjacencyComponents = getRecursiveComponents(pinchGraph, getRecursiveComponents2_excludedEdgesFn);
	hashtable_destroy(getRecursiveComponents2_excludedEdgesHash, FALSE, TRUE);
	return adjacencyComponents;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method for linking the stub components to the
//sink component.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t linkStubComponentsToTheSinkComponent_excludedEdgesFn(void *o) {
	assert(o != NULL);
	return FALSE;
}

void linkStubComponentsToTheSinkComponent(struct PinchGraph *pinchGraph) {
	struct List *components;
	struct List *component;
	struct PinchVertex *vertex;
	struct PinchVertex *sinkVertex;
	int32_t i, j, k;

	//isolate the separate graph components using the components method
	components = getRecursiveComponents(pinchGraph, linkStubComponentsToTheSinkComponent_excludedEdgesFn);

	sinkVertex = pinchGraph->vertices->list[0];
	//for each non-sink component select a random stub to link to the sink vertex.
	k = 0;
	for(i=0; i<components->length; i++) {
		component = components->list[i];
		if(!listContains(component, sinkVertex)) {
			for(j=0; j<component->length; j++) {
				vertex = component->list[j];
				if(isAStubOrCap(getFirstBlackEdge(vertex)) && lengthGreyEdges(vertex) == 0) {
					k++;
					connectVertices(vertex, sinkVertex);
					break;
				}
			}
		}
	}

#ifdef BEN_DEBUG
	assert(k == components->length-1);
#endif

	//clean up
	destructList(components);
}
