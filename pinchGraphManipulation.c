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

void removeTrivialGreyEdge(struct PinchGraph *graph, struct PinchVertex *vertex1, struct PinchVertex *vertex2, Net *net) {
	assert(lengthBlackEdges(vertex1) == lengthBlackEdges(vertex2));
	assert(lengthGreyEdges(vertex1) == 1);
	assert(lengthGreyEdges(vertex2) == 1);
	assert(getFirstGreyEdge(vertex1) == vertex2);
	assert(getFirstGreyEdge(vertex2) == vertex1);

	//if(lengthBlackEdges(vertex1) > 1) {
	//	assert(FALSE);
	//}

	//For each black edge to vertex1 find consecutive edge from vertex2, then join.
	while(lengthBlackEdges(vertex1) > 0) {
		assert(lengthBlackEdges(vertex1) == lengthBlackEdges(vertex2));

		struct PinchEdge *edge1 = getFirstBlackEdge(vertex1);
		assert(edge1 != NULL);
		assert(!isAStub(edge1));
		edge1 = edge1->rEdge;
		assert(edge1->to == vertex1);
		//first find the grey edge to attach to the new vertex we're about to create

		struct PinchEdge *edge2 = getNextEdge(graph, edge1, net);
		assert(edge2 != NULL);
		assert(!isAStub(edge2));
		assert(edge2->from == vertex2);

		struct PinchEdge *edge3 = constructPinchEdge(constructPiece(edge1->piece->contig, edge1->piece->start, edge2->piece->end));
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

void removeTrivialGreyEdgeComponents(struct PinchGraph *graph, struct List *listOfVertices, Net *net) {
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
		if(lengthGreyEdges(vertex1) == 1 && lengthBlackEdges(vertex1) > 0) {
			edge1 = getFirstBlackEdge(vertex1);
			vertex2 = getFirstGreyEdge(vertex1);
			if(lengthGreyEdges(vertex2) == 1 && lengthBlackEdges(vertex2) > 0) {
				edge2 = getFirstBlackEdge(vertex2);
				if(!isAStub(edge1) && !isAStub(edge2)) {
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
		vertex2 = getFirstGreyEdge(vertex1);
		removeTrivialGreyEdge(graph, vertex1, vertex2, net);
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

void removeOverAlignedEdges_P(struct PinchVertex *vertex, int32_t extensionSteps, struct List *list, struct hashtable *hash) {
	void *greyEdgeIterator = getGreyEdgeIterator(vertex);
	struct PinchVertex *vertex2;
	struct PinchVertex *vertex3;

	int32_t distance = *(int32_t *)hashtable_search(hash, vertex);
	if(distance < extensionSteps) {
		while((vertex2 = getNextGreyEdge(vertex, greyEdgeIterator)) != NULL) {
			if(lengthBlackEdges(vertex2) > 0) {
				struct PinchEdge *edge = getFirstBlackEdge(vertex2);
				int32_t length = edge->piece->end - edge->piece->start + 1;
				if(!isAStub(edge)) {
					int32_t *i = hashtable_search(hash, vertex2);
					vertex3 = edge->to;
					int32_t *j = hashtable_search(hash, vertex3);
					if(i == NULL) {
						assert(j == NULL);
						listAppend(list, vertex2->vertexID > vertex3->vertexID ? vertex3 : vertex2);
						hashtable_insert(hash, vertex2, constructInt(distance));
						hashtable_insert(hash, vertex3, constructInt(distance + length));
					}
					else {
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

void removeOverAlignedEdges(struct PinchGraph *pinchGraph, float minimumTreeCoverage, int32_t maxDegree,
		struct List *extraEdgesToUndo,
		int32_t extensionSteps, Net *net) {
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
	for(i=0; i<pinchGraph->vertices->length; i++) {
		vertex = pinchGraph->vertices->list[i];
		if(lengthBlackEdges(vertex) >= 1 &&
			!isAStub(getFirstBlackEdge(vertex))) {
			if (lengthBlackEdges(vertex) > maxDegree ||
				treeCoverage(vertex, net, pinchGraph) < minimumTreeCoverage) { //has a high degree and is not a stub/cap
				vertex2 = getFirstBlackEdge(vertex)->to;
				if(vertex->vertexID < vertex2->vertexID) {
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
	if(extraEdgesToUndo != NULL) {
		for(i=0; i<extraEdgesToUndo->length; i++) {
			edge = extraEdgesToUndo->list[i];
			if(!isAStub(edge)) {
				if(edge->from->vertexID > edge->to->vertexID) {
					edge = edge->rEdge;
				}
				if(hashtable_search(hash, edge->from) == NULL) {
					assert(hashtable_search(hash, edge->to) == NULL);
					hashtable_insert(hash, edge->from, constructInt(0));
					hashtable_insert(hash, edge->to, constructInt(0));
					listAppend(list, edge->from);
				}
				else {
					assert(hashtable_search(hash, edge->to) != NULL);
				}
			}
		}
	}

	logDebug("Got the initial list of over-aligned black edges to undo, total: %i\n", list->length);

	if(extensionSteps > 0) {
		i = 0, k = 10;
		while(list->length != i || k-- > 0) { //k term to ensure distances have been propagated
			assert(list->length >= i);
			i = list->length;
			list2 = listCopy(list); //just use the vertices in the existing list
			for(j=0; j<list2->length; j++) {
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
	for(i=0; i<list->length; i++) {
		vertex = list->list[i];
		if(lengthBlackEdges(vertex) > 1) {
			listAppend(list2, vertex);
		}
	}
	destructList(list);
	list = list2;

	logDebug("Got the list of black edges to undo, total length: %i!\n", list->length);

	list2 = constructEmptyList(0, NULL);
	for(i=0; i<list->length; i++) {
		vertex = list->list[i];
		vertex2 = getFirstBlackEdge(vertex)->to;
		list2->length = 0;
		splitMultipleBlackEdgesFromVertex(pinchGraph, vertex, list2, net);
		splitMultipleBlackEdgesFromVertex(pinchGraph, vertex2, list2, net);
		removeTrivialGreyEdgeComponents(pinchGraph, list2, net); //now get rid of any trivial components
	}

	destructList(list);
	destructList(list2);
	hashtable_destroy(hash, 1, 0);
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

int32_t getRecursiveComponents_excludedEdgesFn(void *o) {
	assert(o != NULL);
	return TRUE;
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

	if(excludedEdgesFn == NULL) {
		excludedEdgesFn = getRecursiveComponents_excludedEdgesFn;
	}

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
	 * Gets the groups, given the list of edges to exlude
	 */
	struct PinchEdge *edge;
	struct List *groups;
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
	groups = getRecursiveComponents(pinchGraph, getRecursiveComponents2_excludedEdgesFn);
	hashtable_destroy(getRecursiveComponents2_excludedEdgesHash, FALSE, TRUE);
	return groups;
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

void linkStubComponentsToTheSinkComponent(struct PinchGraph *pinchGraph, Net *net) {
	struct List *components;
	struct List *component;
	struct PinchVertex *vertex;
	struct PinchVertex *sinkVertex;
	int32_t i, j, k, l;
	struct PinchEdge *edge;
	Sequence *sequence;
	Sequence *longestSequence;

	//isolate the separate graph components using the components method
	components = getRecursiveComponents(pinchGraph, linkStubComponentsToTheSinkComponent_excludedEdgesFn);

	sinkVertex = pinchGraph->vertices->list[0];
	//for each non-sink component select a random stub to link to the sink vertex.
	k = 0;
	l = 0;
	for(i=0; i<components->length; i++) {
		component = components->list[i];
		if(!listContains(component, sinkVertex)) {
			//Get the longest sequence contained in the component and attach
			//its two ends to the source vertex.
			//Make the the two ends attached end_makeAttached(end) / end_makeUnattached(end)
			longestSequence = NULL;
			for(j=0; j<component->length; j++) {
				vertex = component->list[j];
				if(vertex_isDeadEnd(vertex)) {
					assert(lengthGreyEdges(vertex) == 0);
					assert(lengthBlackEdges(vertex) == 1);
					edge = getFirstBlackEdge(getFirstGreyEdge(getFirstBlackEdge(vertex)->to));
					sequence = net_getSequence(net, edge->piece->contig);
					if(longestSequence == NULL ||
							sequence_getLength(sequence) > sequence_getLength(longestSequence)) {
						longestSequence = sequence;
					}
				}
			}
			//assert(0);
			assert(longestSequence != NULL);
			for(j=0; j<component->length; j++) {
				vertex = component->list[j];
				if(vertex_isDeadEnd(vertex)) {
					assert(lengthGreyEdges(vertex) == 0);
					assert(lengthBlackEdges(vertex) == 1);
					edge = getFirstBlackEdge(getFirstGreyEdge(getFirstBlackEdge(vertex)->to));
					sequence = net_getSequence(net, edge->piece->contig);
					if(sequence == longestSequence) {
						connectVertices(vertex, sinkVertex);
						k++;
					}
				}
			}
		}
	}

#ifdef BEN_DEBUG
	assert(k == 2*(components->length-1));
#endif

	//clean up
	destructList(components);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method for assessing how much of the event tree the
//given set of connected pinch edges covers.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


float treeCoverage(struct PinchVertex *vertex, Net *net,
		struct PinchGraph *pinchGraph) {
	/*
	 * Returns the proportion of the tree covered by the block.
	 */
	struct Piece *piece;
	EventTree *eventTree;
	Event *event;
	Event *commonAncestorEvent;
	struct hashtable *hash;
	float treeCoverage;

#ifdef BEN_DEBUG
	assert(lengthBlackEdges(vertex) > 0);
	assert(!isAStub(getFirstBlackEdge(vertex)));
#endif

	eventTree = net_getEventTree(net);
	commonAncestorEvent = NULL;
	void *blackEdgeIterator = getBlackEdgeIterator(vertex);
	struct PinchEdge *edge;
	while((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
		piece = edge->piece;
		event = sequence_getEvent(net_getSequence(net, piece->contig));
		commonAncestorEvent = commonAncestorEvent == NULL ? event : eventTree_getCommonAncestor(event, commonAncestorEvent);
	}
	destructBlackEdgeIterator(blackEdgeIterator);


	treeCoverage = 0.0;
	hash = create_hashtable(eventTree_getEventNumber(eventTree)*2,
					 hashtable_key, hashtable_equalKey,
					 NULL, NULL);

	blackEdgeIterator = getBlackEdgeIterator(vertex);
	while((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
		piece = edge->piece;
		event = sequence_getEvent(net_getSequence(net, piece->contig));
		while(event != commonAncestorEvent && hashtable_search(hash, event) == NULL) {
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
	treeCoverage /= event_getSubTreeBranchLength(event_getChild(eventTree_getRootEvent(eventTree), 0));
	if(treeCoverage < -0.001) {
		uglyf("The tree coverage for this case is: %f\n", treeCoverage);
	}
	assert(treeCoverage >= -0.001);
	assert(treeCoverage <= 1.0001);
	return treeCoverage;
}
