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

	cactusVertex = mallocLocal(sizeof(struct CactusVertex));
	cactusVertex->edges = constructEmptyList(0, (void (*)(void *))destructCactusEdge);

	return cactusVertex;
}

void destructCactusVertex(struct CactusVertex *vertex) {
	destructList(vertex->edges);
	free(vertex);
}

struct CactusEdge *constructCactusEdge(struct List *segments) {
	struct CactusEdge *edge;
	struct CactusEdge *rEdge;
	struct Segment *segment;
	int32_t i;

	edge = mallocLocal(sizeof(struct CactusEdge));
	rEdge = mallocLocal(sizeof(struct CactusEdge));
	edge->rEdge = rEdge;
	rEdge->rEdge = edge;

	edge->segments = constructEmptyList(segments->length, NULL);
	edge->rEdge->segments = constructEmptyList(segments->length, NULL);

	for(i=0; i<segments->length; i++) {
		segment = segments->list[i];
		edge->segments->list[i] = segment;
		edge->rEdge->segments->list[i] = segment->rSegment;
	}

	return edge;
}

void destructCactusEdge(struct CactusEdge *edge) {
	//the rEdge is dealt with by its respective from vertex.
	destructList(edge->segments);
	free(edge);
}

struct CactusGraph *constructCactusGraph(struct PinchGraph *pinchGraph,
										 struct List *extraEdges,
										 struct List *threeEdgeConnectedComponents) {
	struct CactusGraph *cactusGraph;
	struct CactusEdge *cactusEdge;
	struct CactusVertex *cactusVertex;
	struct CactusVertex *cactusVertex2;
	int32_t i, j, k;
	void *greyEdgeIterator;
	struct List *list;
	struct PinchVertex *pinchVertex;
	struct PinchVertex *pinchVertex2;
	struct PinchEdge *pinchEdge;
	struct List *pinchVertexToCactusVertex;
	struct List *emptyList;
	struct List *list2;
	struct List *extras;

	cactusGraph = mallocLocal(sizeof(struct CactusGraph));
	cactusGraph->vertices = constructEmptyList(0, (void (*)(void *))destructCactusVertex);

#ifdef BEN_DEBUG
	list = threeEdgeConnectedComponents->list[0];
	j = FALSE;
	for(i=0; i<list->length; i++) {
		pinchVertex = list->list[0];
		if(pinchVertex->vertexID == 0) {
			j = TRUE;
		}
	}
	assert(j == TRUE);
#endif

	pinchVertexToCactusVertex = constructEmptyList(pinchGraph->vertices->length, NULL);
	assert(pinchVertexToCactusVertex->length == pinchGraph->vertices->length);
	for(i=0; i<pinchVertexToCactusVertex->length; i++) {
		pinchVertexToCactusVertex->list[i] = NULL;
	}

	for(i=0; i<threeEdgeConnectedComponents->length; i++) {
		cactusVertex = constructCactusVertex();
		cactusVertex->vertexID = i;
		listAppend(cactusGraph->vertices, cactusVertex);

		list = threeEdgeConnectedComponents->list[i];
#ifdef BEN_DEBUG
		//checks all components are non empty.
		assert(list->length > 0);
#endif

		for(j=0; j<list->length; j++) {
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
	for(i=0; i<pinchVertexToCactusVertex->length; i++) {
		assert(pinchVertexToCactusVertex->list[i] != NULL);
	}
#endif

	emptyList = constructEmptyList(0, NULL);
	for(i=0; i<threeEdgeConnectedComponents->length; i++) {
		list = threeEdgeConnectedComponents->list[i];
		cactusVertex = cactusGraph->vertices->list[i];

		for(j=0; j<list->length; j++) {
			pinchVertex = list->list[j];

			//black edges
			if(lengthBlackEdges(pinchVertex) > 0) {
				pinchEdge = getFirstBlackEdge(pinchVertex);
				pinchVertex2 = pinchEdge->to;
				cactusVertex2 = pinchVertexToCactusVertex->list[pinchVertex2->vertexID];

				if(pinchVertex->vertexID < pinchVertex2->vertexID) {
					//getting the ordered segments.
					list2 = constructEmptyList(0, NULL);
					void *blackEdgeIterator = getBlackEdgeIterator(pinchVertex);
					while((pinchEdge = getNextBlackEdge(pinchVertex, blackEdgeIterator)) != NULL) {
						listAppend(list2, pinchEdge->segment);
					}
					destructBlackEdgeIterator(blackEdgeIterator);
					cactusEdge = constructCactusEdge(list2);
					destructList(list2);

					listAppend(cactusVertex->edges, cactusEdge);
					cactusEdge->from = cactusVertex;
					cactusEdge->rEdge->to = cactusVertex;

					listAppend(cactusVertex2->edges, cactusEdge->rEdge);
					cactusEdge->to = cactusVertex2;
					cactusEdge->rEdge->from = cactusVertex2;
				}
			}
#ifdef BEN_DEBUG
			else {
				assert(pinchVertex->vertexID == 0);
			}
#endif

			//grey edges
			greyEdgeIterator=getGreyEdgeIterator(pinchVertex);
			while((pinchVertex2 = getNextGreyEdge(pinchVertex, greyEdgeIterator)) != NULL) {
				cactusVertex2 = pinchVertexToCactusVertex->list[pinchVertex2->vertexID];

				if(cactusVertex != cactusVertex2 && cactusVertex < cactusVertex2) {
					cactusEdge = constructCactusEdge(emptyList);

					listAppend(cactusVertex->edges, cactusEdge);
					cactusEdge->from = cactusVertex;
					cactusEdge->rEdge->to = cactusVertex;

					listAppend(cactusVertex2->edges, cactusEdge->rEdge);
					cactusEdge->to = cactusVertex2;
					cactusEdge->rEdge->from = cactusVertex2;
				}
			}
			destructGreyEdgeIterator(greyEdgeIterator);

			//extra edges
			extras = extraEdges->list[pinchVertex->vertexID];
			for(k=0; k<extras->length; k++) {
				pinchVertex2 = extras->list[k];
				cactusVertex2 = pinchVertexToCactusVertex->list[pinchVertex2->vertexID];

				if(cactusVertex != cactusVertex2 && cactusVertex < cactusVertex2) {
					cactusEdge = constructCactusEdge(emptyList);

					listAppend(cactusVertex->edges, cactusEdge);
					cactusEdge->from = cactusVertex;
					cactusEdge->rEdge->to = cactusVertex;

					listAppend(cactusVertex2->edges, cactusEdge->rEdge);
					cactusEdge->to = cactusVertex2;
					cactusEdge->rEdge->from = cactusVertex2;
				}
			}
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

void checkCactusGraph(struct PinchGraph *pinchGraph,
					  struct List *threeEdgeConnectedComponents,
					  struct CactusGraph *cactusGraph) {
#ifdef BEN_DEBUG
	struct hashtable *hashTable;
	int32_t i, j, k, l, m;
	struct CactusVertex *cactusVertex;
	struct CactusVertex *cactusVertex2;
	struct CactusEdge *cactusEdge;
	struct CactusEdge *cactusEdge2;
	struct PinchVertex *pinchVertex;
	struct PinchVertex *pinchVertex2;
	struct PinchEdge *pinchEdge;
	struct Segment *segment;
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
		assert(list->length > 0);

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
					for(l=0; l<cactusEdge->segments->length; l++) {
						segment = cactusEdge->segments->list[l];
						if(segmentComparator(pinchEdge->segment, segment) == 0) {
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
					if(cactusEdge->segments->length == 0) {
						m += 1;
					}
				}
			}
			if(cactusVertex != cactusVertex2) {
				assert(m > 0);
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

		if(list->length > 1) {
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
#ifdef BEN_DEBUG
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
	logDebug("Constructed the biconnected components for making the net\n");


	////////////////////////////////////////////////
	//(2) Check every component is a cycle.
	////////////////////////////////////////////////

	hashTable = create_hashtable(cactusGraph->vertices->length*10,
								 hashtable_key, hashtable_equalKey,
								 NULL, NULL);

	for(i=0; i<biConnectedComponents->length; i++) {
		list = biConnectedComponents->list[i];

		assert(list->length > 0);

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
	logDebug("Checked that all edges in the cactus graph are in a single simple cycle.\n");
#endif
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method to calculate 2-edge connected components
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void computeBiConnectedComponents_P(struct hashtable *flag,
										  int32_t count,
										  int32_t *dFN,
										  int32_t *low,
										  struct List *stack,
										  int32_t *father,
										  struct List *biConnnectedComponents,
										  struct CactusVertex *v) {
	int32_t i;
	struct CactusVertex *w;
	struct CactusEdge *edge;
	struct List *list;

	hashtable_insert(flag, v, v);
	count++;
	dFN[v->vertexID] = count;
	low[v->vertexID] = count;

	for(i=0; i<v->edges->length; i++) {
		edge = v->edges->list[i];
		w = edge->to;

		assert(hashtable_search(flag, edge) == NULL);

		//self edges are special case
		if(w == v) {
			if(hashtable_search(flag, edge->rEdge) == NULL) {
				hashtable_insert(flag, edge, edge);
				list = constructEmptyList(0, NULL);
				listAppend(list, edge);
				listAppend(biConnnectedComponents, list);
			}
			continue;
		}

		if(hashtable_search(flag, edge->rEdge) == NULL) {
			hashtable_insert(flag, edge, edge);
			listAppend(stack, edge);
		}

		if(hashtable_search(flag, w) == NULL) {
			father[w->vertexID] = v->vertexID;
			computeBiConnectedComponents_P(flag, count, dFN, low, stack, father,
											biConnnectedComponents, w);
			if(low[w->vertexID] >= dFN[v->vertexID]) { //is articulation point
				list = constructEmptyList(0, NULL);
				while(stack->length > 0) {
					listAppend(list, stack->list[--stack->length]);
					if(list->list[list->length-1] == edge) {
						break;
					}
				}
				listAppend(biConnnectedComponents, list);
			}
			low[v->vertexID] = low[v->vertexID] < low[w->vertexID] ? low[v->vertexID] : low[w->vertexID];
		}
		else {
			if(w->vertexID != father[v->vertexID]) {
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

	flag = create_hashtable(cactusGraph->vertices->length*10,
							hashtable_key, hashtable_equalKey,
							NULL, NULL);
	count = 0;
	dFN = mallocLocal(sizeof(int32_t)*cactusGraph->vertices->length);
	low = mallocLocal(sizeof(int32_t)*cactusGraph->vertices->length);
	father = mallocLocal(sizeof(int32_t)*cactusGraph->vertices->length);
	for(i=0; i<cactusGraph->vertices->length; i++) {
		dFN[i] = -1;
		low[i] = -1;
		father[i] = -1;
	}
	stack = constructEmptyList(0, NULL);
	biConnectedComponents = constructEmptyList(0, (void (*)(void *))destructList);

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

int sortBiConnectedComponentsP(struct CactusEdge **edge, struct CactusEdge **edge2) {
	return sortBiConnectedComponents_vertexOrdering[(*edge)->from->vertexID] - sortBiConnectedComponents_vertexOrdering[(*edge2)->from->vertexID];
}

struct List *computeSortedBiConnectedComponents(struct CactusGraph *cactusGraph) {
	struct List *biConnectedComponents;
	struct List *biConnectedComponent;
	int32_t i;

#ifdef BEN_DEBUG
	struct CactusEdge *edge;
	struct CactusEdge *edge2;
	struct CactusVertex *vertex;
	int32_t j;
#endif

	////////////////////////////////////////////////
	//(1) Get bi-connected components.
	////////////////////////////////////////////////
	biConnectedComponents = computeBiConnectedComponents(cactusGraph);
	logDebug("Constructed the biconnected components for making the net\n");

	////////////////////////////////////////////////
	//(1) Get DFS ordering on nodes
	////////////////////////////////////////////////
	sortBiConnectedComponents_vertexOrdering = getDFSDiscoveryTimes(cactusGraph);
	logDebug("Got the DFS discovery times\n");

	////////////////////////////////////////////////
	//(2) Now sort each bi connected component
	////////////////////////////////////////////////
	for(i=0; i<biConnectedComponents->length; i++) {
		biConnectedComponent = biConnectedComponents->list[i];
		qsort(biConnectedComponent->list, biConnectedComponent->length, sizeof(void *),
				(int (*)(const void *v, const void *))sortBiConnectedComponentsP);
#ifdef BEN_DEBUG
		edge = biConnectedComponent->list[0];
		vertex = edge->from;
		for(j=0; j+1<biConnectedComponent->length; j++) {
			edge = biConnectedComponent->list[j];
			edge2 = biConnectedComponent->list[j+1];
			assert(edge2->from->vertexID == edge->to->vertexID);
			assert(edge->to != vertex);
			assert(edge2->from != vertex);
			assert(sortBiConnectedComponents_vertexOrdering[edge->to->vertexID] > sortBiConnectedComponents_vertexOrdering[edge->from->vertexID]);
		}
		edge = biConnectedComponent->list[biConnectedComponent->length-1];
		assert(edge->to->vertexID == vertex->vertexID);
#endif
	}
	logDebug("Sorted each bi-connected component\n");
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
		struct CactusVertex *vertex, int32_t counter,
		int32_t *vertexOrdering) {
	int32_t i;
	struct CactusEdge *edge;

	assert(vertexOrdering[vertex->vertexID] == -1);
	vertexOrdering[vertex->vertexID] = counter++;
	for(i=0; i<vertex->edges->length; i++) {
		edge = vertex->edges->list[i];
		if(vertexOrdering[edge->to->vertexID] == -1) {
			counter = getDFSDiscoveryTimesP(cactusGraph, edge->to, counter, vertexOrdering);
		}
	}
	return counter;
}

int32_t *getDFSDiscoveryTimes(struct CactusGraph *cactusGraph) {
	int32_t i;
	int32_t *vertexOrdering;

	vertexOrdering = mallocLocal(sizeof(int32_t)*cactusGraph->vertices->length);
	for(i=0; i<cactusGraph->vertices->length; i++) {
		vertexOrdering[i] = -1;
	}

	getDFSDiscoveryTimesP(cactusGraph, cactusGraph->vertices->list[0], 0, vertexOrdering);

	for(i=0; i<cactusGraph->vertices->length; i++) {
		assert(vertexOrdering[i] >= 0);
		assert(vertexOrdering[i] < cactusGraph->vertices->length);
	}

	return vertexOrdering;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for dealing with  reverse complement matches
//that create self loops and stubs that create free ends.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct List *getEmptyExtraEdges(struct PinchGraph *pinchGraph) {
	int32_t i;
	struct List *extraEdges;

	//do basic allocation of adjacency list like structure
	extraEdges = constructEmptyList(pinchGraph->vertices->length, (void (*)(void *))destructList);
	for(i=0; i<extraEdges->length; i++) {
		extraEdges->list[i] = constructEmptyList(0, NULL);
	}
	return extraEdges;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//I/O Methods to interact with the 3-edge connected component
//code.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void writeOut3EdgeGraph(struct PinchGraph *pinchGraph, struct List *greyEdgeComponents, struct List *extraEdges,
						FILE *fileHandle) {
	/*
	 * Writes a format compatible with the 3-edge connected component algorithm.
	 */
	struct PinchVertex *vertex;
	struct PinchVertex *vertex2;
	struct PinchEdge *edge;
	int32_t i, j, k, l;
	struct List *list;
	struct List *component;
	struct hashtable *vertexHash;

	//setup vertex to grey edge component hash
	vertexHash = create_hashtable(pinchGraph->vertices->length*2,
			hashtable_key, hashtable_equalKey,
			NULL, (void (*)(void *))destructInt);

	for(i=0;i<greyEdgeComponents->length; i++) {
		component = greyEdgeComponents->list[i];
		for(j=0; j<component->length; j++) {
			vertex = component->list[j];
			hashtable_insert(vertexHash, vertex, constructInt(i));
		}
	}
#ifdef BEN_DEBUG
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
	fprintf(fileHandle, "" INT_STRING "\n", greyEdgeComponents->length);

	for(i=0; i<greyEdgeComponents->length;i++) {
		component = greyEdgeComponents->list[i];
		fprintf(fileHandle, "" INT_STRING "", i+1);
		for(j=0; j<component->length; j++) {
			vertex = component->list[j];

			//The black edges
			if(lengthBlackEdges(vertex) > 0) {
				edge = getFirstBlackEdge(vertex);
				k = *((int32_t *)hashtable_search(vertexHash, edge->to));
				fprintf(fileHandle, ">" INT_STRING "", k+1);
			}
#ifdef BEN_DEBUG
			else {
				assert(vertex->vertexID == 0);
			}
#endif

			//Finally, extra edges.
			list = extraEdges->list[vertex->vertexID];
			for(k=0; k<list->length; k++) {
				vertex2 = list->list[k];
				l = *((int32_t *)hashtable_search(vertexHash, vertex2));
				fprintf(fileHandle, ">" INT_STRING "", l+1);
			}
		}
		//Start a newline for the next component
		fprintf(fileHandle, "\n");
	}
	hashtable_destroy(vertexHash, TRUE, FALSE);
}

struct List *readThreeEdgeComponents(struct PinchGraph *pinchGraph, struct List *greyEdgeComponents, char *file) {
	/*
	 * Reads in the three edge connected components written out by the three
	 * edge script.
	 */
	FILE *fileHandle;
	int32_t i, j, k;
	int32_t numberOfComponents;
	int32_t *iA;
	struct List *list;
	struct List *list2;
	struct List *list3;
	struct List *component;
	struct PinchVertex *vertex;

#ifdef BEN_DEBUG
	int32_t l;
#endif

	//construct stuff
	fileHandle = fopen(file, "r");
	iA = mallocLocal(sizeof(int32_t)*pinchGraph->vertices->length);
	list = constructEmptyList(0, (void (*)(void *))destructList);

	readIntegers(fileHandle, 1, iA);
	numberOfComponents = iA[0];
#ifdef BEN_DEBUG
	l = 0;
#endif
	while(numberOfComponents-- > 0) {
		readIntegers(fileHandle, 1, iA);
		i = iA[0];
		readIntegers(fileHandle, i, iA);
		list2 = constructEmptyList(0, NULL);
		listAppend(list, list2);
		for(j=0; j<i; j++) {
			component = greyEdgeComponents->list[iA[j]-1];
			for(k=0; k<component->length; k++) {
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
	for(i=0; i<list->length; i++) {
		list2= list->list[i];
		for(j=0; j<list2->length; j++) {
			vertex = list2->list[j];
			if(vertex->vertexID == 0) {
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
	free(iA);
	fclose(fileHandle);

	return list;
}

void writeOutCactusGraph(struct CactusGraph *cactusGraph, struct PinchGraph *pinchGraph, struct hashtable *names, FILE *fileHandle) {
	/*
	 * Writes out a graph in 'dot' format, compatible with graphviz.
	 *
	 * The associated function 'writeOutEdgeSegments' gives a way of decoding
	 * the segments associated with the black (segment containing) edges
	 * of the graph.
	 */
	int32_t i, j, k;
	struct CactusVertex *vertex;
	struct CactusEdge *edge;
	struct PinchEdge *pinchEdge;
	struct Segment *segment;
	char *name;

	//Write the preliminaries.
	fprintf(fileHandle, "graph G {\n");
	fprintf(fileHandle, "overlap=false\n");
	//Write the vertices.
	for(i=0; i<cactusGraph->vertices->length;i++) {
		vertex = cactusGraph->vertices->list[i];
#ifdef BEN_DEBUG
		assert(vertex->vertexID == i);
#endif
		fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n", vertex->vertexID, vertex->vertexID);
	}

	fprintf(fileHandle, "edge[color=black,len=2.5,weight=100,dir=forward];\n");
	for(i=0; i<cactusGraph->vertices->length;i++) {
		vertex = cactusGraph->vertices->list[i];
		for(j=0; j<vertex->edges->length; j++) {
			edge = vertex->edges->list[j];
#ifdef BEN_DEBUG
			assert(edge != edge->rEdge);
#endif
			if(edge > edge->rEdge) {
				if(edge->segments->length > 0) {
					for(k=0; k<edge->segments->length; k++) {
						segment = edge->segments->list[k];
						pinchEdge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
						name = (char *)hashtable_search(names, pinchEdge);
						fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n [label=\"" INT_STRING "." INT_STRING ":" INT_STRING ":%s\"];\n",
								edge->from->vertexID, edge->to->vertexID, segment->contig, segment->start, segment->end, name);
					}
				}
				else {
					fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n", edge->from->vertexID, edge->to->vertexID);
				}
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

int32_t computeCactusGraph_excludedEdgesFn(void *o) {
	return TRUE;
}

int32_t computeCactusGraph(struct PinchGraph *pinchGraph, struct CactusGraph **cactusGraph, struct List **threeEdgeConnectedComponents, struct List *extraEdges, char *logLevelString) {
	char *threeEdgeOutputFile;
	char *threeEdgeInputFile;
	static char cA[STRING_ARRAY_SIZE];
	struct PinchVertex *vertex;
	struct List *list;
	int32_t i, j;
	FILE *fileHandle;
	struct List *greyEdgeComponents;

	///////////////////////////////////////////////////////////////////////////
	// Run the three-edge connected component algorithm to identify
	// three edge connected components.
	///////////////////////////////////////////////////////////////////////////

	//Get temp files for the three edge stuff.
	threeEdgeOutputFile = getTempFile();
	threeEdgeInputFile = getTempFile();

	//Write out three-edge graph.
	fileHandle = fopen(threeEdgeOutputFile, "w");

	greyEdgeComponents = getRecursiveComponents(pinchGraph, computeCactusGraph_excludedEdgesFn);

	writeOut3EdgeGraph(pinchGraph, greyEdgeComponents, extraEdges, fileHandle);
	fclose(fileHandle);
	logInfo("Output the 3-edge graph description in tmp file: %s\n", threeEdgeOutputFile);

	//Run three edge connected components.
	sprintf(cA, "cactus_3Edge %s %s %s", logLevelString, threeEdgeOutputFile, threeEdgeInputFile);
	assert(strlen(cA) < STRING_ARRAY_SIZE);
	i = system(cA);
	if(i != 0) {
		logInfo("Tried to run the three edge command, but it went wrong: " INT_STRING "\n", i);
		return i;
	}
	logInfo("Seems to have successfully run the three edge command: " INT_STRING "\n", i);
	//Parse results (the three edge connected components).
	*threeEdgeConnectedComponents = readThreeEdgeComponents(pinchGraph, greyEdgeComponents, threeEdgeInputFile);
	logInfo("Read in the three edge components\n");
	for(i=0; i<(*threeEdgeConnectedComponents)->length; i++) {
		list = (*threeEdgeConnectedComponents)->list[i];
		logDebug("3 edge component : " INT_STRING " ", i);
		for(j=0; j<list->length; j++) {
			vertex = list->list[j];
			logDebug(" vertex, " INT_STRING " ", vertex->vertexID);
		}
		logDebug("\n");
	}
	//Cleanup the three edge input/output file
	removeTempFile(threeEdgeInputFile);
	removeTempFile(threeEdgeOutputFile);
	destructList(greyEdgeComponents);

	///////////////////////////////////////////////////////////////////////////
	// Merge vertices in each three edge connected component to convert the main graph
	// into a 'cactus' graph.
	///////////////////////////////////////////////////////////////////////////

	//Collapse the graph to cactus tree
	*cactusGraph = constructCactusGraph(pinchGraph, extraEdges, *threeEdgeConnectedComponents);
	checkCactusGraph(pinchGraph, *threeEdgeConnectedComponents, *cactusGraph);

	return 0;
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

void circulariseStemsP(struct CactusGraph *cactusGraph,
		struct CactusVertex *vertex, struct CactusVertex *pVertex,
		struct hashtable *stemHash, struct hashtable *seen,  struct List *extraEdges,
		struct List *threeEdgeConnectedComponents) {
	int32_t i, j, k;
	struct CactusEdge *edge;
	struct CactusEdge *edge2;
	struct CactusVertex *vertex2;
	struct List *component;
	struct PinchVertex *pinchVertex1;
	struct PinchVertex *pinchVertex2;
	struct List *list;

	if(hashtable_search(seen, vertex) == NULL) {
		hashtable_insert(seen, vertex, vertex);
		for(i=0; i<vertex->edges->length; i++) {
			edge = vertex->edges->list[i];
			if(hashtable_search(stemHash, edge) != NULL) { //is a stem
				vertex2 = edge->to;
#ifdef BEN_DEBUG
			    assert(vertex != vertex2);
#endif
				k = FALSE;
				for(j=0; j<vertex2->edges->length; j++) {
					edge2 = vertex2->edges->list[j];
					if(edge2->to != vertex) {
						if(hashtable_search(stemHash, edge2) != NULL) {
							k = TRUE;
						}
					}
					else {
#ifdef BEN_DEBUG
						assert(edge2 == edge->rEdge);
#endif
					}
				}
				if(k == FALSE) { //is a stem end
					//add edge between vertex and vertex2
#ifdef BEN_DEBUG
					assert(threeEdgeConnectedComponents->length > vertex2->vertexID);
					assert(threeEdgeConnectedComponents->length > pVertex->vertexID);
#endif
					component = threeEdgeConnectedComponents->list[vertex2->vertexID];
#ifdef BEN_DEBUG
					assert(component->length > 0);
#endif
					pinchVertex1 = component->list[0];
					component = threeEdgeConnectedComponents->list[pVertex->vertexID];
#ifdef BEN_DEBUG
					assert(component->length > 0);
#endif
					pinchVertex2 = component->list[0];
#ifdef BEN_DEBUG
					assert(pinchVertex1 != NULL);
					assert(pinchVertex2 != NULL);
					assert(extraEdges != NULL);
#endif
					//give it a double dose to ensure they are in the same 3-edge component
					list = extraEdges->list[pinchVertex1->vertexID];
#ifdef BEN_DEBUG
					assert(list != NULL);
#endif
					listAppend(list, pinchVertex2);
					listAppend(list, pinchVertex2);

					list = extraEdges->list[pinchVertex2->vertexID];
					listAppend(list, pinchVertex1);
					listAppend(list, pinchVertex1);
				}
				circulariseStemsP(cactusGraph, edge->to, pVertex, stemHash, seen,
						extraEdges, threeEdgeConnectedComponents);
			}
			else { //is not a stem
				circulariseStemsP(cactusGraph, edge->to, edge->to, stemHash, seen,
						extraEdges, threeEdgeConnectedComponents);
			}
		}
	}
}

void circulariseStems(struct CactusGraph *cactusGraph, struct List *extraEdges, struct List *threeEdgeConnectedComponents) {
	struct CactusEdge *edge;
	struct List *biConnectedComponent;
	struct List *biConnectedComponents;
	struct hashtable *stemHash;
	struct hashtable *seen;
	int32_t i;

	logDebug("Circularising the stems\n");

	////////////////////////////////////////////////
	//(1) Get bi-connected components.
	////////////////////////////////////////////////
	biConnectedComponents = computeBiConnectedComponents(cactusGraph);
	logDebug("Constructed the biconnected components for making the net\n");

	////////////////////////////////////////////////
	//(2) Put stems in a hash
	////////////////////////////////////////////////

	stemHash = create_hashtable(cactusGraph->vertices->length*2,
								  hashtable_key, hashtable_equalKey,
								  NULL, NULL);
	for(i=0; i<biConnectedComponents->length; i++) {
		biConnectedComponent = biConnectedComponents->list[i];
		if(biConnectedComponent->length == 1) {
			edge = biConnectedComponent->list[0];
			if(edge->from != edge->to) { //must be a stem as not a self loop
				hashtable_insert(stemHash, edge, edge);
				hashtable_insert(stemHash, edge->rEdge, edge->rEdge);
			}
		}
	}
	logDebug("Put the stems in a hash\n");

	////////////////////////////////////////////////
	//(3) Do DFS on graph
	////////////////////////////////////////////////

	seen = create_hashtable(cactusGraph->vertices->length*2,
							hashtable_key, hashtable_equalKey,
						    NULL, NULL);

	circulariseStemsP(cactusGraph, cactusGraph->vertices->list[0], cactusGraph->vertices->list[0],
					  stemHash, seen, extraEdges, threeEdgeConnectedComponents);
	logDebug("Done the DFS\n");

	////////////////////////////////////////////////
	//(4) Cleanup
	////////////////////////////////////////////////

	hashtable_destroy(stemHash, FALSE, FALSE);
	hashtable_destroy(seen, FALSE, FALSE);
	destructList(biConnectedComponents);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method to break loop discontinuities
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct Segment *edgeDiscontinuityP(struct CactusEdge *edge, struct Segment *segment) {
	int32_t i;
	struct Segment *segment2;
	struct Segment *segment3;

	segment2 = NULL;
	for(i=0; i<edge->segments->length; i++) {
		segment3 = edge->segments->list[i];
		if(segment->contig == segment3->contig) {
			if((segment->start < 1 && segment3->start < 1) || (segment->start > 1 && segment3->start > 1)) {
				if(segment2 == NULL || segment3->start < segment2->start) {
					segment2 = segment3;
				}
			}
		}
	}
	return segment2;
}

int32_t edgeDiscontinuity(struct CactusEdge *edge, struct CactusEdge *edge2) {
	struct Segment *segment;
	struct Segment *segment2;
	struct Segment *segment3;
	int32_t i;

	//for each sequence in edge:
	for(i=0; i<edge->segments->length; i++) {
		segment = edge->segments->list[i];
		//get best match in edge2
		segment2 = edgeDiscontinuityP(edge2, segment);
		//get best match in edge->rEdge
		segment3 = edgeDiscontinuityP(edge->rEdge, segment);
		//if no match in either or match only in edge2 then no discontinuity.
		if(segment3 == NULL) {
			continue;
		}
		//if match only in edge->rEdge then discontinuity
		if(segment2 == NULL) {
			return TRUE;
		}
		//if match in both and:
		//match in edge2 is sooner than match in edge->rEdge then no discontunity
		if(segment2->start < segment3->start) {
			continue;
		}
		//else there is a discontinuity.
		return TRUE;
	}
	return FALSE;
}

void breakLoopDiscontinuities(struct CactusGraph *cactusGraph, struct List *extraEdges,
							  struct List *threeEdgeConnectedComponents) {
	struct CactusEdge *edge;
	struct CactusEdge *edge2;
	struct CactusVertex *vertex;
	struct PinchVertex *pinchVertex1;
	struct PinchVertex *pinchVertex2;

	struct List *biConnectedComponent;
	struct List *list;
	struct List *biConnectedComponents;
	struct List *component;
	int32_t i, j;

	////////////////////////////////////////////////
	//(1) Sort each bi-connected component
	////////////////////////////////////////////////

	biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);
	logDebug("Constructed the sorted biconnected components for making the net\n");

	////////////////////////////////////////////////
	//(2) Identify discontinuities and add extra edges
	////////////////////////////////////////////////

	for(i=0; i<biConnectedComponents->length; i++) {
		biConnectedComponent = biConnectedComponents->list[i];
		edge = biConnectedComponent->list[0];
		vertex = edge->from;
		for(j=0; j+1<biConnectedComponent->length; j++) {
			edge = biConnectedComponent->list[j];
			edge2 = biConnectedComponent->list[j+1];
			if(edgeDiscontinuity(edge, edge2) || edgeDiscontinuity(edge2->rEdge, edge->rEdge)) {
				//add edge between vertex and vertex2
				component = threeEdgeConnectedComponents->list[edge->to->vertexID];
				pinchVertex1 = component->list[0];
				component = threeEdgeConnectedComponents->list[vertex->vertexID];
				pinchVertex2 = component->list[0];

				//give it a double dose to ensure they are in the same 3-edge component
				list = extraEdges->list[pinchVertex1->vertexID];
				listAppend(list, pinchVertex2);
				listAppend(list, pinchVertex2);

				list = extraEdges->list[pinchVertex2->vertexID];
				listAppend(list, pinchVertex1);
				listAppend(list, pinchVertex1);
			}
		}
	}
	logDebug("Computed the discontinuities\n");

	////////////////////////////////////////////////
	//(6)Cleanup/return net
	////////////////////////////////////////////////
	destructList(biConnectedComponents);
	free(sortBiConnectedComponents_vertexOrdering);

	logDebug("Cleaned up having calculated the loop discontinuities.\n");
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Choosing which atoms in the chains to keep.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

//first get tree covering score for each atom -
//drop all atoms with score less than X.
//accept chains whose remaining element's combined length is greater than a set length.

float treeCoverage(struct CactusEdge *cactusEdge, struct List *nodeSets,
		struct List *contigIndexToContigStrings, struct PinchGraph *pinchGraph) {
	/*
	 * Returns the proportion of the tree covered by the atom.
	 */
	int32_t i, j, k;
	struct ContigEventSet *contigEventSet;
	struct Segment *segment;
	float treeCoverage = 0.0;
	float totalTreeBranchLength = 0.0;
	char *cA;

#ifdef BEN_DEBUG
	assert(!isAStubOrCapCactusEdge(cactusEdge, pinchGraph));
#endif
	for(i=0; i<nodeSets->length; i++) {
		contigEventSet = nodeSets->list[i];
		totalTreeBranchLength += contigEventSet->branchLength;
		k = 0;
		for(j=0; j<cactusEdge->segments->length; j++) {
			segment = cactusEdge->segments->list[j];
#ifdef BEN_DEBUG
			assert(segment->contig < contigIndexToContigStrings->length);
#endif
			cA = contigIndexToContigStrings->list[segment->contig];
			if(hashtable_search(contigEventSet->contigs, cA) != NULL) {
				k += 1;
			}
		}
		if(k == cactusEdge->segments->length) {
			treeCoverage = 0.0; //reset as all segments must be contained within children of this node
		}
		else if(k > 0) {
			treeCoverage += contigEventSet->branchLength;
		}
	}
#ifdef BEN_DEBUG
	assert(cactusEdge->segments->length > 0);
	assert(treeCoverage >= 0);
	assert(treeCoverage <= totalTreeBranchLength);
#endif
	return treeCoverage/totalTreeBranchLength;
}

float atomScore(struct CactusEdge *cactusEdge, struct List *nodeSets,
				struct List *contigIndexToContigStrings, struct PinchGraph *pinchGraph) {
	struct Segment *segment = cactusEdge->segments->list[0];
	return (segment->end - segment->start + 1) * treeCoverage(cactusEdge, nodeSets, contigIndexToContigStrings, pinchGraph);
}

int32_t chainScore(struct List *biConnectedComponent, struct List *nodeSets,
		struct List *contigIndexToContigStrings, struct PinchGraph *pinchGraph) {
	/*
	 * Gets the 'length' of the chain, including only those atoms in the chain that are
	 * true according to the 'includeAtom' function.
	 *
	 * Length is total length of atoms (all sequences in an atom have the same length).
	 */
	int32_t i, totalScore;
	struct CactusEdge *cactusEdge;

	totalScore = 0;
	for(i=0; i<biConnectedComponent->length; i++) {
		cactusEdge = biConnectedComponent->list[i];
		if(!isAStubOrCapCactusEdge(cactusEdge, pinchGraph)) {
			totalScore += atomScore(cactusEdge, nodeSets, contigIndexToContigStrings, pinchGraph);
		}
	}
	return totalScore;
}

int
compareFloats (const float *a, const float *b)
{
  return (int) (*a - *b);
}

int32_t chainLength(struct List *biConnectedComponent, int32_t includeStubs, struct PinchGraph *pinchGraph) {
	/*
	 * Get the number of links in the chain.
	 */
	int32_t i, j;
	struct CactusEdge *cactusEdge;
	i = 0;
	for(j=0; j<biConnectedComponent->length; j++) {
		cactusEdge = biConnectedComponent->list[j];
		if(includeStubs || !isAStubOrCapCactusEdge(cactusEdge, pinchGraph)) {
			i++;
		}
	}
	return i;
}

int32_t chainBaseLength(struct List *biConnectedComponent, struct PinchGraph *pinchGraph) {
	/*
	 * Get the number of links in the chain.
	 */
	int32_t i, j;
	struct CactusEdge *cactusEdge;
	struct Segment *segment;

	i = 0;
	for(j=0; j<biConnectedComponent->length; j++) {
		cactusEdge = biConnectedComponent->list[j];
		if(!isAStubOrCapCactusEdge(cactusEdge, pinchGraph)) {
			segment = cactusEdge->segments->list[0];
			i += segment->end - segment->start + 1;
		}
	}
	return i;
}

void filterAtomsByTreeCoverageAndLength(struct List *biConnectedComponents,
		struct List *chosenAtoms,
		struct List *nodeSets,
		float proportionToKeep, /*Proportion of all atoms to select to keep*/
		float discardRatio, /*The proportion of an atom's chain's average atom score required to be score to be considered */
		float minimumTreeCoverage, /*Minimum tree coverage to be included */
		int32_t minimumChainLength, /* Minimum chain length to be included */
		struct PinchGraph *pinchGraph,
		struct List *contigIndexToContigStrings) {
	/*
	 * Filters atoms in chains by length and tree coverage score.
	 *
	 * Returns a list of all accepted atoms.
	 */
	int32_t i, j;
	struct CactusEdge *cactusEdge;
	float *fA;
	float minScore;
	float minAtomScore;
	float f;
	struct List *biConnectedComponent;

#ifdef BEN_DEBUG
	assert(proportionToKeep <= 1.0);
	assert(proportionToKeep >= 0.0);
	assert(discardRatio >= 0.0);
#endif

	///////////////
	//Calculate the minimum score to be kept.
	///////////////
	fA = mallocLocal(sizeof(float)*biConnectedComponents->length);
	for(i=0; i<biConnectedComponents->length; i++) {
		biConnectedComponent = biConnectedComponents->list[i];
		fA[i] = chainScore(biConnectedComponent, nodeSets, contigIndexToContigStrings, pinchGraph);
	}
	qsort(fA, biConnectedComponents->length, sizeof(float),
			(int (*)(const void *v, const void *))compareFloats);
	i = (1.0 - proportionToKeep)*biConnectedComponent->length;
	minScore = i < biConnectedComponents->length ? fA[i] : 10000000000000.0;
	logInfo("The minimum chain score to be kept is %f and %f \n", minScore, proportionToKeep);
	free(fA);

	///////////////
	//Gets those atoms whose chain meet the minimum score and have at
	///////////////

	for(i=0; i<biConnectedComponents->length; i++) {
		biConnectedComponent = biConnectedComponents->list[i];
		f = chainScore(biConnectedComponent, nodeSets, contigIndexToContigStrings, pinchGraph);
		if(f >= minScore && chainBaseLength(biConnectedComponent, pinchGraph) >= minimumChainLength) {
			minAtomScore = (minScore / chainLength(biConnectedComponent, FALSE, pinchGraph)) * discardRatio;
			for(j=0; j<biConnectedComponent->length; j++) {
				cactusEdge = biConnectedComponent->list[j];
				if(!isAStubOrCapCactusEdge(cactusEdge, pinchGraph)) {
					if(treeCoverage(cactusEdge, nodeSets, contigIndexToContigStrings, pinchGraph) >= minimumTreeCoverage &&
						atomScore(cactusEdge, nodeSets, contigIndexToContigStrings, pinchGraph) >= minAtomScore) {
						listAppend(chosenAtoms, cactusEdge);
					}
				}
			}
		}
	}
}

void logTheChosenAtomSubset(struct List *biConnectedComponents, struct List *chosenAtoms, struct PinchGraph *pinchGraph,
		struct List *nodeSets, struct List *contigIndexToContigStrings) {
	/*
	 * Produces logging information about the chosen atoms.
	 */
	int32_t i, j;
	struct List *biConnectedComponent;
	struct CactusEdge *cactusEdge;
	struct Segment *segment;
	float totalAtomScore = 0.0;
	float totalAtomLength = 0.0;
	float totalBaseLengthOfAllAtoms = 0.0;
	float totalNumberOfAllAtoms = 0.0;
	float totalNumberOfStubAtoms = 0.0;
	float averageSegmentNumber = 0.0;
	float averageSegmentNumberOfAllAtoms = 0.0;
	for(i=0; i<biConnectedComponents->length; i++) {
		biConnectedComponent = biConnectedComponents->list[i];
		totalBaseLengthOfAllAtoms += chainBaseLength(biConnectedComponent, pinchGraph);
		totalNumberOfAllAtoms += chainLength(biConnectedComponent, FALSE, pinchGraph);
		totalNumberOfStubAtoms += chainLength(biConnectedComponent, TRUE, pinchGraph) - chainLength(biConnectedComponent, FALSE, pinchGraph);
		for(j=0; j<biConnectedComponent->length; j++) {
			cactusEdge = biConnectedComponent->list[j];
			if(!isAStubOrCapCactusEdge(cactusEdge, pinchGraph)) {
				averageSegmentNumberOfAllAtoms += cactusEdge->segments->length;
			}
		}
	}
	j = 0;
	for(i=0; i<chosenAtoms->length; i++) {
		cactusEdge = chosenAtoms->list[i];
		if(!isAStubOrCapCactusEdge(cactusEdge, pinchGraph)) {
			totalAtomScore += treeCoverage(cactusEdge, nodeSets, contigIndexToContigStrings, pinchGraph);
			segment = cactusEdge->segments->list[0];
			totalAtomLength += segment->end - segment->start + 1;
			averageSegmentNumber += cactusEdge->segments->length;
			j++;
		}
	}
	logInfo("Chosen atom subset composed of %i atoms, of average length %f and average tree coverage %f, total base length of all atoms: %f, total number of all atoms %f, average length of all atoms: %f, total number of stub atoms: %f, average segment number of chosen atoms: %f, average segment number of all atoms: %f\n",
				chosenAtoms->length, totalAtomLength/j, totalAtomScore/j, totalBaseLengthOfAllAtoms, totalNumberOfAllAtoms, totalBaseLengthOfAllAtoms/totalNumberOfAllAtoms, totalNumberOfStubAtoms, averageSegmentNumber/j, averageSegmentNumberOfAllAtoms/totalNumberOfAllAtoms);
}
