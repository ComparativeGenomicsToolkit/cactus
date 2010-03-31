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
#include "net.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Chains
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct Chain *constructChain(struct List *segments, struct Chain *pLink) {
	struct Chain *chain;
	int32_t i;
	struct Segment *segment;

	chain = mallocLocal(sizeof(struct Chain));

	chain->subChains = constructEmptyList(0, (void(*)(void *))destructChain);
	chain->adjacencyComponents = constructEmptyList(0, (void (*)(void *))destructList);
	chain->ends = constructEmptyList(0, NULL);
	chain->segments = constructEmptyList(0, NULL);
	chain->pLink = pLink;
	chain->nLink = NULL;
	chain->score = 0.0;

	if(pLink != NULL) {
		pLink->nLink = chain;
	}

	for(i=0; i<segments->length; i++) {
		segment = segments->list[i];
		listAppend(chain->segments, segment);
	}
	return chain;
}

void destructChain(struct Chain *chain) {
	destructList(chain->subChains);
	destructList(chain->adjacencyComponents);
	destructList(chain->ends);
	destructList(chain->segments);
	if(chain->nLink != NULL) {
		destructChain(chain->nLink);
	}
	free(chain);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods on chains
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int sortChains_P(const struct Chain **c1, const struct Chain **c2) {
	//Sort in descending order of score
	return ((*c1)->score > (*c2)->score) ? -1 : ((*c1)->score < (*c2)->score ? 1 : 0);
}

void sortChains(struct List *listOfChains) {
	/*
	 * Sorts a list of alignments by there script (in descending order of score).
	 */
	qsort(listOfChains->list, listOfChains->length, sizeof(void *),
			(int (*)(const void *v, const void *))sortChains_P);
}

void flattenChainList(struct List *chainList, struct List *flatList) {
	int32_t i;
	struct Chain *chain;

	for(i=0; i<chainList->length; i++) {
		chain = chainList->list[i];
		listAppend(flatList, chain);
		while(chain != NULL) {
			flattenChainList(chain->subChains, flatList);
			chain = chain->nLink;
		}
	}
}

int32_t isAStubOrCapChain(struct Chain *chain, struct PinchGraph *pinchGraph) {
	/*
	 * Returns true if this element in the chain represents either a stub or cap.
	 */
	struct Segment *segment;
	struct PinchEdge *pinchEdge;
#ifdef BEN_DEBUG
	assert(chain->segments->length > 0);
#endif
	segment = (struct Segment *)chain->segments->list[0];
	pinchEdge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
	return isAStubOrCap(pinchEdge);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Nets
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct Net *constructNet(struct CactusGraph *cactusGraph) {
	struct Net *net;
	struct Chain *chain;

	struct CactusEdge *edge;
	struct List *chainsAdjacencyList;
	struct List *chains;
	struct List *edges;
	struct List *biConnectedComponent;
	struct List *list;
	struct List *biConnectedComponents;
#ifdef BEN_DEBUG
	struct CactusEdge *edge2;
#endif

	int32_t i, j;

	logDebug("Building the net\n");

	////////////////////////////////////////////////
	//(1) Get sorted bi-connected components.
	////////////////////////////////////////////////

	biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);

	////////////////////////////////////////////////
	//(2) Construct chain for each cycle.
	////////////////////////////////////////////////

	chainsAdjacencyList = constructEmptyList(cactusGraph->vertices->length, NULL);
	for(i=0; i<cactusGraph->vertices->length; i++) {
		chainsAdjacencyList->list[i] = constructEmptyList(0, (void (*)(void *))destructChain);
	}

	chains = constructEmptyList(0, NULL);
	edges = constructEmptyList(0, (void (*)(void *))destructList);

	for(i=0; i<biConnectedComponents->length; i++) {
		biConnectedComponent = biConnectedComponents->list[i];
		edge = biConnectedComponent->list[0];
		chain = constructChain(edge->segments, NULL);
		qsort(chain->segments->list, chain->segments->length, sizeof(void *),
				(int (*)(const void *, const void*))segmentComparatorPointers);
		list = chainsAdjacencyList->list[edge->from->vertexID];
		listAppend(list, chain);
		listAppend(chains, chain);
		list = constructEmptyList(0, NULL);
		listAppend(edges, list);
		listAppend(list, edge);

		for(j=1; j<biConnectedComponent->length; j++) {
			edge = biConnectedComponent->list[j];
			listAppend(list, edge);
			chain = constructChain(edge->segments, chain);
			qsort(chain->segments->list, chain->segments->length, sizeof(void *),
				  (int (*)(const void *, const void*))segmentComparatorPointers);
		}
	}
	logDebug("Constructed chain for each cycle.\n");

	////////////////////////////////////////////////
	//(3) Link sub-chains to parent chains.
	////////////////////////////////////////////////

#ifdef BEN_DEBUG
	assert(edges->length == chains->length);
#endif

	for(i=0; i<chains->length; i++) {

		chain = chains->list[i];
		list = edges->list[i];

		for(j=1; j<list->length; j++) {
#ifdef BEN_DEBUG
			assert(chain != NULL);
#endif
			edge = list->list[j];
#ifdef BEN_DEBUG
			assert(edge->from->vertexID != 0);
			edge2 = list->list[j-1];
			assert(edge2->to == edge->from);
#endif
			destructList(chain->subChains);
			chain->subChains = chainsAdjacencyList->list[edge->from->vertexID];
#ifdef BEN_DEBUG
			assert(chain->subChains != NULL);
			chainsAdjacencyList->list[edge->from->vertexID] = NULL;
#endif
			chain = chain->nLink;
		}
#ifdef BEN_DEBUG
		assert(chain != NULL);
		assert(chain->nLink == NULL);
#endif
	}

	net = mallocLocal(sizeof(struct Net));
	net->adjacencyComponents = constructEmptyList(0, (void (*)(void *))destructList);
	net->ends = constructEmptyList(0, NULL);
	net->chains = chainsAdjacencyList->list[0];
#ifdef BEN_DEBUG
	assert(net->chains != NULL);
	chainsAdjacencyList->list[0] = NULL;
	for(i=0; i<chainsAdjacencyList->length; i++) {
		assert(chainsAdjacencyList->list[i] == NULL);
	}
#endif
	logDebug("Linked sub-chains to parent chains.\n");

	////////////////////////////////////////////////
	//(6)Cleanup/return net
	////////////////////////////////////////////////
	destructList(chainsAdjacencyList);
	destructList(chains);
	destructList(edges);
	destructList(biConnectedComponents);

	//check the constructed net
	checkNet(cactusGraph, net);

	logDebug("Cleaned up and am returning the finished net.\n");

	return net;
}

void destructNet(struct Net *net) {
	destructList(net->adjacencyComponents);
	destructList(net->chains);
	destructList(net->ends);
	free(net);
}

void checkNet(struct CactusGraph *graph, struct Net *net) {
	/*
	 * Check that the net is okay with respect to the pinch graph.
	 *
	 * Currently checks that every base is covered (every aligned base is included).
	 */
#ifdef BEN_DEBUG
	int32_t i, j, k;
	struct Chain *chain;
	struct Segment *segment;
	struct List *flatList;
	struct hashtable *segmentsHash;
	struct CactusVertex *vertex;
	struct CactusEdge *edge;

	flatList = constructEmptyList(0, NULL);
	flattenChainList(net->chains, flatList);

	//check every base is covered (every aligned base is included).
	segmentsHash = create_hashtable(graph->vertices->length*2,
					hashtable_key, hashtable_equalKey,
					NULL, NULL);
	for(i=0; i<flatList->length; i++) {
		chain = flatList->list[i];
		while(chain != NULL) {
			for(j=0; j<chain->segments->length; j++) {
				segment = chain->segments->list[j];

				//check not already in
				assert(hashtable_search(segmentsHash, segment) == NULL);
				assert(hashtable_search(segmentsHash, segment->rSegment) == NULL);

				//now add
				hashtable_insert(segmentsHash, segment, segment);
				hashtable_insert(segmentsHash, segment->rSegment, segment->rSegment);
			}
			chain = chain->nLink;
		}
	}

	for(i=0; i<graph->vertices->length; i++) {
		vertex = graph->vertices->list[i];
		for(j=0; j<vertex->edges->length; j++) {
			edge = vertex->edges->list[j];
			for(k=0; k<edge->segments->length; k++) {
				assert(hashtable_search(segmentsHash, edge->segments->list[k]) != NULL);
			}
		}
	}
	hashtable_destroy(segmentsHash, FALSE, FALSE);

	destructList(flatList);

#endif
}

void writeOutNet(struct Net *net, struct hashtable *names, FILE *fileHandle) {
	/*
	 * Writes out the chains in the net in 'dot' format, compatible with graphviz.
	 */
	int32_t i, j, k, l, m;
	struct Chain *chain;
	struct List *list;
	struct List *adjacencyComponent;
	struct hashtable *chainHash;

	logDebug("Writing a net\n");

	//Write the preliminaries.
	fprintf(fileHandle, "graph G {\n");
	fprintf(fileHandle, "overlap=false\n");

	list = constructEmptyList(0, NULL);
	flattenChainList(net->chains, list);


	chainHash = create_hashtable(list->length*2,
							     hashtable_key, hashtable_equalKey,
								 NULL, (void (*)(void *))destructInt);

	//Write the sink vertex
	fprintf(fileHandle, "node[width=0.3,height=0.3,shape=diamond,fontsize=14];\n");
	fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n", 0, 0);

	//Write the chain vertices.
	j = 1;
	for(i=0; i<list->length; i++) {
		chain = list->list[i];
		fprintf(fileHandle, "node[width=0.3,height=0.3,shape=circle,fontsize=14];\n");
		fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n", j, j);
		hashtable_insert(chainHash, chain, constructInt(j++));
	}

	//Write the organising link vertices.
	for(i=0; i<list->length; i++) {
		chain = list->list[i];
		while(chain->nLink != NULL) {
			hashtable_insert(chainHash, chain->subChains, constructInt(j));
			fprintf(fileHandle, "node[width=0.3,height=0.3,shape=box,fontsize=14];\n");
			fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n", j, j);
			chain = chain->nLink;
			j++;
		}
	}

	//Write the adjacency component vertices.
	for(i=0; i<list->length; i++) {
		chain = list->list[i];
		while(chain->nLink != NULL) {
			for(k=0; k<chain->adjacencyComponents->length; k++) {
				adjacencyComponent = chain->adjacencyComponents->list[k];
				hashtable_insert(chainHash, adjacencyComponent, constructInt(j));
				fprintf(fileHandle, "node[width=0.3,height=0.3,shape=diamond,fontsize=14];\n");
				fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n", j, j);
				j++;
			}
			chain = chain->nLink;
		}
	}

	//Write the adjacency components attached to the sink
	for(k=0; k<net->adjacencyComponents->length; k++) {
		adjacencyComponent = net->adjacencyComponents->list[k];
		hashtable_insert(chainHash, adjacencyComponent, constructInt(j));
		fprintf(fileHandle, "node[width=0.3,height=0.3,shape=diamond,fontsize=14];\n");
		fprintf(fileHandle, "n" INT_STRING "n [label=\"" INT_STRING "\"];\n", j, j);
		j++;
	}

	logDebug("Written the vertices of the net\n");

	//Format of the edges.
	fprintf(fileHandle, "edge[color=black,len=2.5,weight=100,dir=forward];\n");

	//Write the sink edges.

	//edge from sink node to child node
	for(i=0; i<net->chains->length; i++) {
		j = *((int32_t *)hashtable_search(chainHash, net->chains->list[i]));
		fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n", 0, j);
	}

	//edge from sink node to adjacency components
	for(i=0; i<net->adjacencyComponents->length; i++) {
		j = *((int32_t *)hashtable_search(chainHash, net->adjacencyComponents->list[i]));
		fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n", 0, j);
	}

	logDebug("Written the edges from the sink node\n");

	//Write the edges excluding the sink edges
	for(i=0; i<list->length; i++) {
		chain = list->list[i];
		j = *((int32_t *)hashtable_search(chainHash, chain));
		while(chain->nLink != NULL) {
			k = *((int32_t *)hashtable_search(chainHash, chain->subChains));
			//edge from parent chain node to organising node
			fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n", j, k);

			//edge from organising node to child chain.
			for(l=0; l<chain->subChains->length; l++) {
				m = *((int32_t *)hashtable_search(chainHash, chain->subChains->list[l]));
				fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n", k, m);
			}

			//edge from organising node to adjacency component
			for(l=0; l<chain->adjacencyComponents->length; l++) {
				m = *((int32_t *)hashtable_search(chainHash, chain->adjacencyComponents->list[l]));
				fprintf(fileHandle, "n" INT_STRING "n -- n" INT_STRING "n;\n", k, m);
			}

			chain = chain->nLink;
		}
	}

	logDebug("Written the edges of the net\n");

	fprintf(fileHandle, "}\n");

	destructList(list);
	hashtable_destroy(chainHash, TRUE, FALSE);

	logDebug("Written the net\n");
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods to manipulate a net.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void pruneNetP_destructChain(struct Chain *chain) { //safely deletes the chain
	chain->nLink = NULL;
	chain->pLink = NULL;
	chain->subChains->length = 0;
	destructChain(chain);
}

struct List *pruneNetP(struct Chain *chain, struct hashtable *chosenChains) {
	int32_t i, j;
	struct List *list;
	struct List *list2;
	struct List *list3;
	struct Chain *chainToDelete;

	list = constructEmptyList(0, NULL);
	chainToDelete = NULL;
	while(TRUE) {
		list2 = constructEmptyList(0, NULL);
#ifdef BEN_DEBUG
		if(chain->nLink == NULL) {
			assert(chain->subChains->length == 0);
		}
#endif
		for(i=0; i<chain->subChains->length; i++) {
			list3 = pruneNetP(chain->subChains->list[i], chosenChains);
			for(j=0; j<list3->length; j++) {
				listAppend(list2, list3->list[j]);
			}
			destructList(list3);
		}

		if(hashtable_search(chosenChains, chain) == NULL) { //must be removed
			chainToDelete = chain;
			if(chain->nLink != NULL) {
				chain->nLink->pLink = chain->pLink;
			}
			if(chain->pLink != NULL) {
				chain->pLink->nLink = chain->nLink;
				for(j=0; j<list2->length; j++) {
					listAppend(chain->pLink->subChains, list2->list[j]);
				}
				if(chain->nLink == NULL) { //have to remove when last link in chain is removed
					for(j=0; j<chain->pLink->subChains->length; j++) {
						listAppend(list, chain->pLink->subChains->list[j]);
					}
					chain->pLink->subChains->length = 0;
				}
			}
			else {
				for(j=0; j<list2->length; j++) {
					listAppend(list, list2->list[j]);
				}
			}
		}
		else {
			chain->subChains->length = 0;
			for(j=0; j<list2->length; j++) {
				listAppend(chain->subChains, list2->list[j]);
			}
		}
		destructList(list2);
		if(chain->nLink == NULL) {
			while(chain->pLink != NULL) {
				chain = chain->pLink;
			}
			if(hashtable_search(chosenChains, chain) != NULL) {
				listAppend(list, chain);
			}
			if(chainToDelete != NULL) { //now we delete it
				pruneNetP_destructChain(chainToDelete);
				chainToDelete = NULL;
			}
			break;
		}
		chain = chain->nLink;
		if(chainToDelete != NULL) { //now we delete it
			pruneNetP_destructChain(chainToDelete);
			chainToDelete = NULL;
		}
	}
	return list;
}

void pruneNet(struct Net *net, struct List *chosenAtoms) {
	/*
	 * Method removes all atoms not in the chosen atom list from the net.
	 */
	int32_t i, j;
	struct Chain *chain;
	struct hashtable *hash;
	struct List *list;
	struct List *list2;
#ifdef BEN_DEBUG
	int32_t k;
#endif

	hash = create_hashtable(chosenAtoms->length*2,
								hashtable_key, hashtable_equalKey,
								NULL, NULL);
	for(i=0; i<chosenAtoms->length; i++) {
		chain = chosenAtoms->list[i];
		hashtable_insert(hash, chain, chain);
	}

	list = constructEmptyList(0, NULL);
	for(i=0; i<net->chains->length; i++) {
		list2 = pruneNetP(net->chains->list[i], hash);
		for(j=0; j<list2->length; j++) {
			listAppend(list, list2->list[j]);
		}
		destructList(list2);
	}
	net->chains->length = 0;
	for(j=0; j<list->length; j++) {
		listAppend(net->chains, list->list[j]);
	}


#ifdef BEN_DEBUG
	assert(list->destructElement == NULL); //defensive
	list->length = 0;
	flattenChainList(net->chains, list);

	j = 0;
	for(i=0; i<list->length; i++) {
		chain = list->list[i];
		while(chain != NULL) {
			assert(hashtable_search(hash, chain) != NULL);
			j++;
			if(chain->nLink == NULL) {
				assert(chain->subChains->length == 0);
			}
			k = chain->subChains->length;
			listRemoveDuplicates(chain->subChains);
			assert(chain->subChains->length == k);
			chain = chain->nLink;
		}
	}
	assert(j == chosenAtoms->length);
#endif

	destructList(list);
	hashtable_destroy(hash, FALSE, FALSE);
}

struct List *identifyPseudoAdjacencyComponents(struct List *chosenAtoms, struct PinchGraph *pinchGraph) {
	/*
	 * Method identifies pseudo adjacency components.
	 */
	int32_t i, j;
	struct List *list;
	struct List *list2;
	struct Chain *chain;
	struct Segment *segment;
	struct PinchEdge *pinchEdge;
	struct PinchVertex *pinchVertex;
	struct List *adjacencyComponents;
	struct hashtable *capsHash;

	//////////////////
	//get the pseudo adjacency components.
	//////////////////

	list = constructEmptyList(0, NULL);
	for(i=0; i<chosenAtoms->length; i++) { //build the excluded edges hash.
		chain = chosenAtoms->list[i];
#ifdef BEN_DEBUG
		assert(chain->segments->length > 0);
#endif
		segment = chain->segments->list[0];
		pinchEdge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
#ifdef BEN_DEBUG
		assert(pinchEdge != NULL);
#endif
		listAppend(list, pinchEdge);
	}
	adjacencyComponents = getRecursiveComponents2(pinchGraph, list);
	destructList(list);

	/////////////////
	//build the caps hash
	////////////////////
	capsHash = create_hashtable(chosenAtoms->length*2,
								hashtable_key, hashtable_equalKey,
								NULL, NULL);
	hashtable_insert(capsHash, pinchGraph->vertices->list[0], pinchGraph->vertices->list[0]);

	for(i=0; i<chosenAtoms->length; i++) {
		chain = chosenAtoms->list[i];
#ifdef BEN_DEBUG
		assert(chain->segments->length > 0);
#endif
		segment = chain->segments->list[0];
		pinchEdge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
		hashtable_insert(capsHash, pinchEdge->from, pinchEdge->from);
		hashtable_insert(capsHash, pinchEdge->to, pinchEdge->to);
	}

#ifdef BEN_DEBUG
	for(i=0; i<pinchGraph->vertices->length; i++) { //check chosen atoms contains all the caps.
		pinchVertex = pinchGraph->vertices->list[i];
		void *blackEdgeIterator = getBlackEdgeIterator(pinchVertex);
		while((pinchEdge = getNextBlackEdge(pinchVertex, blackEdgeIterator)) != NULL) {
			if(isAStubOrCap(pinchEdge)) {
				assert(hashtable_search(capsHash, pinchEdge->from) != NULL);
				assert(hashtable_search(capsHash, pinchEdge->to) != NULL);
			}
		}
		destructBlackEdgeIterator(blackEdgeIterator);
	}
#endif

	//////////////////////////////
	//get rid of the vertices in the adjacency components not linked to chosen atoms.
	//////////////////////////////
	for(i=0; i<adjacencyComponents->length; i++) {
		list = adjacencyComponents->list[i];
		list2 = constructEmptyList(0, NULL);
		for(j=0; j<list->length; j++) {
			pinchVertex = list->list[j];
			if(hashtable_search(capsHash, pinchVertex) != NULL) {
				listAppend(list2, pinchVertex);
			}
		}
#ifdef BEN_DEBUG
		assert(list2->length > 0);
#endif
		destructList(adjacencyComponents->list[i]);
		adjacencyComponents->list[i] = list2;
	}
	hashtable_destroy(capsHash, FALSE, FALSE);

	return adjacencyComponents;
}

struct List *addAdjacencyComponentsP(struct List *subChains, struct PinchGraph *pinchGraph) {
	struct Chain *chain;
	struct Segment *segment;
	struct PinchEdge *edge;
	struct List *list;
	int32_t i;

	list = constructEmptyList(0, NULL);
	for(i=0; i<subChains->length; i++) {
		chain = subChains->list[i];
		//start of chain
		segment = chain->segments->list[0];
		edge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
		listAppend(list, edge->from);
		while(chain->nLink != NULL) {
			chain = chain->nLink;
		}
		//end of chain
		segment = chain->segments->list[0];
		edge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
		listAppend(list, edge->to);
	}

	return list;
}

void addAdjacencyComponents(struct Net *net, struct PinchGraph *pinchGraph, struct List *adjacencyComponents) {
	/*
	 * Adds the adjacency components to the net.
	 */
	int32_t i, j;
	struct List *list;
	struct List *list2;
	struct List *adjacencyComponent;
	struct Chain *chain;
	struct hashtable *hash;
	struct PinchEdge *edge;
	struct PinchVertex *vertex;
	struct Segment *segment;
#ifdef BEN_DEBUG
	int32_t k, l;
	struct PinchVertex *vertex2;
#endif

	list = constructEmptyList(0, NULL);
	flattenChainList(net->chains, list);
	hash = create_hashtable(list->length*2,
							hashtable_key, hashtable_equalKey,
							NULL, NULL);

	for(i=0; i<adjacencyComponents->length; i++) {
		adjacencyComponent = listCopy(adjacencyComponents->list[i]);
		for(j=0; j<adjacencyComponent->length; j++) {
			vertex = adjacencyComponent->list[j];
			hashtable_insert(hash, vertex, adjacencyComponent);
		}
	}

#ifdef BEN_DEBUG
	k=0;
#endif
	for(i=0; i<list->length; i++) {
		chain = list->list[i];
		while(chain->nLink != NULL) {
			list2 = addAdjacencyComponentsP(chain->subChains, pinchGraph);
			//left side of the chain
			segment = chain->segments->list[0];
			edge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
			listAppend(list2, edge->to);
			//right side of the chain
			segment = chain->nLink->segments->list[0];
			edge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
			listAppend(list2, edge->from);
			//adjacency component
			for(j=0; j<list2->length; j++) {
				vertex = list2->list[j];
				adjacencyComponent = hashtable_search(hash, vertex);
#ifdef BEN_DEBUG
				assert(adjacencyComponent != NULL);
#endif
				if(listContains(chain->adjacencyComponents, adjacencyComponent) == FALSE) {
					listAppend(chain->adjacencyComponents, adjacencyComponent);
#ifdef BEN_DEBUG
					k++;
					for(l=0; l<adjacencyComponent->length; l++) {
						vertex2 = adjacencyComponent->list[l];
						assert(listContains(list2, vertex2) == TRUE);
					}
#endif
				}
			}
			destructList(list2);
			chain = chain->nLink;
		}
#ifdef BEN_DEBUG
		assert(chain->subChains->length == 0);
#endif
	}
	destructList(list);

	list2 = addAdjacencyComponentsP(net->chains, pinchGraph);
	listAppend(list2, pinchGraph->vertices->list[0]);
	for(j=0; j<list2->length; j++) {
		vertex = list2->list[j];
		adjacencyComponent = hashtable_search(hash, vertex);
#ifdef BEN_DEBUG
		assert(adjacencyComponent != NULL);
#endif
		if(listContains(net->adjacencyComponents, adjacencyComponent) == FALSE) {
			listAppend(net->adjacencyComponents, adjacencyComponent);
#ifdef BEN_DEBUG
			k++;
			for(l=0; l<adjacencyComponent->length; l++) {
				vertex2 = adjacencyComponent->list[l];
				assert(listContains(list2, vertex2) == TRUE);
			}
#endif
		}
	}
	destructList(list2);

#ifdef BEN_DEBUG
	assert(k == adjacencyComponents->length);
#endif

	hashtable_destroy(hash, FALSE, FALSE);
}

int32_t isAStubAdjacencyComponent(struct List *adjacencyComponent) {
	int32_t i;
	struct PinchVertex *vertex;

	for(i=0; i<adjacencyComponent->length; i++) {
		vertex = adjacencyComponent->list[i];
		if(vertex->vertexID == 0) {
#ifdef BEN_DEBUG
			int32_t j;
			for(j=0; j<adjacencyComponent->length; j++) {
				vertex = adjacencyComponent->list[j];
				if(vertex->vertexID != 0) {
					assert(lengthBlackEdges(vertex) > 0);
					assert(isAStubOrCap(getFirstBlackEdge(vertex)));
				}
			}
#endif
			return TRUE;
		}
		assert(lengthBlackEdges(vertex) > 0);
		if(isAStubOrCap(getFirstBlackEdge(vertex)) && lengthGreyEdges(vertex) == 0) {
#ifdef BEN_DEBUG
			assert(adjacencyComponent->length == 1);
#endif
			return TRUE;
		}
	}

	return FALSE;
}

struct List *processAdjacencyComponents(struct List *adjacencyComponents) {
	struct List *list;
	int32_t i;

	list = constructEmptyList(0, adjacencyComponents->destructElement);
	for(i=0; i<adjacencyComponents->length; i++) {
		if(isAStubAdjacencyComponent(adjacencyComponents->list[i]) == FALSE) {
			listAppend(list, adjacencyComponents->list[i]);
		}
		else {
			destructList(adjacencyComponents->list[i]);
		}
	}
	adjacencyComponents->destructElement = NULL;
	destructList(adjacencyComponents);
	return list;
}

void removeStubAdjacencyComponents(struct Net *net) {
	/*
	 * Removes adjacency components from the net which contain the invisible ends of stubs/caps.
	 */
	struct List *list;
	int32_t i;
	struct Chain *chain;
#ifdef BEN_DEBUG
	int32_t j;
#endif

	list = constructEmptyList(0, NULL);
	flattenChainList(net->chains, list);
	net->adjacencyComponents = processAdjacencyComponents(net->adjacencyComponents);

#ifdef BEN_DEBUG
	i = net->adjacencyComponents->length;
	listRemoveDuplicates(net->adjacencyComponents);
	assert(i == net->adjacencyComponents->length);
#endif

	for(i=0; i<list->length; i++) {
		chain = list->list[i];
		while(chain->nLink != NULL) {
			chain->adjacencyComponents = processAdjacencyComponents(chain->adjacencyComponents);

#ifdef BEN_DEBUG
			j = chain->adjacencyComponents->length;
			listRemoveDuplicates(chain->adjacencyComponents);
			assert(j == chain->adjacencyComponents->length);
#endif

			chain = chain->nLink;
		}
	}
	destructList(list);
}

void addEndsP(struct List *adjacencyComponents, struct List *chains, struct List *ends, struct List *parentEnds) {
	int32_t i, j;
	struct Chain *chain;
	struct PinchVertex *vertex;
	struct List *list;

	logDebug("Adding ends to this link of the net\n");

	//find the child ends recursively first
	for(i=0; i<chains->length; i++) {
		chain = chains->list[i];
		while(chain->nLink != NULL) {
			addEndsP(chain->adjacencyComponents, chain->subChains, chain->ends, ends);
			chain = chain->nLink;
		}
	}

	logDebug("Got child ends recursively\n");

	//now add the adjacency components
	for(i=0; i<adjacencyComponents->length; i++) {
		list = adjacencyComponents->list[i];
		for(j=0; j<list->length; j++) {
			vertex = list->list[j];
#ifdef BEN_DEBUG
			assert(lengthBlackEdges(vertex) > 0);
#endif
			if(isAStubOrCap(getFirstBlackEdge(vertex))) {
				listAppend(ends, vertex);
			}
		}
	}

#ifdef BEN_DEBUG
	i = ends->length;
	listRemoveDuplicates(ends);
	assert(i == ends->length);
#endif

	logDebug("Added the adjacency components\n");
	//now filter the stubs for the parent
	for(i=0; i<ends->length; i++) {
		vertex = ends->list[i];
#ifdef BEN_DEBUG
		assert(lengthBlackEdges(vertex) > 0);
		assert(isAStubOrCap(getFirstBlackEdge(vertex)) == TRUE);
#endif
		listAppend(parentEnds, vertex);
	}

	logDebug("Added the ends to this bit of the net\n");
}

void addEnds(struct Net *net) {
	/*
	 * Adds the ends pertinent to each link + base of the net.
	 * This allows us to trace the set of strings associated with each link + the base.
	 *
	 */
	struct List *list;

	logDebug("Starting to add ends to net\n");

	list = constructEmptyList(0, NULL);
	addEndsP(net->adjacencyComponents, net->chains, net->ends, list);
	destructList(list);

	logDebug("Finished adding ends to net\n");
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Choosing which atoms in the chains to keep.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct ContigEventSet *constructContigEventSet() {
	struct ContigEventSet *contigEventSet;

	contigEventSet = malloc(sizeof(struct ContigEventSet));
	contigEventSet->contigs = create_hashtable(0,
			hashtable_stringHashKey, hashtable_stringEqualKey,
			free, NULL);
	contigEventSet->event = NULL;
	contigEventSet->branchLength = 0.0;

	return contigEventSet;
}

void destructContigEventSet(struct ContigEventSet *contigEventSet) {
	hashtable_destroy(contigEventSet->contigs, FALSE, TRUE);
	free(contigEventSet->event);
	free(contigEventSet);
}

//first get tree covering score for each atom -
//drop all atoms with score less than X.
//accept chains whose remaining element's combined length is greater than a set length.

float treeCoverage(struct Chain *chain, struct List *nodeSets,
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
	assert(!isAStubOrCapChain(chain, pinchGraph));
#endif
	for(i=0; i<nodeSets->length; i++) {
		contigEventSet = nodeSets->list[i];
		totalTreeBranchLength += contigEventSet->branchLength;
		k = 0;
		for(j=0; j<chain->segments->length; j++) {
			segment = chain->segments->list[j];
#ifdef BEN_DEBUG
			assert(segment->contig < contigIndexToContigStrings->length);
#endif
			cA = contigIndexToContigStrings->list[segment->contig];
			if(hashtable_search(contigEventSet->contigs, cA) != NULL) {
				k += 1;
			}
		}
		if(k == chain->segments->length) {
			treeCoverage = 0.0; //reset as all segments must be contained within children of this node
		}
		else if(k > 0) {
			treeCoverage += contigEventSet->branchLength;
		}
	}
#ifdef BEN_DEBUG
	assert(chain->segments->length > 0);
	assert(treeCoverage >= 0);
	assert(treeCoverage <= totalTreeBranchLength);
#endif
	return treeCoverage/totalTreeBranchLength;
}

float atomScore(struct Chain *chain, struct List *nodeSets,
				struct List *contigIndexToContigStrings, struct PinchGraph *pinchGraph) {
	struct Segment *segment = chain->segments->list[0];
	return (segment->end - segment->start + 1) * treeCoverage(chain, nodeSets, contigIndexToContigStrings, pinchGraph);
}

int32_t chainScore(struct Chain *chain, struct List *nodeSets,
		struct List *contigIndexToContigStrings, struct PinchGraph *pinchGraph) {
	/*
	 * Gets the 'length' of the chain, including only those atoms in the chain that are
	 * true according to the 'includeAtom' function.
	 *
	 * Length is total length of atoms (all sequences in an atom have the same length).
	 */
	int32_t totalScore = 0;

	while(chain != NULL) {
		if(!isAStubOrCapChain(chain, pinchGraph)) {
			totalScore += atomScore(chain, nodeSets, contigIndexToContigStrings, pinchGraph);
		}
		chain = chain->nLink;
	}
	return totalScore;
}

int
compareFloats (const float *a, const float *b)
{
  return (int) (*a - *b);
}

int32_t chainLength(struct Chain *chain, int32_t includeStubs, struct PinchGraph *pinchGraph) {
	/*
	 * Get the number of links in the chain.
	 */
	int32_t i = 0;
	while(chain != NULL) {
		if(includeStubs || !isAStubOrCapChain(chain, pinchGraph)) {
			i++;
		}
		chain = chain->nLink;
	}
	return i;
}

int32_t chainBaseLength(struct Chain *chain, struct PinchGraph *pinchGraph) {
	/*
	 * Get the number of links in the chain.
	 */
	int32_t i = 0;
	struct Segment *segment;
	while(chain != NULL) {
		if(!isAStubOrCapChain(chain, pinchGraph)) {
			segment = chain->segments->list[0];
			i += segment->end - segment->start + 1;
		}
		chain = chain->nLink;
	}
	return i;
}

void filterAtomsByTreeCoverageAndLength(struct List *chainsList,
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
	int32_t i;
	struct Chain *chain;
	float *fA;
	float minScore;
	float minAtomScore;
	float f;

#ifdef BEN_DEBUG
	assert(proportionToKeep <= 1.0);
	assert(proportionToKeep >= 0.0);
	assert(discardRatio >= 0.0);
#endif

	///////////////
	//Calculate the minimum score to be kept.
	///////////////
	fA = mallocLocal(sizeof(float)*chainsList->length);
	for(i=0; i<chainsList->length; i++) {
		chain = chainsList->list[i];
		fA[i] = chainScore(chain, nodeSets, contigIndexToContigStrings, pinchGraph);
	}
	qsort(fA, chainsList->length, sizeof(float),
			(int (*)(const void *v, const void *))compareFloats);
	i = (1.0 - proportionToKeep)*chainsList->length;
	minScore = i < chainsList->length ? fA[i] : 10000000000000.0;
	logInfo("The minimum chain score to be kept is %f and %f \n", minScore, proportionToKeep);
	free(fA);

	///////////////
	//Gets those atoms whose chain meet the minimum score and have at
	///////////////

	for(i=0; i<chainsList->length; i++) {
		chain = chainsList->list[i];
		f = chainScore(chain, nodeSets, contigIndexToContigStrings, pinchGraph);
		if(f >= minScore && chainBaseLength(chain, pinchGraph) >= minimumChainLength) {
			minAtomScore = (minScore / chainLength(chain, FALSE, pinchGraph)) * discardRatio;
			while(chain != NULL) {
				if(!isAStubOrCapChain(chain, pinchGraph)) {
					if(treeCoverage(chain, nodeSets, contigIndexToContigStrings, pinchGraph) >= minimumTreeCoverage &&
						atomScore(chain, nodeSets, contigIndexToContigStrings, pinchGraph) >= minAtomScore) {
						listAppend(chosenAtoms, chain);
					}
				}
				chain = chain->nLink;
			}
		}
	}
}

void filterAtomsByIfStubOrCap(struct List *chainsList,
		struct List *chosenAtoms, struct PinchGraph *pinchGraph) {
	/*
	 * Adds all chains which represent either stubs or caps to the list of chosen atoms.
	 */
	int32_t i;
	struct Chain *chain;

	for(i=0; i<chainsList->length; i++) {
		chain = (struct Chain *)chainsList->list[i];
		while(chain != NULL) {
			if(isAStubOrCapChain(chain, pinchGraph)) {
				listAppend(chosenAtoms, chain);
			}
			chain = chain->nLink;
		}
	}
}

void logTheChosenAtomSubset(struct List *allAtoms, struct List *chosenAtoms, struct PinchGraph *pinchGraph,
		struct List *nodeSets, struct List *contigIndexToContigStrings) {
	/*
	 * Produces logging information about the chosen atoms.
	 */
	int32_t i, j;
	struct Chain *chain;
	struct Segment *segment;
	float totalAtomScore = 0.0;
	float totalAtomLength = 0.0;
	float totalBaseLengthOfAllAtoms = 0.0;
	float totalNumberOfAllAtoms = 0.0;
	float totalNumberOfStubAtoms = 0.0;
	float averageSegmentNumber = 0.0;
	float averageSegmentNumberOfAllAtoms = 0.0;
	for(i=0; i<allAtoms->length; i++) {
		chain = allAtoms->list[i];
		totalBaseLengthOfAllAtoms += chainBaseLength(chain, pinchGraph);
		totalNumberOfAllAtoms += chainLength(chain, FALSE, pinchGraph);
		totalNumberOfStubAtoms += chainLength(chain, TRUE, pinchGraph) - chainLength(chain, FALSE, pinchGraph);
		while(chain != NULL) {
			if(!isAStubOrCapChain(chain, pinchGraph)) {
				averageSegmentNumberOfAllAtoms += chain->segments->length;
			}
			chain = chain->nLink;
		}
	}
	j = 0;
	for(i=0; i<chosenAtoms->length; i++) {
		chain = chosenAtoms->list[i];
		if(!isAStubOrCapChain(chain, pinchGraph)) {
			totalAtomScore += treeCoverage(chain, nodeSets, contigIndexToContigStrings, pinchGraph);
			segment = chain->segments->list[0];
			totalAtomLength += segment->end - segment->start + 1;
			averageSegmentNumber += chain->segments->length;
			j++;
		}
	}
	logInfo("Chosen atom subset composed of %i atoms, of average length %f and average tree coverage %f, total base length of all atoms: %f, total number of all atoms %f, average length of all atoms: %f, total number of stub atoms: %f, average segment number of chosen atoms: %f, average segment number of all atoms: %f\n",
				chosenAtoms->length, totalAtomLength/j, totalAtomScore/j, totalBaseLengthOfAllAtoms, totalNumberOfAllAtoms, totalBaseLengthOfAllAtoms/totalNumberOfAllAtoms, totalNumberOfStubAtoms, averageSegmentNumber/j, averageSegmentNumberOfAllAtoms/totalNumberOfAllAtoms);
}
