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
#include "cactusNetFunctions.h"
#include "avl.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions to construct pinch graphs from nets.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchEdge *hookUpEdge(struct Piece *piece, struct PinchGraph *pinchGraph,
		struct PinchVertex *vertex2, struct PinchVertex *vertex3) {
	struct PinchEdge *edge;

	edge = constructPinchEdge(piece);

	//Connect up each end of the black edge.
	edge->from = vertex2;
	edge->rEdge->to = vertex2;
	insertBlackEdge(vertex2, edge);

	edge->to = vertex3;
	edge->rEdge->from = vertex3;
	insertBlackEdge(vertex3, edge->rEdge);

	//Now add pieces connected to edges to the graph.
	avl_insert(pinchGraph->edges, edge);
	avl_insert(pinchGraph->edges, edge->rEdge);

	return edge;
}

struct PinchGraph *constructPinchGraph(Net *net) {
	struct PinchGraph *graph;
	Net_EndIterator *endIterator;
	End_InstanceIterator *instanceIterator;
	End *end;
	Cap *cap;
	Cap *cap2;
	Sequence *sequence;
	struct PinchVertex *sourceVertex;
	struct PinchVertex *pinchVertex;
	struct PinchVertex *pinchVertex2;
	struct PinchEdge *leftCapEdge;
	struct PinchEdge *edge;
	struct PinchEdge *rightCapEdge;
	struct hashtable *hash;
	struct hashtable *hash2;
	int32_t start;
	int32_t stop;
	int32_t length;

	//make basic object.
	graph = pinchGraph_construct();
	sourceVertex = graph->vertices->list[0];

	//make hashes for ends to vertices
	hash = create_hashtable(net_getEndNumber(net) * 2,
				 hashtable_stringHashKey, hashtable_stringEqualKey,
				 free, NULL);
	hash2 = create_hashtable(net_getEndNumber(net) * 2,
					 hashtable_stringHashKey, hashtable_stringEqualKey,
					 free, NULL);

	//for each cap, build a pair of vertices
	endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		pinchVertex = constructPinchVertex(graph, -1, 0, 1);
		pinchVertex2 = constructPinchVertex(graph, -1, 1, 0);
		//connect to source.
		if(end_isAttached(end)) {
			connectVertices(sourceVertex, pinchVertex);
		}
		hashtable_insert(hash, netMisc_nameToString(end_getName(end)), pinchVertex);
		hashtable_insert(hash2, netMisc_nameToString(end_getName(end)), pinchVertex2);
	}
	net_destructEndIterator(endIterator);

	endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		instanceIterator = end_getInstanceIterator(end);
		while((cap = end_getNext(instanceIterator)) != NULL) {
			cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
			cap2 = cap_getAdjacency(cap);
			sequence = cap_getSequence(cap);

			assert(cap2 != NULL);
			assert(cap_getStrand(cap2));
			assert(sequence == cap_getSequence(cap2));

			//if(length >= 0)  {
			if(!cap_getSide(cap)) {
				assert(cap_getSide(cap2));

				start = cap_getCoordinate(cap);
				stop = cap_getCoordinate(cap2);
				length = stop - start - 1;
				assert(length >= 0);

				//Make black edges for caps/stubs on left end
				leftCapEdge = hookUpEdge(constructPiece(cap_getName(cap), start, start), graph,
						hashtable_search(hash, (void *)netMisc_nameToStringStatic(end_getName(cap_getEnd(cap)))),
						hashtable_search(hash2, (void *)netMisc_nameToStringStatic(end_getName(cap_getEnd(cap)))));


				//Construct the middle sequence, if not zero length.
				if(length > 0) {
					edge = hookUpEdge(constructPiece(sequence_getName(sequence), start+1, stop - 1), graph,
							constructPinchVertex(graph, -1, 0, 0), constructPinchVertex(graph, -1, 0, 0));
				}

				//Construct the right cap/stub
				rightCapEdge = hookUpEdge(constructPiece(cap_getName(cap2), stop, stop), graph,
						hashtable_search(hash2, (void *)netMisc_nameToStringStatic(end_getName(cap_getEnd(cap2)))),
						hashtable_search(hash, (void *)netMisc_nameToStringStatic(end_getName(cap_getEnd(cap2)))));

				//Connect the edges
				if(length > 0) {
					connectVertices(leftCapEdge->to, edge->from);
					connectVertices(edge->to, rightCapEdge->from);
				}
				else {
					connectVertices(leftCapEdge->to, rightCapEdge->from);
				}
			}
		}
		end_destructInstanceIterator(instanceIterator);
	}
	net_destructEndIterator(endIterator);
	//Cleanup the hashes
	hashtable_destroy(hash, FALSE, TRUE);
	hashtable_destroy(hash2, FALSE, TRUE);
	return graph;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions to construct nets.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


struct CactusEdge *getNonDeadEndOfStubCactusEdge(struct CactusEdge *edge, struct PinchGraph *pinchGraph) {
	struct PinchEdge *pinchEdge;
	pinchEdge = cactusEdgeToFirstPinchEdge(edge, pinchGraph);
	assert(isAStubCactusEdge(edge, pinchGraph));
	assert(vertex_isDeadEnd(pinchEdge->from) || vertex_isDeadEnd(pinchEdge->to));
	return vertex_isDeadEnd(pinchEdge->from) ? edge->rEdge : edge;
}

Name cactusEdgeToEndName(struct CactusEdge *edge, struct hashtable *endNamesHash, struct PinchGraph *pinchGraph) {
	struct PinchEdge *pinchEdge = cactusEdgeToFirstPinchEdge(edge, pinchGraph);
	char *cA = (char *)hashtable_search(endNamesHash, pinchEdge->from);
	assert(cA != NULL);
	return netMisc_stringToName(cA);
}

Sequence *copySequence(Net *net, Name name) {
	Sequence *sequence = net_getSequence(net, name);
	if(sequence == NULL) {
		sequence = sequence_construct(netDisk_getMetaSequence(net_getNetDisk(net), name), net);
	}
	return sequence;
}

Block *constructBlockFromCactusEdge(struct CactusEdge *edge, Net *net) {
	/*
	 * Constructs an block and two connected ends.
	 */
	int32_t i;
	Block *block;
	Sequence *sequence;
	struct Piece *piece;
	piece = edge->pieces->list[0];
	block = block_construct(piece->end - piece->start + 1, net);
	for(i=0; i<edge->pieces->length; i++) {
		piece = edge->pieces->list[i];
		sequence = copySequence(net, piece->contig);
		segment_construct2(block, piece->start > 0 ? piece->start : -piece->end, piece->start > 0, sequence);
	}
	return block;
}

struct List *addEnvelopedStubEnds(Net *net, int32_t addToNet) {
	/*
	 * For each net contained within a link in a chain, adds the encompassing ends
	 * of the chain to the nested net.
	 */
	int32_t j;
	End *end, *end2;
	Net *net2;
	struct List *list;
	Net_EndIterator *endIterator;
	Group_EndIterator *adjacencyIterator;
	Group *group;

	adjacencyIterator = net_getGroupIterator(net);
	while((group = net_getNextGroup(adjacencyIterator)) != NULL) {
		net2 = group_getNestedNet(group);
		if(net2 != NULL) {
			list = addEnvelopedStubEnds(net2, 1);
			for(j=0; j<list->length; j++) {
				end = list->list[j];
				if(addToNet && net_getEnd(net, end_getName(end)) == NULL) {
					group_addEnd(group, end_copyConstruct(end, net));
				}
				else {
					end2 = net_getEnd(net, end_getName(end));
					assert(end2 != NULL);
					if(end_getGroup(end2) == NULL) {
						group_addEnd(group, end2);
					}
					else {
						assert(end_getGroup(end2) == group);
					}
				}
			}
			destructList(list);
		}
	}
	net_destructGroupIterator(adjacencyIterator);

	list = constructEmptyList(0, NULL);
	endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(end_isStubEnd(end)) {
			listAppend(list, end);
		}
	}
	net_destructEndIterator(endIterator);
	return list;
}

void addAdjacenciesToEnds(Net *net) {
	netMisc_addAdjacenciesToLeafCaps(net);
	Group_EndIterator *adjacencyIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(adjacencyIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			addAdjacenciesToEnds(group_getNestedNet(group));
		}
	}
	net_destructGroupIterator(adjacencyIterator);
}

static int32_t returnsTrue(Event *event) {
	assert(event != NULL);
	return 1;
}

bool groupIsZeroLength(struct List *endNames, Net *net) {
	int32_t i;

	for(i=0; i<endNames->length; i++) {
		End *end = net_getEnd(net, netMisc_stringToName(endNames->list[i]));
		End_InstanceIterator *iterator = end_getInstanceIterator(end);
		Cap *cap;
		while((cap = end_getNext(iterator)) != NULL) {
			Cap *cap2 = cap_getAdjacency(cap);
			assert(cap2 != NULL);
			uint32_t j = abs(cap_getCoordinate(cap) - cap_getCoordinate(cap2));
			assert(j != 0);
			if(j > 1) {
				return 0;
			}
		}
		end_destructInstanceIterator(iterator);
	}
	return 1;
}

void addGroupsP(Net *net, struct hashtable *groups) {
	/*
	 * Adds the non chain groups to each net.
	 */
	Net_EndIterator *endIterator;
	Net_GroupIterator *adjacencyIterator;
	Group_EndIterator *endIterator2;
	End *end;
	End *end2;
	Group *group;
	struct List *endNames;
	int32_t i;

	//Do for each level of the net.
	adjacencyIterator = net_getGroupIterator(net);
	while((group = net_getNextGroup(adjacencyIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			addGroupsP(group_getNestedNet(group), groups);
		}
		else { //Ends are already in an terminal group!
			endIterator2 = group_getEndIterator(group);
			endNames = NULL;
			while((end = group_getNextEnd(endIterator2)) != NULL) {
				assert((endNames = hashtable_remove(groups, (void *)netMisc_nameToStringStatic(end_getName(end)), 0)) != NULL);
			}
			group_destructEndIterator(endIterator2);
			assert(endNames != NULL);
			destructList(endNames);
		}
	}
	net_destructGroupIterator(adjacencyIterator);

	endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		group = end_getGroup(end);
		if(group == NULL) {
			endNames = hashtable_search(groups, (void *)netMisc_nameToStringStatic(end_getName(end)));
			assert(endNames != NULL);
#ifdef BEN_DEBUG
			for(i=0; i<endNames->length; i++) {
				end2 = net_getEnd(net, netMisc_stringToName(endNames->list[i]));
				assert(end2 != NULL);
				assert(end_getGroup(end2) == NULL);
			}
#endif
			group = group_construct2(net);
			for(i=0; i<endNames->length; i++) {
				end2 = net_getEnd(net, netMisc_stringToName(endNames->list[i]));
				group_addEnd(group, end2);
			}
#ifdef BEN_DEBUG
			for(i=0; i<endNames->length; i++) {
				end2 = net_getEnd(net, netMisc_stringToName(endNames->list[i]));
				assert(end_getGroup(end2) == group);
			}
#endif
			for(i=0; i<endNames->length; i++) {
				assert(hashtable_remove(groups, endNames->list[i], 0) == endNames);
			}
			destructList(endNames);
		}
	}
	net_destructEndIterator(endIterator);

#ifdef BEN_DEBUG
	if(net_getGroupNumber(net) > 0 || net_getBlockNumber(net) > 0) {
		endIterator = net_getEndIterator(net);
		while((end = net_getNextEnd(endIterator)) != NULL) {
			assert(end_getGroup(end) != NULL);
		}
		net_destructEndIterator(endIterator);
	}
#endif
}

void addGroups(Net *net, struct PinchGraph *pinchGraph,
		struct List *chosenBlocks, struct hashtable *endNamesHash) {
	int32_t i, j;
	struct List *chosenPinchEdges = constructEmptyList(0, NULL);
	for(i=0; i<chosenBlocks->length; i++) {
		listAppend(chosenPinchEdges, cactusEdgeToFirstPinchEdge(chosenBlocks->list[i], pinchGraph));
		assert(chosenPinchEdges->list[i] != NULL);
	}
	for(i=0; i<pinchGraph->vertices->length; i++) {
		struct PinchVertex *vertex = pinchGraph->vertices->list[i];
		if(vertex_isDeadEnd(vertex) || vertex_isEnd(vertex)) {
			listAppend(chosenPinchEdges, getFirstBlackEdge(vertex));
		}
	}
	struct List *groupsList = getRecursiveComponents2(pinchGraph, chosenPinchEdges);
	destructList(chosenPinchEdges);
	struct hashtable *groupsHash =
			create_hashtable(groupsList->length * 2,
							 hashtable_stringHashKey, hashtable_stringEqualKey, NULL, NULL);
	for(i=0; i<groupsList->length; i++) {
		struct List *vertices = groupsList->list[i];
		struct List *endNames = constructEmptyList(0, NULL);
		assert(vertices->length > 0);
		struct PinchVertex *vertex = vertices->list[0];
		if(!vertex_isDeadEnd(vertex) && vertex->vertexID != 0) {
			for(j=0; j<vertices->length; j++) {
				vertex = vertices->list[j];
				assert(!vertex_isDeadEnd(vertex));
				assert(vertex->vertexID != 0);
				const char *endNameString = hashtable_search(endNamesHash, vertex);
				//assert(endNameString != NULL);
				if(endNameString != NULL) {
					listAppend(endNames, (void *)endNameString);
					hashtable_insert(groupsHash, (void *)endNameString, endNames);
				}
			}
			assert(endNames->length >= 1);
		}
#ifdef BEN_DEBUG
		else {
			for(j=0; j<vertices->length; j++) {
				vertex = vertices->list[j];
				assert(vertex_isDeadEnd(vertex) || vertex->vertexID == 0);
			}
		}
#endif
	}
	destructList(groupsList);
	addGroupsP(net, groupsHash);
	assert(hashtable_count(groupsHash) == 0);
	hashtable_destroy(groupsHash, FALSE, FALSE);
}


static int32_t *vertexDiscoveryTimes;

#ifdef BEN_DEBUG
void checkBiConnectedComponent(struct List *biConnnectedComponent) {
	int32_t i, j;
	struct CactusEdge *cactusEdge;
	struct CactusEdge *cactusEdge2;
	assert(biConnnectedComponent->length > 0);
	j = INT32_MIN;
	for(i=0; i<biConnnectedComponent->length; i++) {
		cactusEdge = biConnnectedComponent->list[i];
		assert(j < vertexDiscoveryTimes[cactusEdge->from->vertexID]);
		j = vertexDiscoveryTimes[cactusEdge->from->vertexID];
	}
	cactusEdge = biConnnectedComponent->list[0];
	cactusEdge2 = biConnnectedComponent->list[biConnnectedComponent->length-1];
	assert(vertexDiscoveryTimes[cactusEdge->from->vertexID] == vertexDiscoveryTimes[cactusEdge2->to->vertexID]);
}
#endif

int fillOutNetFromInputsP2(struct List **biConnectedComponent1, struct List **biConnectedComponent2) {
	struct CactusEdge *cactusEdge1;
	struct CactusEdge *cactusEdge2;
	int32_t i, j;
	cactusEdge1 = (*biConnectedComponent1)->list[0];
	cactusEdge2 = (*biConnectedComponent2)->list[0];
	i = vertexDiscoveryTimes[cactusEdge1->from->vertexID];
	j = vertexDiscoveryTimes[cactusEdge2->from->vertexID];
	return i - j;
}

void fillOutNetFromInputs(
		Net *parentNet,
		struct CactusGraph *cactusGraph,
		struct PinchGraph *pinchGraph,
		struct List *chosenBlocks) {
	Net *net;
	Net *nestedNet;
	End *end;
	End *end2;
	Block *block;
	Net_EndIterator *endIterator;
	Cap *cap;
	Chain *chain;
	Group *group;
	struct CactusVertex *cactusVertex;
	struct CactusEdge *cactusEdge;
	struct CactusEdge *cactusEdge2;
	struct List *biConnectedComponent;
	struct List *list;
	struct List *biConnectedComponents;
	void **nets;
	void **parentNets;
	int32_t *mergedVertexIDs;
	int32_t i, j, k;
	struct hashtable *chosenBlocksHash;
	struct hashtable *endNamesHash;
	struct PinchEdge *pinchEdge;
	struct Piece *piece;

	logDebug("Building the net\n");

	////////////////////////////////////////////////
	//Get sorted bi-connected components (sorted as in ordered from root of vertex)
	////////////////////////////////////////////////

	biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);

	logDebug("Built the bi-connected components\n");

	////////////////////////////////////////////////
	//Get DFS numbering on cactus vertices
	////////////////////////////////////////////////

	vertexDiscoveryTimes = getDFSDiscoveryTimes(cactusGraph);
#ifdef BEN_DEBUG
	for(i=0; i<biConnectedComponents->length; i++) { //checks the discovery times are as expected.
		checkBiConnectedComponent(biConnectedComponents->list[i]);
	}
#endif
	logDebug("Got the vertex discovery times\n");

	////////////////////////////////////////////////
	//Sort the biconnected components by their start time
	////////////////////////////////////////////////

	qsort(biConnectedComponents->list, biConnectedComponents->length, sizeof(void *),
			(int (*)(const void *v, const void *))fillOutNetFromInputsP2);
	logDebug("Sorted the biconnected components by vertex discovery time");

	////////////////////////////////////////////////
	//Build end names hash
	////////////////////////////////////////////////

	endNamesHash = create_hashtable(chosenBlocks->length*2,
			 hashtable_key, hashtable_equalKey,
			 NULL, free);
	endIterator = net_getEndIterator(parentNet);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		cap = end_getFirst(end);
		pinchEdge = getContainingBlackEdge(pinchGraph, cap_getName(cap), cap_getCoordinate(cap));
		assert(pinchEdge != NULL);
		if(vertex_isEnd(pinchEdge->from)) {
			assert(vertex_isDeadEnd(pinchEdge->to));
			hashtable_insert(endNamesHash, pinchEdge->from, netMisc_nameToString(end_getName(end)));
		}
		else {
			assert(vertex_isEnd(pinchEdge->to));
			assert(vertex_isDeadEnd(pinchEdge->from));
			hashtable_insert(endNamesHash, pinchEdge->to, netMisc_nameToString(end_getName(end)));
		}
	}
	net_destructEndIterator(endIterator);
	logDebug("Built the end names hash\n");

	////////////////////////////////////////////////
	//Prune the cactus graph to include only those edges relevant to the desired net.
	////////////////////////////////////////////////

	chosenBlocksHash = create_hashtable(chosenBlocks->length*2,
							 hashtable_key, hashtable_equalKey,
							 NULL, NULL);
	for(i=0; i<chosenBlocks->length; i++) {
		hashtable_insert(chosenBlocksHash, chosenBlocks->list[i], &i);
	}
	mergedVertexIDs = malloc(sizeof(int32_t)*cactusGraph->vertices->length);
	for(i=0; i<cactusGraph->vertices->length; i++) {
		mergedVertexIDs[i] = ((struct CactusVertex *)cactusGraph->vertices->list[i])->vertexID;
	}
	for(i=0; i<biConnectedComponents->length; i++) {
		biConnectedComponent = biConnectedComponents->list[i];
		list = constructEmptyList(0, NULL);
		for(j=0; j<biConnectedComponent->length; j++) {
			cactusEdge = biConnectedComponent->list[j];
			if((!isAStubCactusEdge(cactusEdge, pinchGraph)) && hashtable_search(chosenBlocksHash, cactusEdge) == NULL) {
				//merge vertices
				if(vertexDiscoveryTimes[cactusEdge->from->vertexID] < vertexDiscoveryTimes[cactusEdge->to->vertexID]) {
					mergedVertexIDs[cactusEdge->to->vertexID] = mergedVertexIDs[cactusEdge->from->vertexID];
				}
				else if (vertexDiscoveryTimes[cactusEdge->from->vertexID] > vertexDiscoveryTimes[cactusEdge->to->vertexID]) {
					assert(j == biConnectedComponent->length-1);
					for(k=0; k <= j; k++) {
						cactusVertex = ((struct CactusEdge *)biConnectedComponent->list[k])->from;
						if(mergedVertexIDs[cactusVertex->vertexID] == mergedVertexIDs[cactusEdge->from->vertexID]) {
							mergedVertexIDs[cactusVertex->vertexID] = mergedVertexIDs[cactusEdge->to->vertexID];
						}
					}
				}
			}
			else {
				listAppend(list, cactusEdge);
			}
		}
		destructList(biConnectedComponent);
		biConnectedComponents->list[i] = list;
	}
	logDebug("Built the chosen blocks hash\n");

	////////////////////////////////////////////////
	//Blocks and ends for each net.
	////////////////////////////////////////////////

	nets = mallocLocal(sizeof(void *) * cactusGraph->vertices->length);
	for(i=1; i<cactusGraph->vertices->length; i++) {
		nets[i] = NULL;
	}
	nets[0] = parentNet;
	parentNets = mallocLocal(sizeof(void *) * biConnectedComponents->length);
	for(i=0; i<biConnectedComponents->length; i++) {
		biConnectedComponent = biConnectedComponents->list[i];
		if(biConnectedComponent->length > 0) {
			cactusEdge = biConnectedComponent->list[0];
			//Get the net.
			net = nets[mergedVertexIDs[cactusEdge->from->vertexID]];
			if(net == NULL) {
				net = net_construct(net_getNetDisk(parentNet));
				eventTree_copyConstruct(net_getEventTree(parentNet), net, returnsTrue);
				nets[mergedVertexIDs[cactusEdge->from->vertexID]] = net;
			}
			parentNets[i] = net;

			//Make the blocks and ends
			for(j=0; j<biConnectedComponent->length; j++) {
				cactusEdge = biConnectedComponent->list[j];
				piece = cactusEdge->pieces->list[0];
				if(!isAStubCactusEdge(cactusEdge, pinchGraph)) {
					block = constructBlockFromCactusEdge(cactusEdge, net);
					pinchEdge = cactusEdgeToFirstPinchEdge(cactusEdge, pinchGraph);
					hashtable_insert(endNamesHash, pinchEdge->from, netMisc_nameToString(end_getName(block_getLeftEnd(block))));
					hashtable_insert(endNamesHash, pinchEdge->to, netMisc_nameToString(end_getName(block_getRightEnd(block))));
					assert(cactusEdgeToEndName(cactusEdge, endNamesHash, pinchGraph) == end_getName(block_getLeftEnd(block)));
					assert(cactusEdgeToEndName(cactusEdge->rEdge, endNamesHash, pinchGraph) == end_getName(block_getRightEnd(block)));
				}
				else {
					assert(j == 0 || j == biConnectedComponent->length-1);
					cactusEdge2 = getNonDeadEndOfStubCactusEdge(cactusEdge, pinchGraph);
					//if(j == 0) { //not using these asserts currently
						//assert(cactusEdge2 == cactusEdge->rEdge);
					//}
					//else {
						//assert(cactusEdge2 == cactusEdge);
					//}
					end = net_getEnd(parentNet, cactusEdgeToEndName(cactusEdge2, endNamesHash, pinchGraph));
					assert(end != NULL);
					if(net != parentNet) {
						end_copyConstruct(end, net);
					}
				}
			}
		}
	}
	logDebug("Constructed blocks and nets for the cycle.\n");

	////////////////////////////////////////////////
	//Link nets to parent nets and construct chains.
	////////////////////////////////////////////////

	for(i=0; i<biConnectedComponents->length; i++) {
		biConnectedComponent = biConnectedComponents->list[i];
		net = parentNets[i];
		if(biConnectedComponent->length > 1) {
			assert(net != NULL);
			chain = chain_construct(net);
			for(j=1; j<biConnectedComponent->length; j++) {
				cactusEdge = biConnectedComponent->list[j-1];
				cactusEdge2 = biConnectedComponent->list[j];
				nestedNet = nets[mergedVertexIDs[cactusEdge->to->vertexID]];
				assert(cactusEdge->to->vertexID != 0);
				if(nestedNet == NULL) { //construct a terminal group.
					group = group_construct2(net);

					end = net_getEnd(net, cactusEdgeToEndName(isAStubCactusEdge(cactusEdge, pinchGraph) ?
							getNonDeadEndOfStubCactusEdge(cactusEdge, pinchGraph) : cactusEdge->rEdge, endNamesHash, pinchGraph));
					assert(end != NULL);
					group_addEnd(group, end);

					end2 = net_getEnd(net, cactusEdgeToEndName(isAStubCactusEdge(cactusEdge2, pinchGraph) ?
							getNonDeadEndOfStubCactusEdge(cactusEdge2, pinchGraph) : cactusEdge2, endNamesHash, pinchGraph));
					assert(end2 != NULL);
					group_addEnd(group, end2);
				}
				else { //construct a link between two existing chains.
					assert(net_getEndNumber(nestedNet) > 0);
					end = net_getEnd(net, cactusEdgeToEndName(isAStubCactusEdge(cactusEdge, pinchGraph) ?
							getNonDeadEndOfStubCactusEdge(cactusEdge, pinchGraph) : cactusEdge->rEdge, endNamesHash, pinchGraph));
					assert(end != NULL);
					end_copyConstruct(end, nestedNet);
					end2 = net_getEnd(net, cactusEdgeToEndName(isAStubCactusEdge(cactusEdge2, pinchGraph) ?
							getNonDeadEndOfStubCactusEdge(cactusEdge2, pinchGraph) : cactusEdge2, endNamesHash, pinchGraph));
					assert(end2 != NULL);
					end_copyConstruct(end2, nestedNet);
					group = group_construct(net, nestedNet);
					nets[mergedVertexIDs[cactusEdge->to->vertexID]] = NULL;
				}
				//Make link chain
				link_construct(end, end2, group, chain);
			}
		}
	}
	net = nets[0];
#ifdef BEN_DEBUG
	for(i=1; i<cactusGraph->vertices->length; i++) {
		assert(nets[i] == NULL);
	}
#endif
	logDebug("Constructed the chains and linked together the nets\n");

	////////////////////////////////////////////////
	//Add nested ends to nets.
	////////////////////////////////////////////////

	destructList(addEnvelopedStubEnds(parentNet, 0));
	logDebug("Added the nested ends to the parent nets\n");

	////////////////////////////////////////////////
	//Add adjacencies between ends.
	////////////////////////////////////////////////

	addAdjacenciesToEnds(net);
	logDebug("Added the adjacencies between the ends\n");

	////////////////////////////////////////////////
	//Add groups.
	////////////////////////////////////////////////

	addGroups(net, pinchGraph, chosenBlocks, endNamesHash);
	logDebug("Added the trivial groups\n");

	////////////////////////////////////////////////
	//Clean up
	////////////////////////////////////////////////

	free(nets);
	free(mergedVertexIDs);
	free(vertexDiscoveryTimes);
	free(parentNets);
	destructList(biConnectedComponents);
	hashtable_destroy(chosenBlocksHash, FALSE, FALSE);
	hashtable_destroy(endNamesHash, TRUE, FALSE);
}

void copyEndTreePhylogenies(Net *parentNet, Net *net) {
	/*
	 * For each end in the parent net that is also in the net, we copy
	 * the the phylogenetic information across.
	 */
	End *end1;
	End *end2;
	Cap *cap1;
	Cap *cap2;
	Cap *cap3;
	Cap *cap4;
	Net_EndIterator *endIterator;
	End_InstanceIterator *instanceIterator;

	endIterator = net_getEndIterator(net);
	while((end1 = net_getNextEnd(endIterator)) != NULL) {
		end2 = net_getEnd(parentNet, end_getName(end1));
		assert(end2 != NULL);
		instanceIterator = end_getInstanceIterator(end1);
		while((cap1 = end_getNext(instanceIterator)) != NULL) {
			assert(cap_getParent(cap1) == NULL);
			cap2 = end_getInstance(end2, cap_getName(cap1));
			assert(cap2 != NULL);
			if((cap3 = cap_getParent(cap2)) != NULL) {
				cap4 = end_getInstance(end1, cap_getName(cap3));
				assert(cap4 != NULL);
				cap_makeParentAndChild(cap4, cap1);
			}
			else {
				assert(end_getRootInstance(end2) == cap2);
			}
		}
		end_destructInstanceIterator(instanceIterator);
	}
	net_destructEndIterator(endIterator);
}
