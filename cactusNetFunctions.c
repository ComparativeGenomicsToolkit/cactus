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

struct PinchEdge *hookUpEdge(struct Segment *segment, struct PinchGraph *pinchGraph,
		struct PinchVertex *vertex2, struct PinchVertex *vertex3) {
	struct PinchEdge *edge;

	edge = constructPinchEdge(segment);

	//Connect up each end of the black edge.
	edge->from = vertex2;
	edge->rEdge->to = vertex2;
	insertBlackEdge(vertex2, edge);

	edge->to = vertex3;
	edge->rEdge->from = vertex3;
	insertBlackEdge(vertex3, edge->rEdge);

	//Now add segments connected to edges to the graph.
	avl_insert(pinchGraph->edges, edge);
	avl_insert(pinchGraph->edges, edge->rEdge);

	return edge;
}

struct PinchGraph *constructPinchGraph(Net *net) {
	struct PinchGraph *graph;
	Net_EndIterator *endIterator;
	End_InstanceIterator *instanceIterator;
	End *end;
	EndInstance *endInstance;
	EndInstance *endInstance2;
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
		if(end_isCap(end)) {
			connectVertices(sourceVertex, pinchVertex);
		}
		hashtable_insert(hash, netMisc_nameToString(end_getName(end)), pinchVertex);
		hashtable_insert(hash2, netMisc_nameToString(end_getName(end)), pinchVertex2);
	}
	net_destructEndIterator(endIterator);

	endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		instanceIterator = end_getInstanceIterator(end);
		while((endInstance = end_getNext(instanceIterator)) != NULL) {
			endInstance = endInstance_getStrand(endInstance) ? endInstance : endInstance_getReverse(endInstance);
			endInstance2 = endInstance_getAdjacency(endInstance);
			sequence = endInstance_getSequence(endInstance);

			assert(endInstance2 != NULL);
			assert(endInstance_getStrand(endInstance2));
			assert(sequence == endInstance_getSequence(endInstance2));

			//if(length >= 0)  {
			if(!endInstance_getSide(endInstance)) {
				assert(endInstance_getSide(endInstance2));

				start = endInstance_getCoordinate(endInstance);
				stop = endInstance_getCoordinate(endInstance2);
				length = stop - start - 1;
				assert(length >= 0);

				//Make black edges for caps/stubs on left end
				leftCapEdge = hookUpEdge(constructSegment(endInstance_getName(endInstance), start, start), graph,
						hashtable_search(hash, (void *)netMisc_nameToStringStatic(end_getName(endInstance_getEnd(endInstance)))),
						hashtable_search(hash2, (void *)netMisc_nameToStringStatic(end_getName(endInstance_getEnd(endInstance)))));


				//Construct the middle sequence, if not zero length.
				if(length > 0) {
					edge = hookUpEdge(constructSegment(sequence_getName(sequence), start+1, stop - 1), graph,
							constructPinchVertex(graph, -1, 0, 0), constructPinchVertex(graph, -1, 0, 0));
				}

				//Construct the right cap/stub
				rightCapEdge = hookUpEdge(constructSegment(endInstance_getName(endInstance2), stop, stop), graph,
						hashtable_search(hash2, (void *)netMisc_nameToStringStatic(end_getName(endInstance_getEnd(endInstance2)))),
						hashtable_search(hash, (void *)netMisc_nameToStringStatic(end_getName(endInstance_getEnd(endInstance2)))));

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


struct CactusEdge *getNonDeadEndOfStubOrCapCactusEdge(struct CactusEdge *edge, struct PinchGraph *pinchGraph) {
	struct PinchEdge *pinchEdge;
	pinchEdge = cactusEdgeToFirstPinchEdge(edge, pinchGraph);
	assert(isAStubOrCapCactusEdge(edge, pinchGraph));
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

Atom *constructAtomFromCactusEdge(struct CactusEdge *edge, Net *net) {
	/*
	 * Constructs an atom and two connected ends.
	 */
	int32_t i;
	Atom *atom;
	Sequence *sequence;
	struct Segment *segment;
	segment = edge->segments->list[0];
	atom = atom_construct(segment->end - segment->start + 1, net);
	for(i=0; i<edge->segments->length; i++) {
		segment = edge->segments->list[i];
		sequence = copySequence(net, segment->contig);
		atomInstance_construct2(atom, segment->start > 0 ? segment->start : -segment->end, segment->start > 0,
				sequence);
	}
	return atom;
}

struct List *addEnvelopedStubEnds(Net *net, int32_t addToNet) {
	/*
	 * For each net contained within a link in a chain, adds the encompassing ends
	 * of the chain to the nested net.
	 */
	int32_t j;
	End *end;
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
					end_copyConstruct(end, net);
				}
				else {
					assert(net_getEnd(net, end_getName(end)) != NULL);
				}
			}
			destructList(list);
			group_updateContainedEnds(group);
		}
	}
	net_destructGroupIterator(adjacencyIterator);

	list = constructEmptyList(0, NULL);
	endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(end_isStub(end)) {
			listAppend(list, end);
		}
	}
	net_destructEndIterator(endIterator);
	return list;
}

void addAdjacenciesToEnds(Net *net) {
	netMisc_addAdjacenciesToLeafEndInstances(net);
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
		EndInstance *endInstance;
		while((endInstance = end_getNext(iterator)) != NULL) {
			EndInstance *endInstance2 = endInstance_getAdjacency(endInstance);
			assert(endInstance2 != NULL);
			uint32_t j = abs(endInstance_getCoordinate(endInstance) - endInstance_getCoordinate(endInstance2));
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
	if(net_getGroupNumber(net) > 0 || net_getAtomNumber(net) > 0) {
		endIterator = net_getEndIterator(net);
		while((end = net_getNextEnd(endIterator)) != NULL) {
			assert(end_getGroup(end) != NULL);
		}
		net_destructEndIterator(endIterator);
	}
#endif
}

void addGroups(Net *net, struct PinchGraph *pinchGraph,
		struct List *chosenAtoms, struct hashtable *endNamesHash) {
	int32_t i, j;
	struct List *chosenPinchEdges = constructEmptyList(0, NULL);
	for(i=0; i<chosenAtoms->length; i++) {
		listAppend(chosenPinchEdges, cactusEdgeToFirstPinchEdge(chosenAtoms->list[i], pinchGraph));
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
		struct List *chosenAtoms) {
	Net *net;
	Net *nestedNet;
	End *end;
	End *end2;
	Atom *atom;
	Net_EndIterator *endIterator;
	EndInstance *endInstance;
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
	struct hashtable *chosenAtomsHash;
	struct hashtable *endNamesHash;
	struct PinchEdge *pinchEdge;
	struct Segment *segment;

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

	endNamesHash = create_hashtable(chosenAtoms->length*2,
			 hashtable_key, hashtable_equalKey,
			 NULL, free);
	endIterator = net_getEndIterator(parentNet);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		endInstance = end_getFirst(end);
		pinchEdge = getContainingBlackEdge(pinchGraph, endInstance_getName(endInstance), endInstance_getCoordinate(endInstance));
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

	chosenAtomsHash = create_hashtable(chosenAtoms->length*2,
							 hashtable_key, hashtable_equalKey,
							 NULL, NULL);
	for(i=0; i<chosenAtoms->length; i++) {
		hashtable_insert(chosenAtomsHash, chosenAtoms->list[i], &i);
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
			if((!isAStubOrCapCactusEdge(cactusEdge, pinchGraph)) && hashtable_search(chosenAtomsHash, cactusEdge) == NULL) {
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
	logDebug("Built the chosen atoms hash\n");

	////////////////////////////////////////////////
	//Atoms and ends for each net.
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
			//get the net.
			net = nets[mergedVertexIDs[cactusEdge->from->vertexID]];
			if(net == NULL) {
				net = net_construct(net_getNetDisk(parentNet));
				eventTree_copyConstruct(net_getEventTree(parentNet), net, returnsTrue);
				nets[mergedVertexIDs[cactusEdge->from->vertexID]] = net;
			}
			parentNets[i] = net;

			//Make the atoms and ends
			for(j=0; j<biConnectedComponent->length; j++) {
				cactusEdge = biConnectedComponent->list[j];
				segment = cactusEdge->segments->list[0];
				if(!isAStubOrCapCactusEdge(cactusEdge, pinchGraph)) {
					atom = constructAtomFromCactusEdge(cactusEdge, net);
					pinchEdge = cactusEdgeToFirstPinchEdge(cactusEdge, pinchGraph);
					hashtable_insert(endNamesHash, pinchEdge->from, netMisc_nameToString(end_getName(atom_getLeftEnd(atom))));
					hashtable_insert(endNamesHash, pinchEdge->to, netMisc_nameToString(end_getName(atom_getRightEnd(atom))));
					assert(cactusEdgeToEndName(cactusEdge, endNamesHash, pinchGraph) == end_getName(atom_getLeftEnd(atom)));
					assert(cactusEdgeToEndName(cactusEdge->rEdge, endNamesHash, pinchGraph) == end_getName(atom_getRightEnd(atom)));
				}
				else {
					assert(j == 0 || j == biConnectedComponent->length-1);
					cactusEdge2 = getNonDeadEndOfStubOrCapCactusEdge(cactusEdge, pinchGraph);
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
	logDebug("Constructed atoms and nets for the cycle.\n");

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

					end = net_getEnd(net, cactusEdgeToEndName(isAStubOrCapCactusEdge(cactusEdge, pinchGraph) ?
							getNonDeadEndOfStubOrCapCactusEdge(cactusEdge, pinchGraph) : cactusEdge->rEdge, endNamesHash, pinchGraph));
					assert(end != NULL);
					group_addEnd(group, end);

					end2 = net_getEnd(net, cactusEdgeToEndName(isAStubOrCapCactusEdge(cactusEdge2, pinchGraph) ?
							getNonDeadEndOfStubOrCapCactusEdge(cactusEdge2, pinchGraph) : cactusEdge2, endNamesHash, pinchGraph));
					assert(end2 != NULL);
					group_addEnd(group, end2);
				}
				else { //construct a link between two existing chains.
					assert(net_getEndNumber(nestedNet) > 0);
					end = net_getEnd(net, cactusEdgeToEndName(isAStubOrCapCactusEdge(cactusEdge, pinchGraph) ?
							getNonDeadEndOfStubOrCapCactusEdge(cactusEdge, pinchGraph) : cactusEdge->rEdge, endNamesHash, pinchGraph));
					assert(end != NULL);
					end_copyConstruct(end, nestedNet);
					end2 = net_getEnd(net, cactusEdgeToEndName(isAStubOrCapCactusEdge(cactusEdge2, pinchGraph) ?
							getNonDeadEndOfStubOrCapCactusEdge(cactusEdge2, pinchGraph) : cactusEdge2, endNamesHash, pinchGraph));
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

	addGroups(net, pinchGraph, chosenAtoms, endNamesHash);
	logDebug("Added the trivial groups\n");

	////////////////////////////////////////////////
	//Clean up
	////////////////////////////////////////////////

	free(nets);
	free(mergedVertexIDs);
	free(vertexDiscoveryTimes);
	free(parentNets);
	destructList(biConnectedComponents);
	hashtable_destroy(chosenAtomsHash, FALSE, FALSE);
	hashtable_destroy(endNamesHash, TRUE, FALSE);
}

void copyEndTreePhylogenies(Net *parentNet, Net *net) {
	/*
	 * For each end in the parent net that is also in the net, we copy
	 * the the phylogenetic information across.
	 */
	End *end1;
	End *end2;
	EndInstance *endInstance1;
	EndInstance *endInstance2;
	EndInstance *endInstance3;
	EndInstance *endInstance4;
	Net_EndIterator *endIterator;
	End_InstanceIterator *instanceIterator;

	endIterator = net_getEndIterator(net);
	while((end1 = net_getNextEnd(endIterator)) != NULL) {
		end2 = net_getEnd(parentNet, end_getName(end1));
		assert(end2 != NULL);
		instanceIterator = end_getInstanceIterator(end1);
		while((endInstance1 = end_getNext(instanceIterator)) != NULL) {
			assert(endInstance_getParent(endInstance1) == NULL);
			endInstance2 = end_getInstance(end2, endInstance_getName(endInstance1));
			assert(endInstance2 != NULL);
			if((endInstance3 = endInstance_getParent(endInstance2)) != NULL) {
				endInstance4 = end_getInstance(end1, endInstance_getName(endInstance3));
				assert(endInstance4 != NULL);
				endInstance_makeParentAndChild(endInstance4, endInstance1);
			}
			else {
				assert(end_getRootInstance(end2) == endInstance2);
			}
		}
		end_destructInstanceIterator(instanceIterator);
	}
	net_destructEndIterator(endIterator);
}
