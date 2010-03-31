#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <iostream>

#include "XMLTools.h"
#include "xmlParser.h"
#include "reconstructionTree.h"

extern "C" {
	#include "pinchGraph.h"
	#include "commonC.h"
	#include "fastCMaths.h"
	#include "bioioC.h"
	#include "hashTableC.h"
	#include "net.h"
	#include "pairwiseAlignment.h"
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods to extract a pinch graph from a reconstruction problem
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchEdge *hookUpEdge(struct Segment *segment, struct PinchGraph *pinchGraph) {
	struct PinchEdge *edge;
	struct PinchVertex *vertex2;
	struct PinchVertex *vertex3;

	edge = constructPinchEdge(segment);

	//Connect up each end of the black edge.
	vertex2 = constructPinchVertex(pinchGraph, -1);
	edge->from = vertex2;
	edge->rEdge->to = vertex2;
	insertBlackEdge(vertex2, edge);

	vertex3 = constructPinchVertex(pinchGraph, -1);
	edge->to = vertex3;
	edge->rEdge->from = vertex3;
	insertBlackEdge(vertex3, edge->rEdge);

	//Now add segments connected to edges to the graph.
	avl_insert(pinchGraph->edges, edge);
	avl_insert(pinchGraph->edges, edge->rEdge);

	return edge;
}

struct PinchGraph *constructPinchGraph(XMLNode xMainNode,
									   struct List *contigIndexToContigStrings,
									   struct IntList *contigIndexToContigStart) {
	/*
	 * Method to construct the pinch graph from a reconstruction problem node.
	 */
	struct PinchGraph *pinchGraph;
	struct List *capEdges;
	struct List *capStrings;

	struct PinchEdge *leftCapEdge;
	struct PinchEdge *edge;
	struct PinchEdge *edge2;
	struct PinchEdge *rightCapEdge;
	struct PinchVertex *vertex1;

	int32_t i, start, length;
	const char *cA;
	char *cA2;
	struct Segment *segment;
	struct hashtable *caps;

	pinchGraph = (struct PinchGraph *)mallocLocal(sizeof(struct PinchGraph));
	pinchGraph->edges = avl_create((int32_t (*)(const void *, const void *, void *a))edgeComparator, NULL, NULL);
	pinchGraph->vertices = constructEmptyList(0, (void (*)(void *))destructPinchVertex);
	vertex1 = constructPinchVertex(pinchGraph, -1);
	logDebug("Setup the basic graph data structures ready to add the details.\n");

	XMLNode stringsNode = xMainNode.getChildNode("strings");

	capEdges = constructEmptyList(0, NULL);
	capStrings = constructEmptyList(0, NULL);

	logDebug("Got the strings node for the reconstruction problem, there are: %i strings\n", stringsNode.nChildNode("string"));

	//construct the black edges
	for(i=0; i<stringsNode.nChildNode("string"); i++) {
		XMLNode stringNode = stringsNode.getChildNode("string", i);

		assert(sscanf(stringNode.getAttribute("start"), "%i", &start) == 1);
		assert(sscanf(stringNode.getAttribute("length"), "%i", &length) == 1);

		//Construct left cap/stub
		cA = stringNode.getAttribute("left");
		listAppend(contigIndexToContigStrings, stringCopy(cA));
		intListAppend(contigIndexToContigStart, start+1);
		listAppend(capStrings, contigIndexToContigStrings->list[contigIndexToContigStrings->length-1]);
		segment = constructSegment(contigIndexToContigStrings->length-1, start+1, start+1);
		leftCapEdge = hookUpEdge(segment, pinchGraph);
		listAppend(capEdges, leftCapEdge);
		if(cA[0] != '$') {
			connectVertices(vertex1, leftCapEdge->from);
		}

		//Construct the middle sequence.
		//This order is important. It ensures that the contig index of the stub/cap instance
		//to the left of the sequence is always one less than
		//contig index of the actual sequence.
		//Similarly it ensures that the contig index of the stub/cap instance
		//to the right of the sequence is always one greater than the contig index of the actual sequence!
		//This allows us to unambiguously identify which stub/cap instances lie to the left and right of which sequences,
		//useful if we unalign edges and need to work out there original (grey) adjacencies with a stub.
		cA = stringNode.getAttribute("contig");
		listAppend(contigIndexToContigStrings, stringCopy(cA));
		intListAppend(contigIndexToContigStart, start+2);
		if(length > 0) {
			segment = constructSegment(contigIndexToContigStrings->length-1, start+2, start+length+1);
			edge = hookUpEdge(segment, pinchGraph);
		}
		else {
			edge = NULL;
		}

		//Construct the right cap/stub
		cA = stringNode.getAttribute("right");
		listAppend(contigIndexToContigStrings, stringCopy(cA));
		intListAppend(contigIndexToContigStart, start+length+2);
		listAppend(capStrings, contigIndexToContigStrings->list[contigIndexToContigStrings->length-1]);
		segment = constructSegment(contigIndexToContigStrings->length-1, start+length+2, start+length+2);
		rightCapEdge = hookUpEdge(segment, pinchGraph);
		listAppend(capEdges, rightCapEdge->rEdge);
		if(cA[0] != '$') {
			connectVertices(vertex1, rightCapEdge->to);
		}

		//Connect the edges
		if(length > 0) {
			assert(edge != NULL);
			connectVertices(leftCapEdge->to, edge->from);
			connectVertices(edge->to, rightCapEdge->from);
		}
		else {
			connectVertices(leftCapEdge->to, rightCapEdge->from);
		}
	}

	logDebug("Constructed the black edges from the pinch graph\n");

	caps = create_hashtable(pinchGraph->vertices->length*10,
							hashtable_stringHashKey, hashtable_stringEqualKey,
							free, NULL);

	//merge the telomere edges
	for(i=0; i<capStrings->length; i++) {
		cA = (char *)capStrings->list[i];
		edge = (PinchEdge *)capEdges->list[i];
#ifdef BEN_DEBUG
		isAStubOrCap(edge);
#endif

		cA2 = removeInstance((char *)cA);
		edge2 = (struct PinchEdge *)hashtable_search(caps, cA2);
		if(edge2 != NULL) {
			mergeVertices(pinchGraph, edge->from, edge2->from);
			mergeVertices(pinchGraph, edge->to, edge2->to);
			free(cA2);
		}
		else {
			hashtable_insert(caps, (void *)cA2, edge);
		}
	}
	hashtable_destroy(caps, FALSE, TRUE);

	destructList(capEdges);
	destructList(capStrings);
	return pinchGraph;
}

struct Sequence *constructSequence(const char *contig, const char *event, const char *file, int32_t length) {
	struct Sequence *sequence;

	sequence = (struct Sequence *)mallocLocal(sizeof(struct Sequence));
	sequence->contig = stringCopy(contig);
	sequence->event = stringCopy(event);
	sequence->file = stringCopy(file);
	sequence->length = length;

	return sequence;
}

void destructSequence(struct Sequence *sequence) {
	free(sequence->contig);
	free(sequence->event);
	free(sequence->file);
	free(sequence);
}

struct hashtable *parseSequences(XMLNode xMainNode) {
	int32_t i, j;
	struct hashtable *contigNamesToSequences;
	struct Sequence *sequence;

	XMLNode sequencesNode = xMainNode.getChildNode("sequences");
	contigNamesToSequences = create_hashtable(MEDIUM_CHUNK_SIZE,
				hashtable_stringHashKey, hashtable_stringEqualKey,
				NULL, (void (*)(void *))destructSequence);

	for(i=0; i<sequencesNode.nChildNode("sequence"); i++) {
		XMLNode sequenceNode = sequencesNode.getChildNode("sequence", i);
		assert(sscanf(sequenceNode.getAttribute("length"), "%i", &j) == 1);
		sequence = constructSequence(sequenceNode.getAttribute("contig"),
									 sequenceNode.getAttribute("event"),
									 sequenceNode.getAttribute("sequence_file"), j);
		hashtable_insert(contigNamesToSequences, sequence->contig, sequence);
	}

	return contigNamesToSequences;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for writing out a reconstruction tree.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

char *getAbsoluteFilePath(const char *absoluteFilePrefix, const char *relativeReconstructionFilePath) {
	/*
	 * Creates an absolute file path by concatenating an absolute and relative file path.
	 */
	char *cA;
	cA = (char *)mallocLocal(sizeof(char) * (strlen(absoluteFilePrefix) + strlen(relativeReconstructionFilePath) + 2));
	sprintf(cA, "%s/%s", absoluteFilePrefix, relativeReconstructionFilePath);
	return cA;
}

void convertToXMLCoordinate(int32_t i, char *coordinate, char *strand) {
	/*
	 * Takes a coordinate from the internal cactus graph coordinates and makes into a coordinate for the XML file tree.
	 */
	int32_t j;

#ifdef BEN_DEBUG
	assert(i >= 2 || i <= -2); //we can't deal with caps
#endif

	if(i >= 1) { //on positive strand
		sprintf(coordinate, "%i", i-2);
		sprintf(strand, "+");
	}
	else {
		j = (-i)-2;
#ifdef BEN_DEBUG
		assert(j >= 0);
#endif
		sprintf(coordinate, "%i", j);
		sprintf(strand, "-");
	}
}

char *getXMLContigString(int32_t i, struct List *contigIndexToContigStrings) {
	/*
	 * Gets the contig names associated with a contig index i.
	 */
	char *cA;

	cA = (char *)contigIndexToContigStrings->list[i - (i%3) + 1];
#ifdef BEN_DEBUG
	assert(cA != NULL);
#endif
	return cA;
}

char *getCapInstanceName(struct PinchEdge *edge, struct PinchVertex *vertex, struct hashtable *names) {
	/*
	 * Gets the name of a cap instance.
	 */
	char *cA;
	char *cA2;
	char *cA3;

#ifdef BEN_DEBUG
	assert(containsBlackEdge(vertex, edge));
#endif

	if(isAStubOrCap(edge) == TRUE) {
		cA = (char *)hashtable_search(names, edge);
		cA = stringCopy(cA);
	}
	else {
		//copy edge name
		cA2 = (char *)hashtable_search(names, vertex);
		cA3 = getInstance((char *)hashtable_search(names, edge));
		cA = (char *)mallocLocal(sizeof(char)*(2+strlen(cA2)+strlen(cA3)));
		sprintf(cA, "%s.%s", cA2, cA3);
		free(cA3);
	}

	return cA;
}

struct List *getOtherCap(struct PinchGraph *graph, struct List *adjacencyComponent, struct PinchEdge *edge, struct PinchEdge **edge2) {
	/*
	 * Gets the other cap on a sequence starting from the sequence given by the edge.
	 */
	struct List *list;

	list = constructEmptyList(0, NULL);

#ifdef BEN_DEBUG
	assert(listContains(adjacencyComponent, edge->to) == TRUE);
	assert((edge->segment->contig % 3) != 2);
	assert(edge->segment->start >= 1);
#endif

	while(TRUE) {
		*edge2 = getNextEdge(graph, edge);
		if(listContains(adjacencyComponent, (*edge2)->from)) {
			return list;
		}
		listAppend(list, *edge2);
		edge = *edge2;
	}
}

XMLNode getNode(XMLNode xMainNode, char *name) {
	/*
	 * Creates an empty node for a given name, if a node with the name already exists then the existing node is deleted.
	 * Can only be one or zero existing nodes with given name.
	 */

#ifdef BEN_DEBUG
	assert(xMainNode.nChildNode(name) <= 1);
#endif

	if(xMainNode.nChildNode(name) == 1) {
		XMLNode oldNode = xMainNode.getChildNode(name);
		oldNode.deleteNodeContent();
	}
	return xMainNode.addChild(name);
}

int32_t createReconstructionFilePath(const char *absolutePathPrefix, const char *relativeReconstructionPath, int32_t index, char **relativeFile) {
	/*
	 * Creates a nested path for a new reconstruction problem file. The new file is placed inside of a newly created directory contained with the directory
	 * pointed to by the absolute+ relative reconstruction path. The index is used to ensure a unique name.
	 */
	char *cA;
	int32_t i;
	cA = (char *)mallocLocal(sizeof(char)*(100+strlen(absolutePathPrefix)+strlen(relativeReconstructionPath)));
	sprintf(cA, "mkdir %s/%s/%i", absolutePathPrefix, relativeReconstructionPath, index);
	i = system(cA);
	if(i != 0) {
		logDebug("Tried to create a directory to hold the child reconstruction problem: " INT_STRING "\n", i);
		return i;
	}
	if(strlen(relativeReconstructionPath) > 0) {
		sprintf(cA, "%s/%i/reconstructionProblem.xml", relativeReconstructionPath, index);
	}
	else {
		sprintf(cA, "%i/reconstructionProblem.xml", index);
	}
	*relativeFile = cA;
	return 0;
}

struct hashtable *createValidMembersHash(struct PinchGraph *pinchGraph,
									     struct List *subChains, struct List *ends) {
	/*
	 * Create a hash to contain the valid members (atom instances and cap instances) that are contained in the considered reconstruction problem
	 * (defined by the input subchains and ends).
	 */
	int32_t k, l;
	struct hashtable *validMembers;
	struct Chain *chain;
	struct Segment *segment;
	struct PinchEdge *edge;

	logDebug("Creating the valid members for the reconstruction problem\n");

	validMembers = create_hashtable(subChains->length,
							hashtable_key, hashtable_equalKey,
							NULL, NULL);
	for(k=0; k<subChains->length; k++) {
		 chain = (struct Chain *)subChains->list[k];
		 while(chain != NULL) {
			 for(l=0; l<chain->segments->length; l++) {
				 segment = (struct Segment *)chain->segments->list[l];
				 edge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
#ifdef BEN_DEBUG
				 assert(edge != NULL);
				 assert(hashtable_search(validMembers, edge) == NULL);
				 assert(hashtable_search(validMembers, edge->rEdge) == NULL);
#endif
				 hashtable_insert(validMembers, edge, edge);
			 }
			 chain = chain->nLink;
		 }
	}
	for(k=0; k<ends->length; k++) {
		hashtable_insert(validMembers, ends->list[k], &k);
	}

	logDebug("Created the valid members for the considered reconstruction problem\n");

	return validMembers;
}

void createLeafEventsP(XMLNode cladeNode, struct hashtable *leafEvents) {
	int32_t i;

	if(cladeNode.nChildNode("clade") == 0) { //is leaf
		hashtable_insert(leafEvents, stringCopy(cladeNode.getAttribute("event")), &i);
	}
	else {
		for(i=0; i<cladeNode.nChildNode("clade"); i++) {
			createLeafEventsP(cladeNode.getChildNode("clade", i), leafEvents);
		}
	}
}

struct hashtable *createLeafEvents(XMLNode xMainNode) {
	struct hashtable *leafEvents;

	leafEvents = create_hashtable(LARGE_CHUNK_SIZE,
			hashtable_stringHashKey, hashtable_stringEqualKey,
			free, NULL);

	createLeafEventsP(xMainNode.getChildNode("event_tree").getChildNode("phylogeny").getChildNode("clade"), leafEvents);

	return leafEvents;
}

void createReconstructionProblem_updateCapAdjacenciesP(XMLNode cladeNode,
		struct PinchGraph *pinchGraph, struct hashtable *hash,
		struct hashtable *names, struct hashtable *validMembers) {
	int32_t i;
	struct PinchEdge *edge;
	struct PinchEdge *edge2;
	char *cA;

	if(cladeNode.nChildNode("clade") == 0) {
		edge = (struct PinchEdge *)hashtable_search(hash, (void *)cladeNode.getAttribute("instance"));
#ifdef BEN_DEBUG
		assert(edge != NULL);
#endif
		//add the adjacency
		edge2 = getNextEdge(pinchGraph, edge->rEdge);
		while(hashtable_search(validMembers, edge2) == NULL && hashtable_search(validMembers, edge2->rEdge) == NULL //ensure element is part of the problem at this level
				&& hashtable_search(validMembers, edge2->from) == NULL) {
			edge2 = getNextEdge(pinchGraph, edge2);
		}
		cA = getCapInstanceName(edge2, edge2->from, names);
		cladeNode.deleteAttribute("adjacency");
		cladeNode.addAttribute("adjacency", cA);
		free(cA);
	}
	else {
		for(i=0; i<cladeNode.nChildNode("clade"); i++) {
			createReconstructionProblem_updateCapAdjacenciesP(cladeNode.getChildNode("clade", i), pinchGraph, hash,
					names, validMembers);
		}
	}
}

void createReconstructionProblem_updateCapAdjacencies(XMLNode xMainNode, struct PinchGraph *pinchGraph,
													  struct List *subChains, struct List *ends, struct hashtable *names,
													  struct hashtable *validMembers) {
	/*
	 * Updates the adjacencies for the caps of a reconstruction problem.
	 */
	int32_t i;
	struct hashtable *hash;
	struct PinchEdge *edge;
	struct PinchVertex *vertex;
	char *cA;

	logDebug("Creating the cap adjacencies\n");

	hash = create_hashtable(subChains->length,
							hashtable_stringHashKey, hashtable_stringEqualKey,
							free, NULL);
	for(i=0; i<pinchGraph->vertices->length; i++) { //put the edges in a hash with there end instance name
		vertex = (struct PinchVertex *)pinchGraph->vertices->list[i];
		if(isADeadEnd(vertex) == FALSE) {
			void *blackEdgeIterator = getBlackEdgeIterator(vertex);
			while((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
				cA = getCapInstanceName(edge, edge->from, names);
				hashtable_insert(hash, cA, edge);
			}
			destructBlackEdgeIterator(blackEdgeIterator);
		}
	}

	XMLNode capsNode = xMainNode.getChildNode("caps");
	for(i=0; i<capsNode.nChildNode("cap"); i++) {
		XMLNode capNode = capsNode.getChildNode("cap", i);
		createReconstructionProblem_updateCapAdjacenciesP(capNode.getChildNode("phylogeny").getChildNode("clade"),
				pinchGraph, hash, names, validMembers);
	}
	logDebug("Updated the cap adjacencies\n");

	XMLNode atomsNode = xMainNode.getChildNode("atoms");
	for(i=0; i<atomsNode.nChildNode("atom"); i++) {
		XMLNode atomNode = atomsNode.getChildNode("atom", i);
#ifdef BEN_DEBUG
		assert(atomNode.nChildNode("cap") == 2);
#endif
		createReconstructionProblem_updateCapAdjacenciesP(atomNode.getChildNode("cap", 0).getChildNode("phylogeny").getChildNode("clade"),
														  pinchGraph, hash, names, validMembers);
		createReconstructionProblem_updateCapAdjacenciesP(atomNode.getChildNode("cap", 1).getChildNode("phylogeny").getChildNode("clade"),
																  pinchGraph, hash, names, validMembers);
	}
	logDebug("Updated the atom cap adjacencies\n");

	hashtable_destroy(hash, FALSE, TRUE);

}

void createReconstructionProblem_AddAtoms(XMLNode xMainNode, struct PinchGraph *pinchGraph, struct List *subChains,
		struct hashtable *names, struct List *contigIndexToContigStrings) {
	/*
	 * Adds the atoms to a reconstruction problem node.
	 */

	int32_t i, j;
	struct Chain *chain;
	struct Segment *segment;
	struct PinchEdge *edge;
	char *cA;
	static char scratch[MEDIUM_CHUNK_SIZE];
	static char scratch2[MEDIUM_CHUNK_SIZE];
	struct List *list;

	logDebug("Adding the atoms to the reconstruction problem\n");

	//create a list of atoms
	XMLNode atomsNode = getNode(xMainNode, "atoms");

	//for each atom in each chain create an atom node.
	for(i=0; i<subChains->length; i++) {
		chain = (struct Chain *)subChains->list[i];
		do {
#ifdef BEN_DEBUG
			assert(chain->segments->length > 0);
#endif
			segment = (struct Segment *)chain->segments->list[0];
			edge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
			if(isAStubOrCap(edge) == FALSE) {
				XMLNode atomNode = atomsNode.addChild("atom");
				segment = (struct Segment *)chain->segments->list[0];

				sprintf(scratch, "%i", segment->end - segment->start + 1);
				atomNode.addAttribute("length", scratch);
				sprintf(scratch, "%f", chain->score);
				atomNode.addAttribute("score", scratch);

				//get atom name
				cA = (char *)hashtable_search(names, edge);
#ifdef BEN_DEBUG
				assert(cA != NULL);
#endif
				cA = removeInstance(cA);
				list = constructEmptyList(0, NULL);
				if(cA[0] == '-') { //is reversed, we want in the positive orientation
					atomNode.updateText(cA+1);
					for(j=0; j<chain->segments->length; j++) {
						segment = (struct Segment *)chain->segments->list[j];
						listAppend(list, segment->rSegment);
					}
				}
				else {
					atomNode.updateText(cA);
					for(j=0; j<chain->segments->length; j++) {
						listAppend(list, chain->segments->list[j]);
					}
				}
				free(cA);

				XMLNode atomInstancesNode = atomNode.addChild("atom_instances");
				for(j=0; j<list->length; j++) {
					XMLNode atomInstanceNode = atomInstancesNode.addChild("atom_instance");
					//set the name
					segment = (struct Segment *)list->list[j];
					edge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
					cA = (char *)hashtable_search(names, edge);
#ifdef BEN_DEBUG
					assert(cA != NULL);
					assert(cA[0] != '-');
#endif
					atomInstanceNode.updateText(cA);
					//set the instance coordinates
					XMLNode coordinatesNode = atomInstanceNode.addChild("coordinates");
					coordinatesNode.addAttribute("contig", getXMLContigString(segment->contig, contigIndexToContigStrings));
					if(segment->start >= 1) {
						convertToXMLCoordinate(segment->start, scratch, scratch2);
					}
					else {
						convertToXMLCoordinate(segment->end, scratch, scratch2);
					}
					coordinatesNode.addAttribute("start", scratch);
					coordinatesNode.addAttribute("strand", scratch2);
				}
				//cleanup
				destructList(list);
			}
			else {
#ifdef BEN_DEBUG
				assert(chain->nLink == NULL || chain->pLink == NULL);
#endif
			}
			chain = chain->nLink;
		} while(chain != NULL);
	}

	logDebug("Created atoms\n");
}

void createReconstructionProblem_AddStrings(XMLNode xMainNode, struct PinchGraph *pinchGraph,
		struct List *ends, struct hashtable *names, struct List *contigIndexToContigStrings, struct hashtable *validMembers) {
	/*
	 * Adds strings to the reconstruction problem node, replacing any existing strings node.
	 */

	int32_t i, k;
	struct PinchEdge *edge;
	struct PinchEdge *edge2;
	struct PinchVertex *vertex;
	char *cA;
	static char scratch[MEDIUM_CHUNK_SIZE];
	static char scratch2[MEDIUM_CHUNK_SIZE];
	struct List *list;
	char *capInstanceName;

	logDebug("Creating strings\n");

	XMLNode stringsNode = getNode(xMainNode, "strings");

	//create the strings
	for(i=0; i<ends->length; i++) {
		vertex = (struct PinchVertex *)ends->list[i];
		//get the black edges for vertex, each of which is an instance of the cap.
		void *blackEdgeIterator = getBlackEdgeIterator(vertex);
		while((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
			//add its position attribute
			if(edge->segment->start <= 0) { //on the negative strand in the other direction, so generate the sequence also
				//now create a string to hold the atom instances
				XMLNode stringNode = stringsNode.addChild("string");
				stringNode.addAttribute("contig", getXMLContigString(edge->segment->contig, contigIndexToContigStrings));

				convertToXMLCoordinate(edge->rEdge->segment->end+1, scratch, scratch2);
				stringNode.addAttribute("start", scratch);

				capInstanceName = getCapInstanceName(edge, vertex, names);
				stringNode.addAttribute("left", capInstanceName);
				free(capInstanceName);

				//now do the rest of the string node
				list = getOtherCap(pinchGraph, ends, edge->rEdge, &edge2);
				capInstanceName = getCapInstanceName(edge2, edge2->from, names);
				stringNode.addAttribute("right", capInstanceName);
				free(capInstanceName);

				//add length attribute
				sprintf(scratch, "%i", edge2->segment->start - edge->rEdge->segment->end - 1);
				stringNode.addAttribute("length", scratch);

				using std::ostringstream; // stream insertion operators
				ostringstream cA2;

				cA2 << " ";
				for(k=0; k<list->length; k++) {
					edge = (struct PinchEdge *)list->list[k];
#ifdef BEN_DEBUG
					assert(edge != NULL);
#endif
					if(hashtable_search(validMembers, edge) != NULL || hashtable_search(validMembers, edge->rEdge) != NULL) {
						//check is edge is part of list
						cA = (char *)hashtable_search(names, edge);
#ifdef BEN_DEBUG
						assert(cA != NULL);
#endif
						cA2 << cA << " ";
					}
				}
				stringNode.updateText(cA2.str().c_str());
				destructList(list);
			}
		}
		destructBlackEdgeIterator(blackEdgeIterator);
	}

	logDebug("Created strings\n");
}

void createReconstructionProblem_AddSequences(XMLNode xMainNode, struct List *subChains,
		struct List *ends, struct List *contigIndexToContigStrings,
		struct hashtable *contigNamesToSequences) {
	/*
	 * Creates a description of the sequences for the given reconstruction problem.
	 */
	struct hashtable *hash;
	struct PinchVertex *vertex;
	struct PinchEdge *edge;
	char *cA;
	struct Sequence *sequence;
	int32_t i;
	static char scratch[MEDIUM_CHUNK_SIZE];

	logDebug("Creating sequences\n");

	XMLNode sequencesNode = getNode(xMainNode, "sequences");

	hash = create_hashtable(subChains->length,
							hashtable_stringHashKey, hashtable_stringEqualKey,
							NULL, NULL);

	//create the sequences
	for(i=0; i<ends->length; i++) {
		vertex = (struct PinchVertex *)ends->list[i];
		void *blackEdgeIterator = getBlackEdgeIterator(vertex);
		while((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
			cA = getXMLContigString(edge->segment->contig, contigIndexToContigStrings);
			sequence = (struct Sequence *)hashtable_search(contigNamesToSequences, cA);
#ifdef BEN_DEBUG
			assert(sequence != NULL);
			assert(strcmp(sequence->contig, cA) == 0);
#endif
			if(hashtable_search(hash, sequence->contig) == NULL) { //not yet created
				hashtable_insert(hash, sequence->contig, sequence->contig);

				XMLNode sequenceNode = sequencesNode.addChild("sequence");

				//Add contig name
				sequenceNode.addAttribute("contig", sequence->contig);

				//Add file path
				sequenceNode.addAttribute("sequence_file", sequence->file);

				//Add sequence length
				sprintf(scratch, "%i", sequence->length);
				sequenceNode.addAttribute("length", scratch);

				//Add event?
				sequenceNode.addAttribute("event", sequence->event);
			}
		}
		destructBlackEdgeIterator(blackEdgeIterator);
	}

	hashtable_destroy(hash, FALSE, FALSE);

	logDebug("Created sequences\n");
}

struct ReconstructionProblemGraphBits {
	struct List *subChains;
	struct List *adjacencyComponents;
	struct List *ends;
};

void destructReconstructionProblemGraphBits(struct ReconstructionProblemGraphBits *reconstructionProblemGraphBits) {
	//Destructs the basic lists, but not the stuff within them.
	reconstructionProblemGraphBits->subChains->destructElement = NULL;
	destructList(reconstructionProblemGraphBits->subChains);
	reconstructionProblemGraphBits->adjacencyComponents->destructElement = NULL;
	destructList(reconstructionProblemGraphBits->adjacencyComponents);
	reconstructionProblemGraphBits->ends->destructElement = NULL;
	destructList(reconstructionProblemGraphBits->ends);
	free(reconstructionProblemGraphBits);
}

int32_t createReconstructionProblem_AddComplexAdjacencyComponents(XMLNode xMainNode, struct PinchGraph *pinchGraph,
		struct List *subChains, struct List *adjacencyComponents, struct hashtable *names,
		const char *absolutePathPrefix, const char *relativeReconstructionPath, struct hashtable *reconstructionProblemGraphBitsHash) {
	/*
	 * Adds the complex (non chain) adjacency components to the graph.
	 */
	struct PinchEdge *edge;
	char *cA;
	int32_t i, j, l;
	static char scratch[MEDIUM_CHUNK_SIZE];
	struct Chain *chain;
	struct Segment *segment;
	struct List *list;
	struct ReconstructionProblemGraphBits *reconstructionProblemGraphBits;

	logDebug("Creating adjacency components\n");

	//for each adjacency component
	//create an adjacency component node.
	//create a file to contain the child problem associated with the adjacency component.
	//create the child node (by call to createSkeletalReconstructionProblem)
	//write the child node into the file and link to it in the parent

	//we only create a sub adjacency problem if the current problem is non-terminal.
	//a terminal problem contains no atoms and has only one effective adjacency component.
	j = FALSE;
	for(i=0; i<subChains->length; i++) {
		chain = (struct Chain *)subChains->list[i];
		segment = (struct Segment *)chain->segments->list[0];
		edge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
		if(chain->nLink != NULL || isAStubOrCap(edge) == FALSE) {
			j = TRUE;
			break;
		}
	}
	XMLNode adjacencyComponentsNode = getNode(xMainNode, "adjacency_components");
	if(j == TRUE || adjacencyComponents->length > 1) {
		for(i=0; i<adjacencyComponents->length; i++) {
			list = (struct List *)adjacencyComponents->list[i];
			XMLNode adjacencyComponentNode = adjacencyComponentsNode.addChild("adjacency_component");
			sprintf(scratch, "%i", i);
			adjacencyComponentNode.addAttribute("adjacency_component_index", scratch);
			j = createReconstructionFilePath(absolutePathPrefix, relativeReconstructionPath, i, &cA);
			if(j != 0) {
				return j;
			}
			adjacencyComponentNode.addAttribute("child_file", cA);

			reconstructionProblemGraphBits = (struct ReconstructionProblemGraphBits *)mallocLocal(sizeof(struct ReconstructionProblemGraphBits));
			reconstructionProblemGraphBits->subChains = constructEmptyList(0, NULL);
			reconstructionProblemGraphBits->adjacencyComponents = constructEmptyList(0, NULL);
			reconstructionProblemGraphBits->ends = listCopy(list);
			hashtable_insert(reconstructionProblemGraphBitsHash, cA, reconstructionProblemGraphBits);

			//add text string to adjacency component containing cap names
			using std::ostringstream; // stream insertion operators
			ostringstream cA2;
			for(l=0; l<list->length; l++) {
	#ifdef BEN_DEBUG
				assert(hashtable_search(names, list->list[l]) != NULL);
	#endif
				cA2 << (char *)hashtable_search(names, list->list[l]) << " ";
			}
			adjacencyComponentNode.updateText(cA2.str().c_str());
		}
	}

	logDebug("Created adjacency components\n");
	return 0;
}

int32_t createReconstructionProblem_AddChainsAndChainAdjacencyComponents(XMLNode xMainNode, struct PinchGraph *pinchGraph,
		struct List *subChains, struct List *adjacencyComponents, struct hashtable *names,
		const char *absolutePathPrefix, const char *relativeReconstructionPath, struct hashtable *reconstructionProblemGraphBitsHash) {
	/*
	 * Adds the chains and chain adjacency components to the reconstruction problem node.
	 */
	struct PinchEdge *edge;
	char *cA;
	int32_t i, j, l, k, m;
	static char scratch[MEDIUM_CHUNK_SIZE];
	struct PinchVertex *vertex;
	struct PinchVertex *vertex2;
	struct Chain *chain;
	//struct Chain *chain2;
	struct Segment *segment;
	struct List *list;
	struct ReconstructionProblemGraphBits * reconstructionProblemGraphBits;
	const char *rootEventName;

	logDebug("Creating chains and chain adjacency components\n");

	//create a chains node.
	XMLNode chainsNode = getNode(xMainNode, "chains");

	//for each link in a chain, create an adjacency component node and
	//a chain node, this requires making a list with the two caps of the chain link in.
	//again, create a file to contain the child problem.
	//create the child node  (by call to createSkeletalReconstructionProblem)
	//call createSkeletalReconstructionTree on the link in the chains subChains and adjacency components.
	j = adjacencyComponents->length;
	m = 0;
	rootEventName = xMainNode.getChildNode("event_tree").getChildNode("phylogeny").getChildNode("clade").getAttribute("event");
	XMLNode adjacencyComponentsNode = xMainNode.getChildNode("adjacency_components");
	for(i=0; i<subChains->length; i++) {
		chain = (struct Chain *)subChains->list[i];
		if(chain->nLink != NULL) {
			XMLNode chainNode = chainsNode.addChild("chain");

			//calculate the number of links in the chain
			/*k = 0;
			chain2 = chain;
			while(chain2->nLink != NULL) {
				k++;
				chain2 = chain2->nLink;
			}
			sprintf(scratch, "%i", k);
			//add length attribute
			chainNode.addAttribute("length", scratch);*/

			//add chain index attribute
			sprintf(scratch, "%i", m++);
			chainNode.addAttribute("chain_index", scratch);

			XMLNode configurationNode = chainNode.addChild("configuration");
			configurationNode.addAttribute("event", rootEventName); //all chains currently span to the root.
			XMLNode adjacencyPairsNode = configurationNode.addChild("adjacency_pairs");
			k = 0;
			while(chain->nLink != NULL) {
				//make the adjacency node.
				XMLNode adjacencyComponentNode = adjacencyComponentsNode.addChild("adjacency_component");
				sprintf(scratch, "%i", j);
				adjacencyComponentNode.addAttribute("adjacency_component_index", scratch);
				l = createReconstructionFilePath(absolutePathPrefix, relativeReconstructionPath, j, &cA);
				if(l != 0) {
					return l;
				}
				adjacencyComponentNode.addAttribute("child_file", cA); //cA is added to the hash below, so don't clean it up now.

				list = constructEmptyList(0, NULL);
				for(l=0; l<chain->ends->length; l++) {
					listAppend(list, chain->ends->list[l]);
				}
				//the two ends of the chain
				segment = (struct Segment *)chain->segments->list[0];
				edge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
				vertex = edge->to;
				if(listContains(list, vertex) == FALSE) {
					listAppend(list, vertex);
				}

				segment = (struct Segment *)chain->nLink->segments->list[0];
				edge = getContainingBlackEdge(pinchGraph, segment->contig, segment->start);
				vertex2 = edge->from;
				if(listContains(list, vertex2) == FALSE) {
					listAppend(list, vertex2);
				}

				reconstructionProblemGraphBits = (struct ReconstructionProblemGraphBits *)mallocLocal(sizeof(struct ReconstructionProblemGraphBits));
				reconstructionProblemGraphBits->subChains = listCopy(chain->subChains);
				reconstructionProblemGraphBits->adjacencyComponents = listCopy(chain->adjacencyComponents);
				reconstructionProblemGraphBits->ends = list;
				hashtable_insert(reconstructionProblemGraphBitsHash, cA, reconstructionProblemGraphBits);

				//add text string to adjacency component containing cap names
				using std::ostringstream; // stream insertion operators
				ostringstream cA2;
				for(l=0; l<list->length; l++) {
#ifdef BEN_DEBUG
					assert(hashtable_search(names, list->list[l]) != NULL);
#endif
					cA2 << (char *)hashtable_search(names, list->list[l]) << " ";
				}
				adjacencyComponentNode.updateText(cA2.str().c_str());

				//make the link node
				XMLNode adjacencyPairNode = adjacencyPairsNode.addChild("adjacency_pair");

				//index of the link in the chain
				//sprintf(scratch, "%i", k);
				//linkNode.addAttribute("link_index", scratch);

				//left side of chain link
				cA = (char *)hashtable_search(names, vertex);
#ifdef BEN_DEBUG
				assert(cA != NULL);
#endif
				adjacencyPairNode.addAttribute("left", cA);

				//right side of the chain link
				cA = (char *)hashtable_search(names, vertex2);
#ifdef BEN_DEBUG
				assert(cA != NULL);
#endif
				adjacencyPairNode.addAttribute("right", cA);

				sprintf(scratch, "%i", j);
				adjacencyPairNode.addAttribute("adjacency_component_index", scratch);

				chain = chain->nLink;
				j++;
				k++;
			}
		}
	}

	logDebug("Created chains and chain adjacency components\n");
	return 0;
}

int32_t addInternalInstancesToTree(const char *atomName, int32_t internalInstanceNumber,
										XMLNode atomInstancesNode, XMLNode cladeNode) {
	/*
	 * Adds internal instances (in an atom tree) to an atom XML.
	 */
	int32_t i;
	char *cA;

	if(cladeNode.nChildNode("clade") > 0) { //only for internal nodes.
		//call recursively
		for(i=0; i<cladeNode.nChildNode("clade"); i++) {
			internalInstanceNumber = addInternalInstancesToTree(atomName, internalInstanceNumber,
					atomInstancesNode, cladeNode.getChildNode("clade", i));
		}

		//add atom instance
		XMLNode atomInstanceNode = atomInstancesNode.addChild("atom_instance");
		cA = (char *)mallocLocal(sizeof(char)*(20+strlen(atomName)));
		sprintf(cA, "%s.%i", atomName, internalInstanceNumber++);
		atomInstanceNode.updateText(cA);
		free(cA);
		cladeNode.addAttribute("instance", atomInstanceNode.getText());
	}
	return internalInstanceNumber;
}

void addEventsToAtomInstancesP(XMLNode cladeNode, struct hashtable *hash) {
	int32_t i;

	hashtable_insert(hash, (void *)stringCopy(cladeNode.getAttribute("instance")),
			(void *)stringCopy(cladeNode.getAttribute("event")));
	for(i=0; i<cladeNode.nChildNode("clade"); i++) {
		addEventsToAtomInstancesP(cladeNode.getChildNode("clade", i), hash);
	}
}

void addEventsToAtomInstances(XMLNode atomNode) {
	/*
	 * Labels each atom instance with its corresponding event from the event tree.
	 */
	struct hashtable *hash;
	int32_t i;
	const char *cA;

	XMLNode atomInstancesNode = atomNode.getChildNode("atom_instances");

	hash = create_hashtable(atomInstancesNode.nChildNode("atom_instance")*2+1,
							hashtable_stringHashKey, hashtable_stringEqualKey,
							free, (void (*)(void *))freeXMLString);

	addEventsToAtomInstancesP(atomNode.getChildNode("phylogeny").getChildNode("clade"), hash);

	for(i=0; i<atomInstancesNode.nChildNode("atom_instance"); i++) {
		XMLNode atomInstanceNode = atomInstancesNode.getChildNode("atom_instance", i);
		cA = (char *)hashtable_search(hash, (void *)atomInstanceNode.getText());
#ifdef BEN_DEBUG
		assert(cA != NULL);
#endif
		atomInstanceNode.addAttribute("event", cA);
	}
	hashtable_destroy(hash, TRUE, TRUE);
}

void addBracketToAtomInstanceNamesInTree(XMLNode treeCopy, const char *name) {
	/*
	 * Used to convert the atom tree into a cap tree.
	 */
	int32_t i;
	const char *cA;
	char *cA2;
	char *cA3;

	cA = treeCopy.getAttribute("instance");
	if(cA != NULL) {
		cA2 = getInstance(cA);
		cA3 = (char *)mallocLocal(sizeof(char)*strlen(cA)+2);
		sprintf(cA3, "%s.%s", name, cA2);
		treeCopy.deleteAttribute("instance");
#ifdef BEN_DEBUG
		assert(treeCopy.getAttribute("instance") == NULL);
#endif
		treeCopy.addAttribute("instance", cA3);
		free(cA2);
		free(cA3);
	}

	for(i=0; i<treeCopy.nChildNode("clade"); i++) {
		addBracketToAtomInstanceNamesInTree(treeCopy.getChildNode("clade", i), name);
	}
}

void addAtomCap(XMLNode atomNode, const char *name) {
	/*
	 * Creates a cap in an atom using the contained 'atom' tree (produced by the subprogram).
	 */
	char *cA;

	XMLNode capNode = atomNode.addChild("cap");
	capNode.updateText(name);
	cA = atomNode.getChildNode("phylogeny").createXMLString();
	XMLNode capTreeNode = XMLNode::parseString(cA);
	freeXMLString(cA);
	addBracketToAtomInstanceNamesInTree(capTreeNode, name);
	capNode.addChild(capTreeNode);
}


char *getUniqueNamePrefix_uniqueNamePrefix = NULL;
int32_t getUniqueNamePrefix_counter = 0;

void setUniqueNamePrefix(const char *uniqueNamePrefix) {
	getUniqueNamePrefix_uniqueNamePrefix = stringCopy(uniqueNamePrefix);
}

char *getUniqueNamePrefix() {
	static char cA[100];

	assert(getUniqueNamePrefix_uniqueNamePrefix != NULL);
	sprintf(cA, "%s%i", getUniqueNamePrefix_uniqueNamePrefix, getUniqueNamePrefix_counter++);
	assert(strlen(cA) < 100);
	return cA;
}

int32_t createReconstructionProblem_RunExternalProgram(XMLNode *xMainNode, const char *absolutePathPrefix, const char *relativeReconstructionFilePath,
		const char *program, const char *tempFilePath) {
	/*
	 * Runs the adjacencies building program.
	 */
	char *cA;
	char *cA2;
	char *absoluteFilePath;
	char *uniqueNamePrefix;
	int32_t i;

	logDebug("Starting to call the external program: %s/%s\n", absolutePathPrefix, relativeReconstructionFilePath);

	//now writing out the file for the adjacencies building script
	absoluteFilePath = getAbsoluteFilePath(absolutePathPrefix, relativeReconstructionFilePath);
	logDebug("Writing the reconstruction problem to the file: %s\n", absoluteFilePath);
	writeXMLFile(absoluteFilePath, *xMainNode);

	//call the tree program with reconstruction program.
	cA = (char *)mallocLocal(sizeof(char)*(strlen(program) + strlen(absolutePathPrefix) + strlen(relativeReconstructionFilePath) +
										   + strlen(tempFilePath) + 300));
	i = constructRandomDir(tempFilePath, &cA2);
	if(i != 0) {
		logDebug("Tried to make a recursive directory of temp files but failed\n");
		return i;
	}
	uniqueNamePrefix = getUniqueNamePrefix(); //no need to free.
	sprintf(cA, "%s --absolutePathPrefix %s --reconstructionProblem %s --tempDirRoot %s --uniqueNamePrefix %s",
			program, absolutePathPrefix, relativeReconstructionFilePath, cA2, uniqueNamePrefix);
	i = system(cA);
	if(i != 0) {
		logDebug("Something went wrong calling the tree reconstruction program 1149\n");
		return i;
	}
	free(cA);
	i = destructRandomDir(cA2);
	if(i != 0) {
		logDebug("Tried to destroy a recursive directory of temp files but failed\n");
		return i;
	}

	free(absoluteFilePath);

	logDebug("Ran the external program script okay, apparently\n");
	return 0;
}

int32_t createReconstructionProblem_AddTrees(XMLNode *xMainNode, const char *absolutePathPrefix, const char *relativeReconstructionFilePath,
		const char *treeProgram, const char *tempFilePath,
		struct PinchGraph *pinchGraph, struct List *subChains, struct List *ends, struct hashtable *names, struct hashtable *validMembers) {
	/*
	 * In overview:
	 *
	 * run external code with current reconstruction problem + event dag, adding in atom trees.
	 * Then add internal instances of the atoms, and convert the atom tree into two contained caps.
	 */
	char *cA;
	char *absoluteFilePath;
	int32_t i, startTime;

	startTime = time(NULL);
	i = createReconstructionProblem_RunExternalProgram(xMainNode, absolutePathPrefix, relativeReconstructionFilePath,
													   treeProgram, tempFilePath);
	if(i != 0) {
		logDebug("Something went wrong trying to create the reconstruction trees (1173)\n");
		return i;
	}
	logInfo("Ran the treebuilding program in: %i seconds\n", time(NULL) - startTime);

	//read reconstruction node.
	absoluteFilePath = getAbsoluteFilePath(absolutePathPrefix, relativeReconstructionFilePath);
	*xMainNode=XMLNode::openFileHelper(absoluteFilePath, "reconstruction_problem");
	logDebug("Parsed the reconstruction problem from the XML file: %s\n", absoluteFilePath);

	//get the parent event tree and parent caps + atoms
	XMLNode eventTreeNode = (*xMainNode).getChildNode("event_tree");
	XMLNode atomsNode = (*xMainNode).getChildNode("atoms");

	//Now add the instances to each of the atom trees.
	//traverse through each atom tree and find internal nodes, adding instances
	for(i=0; i<atomsNode.nChildNode("atom"); i++) {
		XMLNode atomNode = atomsNode.getChildNode("atom", i);
		XMLNode cladeNode = atomNode.getChildNode("phylogeny").getChildNode("clade");
		XMLNode atomInstancesNode = atomNode.getChildNode("atom_instances");

#ifdef BEN_DEBUG
		assert(strcmp(cladeNode.getAttribute("event"), (*xMainNode).getChildNode("event_tree").getChildNode("phylogeny").getChildNode("clade").getAttribute("event")) == 0);
		assert(cladeNode.nChildNode("clade") == 1);
		int32_t k;
		for(j=0; j<atomInstancesNode.nChildNode("atom_instance"); j++) {
			XMLNode atomInstanceNode = atomInstancesNode.getChildNode("atom_instance", j);
			cA = getInstance(atomInstanceNode.getText());
			assert(sscanf(cA, "%i", &k) == 1);
			assert(k < atomInstancesNode.nChildNode("atom_instance"));
			free(cA);
		}
#endif
		//make the internal events
		addInternalInstancesToTree(atomNode.getText(), atomInstancesNode.nChildNode("atom_instance"),
								   atomInstancesNode, cladeNode);
		addEventsToAtomInstances(atomNode); //label the atom instances with the events

		//now modify the tree into two cap trees.
		cA = (char *)mallocLocal(sizeof(char)*(strlen(atomNode.getText()) + 2));
		sprintf(cA, "[%s", atomNode.getText());
		addAtomCap(atomNode, cA);
		sprintf(cA, "%s]", atomNode.getText());
		addAtomCap(atomNode, cA);
		free(cA);

		//now remove the old tree
		atomNode.getChildNode("phylogeny").deleteNodeContent();
	}

	writeXMLFile(absoluteFilePath, *xMainNode);

	//now we add in the adjacencies for the leaves of the cap
	createReconstructionProblem_updateCapAdjacencies(*xMainNode, pinchGraph,
													 subChains, ends, names, validMembers);

	free(absoluteFilePath);

	logDebug("Created trees\n");

	return 0;
}

void pruneEventTree_1(struct hashtable *hash, XMLNode cladeNode) {
	/*
	 * Puts events in tree into a hash.
	 */
	char *cA;
	int32_t i;

	cA = stringCopy(cladeNode.getAttribute("event"));
	hashtable_insert(hash, cA, cA);

	for(i=0; i<cladeNode.nChildNode("clade"); i++) {
		pruneEventTree_1(hash, cladeNode.getChildNode("clade", i));
	}
}

void pruneEventTree_2(struct hashtable *hash, XMLNode cladeNode, XMLNode parentCladeNode) {
	/*
	 * Removes the unary events in the tree which are not in the hash.
	 */
	int32_t i;

	if(hashtable_search(hash, (void *)cladeNode.getAttribute("event")) != NULL || cladeNode.nChildNode("clade") != 1) {
		for(i=0; i<cladeNode.nChildNode("clade"); i++) {
			pruneEventTree_2(hash, cladeNode.getChildNode("clade", i), cladeNode);
		}
	}
	else {
#ifdef BEN_DEBUG
		assert(cladeNode.nChildNode("clade") == 1);
		assert(strcmp(parentCladeNode.getAttribute("event"), cladeNode.getAttribute("event")) != 0);
#endif

		while(cladeNode.nChildNode() > 0) {
			//edit out this node
			parentCladeNode.addChild(cladeNode.getChildNode(0));
		}
		cladeNode.deleteNodeContent();

		//now recurse
		for(i=0; i<parentCladeNode.nChildNode("clade"); i++) {
			pruneEventTree_2(hash, parentCladeNode.getChildNode("clade", i), parentCladeNode);
		}
	}
}

XMLNode pruneEventTree(XMLNode eventTreeNode, XMLNode childCapsNode) {
	/*
	 * Copies the event tree and prunes out all the unary (duplication) nodes not referenced by
	 *  the caps of the child problem.
	 */
	int32_t i;
	struct hashtable *hash;
	char *cA;

	//Copies the node
	cA = eventTreeNode.createXMLString();
	eventTreeNode = XMLNode::parseString(cA);
	freeXMLString(cA);
	logDebug("Copied the parent tree\n");

	hash = create_hashtable((1+childCapsNode.nChildNode("cap"))*100,
							hashtable_stringHashKey, hashtable_stringEqualKey,
							free, NULL);
	//put child events in a hash
	for(i=0; i<childCapsNode.nChildNode("cap"); i++) {
#ifdef BEN_DEBUG
		assert(childCapsNode.getChildNode("cap", i).nChildNode("phylogeny") == 1);
#endif
		XMLNode childCapTreeNode = childCapsNode.getChildNode("cap", i).getChildNode("phylogeny");
#ifdef BEN_DEBUG
		assert(childCapTreeNode.nChildNode("clade") == 1);
#endif
		pruneEventTree_1(hash, childCapTreeNode.getChildNode("clade"));
	}

#ifdef BEN_DEBUG
		assert(hashtable_count(hash) >= 2);
#endif

	logDebug("Put the child events in a hash\n");

	//now prune the tree
	XMLNode topCladeNode = eventTreeNode.getChildNode("phylogeny").getChildNode("clade");
	pruneEventTree_2(hash, topCladeNode, topCladeNode); //top level node can not be pruned out

	logDebug("Pruned the parent event tree\n");

	//cleanup
	hashtable_destroy(hash, FALSE, TRUE);

	logDebug("Have finished pruning the child event tree\n");
	return eventTreeNode;
}

int32_t createReconstructionProblem_createSkeletalChildReconstructionProblems(XMLNode xMainNode,
								const char *absolutePathPrefix, const char *relativeReconstructionFilePath,
								struct hashtable *reconstructionProblemGraphBitsHash) {
	/*
	 * Creates the basic child reconstruction problem nodes for a given reconstruction problem node.
	 */
	int32_t i;
	char *cA;
	char *cA2;
	struct hashtable *hash;

	logDebug("Creating skeletal child reconstruction problems\n");

	XMLNode atomsNode = xMainNode.getChildNode("atoms");
	XMLNode capsNode = xMainNode.getChildNode("caps");
	XMLNode eventTreeNode = xMainNode.getChildNode("event_tree");

	//make hash of relevant cap+atom trees.
	hash = create_hashtable((capsNode.nChildNode("cap")+atomsNode.nChildNode("atom"))*2,
												 hashtable_stringHashKey, hashtable_stringEqualKey,
												 free, (void (*)(void *))freeXMLString);
#ifdef BEN_DEBUG
	int32_t k;
	k = 0;
#endif

	for(i=0; i<capsNode.nChildNode("cap"); i++) {
		XMLNode capNode = capsNode.getChildNode("cap", i);
#ifdef BEN_DEBUG
		assert(capNode.nChildNode("phylogeny") == 1);
		k++;
#endif
		hashtable_insert(hash, stringCopy(capNode.getText()), capNode.createXMLString());
	}
	logDebug("Pushed parent caps into the hash\n");

	for(i=0; i<atomsNode.nChildNode("atom"); i++) {
		XMLNode atomNode = atomsNode.getChildNode("atom", i);
#ifdef BEN_DEBUG
		k += 2;
#endif
		XMLNode capNode = atomNode.getChildNode("cap", 0);
		hashtable_insert(hash, stringCopy(capNode.getText()), capNode.createXMLString());
		capNode = atomNode.getChildNode("cap", 1);
		hashtable_insert(hash, stringCopy(capNode.getText()), capNode.createXMLString());
	}
	logDebug("Pushed the atom caps into the hash\n");

	//for each adjacency component, get the relevant cap trees + pruned event tree and put in the child reconstruction problem.
	//finally call the child reconstruction program.
	XMLNode adjacencyComponentsNode = xMainNode.getChildNode("adjacency_components");
	for(i=0; i<adjacencyComponentsNode.nChildNode("adjacency_component"); i++) {
		XMLNode adjacencyComponentNode = adjacencyComponentsNode.getChildNode("adjacency_component", i);
		//build the basic child file
		logDebug("Starting to create skeletal reconstruction problem node with caps: %s\n", adjacencyComponentNode.getText());

		//create the reconstruction problem node.
		XMLNode childReconstructionProblemNode = XMLNode::createXMLTopNode("reconstruction_problem",FALSE);

		childReconstructionProblemNode.addAttribute("parent_file", relativeReconstructionFilePath);

		//create the operations node
		childReconstructionProblemNode.addChild("operations");

		//create the caps node
		XMLNode childCapsNode = childReconstructionProblemNode.addChild("caps");
		logDebug("Created the basic reconstruction problem node\n");

		cA2 = stringCopy(adjacencyComponentNode.getText());
		cA = strtok(cA2, " ");
		while(cA != NULL) {
			cA = (char *)hashtable_search(hash, cA);
#ifdef BEN_DEBUG
			assert(cA != NULL);
			k--;
#endif
			childCapsNode.addChild(XMLNode::parseString(cA));
			cA = strtok(NULL, " ");
		}
		free(cA2);
		logDebug("Pushed the caps into the child\n");

		XMLNode childEventTreeNode = pruneEventTree(eventTreeNode, childCapsNode);
		childReconstructionProblemNode.addChild(childEventTreeNode);

		logDebug("Pushed the pruned event tree into the child\n");

		//write the problem to disk
		cA = getAbsoluteFilePath(absolutePathPrefix, adjacencyComponentNode.getAttribute("child_file"));
		writeXMLFile(cA, childReconstructionProblemNode);
		free(cA);
		logDebug("Written an updated child file: %s\n", adjacencyComponentNode.getAttribute("child_file"));
	}
#ifdef BEN_DEBUG
	if(adjacencyComponentsNode.nChildNode("adjacency_component") != 0) {
		assert(k == 0);
	}
#endif

	hashtable_destroy(hash, TRUE, TRUE);
	logDebug("Finished creating the child reconstruction problem nodes\n");

	return 0;
}

int32_t createReconstructionProblem(
		const char *absolutePathPrefix, const char *relativeReconstructionFilePath,
		struct List *subChains, struct List *adjacencyComponents, struct List *ends,
		struct PinchGraph *pinchGraph, struct hashtable *names,
		struct List *contigIndexToContigStrings, struct hashtable *contigNamesToSequences,
		const char *treeProgram, const char *tempFilePath) {
	/*
	 * The function fills out a reconstruction problem. The input problem must contain valid caps with trees etc.
	 */

	int32_t i, j;
	struct hashtable *validMembers;
	struct hashtable *reconstructionProblemGraphBitsHash;
	char *relativeReconstructionPath;
	char *absoluteFilePath;
	struct ReconstructionProblemGraphBits *reconstructionProblemGraphBits;

	absoluteFilePath = getAbsoluteFilePath(absolutePathPrefix, relativeReconstructionFilePath);
	XMLNode xMainNode=XMLNode::openFileHelper(absoluteFilePath);
	logDebug("Parsed the reconstruction problem from the XML file: %s\n", absoluteFilePath);

	//build the relative path to reconstruction problem file (this is used for building child reconstruction problems.
	relativeReconstructionPath = stringCopy(relativeReconstructionFilePath);
	i = strlen(relativeReconstructionPath);
	while(relativeReconstructionPath[i] != '/' && i > 0) {
		i--;
	}
#ifdef BEN_DEBUG
	assert(i >= 0);
#endif
	relativeReconstructionPath[i] = '\0';

	//build a hash used for getting the bits for building the child reconstruction problems from the graph
	reconstructionProblemGraphBitsHash = create_hashtable(adjacencyComponents->length*2 + SMALL_CHUNK_SIZE,
														  hashtable_stringHashKey, hashtable_stringEqualKey,
														  free, (void (*)(void *))destructReconstructionProblemGraphBits);

	logDebug("Setup preliminaries: %s %s\n", absolutePathPrefix, relativeReconstructionPath);

	//construct a hash of the 'valid members', that is the caps and atoms in this reconstruction problem node.
	validMembers = createValidMembersHash(pinchGraph, subChains, ends);
	//reset the cap adjacencies.
	createReconstructionProblem_updateCapAdjacencies(xMainNode, pinchGraph, subChains, ends, names, validMembers);
	//add the atoms to the node
	createReconstructionProblem_AddAtoms(xMainNode, pinchGraph, subChains, names, contigIndexToContigStrings);
	//add strings (replaces any existing)
	createReconstructionProblem_AddStrings(xMainNode, pinchGraph, ends, names, contigIndexToContigStrings, validMembers);
	//add sequences (replaces any existing)
	createReconstructionProblem_AddSequences(xMainNode, subChains, ends, contigIndexToContigStrings, contigNamesToSequences);
	//add the complex (non chain) adjacency components.
	i = createReconstructionProblem_AddComplexAdjacencyComponents(xMainNode, pinchGraph, subChains, adjacencyComponents,
			names, absolutePathPrefix, relativeReconstructionPath, reconstructionProblemGraphBitsHash);
	if(i != 0) {
		logInfo("Something went wrong constructing the adjacency components 1189\n");
		return i;
	}
	//now add the chains and chain adjacency components.
	i = createReconstructionProblem_AddChainsAndChainAdjacencyComponents(xMainNode, pinchGraph, subChains, adjacencyComponents,
			names, absolutePathPrefix, relativeReconstructionPath, reconstructionProblemGraphBitsHash);
	if(i != 0) {
		logInfo("Something went wrong constructing the chains and chain adjacency components 1196\n");
		return i;
	}

	//now add the trees to the problem.
	i = createReconstructionProblem_AddTrees(&xMainNode, absolutePathPrefix, relativeReconstructionFilePath,
			treeProgram, tempFilePath,
			pinchGraph, subChains, ends, names, validMembers);
	if(i != 0) {
		logInfo("Something went wrong constructing the trees 1202\n");
		return i;
	}

	//now create recursively the child reconstruction problems.
	i = createReconstructionProblem_createSkeletalChildReconstructionProblems(xMainNode,
									absolutePathPrefix, relativeReconstructionFilePath,
									reconstructionProblemGraphBitsHash);
	if(i != 0) {
		logInfo("Something went wrong constructing a child reconstruction problem 1209\n");
		return i;
	}

	logDebug("Writing file for main node (completed): %s\n", absoluteFilePath);
	writeXMLFile(absoluteFilePath, xMainNode);

	XMLNode adjacencyComponentsNode = xMainNode.getChildNode("adjacency_components");
	for(i=0; i<adjacencyComponentsNode.nChildNode("adjacency_component"); i++) {
		XMLNode adjacencyComponentNode = adjacencyComponentsNode.getChildNode("adjacency_component", i);

		//Now call the function to do everything recursively
		reconstructionProblemGraphBits = (struct ReconstructionProblemGraphBits *)hashtable_search(reconstructionProblemGraphBitsHash,
				(char *)adjacencyComponentNode.getAttribute("child_file"));
#ifdef BEN_DEBUG
		assert(reconstructionProblemGraphBits != NULL);
#endif
		j = createReconstructionProblem(absolutePathPrefix, adjacencyComponentNode.getAttribute("child_file"),
				reconstructionProblemGraphBits->subChains, reconstructionProblemGraphBits->adjacencyComponents, reconstructionProblemGraphBits->ends,
				pinchGraph, names, contigIndexToContigStrings, contigNamesToSequences,
				treeProgram, tempFilePath);
		if(j != 0) {
			logInfo("A problem creating the child problem file: %s 1331\n", adjacencyComponentNode.getAttribute("child_file"));
			return j;
		}
	}

	//cleanup
	free(absoluteFilePath);
	free(relativeReconstructionPath);
	hashtable_destroy(reconstructionProblemGraphBitsHash, TRUE, TRUE);
	hashtable_destroy(validMembers, FALSE, FALSE);

	logDebug("Done cleaning up, and have finished building the reconstruction tree!\n");

	return 0;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods to construct contig-event sets.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct hashtable *partitionSequencesByEvent(XMLNode sequences) {
	struct hashtable *eventsToSequences;
	struct List *list;
	int32_t i;
	char *cA;

	eventsToSequences = create_hashtable(0,
								hashtable_stringHashKey, hashtable_stringEqualKey,
								free, (void (*)(void *))destructList);

	for(i=0; i<sequences.nChildNode("sequence"); i++) {
		XMLNode sequence = sequences.getChildNode("sequence", i);
#ifdef BEN_DEBUG
		assert(sequence.getAttribute("contig") != NULL);
		assert(sequence.getAttribute("event") != NULL);
#endif
		cA = stringCopy(sequence.getAttribute("contig"));
		if(hashtable_search(eventsToSequences, (void *)sequence.getAttribute("event")) == NULL) {
			hashtable_insert(eventsToSequences, stringCopy(sequence.getAttribute("event")),
					constructEmptyList(0, free));
		}
		list = (struct List *)hashtable_search(eventsToSequences, (void *)sequence.getAttribute("event"));
#ifdef BEN_DEBUG
		assert(list != NULL);
#endif
		listAppend(list, cA);
	}
	return eventsToSequences;
}

struct List *constructContigEventSets_P(XMLNode cladeNode,
		struct List *contigEventSets,
		struct hashtable *eventsToSequences, float parentBranchLength) {
	/*
	 * Gets contigs in subtree of given clade and build a contig event set. Nodes are ordered by depth first discovery time.
	 */
	struct List *list;
	struct List *list2;
	struct ContigEventSet *contigEventSet;
	int32_t i, j;
	char *cA;

	contigEventSet = constructContigEventSet();
	contigEventSet->event = stringCopy(cladeNode.getAttribute("event"));
	if(cladeNode.getAttribute("branch_length") != NULL) {
		assert(sscanf(cladeNode.getAttribute("branch_length"), "%f", &contigEventSet->branchLength) == 1);
	}
	else {
		contigEventSet->branchLength = 0.0;
	}
	contigEventSet->branchLength += parentBranchLength;

	list = constructEmptyList(0, NULL);

	if(cladeNode.nChildNode("clade") == 0) { //is leaf
		listAppend(contigEventSets, contigEventSet);
		list2 = (struct List *)hashtable_search(eventsToSequences, (void *)cladeNode.getAttribute("event"));
		if(list2 != NULL) {
			for(i=0; i<list2->length; i++) {
				cA = stringCopy((char *)list2->list[i]);
#ifdef BEN_DEBUG
				assert(hashtable_search(contigEventSet->contigs, cA) == NULL);
#endif
				hashtable_insert(contigEventSet->contigs, cA, cA);
				listAppend(list, cA);
			}
		}
	}
	else if(cladeNode.nChildNode("clade") == 1) { //is unary node, so we factor it out.
		list2 = constructContigEventSets_P(cladeNode.getChildNode("clade"),
				contigEventSets, eventsToSequences, contigEventSet->branchLength);
		destructContigEventSet(contigEventSet);
		destructList(list);
		return list2;
	}
	else {
		listAppend(contigEventSets, contigEventSet);
		for(i=0; i<cladeNode.nChildNode("clade"); i++) {
			list2 = constructContigEventSets_P(cladeNode.getChildNode("clade", i),
					contigEventSets, eventsToSequences, 0.0);
			for(j=0; j<list2->length; j++) {
				cA = stringCopy((char *)list2->list[j]);
				listAppend(list, cA);
#ifdef BEN_DEBUG
				assert(hashtable_search(contigEventSet->contigs, cA) == NULL);
#endif
				hashtable_insert(contigEventSet->contigs, cA, cA);
			}
			destructList(list2);
		}
	}
	return list;
}

struct List *constructContigEventSets(XMLNode xMainNode) {
	/*
	 * Construct a contig-to-events map for each node in the event tree (excluding the root).
	 */
	struct List *contigEventSets;
	struct List *list;
	struct hashtable *eventsToSequences;

	//get preliminaries
	eventsToSequences = partitionSequencesByEvent(xMainNode.getChildNode("sequences"));
	XMLNode cladeNode = xMainNode.getChildNode("event_tree").getChildNode("phylogeny").getChildNode("clade").getChildNode("clade");
	contigEventSets = constructEmptyList(0, (void (*)(void *))destructContigEventSet);
	list = constructContigEventSets_P(cladeNode, contigEventSets, eventsToSequences, 0.0);
	destructList(list);

	//cleanup
	hashtable_destroy(eventsToSequences, TRUE, TRUE);

	return contigEventSets;
}
