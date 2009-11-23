/*
 * The script builds a cactus tree representation of the chains and nets.
 * The format of the output graph is dot format.
 */
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"

/*
 * Global variables.
 */
static bool edgeColours = 1;
static bool nameLabels = 0;

static void usage() {
	fprintf(stderr, "cactus_graphViewer, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-e --outputFile : The file to write the dot graph file in.\n");
	fprintf(stderr, "-f --chainColours : Do not give chains distinct colours (instead of just black)\n");
	fprintf(stderr, "-g --nameLabels : Give chain and net nodes name labels.\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

void addEndNodeToGraph(End *end, FILE *fileHandle) {
	const char *nameString = netMisc_nameToStringStatic(end_getName(end));
	graphViz_addNodeToGraph(nameString, fileHandle, nameString, 0.5, 0.5, "circle", "black", 14);
}

void addEdgeToGraph(End *end1, End *end2, const char *colour, const char *label, double length, double weight, const char *direction, FILE *fileHandle) {
	char *nameString1 = netMisc_nameToString(end_getName(end1));
	char *nameString2 = netMisc_nameToString(end_getName(end2));
	graphViz_addEdgeToGraph(nameString1, nameString2, fileHandle, nameLabels ? label : "", colour, length, weight, direction);
	free(nameString1);
	free(nameString2);
}

void addAtomToGraph(Atom *atom, const char *colour, FILE *fileHandle) {
	static char label[100000];
	End *leftEnd = atom_getLeftEnd(atom);
	End *rightEnd = atom_getRightEnd(atom);
	addEndNodeToGraph(leftEnd, fileHandle);
	addEndNodeToGraph(rightEnd, fileHandle);
	Atom_InstanceIterator *iterator = atom_getInstanceIterator(atom);
	AtomInstance *atomInstance;
	while((atomInstance = atom_getNext(iterator)) != NULL) {
		atomInstance = atomInstance_getStrand(atomInstance) ? atomInstance : atomInstance_getReverse(atomInstance);
		if(atomInstance_getSequence(atomInstance) != NULL) {
			sprintf(label, "%s:%i:%i", netMisc_nameToStringStatic(sequence_getName(atomInstance_getSequence(atomInstance))),
					atomInstance_getStart(atomInstance), atomInstance_getStart(atomInstance)+atomInstance_getLength(atomInstance));
			addEdgeToGraph(endInstance_getEnd(atomInstance_get5End(atomInstance)),
						   endInstance_getEnd(atomInstance_get3End(atomInstance)),
						   edgeColours ? colour : "black", label, 1.5, 100, "forward", fileHandle);
		}
	}
	atom_destructInstanceIterator(iterator);
}

void addTrivialChainsToGraph(Net *net, FILE *fileHandle) {
	/*
	 * Add atoms not part of chain to the graph
	 */
	Net_AtomIterator *atomIterator = net_getAtomIterator(net);
	Atom *atom;
	while((atom = net_getNextAtom(atomIterator)) != NULL) {
		if(atom_getChain(atom) == NULL) {
			addAtomToGraph(atom, "black", fileHandle);
		}
	}
	net_destructAtomIterator(atomIterator);
}

void addChainsToGraph(Net *net, FILE *fileHandle) {
	/*
	 * Add atoms part of a chain to the graph.
	 */
	Net_ChainIterator *chainIterator = net_getChainIterator(net);
	Chain *chain;
	while((chain = net_getNextChain(chainIterator)) != NULL) {
		int32_t i, j;
		const char *chainColour;
		while((chainColour = graphViz_getColour(chainColour)) != NULL) { //ensure the chain colours don't match the trivial atom chains and the adjacencies.
			if(strcmp(chainColour, "black") != 0 && strcmp(chainColour, "grey") != 0) {
				break;
			}
		}
		Atom **atoms = chain_getAtomChain(chain, &i);
		for(j=0; j<i; j++) {
			addAtomToGraph(atoms[j], chainColour, fileHandle);
		}
		free(atoms);
	}
	net_destructChainIterator(chainIterator);
}

void addAdjacencies(Net *net, FILE *fileHandle) {
	/*
	 * Adds adjacency edges to the graph.
	 */
	static char label[10000];
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
		EndInstance *endInstance;
		char *netName = netMisc_nameToString(net_getName(net));
		while((endInstance = end_getNext(instanceIterator)) != NULL) {
			if(endInstance_getSequence(endInstance) != NULL) {
				endInstance = endInstance_getStrand(endInstance) ? endInstance : endInstance_getReverse(endInstance);
				EndInstance *endInstance2 = endInstance_getAdjacency(endInstance);
				if(!endInstance_getSide(endInstance)) {
					assert(endInstance_getCoordinate(endInstance) < endInstance_getCoordinate(endInstance2));
					sprintf(label, "%s:%i:%i:%s:%i", netMisc_nameToStringStatic(sequence_getName(endInstance_getSequence(endInstance))),
							endInstance_getCoordinate(endInstance), endInstance_getCoordinate(endInstance2),
							netName,
							net_getEndNumber(net));
					//sprintf(label, "%s:%i",
					//		netName,
					//		net_getEndNumber(net));
					addEdgeToGraph(endInstance_getEnd(endInstance), endInstance_getEnd(endInstance2), "grey", label, 1.5, 1, "forward", fileHandle);
				}
			}
		}
		free(netName);
		end_destructInstanceIterator(instanceIterator);
	}
	net_destructEndIterator(endIterator);
}

void addStubAndCapEndsToGraph(Net *net, FILE *fileHandle) {
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(!end_isAtomEnd(end)) {
			addEndNodeToGraph(end, fileHandle);
		}
	}
	net_destructEndIterator(endIterator);
}

void makeCactusGraph(Net *net, FILE *fileHandle) {
	if(net_getParentAdjacencyComponent(net) == NULL) {
		addStubAndCapEndsToGraph(net, fileHandle);
	}
	addTrivialChainsToGraph(net, fileHandle);
	addChainsToGraph(net, fileHandle);
	if(net_getAdjacencyComponentNumber(net) == 0) {
		addAdjacencies(net, fileHandle);
	}
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		makeCactusGraph(adjacencyComponent_getNestedNet(adjacencyComponent), fileHandle);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);
}

int main(int argc, char *argv[]) {
	NetDisk *netDisk;
	Net *net;
	FILE *fileHandle;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * netDiskName = NULL;
	char * netName = NULL;
	char * outputFile = NULL;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "netDisk", required_argument, 0, 'c' },
			{ "netName", required_argument, 0, 'd' },
			{ "outputFile", required_argument, 0, 'e' },
			{ "edgeColours", no_argument, 0, 'f' },
			{ "nameLabels", no_argument, 0, 'g' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:c:d:e:fgh", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
			case 'a':
				logLevelString = stringCopy(optarg);
				break;
			case 'c':
				netDiskName = stringCopy(optarg);
				break;
			case 'd':
				netName = stringCopy(optarg);
				break;
			case 'e':
				outputFile = stringCopy(optarg);
				break;
			case 'f':
				edgeColours = !edgeColours;
				break;
			case 'g':
				nameLabels = !nameLabels;
				break;
			case 'h':
				usage();
				return 0;
			default:
				usage();
				return 1;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// (0) Check the inputs.
	///////////////////////////////////////////////////////////////////////////

	assert(netDiskName != NULL);
	assert(netName != NULL);
	assert(outputFile != NULL);

	//////////////////////////////////////////////
	//Set up logging
	//////////////////////////////////////////////

	if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
		setLogLevel(LOGGING_INFO);
	}
	if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
		setLogLevel(LOGGING_DEBUG);
	}

	//////////////////////////////////////////////
	//Log (some of) the inputs
	//////////////////////////////////////////////

	logInfo("Net disk name : %s\n", netDiskName);
	logInfo("Net name : %s\n", netName);
	logInfo("Output graph file : %s\n", outputFile);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Parse the basic reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
	logInfo("Parsed the top level net of the cactus tree to build\n");

	///////////////////////////////////////////////////////////////////////////
	// Build the graph.
	///////////////////////////////////////////////////////////////////////////

	fileHandle = fopen(outputFile, "w");
	graphViz_setupGraphFile(fileHandle);
	makeCactusGraph(net, fileHandle);
	graphViz_finishGraphFile(fileHandle);
	fclose(fileHandle);
	logInfo("Written the tree to file\n");

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	netDisk_destruct(netDisk);

	return 0;
}
