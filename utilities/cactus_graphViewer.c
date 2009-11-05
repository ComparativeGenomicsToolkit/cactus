/*
 * The script builds a cactus graph representation of the chains and nets.
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

#include "utilitiesShared.h"

/*
 * Global variables.
 */
static double totalProblemSize;
static bool scaleNodeSizes = 1;
static bool nameLabels = 0;

static void usage() {
	fprintf(stderr, "cactus_graphViewer, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-e --outputFile : The file to write the dot graph file in.\n");
	fprintf(stderr, "-f --scaleNodeSizes : Scale the node sizes according to the volume of sequence they contained.\n");
	fprintf(stderr, "-g --nameLabels : Give chain and net nodes name labels.\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

void makeCactusTree_net(Net *net, FILE *fileHandle, const char *parentNodeName, const char *parentEdgeColour) {
	//Write the net nodes.
	char *netNameString = netMisc_nameToString(net_getName(net));
	const char *edgeColour = graphViz_getColour();
	addNodeToGraph(netNameString, fileHandle, calculateTotalContainedSequence(net)/totalProblemSize, "ellipse", netNameString);
	//Write in the parent edge.
	if(parentNodeName != NULL) {
		graphViz_addEdgeToGraph(parentNodeName, netNameString, fileHandle, "", parentEdgeColour, 10, 1, "forward");
	}
	//Create the chains.
	Net_ChainIterator *chainIterator = net_getChainIterator(net);
	Chain *chain;
	while((chain = net_getNextChain(chainIterator)) != NULL) {
		makeCactusTree_chain(chain, fileHandle, netNameString, edgeColour);
	}
	net_destructChainIterator(chainIterator);

	//Create the diamond node
	char *diamondNodeNameString = malloc(sizeof(char)*(strlen(netNameString) + 2));
	sprintf(diamondNodeNameString, "z%s", netNameString);
	const char *diamondEdgeColour = graphViz_getColour();
	//Create all the adjacency components linked to the diamond.
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	double size = 0.0; //get the size of the adjacency component organising node..
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		if(adjacencyComponent_getLink(adjacencyComponent) == NULL) { //linked to the diamond node.
			size += calculateTotalContainedSequence(adjacencyComponent_getNestedNet(adjacencyComponent));
		}
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);
	if(size > 0.0) {
		addNodeToGraph(diamondNodeNameString, fileHandle, size/totalProblemSize, "diamond", "");
		graphViz_addEdgeToGraph(netNameString, diamondNodeNameString, fileHandle, "", edgeColour, 10, 1, "forward");
		adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
		while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
			if(adjacencyComponent_getLink(adjacencyComponent) == NULL) { //linked to the diamond node.
				makeCactusTree_net(adjacencyComponent_getNestedNet(adjacencyComponent), fileHandle, diamondNodeNameString, diamondEdgeColour);
			}
		}
		net_destructAdjacencyComponentIterator(adjacencyComponentIterator);
	}
	free(netNameString);
	free(diamondNodeNameString);
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
			{ "scaleNodeSizes", no_argument, 0, 'f' },
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
				scaleNodeSizes = !scaleNodeSizes;
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

	assert(logLevelString == NULL || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
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

	totalProblemSize = calculateTotalContainedSequence(net);
	fileHandle = fopen(outputFile, "w");
	graphViz_setupGraphFile(fileHandle);
	makeCactusTree_net(net, fileHandle, NULL, NULL);
	graphViz_finishGraphFile(fileHandle);
	fclose(fileHandle);
	logInfo("Written the tree to file\n");

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	netDisk_destruct(netDisk);

	return 0;
}
