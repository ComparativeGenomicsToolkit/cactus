#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "pinchGraph.h"
#include "cactusGraph.h"
#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "cactus.h"
#include "pairwiseAlignment.h"
#include "cactusNetFunctions.h"

void writePinchGraph(char *name, struct PinchGraph *pinchGraph,
						struct List *biConnectedComponents, struct List *adjacencyComponents) {
	FILE *fileHandle;
	fileHandle = fopen(name, "w");
	writeOutPinchGraphWithChains(pinchGraph, biConnectedComponents, adjacencyComponents, fileHandle);
	fclose(fileHandle);
}

void writeCactusGraph(char *name, struct PinchGraph *pinchGraph,
						 struct CactusGraph *cactusGraph) {
	FILE *fileHandle;
	fileHandle = fopen(name, "w");
	writeOutCactusGraph(cactusGraph, pinchGraph, fileHandle);
	fclose(fileHandle);
}

void usage() {
	fprintf(stderr, "cactus_core, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --alignments : The input alignments file\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-e --tempDirRoot : The temp file root directory\n");
	fprintf(stderr, "-f --maxEdgeDegree : Maximum degree of aligned edges\n");
	fprintf(stderr, "-g --writeDebugFiles : Write the debug files\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
	fprintf(stderr, "-k --proportionToKeep : The proportion of the highest scoring atoms to keep\n");
	fprintf(stderr, "-l --discardRatio : The proportion of the average atom score in an atom's chain an atom must score to be kept\n");
	fprintf(stderr, "-m --minimumTreeCoverage : Minimum tree coverage proportion to be included in the problem\n");
	fprintf(stderr, "-n --minimumChainLength : The minimum chain length required to be included in the problem\n");
}

int main(int argc, char *argv[]) {
	/*
	 * The script has a number of reasonably distinct stages (see the accompanying paper for a description)
	 *
	 * (1) Constructing the basic pinch graph.
	 *
	 * (2) Adding alignments to the pinch graph
	 *
	 * (3) Removing over aligned stuff.
	 *
	 * (4) Linking stub components to the sink component.
	 *
	 * (5) Constructing the basic cactus.
	 *
	 * (6) Circularising the stems in the cactus.
	 *
	 * (7) Eliminating chain discontinuities.
	 *
	 * (8) Constructing the net.
	 *
	 * (9) Choosing an atom subset.
	 *
	 * (10) Pruning the net.
	 *
	 * (11) Identifying pseudo adjacency components.
	 *
	 * (12) Adding the pseudo adjacency components to the net.
	 *
	 * (13) Adding the ends (caps/stubs) to the net
	 *
	 * (14) Constructing the recursion tree + build trees/events trees for each reconstruction problem.
	 *
	 * Additionally the script reads in the set of caps, inputs and alignments and
	 * outputs the recursion tree.
	 *
	 */
	struct PinchGraph *pinchGraph;
	struct CactusGraph *cactusGraph;
	struct List *extraEdges;
	int32_t i, startTime;
	FILE *fileHandle;
	struct PairwiseAlignment *pairwiseAlignment;
	struct List *threeEdgeConnectedComponents;
	struct List *chosenAtoms;
	struct List *list;
	NetDisk *netDisk;
	Net *net;
	struct List *biConnectedComponents;
	int key;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * alignmentsFile = NULL;
	char * netDiskName = NULL;
	char * netName = NULL;
	char * tempFileRootDirectory = NULL;
	int32_t maxEdgeDegree = 50;
	bool writeDebugFiles = 0;
	float proportionToKeep = 1.0;
	float discardRatio = 0.0;
	float minimumTreeCoverage = 0.8;
	int32_t minimumChainLength = 10;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "alignments", required_argument, 0, 'b' },
			{ "netDisk", required_argument, 0, 'c' },
			{ "netName", required_argument, 0, 'd' },
			{ "tempDirRoot", required_argument, 0, 'e' },
			{ "maxEdgeDegree", required_argument, 0, 'f' },
			{ "writeDebugFiles", no_argument, 0, 'g' },
			{ "help", no_argument, 0, 'h' },
			{ "proportionToKeep", required_argument, 0, 'k' },
			{ "discardRatio", required_argument, 0, 'l' },
			{ "minimumTreeCoverage", required_argument, 0, 'm' },
			{ "minimumChainLength", required_argument, 0, 'n' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		key = getopt_long(argc, argv, "a:b:c:d:e:f:ghk:l:m:n:", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
			case 'a':
				logLevelString = stringCopy(optarg);
				break;
			case 'b':
				alignmentsFile = stringCopy(optarg);
				break;
			case 'c':
				netDiskName = stringCopy(optarg);
				break;
			case 'd':
				netName = stringCopy(optarg);
				break;
			case 'e':
				tempFileRootDirectory = stringCopy(optarg);
				break;
			case 'f':
				assert(sscanf(optarg, "%i", &maxEdgeDegree) == 1);
				break;
			case 'g':
				writeDebugFiles = 1;
				break;
			case 'h':
				usage();
				return 0;
			case 'k':
				assert(sscanf(optarg, "%f", &proportionToKeep) == 1);
				break;
			case 'l':
				assert(sscanf(optarg, "%f", &discardRatio) == 1);
				break;
			case 'm':
				assert(sscanf(optarg, "%f", &minimumTreeCoverage) == 1);
				break;
			case 'n':
				assert(sscanf(optarg, "%i", &minimumChainLength) == 1);
				break;
			default:
				usage();
				return 1;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// (0) Check the inputs.
	///////////////////////////////////////////////////////////////////////////

	assert(logLevelString == NULL || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
	assert(alignmentsFile != NULL);
	assert(netDiskName != NULL);
	assert(netName != NULL);
	assert(tempFileRootDirectory != NULL);
	assert(maxEdgeDegree > 0);
	assert(proportionToKeep >= 0.0 && proportionToKeep <= 1.0);
	assert(discardRatio >= 0.0);
	assert(minimumTreeCoverage >= 0.0);
	assert(minimumChainLength >= 0);

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

	logInfo("Pairwise alignments file : %s\n", alignmentsFile);
	logInfo("Net disk name : %s\n", netDiskName);
	logInfo("Net name : %s\n", netName);
	logInfo("Temp file root directory : %s\n", tempFileRootDirectory);
	logInfo("Max edge degree : %i\n", maxEdgeDegree);

	//////////////////////////////////////////////
	//Set up the temp file root directory
	//////////////////////////////////////////////

	initialiseTempFileTree(tempFileRootDirectory, 100, 4);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Parse the basic reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
	logInfo("Parsed the net to be refined\n");

	startTime = time(NULL);

	///////////////////////////////////////////////////////////////////////////
	//Setup the basic pinch graph
	///////////////////////////////////////////////////////////////////////////

	pinchGraph = constructPinchGraph(net);

	if(writeDebugFiles) {
		writePinchGraph("pinchGraph1.dot", pinchGraph, NULL, NULL);
		logDebug("Finished writing out dot formatted version of initial pinch graph\n");
	}
	//check the graph is consistent
	checkPinchGraph(pinchGraph);

	logInfo("Constructed the graph in: %i seconds\n", time(NULL) - startTime);
	logInfo("Vertex number %i \n", pinchGraph->vertices->length);

	///////////////////////////////////////////////////////////////////////////
	//  (2) Adding alignments to the pinch graph
	///////////////////////////////////////////////////////////////////////////

	//Now run through all the alignments.
	startTime = time(NULL);
	fileHandle = fopen(alignmentsFile, "r");
	pairwiseAlignment = cigarRead(fileHandle);
	logInfo("Now doing the pinch merges:\n");
	i = 0;
	while(pairwiseAlignment != NULL) {
		logDebug("Alignment : %i , score %f\n", i++, pairwiseAlignment->score);
		logPairwiseAlignment(pairwiseAlignment);
		pinchMerge(pinchGraph, pairwiseAlignment);
		destructPairwiseAlignment(pairwiseAlignment);
		pairwiseAlignment = cigarRead(fileHandle);
	}
	logInfo("Finished pinch merges\n");

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of pinch graph with alignments added\n");
		writePinchGraph("pinchGraph2.dot", pinchGraph, NULL, NULL);
		logDebug("Finished writing out dot formatted version of pinch graph with alignments added\n");
	}

	checkPinchGraph(pinchGraph);
	logInfo("Pinched the graph in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// (3) Removing over aligned stuff.
	///////////////////////////////////////////////////////////////////////////

	startTime = time(NULL);
	assert(maxEdgeDegree >= 1);
	logInfo("Before removing over aligned edges the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));
	removeOverAlignedEdges(pinchGraph, maxEdgeDegree, net);
	logInfo("After removing over aligned edges (degree %i) the graph has %i vertices and %i black edges\n", maxEdgeDegree, pinchGraph->vertices->length, avl_count(pinchGraph->edges));
	removeTrivialGreyEdgeComponents(pinchGraph, pinchGraph->vertices);
	logInfo("After removing the trivial graph components the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));
	checkPinchGraphDegree(pinchGraph, maxEdgeDegree);

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of pinch graph with over aligned edges removed\n");
		writePinchGraph("pinchGraph3.dot", pinchGraph, NULL, NULL);
		logDebug("Finished writing out dot formatted version of pinch graph with over aligned edges removed\n");
	}

	checkPinchGraph(pinchGraph);
	logInfo("Removed the over aligned edges in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// (4) Linking stub components to the sink component.
	///////////////////////////////////////////////////////////////////////////

	startTime = time(NULL);
	linkStubComponentsToTheSinkComponent(pinchGraph);

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of pinch graph stub components linked to the sink vertex\n");
		writePinchGraph("pinchGraph4.dot", pinchGraph, NULL, NULL);
		logDebug("Finished writing out dot formatted version of pinch graph with stub components linked to the sink vertex\n");
	}

	checkPinchGraph(pinchGraph);
	logInfo("Linked stub components to the sink component in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// (5) Constructing the basic cactus.
	///////////////////////////////////////////////////////////////////////////

	startTime = time(NULL);
	extraEdges = getEmptyExtraEdges(pinchGraph);
	i = computeCactusGraph(pinchGraph, &cactusGraph, &threeEdgeConnectedComponents, extraEdges, (char *)logLevelString);

	if(i != 0) {
		logInfo("Something went wrong constructing the initial cactus graph\n");
		return i;
	}

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of initial cactus graph\n");
		writeCactusGraph("cactusGraph1.dot", pinchGraph, cactusGraph);
		logDebug("Finished writing out dot formatted version of initial cactus graph\n");
	}

	logInfo("Constructed the initial cactus graph in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// (6) Circularising the stems in the cactus.
	///////////////////////////////////////////////////////////////////////////

	startTime = time(NULL);
	circulariseStems(cactusGraph, extraEdges, threeEdgeConnectedComponents);
	destructCactusGraph(cactusGraph); //clean up the initial cactus graph.
	destructList(threeEdgeConnectedComponents);
	i = computeCactusGraph(pinchGraph, &cactusGraph, &threeEdgeConnectedComponents, extraEdges, (char *)logLevelString);

	if(i != 0) {
		logInfo("Something went wrong constructing the cactus with circularised stems\n");
		return i;
	}

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of 2-edge component only cactus graph\n");
		writeCactusGraph("cactusGraph2.dot", pinchGraph, cactusGraph);
		logDebug("Finished writing out dot formatted version of 2-edge component only cactus graph\n");
	}

	logInfo("Constructed the 2-edge component only cactus graph\n");

	checkCactusContainsOnly2EdgeConnectedComponents(cactusGraph);
	logInfo("Checked the cactus contains only 2-edge connected components in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// (7) Eliminating chain discontinuities.
	///////////////////////////////////////////////////////////////////////////

	/*startTime = time(NULL);
	breakLoopDiscontinuities(cactusGraph, extraEdges, threeEdgeConnectedComponents);
	destructCactusGraph(cactusGraph); //clean up the initial cactus graph.
	destructList(threeEdgeConnectedComponents);
	i = computeCactusGraph(pinchGraph, &cactusGraph, &threeEdgeConnectedComponents, extraEdges, (char *)logLevelString);

	if(i != 0) {
		logInfo("Something went wrong constructing the cactus without loop discontinuities\n");
		return i;
	}

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of the final cactus graph\n");
		writeCactusGraph("cactusGraph3.dot", pinchGraph, cactusGraph);
		logDebug("Finished writing out dot formatted version of the final cactus graph\n");
	}

	logInfo("Constructed the final cactus graph in: %i seconds\n", time(NULL) - startTime);*/

	////////////////////////////////////////////////
	//Get sorted bi-connected components.
	////////////////////////////////////////////////

	biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted final pinch graph showing chains prior to pruning\n");
		writePinchGraph("pinchGraph5.dot", pinchGraph, biConnectedComponents, NULL);
		logDebug("Finished writing out final pinch graph showing chains prior to pruning\n");
	}

	///////////////////////////////////////////////////////////////////////////
	// (9) Choosing an atom subset.
	///////////////////////////////////////////////////////////////////////////

	//first get tree covering score for each atom -
	//drop all atoms with score less than X.
	//accept chains whose remaining element's combined length is greater than a set length.

	startTime = time(NULL);
	chosenAtoms = filterAtomsByTreeCoverageAndLength(biConnectedComponents,
			net, proportionToKeep,
			discardRatio, minimumTreeCoverage, minimumChainLength,
			pinchGraph);
	//now report the results
	logTheChosenAtomSubset(biConnectedComponents, chosenAtoms, pinchGraph, net);

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted final pinch graph showing chains after pruning\n");
		list = constructEmptyList(0, NULL);
		listAppend(list, chosenAtoms);
		writePinchGraph("pinchGraph6.dot", pinchGraph, list, NULL);
		destructList(list);
		logDebug("Finished writing out final pinch graph showing chains prior to pruning\n");
	}

	///////////////////////////////////////////////////////////////////////////
	// (8) Constructing the net.
	///////////////////////////////////////////////////////////////////////////

	fillOutNetFromInputs(net, cactusGraph, pinchGraph, chosenAtoms);

	///////////////////////////////////////////////////////////////////////////
	// (9) Write the net to disk.
	///////////////////////////////////////////////////////////////////////////

	netDisk_write(netDisk);
	logInfo("Updated the net on disk\n");

	///////////////////////////////////////////////////////////////////////////
	//(15) Clean up.
	///////////////////////////////////////////////////////////////////////////

	//Destruct stuff
	startTime = time(NULL);
	destructList(biConnectedComponents);
	destructPinchGraph(pinchGraph);
	destructList(extraEdges);
	destructList(threeEdgeConnectedComponents);
	destructCactusGraph(cactusGraph);
	destructList(chosenAtoms);
	removeAllTempFiles();
	netDisk_destruct(netDisk);

	logInfo("Cleaned stuff up and am finished in: %i seconds\n", time(NULL) - startTime);
	return 0;
}
