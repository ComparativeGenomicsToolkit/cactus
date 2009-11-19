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
	fprintf(stderr, "-i --minimumTreeCoverage : Minimum tree coverage proportion of an atom to be included in the problem\n");
	fprintf(stderr, "-j --minimumAtomLength : The minimum length of an atom required to be included in the problem\n");
	fprintf(stderr, "-k --minimumChainLength : The minimum chain length required to be included in the problem\n");
	fprintf(stderr, "-l --trim : The length of bases to remove from the end of each alignment\n");
	fprintf(stderr, "-m --alignRepeats : Allow bases marked as repeats to be aligned (else alignments to these bases to be excluded)\n");
}

char *segment_getString(struct Segment *segment, Net *net) {
	Sequence *sequence = net_getSequence(net, segment->contig);
	if(segment->start >= 1) {
		return sequence_getString(sequence, segment->start, segment->end - segment->start + 1, 1);
	}
	else {
		return sequence_getString(sequence, -segment->end, segment->end - segment->start + 1, 0);
	}
}

bool containsRepeatBases(char *string) {
	/*
	 * Function returns non zero if the string contains lower case bases or a base of type 'N'
	 */
	int32_t i, j;
	j = strlen(string);
	for(i=0; i<j; i++) {
		char c = string[i];
		if(c != '-') {
			assert((c >= 65 && c <= 90) || (c >= 97 && c <= 122));
			if((c >= 97 && c <= 122) || c == 'N') {
				return 1;
			}
		}
	}
	return 0;
}

struct FilterAlignmentParameters {
	int32_t alignRepeats;
	int32_t trim;
	Net *net;
};

void filterSegmentAndThenAddToGraph(struct PinchGraph *pinchGraph, struct Segment *segment, struct Segment *segment2,
		struct FilterAlignmentParameters *filterParameters) {
	/*
	 * Function is used to filter the alignments added to the graph to optionally exclude alignments to repeats and to trim the edges of matches
	 * to avoid misalignments due to edge wander effects.
	 */
	assert(segment->end - segment->start == segment2->end - segment2->start);
	if(segment->end - segment->start + 1 > 2 * filterParameters->trim) { //only add to graph if non trivial in length.
		//Do the trim.
		segment->end -= filterParameters->trim;
		segment->start += filterParameters->trim;
		segment2->end -= filterParameters->trim;
		segment2->start += filterParameters->trim;
		assert(segment->end - segment->start == segment2->end - segment2->start);
		assert(segment->end - segment->start >= 0);

		//Now filter by repeat content.
		if(!filterParameters->alignRepeats) {
			char *string1 = segment_getString(segment, filterParameters->net);
			char *string2 = segment_getString(segment2, filterParameters->net);
			if(!containsRepeatBases(string1) && !containsRepeatBases(string2)) {
				pinchMergeSegment(pinchGraph, segment, segment2);
			}
			free(string1);
			free(string2);
		}
		else {
			pinchMergeSegment(pinchGraph, segment, segment2);
		}
	}
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
	float minimumTreeCoverage = 0.7;
	int32_t minimumAtomLength = 4;
	int32_t minimumChainLength = 12;
	int32_t trim = 3;
	int32_t alignRepeats = 0;

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
			{ "minimumTreeCoverage", required_argument, 0, 'i' },
			{ "minimumAtomLength", required_argument, 0, 'j' },
			{ "minimumChainLength", required_argument, 0, 'k' },
			{ "trim", required_argument, 0, 'l' },
			{ "alignRepeats", no_argument, 0, 'm' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		key = getopt_long(argc, argv, "a:b:c:d:e:f:ghi:j:k:l:m", long_options, &option_index);

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
			case 'i':
				assert(sscanf(optarg, "%f", &minimumTreeCoverage) == 1);
				break;
			case 'j':
				assert(sscanf(optarg, "%i", &minimumAtomLength) == 1);
				break;
			case 'k':
				assert(sscanf(optarg, "%i", &minimumChainLength) == 1);
				break;
			case 'l':
				assert(sscanf(optarg, "%i", &trim) == 1);
				break;
			case 'm':
				alignRepeats = !alignRepeats;
				break;
			default:
				usage();
				return 1;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// (0) Check the inputs.
	///////////////////////////////////////////////////////////////////////////

	assert(logLevelString == NULL || strcmp(logLevelString, "CRITICAL") == 0 || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
	assert(alignmentsFile != NULL);
	assert(netDiskName != NULL);
	assert(netName != NULL);
	assert(tempFileRootDirectory != NULL);
	assert(maxEdgeDegree > 0);
	assert(minimumTreeCoverage >= 0.0);
	assert(minimumAtomLength >= 0.0);
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

	struct FilterAlignmentParameters *filterParameters = malloc(sizeof(struct FilterAlignmentParameters));
	filterParameters->trim = trim;
	filterParameters->alignRepeats = alignRepeats;
	filterParameters->net = net;
	while(pairwiseAlignment != NULL) {
		logDebug("Alignment : %i , score %f\n", i++, pairwiseAlignment->score);
		logPairwiseAlignment(pairwiseAlignment);
		pinchMerge(pinchGraph, pairwiseAlignment,
				(void (*)(struct PinchGraph *pinchGraph, struct Segment *, struct Segment *, void *))filterSegmentAndThenAddToGraph, filterParameters);
		destructPairwiseAlignment(pairwiseAlignment);
		pairwiseAlignment = cigarRead(fileHandle);
	}
	free(filterParameters);
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
	i = computeCactusGraph(pinchGraph, &cactusGraph, &threeEdgeConnectedComponents, (char *)logLevelString);

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
	circulariseStems(cactusGraph);

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of 2-edge component only cactus graph\n");
		writeCactusGraph("cactusGraph2.dot", pinchGraph, cactusGraph);
		logDebug("Finished writing out dot formatted version of 2-edge component only cactus graph\n");
	}

	logInfo("Constructed the 2-edge component only cactus graph\n");

	checkCactusContainsOnly2EdgeConnectedComponents(cactusGraph);
	logInfo("Checked the cactus contains only 2-edge connected components in: %i seconds\n", time(NULL) - startTime);

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
			net, minimumTreeCoverage, minimumAtomLength, minimumChainLength,
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
	destructList(threeEdgeConnectedComponents);
	destructCactusGraph(cactusGraph);
	destructList(chosenAtoms);
	removeAllTempFiles();
	netDisk_destruct(netDisk);

	logInfo("Cleaned stuff up and am finished in: %i seconds\n", time(NULL) - startTime);
	return 0;
}
