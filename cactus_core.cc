#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <time.h>

#include "xmlParser.h"
#include "Argument_helper.h"

extern "C" {
	#include "pinchGraph.h"
	#include "cactusGraph.h"
	#include "commonC.h"
	#include "fastCMaths.h"
	#include "bioioC.h"
	#include "hashTableC.h"
	#include "net.h"
	#include "pairwiseAlignment.h"
};

void writePinchGraph(char *name, struct PinchGraph *pinchGraph, const char *namePrefix,
						struct List *biConnectedComponents, struct List *adjacencyComponents,
						struct List *contigIndexToContigStrings) {
	struct hashtable *names;
	FILE *fileHandle;

	fileHandle = fopen(name, "w");
	names = getNames(pinchGraph, contigIndexToContigStrings, namePrefix);
	writeOutPinchGraphWithChains(pinchGraph, biConnectedComponents, adjacencyComponents, names, fileHandle);
	fclose(fileHandle);
	hashtable_destroy(names, TRUE, FALSE);
}

void writeCactusGraph(char *name, struct PinchGraph *pinchGraph,
						 struct CactusGraph *cactusGraph, const char *namePrefix,
						 struct List *contigIndexToContigStrings) {
	struct hashtable *names;
	FILE *fileHandle;

	fileHandle = fopen(name, "w");
	names = getNames(pinchGraph, contigIndexToContigStrings, namePrefix);
	writeOutCactusGraph(cactusGraph, pinchGraph, names, fileHandle);
	fclose(fileHandle);
	hashtable_destroy(names, TRUE, FALSE);
}

void writeNet(char *name, struct PinchGraph *pinchGraph,
				 struct Net *net, const char *namePrefix,
				 struct List *contigIndexToContigStrings) {
	struct hashtable *names;
	FILE *fileHandle;

	fileHandle = fopen(name, "w");
	names = getNames(pinchGraph, contigIndexToContigStrings, namePrefix);
	writeOutNet(net, names, fileHandle);
	fclose(fileHandle);
	hashtable_destroy(names, TRUE, FALSE);
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
	struct List *contigIndexToContigStrings;
	struct IntList *contigIndexToContigStart;
	struct hashtable *contigStringToContigIndex;
	struct List *threeEdgeConnectedComponents;
	struct Net *net;
	struct List *chosenAtoms;
	struct List *list;
	struct List *adjacencyComponents;
	struct hashtable *contigStringsToSequences;
	struct hashtable *names;
	struct List *contigEventSets;
	NetDisk *netDisk;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	std::string logLevelString = "None";
	std::string alignmentsFile = "None";
	std::string netDiskName = "None";
	std::string netName = "None";
	std::string tempFileRootDirectory = "None";
	std::string treeProgram = "None";
	std::string uniqueNamePrefix = "None";
	int32_t maxEdgeDegree = 50;
	bool writeDebugFiles = false;
	double proportionToKeep = 1.0;
	double discardRatio = 0.0;
	double minimumTreeCoverage = 0.8;
	int32_t minimumChainLength = 10;
	struct List *biConnectedComponents;

	dsr::Argument_helper ah;

	ah.new_named_string('a', "logLevel", "", "Set the log level", logLevelString);
	ah.new_named_string('b', "alignments", "", "The input alignments file", alignmentsFile);
	ah.new_named_string('c', "netDisk", "", "The location of the net disk", netDiskName);
	ah.new_named_string('d', "netName", "", "The name of the net to be aligned", netName);
	ah.new_named_string('e', "tempDirRoot", "", "The temp file root directory", tempFileRootDirectory);
	ah.new_named_int('f', "maxEdgeDegree", "", "Maximum degree of aligned edges", maxEdgeDegree);
	ah.new_flag('g', "writeDebugFiles", "Write the debug files", writeDebugFiles);
	ah.new_named_string('h', "treeProgram", "", "Program to add trees to a reconstruction problem", treeProgram);
	ah.new_named_string('i', "uniqueNamePrefix", "", "An alpha-numeric prefix which, when appended with any alpha-numeric characters is guaranteed to produce a unique name", uniqueNamePrefix);
	ah.new_named_double('j', "proportionToKeep", "", "The proportion of the highest scoring atoms to keep", proportionToKeep);
	ah.new_named_double('k', "discardRatio", "", "The proportion of the average atom score in an atom's chain an atom must score to be kept", discardRatio);
	ah.new_named_double('l', "minimumTreeCoverage", "", "Minimum tree coverage proportion to be included in the problem", minimumTreeCoverage);
	ah.new_named_int('m', "minimumChainLength", "", "The minimum chain length required to be included in the problem", minimumChainLength);

	ah.set_description("cactus");
	ah.set_author("Benedict Paten, benedict@soe.ucsc.edu");
	ah.set_version("0.1");

	ah.process(argc, argv);

	assert(alignmentsFile != "None");
	assert(absolutePathPrefix != "None");
	assert(relativeReconstructionProblemFile != "None");
	assert(tempFileRootDirectory != "None");

	std::string absoluteReconstructionProblemFile = (absolutePathPrefix + "/" + relativeReconstructionProblemFile).c_str();

	//////////////////////////////////////////////
	//Set up logging/log inputs
	//////////////////////////////////////////////

	if(strcmp(logLevelString.c_str(), "INFO") == 0) {
		setLogLevel(LOGGING_INFO);
	}
	if(strcmp(logLevelString.c_str(), "DEBUG") == 0) {
		setLogLevel(LOGGING_DEBUG);
	}

	logInfo("Pairwise alignments file : %s\n", alignmentsFile.c_str());
	logInfo("The reconstruction problem file (top of hierarchy) : %s\n", absoluteReconstructionProblemFile.c_str());
	logInfo("Root directory of temporary files: %s\n", tempFileRootDirectory.c_str());

	//////////////////////////////////////////////
	//Set up the temp file root directory
	//////////////////////////////////////////////

	initialiseTempFileTree((char *)tempFileRootDirectory.c_str(), 100, 4);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName.c_str());
	logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Parse the basic reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	net = netDisk_getNet(netName.c_str());
	logInfo("Parsed the net to be refined\n");

	startTime = time(NULL);
	setUniqueNamePrefix(uniqueNamePrefix.c_str());

	///////////////////////////////////////////////////////////////////////////
	//Setup the basic pinch graph
	///////////////////////////////////////////////////////////////////////////

	//Construct graph
	contigIndexToContigStrings = constructEmptyList(0, free);
	contigIndexToContigStart = constructEmptyIntList(0);

	pinchGraph = constructPinchGraph(net, contigIndexToContigStrings, contigIndexToContigStart);
	contigStringsToSequences = parseSequences(xMainNode);

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of initial pinch graph\n");
		writePinchGraph("pinchGraph1.dot", pinchGraph, uniqueNamePrefix.c_str(), NULL, NULL, contigIndexToContigStrings);
		logDebug("Finished writing out dot formatted version of initial pinch graph\n");
	}
	//check the graph is consistent
	checkPinchGraph(pinchGraph);

	//create a reverse hash of the contig strings.
	contigStringToContigIndex = create_hashtable(pinchGraph->vertices->length*2,
												hashtable_stringHashKey, hashtable_stringEqualKey,
												NULL, (void (*)(void *))destructList);
#ifdef BEN_DEBUG
	assert(contigIndexToContigStart->length == contigIndexToContigStrings->length);
#endif
	for(i=0; i<contigIndexToContigStrings->length; i++) {
#ifdef BEN_DEBUG
		assert(contigIndexToContigStrings->list[i] != NULL);
#endif
		list = (struct List *)hashtable_search(contigStringToContigIndex, contigIndexToContigStrings->list[i]);
		if(list == NULL) {
			list = constructEmptyList(0, (void (*)(void *))destructIntPair);
			hashtable_insert(contigStringToContigIndex, contigIndexToContigStrings->list[i], list);
		}
		listAppend(list, constructIntPair(contigIndexToContigStart->list[i], i));
	}

	logInfo("Constructed the graph in: %i seconds\n", time(NULL) - startTime);
	logInfo("Vertex number %i \n", pinchGraph->vertices->length);

	///////////////////////////////////////////////////////////////////////////
	//  (2) Adding alignments to the pinch graph
	///////////////////////////////////////////////////////////////////////////

	//Now run through all the alignments.
	startTime = time(NULL);
	fileHandle = fopen(alignmentsFile.c_str(), "r");
	pairwiseAlignment = cigarRead(fileHandle);
	logInfo("Now doing the pinch merges:\n");
	i = 0;
	while(pairwiseAlignment != NULL) {
		logDebug("Alignment : %i , score %f\n", i++, pairwiseAlignment->score);
		logPairwiseAlignment(pairwiseAlignment);
		pinchMerge(pinchGraph, pairwiseAlignment, contigStringToContigIndex);
		destructPairwiseAlignment(pairwiseAlignment);
		pairwiseAlignment = cigarRead(fileHandle);
	}
	logInfo("Finished pinch merges\n");

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of pinch graph with alignments added\n");
		writePinchGraph("pinchGraph2.dot", pinchGraph, uniqueNamePrefix.c_str(), NULL, NULL, contigIndexToContigStrings);
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
	removeOverAlignedEdges(pinchGraph, maxEdgeDegree);
	logInfo("After removing over aligned edges (degree %i) the graph has %i vertices and %i black edges\n", maxEdgeDegree, pinchGraph->vertices->length, avl_count(pinchGraph->edges));
	removeTrivialGreyEdgeComponents(pinchGraph, pinchGraph->vertices);
	logInfo("After removing the trivial graph components the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));
	checkPinchGraphDegree(pinchGraph, maxEdgeDegree);

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of pinch graph with over aligned edges removed\n");
		writePinchGraph("pinchGraph3.dot", pinchGraph, uniqueNamePrefix.c_str(), NULL, NULL, contigIndexToContigStrings);
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
		writePinchGraph("pinchGraph4.dot", pinchGraph, uniqueNamePrefix.c_str(), NULL, NULL, contigIndexToContigStrings);
		logDebug("Finished writing out dot formatted version of pinch graph with stub components linked to the sink vertex\n");
	}

	checkPinchGraph(pinchGraph);
	logInfo("Linked stub components to the sink component in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// (5) Constructing the basic cactus.
	///////////////////////////////////////////////////////////////////////////

	startTime = time(NULL);
	extraEdges = getEmptyExtraEdges(pinchGraph);
	i = computeCactusGraph(pinchGraph, &cactusGraph, &threeEdgeConnectedComponents, extraEdges, (char *)logLevelString.c_str());

	if(i != 0) {
		logInfo("Something went wrong constructing the initial cactus graph\n");
		return i;
	}

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of initial cactus graph\n");
		writeCactusGraph("cactusGraph1.dot", pinchGraph, cactusGraph, uniqueNamePrefix.c_str(),
						 contigIndexToContigStrings);
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
	i = computeCactusGraph(pinchGraph, &cactusGraph, &threeEdgeConnectedComponents, extraEdges, (char *)logLevelString.c_str());

	if(i != 0) {
		logInfo("Something went wrong constructing the cactus with circularised stems\n");
		return i;
	}

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of 2-edge component only cactus graph\n");
		writeCactusGraph("cactusGraph2.dot", pinchGraph, cactusGraph, uniqueNamePrefix.c_str(),
						 contigIndexToContigStrings);
		logDebug("Finished writing out dot formatted version of 2-edge component only cactus graph\n");
	}

	logInfo("Constructed the 2-edge component only cactus graph\n");

	checkCactusContainsOnly2EdgeConnectedComponents(cactusGraph);
	logInfo("Checked the cactus contains only 2-edge connected components in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// (7) Eliminating chain discontinuities.
	///////////////////////////////////////////////////////////////////////////

	startTime = time(NULL);
	breakLoopDiscontinuities(cactusGraph, extraEdges, threeEdgeConnectedComponents);
	destructCactusGraph(cactusGraph); //clean up the initial cactus graph.
	destructList(threeEdgeConnectedComponents);
	i = computeCactusGraph(pinchGraph, &cactusGraph, &threeEdgeConnectedComponents, extraEdges, (char *)logLevelString.c_str());

	if(i != 0) {
		logInfo("Something went wrong constructing the cactus without loop discontinuities\n");
		return i;
	}

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted version of the final cactus graph\n");
		writeCactusGraph("cactusGraph3.dot", pinchGraph, cactusGraph, uniqueNamePrefix.c_str(),
											contigIndexToContigStrings);
		logDebug("Finished writing out dot formatted version of the final cactus graph\n");
	}

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted final pinch graph showing chains prior to pruning\n");
		writePinchGraph("pinchGraph5.dot", pinchGraph, uniqueNamePrefix.c_str(), biConnectedComponents, NULL, contigIndexToContigStrings);
		logDebug("Finished writing out final pinch graph showing chains prior to pruning\n");
	}

	logInfo("Constructed the final cactus graph in: %i seconds\n", time(NULL) - startTime);

	////////////////////////////////////////////////
	//Get sorted bi-connected components.
	////////////////////////////////////////////////

	biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);

	///////////////////////////////////////////////////////////////////////////
	// (9) Choosing an atom subset.
	///////////////////////////////////////////////////////////////////////////

	//first get tree covering score for each atom -
	//drop all atoms with score less than X.
	//accept chains whose remaining element's combined length is greater than a set length.

	startTime = time(NULL);
	chosenAtoms = constructEmptyList(0, NULL);
	contigEventSets = constructContigEventSets(xMainNode);
	filterAtomsByTreeCoverageAndLength(biConnectedComponents, chosenAtoms,
			contigEventSets, proportionToKeep,
			discardRatio, minimumTreeCoverage, minimumChainLength,
			pinchGraph, contigIndexToContigStrings);
	filterAtomsByIfStubOrCap(biConnectedComponents, chosenAtoms, pinchGraph);
	//now report the results
	logTheChosenAtomSubset(biConnectedComponents, chosenAtoms, pinchGraph, contigEventSets, contigIndexToContigStrings);

	//cleanup the selection.
	destructList(contigEventSets);
	destructList(list);

	if(writeDebugFiles) {
		logDebug("Writing out dot formatted final pinch graph showing chains after pruning\n");
		writePinchGraph("pinchGraph5.dot", pinchGraph, uniqueNamePrefix.c_str(), net->chains, NULL, contigIndexToContigStrings);
		logDebug("Finished writing out final pinch graph showing chains prior to pruning\n");
	}

	///////////////////////////////////////////////////////////////////////////
	// (8) Constructing the net.
	///////////////////////////////////////////////////////////////////////////

	constructNetFromInputs(cactusGraph,
			pinchGraph, names, includeFn,
			contigStringsToSequences, contigIndexToContigStrings,
			netDisk, getUniqueName);

	///////////////////////////////////////////////////////////////////////////
	// (9) Write the net to disk.
	///////////////////////////////////////////////////////////////////////////

	netDisk_write(netDisk);

	logInfo("Updated the net on disk\m");

	///////////////////////////////////////////////////////////////////////////
	//(15) Clean up.
	///////////////////////////////////////////////////////////////////////////

	//Destruct stuff
	startTime = time(NULL);
	destructList(biConnectedComponents);
	destructPinchGraph(pinchGraph);
	destructList(extraEdges);
	destructList(contigIndexToContigStrings);
	destructIntList(contigIndexToContigStart);
	hashtable_destroy(contigStringToContigIndex, TRUE, FALSE);
	destructList(threeEdgeConnectedComponents);
	destructCactusGraph(cactusGraph);
	destructNet(net);
	destructList(chosenAtoms);
	hashtable_destroy(contigStringsToSequences, TRUE, FALSE);
	removeAllTempFiles();

	logInfo("Cleaned stuff up and am finished in: %i seconds\n", time(NULL) - startTime);
	return 0;
}
