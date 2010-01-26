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
#include "cactus_core.h"

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

CactusCoreInputParameters *constructCactusCoreInputParameters() {
	CactusCoreInputParameters *cCIP = (CactusCoreInputParameters *)malloc(sizeof(CactusCoreInputParameters));
	cCIP->extensionSteps = 3;
	cCIP->maxEdgeDegree = 50;
	cCIP->writeDebugFiles = 0;
	cCIP->minimumTreeCoverage = 0.5;
	cCIP->minimumTreeCoverageForAtoms = 0.9;
	cCIP->minimumAtomLength = 4;
	cCIP->minimumChainLength = 12;
	cCIP->trim = 3;
	cCIP->alignRepeats = 0;
	return cCIP;
}

void destructCactusCoreInputParameters(CactusCoreInputParameters *cCIP) {
	free(cCIP);
}

int32_t cactusCorePipeline(Net *net,
		CactusCoreInputParameters *cCIP,
		struct PairwiseAlignment *(*getNextAlignment)()) {
	struct PinchGraph *pinchGraph;
	struct CactusGraph *cactusGraph;
	int32_t i, startTime;
	struct PairwiseAlignment *pairwiseAlignment;
	struct List *threeEdgeConnectedComponents;
	struct List *chosenAtoms;
	struct List *list;
	struct List *biConnectedComponents;

	///////////////////////////////////////////////////////////////////////////
	//Setup the basic pinch graph
	///////////////////////////////////////////////////////////////////////////

	pinchGraph = constructPinchGraph(net);

	if(cCIP->writeDebugFiles) {
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
	pairwiseAlignment = getNextAlignment();
	logInfo("Now doing the pinch merges:\n");
	i = 0;

	struct FilterAlignmentParameters *filterParameters = (struct FilterAlignmentParameters *)malloc(sizeof(struct FilterAlignmentParameters));
	filterParameters->trim = cCIP->trim;
	filterParameters->alignRepeats = cCIP->alignRepeats;
	filterParameters->net = net;
	while(pairwiseAlignment != NULL) {
		logDebug("Alignment : %i , score %f\n", i++, pairwiseAlignment->score);
		logPairwiseAlignment(pairwiseAlignment);
		pinchMerge(pinchGraph, pairwiseAlignment,
				(void (*)(struct PinchGraph *pinchGraph, struct Segment *, struct Segment *, void *))filterSegmentAndThenAddToGraph, filterParameters);
		destructPairwiseAlignment(pairwiseAlignment);
		pairwiseAlignment = getNextAlignment();
	}
	free(filterParameters);
	logInfo("Finished pinch merges\n");

	if(cCIP->writeDebugFiles) {
		logDebug("Writing out dot formatted version of pinch graph with alignments added\n");
		writePinchGraph("pinchGraph2.dot", pinchGraph, NULL, NULL);
		logDebug("Finished writing out dot formatted version of pinch graph with alignments added\n");
	}

	checkPinchGraph(pinchGraph);
	logInfo("Pinched the graph in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// (3) Removing over aligned stuff.
	///////////////////////////////////////////////////////////////////////////

	//calculate the optimum undo threshold.
	//for steadily increasing extension threshold
	//calculate edges to ignore.
	//calculate edges to accept.. (does not include any ignored edges)..
	//calculate grey edge adjacency components, ignoring the ignored edges..
	//calculate sizes of adjacency components, record the largest.

	//choose the extension threshold which reduces the size of the

	startTime = time(NULL);
	assert(cCIP->maxEdgeDegree >= 1);
	logInfo("Before removing over aligned edges the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));
	removeOverAlignedEdges(pinchGraph, 0.0, cCIP->maxEdgeDegree, cCIP->extensionSteps, net);
	logInfo("After removing over aligned edges (degree %i) the graph has %i vertices and %i black edges\n", cCIP->maxEdgeDegree, pinchGraph->vertices->length, avl_count(pinchGraph->edges));
	removeOverAlignedEdges(pinchGraph, cCIP->minimumTreeCoverage, INT32_MAX, 0, net);
	logInfo("After removing atoms with less than the minimum tree coverage (%f) the graph has %i vertices and %i black edges\n", cCIP->minimumTreeCoverage, pinchGraph->vertices->length, avl_count(pinchGraph->edges));
	removeTrivialGreyEdgeComponents(pinchGraph, pinchGraph->vertices, net);
	logInfo("After removing the trivial graph components the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));
	checkPinchGraphDegree(pinchGraph, cCIP->maxEdgeDegree);

	if(cCIP->writeDebugFiles) {
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

	if(cCIP->writeDebugFiles) {
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
	computeCactusGraph(pinchGraph, &cactusGraph, &threeEdgeConnectedComponents);

	if(cCIP->writeDebugFiles) {
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

	if(cCIP->writeDebugFiles) {
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

	if(cCIP->writeDebugFiles) {
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
			net, cCIP->minimumTreeCoverageForAtoms, cCIP->minimumAtomLength, cCIP->minimumChainLength,
			pinchGraph);
	//now report the results
	logTheChosenAtomSubset(biConnectedComponents, chosenAtoms, pinchGraph, net);

	if(cCIP->writeDebugFiles) {
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
	//(15) Clean up.
	///////////////////////////////////////////////////////////////////////////

	//Destruct stuff
	startTime = time(NULL);
	destructList(biConnectedComponents);
	destructPinchGraph(pinchGraph);
	destructList(threeEdgeConnectedComponents);
	destructCactusGraph(cactusGraph);
	destructList(chosenAtoms);

	logInfo("Ran the core pipeline script in: %i seconds\n", time(NULL) - startTime);
	return 0;
}
