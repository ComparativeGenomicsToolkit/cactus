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
						struct List *biConnectedComponents, struct List *groups) {
	FILE *fileHandle;
	fileHandle = fopen(name, "w");
	writeOutPinchGraphWithChains(pinchGraph, biConnectedComponents, groups, fileHandle);
	fclose(fileHandle);
}

void writeCactusGraph(char *name, struct PinchGraph *pinchGraph,
						 struct CactusGraph *cactusGraph) {
	FILE *fileHandle;
	fileHandle = fopen(name, "w");
	writeOutCactusGraph(cactusGraph, pinchGraph, fileHandle);
	fclose(fileHandle);
}

char *piece_getString(struct Piece *piece, Net *net) {
	Sequence *sequence = net_getSequence(net, piece->contig);
	if(piece->start >= 1) {
		return sequence_getString(sequence, piece->start, piece->end - piece->start + 1, 1);
	}
	else {
		return sequence_getString(sequence, -piece->end, piece->end - piece->start + 1, 0);
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

void filterPieceAndThenAddToGraph(struct PinchGraph *pinchGraph,
		struct Piece *piece, struct Piece *piece2,
		struct hashtable *vertexAdjacencyComponents,
		struct FilterAlignmentParameters *filterParameters) {
	/*
	 * Function is used to filter the alignments added to the graph to optionally exclude alignments to repeats and to trim the edges of matches
	 * to avoid misalignments due to edge wander effects.
	 */
	assert(piece->end - piece->start == piece2->end - piece2->start);
	if(piece->end - piece->start + 1 > 2 * filterParameters->trim) { //only add to graph if non trivial in length.
		//Do the trim.
		piece->end -= filterParameters->trim;
		piece->start += filterParameters->trim;
		piece2->end -= filterParameters->trim;
		piece2->start += filterParameters->trim;
		assert(piece->end - piece->start == piece2->end - piece2->start);
		assert(piece->end - piece->start >= 0);

		//Now filter by repeat content.
		if(!filterParameters->alignRepeats) {
			char *string1 = piece_getString(piece, filterParameters->net);
			char *string2 = piece_getString(piece2, filterParameters->net);
			if(!containsRepeatBases(string1) && !containsRepeatBases(string2)) {
				pinchMergePiece(pinchGraph, piece, piece2, vertexAdjacencyComponents);
			}
			free(string1);
			free(string2);
		}
		else {
			pinchMergePiece(pinchGraph, piece, piece2, vertexAdjacencyComponents);
		}
	}
}

CactusCoreInputParameters *constructCactusCoreInputParameters() {
	CactusCoreInputParameters *cCIP = (CactusCoreInputParameters *)malloc(sizeof(CactusCoreInputParameters));
	//Everything is essentially 'turned off' by default.
	cCIP->writeDebugFiles = 0;
	cCIP->alignRepeats = 0;
	cCIP->alignUndoLoops = 1;
	cCIP->maxEdgeDegree = 50;
	cCIP->extensionSteps = 0;
	cCIP->extensionStepsReduction = 0;
	cCIP->trim = 0;
	cCIP->trimReduction = 0;
	cCIP->minimumTreeCoverage = 0.0;
	cCIP->minimumTreeCoverageReduction = 0.0;
	cCIP->minimumBlockLength = 0;
	cCIP->minimumChainLength = 0;
	cCIP->minimumChainLengthReduction = 0;
	return cCIP;
}

void destructCactusCoreInputParameters(CactusCoreInputParameters *cCIP) {
	free(cCIP);
}

int32_t cactusCorePipeline(Net *net,
		CactusCoreInputParameters *cCIP,
		struct PairwiseAlignment *(*getNextAlignment)(),
		void (*startAlignmentStack)()) {
	struct PinchGraph *pinchGraph;
	struct PinchVertex *vertex;
	struct PinchEdge *pinchEdge;
	struct CactusGraph *cactusGraph;
	int32_t i, k, startTime;
	struct PairwiseAlignment *pairwiseAlignment;
	struct List *threeEdgeConnectedComponents;
	struct List *chosenBlocks;
	struct List *chosenPinchEdges;
	struct List *list;
	struct List *biConnectedComponents;
	struct hashtable *vertexAdjacencyComponents;

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
	//  Loop between adding and undoing pairwise alignments
	///////////////////////////////////////////////////////////////////////////

	/*
	 * These parameters are altered during the loops to push the sequences together.
	 */
	double minimumTreeCoverage = cCIP->minimumTreeCoverage;
	int32_t minimumChainLength = cCIP->minimumChainLength;
	int32_t trim = cCIP->trim;
	int32_t extensionSteps = cCIP->extensionSteps;

	vertexAdjacencyComponents = create_hashtable(pinchGraph->vertices->length*2, hashtable_intHashKey, hashtable_intEqualKey, NULL, free);

	//Build a hash putting the vertices all in the same adjacency component.
	for(i=0; i<pinchGraph->vertices->length; i++) {
		hashtable_insert(vertexAdjacencyComponents, pinchGraph->vertices->list[i], constructInt(0));
	}

	do {

#ifdef BEN_DEBUG
		///////////////////////////////////////////////////////////////////////////
		//  Check the adjacency vertex components.
		///////////////////////////////////////////////////////////////////////////

		assert(hashtable_count(vertexAdjacencyComponents) == pinchGraph->vertices->length);
		for(i=0; i<pinchGraph->vertices->length; i++) {
			vertex = pinchGraph->vertices->list[i];
			assert(hashtable_search(vertexAdjacencyComponents, vertex) != NULL);
		}
#endif

		///////////////////////////////////////////////////////////////////////////
		//  Adding alignments to the pinch graph
		///////////////////////////////////////////////////////////////////////////

		startTime = time(NULL);
		//Must be called to initialise the alignment stack..
		startAlignmentStack();

		//Now run through all the alignments.
		pairwiseAlignment = getNextAlignment();
		logInfo("Now doing the pinch merges:\n");
		i = 0;

		struct FilterAlignmentParameters *filterParameters = (struct FilterAlignmentParameters *)malloc(sizeof(struct FilterAlignmentParameters));
		assert(cCIP->trim >= 0);
		filterParameters->trim = trim;
		filterParameters->alignRepeats = cCIP->alignRepeats;
		filterParameters->net = net;
		while(pairwiseAlignment != NULL) {
			logDebug("Alignment : %i , score %f\n", i++, pairwiseAlignment->score);
			logPairwiseAlignment(pairwiseAlignment);
			pinchMerge(pinchGraph, pairwiseAlignment,
					(void (*)(struct PinchGraph *pinchGraph, struct Piece *, struct Piece *, struct hashtable *, void *))filterPieceAndThenAddToGraph,
					filterParameters, vertexAdjacencyComponents);
			pairwiseAlignment = getNextAlignment();
		}
		free(filterParameters);
		logInfo("Finished pinch merges\n");

#ifdef BEN_DEBUG
		for(i=0; i<pinchGraph->vertices->length; i++) {
			assert(hashtable_search(vertexAdjacencyComponents, pinchGraph->vertices->list[i]) != NULL);
		}
		assert(hashtable_count(vertexAdjacencyComponents) == pinchGraph->vertices->length);
#endif

		if(cCIP->writeDebugFiles) {
			logDebug("Writing out dot formatted version of pinch graph with alignments added\n");
			writePinchGraph("pinchGraph2.dot", pinchGraph, NULL, NULL);
			logDebug("Finished writing out dot formatted version of pinch graph with alignments added\n");
		}

		//Cleanup the adjacency component vertex hash.
		hashtable_destroy(vertexAdjacencyComponents, 1, 0);

		checkPinchGraph(pinchGraph);
		logInfo("Pinched the graph in: %i seconds\n", time(NULL) - startTime);

		///////////////////////////////////////////////////////////////////////////
		// Removing over aligned stuff.
		///////////////////////////////////////////////////////////////////////////

		startTime = time(NULL);
		assert(cCIP->maxEdgeDegree >= 1);
		logInfo("Before removing over aligned edges the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));
		assert(cCIP->extensionSteps >= 0);
		removeOverAlignedEdges(pinchGraph, 0.0, cCIP->maxEdgeDegree, NULL, extensionSteps, net);
		logInfo("After removing over aligned edges (degree %i) the graph has %i vertices and %i black edges\n", cCIP->maxEdgeDegree, pinchGraph->vertices->length, avl_count(pinchGraph->edges));
		//removeOverAlignedEdges(pinchGraph, cCIP->minimumTreeCoverage, INT32_MAX, 0, net);
		//logInfo("After removing blocks with less than the minimum tree coverage (%f) the graph has %i vertices and %i black edges\n", cCIP->minimumTreeCoverage, pinchGraph->vertices->length, avl_count(pinchGraph->edges));
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
		// Linking stub components to the sink component (if they haven't been already.
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
		// Constructing the basic cactus.
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
		// Circularising the stems in the cactus.
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
		// Get sorted bi-connected components.
		////////////////////////////////////////////////

		biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);

		if(cCIP->writeDebugFiles) {
			logDebug("Writing out dot formatted final pinch graph showing chains prior to pruning\n");
			writePinchGraph("pinchGraph5.dot", pinchGraph, biConnectedComponents, NULL);
			logDebug("Finished writing out final pinch graph showing chains prior to pruning\n");
		}

		///////////////////////////////////////////////////////////////////////////
		// Choosing a block subset.
		///////////////////////////////////////////////////////////////////////////

		//first get tree covering score for each block -
		//drop all blocks with score less than X.
		//accept chains whose remaining element's combined length is greater than a set length.

		startTime = time(NULL);
		chosenBlocks = filterBlocksByTreeCoverageAndLength(biConnectedComponents,
				net, minimumTreeCoverage, cCIP->minimumBlockLength, minimumChainLength,
				pinchGraph);
		//now report the results
		logTheChosenBlockSubset(biConnectedComponents, chosenBlocks, pinchGraph, net);

		chosenPinchEdges = constructEmptyList(0, NULL);
		for(i=0; i<chosenBlocks->length; i++) {
			pinchEdge = cactusEdgeToFirstPinchEdge(chosenBlocks->list[i], pinchGraph);
			if(!isAStubOrCap(pinchEdge)) {
				listAppend(chosenPinchEdges, pinchEdge);
			}
		}

		if(cCIP->writeDebugFiles) {
			logDebug("Writing out dot formatted final pinch graph showing chains after pruning\n");
			list = constructEmptyList(0, NULL);
			listAppend(list, chosenBlocks);
			writePinchGraph("pinchGraph6.dot", pinchGraph, list, NULL);
			destructList(list);
			logDebug("Finished writing out final pinch graph showing chains prior to pruning\n");
		}

		logInfo("I have chosen %i blocks which meet the requirements\n", chosenBlocks->length);


		if(cCIP->alignUndoLoops > 1) {
			///////////////////////////////////////////////////////////////////////////
			// Cleanup the cactus graph.
			///////////////////////////////////////////////////////////////////////////

			destructList(biConnectedComponents);
			destructList(threeEdgeConnectedComponents);
			destructCactusGraph(cactusGraph);
			destructList(chosenBlocks);

			///////////////////////////////////////////////////////////////////////////
			// Removing stuff not in the chosen list of blocks
			///////////////////////////////////////////////////////////////////////////

			if(1) {
				startTime = time(NULL);
				/*
				 * Make not chosen pinch edges.
				 */
				struct hashtable *chosenVerticesHash = create_hashtable(pinchGraph->vertices->length*2, hashtable_key, hashtable_equalKey, NULL, NULL);;
				for(i=0; i<chosenPinchEdges->length; i++) {
					pinchEdge = chosenPinchEdges->list[i];
					hashtable_insert(chosenVerticesHash, pinchEdge->from, pinchEdge->from);
					hashtable_insert(chosenVerticesHash, pinchEdge->to, pinchEdge->to);
				}
				struct List *notChosenPinchEdges = constructEmptyList(0, NULL);
				for(i=0; i<pinchGraph->vertices->length; i++) {
					vertex = pinchGraph->vertices->list[i];
					if(lengthBlackEdges(vertex) > 1 && !isAStubOrCap(getFirstBlackEdge(vertex)) &&
							hashtable_search(chosenVerticesHash, vertex) == NULL) {
						listAppend(notChosenPinchEdges, getFirstBlackEdge(vertex));
					}
					//if(hashtable_search(hash, vertex) == NULL && isAS
				}
				removeOverAlignedEdges(pinchGraph, 0.0, INT32_MAX, notChosenPinchEdges, 0, net);
				destructList(notChosenPinchEdges);
				hashtable_destroy(chosenVerticesHash, 0, 0);
				logInfo("After removing edges which were not chosen, the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));
				removeTrivialGreyEdgeComponents(pinchGraph, pinchGraph->vertices, net);
				logInfo("After removing the trivial graph components the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));
				checkPinchGraphDegree(pinchGraph, cCIP->maxEdgeDegree);

				if(cCIP->writeDebugFiles) {
					logDebug("Writing out dot formatted version of pinch graph with non chosen blocks removed.\n");
					writePinchGraph("pinchGraph7.dot", pinchGraph, NULL, NULL);
					logDebug("Finished writing out dot formatted version of pinch graph with with non chosen blocks removed.\n");
				}

				checkPinchGraph(pinchGraph);
				destructList(chosenPinchEdges);
				logInfo("Removed the non-chosen block edges in: %i seconds\n", time(NULL) - startTime);

				///////////////////////////////////////////////////////////////////////////
				// Calculate the adjacency components for the next loop.
				///////////////////////////////////////////////////////////////////////////

				//Build a hash putting each vertex in its own adjacency component.
				//Iterate through all the edges, keeping only those whose degree is greater than zero.
				vertexAdjacencyComponents = create_hashtable(pinchGraph->vertices->length*2, hashtable_intHashKey, hashtable_intEqualKey, NULL, free);
				struct List *groupsList = getRecursiveComponents(pinchGraph, NULL);
				for(i=0; i<groupsList->length; i++) {
					struct List *vertices = groupsList->list[i];
					for(k=0; k<vertices->length; k++) {
						hashtable_insert(vertexAdjacencyComponents, vertices->list[k], constructInt(i));
					}
				}
				destructList(groupsList);
			}
			else {
				///////////////////////////////////////////////////////////////////////////
				// Calculate the adjacency components for the next loop.
				///////////////////////////////////////////////////////////////////////////

				//Build a hash putting each vertex in its own adjacency component.
				//Iterate through all the edges, keeping only those whose degree is greater than zero.
				vertexAdjacencyComponents = create_hashtable(pinchGraph->vertices->length*2, hashtable_intHashKey, hashtable_intEqualKey, NULL, free);
				struct List *groupsList = getRecursiveComponents2(pinchGraph, chosenPinchEdges);
				for(i=0; i<groupsList->length; i++) {
					struct List *vertices = groupsList->list[i];
					for(k=0; k<vertices->length; k++) {
						hashtable_insert(vertexAdjacencyComponents, vertices->list[k], constructInt(i));
					}
				}
				destructList(groupsList);
				destructList(chosenPinchEdges);
			}

			///////////////////////////////////////////////////////////////////////////
			// Modify parameters for next loop
			///////////////////////////////////////////////////////////////////////////

			trim = trim - cCIP->trimReduction > 0 ? trim - cCIP->trimReduction : 0;
			extensionSteps = extensionSteps - cCIP->extensionStepsReduction > 0 ? extensionSteps - cCIP->extensionStepsReduction : 0;
			minimumTreeCoverage = minimumTreeCoverage - cCIP->minimumTreeCoverageReduction > 0 ? minimumTreeCoverage - cCIP->minimumTreeCoverageReduction : 0.0;
			minimumChainLength = minimumChainLength - cCIP->minimumChainLengthReduction > 0 ? minimumChainLength - cCIP->minimumChainLength : 0;


		}
		else {
			///////////////////////////////////////////////////////////////////////////
			// Constructing the net.
			///////////////////////////////////////////////////////////////////////////

			fillOutNetFromInputs(net, cactusGraph, pinchGraph, chosenBlocks);
		}
	} while(--cCIP->alignUndoLoops > 0);

	///////////////////////////////////////////////////////////////////////////
	//Clean up remaining stuff.
	///////////////////////////////////////////////////////////////////////////

	destructList(biConnectedComponents);
	destructList(threeEdgeConnectedComponents);
	destructCactusGraph(cactusGraph);
	destructList(chosenBlocks);
	destructList(chosenPinchEdges);
	destructPinchGraph(pinchGraph);

	logInfo("Ran the core pipeline script\n");
	return 0;
}
