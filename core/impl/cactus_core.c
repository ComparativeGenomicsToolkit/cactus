/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "pinchGraph.h"
#include "pinchGraphManipulation.h"
#include "cactusGraph.h"
#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "cactus.h"
#include "pairwiseAlignment.h"
#include "cactusFlowerFunctions.h"
#include "cactus_core.h"
#include "sonLib.h"
#include "adjacencyComponents.h"

void parseRequiredSpeciesTree(const char *string, CactusCoreInputParameters *cCIP) {
    cCIP->listOfSetsOfRequiredSpecies = stList_construct3(0, (void (*)(void *))stSortedSet_destruct);
    cCIP->listOfRequiredSpeciesCoverages = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    stTree *tree = stTree_parseNewickString(string);
    for (int32_t i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree *child = stTree_getChild(tree, i);
        int32_t j;
        assert(stTree_getChildNumber(child) > 0);
        assert(stTree_getLabel(stTree_getChild(child, 0)) != NULL);
        int32_t k =
                sscanf(stTree_getLabel(stTree_getChild(child, 0)), "%i", &j);
        assert(k == 1);
        assert(j >= 1);
        stIntTuple *coverage = stIntTuple_construct(1, j);
        stSortedSet *requiredSpeciesSet = stSortedSet_construct3(
                (int(*)(const void *, const void *)) strcmp, free);
        for (j = 1; j < stTree_getChildNumber(child); j++) {
            stTree *grandChild = stTree_getChild(child, j);
            assert(stTree_getLabel(grandChild) != NULL);
            assert(
                    stSortedSet_search(requiredSpeciesSet,
                            (void *)stTree_getLabel(grandChild)) == NULL);
            stSortedSet_insert(requiredSpeciesSet,
                    stString_copy(stTree_getLabel(grandChild)));
            assert(
                                stSortedSet_search(requiredSpeciesSet,
                                        (void *)stTree_getLabel(grandChild)) != NULL);
        }
        assert(
                stIntTuple_getPosition(coverage, 0) <= stSortedSet_size(
                        requiredSpeciesSet));
        stList_append(cCIP->listOfSetsOfRequiredSpecies, requiredSpeciesSet);
        stList_append(cCIP->listOfRequiredSpeciesCoverages, coverage);
    }
    stTree_destruct(tree);
    assert(stList_length(cCIP->listOfRequiredSpeciesCoverages) == stList_length(cCIP->listOfSetsOfRequiredSpecies));
    st_logDebug("Parsed %i required species sets\n", stList_length(cCIP->listOfRequiredSpeciesCoverages));
}

void writePinchGraph(char *name, struct PinchGraph *pinchGraph, struct List *biConnectedComponents, struct List *groups) {
    FILE *fileHandle = fopen(name, "w");
    struct hashtable *hash = createHashColouringPinchEdgesByChains(pinchGraph, biConnectedComponents);
    writeOutPinchGraphWithChains(pinchGraph, hash, groups, fileHandle);
    fclose(fileHandle);
    hashtable_destroy(hash, TRUE, FALSE);
}

void writeCactusGraph(char *name, struct PinchGraph *pinchGraph, struct CactusGraph *cactusGraph) {
    FILE *fileHandle = fopen(name, "w");
    writeOutCactusGraph(cactusGraph, pinchGraph, fileHandle);
    fclose(fileHandle);
}

char *piece_getString(struct Piece *piece, Flower *flower) {
    Sequence *sequence = flower_getSequence(flower, piece->contig);
    if (piece->start >= 1) {
        return sequence_getString(sequence, piece->start, piece->end - piece->start + 1, 1);
    } else {
        return sequence_getString(sequence, -piece->end, piece->end - piece->start + 1, 0);
    }
}

bool containsRepeatBases(char *string) {
    /*
     * Function returns non zero if the string contains lower case bases or a base of type 'N'
     */
    int32_t i, j;
    j = strlen(string);
    for (i = 0; i < j; i++) {
        char c = string[i];
        if (c != '-') {
            assert((c >= 65 && c <= 90) || (c >= 97 && c <= 122));
            if ((c >= 97 && c <= 122) || c == 'N') {
                return 1;
            }
        }
    }
    return 0;
}

struct FilterAlignmentParameters {
    int32_t alignRepeats;
    int32_t trim;
    Flower *flower;
};

void filterPieceAndThenAddToGraph(struct PinchGraph *pinchGraph, struct Piece *piece, struct Piece *piece2,
        stHash *vertexToAdjacencyComponent, struct FilterAlignmentParameters *filterParameters) {
    /*
     * Function is used to filter the alignments added to the graph to optionally exclude alignments to repeats and to trim the edges of matches
     * to avoid misalignments due to edge wander effects.
     */
    assert(piece->end - piece->start == piece2->end - piece2->start);
    if (piece->end - piece->start + 1 > 2 * filterParameters->trim) { //only add to graph if non trivial in length.
        //Do the trim.
        piece->end -= filterParameters->trim;
        piece->start += filterParameters->trim;
        piece2->end -= filterParameters->trim;
        piece2->start += filterParameters->trim;
#ifdef BEN_DEBUG
        assert(piece->end - piece->start == piece2->end - piece2->start);
        assert(piece->end - piece->start >= 0);
#endif

        //Now filter by repeat content.
        if (!filterParameters->alignRepeats) {
            char *string1 = piece_getString(piece, filterParameters->flower);
            char *string2 = piece_getString(piece2, filterParameters->flower);
            if (!containsRepeatBases(string1) && !containsRepeatBases(string2)) {
                pinchMergePiece(pinchGraph, piece, piece2, vertexToAdjacencyComponent);
            }
            free(string1);
            free(string2);
        } else {
            pinchMergePiece(pinchGraph, piece, piece2, vertexToAdjacencyComponent);
        }
    }
}

CactusCoreInputParameters *constructCactusCoreInputParameters() {
    CactusCoreInputParameters *cCIP = (CactusCoreInputParameters *) st_malloc(sizeof(CactusCoreInputParameters));
    //Everything is essentially 'turned off' by default.
    cCIP->writeDebugFiles = 0;

    cCIP->annealingRoundsLength = 1;
    cCIP->annealingRounds = st_malloc(sizeof(int32_t) * 1);
    cCIP->annealingRounds[0] = 0;
    cCIP->minimumChainLength = 0;
    cCIP->maximumAdjacencyComponentSize = INT64_MAX;

    cCIP->deannealingRoundsLength = 0;
    cCIP->deannealingRounds = st_malloc(0);

    cCIP->alignRepeatsAtRound = 0;

    cCIP->trim = st_malloc(0);
    cCIP->trimLength = 0;

    cCIP->minimumTreeCoverage = 0.0;
    cCIP->blockTrim = 0;
    cCIP->listOfSetsOfRequiredSpecies = NULL;
    cCIP->listOfRequiredSpeciesCoverages = NULL;
    cCIP->singleCopySpecies = NULL;
    cCIP->minimumDegree = 2;
    return cCIP;
}

void destructCactusCoreInputParameters(CactusCoreInputParameters *cCIP) {
    free(cCIP->annealingRounds);
    free(cCIP->deannealingRounds);
    free(cCIP->trim);
    stList_destruct(cCIP->listOfRequiredSpeciesCoverages);
    stList_destruct(cCIP->listOfSetsOfRequiredSpecies);
    if(cCIP->singleCopySpecies != NULL) {
        stSortedSet_destruct(cCIP->singleCopySpecies);
    }
    free(cCIP);
}

/*static bool computeCactusGraphP(struct PinchEdge *edge) {
 assert(edge != NULL);
 return 0;
 }*/

static struct CactusGraph *cactusCorePipeline_2(struct PinchGraph *pinchGraph, Flower *flower,
        bool(*passThroughEdgeFn)(struct PinchEdge *), int32_t attachEnds) {
    ///////////////////////////////////////////////////////////////////////////
    // Linking stub components to the sink component (if they haven't been already been).
    ///////////////////////////////////////////////////////////////////////////

    int32_t startTime = time(NULL);
    linkStubComponentsToTheSinkComponent(pinchGraph, flower, attachEnds);
    checkPinchGraph(pinchGraph);
    st_logInfo("Linked stub components to the sink component in: %i seconds\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Constructing the basic cactus.
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);
    struct CactusGraph *cactusGraph = computeCactusGraph(pinchGraph, passThroughEdgeFn);
    st_logInfo("Constructed the initial cactus graph in: %i seconds\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Circularising the stems in the cactus.
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);
    circulariseStems(cactusGraph, pinchGraph, flower);
    st_logInfo("Constructed the 2-edge component only cactus graph\n");
    checkCactusContainsOnly2EdgeConnectedComponents(cactusGraph);
    st_logInfo("Checked the cactus contains only 2-edge connected components in: %i seconds\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup.
    ///////////////////////////////////////////////////////////////////////////

    return cactusGraph;
}

struct List *getChosenBlockPinchEdges(stSortedSet *chosenBlocks, struct PinchGraph *pinchGraph) {
    struct CactusEdge *cactusEdge;
    struct List *chosenPinchEdges = constructEmptyList(0, NULL);
    stSortedSetIterator *it = stSortedSet_getIterator(chosenBlocks);
    while ((cactusEdge = stSortedSet_getNext(it)) != NULL) {
        struct PinchEdge *pinchEdge = cactusEdgeToFirstPinchEdge(cactusEdge, pinchGraph);
        if (!isAStub(pinchEdge)) {
            listAppend(chosenPinchEdges, pinchEdge);
        }
    }
    stSortedSet_destructIterator(it);
    return chosenPinchEdges;
}

struct CactusGraph *deanneal(Flower *flower, struct PinchGraph *pinchGraph, struct CactusGraph *cactusGraph,
        struct List **biConnectedComponents, int32_t minimumChainLengthInGraph, double minimumTreeCoverage,
        int32_t minimumBlockDegree, stList *listOfSetsOfRequiredSpecies, stList *listOfRequiredSpeciesCoverages,
        stSortedSet *singleCopySpecies) {
    ///////////////////////////////////////////////////////////////////////////
    // Choosing a block subset to undo.
    ///////////////////////////////////////////////////////////////////////////

    //Get all the blocks.
    stSortedSet *allBlocksOfDegree2OrHigher = filterBlocksByTreeCoverageAndLength(*biConnectedComponents, flower, 0.0,
            2, 0, 0, NULL, NULL, NULL, pinchGraph);
    //Get the blocks we want to keep
    stSortedSet *chosenBlocksToKeep = filterBlocksByTreeCoverageAndLength(*biConnectedComponents, flower,
            minimumTreeCoverage, minimumBlockDegree, 0, minimumChainLengthInGraph + 1, listOfSetsOfRequiredSpecies,
            listOfRequiredSpeciesCoverages,
            singleCopySpecies, pinchGraph);
    //Now get the blocks to undo by computing the difference.
    stSortedSet *blocksToUndo = stSortedSet_getDifference(allBlocksOfDegree2OrHigher, chosenBlocksToKeep);
    stSortedSet_destruct(chosenBlocksToKeep);
    stSortedSet_destruct(allBlocksOfDegree2OrHigher);

    //assert(stSortedSet_size(blocksToUndo) > 0);
    //now report the results
    //logTheChosenBlockSubset(biConnectedComponents, //We don't call this as it burns compute.
    //       blocksToUndo, pinchGraph, flower);
    st_logInfo("I have chosen %i blocks which meet the requirements to be undone\n", stSortedSet_size(blocksToUndo));

    ///////////////////////////////////////////////////////////////////////////
    // Undo the blocks.
    ///////////////////////////////////////////////////////////////////////////

    struct List *list = getChosenBlockPinchEdges(blocksToUndo, pinchGraph);
    removeOverAlignedEdges(pinchGraph, 0.0, INT32_MAX, list, 0, flower);
    destructList(list);
    st_logInfo("After removing edges which were not chosen, the graph has %i vertices and %i black edges\n",
            pinchGraph->vertices->length, avl_count(pinchGraph->edges));
    removeTrivialGreyEdgeComponents(pinchGraph, pinchGraph->vertices, flower);
    st_logInfo("After removing the trivial graph components the graph has %i vertices and %i black edges\n",
            pinchGraph->vertices->length, avl_count(pinchGraph->edges));
    stSortedSet_destruct(blocksToUndo);

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup the old cactus graph (it is now out of sync with the pinch graph,
    // after undoing the selected edges).
    ///////////////////////////////////////////////////////////////////////////

    destructList(*biConnectedComponents);
    destructCactusGraph(cactusGraph);

    ////////////////////////////////////////////////
    // Re-compute the cactus graph
    ////////////////////////////////////////////////

    cactusGraph = cactusCorePipeline_2(pinchGraph, flower,
            minimumBlockDegree == 0 ? doNotPassThroughDegree1EdgesFn : passThroughDegree1EdgesFn, 0);

    ////////////////////////////////////////////////
    // Get the sorted bi-connected components, again
    ////////////////////////////////////////////////

    *biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);

    return cactusGraph;
}

int32_t getMinimumChainLengthInGraph(struct List *biConnectedComponents, struct PinchGraph *pinchGraph) {
    /*
     * Gets the length of the smallest non-zero length chain in the graph (or INT32_MAX if the chain is empty).
     */
    int32_t minimumChainLengthInGraph = INT32_MAX;
    for (int32_t i = 0; i < biConnectedComponents->length; i++) {
        struct List *biConnectedComponent = biConnectedComponents->list[i];
        int32_t k = maxChainDegree(biConnectedComponent, pinchGraph);
        if (k > 1) {
            int32_t j = chainBaseLength(biConnectedComponent, pinchGraph);
            if (j >= 1 && j < minimumChainLengthInGraph) { //The greater than 1 is to avoid trying to undo chains consisting only of stubs or unaligned segments
                minimumChainLengthInGraph = j;
            }
        }
    }
    return minimumChainLengthInGraph;
}

void buildOutPinchGraph(struct PinchGraph *pinchGraph, stList *adjacencyComponents, Flower *flower,
        CactusCoreInputParameters *cCIP, struct PairwiseAlignment *(*getNextAlignment)(), void(*startAlignmentStack)(),
        int32_t minimumChainLength, int32_t trim, int32_t alignRepeats) {
    //struct PinchVertex *vertex;
    struct CactusGraph *cactusGraph;
    int32_t i, startTime;
    struct List *biConnectedComponents;
    struct PairwiseAlignment *pairwiseAlignment;

    ///////////////////////////////////////////////////////////////////////////
    //  Construct the extra adjacency components hash
    ///////////////////////////////////////////////////////////////////////////

    stHash *vertexToAdjacencyComponents = stHash_construct();
    int32_t maxAdjacencyComponentSize = 0;
    for (i = 0; i < stList_length(adjacencyComponents); i++) {
        stSortedSet *adjacencyComponent = stList_get(adjacencyComponents, i);
        stSortedSetIterator *it = stSortedSet_getIterator(adjacencyComponent);
        struct PinchVertex *vertex;
        while ((vertex = stSortedSet_getNext(it)) != NULL) {
            stHash_insert(vertexToAdjacencyComponents, vertex, adjacencyComponent);
        }
        stSortedSet_destructIterator(it);
        if (stSortedSet_size(adjacencyComponent) > maxAdjacencyComponentSize) {
            maxAdjacencyComponentSize = stSortedSet_size(adjacencyComponent);
        }
    }
    st_logInfo(
            "For min chain length %i we have %i adjacency components, the largest is %i vertices and the total vertices is %i\n",
            minimumChainLength, stList_length(adjacencyComponents), maxAdjacencyComponentSize,
            pinchGraph->vertices->length);

#ifdef BEN_DEBUG
    ///////////////////////////////////////////////////////////////////////////
    //  Check the adjacency vertex components.
    ///////////////////////////////////////////////////////////////////////////

    assert((int32_t)stHash_size(vertexToAdjacencyComponents) == pinchGraph->vertices->length);
    for (i = 0; i < pinchGraph->vertices->length; i++) {
        struct PinchVertex *vertex = pinchGraph->vertices->list[i];
        assert(stHash_search(vertexToAdjacencyComponents, vertex) != NULL);
    }
#endif

    ///////////////////////////////////////////////////////////////////////////
    //  Adding alignments to the pinch graph
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);
    //Must be called to initialise the alignment stack..
    startAlignmentStack();

    //Now run through all the alignments.
    pairwiseAlignment = getNextAlignment(); //we assume we own this memory, and will clean it up.
    st_logInfo("Now doing the pinch merges:\n");
    i = 0;

    struct FilterAlignmentParameters *filterParameters = (struct FilterAlignmentParameters *) st_malloc(
            sizeof(struct FilterAlignmentParameters));
    filterParameters->trim = trim;
    filterParameters->alignRepeats = alignRepeats; //cCIP->alignRepeats;
    filterParameters->flower = flower;

    while (pairwiseAlignment != NULL) {
        pinchMerge(
                pinchGraph,
                pairwiseAlignment,
                (void(*)(struct PinchGraph *pinchGraph, struct Piece *, struct Piece *, stHash *, void *)) filterPieceAndThenAddToGraph,
                filterParameters, vertexToAdjacencyComponents);
        destructPairwiseAlignment(pairwiseAlignment); //cleanup the previous alignment
        pairwiseAlignment = getNextAlignment();
    }
    free(filterParameters);
    st_logInfo("Finished pinch merges\n");

    ////////////////////////////////////////////////
    // Cleanup the residual bits of the graph.
    ////////////////////////////////////////////////

    //Cleanup the adjacency component vertex hash.
    stList_destruct(adjacencyComponents);
    stHash_destruct(vertexToAdjacencyComponents);

    checkPinchGraph(pinchGraph); //check the graph is all good.
    st_logInfo("Pinched the graph in: %i seconds\n", time(NULL) - startTime);

    removeTrivialGreyEdgeComponents(pinchGraph, pinchGraph->vertices, flower); //remove any pointless adjacencies.
    st_logInfo("After removing the trivial graph components the graph has %i vertices and %i black edges\n",
            pinchGraph->vertices->length, avl_count(pinchGraph->edges));
    checkPinchGraph(pinchGraph);

    ////////////////////////////////////////////////
    // Compute the cactus graph
    ////////////////////////////////////////////////

    cactusGraph = cactusCorePipeline_2(pinchGraph, flower,
            cCIP->minimumDegree <= 1 ? doNotPassThroughDegree1EdgesFn : passThroughDegree1EdgesFn, 0);

    ////////////////////////////////////////////////
    // Get sorted bi-connected components.
    ////////////////////////////////////////////////

    biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);

    ////////////////////////////////////////////////
    // Do the first deanneal of bad blocks, not worrying about minimum chain length.
    ////////////////////////////////////////////////

    if(cCIP->minimumTreeCoverage > 0.0 || cCIP->minimumDegree > 1 || cCIP->listOfRequiredSpeciesCoverages != NULL || cCIP->singleCopySpecies != NULL) {
        cactusGraph = deanneal(flower, pinchGraph, cactusGraph, &biConnectedComponents, 0,
                        cCIP->minimumTreeCoverage,
                        cCIP->minimumDegree, cCIP->listOfSetsOfRequiredSpecies, cCIP->listOfRequiredSpeciesCoverages, cCIP->singleCopySpecies);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Setup the deannealing rounds.
    ///////////////////////////////////////////////////////////////////////////

    int32_t minimumChainLengthInGraph = getMinimumChainLengthInGraph(biConnectedComponents, pinchGraph);
    assert(minimumChainLengthInGraph > 0);
    int32_t deannealingRound = 0;
    while (minimumChainLengthInGraph < minimumChainLength) {
        ///////////////////////////////////////////////////////////////////////////
        // Get the length of chains to remove in this deannealing round.
        ///////////////////////////////////////////////////////////////////////////

        int32_t minimumChainLengthToRemove = minimumChainLength - 1; //We will remove all chains less the minimum chain length
        if (deannealingRound < cCIP->deannealingRoundsLength) { //We have deannealing rounds to perform.
            if (cCIP->deannealingRounds[deannealingRound] < minimumChainLength) { //We can deanneal with a smaller value first
                minimumChainLengthToRemove = cCIP->deannealingRounds[deannealingRound++];
            }
        }

        //Start another loop if the minimum chain length in the graph is greater than the minimum chain length to remove.
        if (minimumChainLengthInGraph > minimumChainLengthToRemove) {
            continue;
        }

        ///////////////////////////////////////////////////////////////////////////
        // Do the actual deannealing of the blocks.
        ///////////////////////////////////////////////////////////////////////////

        cactusGraph = deanneal(flower, pinchGraph, cactusGraph, &biConnectedComponents, minimumChainLengthToRemove,
                0.0,
                0, NULL, NULL, NULL);

        ///////////////////////////////////////////////////////////////////////////
        // Recalculate the minimum length of chains in the graph
        ///////////////////////////////////////////////////////////////////////////

        minimumChainLengthInGraph = getMinimumChainLengthInGraph(biConnectedComponents, pinchGraph);

        st_logDebug(
                "The longest non-empty chain in the graph is %i bases, we removed chains less than or equal to %i bases and the required minimum length chain is %i bases\n",
                minimumChainLengthInGraph, minimumChainLengthToRemove, minimumChainLength);
    }

    ////////////////////////////////////////////////
    // Destruct the cactus graph for the loop
    ////////////////////////////////////////////////

    destructCactusGraph(cactusGraph);
    destructList(biConnectedComponents);

    ////////////////////////////////////////////////
    // Trim the edges of the pinch graph.
    ////////////////////////////////////////////////

    trimEdges(pinchGraph, cCIP->blockTrim, flower);
    st_logInfo("Trimmed %i from the end of edges\n", cCIP->blockTrim);
}

int64_t getMaximumAdjacencyComponentSize(struct PinchGraph *pinchGraph) {
    stList *adjacencyComponents = getAdjacencyComponents(pinchGraph);
    int64_t maxAdjacencyComponentSize = 0;
    for (int32_t i = 0; i < stList_length(adjacencyComponents); i++) {
        int64_t adjacencyComponentSize = getAdjacencyComponentSize(stList_get(adjacencyComponents, i));
        if (adjacencyComponentSize > maxAdjacencyComponentSize) {
            maxAdjacencyComponentSize = adjacencyComponentSize;
        }
    }
    stList_destruct(adjacencyComponents);
    return maxAdjacencyComponentSize;
}

int32_t cactusCorePipelineP(Flower *flower, CactusCoreInputParameters *cCIP,
        struct PairwiseAlignment *(*getNextAlignment)(), void(*startAlignmentStack)(),
        int32_t maxAnnealingRounds, bool recurse) {
    ///////////////////////////////////////////////////////////////////////////
    //Setup the basic pinch graph
    ///////////////////////////////////////////////////////////////////////////

    int32_t startTime = time(NULL);
    struct PinchGraph *pinchGraph = constructPinchGraph(flower);

    if (cCIP->writeDebugFiles) {
        writePinchGraph("pinchGraph1.dot", pinchGraph, NULL, NULL);
        st_logDebug("Finished writing out dot formatted version of initial pinch graph\n");
    }

    //check the graph is consistent
    checkPinchGraph(pinchGraph);

    st_logInfo("Constructed the graph in: %i seconds\n", time(NULL) - startTime);
    st_logInfo("Vertex number %i \n", pinchGraph->vertices->length);

    ///////////////////////////////////////////////////////////////////////////
    //  Loop between adding and undoing pairwise alignments
    ///////////////////////////////////////////////////////////////////////////

    /*
     * These parameters are altered during the loops to push/pull the sequences together/apart.
     */

    //Construct an initial adjacency component containing all the vertices
    stList *adjacencyComponents = stList_construct3(0, (void(*)(void *)) stSortedSet_destruct);
    stSortedSet *adjacencyComponent = stSortedSet_construct();
    stList_append(adjacencyComponents, adjacencyComponent);
    for (int32_t i = 0; i < pinchGraph->vertices->length; i++) {
        stSortedSet_insert(adjacencyComponent, pinchGraph->vertices->list[i]);
    }

    int64_t maxComponentSizes[cCIP->annealingRoundsLength]; // = st_malloc(sizeof(int64_t) * cCIP->annealingRoundsLength);

    for (int32_t annealingRound = 0; annealingRound < maxAnnealingRounds; annealingRound++) {
        /*
         * Need to record max adjacency component size.
         * If less than min and chain size is less than min then break.
         * Else record the max adjacency component size and continue up to max.
         * Repeat with best value.
         */

        int32_t trim = 0;
        bool alignRepeats = annealingRound >= cCIP->alignRepeatsAtRound;
        assert(annealingRound < cCIP->annealingRoundsLength);
        int32_t minimumChainLength = cCIP->annealingRounds[annealingRound];
        if (annealingRound < cCIP->trimLength) {
            trim = cCIP->trim[0];
            assert(trim >= 0);
        }
        st_logInfo("Starting annealing round %i, with minimum chain length %i, aiming at at overall minimum chain length of %i and maximum group size of %lli\n", annealingRound,
                minimumChainLength, cCIP->minimumChainLength, cCIP->maximumAdjacencyComponentSize);

        buildOutPinchGraph(pinchGraph, adjacencyComponents, flower, cCIP, getNextAlignment,
                startAlignmentStack, cCIP->annealingRounds[annealingRound], trim, alignRepeats);

        int64_t maxAdjacencyComponentSize = getMaximumAdjacencyComponentSize(pinchGraph);
        maxComponentSizes[annealingRound] = maxAdjacencyComponentSize;

        st_logInfo("The group size is %lli, the max group size allowed is %lli\n", maxAdjacencyComponentSize, cCIP->maximumAdjacencyComponentSize);

        ///////////////////////////////////////////////////////////////////////////
        // Un-link stub components from the sink component
        ///////////////////////////////////////////////////////////////////////////

        unlinkStubComponentsFromTheSinkComponent(pinchGraph, flower);

        if (minimumChainLength < cCIP->minimumChainLength || (maxAdjacencyComponentSize
                > cCIP->maximumAdjacencyComponentSize && annealingRound + 1 < maxAnnealingRounds)) { //We will loop around again.
            assert(annealingRound+1 < maxAnnealingRounds);
            st_logDebug("We will loop again\n");

            ///////////////////////////////////////////////////////////////////////////
            // Calculate the adjacency components for the next loop.
            ///////////////////////////////////////////////////////////////////////////

            adjacencyComponents = getAdjacencyComponents(pinchGraph);

        } else {
            if(recurse && maxAdjacencyComponentSize
                > cCIP->maximumAdjacencyComponentSize) {
                st_logInfo("The maximum size group is too big %lli, should be: %i\n", maxAdjacencyComponentSize, cCIP->maximumAdjacencyComponentSize);
                int32_t j=annealingRound;
                for(int32_t i=0; i<cCIP->annealingRoundsLength; i++) {
                    if(cCIP->annealingRounds[i] >= cCIP->minimumChainLength) {
                        if(maxComponentSizes[i] < maxAdjacencyComponentSize) {
                            maxAdjacencyComponentSize = maxComponentSizes[i];
                            j=i;
                        }
                    }
                }
                assert(j <= annealingRound);
                assert(j >= 0);
                if(j != annealingRound) {
                    /*
                     * Recurse to get back to the cCIP->annealingRounds[j] as the minimum chain length
                     */
                    st_logInfo("Recursing to get the optimal annealing round of %i with group size of %lli\n", j, maxAdjacencyComponentSize);
                    destructPinchGraph(pinchGraph);
                    assert(j < annealingRound);
                    return cactusCorePipelineP(flower, cCIP,
                            getNextAlignment, startAlignmentStack, j+1, 0);
                }
            }
            st_logDebug("We have finished iterating and will now fill out the net.\n");

            ////////////////////////////////////////////////
            // Recompute the cactus graph
            ////////////////////////////////////////////////

            struct CactusGraph *cactusGraph = cactusCorePipeline_2(pinchGraph, flower,
                    cCIP->minimumDegree <= 1 ? doNotPassThroughDegree1EdgesFn : passThroughDegree1EdgesFn, 1);

            adjacencyComponents = getAdjacencyComponents2(pinchGraph, cCIP->minimumDegree <= 1 ? doNotPassThroughDegree1EdgesFn : passThroughDegree1EdgesFn);

            ///////////////////////////////////////////////////////////////////////////
            // Constructing the flower.
            ///////////////////////////////////////////////////////////////////////////

            fillOutFlowerFromInputs(flower, cactusGraph, pinchGraph, adjacencyComponents);

#ifdef BEN_DEBUG
            flower_checkRecursive(flower);
            flower_checkNotEmpty(flower, 1);
#endif

            ///////////////////////////////////////////////////////////////////////////
            //Clean up remaining stuff.
            ///////////////////////////////////////////////////////////////////////////

            destructCactusGraph(cactusGraph);
            stList_destruct(adjacencyComponents);
            destructPinchGraph(pinchGraph);

            //assert(0);
            st_logInfo("Ran the core pipeline script\n");
            return 0;
        }
    }
    st_logCritical("We've run out of iterations of annealing before finishing!\n");
    assert(0);
}

int32_t cactusCorePipeline(Flower *flower, CactusCoreInputParameters *cCIP,
        struct PairwiseAlignment *(*getNextAlignment)(), void(*startAlignmentStack)()) {
    ////////////////////////////////////////////////
    //Check the flower to fill in terminal, and get rid of the group it contains and any terminal chain.
    ////////////////////////////////////////////////

#ifdef BEN_DEBUG
    assert(!flower_builtBlocks(flower)); //We can't do this if we've already built blocks for the flower!.
    flower_check(flower);
    assert(flower_isTerminal(flower));
    assert(flower_getGroupNumber(flower) == 1);
    assert(group_isLeaf(flower_getFirstGroup(flower))); //this should be true by the previous assert
    //Destruct any chain
    assert(flower_getChainNumber(flower) <= 1);
#endif
    if (flower_getChainNumber(flower) == 1) {
        Chain *chain = flower_getFirstChain(flower);
        chain_destruct(chain);
    }
    group_destruct(flower_getFirstGroup(flower));

    return cactusCorePipelineP(flower, cCIP, getNextAlignment, startAlignmentStack, cCIP->annealingRoundsLength, 1);
}
