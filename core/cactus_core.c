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
#include "sonLib.h"

void writePinchGraph(char *name, struct PinchGraph *pinchGraph,
        struct List *biConnectedComponents, struct List *groups) {
    FILE *fileHandle;
    fileHandle = fopen(name, "w");
    struct hashtable *hash = createHashColouringPinchEdgesByChains(pinchGraph,
            biConnectedComponents);
    writeOutPinchGraphWithChains(pinchGraph, hash, groups, fileHandle);
    fclose(fileHandle);
    hashtable_destroy(hash, TRUE, FALSE);
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
    if (piece->start >= 1) {
        return sequence_getString(sequence, piece->start, piece->end
                - piece->start + 1, 1);
    } else {
        return sequence_getString(sequence, -piece->end, piece->end
                - piece->start + 1, 0);
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
    if (piece->end - piece->start + 1 > 2 * filterParameters->trim) { //only add to graph if non trivial in length.
        //Do the trim.
        piece->end -= filterParameters->trim;
        piece->start += filterParameters->trim;
        piece2->end -= filterParameters->trim;
        piece2->start += filterParameters->trim;
        assert(piece->end - piece->start == piece2->end - piece2->start);
        assert(piece->end - piece->start >= 0);

        //Now filter by repeat content.
        if (!filterParameters->alignRepeats) {
            char *string1 = piece_getString(piece, filterParameters->net);
            char *string2 = piece_getString(piece2, filterParameters->net);
            if (!containsRepeatBases(string1) && !containsRepeatBases(string2)) {
                pinchMergePiece(pinchGraph, piece, piece2,
                        vertexAdjacencyComponents);
            }
            free(string1);
            free(string2);
        } else {
            pinchMergePiece(pinchGraph, piece, piece2,
                    vertexAdjacencyComponents);
        }
    }
}

CactusCoreInputParameters *constructCactusCoreInputParameters() {
    CactusCoreInputParameters *cCIP = (CactusCoreInputParameters *) st_malloc(
            sizeof(CactusCoreInputParameters));
    //Everything is essentially 'turned off' by default.
    cCIP->writeDebugFiles = 0;
    cCIP->annealingRounds = 1;
    cCIP->alignRepeatsAtRound = 0;
    cCIP->trim = 0;
    cCIP->trimChange = 0.0;
    cCIP->minimumTreeCoverage = 0.0;
    cCIP->minimumBlockLength = 0;
    cCIP->minimumBlockLengthChange = 0.0;
    cCIP->minimumChainLength = 0;
    cCIP->minimumChainLengthChange = 0.0;
    cCIP->deannealingRounds = 1.0;
    return cCIP;
}

void destructCactusCoreInputParameters(CactusCoreInputParameters *cCIP) {
    free(cCIP);
}

struct CactusGraph *cactusCorePipeline_2(struct PinchGraph *pinchGraph,
        Net *net) {
    struct CactusGraph *cactusGraph;
    struct List *threeEdgeConnectedComponents;
    int32_t startTime;

    ///////////////////////////////////////////////////////////////////////////
    // Linking stub components to the sink component (if they haven't been already been).
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);
    linkStubComponentsToTheSinkComponent(pinchGraph, net);
    checkPinchGraph(pinchGraph);
    st_logInfo("Linked stub components to the sink component in: %i seconds\n",
            time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Constructing the basic cactus.
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);
    computeCactusGraph(pinchGraph, &cactusGraph, &threeEdgeConnectedComponents);
    st_logInfo("Constructed the initial cactus graph in: %i seconds\n", time(
            NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Circularising the stems in the cactus.
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);
    circulariseStems(cactusGraph);
    st_logInfo("Constructed the 2-edge component only cactus graph\n");
    checkCactusContainsOnly2EdgeConnectedComponents(cactusGraph);
    st_logInfo(
            "Checked the cactus contains only 2-edge connected components in: %i seconds\n",
            time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup.
    ///////////////////////////////////////////////////////////////////////////

    destructList(threeEdgeConnectedComponents);
    return cactusGraph;
}

struct List *getChosenBlockPinchEdges(stSortedSet *chosenBlocks,
        struct PinchGraph *pinchGraph) {
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

int32_t cactusCorePipeline(Net *net, CactusCoreInputParameters *cCIP,
        struct PairwiseAlignment *(*getNextAlignment)(),
        void(*startAlignmentStack)(), int32_t terminateRecursion) {
    struct PinchGraph *pinchGraph;
    struct PinchVertex *vertex;
    struct CactusGraph *cactusGraph;
    int32_t i, k, startTime;
    struct List *biConnectedComponents;
    struct PairwiseAlignment *pairwiseAlignment;
    struct List *list;
    struct hashtable *vertexAdjacencyComponents;

    assert(!net_builtBlocks(net)); //We can't do this if we've already built blocks for the net!.

    ///////////////////////////////////////////////////////////////////////////
    //Setup the basic pinch graph
    ///////////////////////////////////////////////////////////////////////////

    pinchGraph = constructPinchGraph(net);

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

    float trim = cCIP->trim;
    float minimumChainLength = cCIP->minimumChainLength;
    float minimumBlockLength = cCIP->minimumBlockLength;

    vertexAdjacencyComponents = create_hashtable(pinchGraph->vertices->length*2, hashtable_key, hashtable_equalKey, NULL, free);

    //Build a hash putting the vertices all in the same adjacency component.
    for(i=0; i<pinchGraph->vertices->length; i++) {
        hashtable_insert(vertexAdjacencyComponents, pinchGraph->vertices->list[i], constructInt(0));
    }

    int32_t loop = 0;
    while(1) {

#ifdef BEN_DEBUG
        ///////////////////////////////////////////////////////////////////////////
        //  Check the adjacency vertex components.
        ///////////////////////////////////////////////////////////////////////////

        assert((int32_t)hashtable_count(vertexAdjacencyComponents) == pinchGraph->vertices->length);
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
        pairwiseAlignment = getNextAlignment(); //we assume we own this memory, and will clean it up.
        st_logInfo("Now doing the pinch merges:\n");
        i = 0;

        struct FilterAlignmentParameters *filterParameters = (struct FilterAlignmentParameters *)st_malloc(sizeof(struct FilterAlignmentParameters));
        assert(trim >= 0);
        filterParameters->trim = trim;
        filterParameters->alignRepeats = loop >= cCIP->alignRepeatsAtRound; //cCIP->alignRepeats;
        filterParameters->net = net;

        while(pairwiseAlignment != NULL) {
            st_logDebug("Alignment : %i , score %f\n", i++, pairwiseAlignment->score);
            logPairwiseAlignment(pairwiseAlignment);
            pinchMerge(pinchGraph, pairwiseAlignment,
                    (void (*)(struct PinchGraph *pinchGraph, struct Piece *, struct Piece *, struct hashtable *, void *))filterPieceAndThenAddToGraph,
                    filterParameters, vertexAdjacencyComponents);
            destructPairwiseAlignment(pairwiseAlignment); //cleanup the previous alignment
            pairwiseAlignment = getNextAlignment();
        }
        free(filterParameters);
        st_logInfo("Finished pinch merges\n");

#ifdef BEN_DEBUG
        for(i=0; i<pinchGraph->vertices->length; i++) {
            assert(hashtable_search(vertexAdjacencyComponents, pinchGraph->vertices->list[i]) != NULL);
        }
        assert((int32_t)hashtable_count(vertexAdjacencyComponents) == pinchGraph->vertices->length);
#endif

        //Cleanup the adjacency component vertex hash.
        hashtable_destroy(vertexAdjacencyComponents, 1, 0);

        checkPinchGraph(pinchGraph); //check the graph is all good.
        st_logInfo("Pinched the graph in: %i seconds\n", time(NULL) - startTime);

        removeTrivialGreyEdgeComponents(pinchGraph, pinchGraph->vertices, net); //remove any pointless adjacencies.
        st_logInfo("After removing the trivial graph components the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));
        checkPinchGraph(pinchGraph);

        ////////////////////////////////////////////////
        // Compute the cactus graph
        ////////////////////////////////////////////////

        cactusGraph = cactusCorePipeline_2(pinchGraph, net);

        ////////////////////////////////////////////////
        // Get sorted bi-connected components.
        ////////////////////////////////////////////////

        biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);

        ////////////////////////////////////////////////
        // Loop a bunch of times to progressively remove longer and longer (upto minimum chain length) chains.
        ////////////////////////////////////////////////

        //assert(cCIP->deannealingRounds >= 1.0);
        if(cCIP->deannealingRounds >= 1) {
            float deannealingChainLengthStepSize = ((float)minimumChainLength) / cCIP->deannealingRounds;
            /*
            if(cCIP->deannealingRounds >= 1.0) {
                deannealingChainLengthStepSize = ((float)minimumChainLength) / cCIP->deannealingRounds;
            }
            else {
                deannealingChainLengthStepSize = minimumChainLength;
            }*/
            float deannealingChainLength = deannealingChainLengthStepSize;
    //if(loop+1 < cCIP->annealingRounds) {
            while(1) {
                ///////////////////////////////////////////////////////////////////////////
                // Choosing a block subset to undo.
                ///////////////////////////////////////////////////////////////////////////

                startTime = time(NULL);
                //Get all the blocks.
                stSortedSet *allBlocksOfDegree2OrHigher = filterBlocksByTreeCoverageAndLength(biConnectedComponents, net, 0.0, 2, 0, 0, pinchGraph);
                //Get the blocks we want to keep
                stSortedSet *chosenBlocksToKeep = filterBlocksByTreeCoverageAndLength(biConnectedComponents,
                        net, cCIP->minimumTreeCoverage, 0, minimumBlockLength, deannealingChainLength, pinchGraph);
                //Now get the blocks to undo by computing the difference.
                stSortedSet *blocksToUndo = stSortedSet_getDifference(allBlocksOfDegree2OrHigher, chosenBlocksToKeep);
                stSortedSet_destruct(chosenBlocksToKeep);
                stSortedSet_destruct(allBlocksOfDegree2OrHigher);

                if(stSortedSet_size(blocksToUndo) > 0) {
                    //now report the results
                    logTheChosenBlockSubset(biConnectedComponents, blocksToUndo, pinchGraph, net);
                    st_logInfo("I have chosen %i blocks which meet the requirements to be undone\n", stSortedSet_size(blocksToUndo));

                    ///////////////////////////////////////////////////////////////////////////
                    // Undo the blocks.
                    ///////////////////////////////////////////////////////////////////////////

                    list = getChosenBlockPinchEdges(blocksToUndo, pinchGraph);
                    removeOverAlignedEdges(pinchGraph, 0.0, INT32_MAX, list, 0, net);
                    destructList(list);
                    st_logInfo("After removing edges which were not chosen, the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));
                    removeTrivialGreyEdgeComponents(pinchGraph, pinchGraph->vertices, net);
                    st_logInfo("After removing the trivial graph components the graph has %i vertices and %i black edges\n", pinchGraph->vertices->length, avl_count(pinchGraph->edges));

                    ///////////////////////////////////////////////////////////////////////////
                    // Cleanup the old cactus graph (it is now out of sync with the pinch graph,
                    // after undoing the selected edges).
                    ///////////////////////////////////////////////////////////////////////////

                    destructList(biConnectedComponents);
                    destructCactusGraph(cactusGraph);

                    ////////////////////////////////////////////////
                    // Re-compute the cactus graph
                    ////////////////////////////////////////////////

                    cactusGraph = cactusCorePipeline_2(pinchGraph, net);

                    ////////////////////////////////////////////////
                    // Get the sorted bi-connected components, again
                    ////////////////////////////////////////////////

                    biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);
                }
                stSortedSet_destruct(blocksToUndo);

                if(deannealingChainLength >= minimumChainLength) {
                    break;
                }
                deannealingChainLength += deannealingChainLengthStepSize;
            }
        }
//}

        ///////////////////////////////////////////////////////////////////////////
        // Choosing a block subset to keep in the final set of chains.
        ///////////////////////////////////////////////////////////////////////////

        if(++loop < cCIP->annealingRounds) {
            ///////////////////////////////////////////////////////////////////////////
            // Calculate the blocks used to partition the graph for the next loop.
            ///////////////////////////////////////////////////////////////////////////

            stSortedSet *chosenBlocks = filterBlocksByTreeCoverageAndLength(biConnectedComponents,
                    net, 0.0, 2, 0, 0, pinchGraph); //this gets all blocks that are aligned to something else..

#ifdef BEN_DEBUG
            /*
             * This checks that the chosen blocks all have more than one sequence in them.
             */
            stSortedSetIterator *chosenBlocksIterator = stSortedSet_getIterator(chosenBlocks);
            struct CactusEdge *cactusEdge;
            while((cactusEdge = stSortedSet_getNext(chosenBlocksIterator)) != NULL) {
                assert(cactusEdge->pieces->length > 1);
            }
            stSortedSet_destructIterator(chosenBlocksIterator);

            /*
             * This test checks all blocks not in the chosen blocks have degree 1 or are stubs.
             */
            stSortedSet *allBlocks = filterBlocksByTreeCoverageAndLength(biConnectedComponents, net, 0.0, 0, 0, 0, pinchGraph);
            assert(stSortedSet_size(allBlocks) == cactusGraph_getEdgeNumber(cactusGraph) - net_getStubEndNumber(net)); //check that this does slurp up all the block edges in the graph except those representing stub ends.
            //Now get the blocks to undo by computing the difference.
            stSortedSet *otherBlocks = stSortedSet_getDifference(allBlocks, chosenBlocks);
            assert(stSortedSet_size(allBlocks) == stSortedSet_size(chosenBlocks) + stSortedSet_size(otherBlocks)); //just to be sure.
            stSortedSetIterator *otherBlocksIterator = stSortedSet_getIterator(otherBlocks);
            while((cactusEdge = stSortedSet_getNext(otherBlocksIterator)) != NULL) {
                assert(cactusEdge->pieces->length == 1);
            }
            stSortedSet_destructIterator(otherBlocksIterator);
            stSortedSet_destruct(allBlocks);
            stSortedSet_destruct(otherBlocks);
#endif


            //now report the results
            logTheChosenBlockSubset(biConnectedComponents, chosenBlocks, pinchGraph, net);
            st_logInfo("I have chosen %i blocks which meet the requirements to keep\n", stSortedSet_size(chosenBlocks));

            ///////////////////////////////////////////////////////////////////////////
            // Calculate the adjacency components for the next loop.
            ///////////////////////////////////////////////////////////////////////////

            //Build a hash putting each vertex in its own adjacency component.
            //Iterate through all the edges, keeping only those whose degree is greater than zero.
            vertexAdjacencyComponents = create_hashtable(pinchGraph->vertices->length*2, hashtable_key, hashtable_equalKey, NULL, free);
            list = getChosenBlockPinchEdges(chosenBlocks, pinchGraph);
            struct List *groupsList = getRecursiveComponents2(pinchGraph, list);
            for(i=0; i<groupsList->length; i++) {
                struct List *vertices = groupsList->list[i];
                for(k=0; k<vertices->length; k++) {
                    hashtable_insert(vertexAdjacencyComponents, vertices->list[k], constructInt(i));
                }
            }
            destructList(list);
            destructList(groupsList);

            ///////////////////////////////////////////////////////////////////////////
            // Modify parameters for next loop
            ///////////////////////////////////////////////////////////////////////////

            minimumBlockLength += cCIP->minimumBlockLengthChange;
            minimumBlockLength = minimumBlockLength < 0.0 ? 0.0 : minimumBlockLength;
            minimumChainLength += cCIP->minimumChainLengthChange;
            minimumChainLength = minimumChainLength < 0.0 ? 0.0 : minimumChainLength;
            trim += cCIP->trimChange;
            trim = trim < 0.0 ? 0.0 : trim;


            ///////////////////////////////////////////////////////////////////////////
            // Cleanup the loop.
            ///////////////////////////////////////////////////////////////////////////

            destructCactusGraph(cactusGraph);
            stSortedSet_destruct(chosenBlocks);
            destructList(biConnectedComponents);
        }
        else {
            ///////////////////////////////////////////////////////////////////////////
            // Constructing the net.
            ///////////////////////////////////////////////////////////////////////////

            stSortedSet *chosenBlocks = filterBlocksByTreeCoverageAndLength(biConnectedComponents,
                                net, 0.0, terminateRecursion ? 0 : 2, 0, 0, pinchGraph);
            //assert(stSortedSet_size(chosenBlocks) == cactusGraph_getEdgeNumber(cactusGraph) - net_getStubEndNumber(net)); //check that this does slurp up all the block edges in the graph except those representing stub ends.
            fillOutNetFromInputs(net, cactusGraph, pinchGraph, chosenBlocks);

            if(cCIP->writeDebugFiles) {
                ///////////////////////////////////////////////////////////////////////////
                //Write out the graphs.
                ///////////////////////////////////////////////////////////////////////////

                st_logDebug("Writing out dot formatted final pinch graph showing all chains\n");
                writePinchGraph("pinchGraph2.dot", pinchGraph, biConnectedComponents, NULL);
                st_logDebug("Finished writing out final pinch graph showing all chains\n");

                st_logDebug("Writing out dot formatted final pinch graph showing chosen blocks\n");
                list = constructEmptyList(0, NULL);
                listAppend(list, chosenBlocks);
                writePinchGraph("pinchGraph3.dot", pinchGraph, list, NULL);
                destructList(list);
                st_logDebug("Finished writing out final pinch graph\n");

                st_logDebug("Writing out dot formatted version of final cactus graph\n");
                writeCactusGraph("cactusGraph.dot", pinchGraph, cactusGraph);
                st_logDebug("Finished writing out dot formatted version of cactus graph\n");
            }

            stSortedSet_destruct(chosenBlocks);
            break;
        }

    }

    ///////////////////////////////////////////////////////////////////////////
    //Clean up remaining stuff.
    ///////////////////////////////////////////////////////////////////////////

    destructCactusGraph(cactusGraph);
    destructList(biConnectedComponents);
    destructPinchGraph(pinchGraph);

    st_logInfo("Ran the core pipeline script\n");
    return 0;
}
