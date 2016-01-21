/*
 * stCore.h
 *
 *  Created on: 28 Apr 2012
 *      Author: benedictpaten
 */

#ifndef STCAF_H_
#define STCAF_H_

#include "sonLib.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stCactusGraphs.h"
#include "stOnlineCactus.h"
#include "cactus.h"

///////////////////////////////////////////////////////////////////////////
// Setup the pinch graph from a cactus graph
///////////////////////////////////////////////////////////////////////////

/*
 * Add the proper caps (grouped into blocks representing ends) to a
 * pinch graph from a flower containing no blocks. Stubs are
 * represented by length 1 blocks.
 */
void stCaf_addCapsToPinchGraph(Flower *flower, stPinchThreadSet *threadSet);

/*
 * Add the proper threads to an empty pinch graph. Each thread is
 * given the name of its 5' cap.
 */
void stCaf_addThreadsToPinchGraph(Flower *flower, stPinchThreadSet *threadSet);

/*
 * Initialises the flower for filling out and creates a pinch graph.
 */
stPinchThreadSet *stCaf_setup(Flower *flower);

/*
 * Initializes the flower for filling out, and creates an online pinch to cactus mapping.
 * A pointer to the online cactus is stored in the "cactus" parameter.
 */
stPinchThreadSet *stCaf_setupForOnlineCactus(Flower *flower, stOnlineCactus **cactus);

void stCaf_disableAndCleanupOnlineCactus(stPinchThreadSet *threadSet, stOnlineCactus *cactus);

///////////////////////////////////////////////////////////////////////////
// Annealing fuctions -- adding alignments to pinch graph
///////////////////////////////////////////////////////////////////////////

/*
 * Add the set of alignments, represented as pinches, to the graph.
 */
void stCaf_anneal(stPinchThreadSet *threadSet, stPinchIterator *pinchIterator, bool (*filterFn)(stPinchSegment *, stPinchSegment *));

/*
 * Add the set of alignments, represented as pinches, to the graph, allowing alignments only between segments in the same component.
 */
void stCaf_annealBetweenAdjacencyComponents(stPinchThreadSet *threadSet, stPinchIterator *pinchIterator, bool (*filterFn)(stPinchSegment *, stPinchSegment *));

/*
 * Joins all trivial boundaries, but not joining stub boundaries.
 */
void stCaf_joinTrivialBoundaries(stPinchThreadSet *threadSet);

typedef enum {
    PRESERVE_NON_UNDOABLE_CHAINS,
    REMOVE_NON_UNDOABLE_CHAINS,
    ONLY_UNDO,
    ONLY_REMOVE,
    NONE
} stCaf_meltingMethod;

/*
 * Anneals from the pinch iterator, but checks after every pinch to
 * prevent chains with length less than minimumChainLength from
 * forming.
 */
void stCaf_annealPreventingSmallChains(Flower *flower, stPinchThreadSet *threadSet,
                                       stOnlineCactus *cactus,
                                       const char *alignmentsFile,
                                       stList *alignments,
                                       int64_t alignmentTrim,
                                       bool (*filterFn)(stPinchSegment *, stPinchSegment *),
                                       stList *minimumChainLengths,
                                       stCaf_meltingMethod meltingMethod,
                                       int64_t numAlignmentsPerBatch,
                                       int64_t maxNumAlignments,
                                       double maxRedundantFraction,
                                       const char *dumpPath);

/*
 * Dump statistics on the maximum block degree and total aligned bases to the dump file.
 */
void dumpMaxBlockDegreeAndTAB(stPinchThreadSet *threadSet, FILE *dumpFile);


///////////////////////////////////////////////////////////////////////////
// Melting fuctions -- removing alignments from the pinch graph
///////////////////////////////////////////////////////////////////////////

void stCaf_undoChainsSmallerThanThis_onlyUndo(stOnlineCactus *cactus, stPinchThreadSet *threadSet,
                                              stList *pinches, stPinchUndo *undo,
                                              int64_t minimumChainLength);

void stCaf_undoChainsSmallerThanThis_onlyRemove(stOnlineCactus *cactus, stPinchThreadSet *threadSet,
                                                stList *pinches, stPinchUndo *undo,
                                                int64_t minimumChainLength);

void stCaf_undoChainsSmallerThanThis_preserveNonUndoableChains(stOnlineCactus *cactus, stPinchThreadSet *threadSet,
                                                               stList *pinches, stPinchUndo *undo,
                                                               int64_t minimumChainLength);

void stCaf_undoChainsSmallerThanThis_removeNonUndoableChains(stOnlineCactus *cactus, stPinchThreadSet *threadSet,
                                                             stList *pinches, stPinchUndo *undo,
                                                             int64_t minimumChainLength);

/*
 * Get a list of blocks that participate in a chain of length less
 * than minimumChainLength.
 */
stList *stCaf_getBlocksInChainsLessThanGivenLength(stCactusGraph *cactusGraph, int64_t minimumChainLength);

/*
 * Removes homologies from the graph.
 */
void stCaf_melt(Flower *flower, stPinchThreadSet *threadSet, bool blockFilterfn(stPinchBlock *), int64_t blockEndTrim,
        int64_t minimumChainLength, bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds);

/*
 * Function used to determine if blocks contains sufficient numbers of sequences of ingroup/outgroup species.
 */
bool stCaf_containsRequiredSpecies(stPinchBlock *pinchBlock, Flower *flower, int64_t minimumIngroupDegree,
        int64_t minimumOutgroupDegree, int64_t minimumAllDegree);

/*
 * Returns the proportion of the tree covered by the block.
 */
bool stCaf_treeCoverage(stPinchBlock *pinchBlock, Flower *flower);

/*
 * Returns 1 if any ingroup species is present in multiple copies.
 */
bool stCaf_containsMultipleCopiesOfIngroupSpecies(stPinchBlock *pinchBlock, Flower *flower);

/*
 * Returns 1 if any outgroup species is present in multiple copies.
 */
bool stCaf_containsMultipleCopiesOfOutgroupSpecies(stPinchBlock *pinchBlock, Flower *flower);

/*
 * Returns 1 if any species is present in multiple copies.
 */
bool stCaf_containsMultipleCopiesOfAnySpecies(stPinchBlock *pinchBlock, Flower *flower);

Event *getEvent(stPinchSegment *segment, Flower *flower);

/*
 * Simply returns the average degree of the blocks in the list.
 */
double stCaf_averageBlockDegree(stList *blocks);

/*
 * Returns the number of aligned bases represented in the blocks in the list.
 */
uint64_t stCaf_totalAlignedBases(stList *blocks);

///////////////////////////////////////////////////////////////////////////
// Pinch graph to cactus graph
///////////////////////////////////////////////////////////////////////////

/*
 * Function that merges together adjacency components.
 */
void *stCaf_mergeNodeObjects(void *a, void *b);

/*
 * Functions which converts a pinch graph into a cactus graph.
 */
stCactusGraph *stCaf_getCactusGraphForThreadSet(Flower *flower, stPinchThreadSet *threadSet, stCactusNode **startCactusNode,
        stList **deadEndComponent, bool attachEndsInFlower, int64_t minLengthForChromosome, double proportionOfUnalignedBasesForNewChromosome,
        bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds);

///////////////////////////////////////////////////////////////////////////
// Finishing: Converting a pinch graph into the flower hierarchy
///////////////////////////////////////////////////////////////////////////

/*
 * Create blocks for unaligned segments.
 */
void stCaf_makeDegreeOneBlocks(stPinchThreadSet *threadSet);

/*
 * Add the adjacencies between caps of a flower.
 * All ends must be in the flower before this function is called.
 */
void stCaf_addAdjacencies(Flower *flower);

/*
 * Takes a pinch graph for a flower and adds the alignments it contains back to the flower as a cactus.
 */
void stCaf_finish(Flower *flower, stPinchThreadSet *threadSet, int64_t chainLengthForBigFlower,
        int64_t longChain, int64_t minLengthForChromosome,
        double proportionOfUnalignedBasesForNewChromosome);

#endif /* STCAF_H_ */
