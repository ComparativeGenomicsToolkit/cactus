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
#include "cactus.h"

///////////////////////////////////////////////////////////////////////////
// Setup the pinch graph from a cactus graph
///////////////////////////////////////////////////////////////////////////

/*
 * Create a pinch graph from a flower containing no blocks. Each thread is given the name of its 5' cap,
 * and stubs are represented by length 1 blocks.
 */
stPinchThreadSet *stCaf_constructEmptyPinchGraph(Flower *flower);

/*
 * Initialises the flower for filling out and creates a pinch graph.
 */
stPinchThreadSet *stCaf_setup(Flower *flower);

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

///////////////////////////////////////////////////////////////////////////
// Melting fuctions -- removing alignments from the pinch graph
///////////////////////////////////////////////////////////////////////////

/*
 * Removes homologies from the graph.
 */
void stCaf_melt(Flower *flower, stPinchThreadSet *threadSet, bool blockFilterfn(stPinchBlock *), int64_t blockEndTrim,
        int64_t minimumChainLength, bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds);

/*
 * Removes any recoverable chains (those expected to be picked up by
 * bar phase) from the graph. Only chains that are recoverable *and*
 * where the "recoverabilityFilter" function returns 1 are removed.
 */
void stCaf_meltRecoverableChains(Flower *flower, stPinchThreadSet *threadSet, bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds, bool (*recoverabilityFilter)(stCactusEdgeEnd *, Flower *), int64_t maxNumIterations, int64_t maxRecoverableChainLength);

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

///////////////////////////////////////////////////////////////////////////
// Filtering fuctions -- filtering incoming alignments or entire blocks
///////////////////////////////////////////////////////////////////////////

/*
 * Must be used before any of stCaf_filterByOutgroup,
 * stCaf_relaxedFilterByOutgroup, or stCaf_filterByRepeatSpecies are
 * used.
 */
void stCaf_setFlowerForAlignmentFiltering(Flower *input);

/*
 * Filters incoming alignments by presence of outgroup, to ensure at
 * most one outgroup segment is in any block.
 */

bool stCaf_filterByOutgroup(stPinchSegment *segment1,
                            stPinchSegment *segment2);

/*
 * A "relaxed" version of the above which allows the addition of more
 * than one outgroup *segment* to a block, but never allows two
 * different blocks containing outgroups to be pinched.
 */
bool stCaf_relaxedFilterByOutgroup(stPinchSegment *segment1,
                                   stPinchSegment *segment2);

/*
 * Filters incoming alignments by presence of repeat species in
 * block. This code is inefficient and does not scale.
 */
bool stCaf_filterByRepeatSpecies(stPinchSegment *segment1,
                                 stPinchSegment *segment2);

/*
 * As above, but allows pinching a duplicated segment to a block, but
 * not pinching two duplicated blocks together.
 */
bool stCaf_relaxedFilterByRepeatSpecies(stPinchSegment *segment1,
                                        stPinchSegment *segment2);

bool stCaf_singleCopyIngroup(stPinchSegment *segment1,
                             stPinchSegment *segment2);

bool stCaf_relaxedSingleCopyIngroup(stPinchSegment *segment1,
                                    stPinchSegment *segment2);

/*
 * Filters block alignments that would merge blocks that each already contain
 * sequences from multiple species. The rationale of this filter is to avoid
 * aligning together paralogous alignments that predate the speciation event.
 */
bool stCaf_filterByMultipleSpecies(stPinchSegment *segment1, stPinchSegment *segment2);

/*
 * Forbids pinching together two copies within the same sequence.
 */
bool stCaf_singleCopyChr(stPinchSegment *segment1, stPinchSegment *segment2);

/*
 * Run stCaf_filterToEnsureCycleFreeIsolatedComponents so that every
 * non-alt thread (determined by its header not ending in "_alt") in
 * the given event is in a separate thread component, and has no cycles.
 */
void stCaf_setupHGVMFiltering(Flower *flower, stPinchThreadSet *threadSet,
                              char *eventName);

/*
 * Set up some additional information for
 * stCaf_filterToEnsureCycleFreeIsolatedComponents. Required to run
 * this first!
 */
void stCaf_setThreadsToBeCycleFreeIsolatedComponents(stPinchThreadSet *threadSet,
                                                     stSet *threads);

/*
 * Special filtering for a draft HGVM, where every non-alt chromosome
 * should be in its own component and should have no within-component
 * cycles.
 */
bool stCaf_filterToEnsureCycleFreeIsolatedComponents(stPinchSegment *segment1,
                                                     stPinchSegment *segment2);

/*
 * Returns true for chains that have an unequal number of ingroup
 * copies, e.g. 0 in ingroup 1, 1 in ingroup 2; or 2 in ingroup 1, 2
 * in ingroup 2, 3 in ingroup 3.
 */
bool stCaf_chainHasUnequalNumberOfIngroupCopies(stCactusEdgeEnd *chainEnd,
                                                Flower *flower);

/*
 * Returns true for chains that have an unequal number of copies
 * distributed among each of the ingroups, or has no copies among the
 * outgroups collectively.
 */
bool stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup(stCactusEdgeEnd *chainEnd,
                                                            Flower *flower);

/*
 * Function used to determine if blocks contains sufficient numbers of sequences of ingroup/outgroup species.
 */
bool stCaf_containsRequiredSpecies(stPinchBlock *pinchBlock, Flower *flower, int64_t minimumIngroupDegree,
        int64_t minimumOutgroupDegree, int64_t minimumDegree, int64_t minimumNumberOfSpecies);

/*
 * Returns the proportion of the tree covered by the block.
 */
bool stCaf_treeCoverage(stPinchBlock *pinchBlock, Flower *flower);

/*
 * Short way to get the event corresponding to a given segment.
 */
Event *stCaf_getEvent(stPinchSegment *segment, Flower *flower);

#endif /* STCAF_H_ */
