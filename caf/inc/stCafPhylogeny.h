/*
 * stCafPhylogeny.h
 *
 *  Created on: Jun 2, 2014
 *      Author: benedictpaten
 */

#ifndef STCAFPHYLOGENY_H_
#define STCAFPHYLOGENY_H_

#include "sonLib.h"
#include "stPinchPhylogeny.h"

enum stCaf_TreeBuildingMethod {
    NEIGHBOR_JOINING,        // Traditional neighbor-joining algorithm.
    GUIDED_NEIGHBOR_JOINING  // Neighbor-joining guided by a species
                             // tree, which prefers to make trees that
                             // have a low amount of dups/losses if
                             // there is not a substantial amount of
                             // information to the contrary.
};

enum stCaf_RootingMethod {
    OUTGROUP_BRANCH, // Root all trees on their longest outgroup branch.
    LONGEST_BRANCH,  // Root all trees on their longest branch.
    BEST_RECON       // Root all trees so that they have minimal number of dups.
};

// Scoring method (how to choose the "best" tree from the canonical
// tree plus the bootstraps).
enum stCaf_ScoringMethod {
    RECON_COST,            // Minimize (dups + losses) according to
                           // the reconciliation

    NUCLEOTIDE_LIKELIHOOD, // Maximize likelihood of the bases in the columns

    RECON_LIKELIHOOD,      // Maximize likelihood according to a
                           // simple method where duplications are
                           // assumed to follow a Poisson distribution
                           // parameterized by the branch length

    COMBINED_LIKELIHOOD    // Maximize reconLikelihood * nucLikelihood
};

// Parameters for a phylogeny run on a particular flower.
typedef struct {
    // See above for definition and options.
    enum stCaf_TreeBuildingMethod treeBuildingMethod;
    // See above for definition and options.
    enum stCaf_RootingMethod rootingMethod;
    // See above for definition and options.
    enum stCaf_ScoringMethod scoringMethod;
    // The amount to weight the breakpoint information in the graph
    // by. Importantly, this is not relative to the substitution
    // information, but just scales the raw breakpoint
    // similarities/differences.
    double breakpointScalingFactor;
    // Skip building trees for single-copy blocks: blocks that have at
    // most 1 segment per species.
    bool skipSingleCopyBlocks;
    // Some trees will partition a block in a way that leaves a
    // segment by itself. This controls whether to keep those segments
    // in a single-degree block or to leave it without a block. It's
    // probably best to leave this off, as single-degree blocks will
    // block further alignment during bar phase.
    bool keepSingleDegreeBlocks;
    // Parameter for guided neighbor-joining. The penalty (in raw # of
    // differences) for choosing a join that implies a dup. This value
    // is multiplied by maxBaseDistance.
    // Ignored if not using guided neighbor-joining.
    double costPerDupPerBase;
    // Parameter for guided neighbor-joining. Similar to
    // costPerDupPerBase, except it is applied n times for a join that
    // implies n losses.
    // Ignored if not using guided neighbor-joining.
    double costPerLossPerBase;
    // The maximum distance (in bases) that will be traversed out on
    // each side looking for breakpoint/substitution information to
    // build trees with.
    int64_t maxBaseDistance;
    // The maximum distance (in blocks) that will be traversed out on
    // each side looking for breakpoint/substitution information to
    // build trees with.
    int64_t maxBlockDistance;
    // The total number of trees to build. A canonical tree with no
    // resampling is always built, plus (numTrees - 1) bootstraps.
    int64_t numTrees;
    // Don't count unaligned bases toward the
    // maxBaseDistance/maxBlockDistance limit, since they will never
    // be included in feature columns.
    bool ignoreUnalignedBases;
    // Only include features which include information for every segment
    // in the block.
    bool onlyIncludeCompleteFeatureBlocks;
    // To save a lot of time, do all splits with support higher than
    // this value at once, and then recompute the necessary trees,
    // instead of taking one at a time and recomputing after every
    // split. This assumes that very well supported splits will likely
    // be good no matter what the breakpoint information around them
    // is, which should usually be correct.
    // Any value greater than 1.0 disables this.
    bool doSplitsWithSupportHigherThanThisAllAtOnce;
} stCaf_PhylogenyParameters;

// Split a block according to a partition.
void splitBlock(stPinchBlock *block, stList *partitions,
                bool allowSingleDegreeBlocks);

/*
 * Build tree for each block and then use it to partition homologies in the block into
 * those which occur before and after the speciation.
 */
void stCaf_buildTreesToRemoveAncientHomologies(stPinchThreadSet *threadSet,
                                               stHash *threadStrings,
                                               stSet *outgroupThreads,
                                               Flower *flower,
                                               stCaf_PhylogenyParameters *params,
                                               FILE *debugFile,
                                               const char *referenceEventHeader);

/*
 * Gets the string for each pinch thread in a set.
 */
stHash *stCaf_getThreadStrings(Flower *flower, stPinchThreadSet *threadSet);

/*
 * Gets the sub-set of threads that are part of outgroup events.
 */
stSet *stCaf_getOutgroupThreads(Flower *flower, stPinchThreadSet *threadSet);


#endif /* STCAFPHYLOGENY_H_ */
