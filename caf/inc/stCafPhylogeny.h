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

// A "split branch": a branch in a block tree that, if removed, would
// produce a partition of the leaf set that would remove an ancient
// homology. In practice, this means that split branches have a
// duplication node as a parent, which should be reconciled at or
// above the reference event in order for the splits to work out
// properly.
typedef struct {
    stTree *child; // Child of the branch.
    stPinchBlock *block; // Block the tree refers to (can and should
                         // refer to more than the child subtree).
    double support; // Bootstrap support for this branch.
} stCaf_SplitBranch;

// A choice of tree-building methods
enum stCaf_TreeBuildingMethod {
    NEIGHBOR_JOINING,        // Traditional neighbor-joining algorithm.
    GUIDED_NEIGHBOR_JOINING  // Neighbor-joining guided by a species
                             // tree, which prefers to make trees that
                             // have a low amount of dups/losses if
                             // there is not a substantial amount of
                             // information to the contrary.
};

// A choice of rooting methods for each tree built
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
    // The amount to weight the nucleotide information in the graph
    // by.
    double nucleotideScalingFactor;
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
    // Number of additional threads to spawn to do tree-building with
    // (has to be more than 0). The master thread is almost always
    // stalled while tree-building is running, so you should expect at
    // most numTreeBuildingThreads cpus to be occupied.
    int64_t numTreeBuildingThreads;
} stCaf_PhylogenyParameters;

// Split a block according to a partition (a list of lists of
// stIntTuples representing the segment indices in the block).
void stCaf_splitBlock(stPinchBlock *block, stList *partitions,
                      bool allowSingleDegreeBlocks);

// Compare two bootstrap scores of split branches. Use the pointer
// value of the branches as a tiebreaker since we are using a sorted
// set and don't want to merge together all branches with the same
// support value.
int stCaf_SplitBranch_cmp(stCaf_SplitBranch *branch1,
                          stCaf_SplitBranch *branch2);

// Find new split branches from the block and add them to the sorted set.
// speciesToSplitOn is just the species that are on the path from the
// reference node to the root.
void stCaf_findSplitBranches(stPinchBlock *block, stTree *tree,
                             stSortedSet *splitBranches,
                             stSet *speciesToSplitOn);

// Remove any split branches that appear in this tree from the
// set of split branches.
void stCaf_removeSplitBranches(stPinchBlock *block, stTree *tree,
                               stSet *speciesToSplitOn,
                               stSortedSet *splitBranches);

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
