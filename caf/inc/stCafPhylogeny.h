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

enum stCaf_RootingMethod {
    OUTGROUP_BRANCH,
    LONGEST_BRANCH,
    BEST_RECON
};

enum stCaf_ScoringMethod {
    RECON_COST,
    NUCLEOTIDE_LIKELIHOOD
};

/*
 * Build tree for each block and then use it to partition homologies in the block into
 * those which occur before and after the speciation.
 */
void stCaf_buildTreesToRemoveAncientHomologies(stPinchThreadSet *threadSet, stHash *threadStrings, stSet *outgroupThreads, Flower *flower, int64_t maxBaseDistance, int64_t maxBlockDistance, int64_t numTrees, enum stCaf_RootingMethod rootingMethod, enum stCaf_ScoringMethod scoringMethod, double breakPointScalingFactor, bool skipSingleCopyBlocks);

/*
 * Gets the string for each pinch thread in a set.
 */
stHash *stCaf_getThreadStrings(Flower *flower, stPinchThreadSet *threadSet);

/*
 * Gets the sub-set of threads that are part of outgroup events.
 */
stSet *stCaf_getOutgroupThreads(Flower *flower, stPinchThreadSet *threadSet);


#endif /* STCAFPHYLOGENY_H_ */
