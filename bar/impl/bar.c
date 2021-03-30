#include "cactus.h"
#include "sonLib.h"
#include "endAligner.h"
#include "poaBarAligner.h"
#include "flowerAligner.h"
#include "rescue.h"
#include "commonC.h"
#include "stCaf.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stateMachine.h"
#include "pairwiseAligner.h"
#include "../../caf/inc/stCaf.h"

// OpenMP
#if defined(_OPENMP)
#include <omp.h>
#endif

PairwiseAlignmentParameters *pairwiseAlignmentParameters_constructFromCactusParams(CactusParams *params) {
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->gapGamma = cactusParams_get_float(params, 2, "bar", "gapGamma");
    p->splitMatrixBiggerThanThis = cactusParams_get_int(params, 2, "bar", "splitMatrixBiggerThanThis");
    p->anchorMatrixBiggerThanThis = cactusParams_get_int(params, 2, "bar", "anchorMatrixBiggerThanThis");
    p->repeatMaskMatrixBiggerThanThis = cactusParams_get_int(params, 2, "bar", "repeatMaskMatrixBiggerThanThis");
    p->diagonalExpansion = cactusParams_get_int(params, 2, "bar", "diagonalExpansion");
    p->constraintDiagonalTrim = cactusParams_get_int(params, 2, "bar", "constraintDiagonalTrim");
    p->alignAmbiguityCharacters = cactusParams_get_int(params, 2, "bar", "alignAmbiguityCharacters");
    return p;
}

stPinch *getNextAlignedPairAlignment(stSortedSetIterator *it, stPinch *pinchToFillOut) {
    AlignedPair *alignedPair = stSortedSet_getNext(it);
    if (alignedPair == NULL) {
        return NULL;
    }
    stPinch_fillOut(pinchToFillOut, alignedPair->subsequenceIdentifier, alignedPair->reverse->subsequenceIdentifier, alignedPair->position,
                    alignedPair->reverse->position, 1, alignedPair->strand == alignedPair->reverse->strand);
    return pinchToFillOut;
}

bool blockFilterFn(stPinchBlock *pinchBlock, void *extraArg) {
    FilterArgs *f = extraArg;
    return !stCaf_containsRequiredSpecies(pinchBlock, f->flower, f->minimumIngroupDegree, f->minimumOutgroupDegree, f->minimumDegree, f->minimumNumberOfSpecies);
}

void bar(stList *flowers, CactusParams *params, CactusDisk *cactusDisk, stList *listOfEndAlignmentFiles, bool cleanupMemory) {
    //////////////////////////////////////////////
    //Parse the many, many necessary parameters from the params file
    //////////////////////////////////////////////

    // Hardcoded parameters
    int64_t chainLengthForBigFlower = 1000000;
    int64_t longChain = 2;
    int64_t maximumLength = 1500;

    int64_t spanningTrees = cactusParams_get_int(params, 2, "bar", "spanningTrees");
    bool useProgressiveMerging = cactusParams_get_int(params, 2, "bar", "useProgressiveMerging");
    float matchGamma = cactusParams_get_float(params, 2, "bar", "matchGamma");

    //bool useBanding = cactusParams_get_int(params, 2, "bar", "useBanding");
    //int64_t minimumSizeToRescue = cactusParams_get_int(params, 2, "bar", "minimumSizeToRescue");
    //double minimumCoverageToRescue = cactusParams_get_float(params, 2, "bar", "minimumCoverageToRescue");

    // toggle from pecan to abpoa for multiple alignment, by setting to non-zero
    // Note that poa uses about N^2 memory, so maximum value is generally in 10s of kb
    int64_t usePoa = cactusParams_get_int(params, 2, "bar", "partialOrderAlignment");
    int64_t poaWindow = cactusParams_get_int(params, 2, "bar", "partialOrderAlignmentWindow");
    int64_t maskFilter = cactusParams_get_int(params, 2, "bar", "partialOrderAlignmentMaskFilter");
    int64_t poaBandConstant = cactusParams_get_int(params, 2, "bar", "partialOrderAlignmentBandConstant");
    //defaults from abpoa
    double poaBandFraction = cactusParams_get_float(params, 2, "bar", "partialOrderAlignmentBandFraction");

    PairwiseAlignmentParameters *pairwiseAlignmentParameters = pairwiseAlignmentParameters_constructFromCactusParams(params);
    bool pruneOutStubAlignments = cactusParams_get_int(params, 2, "bar", "pruneOutStubAlignments");

    //////////////////////////////////////////////
    //Run the bar algorithm
    //////////////////////////////////////////////

    if (listOfEndAlignmentFiles != NULL && stList_length(flowers) != 1) {
        st_errAbort("We have precomputed alignments but %" PRIi64 " flowers to align.\n", stList_length(flowers));
    }
    cactusDisk_preCacheStrings(cactusDisk, flowers);

#if defined(_OPENMP)
#pragma omp parallel for
    for (int64_t j = 0; j<stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);

        // These are all variables used by the filter fns
        FilterArgs *fa = st_calloc(1, sizeof(FilterArgs));
        fa->minimumIngroupDegree = cactusParams_get_int(params, 2, "bar", "minimumIngroupDegree");
        fa->minimumOutgroupDegree = cactusParams_get_int(params, 2, "bar", "minimumOutgroupDegree");
        fa->minimumDegree = cactusParams_get_int(params, 2, "bar", "minimumBlockDegree");
        fa->minimumNumberOfSpecies = cactusParams_get_int(params, 2, "bar", "minimumNumberOfSpecies");
        fa->flower = flower;

        void *alignments;
        if (usePoa) {
            /*
             * This makes a consistent set of alignments using abPoa.
             *
             * It does not use any precomputed alignments, if they are provided they will be ignored
             */
            alignments = make_flower_alignment_poa(flower, maximumLength, poaWindow, maskFilter, poaBandConstant, poaBandFraction);
            st_logDebug("Created the poa alignments: %" PRIi64 " poa alignment blocks for flower\n", stList_length(alignments));
        } else {
            StateMachine *sM = stateMachine5_construct(fiveState);
            alignments = makeFlowerAlignment3(sM, flower, listOfEndAlignmentFiles, spanningTrees, maximumLength,
                                              useProgressiveMerging, matchGamma,
                                              pairwiseAlignmentParameters,
                                              pruneOutStubAlignments);
            stateMachine_destruct(sM);
            st_logDebug("Created the alignment: %" PRIi64 " pairs for flower\n", stSortedSet_size(alignments));
        }

        stPinchIterator *pinchIterator = NULL;
        if(usePoa) {
            pinchIterator = stPinchIterator_constructFromAlignedBlocks(alignments);
        }
        else {
            pinchIterator = stPinchIterator_constructFromAlignedPairs(alignments, getNextAlignedPairAlignment);
        }
        /*
         * Run the cactus caf functions to build cactus.
         */

        stPinchThreadSet *threadSet = stCaf_setup(flower);

        stCaf_anneal(threadSet, pinchIterator, NULL, flower);

        if (fa->minimumDegree < 2) {
            stCaf_makeDegreeOneBlocks(threadSet);
        }

        if (fa->minimumIngroupDegree > 0 || fa->minimumOutgroupDegree > 0 || fa->minimumDegree > 1) {
            stCaf_melt(flower, threadSet, blockFilterFn, fa, 0, 0, 0, INT64_MAX);
        }

        stCaf_finish(flower, threadSet, chainLengthForBigFlower, longChain, INT64_MAX, INT64_MAX, cleanupMemory); //Flower now destroyed.

        stPinchThreadSet_destruct(threadSet);
        st_logDebug("Ran the cactus core script.\n");

        /*
         * Cleanup
         */
        //Clean up the sorted set after cleaning up the iterator
        stPinchIterator_destruct(pinchIterator);
        if(poaWindow != 0) {
            stList_destruct(alignments);
        }
        else {
            stSortedSet_destruct(alignments);
        }
        free(fa);

        st_logDebug("Finished filling in the alignments for the flower\n");
    }

    //st_logDebug("Created the alignment: %" PRIi64 " pairs\n", stSortedSet_size(alignedPairs));
    //pinchIterator = stPinchIterator_constructFromAlignedPairs(alignedPairs, getNextAlignedPairAlignment);
    //}
#else
    assert(0);
    for (int64_t j = 0; j < stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);
        //fa->flower = flower;
        st_logDebug("Processing a flower\n");

        stPinchIterator *pinchIterator = NULL;
        stSortedSet *alignedPairs = NULL;
        stList *alignment_blocks = NULL;

        if(usePoa) {
            /*
             * This makes a consistent set of alignments using abPoa.
             *
             * It does not use any precomputed alignments, if they are provided they will be ignored
             */
            alignment_blocks = make_flower_alignment_poa(flower, maximumLength, poaWindow, maskFilter, poaBandConstant, poaBandFraction);
            st_logDebug("Created the poa alignments: %" PRIi64 " poa alignment blocks\n", stList_length(alignment_blocks));
            pinchIterator = stPinchIterator_constructFromAlignedBlocks(alignment_blocks);
        }
        else {
            StateMachine *sM = stateMachine5_construct(fiveState);
            alignedPairs = makeFlowerAlignment3(sM, flower, listOfEndAlignmentFiles, spanningTrees, maximumLength,
                                                useProgressiveMerging, matchGamma,
                                                pairwiseAlignmentParameters,
                                                pruneOutStubAlignments);
            stateMachine_destruct(sM);
            st_logDebug("Created the alignment: %" PRIi64 " pairs\n", stSortedSet_size(alignedPairs));
            pinchIterator = stPinchIterator_constructFromAlignedPairs(alignedPairs, getNextAlignedPairAlignment);
        }
        /*
         * Run the cactus caf functions to build cactus.
         */
        stPinchThreadSet *threadSet = stCaf_setup(flower);
        stCaf_anneal(threadSet, pinchIterator, NULL);
        if (fa->minimumDegree < 2) {
            stCaf_makeDegreeOneBlocks(threadSet);
        }
        if (fa->minimumIngroupDegree > 0 || fa->minimumOutgroupDegree > 0 || fa->minimumDegree > 1) {
            stCaf_melt(flower, threadSet, blockFilterFn, fa, 0, 0, 0, INT64_MAX);
        }

        stCaf_finish(flower, threadSet, chainLengthForBigFlower, longChain, INT64_MAX, INT64_MAX, cleanupMemory); //Flower now destroyed.
        stPinchThreadSet_destruct(threadSet);
        st_logDebug("Ran the cactus core script.\n");

        /*
         * Cleanup
         */
        //Clean up the sorted set after cleaning up the iterator
        stPinchIterator_destruct(pinchIterator);
        if(poaWindow != 0) {
            stList_destruct(alignment_blocks);
        }
        else {
            stSortedSet_destruct(alignedPairs);
        }

        st_logDebug("Finished filling in the alignments for the flower\n");
    }
#endif

    //////////////////////////////////////////////
    //Clean up
    //////////////////////////////////////////////

    //stList_destruct(flowers);
    pairwiseAlignmentBandingParameters_destruct(pairwiseAlignmentParameters);

    /*if (bedRegions != NULL) {
        // Clean up our mapping.
        munmap(bedRegions, numBeds * sizeof(bedRegion));
    }*/
}
