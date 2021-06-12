#include "cactus.h"
#include "sonLib.h"
#include "stCaf.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stLastzAlignments.h"
#include "stGiantComponent.h"
#include "stCafPhylogeny.h"

static bool blockFilterFn(stPinchBlock *pinchBlock, void *extraArg) {
    FilterArgs *f = extraArg;
    if (!stCaf_containsRequiredSpecies(pinchBlock, f->flower, f->minimumIngroupDegree,
                                       f->minimumOutgroupDegree, f->minimumDegree,
                                       f->minimumNumberOfSpecies)) {
        return 1;
    }
    if (f->minimumTreeCoverage > 0.0 && stCaf_treeCoverage(pinchBlock, f->flower) < f->minimumTreeCoverage) { //Tree coverage
        return 1;
    }
    return 0;
}

static uint64_t choose2(uint64_t n) {
#define CHOOSE_TWO_CACHE_LEN 256
    int64_t chooseTwoCache[256] = {0}; // removing static variable, even if fast
    if (n <= 1) {
        return 0;
    } else if (n >= CHOOSE_TWO_CACHE_LEN) {
        return n * (n - 1) / 2;
    } else {
        if (chooseTwoCache[n] != 0) {
            return chooseTwoCache[n];
        } else {
            chooseTwoCache[n] = n * (n - 1) / 2;
            return chooseTwoCache[n];
        }
    }
}

// Get the number of possible pairwise alignments that could support
// this block. Ordinarily this is (degree choose 2), but since we
// don't do outgroup self-alignment, it's a bit smaller.
static uint64_t numPossibleSupportingHomologies(stPinchBlock *block, Flower *flower) {
    uint64_t outgroupDegree = 0, ingroupDegree = 0;
    stPinchBlockIt segIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&segIt)) != NULL) {
        Name capName = stPinchSegment_getName(segment);
        Cap *cap = flower_getCap(flower, capName);
        Event *event = cap_getEvent(cap);
        if (event_isOutgroup(event)) {
            outgroupDegree++;
        } else {
            ingroupDegree++;
        }
    }
    assert(outgroupDegree + ingroupDegree == stPinchBlock_getDegree(block));
    // We do the ingroup-ingroup alignments as an all-against-all
    // alignment, so we can see each ingroup-ingroup homology up to
    // twice.
    return choose2(ingroupDegree) * 2 + ingroupDegree * outgroupDegree;
}

// for printThreadSetStatistics
static int uint64_cmp(const uint64_t *x, const uint64_t *y) {
    if (*x < *y) {
        return -1;
    } else if (*x == *y) {
        return 0;
    } else {
        return 1;
    }
}

// for printThreadSetStatistics
static int double_cmp(const double *x, const double *y) {
    if (*x < *y) {
        return -1;
    } else if (*x == *y) {
        return 0;
    } else {
        return 1;
    }
}

// Print a set of statistics (avg, median, max, min) for degree and
// support percentage in the pinch graph.
static void printThreadSetStatistics(stPinchThreadSet *threadSet, Flower *flower, FILE *f)
{
    // Naively finds the max, median, and min by sorting: the lists
    // will have "only" millions of elements, so they should fit
    // comfortably into tens of MB of memory.

    uint64_t numBlocks = stPinchThreadSet_getTotalBlockNumber(threadSet);

    uint64_t *blockDegrees = malloc(numBlocks * sizeof(uint64_t));
    double totalDegree = 0.0;
    double *blockSupports = malloc(numBlocks * sizeof(double));
    double totalSupport = 0.0;

    uint64_t totalAlignedBases = 0;

    stPinchThreadSetBlockIt it = stPinchThreadSet_getBlockIt(threadSet);
    uint64_t i = 0;
    stPinchBlock *block;
    while ((block = stPinchThreadSetBlockIt_getNext(&it)) != NULL) {
        blockDegrees[i] = stPinchBlock_getDegree(block);
        totalDegree += stPinchBlock_getDegree(block);
        uint64_t supportingHomologies = stPinchBlock_getNumSupportingHomologies(block);
        uint64_t possibleSupportingHomologies = numPossibleSupportingHomologies(block, flower);
        double support = 0.0;
        if (possibleSupportingHomologies != 0) {
            support = ((double) supportingHomologies) / possibleSupportingHomologies;
        }
        blockSupports[i] = support;
        totalSupport += support;

        totalAlignedBases += stPinchBlock_getLength(block) * stPinchBlock_getDegree(block);

        i++;
    }

    fprintf(f, "There were %" PRIu64 " blocks in the sequence graph, representing %" PRIi64
    " total aligned bases\n", numBlocks, totalAlignedBases);

    qsort(blockDegrees, numBlocks, sizeof(uint64_t),
          (int (*)(const void *, const void *)) uint64_cmp);
    qsort(blockSupports, numBlocks, sizeof(double),
          (int (*)(const void *, const void *)) double_cmp);
    fprintf(f, "Block degree stats: min %" PRIu64 ", avg %lf, median %" PRIu64 ", max %" PRIu64 "\n",
            blockDegrees[0], totalDegree/numBlocks, blockDegrees[(numBlocks - 1) / 2],
            blockDegrees[numBlocks - 1]);
    fprintf(f, "Block support stats: min %lf, avg %lf, median %lf, max %lf\n",
           blockSupports[0], totalSupport/numBlocks, blockSupports[(numBlocks - 1) / 2],
           blockSupports[numBlocks - 1]);
    free(blockDegrees);
    free(blockSupports);
}

void caf(Flower *flower, CactusParams *params, char *alignmentsFile, char *secondaryAlignmentsFile, char *constraintsFile) {
    //////////////////////////////////////////////
    //Parse the many, many necessary parameters from the params file
    //////////////////////////////////////////////

    // Fixed
    bool breakChainsAtReverseTandems = 1;

    // These are all variables used by the filter fns
    FilterArgs *fa = st_malloc(sizeof(FilterArgs));
    fa->flower = flower;
    fa->minimumIngroupDegree = cactusParams_get_int(params, 2, "caf", "minimumIngroupDegree");
    fa->minimumOutgroupDegree = cactusParams_get_int(params, 2, "caf", "minimumOutgroupDegree");
    fa->minimumDegree = cactusParams_get_int(params, 2, "caf", "minimumBlockDegree");
    fa->minimumNumberOfSpecies = cactusParams_get_int(params, 2, "caf", "minimumNumberOfSpecies");
    fa->minimumTreeCoverage = cactusParams_get_float(params, 2, "caf", "minimumTreeCoverage");

    //Parameters for annealing/melting rounds
    int64_t annealingRoundsLength;
    int64_t *annealingRounds = cactusParams_get_ints(params, &annealingRoundsLength, 2, "caf", "annealingRounds");

    int64_t meltingRoundsLength;
    int64_t *meltingRounds = cactusParams_get_ints(params, &meltingRoundsLength, 2, "caf", "deannealingRounds");

    //Parameters for melting
    float maximumAdjacencyComponentSizeRatio = cactusParams_get_int(params, 2, "caf", "maxAdjacencyComponentSizeRatio");
    int64_t blockTrim = cactusParams_get_int(params, 2, "caf", "blockTrim");

    int64_t alignmentTrimLength = 0;
    int64_t *alignmentTrims = cactusParams_get_ints(params, &alignmentTrimLength, 2, "caf", "trim");

    int64_t minLengthForChromosome = cactusParams_get_int(params, 2, "caf", "minLengthForChromosome");
    float proportionOfUnalignedBasesForNewChromosome = cactusParams_get_float(params, 2, "caf", "proportionOfUnalignedBasesForNewChromosome");
    int64_t maximumMedianSequenceLengthBetweenLinkedEnds = cactusParams_get_int(params, 2, "caf", "maximumMedianSequenceLengthBetweenLinkedEnds");

    char *removeRecoverableChainsStr = (char *)cactusParams_get_string(params, 2, "caf", "removeRecoverableChains");
    bool removeRecoverableChains = false;
    bool (*recoverableChainsFilter)(stCactusEdgeEnd *, Flower *) = NULL;
    if (strcmp(removeRecoverableChainsStr, "1") == 0) {
        removeRecoverableChains = true;
        recoverableChainsFilter = NULL;
    } else if (strcmp(removeRecoverableChainsStr, "unequalNumberOfIngroupCopies") == 0) {
        removeRecoverableChains = true;
        recoverableChainsFilter = stCaf_chainHasUnequalNumberOfIngroupCopies;
    } else if (strcmp(removeRecoverableChainsStr, "unequalNumberOfIngroupCopiesOrNoOutgroup") == 0) {
        removeRecoverableChains = true;
        recoverableChainsFilter = stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup;
    } else if (strcmp(removeRecoverableChainsStr, "0") == 0) {
        removeRecoverableChains = false;
    } else {
        st_errAbort("Could not parse removeRecoverableChains argument");
    }
    free(removeRecoverableChainsStr);

    int64_t maxRecoverableChainsIterations = cactusParams_get_int(params, 2, "caf", "maxRecoverableChainsIterations");
    int64_t maxRecoverableChainLength = cactusParams_get_int(params, 2, "caf", "maxRecoverableChainLength");

    int64_t minimumBlockDegreeToCheckSupport = cactusParams_get_int(params, 2, "caf", "minimumBlockDegreeToCheckSupport");
    double minimumBlockHomologySupport = cactusParams_get_float(params, 2, "caf", "minimumBlockHomologySupport");

    // Setting the alignment filters
    char *alignmentFilter = (char *)cactusParams_get_string(params, 2, "caf", "alignmentFilter");
    bool sortAlignments = false;
    bool (*filterFn)(stPinchSegment *, stPinchSegment *, Flower *) = NULL;
    bool (*secondaryFilterFn)(stPinchSegment *, stPinchSegment *, Flower *) = NULL;
    char * singleCopyEventName = NULL;
    bool sortSecondaryAlignments = false;
    char *hgvmEventName = NULL;
    if (strcmp(alignmentFilter, "singleCopyOutgroup") == 0) {
        sortAlignments = true;
        filterFn = stCaf_filterByOutgroup;
    } else if (strcmp(alignmentFilter, "filterSecondariesByMultipleSpecies") == 0) {
        sortAlignments = false;
        filterFn = NULL;
        secondaryFilterFn = stCaf_filterByMultipleSpecies;
    } else if (strcmp(alignmentFilter, "relaxedSingleCopyOutgroup") == 0) {
        sortAlignments = true;
        filterFn = stCaf_relaxedFilterByOutgroup;
    } else if (strcmp(alignmentFilter, "singleCopy") == 0) {
        sortAlignments = true;
        filterFn = stCaf_filterByRepeatSpecies;
    } else if (strcmp(alignmentFilter, "relaxedSingleCopy") == 0) {
        sortAlignments = true;
        filterFn = stCaf_relaxedFilterByRepeatSpecies;
    } else if (strncmp(alignmentFilter, "singleCopyEvent:", 16) == 0) {
        singleCopyEventName = stString_copy(alignmentFilter + 16);
        filterFn = stCaf_filterBySingleCopyEvent;
    } else if (strcmp(alignmentFilter, "singleCopyChr") == 0) {
        sortAlignments = true;
        filterFn = stCaf_singleCopyChr;
    } else if (strcmp(alignmentFilter, "singleCopyIngroup") == 0) {
        sortAlignments = true;
        filterFn = stCaf_singleCopyIngroup;
    } else if (strcmp(alignmentFilter, "relaxedSingleCopyIngroup") == 0) {
        sortAlignments = true;
        filterFn = stCaf_relaxedSingleCopyIngroup;
    } else if (strncmp(alignmentFilter, "hgvm:", 5) == 0) {
        sortAlignments = true;
        size_t argLen = strlen(alignmentFilter);
        if (argLen < 6) {
            st_errAbort("alignmentFilter option \"hgvm\" needs an additional argument: "
                        "the event name to filter on. E.g. \"hgvm:human\"");
        }
        hgvmEventName = stString_copy(alignmentFilter + 5);
        filterFn = stCaf_filterToEnsureCycleFreeIsolatedComponents;
    } else if (strcmp(alignmentFilter, "none") == 0) {
        sortAlignments = false;
        filterFn = NULL;
    } else {
        st_errAbort("Could not recognize alignmentFilter option %s", alignmentFilter);
    }
    free(alignmentFilter);
    // by default we apply all primary filtering to secondary alignments too
    if (secondaryFilterFn == NULL && filterFn != NULL) {
        secondaryFilterFn = filterFn;
        sortSecondaryAlignments = sortAlignments;
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    //TODO: Add more here.
    assert(fa->minimumTreeCoverage >= 0.0);
    assert(fa->minimumTreeCoverage <= 1.0);
    assert(blockTrim >= 0);
    assert(annealingRoundsLength >= 0);
    for (int64_t i = 0; i < annealingRoundsLength; i++) {
        assert(annealingRounds[i] >= 0);
    }
    assert(meltingRoundsLength >= 0);
    for (int64_t i = 1; i < meltingRoundsLength; i++) {
        assert(meltingRounds[i - 1] < meltingRounds[i]);
        assert(meltingRounds[i - 1] >= 1);
    }
    assert(alignmentTrimLength >= 0);
    for (int64_t i = 0; i < alignmentTrimLength; i++) {
        assert(alignmentTrims[i] >= 0);
    }
    assert(fa->minimumOutgroupDegree >= 0);
    assert(fa->minimumIngroupDegree >= 0);

    ///////////////////////////////////////////////////////////////////////////
    // Get the constraints
    ///////////////////////////////////////////////////////////////////////////

    stPinchIterator *pinchIteratorForConstraints = NULL;
    if (constraintsFile != NULL) {
        pinchIteratorForConstraints = stPinchIterator_constructFromFile(constraintsFile);
        st_logDebug("Created an iterator for the alignment constaints from file: %s\n", constraintsFile);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Do the alignment
    ///////////////////////////////////////////////////////////////////////////

    char *tempFile1 = NULL;
    char *tempFile2 = NULL;

    if (!flower_builtBlocks(flower)) { // Do nothing if the flower already has defined blocks
        st_logDebug("Processing flower: %lli\n", flower_getName(flower));

        //Set up the graph and add the initial alignments
        stPinchThreadSet *threadSet = stCaf_setup(flower);

        //Build the set of outgroup threads
        stSet *outgroupThreads = stCaf_getOutgroupThreads(flower, threadSet);

        // Set the single copy event
        if (singleCopyEventName != NULL) {
            stCaf_setSingleCopyEvent(flower, singleCopyEventName);
        }

        if (filterFn == stCaf_filterToEnsureCycleFreeIsolatedComponents) {
            stCaf_setupHGVMFiltering(flower, threadSet, hgvmEventName);
        }

        //Setup the alignments
        stPinchIterator *pinchIterator = NULL;
        stPinchIterator *secondaryPinchIterator = NULL;
        stList *alignmentsList = NULL;
        assert(alignmentsFile != NULL);

        if (sortAlignments) {
            tempFile1 = getTempFile();
            stCaf_sortCigarsFileByScoreInDescendingOrder(alignmentsFile, tempFile1);
            pinchIterator = stPinchIterator_constructFromFile(tempFile1);
        } else {
            pinchIterator = stPinchIterator_constructFromFile(alignmentsFile);
        }

        if(secondaryAlignmentsFile != NULL) {
            if (sortSecondaryAlignments) {
                tempFile2 = getTempFile();
                stCaf_sortCigarsFileByScoreInDescendingOrder(secondaryAlignmentsFile, tempFile2);
                secondaryPinchIterator = stPinchIterator_constructFromFile(tempFile2);
            } else {
                secondaryPinchIterator = stPinchIterator_constructFromFile(secondaryAlignmentsFile);
            }
        }

        for (int64_t annealingRound = 0; annealingRound < annealingRoundsLength; annealingRound++) {
            int64_t minimumChainLength = annealingRounds[annealingRound];
            int64_t alignmentTrim = annealingRound < alignmentTrimLength ? alignmentTrims[annealingRound] : 0;
            st_logDebug("Starting annealing round with a minimum chain length of %" PRIi64 " and an alignment trim of %" PRIi64 "\n", minimumChainLength, alignmentTrim);

            stPinchIterator_setTrim(pinchIterator, alignmentTrim);
            if(secondaryPinchIterator != NULL) {
                stPinchIterator_setTrim(secondaryPinchIterator, alignmentTrim);
            }

            //Add back in the constraints
            if (pinchIteratorForConstraints != NULL) {
                stCaf_anneal(threadSet, pinchIteratorForConstraints, NULL, flower);
            }

            //Do the annealing
            if (annealingRound == 0) {
                stCaf_anneal(threadSet, pinchIterator, filterFn, flower);
            } else {
                stCaf_annealBetweenAdjacencyComponents(threadSet, pinchIterator, filterFn, flower);
            }

            // Do the secondary annealing
            if(secondaryPinchIterator != NULL) {
                if (annealingRound == 0) {
                    stCaf_anneal(threadSet, secondaryPinchIterator, secondaryFilterFn, flower);
                } else {
                    stCaf_annealBetweenAdjacencyComponents(threadSet, secondaryPinchIterator, secondaryFilterFn, flower);
                }
            }

            st_logDebug("Sequence graph statistics after annealing:\n");
            printThreadSetStatistics(threadSet, flower, stderr);

            if (minimumBlockHomologySupport > 0) {
                // Check for poorly-supported blocks--those that have
                // been transitively aligned together but with very
                // few homologies supporting the transitive
                // alignment. These "megablocks" can snarl up the
                // graph so that a lot of extra gets thrown away in
                // the first melting step.
                stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
                stPinchBlock *block;
                while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
                    if (stPinchBlock_getDegree(block) > minimumBlockDegreeToCheckSupport) {
                        uint64_t supportingHomologies = stPinchBlock_getNumSupportingHomologies(block);
                        uint64_t possibleSupportingHomologies = numPossibleSupportingHomologies(block, flower);
                        double support = ((double) supportingHomologies) / possibleSupportingHomologies;
                        if (support < minimumBlockHomologySupport) {
                            st_logDebug("Destroyed a megablock with degree %" PRIi64
                            " and %" PRIi64 " supporting homologies out of a maximum "
                                            "of %" PRIi64 " (%lf%%).\n", stPinchBlock_getDegree(block),
                                    supportingHomologies, possibleSupportingHomologies, support);
                            stPinchBlock_destruct(block);
                        }
                    }
                }
            }

            //Do the melting rounds
            for (int64_t meltingRound = 0; meltingRound < meltingRoundsLength; meltingRound++) {
                int64_t minimumChainLengthForMeltingRound = meltingRounds[meltingRound];
                st_logDebug("Starting melting round with a minimum chain length of %" PRIi64 " \n", minimumChainLengthForMeltingRound);
                if (minimumChainLengthForMeltingRound >= minimumChainLength) {
                    break;
                }
                stCaf_melt(flower, threadSet, NULL, NULL, 0, minimumChainLengthForMeltingRound, 0, INT64_MAX);
            } st_logDebug("Last melting round of cycle with a minimum chain length of %" PRIi64 " \n", minimumChainLength);
            stCaf_melt(flower, threadSet, NULL, NULL, 0, minimumChainLength, breakChainsAtReverseTandems, maximumMedianSequenceLengthBetweenLinkedEnds);
            //This does the filtering of blocks that do not have the required species/tree-coverage/degree.
            stCaf_melt(flower, threadSet, blockFilterFn, fa, blockTrim, 0, 0, INT64_MAX);
        }

        if (removeRecoverableChains) {
            stCaf_meltRecoverableChains(flower, threadSet, breakChainsAtReverseTandems, maximumMedianSequenceLengthBetweenLinkedEnds, recoverableChainsFilter, maxRecoverableChainsIterations, maxRecoverableChainLength);
        }

        st_logDebug("Sequence graph statistics after melting:\n");
        if(st_getLogLevel() == debug) {
            printThreadSetStatistics(threadSet, flower, stderr);
        }

        //Sort out case when we allow blocks of degree 1
        if (fa->minimumDegree < 2) {
            st_logDebug("Creating degree 1 blocks\n");
            stCaf_makeDegreeOneBlocks(threadSet);
            stCaf_melt(flower, threadSet, blockFilterFn, fa, blockTrim, 0, 0, INT64_MAX);
        } else if (maximumAdjacencyComponentSizeRatio < INT64_MAX) { //Deal with giant components
            st_logDebug("Breaking up components greedily\n");
            stCaf_breakupComponentsGreedily(threadSet, maximumAdjacencyComponentSizeRatio);
        }

        //Finish up
        stCaf_finish(flower, threadSet, minLengthForChromosome, proportionOfUnalignedBasesForNewChromosome);
        st_logDebug("Ran the cactus core script\n");

        //Cleanup
        stPinchThreadSet_destruct(threadSet);
        stPinchIterator_destruct(pinchIterator);
        if(secondaryPinchIterator != NULL) {
            stPinchIterator_destruct(secondaryPinchIterator);
        }
        stSet_destruct(outgroupThreads);

        if (alignmentsList != NULL) {
            stList_destruct(alignmentsList);
        }
        st_logDebug("Cleaned up from main loop\n");
    } else {
        st_logDebug("We've already built blocks / alignments for this flower\n");
    }

    // Cleanup
    free(annealingRounds);
    free(meltingRounds);
    free(alignmentTrims);
    free(fa);

    if (constraintsFile != NULL) {
        stPinchIterator_destruct(pinchIteratorForConstraints);
    }

    if (tempFile1 != NULL) {
        st_system("rm %s", tempFile1);
    }
    if (tempFile2 != NULL) {
        st_system("rm %s", tempFile2);
    }
}
