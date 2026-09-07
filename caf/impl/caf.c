#include "cactus.h"
#include "sonLib.h"
#include "stCaf.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stGiantComponent.h"
#include "stCafPhylogeny.h"
#include <string.h>
#include <time.h>
#if defined(_OPENMP)
#include <omp.h>
#endif

double stCaf_now(void) {
#if defined(_OPENMP)
    return omp_get_wtime();
#else
    return ((double) clock()) / CLOCKS_PER_SEC;
#endif
}

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
    return n <= 1 ? 0 : n * (n - 1) / 2;
}

// Get the number of possible pairwise alignments that could support
// this block. Ordinarily this is (degree choose 2), but since we
// don't do outgroup self-alignment, it's a bit smaller.
static uint64_t numPossibleSupportingHomologies(stPinchBlock *block, Flower *flower) {
    uint64_t outgroupDegree = 0, ingroupDegree = 0;
    stPinchBlockIt segIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&segIt)) != NULL) {
        Event *event = stCaf_getEvent(segment, flower);
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

/*
 * The k-th smallest (counting from 0) of n doubles, i.e. what an ascending sort would leave at index k,
 * found by quickselect in expected O(n). The array is reordered.
 */
static double selectKthSmallest(double *a, int64_t n, int64_t k) {
    int64_t left = 0, right = n - 1;
    while (right > left) {
        int64_t mid = left + (right - left) / 2;
        double t;
        if (a[mid] < a[left]) { t = a[mid]; a[mid] = a[left]; a[left] = t; }
        if (a[right] < a[left]) { t = a[right]; a[right] = a[left]; a[left] = t; }
        if (a[right] < a[mid]) { t = a[right]; a[right] = a[mid]; a[mid] = t; }
        double pivot = a[mid];
        int64_t i = left, j = right;
        while (i <= j) {
            while (a[i] < pivot) {
                i++;
            }
            while (a[j] > pivot) {
                j--;
            }
            if (i <= j) {
                t = a[i]; a[i] = a[j]; a[j] = t;
                i++;
                j--;
            }
        }
        if (k <= j) {
            right = j;
        } else if (k >= i) {
            left = i;
        } else {
            return a[k]; //everything between j and i equals the pivot
        }
    }
    return a[k];
}

// Print a set of statistics (avg, median, max, min) for degree and
// support percentage in the pinch graph.
static void printThreadSetStatistics(stPinchThreadSet *threadSet, Flower *flower, FILE *f)
{
    // One pass over the blocks. The degree median comes from a histogram and the support median from
    // an O(n) selection; sorting two arrays with an entry per block, as this used to, was most of the
    // cost of the diagnostic once the graph reached hundreds of millions of blocks.
    uint64_t numBlocks = 0;
    double totalDegree = 0.0, totalSupport = 0.0;
    uint64_t totalAlignedBases = 0;
    uint64_t minDegree = UINT64_MAX, maxDegree = 0;
    double minSupport = 0.0, maxSupport = 0.0;
    uint64_t *degreeHistogram = NULL;
    uint64_t degreeHistogramSize = 0;
    double *supports = NULL;
    uint64_t supportsCapacity = 0;

    stPinchThreadSetBlockIt it = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *block;
    while ((block = stPinchThreadSetBlockIt_getNext(&it)) != NULL) {
        uint64_t degree = stPinchBlock_getDegree(block);
        totalDegree += degree;
        if (degree >= degreeHistogramSize) {
            uint64_t newSize = degreeHistogramSize == 0 ? 256 : degreeHistogramSize;
            while (newSize <= degree) {
                newSize *= 2;
            }
            degreeHistogram = st_realloc(degreeHistogram, newSize * sizeof(uint64_t));
            memset(degreeHistogram + degreeHistogramSize, 0, (newSize - degreeHistogramSize) * sizeof(uint64_t));
            degreeHistogramSize = newSize;
        }
        degreeHistogram[degree]++;
        if (degree < minDegree) {
            minDegree = degree;
        }
        if (degree > maxDegree) {
            maxDegree = degree;
        }
        uint64_t supportingHomologies = stPinchBlock_getNumSupportingHomologies(block);
        uint64_t possibleSupportingHomologies = numPossibleSupportingHomologies(block, flower);
        double support = 0.0;
        if (possibleSupportingHomologies != 0) {
            support = ((double) supportingHomologies) / possibleSupportingHomologies;
        }
        if (numBlocks == supportsCapacity) {
            supportsCapacity = supportsCapacity == 0 ? 1024 : supportsCapacity * 2;
            supports = st_realloc(supports, supportsCapacity * sizeof(double));
        }
        supports[numBlocks] = support;
        totalSupport += support;
        if (numBlocks == 0 || support < minSupport) {
            minSupport = support;
        }
        if (numBlocks == 0 || support > maxSupport) {
            maxSupport = support;
        }

        totalAlignedBases += stPinchBlock_getLength(block) * degree;

        numBlocks++;
    }

    fprintf(f, "There were %" PRIu64 " blocks in the sequence graph, representing %" PRIi64
    " total aligned bases\n", numBlocks, totalAlignedBases);

    if (numBlocks > 0) {
        // The medians are what the ascending sorts used to leave at index (numBlocks - 1) / 2
        uint64_t medianIndex = (numBlocks - 1) / 2, seen = 0, medianDegree = 0;
        for (uint64_t degree = 0; degree < degreeHistogramSize; degree++) {
            seen += degreeHistogram[degree];
            if (seen > medianIndex) {
                medianDegree = degree;
                break;
            }
        }
        double medianSupport = selectKthSmallest(supports, numBlocks, medianIndex);
        fprintf(f, "Block degree stats: min %" PRIu64 ", avg %lf, median %" PRIu64 ", max %" PRIu64 "\n",
                minDegree, totalDegree/numBlocks, medianDegree, maxDegree);
        fprintf(f, "Block support stats: min %lf, avg %lf, median %lf, max %lf\n",
               minSupport, totalSupport/numBlocks, medianSupport, maxSupport);
    }
    free(degreeHistogram);
    free(supports);
}

void caf(Flower *flower, CactusParams *params, char *alignmentsFile, char *secondaryAlignmentsFile, char *constraintsFile,
         Event *referenceEvent) {
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

    // As these parameters depend on the ingroup subtree we parse it here
    stTree *tree = event_getStTree(referenceEvent);
    double max_path_distance = stTree_getLongestPathLength(tree); // This is the longest path distance between ingroups, we use
    // this distance to choose the min chain length (the annealingRounds parameter)
    stTree_destruct(tree); // Cleanup the tree

    // Pick the annealing round parameter based on the distance
    int64_t annealingRoundsLength;
    int64_t *annealingRounds = NULL;
    if (cactusParams_get_float(params, 3, "constants", "divergences", "useDefault") == 0) {
        if(max_path_distance < cactusParams_get_float(params, 3, "constants", "divergences", "one")) {
            annealingRounds = cactusParams_get_ints(params, &annealingRoundsLength, 3, "caf", "annealingRounds", "one");
        } else if(max_path_distance < cactusParams_get_float(params, 3, "constants", "divergences", "two")) {
            annealingRounds = cactusParams_get_ints(params, &annealingRoundsLength, 3, "caf", "annealingRounds", "two");
        } else if(max_path_distance < cactusParams_get_float(params, 3, "constants", "divergences", "three")) {
            annealingRounds = cactusParams_get_ints(params, &annealingRoundsLength, 3, "caf", "annealingRounds", "three");
        } else if(max_path_distance < cactusParams_get_float(params, 3, "constants", "divergences", "four")) {
            annealingRounds = cactusParams_get_ints(params, &annealingRoundsLength, 3, "caf", "annealingRounds", "four");
        } else if(max_path_distance < cactusParams_get_float(params, 3, "constants", "divergences", "five")) {
            annealingRounds = cactusParams_get_ints(params, &annealingRoundsLength, 3, "caf", "annealingRounds", "five");
        }
    }
    if (annealingRounds == NULL) {
        annealingRounds = cactusParams_get_ints(params, &annealingRoundsLength, 3, "caf", "annealingRounds", "default");
    }

    // Log the annealing round parameters
    char *tree_string = eventTree_makeNewickString(flower_getEventTree(flower));
    st_logInfo("We found a max path distance between ingroups in the tree (%s) of %f, giving us and min final chain length of: %" PRIi64 "\n",
               tree_string, max_path_distance, annealingRounds[annealingRoundsLength-1]);
    free(tree_string);

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
    } else if (strcmp(alignmentFilter, "filterSecondariesByMultipleSequences") == 0) {
        sortAlignments = false;
        filterFn = NULL;
        secondaryFilterFn = stCaf_filterByMultipleSequences;
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
        double cafStartTime = stCaf_now(), t = cafStartTime;

        //Set up the graph and add the initial alignments
        stPinchThreadSet *threadSet = stCaf_setup(flower);

        //Build the set of outgroup threads
        stSet *outgroupThreads = stCaf_getOutgroupThreads(flower, threadSet);
        st_logInfo("caf-timing: setup %.3fs\n", stCaf_now() - t);

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
            //stCaf_sortCigarsFileByScoreInDescendingOrder(alignmentsFile, tempFile1);
            pinchIterator = stPinchIterator_constructFromFile(tempFile1);
        } else {
            pinchIterator = stPinchIterator_constructFromFile(alignmentsFile);
        }

        if(secondaryAlignmentsFile != NULL) {
            if (sortSecondaryAlignments) {
                tempFile2 = getTempFile();
                //stCaf_sortCigarsFileByScoreInDescendingOrder(secondaryAlignmentsFile, tempFile2);
                secondaryPinchIterator = stPinchIterator_constructFromFile(tempFile2);
            } else {
                secondaryPinchIterator = stPinchIterator_constructFromFile(secondaryAlignmentsFile);
            }
        }

        for (int64_t annealingRound = 0; annealingRound < annealingRoundsLength; annealingRound++) {
            int64_t minimumChainLength = annealingRounds[annealingRound];
            int64_t alignmentTrim = annealingRound < alignmentTrimLength ? alignmentTrims[annealingRound] : 0;
            st_logInfo("Starting annealing round with a minimum chain length of %" PRIi64 " and an alignment trim of %" PRIi64 "\n", minimumChainLength, alignmentTrim);

            stPinchIterator_setTrim(pinchIterator, alignmentTrim);
            if(secondaryPinchIterator != NULL) {
                stPinchIterator_setTrim(secondaryPinchIterator, alignmentTrim);
            }

            //Add back in the constraints
            if (pinchIteratorForConstraints != NULL) {
                stCaf_anneal(threadSet, pinchIteratorForConstraints, NULL, flower);
            }

            //Do the annealing
            t = stCaf_now();
            if (annealingRound == 0) {
                stCaf_anneal(threadSet, pinchIterator, filterFn, flower);
            } else {
                stCaf_annealBetweenAdjacencyComponents(threadSet, pinchIterator, filterFn, flower);
            }
            double primaryAnnealTime = stCaf_now() - t;

            // Do the secondary annealing
            t = stCaf_now();
            if(secondaryPinchIterator != NULL) {
                if (annealingRound == 0) {
                    stCaf_anneal(threadSet, secondaryPinchIterator, secondaryFilterFn, flower);
                } else {
                    stCaf_annealBetweenAdjacencyComponents(threadSet, secondaryPinchIterator, secondaryFilterFn, flower);
                }
            }
            double secondaryAnnealTime = stCaf_now() - t;

            t = stCaf_now();
            st_logInfo("Sequence graph statistics after annealing:\n");
            printThreadSetStatistics(threadSet, flower, stderr);
            double statsTime = stCaf_now() - t;
            t = stCaf_now();

            // The support test below is only ever true for blocks with a degree above
            // minimumBlockDegreeToCheckSupport, so with that unset there is nothing to look for
            if (minimumBlockHomologySupport > 0 && minimumBlockDegreeToCheckSupport > 0) {
                // Check for poorly-supported blocks--those that have
                // been transitively aligned together but with very
                // few homologies supporting the transitive
                // alignment. These "megablocks" can snarl up the
                // graph so that a lot of extra gets thrown away in
                // the first melting step.
                stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
                stPinchBlock *block;
                int64_t num_megablocks_destroyed = 0;
                int64_t num_homologies_destroyed = 0;
                while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
                    if (minimumBlockDegreeToCheckSupport > 0 && stPinchBlock_getDegree(block) > minimumBlockDegreeToCheckSupport) {
                        uint64_t supportingHomologies = stPinchBlock_getNumSupportingHomologies(block);
                        uint64_t possibleSupportingHomologies = numPossibleSupportingHomologies(block, flower);
                        double support = ((double) supportingHomologies) / possibleSupportingHomologies;
                        if (support < minimumBlockHomologySupport) {
                            st_logDebug("Destroyed a megablock with degree %" PRIi64
                            " and %" PRIi64 " supporting homologies out of a maximum "
                                            "of %" PRIi64 " (%lf%%).\n", stPinchBlock_getDegree(block),
                                    supportingHomologies, possibleSupportingHomologies, support);
                            stPinchBlock_destruct(block);
                            ++num_megablocks_destroyed;
                            num_homologies_destroyed += supportingHomologies;
                        }
                    }
                }
                if (num_megablocks_destroyed > 0) {
                  st_logInfo("Destroyed %" PRIi64 " megablocks with a total of %" PRIi64 " supporting homologies\n",
                             num_megablocks_destroyed, num_homologies_destroyed);
                }
            }
            st_logInfo("caf-timing: anneal round %" PRIi64 " minChain=%" PRIi64 " primary %.3fs secondary %.3fs stats %.3fs megablocks %.3fs\n",
                       annealingRound, minimumChainLength, primaryAnnealTime, secondaryAnnealTime, statsTime, stCaf_now() - t);

            //Do the melting rounds
            for (int64_t meltingRound = 0; meltingRound < meltingRoundsLength; meltingRound++) {
                int64_t minimumChainLengthForMeltingRound = meltingRounds[meltingRound];
                st_logInfo("Starting melting round with a minimum chain length of %" PRIi64 " \n", minimumChainLengthForMeltingRound);
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
            t = stCaf_now();
            stCaf_meltRecoverableChains(flower, threadSet, breakChainsAtReverseTandems, maximumMedianSequenceLengthBetweenLinkedEnds, recoverableChainsFilter, maxRecoverableChainsIterations, maxRecoverableChainLength);
            st_logInfo("caf-timing: recoverable total %.3fs\n", stCaf_now() - t);
        }

        t = stCaf_now();
        st_logInfo("Sequence graph statistics after melting:\n");
        printThreadSetStatistics(threadSet, flower, stderr);
        st_logInfo("caf-timing: stats %.3fs\n", stCaf_now() - t);

        //Sort out case when we allow blocks of degree 1
        t = stCaf_now();
        if (fa->minimumDegree < 2) {
            st_logDebug("Creating degree 1 blocks\n");
            stCaf_makeDegreeOneBlocks(threadSet);
            stCaf_melt(flower, threadSet, blockFilterFn, fa, blockTrim, 0, 0, INT64_MAX);
            st_logInfo("caf-timing: degree-one %.3fs\n", stCaf_now() - t);
        } else if (maximumAdjacencyComponentSizeRatio < INT64_MAX) { //Deal with giant components
            st_logDebug("Breaking up components greedily\n");
            stCaf_breakupComponentsGreedily(threadSet, maximumAdjacencyComponentSizeRatio);
            st_logInfo("caf-timing: breakup-components %.3fs\n", stCaf_now() - t);
        }

        //Finish up
        t = stCaf_now();
        stCaf_finish(flower, threadSet, minLengthForChromosome, proportionOfUnalignedBasesForNewChromosome);
        st_logInfo("caf-timing: finish %.3fs\n", stCaf_now() - t);
        st_logInfo("caf-timing: caf total %.3fs\n", stCaf_now() - cafStartTime);
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
