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

#include <sys/resource.h>
#include <time.h>

// How often to report progress through the flower list, in seconds.
#define BAR_PROGRESS_INTERVAL 600

PairwiseAlignmentParameters *pairwiseAlignmentParameters_constructFromCactusParams(CactusParams *params) {
    PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
    p->gapGamma = cactusParams_get_float(params, 3, "bar", "pecan", "gapGamma");
    p->splitMatrixBiggerThanThis = cactusParams_get_int(params, 3, "bar", "pecan", "splitMatrixBiggerThanThis");
    p->splitMatrixBiggerThanThis *= p->splitMatrixBiggerThanThis;
    p->anchorMatrixBiggerThanThis = cactusParams_get_int(params, 3, "bar", "pecan", "anchorMatrixBiggerThanThis");
    p->anchorMatrixBiggerThanThis *= p->anchorMatrixBiggerThanThis;
    p->repeatMaskMatrixBiggerThanThis = cactusParams_get_int(params, 3, "bar", "pecan", "repeatMaskMatrixBiggerThanThis");
    p->repeatMaskMatrixBiggerThanThis *= p->repeatMaskMatrixBiggerThanThis;
    p->diagonalExpansion = cactusParams_get_int(params, 3, "bar", "pecan", "diagonalExpansion");
    p->constraintDiagonalTrim = cactusParams_get_int(params, 3, "bar", "pecan", "constraintDiagonalTrim");
    p->alignAmbiguityCharacters = cactusParams_get_int(params, 3, "bar", "pecan", "alignAmbiguityCharacters");
    p->useMumAnchors = cactusParams_get_int(params, 3, "bar", "pecan", "useMumAnchors");
    p->recursiveMums = cactusParams_get_int(params, 3, "bar", "pecan", "recursiveMums");
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

void bar(stList *flowers, CactusParams *params, CactusDisk *cactusDisk, stList *listOfEndAlignmentFiles) {
    //////////////////////////////////////////////
    //Parse the many, many necessary parameters from the params file
    //////////////////////////////////////////////

    int64_t maximumLength = cactusParams_get_int(params, 2, "bar", "bandingLimit");
    int64_t usePoa = cactusParams_get_int(params, 2, "bar", "partialOrderAlignment");

    // Pecan prams
    int64_t spanningTrees = cactusParams_get_int(params, 3, "bar", "pecan", "spanningTrees");
    bool useProgressiveMerging = cactusParams_get_int(params, 3, "bar", "pecan", "useProgressiveMerging");
    float matchGamma = cactusParams_get_float(params, 3, "bar", "pecan", "matchGamma");
    PairwiseAlignmentParameters *pairwiseAlignmentParameters = pairwiseAlignmentParameters_constructFromCactusParams(params);
    StateMachine *sM = stateMachine5_construct(fiveState);
    bool pruneOutStubAlignments = cactusParams_get_int(params, 3, "bar", "pecan", "pruneOutStubAlignments");

    // Poa params
    // toggle from pecan to abpoa for multiple alignment, by setting to non-zero
    // Note that poa uses about N^2 memory, so maximum value is generally in 10s of kb
    int64_t poaWindow = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentWindow");
    int64_t maskFilter = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentMaskFilter");
    int64_t poaMaxProgRows = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentProgressiveMaxRows");
    double poaMaxLenDiff = cactusParams_get_float(params, 3, "bar", "poa", "partialOrderAlignmentProgressiveMaxLengthDiff");
    abpoa_para_t *poaParameters = usePoa ? abpoaParamaters_constructFromCactusParams(params) : NULL;

    //////////////////////////////////////////////
    //Run the bar algorithm
    //////////////////////////////////////////////

    if (listOfEndAlignmentFiles != NULL && stList_length(flowers) != 1) {
        st_errAbort("We have precomputed alignments but %" PRIi64 " flowers to align.\n", stList_length(flowers));
    }

    /*
     * Give each flower its own interval of names, so that the names of the blocks,
     * segments and so on that get built below depend on the flower's position in
     * the list and not on the order the threads happen to reach it.
     *
     * The bound: every base in the flower can end up in a pinch segment of its own,
     * each segment costs three names (itself and its two caps), and there is at most
     * one block per segment, at three names again.  Doubling that covers the ends,
     * groups, chains and nested flowers, which are all bounded by the number of
     * blocks.  A flower's bases are its unaligned length plus two per adjacency.
     *
     * The same pass totals the bases, for the progress report below.  The list is
     * sorted largest-first, so the fraction of flowers done is a poor guide to the
     * fraction of work done -- the first flowers are enormous and the last millions
     * are tiny -- and weighting by sequence tracks it far better.  The length is
     * already being read here, so the total costs nothing.
     */
    int64_t flowerNumber = stList_length(flowers);
    int64_t *nameOffsets = st_malloc(sizeof(int64_t) * (flowerNumber + 1));
    int64_t totalBases = 0;
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 1) reduction(+:totalBases)
#endif
    for (int64_t j = 0; j < flowerNumber; j++) {
        Flower *flower = stList_get(flowers, j);
        int64_t flowerBaseLength = flower_getTotalBaseLength(flower);
        nameOffsets[j] = 12 * (flowerBaseLength + flower_getCapNumber(flower)) + 4096;
        totalBases += flowerBaseLength;
    }
    int64_t nameTotal = 0; // Turn the sizes into offsets
    for (int64_t j = 0; j < flowerNumber; j++) {
        int64_t size = nameOffsets[j];
        nameOffsets[j] = nameTotal;
        nameTotal += size;
    }
    nameOffsets[flowerNumber] = nameTotal;
    Name nameBase = cactusDisk_reserveNames(cactusDisk, nameTotal);

    const bool reportProgress = st_getLogLevel() >= info && flowerNumber > 1;
    const time_t barStartTime = time(NULL);
    time_t lastReportTime = barStartTime;
    int64_t lastReportBases = 0;
    int64_t flowersDone = 0, basesDone = 0;

#if defined(_OPENMP)
    // Enable nested parallelism for large flowers with many ends
    // Level 0: serial (main thread)
    // Level 1: this parallel region (flower-level parallelism)
    // Level 2: nested parallel regions in make_consistent_partial_order_alignments
    omp_set_max_active_levels(2);
    
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int64_t j = 0; j<stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);

        cactusDisk_pushNameInterval(cactusDisk, nameBase + nameOffsets[j], nameOffsets[j+1] - nameOffsets[j]);
        st_randomSeed(flower_getName(flower)); // Any random choices below are the flower's own, not the thread's

        // Must be read before stCaf_finish, which adds block ends to the flower
        const int64_t flowerBases = reportProgress ? flower_getTotalBaseLength(flower) : 0;

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
            alignments = make_flower_alignment_poa(flower, maximumLength, poaWindow, maskFilter,
                                                   poaMaxProgRows, poaMaxLenDiff, poaParameters);
            st_logDebug("Created the poa alignments: %" PRIi64 " poa alignment blocks for flower\n", stList_length(alignments));
        } else {
            alignments = makeFlowerAlignment3(sM, flower, listOfEndAlignmentFiles, spanningTrees, maximumLength,
                                              useProgressiveMerging, matchGamma, pairwiseAlignmentParameters,
                                              pruneOutStubAlignments);
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

        stCaf_finish(flower, threadSet, INT64_MAX, INT64_MAX); //Flower now destroyed.

        stPinchThreadSet_destruct(threadSet);
        st_logDebug("Ran the cactus core script.\n");

        /*
         * Cleanup
         */
        //Clean up the sorted set after cleaning up the iterator
        stPinchIterator_destruct(pinchIterator);
        if(usePoa) {
            stList_destruct(alignments);
        }
        else {
            stSortedSet_destruct(alignments);
        }
        free(fa);

        cactusDisk_popNameInterval(cactusDisk);

        st_logDebug("Finished filling in the alignments for the flower\n");

        if (reportProgress) {
            int64_t done, bases;
#pragma omp atomic capture
            done = ++flowersDone;
#pragma omp atomic capture
            { basesDone += flowerBases; bases = basesDone; }

            time_t now = time(NULL), lastSeen;
#pragma omp atomic read
            lastSeen = lastReportTime;
            if (now - lastSeen >= BAR_PROGRESS_INTERVAL) {
#if defined(_OPENMP)
#pragma omp critical(barProgress)
#endif
                {
                    // re-check now that we hold the lock, so only one thread reports
                    if (now - lastReportTime >= BAR_PROGRESS_INTERVAL) {
                        // capture the window before advancing the marks
                        int64_t windowSeconds = now - lastReportTime;
                        int64_t windowBases = bases - lastReportBases;
                        lastReportTime = now;
                        lastReportBases = bases;
                        struct rusage ru;
                        getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
                        int64_t peakMemMB = ru.ru_maxrss / (1024 * 1024); // bytes on macOS
#else
                        int64_t peakMemMB = ru.ru_maxrss / 1024;          // kilobytes on Linux
#endif
                        int64_t elapsed = now - barStartTime;
                        double baseFraction = totalBases > 0 ? (double)bases / (double)totalBases : 0.0;
                        /*
                         * Estimate from the rate over the last interval, not from the average
                         * since the phase began.  The list is sorted largest-first, so the
                         * average is dominated by the huge flowers at the head and stays
                         * pessimistic for hours -- and while one of those is in flight nothing
                         * completes at all, which an average reads as "slow" rather than as
                         * "no information".  A zero-progress interval gives no estimate.
                         */
                        int64_t eta = (windowSeconds > 0 && windowBases > 0) ?
                            (int64_t)((double)(totalBases - bases) * windowSeconds / windowBases) : -1;
                        int64_t maxEnds, nestedFlowers, maxNestedThreads;
                        poa_get_nesting_stats(&maxEnds, &nestedFlowers, &maxNestedThreads);
                        st_logInfo("Bar progress: %" PRIi64 "/%" PRIi64 " flowers (%.2f%%), "
                                   "%" PRIi64 "/%" PRIi64 " bases (%.2f%%), %" PRIi64 " seconds in bar, "
                                   "eta %" PRIi64 " seconds, peak memory %" PRIi64 " MB, "
                                   "max ends %" PRIi64 ", nested on %" PRIi64 " flowers with up to %" PRIi64 " threads\n",
                                   done, flowerNumber, 100.0 * (double)done / (double)flowerNumber,
                                   bases, totalBases, 100.0 * baseFraction,
                                   elapsed, eta, peakMemMB,
                                   maxEnds, nestedFlowers, maxNestedThreads);
                    }
                }
            }
        }
    }
    free(nameOffsets);

    //////////////////////////////////////////////
    //Clean up
    //////////////////////////////////////////////

    pairwiseAlignmentBandingParameters_destruct(pairwiseAlignmentParameters);
    stateMachine_destruct(sM);

    if (poaParameters) {
        abpoa_free_para(poaParameters);
    }
}
