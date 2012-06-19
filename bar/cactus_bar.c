/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <getopt.h>

#include "cactus.h"
#include "sonLib.h"
#include "endAligner.h"
#include "flowerAligner.h"
#include "commonC.h"
#include "stCaf.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"

void usage() {
    fprintf(stderr, "cactus_baseAligner [flower-names, ordered by order they should be processed], version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --cactusDisk : The location of the flower disk directory\n");
    fprintf(
            stderr,
            "-i --spanningTrees (int >= 0) : The number of spanning trees construct in forming the set of pairwise alignments to include. If the number of pairwise alignments is less than the product of the total number of sequences and the number of spanning trees then all pairwise alignments will be included.\n");
    fprintf(
            stderr,
            "-j --maximumLength (int  >= 0 ) : The maximum length of a sequence to align, only the prefix and suffix maximum length bases are aligned\n");
    fprintf(stderr, "-k --useBanding : Use banding to speed up the alignments\n");
    fprintf(stderr, "-l --gapGamma : (float [0, 1]) The gap gamma (as in the AMAP function)\n");

    fprintf(stderr,
            "-o --splitMatrixBiggerThanThis : (int >= 0)  No dp matrix bigger than this number squared will be computed.\n");
    fprintf(
            stderr,
            "-p --anchorMatrixBiggerThanThis : (int >= 0)  Any matrix bigger than this number squared will be broken apart with banding.\n");
    fprintf(
            stderr,
            "-q --repeatMaskMatrixBiggerThanThis : (int >= 0) Any matrix bigger than this after initial banding will be broken apart without repeat masking of the sequences\n");
    fprintf(stderr, "-r --digaonalExpansion : (int >= 0 and even) Number of x-y diagonals to expand around anchors\n");
    fprintf(stderr, "-t --constraintDiagonalTrim : (int >= 0) Amount to trim from ends of each anchor\n");

    fprintf(stderr,
            "-u --minimumDegree : (int >= 0) Minimum number of sequences in a block to be included in the output graph\n");

    fprintf(stderr,
            "-w --alignAmbiguityCharacters : Align ambiguity characters (anything not ACTGactg) as a wildcard\n");

    fprintf(stderr,
            "-y --pruneOutStubAlignments : Prune out alignments of sequences that terminates in free stubs stubs\n");

    fprintf(stderr, "-n --numThreads : (int > 0) Maximum number of threads to use (OpenMP compilation required)\n");

    fprintf(stderr, "-A --requiredIngroupFraction : Fraction of ingroup events required in a block.\n");

    fprintf(stderr, "-B --requiredOutgroupFraction : Fraction of outgroup events required in a block.\n");

    fprintf(stderr, "-C --requiredAllFraction : Fraction of all events required in a block.\n");

    fprintf(stderr, "-h --help : Print this help screen\n");
}

stPinch *getNextAlignedPairAlignment(stSortedSetIterator *it) {
    AlignedPair *alignedPair = stSortedSet_getNext(it);
    if (alignedPair == NULL) {
        return NULL;
    }
    static stPinch pinch;
    stPinch_fillOut(&pinch, alignedPair->subsequenceIdentifier, alignedPair->reverse->subsequenceIdentifier,
            alignedPair->position, alignedPair->reverse->position, 1,
            alignedPair->strand == alignedPair->reverse->strand);
    return &pinch;
}

static int32_t requiredIngroupSpecies, requiredOutgroupSpecies, requiredAllSpecies;
static int32_t minimumDegree = 0;
static Flower *flower;

bool blockFilterFn(stPinchBlock *pinchBlock) {
    if (stPinchBlock_getDegree(pinchBlock) < minimumDegree) {
        return 1;
    }
    if (!stCaf_containsRequiredSpecies(pinchBlock, flower, requiredIngroupSpecies, requiredOutgroupSpecies,
            requiredAllSpecies)) {
        return 1;
    }
    return 0;
}

int main(int argc, char *argv[]) {

    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    int32_t i, j;
    int32_t spanningTrees = 10;
    int32_t maximumLength = 1500;
    float gapGamma = 0.5;
    bool useBanding = 0;
    int32_t numThreads = 1;
    int32_t k;

    PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters = pairwiseAlignmentBandingParameters_construct();

    /*
     * Setup the input parameters for cactus core.
     */
    bool pruneOutStubAlignments = 0;
    float requiredIngroupFraction = 0.0, requiredOutgroupFraction = 0.0, requiredAllFraction = 0.0;

    /*
     * Parse the options.
     */
    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "cactusDisk",
                required_argument, 0, 'b' }, { "help", no_argument, 0, 'h' }, { "spanningTrees", required_argument, 0,
                'i' }, { "maximumLength", required_argument, 0, 'j' }, { "useBanding", no_argument, 0, 'k' }, {
                "gapGamma", required_argument, 0, 'l' }, { "splitMatrixBiggerThanThis", required_argument, 0, 'o' }, {
                "anchorMatrixBiggerThanThis", required_argument, 0, 'p' }, { "repeatMaskMatrixBiggerThanThis",
                required_argument, 0, 'q' }, { "diagonalExpansion", required_argument, 0, 'r' }, {
                "constraintDiagonalTrim", required_argument, 0, 't' }, { "minimumDegree", required_argument, 0, 'u' },
                { "alignAmbiguityCharacters", no_argument, 0, 'w' }, { "pruneOutStubAlignments", no_argument, 0, 'y' },
                { "numThreads", required_argument, 0, 'n' }, { "requiredIngroupFraction", required_argument, 0, 'A' },
                { "requiredOutgroupFraction", required_argument, 0, 'B' }, { "requiredAllFraction", required_argument,
                        0, 'C' },

                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:hi:j:kl:o:p:q:r:t:u:wy:n:A:B:C:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                st_setLogLevelFromString(logLevelString);
                break;
            case 'b':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case 'i':
                i = sscanf(optarg, "%i", &spanningTrees);
                assert(i == 1);
                assert(spanningTrees >= 0);
                break;
            case 'j':
                i = sscanf(optarg, "%i", &maximumLength);
                assert(i == 1);
                assert(maximumLength >= 0);
                break;
            case 'k':
                useBanding = !useBanding;
                break;
            case 'l':
                i = sscanf(optarg, "%f", &gapGamma);
                assert(i == 1);
                assert(gapGamma >= 0.0);
                break;
            case 'o':
                i = sscanf(optarg, "%i", &k);
                assert(i == 1);
                assert(k >= 0);
                pairwiseAlignmentBandingParameters->splitMatrixBiggerThanThis = (int64_t) k * k;
                break;
            case 'p':
                i = sscanf(optarg, "%i", &k);
                assert(i == 1);
                assert(k >= 0);
                pairwiseAlignmentBandingParameters->anchorMatrixBiggerThanThis = (int64_t) k * k;
                break;
            case 'q':
                i = sscanf(optarg, "%i", &k);
                assert(i == 1);
                assert(k >= 0);
                pairwiseAlignmentBandingParameters->repeatMaskMatrixBiggerThanThis = (int64_t) k * k;
                break;
            case 'r':
                i = sscanf(optarg, "%i", &pairwiseAlignmentBandingParameters->diagonalExpansion);
                assert(i == 1);
                assert(pairwiseAlignmentBandingParameters->diagonalExpansion >= 0);
                assert(pairwiseAlignmentBandingParameters->diagonalExpansion % 2 == 0);
            case 't':
                i = sscanf(optarg, "%i", &pairwiseAlignmentBandingParameters->constraintDiagonalTrim);
                assert(i == 1);
                assert(pairwiseAlignmentBandingParameters->constraintDiagonalTrim >= 0);
                break;
            case 'u':
                i = sscanf(optarg, "%i", &minimumDegree);
                assert(i == 1);
                break;
            case 'w':
                pairwiseAlignmentBandingParameters->alignAmbiguityCharacters = 1;
                break;
            case 'y':
                pruneOutStubAlignments = 1;
                break;
            case 'n':
                i = sscanf(optarg, "%i", &numThreads);
                assert(i == 1);
                assert(numThreads > 0);
                break;
            case 'A':
                i = sscanf(optarg, "%f", &requiredIngroupFraction);
                assert(i == 1);
                break;
            case 'B':
                i = sscanf(optarg, "%f", &requiredOutgroupFraction);
                assert(i == 1);
                break;
            case 'C':
                i = sscanf(optarg, "%f", &requiredAllFraction);
                assert(i == 1);
                break;
            default:
                usage();
                return 1;
        }
    }

    st_setLogLevelFromString(logLevelString);

#ifdef _OPENMP
    omp_set_nested(0);
    omp_set_num_threads(numThreads);
#else
    if (numThreads > 1) {
        st_logCritical("numThreads option ignored because not compiled with -fopenmp");
    }
#endif

    /*
     * Load the flowerdisk
     */
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0); //We precache the sequences
    st_logInfo("Set up the flower disk\n");

    /*
     * For each flower.
     */
    stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    cactusDisk_preCacheStrings(cactusDisk, flowers);
    for (j = 0; j < stList_length(flowers); j++) {
        flower = stList_get(flowers, j);
        st_logInfo("Processing a flower\n");

        stSortedSet *alignedPairs = makeFlowerAlignment(flower, spanningTrees, maximumLength, gapGamma,
                pairwiseAlignmentBandingParameters, pruneOutStubAlignments);
        st_logInfo("Created the alignment: %i pairs\n", stSortedSet_size(alignedPairs));
        stPinchIterator *pinchIterator = stPinchIterator_constructFromAlignedPairs(alignedPairs,
                getNextAlignedPairAlignment);

        /*
         * Run the cactus caf functions to build cactus.
         */
        stPinchThreadSet *threadSet = stCaf_setup(flower);
        stCaf_anneal(threadSet, pinchIterator);
        if (minimumDegree < 2) {
            stCaf_makeDegreeOneBlocks(threadSet);
        }
        if (requiredIngroupFraction > 0.0 || requiredOutgroupFraction > 0.0 || requiredAllFraction > 0.0
                || minimumDegree > 1) {
            stCaf_calculateRequiredFractionsOfSpecies(flower,
                    requiredIngroupFraction, requiredOutgroupFraction, requiredAllFraction,
                    &requiredIngroupSpecies, &requiredOutgroupSpecies, &requiredAllSpecies);
            stCaf_melt(flower, threadSet, blockFilterFn, 0, 0);
        }
        stCaf_finish(flower, threadSet);
        stPinchThreadSet_destruct(threadSet);
        st_logInfo("Ran the cactus core script.\n");
        assert(!flower_isParentLoaded(flower));

        /*
         * Cleanup
         */
        //Clean up the sorted set after cleaning up the iterator
        stPinchIterator_destruct(pinchIterator);
        stSortedSet_destruct(alignedPairs);

        st_logInfo("Finished filling in the alignments for the flower\n");
    }
    stList_destruct(flowers);

    /*
     * Write and close the cactusdisk.
     */
    cactusDisk_write(cactusDisk);

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup
    ///////////////////////////////////////////////////////////////////////////

    return 0; //Exit without clean up is quicker, enable cleanup when doing memory leak detection.

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    //destructCactusCoreInputParameters(cCIP);
    free(cactusDiskDatabaseString);
    if (logLevelString != NULL) {
        free(logLevelString);
    }
    st_logInfo("Finished with the flower disk for this flower.\n");

    //while(1);

    return 0;
}

