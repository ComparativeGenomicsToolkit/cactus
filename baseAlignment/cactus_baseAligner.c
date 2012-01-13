/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <getopt.h>
#include <omp.h>

#include "cactus.h"
#include "sonLib.h"
#include "endAligner.h"
#include "flowerAligner.h"
#include "cactus_core.h"
#include "commonC.h"
#include "pairwiseAlignment.h"

void usage() {
    fprintf(
            stderr,
            "cactus_baseAligner [flower-names, ordered by order they should be processed], version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-b --cactusDisk : The location of the flower disk directory\n");
    fprintf(
            stderr,
            "-i --spanningTrees (int >= 0) : The number of spanning trees construct in forming the set of pairwise alignments to include. If the number of pairwise alignments is less than the product of the total number of sequences and the number of spanning trees then all pairwise alignments will be included.\n");
    fprintf(
            stderr,
            "-j --maximumLength (int  >= 0 ) : The maximum length of a sequence to align, only the prefix and suffix maximum length bases are aligned\n");
    fprintf(stderr,
            "-k --useBanding : Use banding to speed up the alignments\n");
    fprintf(stderr,
            "-l --gapGamma : (float [0, 1]) The gap gamma (as in the AMAP function)\n");

    fprintf(
            stderr,
            "-o --maxBandingSize : (int >= 0)  No dp matrix biggger than this number squared will be computed.\n");
    fprintf(
            stderr,
            "-p --minBandingSize : (int >= 0)  Any matrix bigger than this number squared will be broken apart with banding.\n");
    fprintf(
            stderr,
            "-q --minBandingConstraintDistance : (int >= 0)  The minimum size of a dp matrix between banding constraints.\n");
    fprintf(
            stderr,
            "-r --minTraceBackDiag : (int >= 0)  The x+y diagonal to leave between the cut point and the place we choose new cutpoints.\n");
    fprintf(
            stderr,
            "-s --minTraceGapDiags : (int >= 0)  The x+y diagonal distance to leave between a cutpoint and the traceback.\n");
    fprintf(
            stderr,
            "-t --constraintDiagonalTrim : (int >= 0)  The amount to be removed from each end of a diagonal to be considered a banding constraint.\n");

    fprintf(
            stderr,
            "-u --minimumDegree : (int >= 0) Minimum number of sequences in a block to be included in the output graph\n");

    fprintf(
            stderr,
            "-w --alignAmbiguityCharacters : Align ambiguity characters (anything not ACTGactg) as a wildcard\n");

    fprintf(
            stderr,
            "-y --pruneOutStubAlignments : Prune out alignments of sequences that terminates in free stubs stubs\n");

    fprintf(
            stderr,
            "-n --numThreads : (int > 0) Maximum number of threads to use (OpenMP compilation required)\n");

    fprintf(
            stderr,
            "-A --requiredIngroupFraction : Fraction of ingroup events required in a block.\n");

    fprintf(
            stderr,
            "-B --requiredOutgroupFraction : Fraction of outgroup events required in a block.\n");

    fprintf(stderr,
            "-C --requiredAllFraction : Fraction of all events required in a block.\n");

    fprintf(stderr, "-h --help : Print this help screen\n");
}

static stSortedSet *getAlignment_alignedPairs;
static stSortedSetIterator *getAlignment_iterator = NULL;

static struct PairwiseAlignment *getAlignments() {
    AlignedPair *alignedPair = stSortedSet_getNext(getAlignment_iterator);
    if (alignedPair == NULL) {
        return NULL;
    }
    struct List *opList = constructEmptyList(0,
            (void(*)(void *)) destructAlignmentOperation);
    listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, 1, 0.0));
    char *cA = cactusMisc_nameToString(alignedPair->sequence);
    char *cA2 = cactusMisc_nameToString(alignedPair->reverse->sequence);

    int32_t i = alignedPair->position;
    int32_t j = i + 1;
    if (!alignedPair->strand) {
        j = i;
        i = alignedPair->position + 1;
    }

    int32_t k = alignedPair->reverse->position;
    int32_t l = k + 1;
    if (!alignedPair->reverse->strand) {
        l = k;
        k = alignedPair->reverse->position + 1;
    }

    struct PairwiseAlignment *pairwiseAlignment = constructPairwiseAlignment(
            cA, i, j, alignedPair->strand, cA2, k, l,
            alignedPair->reverse->strand, 1.0, opList);
    free(cA);
    free(cA2);
    return pairwiseAlignment;
}

static void startAlignmentStack() {
    if (getAlignment_iterator != NULL) {
        stSortedSet_destructIterator(getAlignment_iterator);
        getAlignment_iterator = NULL;
    }
    getAlignment_iterator = stSortedSet_getIterator(getAlignment_alignedPairs);
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

    PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters =
            pairwiseAlignmentBandingParameters_construct();

    /*
     * Setup the input parameters for cactus core.
     */
    CactusCoreInputParameters *cCIP = constructCactusCoreInputParameters();
    free(cCIP->annealingRounds);
    cCIP->annealingRounds = st_malloc(sizeof(int32_t) * 1); //Set to one annealing round.
    cCIP->annealingRounds[0] = 0;
    cCIP->annealingRoundsLength = 1;
    cCIP->minimumDegree = 0;
    cCIP->minimumTreeCoverage = 0;
    bool pruneOutStubAlignments = 0;

    /*
     * Parse the options.
     */
    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'b' }, { "help", no_argument, 0, 'h' }, { "spanningTrees",
                required_argument, 0, 'i' }, { "maximumLength",
                required_argument, 0, 'j' }, { "useBanding", no_argument, 0,
                'k' }, { "gapGamma", required_argument, 0, 'l' }, {
                "maxBandingSize", required_argument, 0, 'o' }, {
                "minBandingSize", required_argument, 0, 'p' }, {
                "minBandingConstraintDistance", required_argument, 0, 'q' }, {
                "minTraceBackDiag", required_argument, 0, 'r' }, {
                "minTraceGapDiags", required_argument, 0, 's' }, {
                "constraintDiagonalTrim", required_argument, 0, 't' }, {
                "minimumDegree", required_argument, 0, 'u' }, {
                "alignAmbiguityCharacters", no_argument, 0, 'w' }, {
                "pruneOutStubAlignments", no_argument, 0, 'y' }, {
                "numThreads", required_argument, 0, 'n' }, {
                "requiredIngroupFraction", required_argument, 0, 'A' }, {
                "requiredOutgroupFraction", required_argument, 0, 'B' }, {
                "requiredAllFraction", required_argument, 0, 'C' },

        { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv,
                "a:b:hi:j:kl:o:p:q:r:s:t:u:wy:n:A:B:C:", long_options,
                &option_index);

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
                i = sscanf(optarg, "%i",
                        &(pairwiseAlignmentBandingParameters->maxBandingSize));
                assert(i == 1);
                assert(pairwiseAlignmentBandingParameters->maxBandingSize >= 0);
                break;
            case 'p':
                i = sscanf(optarg, "%i",
                        &pairwiseAlignmentBandingParameters->minBandingSize);
                assert(i == 1);
                assert(pairwiseAlignmentBandingParameters->minBandingSize >= 0);
                break;
            case 'q':
                i
                        = sscanf(
                                optarg,
                                "%i",
                                &pairwiseAlignmentBandingParameters->minBandingConstraintDistance);
                assert(i == 1);
                assert(
                        pairwiseAlignmentBandingParameters->minBandingConstraintDistance
                                >= 0);
                break;
            case 'r':
                i = sscanf(optarg, "%i",
                        &pairwiseAlignmentBandingParameters->minTraceBackDiag);
                assert(i == 1);
                assert(
                        pairwiseAlignmentBandingParameters->minTraceBackDiag
                                >= 0);
                break;
            case 's':
                i = sscanf(optarg, "%i",
                        &pairwiseAlignmentBandingParameters->minTraceGapDiags);
                assert(i == 1);
                assert(
                        pairwiseAlignmentBandingParameters->minTraceGapDiags
                                >= 0);
                break;
            case 't':
                i
                        = sscanf(
                                optarg,
                                "%i",
                                &pairwiseAlignmentBandingParameters->constraintDiagonalTrim);
                assert(i == 1);
                assert(
                        pairwiseAlignmentBandingParameters->constraintDiagonalTrim
                                >= 0);
                break;
            case 'u':
                i = sscanf(optarg, "%i", &cCIP->minimumDegree);
                assert(i == 1);
                break;
            case 'w':
                pairwiseAlignmentBandingParameters->alignAmbiguityCharacters
                        = 1;
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
                i = sscanf(optarg, "%f", &cCIP->requiredIngroupFraction);
                assert(i == 1);
                break;
            case 'B':
                i = sscanf(optarg, "%f", &cCIP->requiredOutgroupFraction);
                assert(i == 1);
                break;
            case 'C':
                i = sscanf(optarg, "%f", &cCIP->requiredAllFraction);
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
        st_logCritical(
                "numThreads option ignored because not compiled with -fopenmp");
    }
#endif

    /*
     * Load the flowerdisk
     */
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(
            cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0); //We precache the sequences
    st_logInfo("Set up the flower disk\n");

    /*
     * For each flower.
     */
    stList *flowers = parseFlowers(argv + optind, argc - optind, cactusDisk);
    for(j = 0; j < stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);
        st_logInfo("Processing a flower\n");

        getAlignment_alignedPairs = makeFlowerAlignment(flower, spanningTrees,
                maximumLength, gapGamma, useBanding,
                pairwiseAlignmentBandingParameters, pruneOutStubAlignments);
        st_logInfo("Created the alignment: %i pairs\n",
                stSortedSet_size(getAlignment_alignedPairs));
        //getAlignment_alignedPairs = stSortedSet_construct();

        /*
         * Run the cactus core script.
         */
        cactusCorePipeline(flower, cCIP, getAlignments, startAlignmentStack);
        st_logInfo("Ran the cactus core script.\n");

        /*
         * Cleanup
         */
        assert(getAlignment_iterator != NULL);
        stSortedSet_destructIterator(getAlignment_iterator);
        getAlignment_iterator = NULL;

        //Clean up the sorted set after cleaning up the iterator
        stSortedSet_destruct(getAlignment_alignedPairs);

        st_logInfo("Finished filling in the alignments for the flower\n");
    }

    ///////////////////////////////////////////////////////////////////////////
    // Unload the parent flowers
    ///////////////////////////////////////////////////////////////////////////

    for(j = 0; j < stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);
        assert(flower != NULL);
        //assert(!flower_isParentLoaded(flower));
        flower_unloadParent(flower); //We have this line just in case we are loading the parent..
    }
    stList_destruct(flowers);

    //assert(0);

    /*
     * Write and close the cactusdisk.
     */
    cactusDisk_write(cactusDisk);
    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    destructCactusCoreInputParameters(cCIP);
    free(cactusDiskDatabaseString);
    if (logLevelString != NULL) {
        free(logLevelString);
    }
    st_logInfo("Finished with the flower disk for this flower.\n");

    return 0;
}

