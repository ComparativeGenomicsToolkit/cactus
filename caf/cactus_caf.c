/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <getopt.h>
#include <time.h>

#include "cactus.h"
#include "sonLib.h"
#include "commonC.h"
#include "stCaf.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stLastzAlignments.h"
#include "stGiantComponent.h"

static void usage() {
    fprintf(stderr, "cactus_caf, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --alignments : The input alignments file\n");
    fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-d --lastzArguments : Lastz arguments\n");
    fprintf(stderr, "-h --help : Print this help screen\n");

    fprintf(stderr, "-i --annealingRounds (array of ints, each greater than or equal to 1) : The rounds of annealing\n");
    fprintf(stderr,
            "-o --deannealingRounds (array of ints, each greater than or equal to 1 and each greater than the last) : The rounds of deannealing\n");

    fprintf(
            stderr,
            "-k --trim (array of integers, each greater or equal to zero) : An array giving the trim for each annealing round. If the array is shorter than the annealing rounds then a trim value of 0 is assumed for annealing rounds greater than the length of the trim array\n");

    fprintf(stderr,
            "-m --minimumTreeCoverage : (float [0.0, 1.0]) Minimum tree coverage proportion of a block to be included in the graph\n");

    fprintf(
            stderr,
            "-n --blockTrim : (int >= 0) The number of bases to trim from the ends of each block in a chain before accepting, this filtering is done after choosing the length of chains\n");

    fprintf(stderr, "-p --minimumDegree : (int >= 0) Minimum number of sequences in a block to be included in the output graph\n");

    fprintf(stderr, "-q --minimumIngroupDegree : Number of ingroup sequences required in a block.\n");

    fprintf(stderr, "-r --minimumOutgroupDegree : Number of outgroup sequences required in a block.\n");

    fprintf(stderr, "-s --singleCopyIngroup : Require that in-group sequences have only single coverage\n");

    fprintf(stderr, "-t --singleCopyOutgroup : Require that out-group sequences have only single coverage\n");

    fprintf(stderr, "-v --minimumSequenceLengthForBlast : The minimum length of a sequence to include when blasting\n");

    fprintf(
            stderr,
            "-w --maxAdjacencyComponentSizeRatio : The components equal or less than log(n) * of this size will be allowed in the cactus. Used to fight giant components.\n");

    fprintf(stderr, "-x --constraints : A file of alignments that will be enforced upon the cactus\n");

    fprintf(stderr,
            "-y --minLengthForChromosome : The minimum length required for a sequence to be considered as a candidate to be chromosome.\n");

    fprintf(
            stderr,
            "-z --proportionOfUnalignedBasesForNewChromosome : Proportion of aligned bases to be not contained in an existing chromosome to cause generation of a new chromosome.\n");
    fprintf(
                stderr,
                "-A --maximumMedianSequenceLengthBetweenLinkedEnds : Maximum nedian length of sequences between linked ends to allow before breaking chains.\n");
    fprintf(stderr, "-B --realign : Realign the lastz hits.\n");
    fprintf(stderr, "-C --realignArguments : Arguments for realignment.\n");
}

static int64_t *getInts(const char *string, int64_t *arrayLength) {
    int64_t *iA = st_malloc(sizeof(int64_t) * strlen(string));
    char *cA = stString_copy(string);
    char *cA2 = cA;
    char *cA3;
    *arrayLength = 0;
    while ((cA3 = stString_getNextWord(&cA)) != NULL) {
        int64_t i = sscanf(cA3, "%" PRIi64 "", &iA[(*arrayLength)++]);
        (void) i;
        assert(i == 1);
        free(cA3);
    }
    free(cA2);
    return iA;
}

static int64_t minimumIngroupDegree = 0, minimumOutgroupDegree = 0, minimumDegree = 0, minimumNumberOfSpecies = 0;
static float minimumTreeCoverage = 0.0;
static Flower *flower = NULL;

static bool blockFilterFn(stPinchBlock *pinchBlock) {
    if (!stCaf_containsRequiredSpecies(pinchBlock, flower, minimumIngroupDegree,
                                       minimumOutgroupDegree, minimumDegree,
                                       minimumNumberOfSpecies)) {
        return 1;
    }
    if (minimumTreeCoverage > 0.0 && stCaf_treeCoverage(pinchBlock, flower) < minimumTreeCoverage) { //Tree coverage
        return 1;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    /*
     * Script for adding alignments to cactus tree.
     */
    int64_t startTime;
    stKVDatabaseConf *kvDatabaseConf;
    CactusDisk *cactusDisk;
    int key, k;

    bool (*filterFn)(stPinchSegment *, stPinchSegment *) = NULL;
    stSet *outgroupThreads = NULL;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * alignmentsFile = NULL;
    char * constraintsFile = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * lastzArguments = "";
    int64_t minimumSequenceLengthForBlast = 1;

    //Parameters for annealing/melting rounds
    int64_t *annealingRounds = NULL;
    int64_t annealingRoundsLength = 0;
    int64_t *meltingRounds = NULL;
    int64_t meltingRoundsLength = 0;

    //Parameters for melting
    float maximumAdjacencyComponentSizeRatio = 10;
    int64_t blockTrim = 0;
    int64_t alignmentTrimLength = 0;
    int64_t *alignmentTrims = NULL;
    bool singleCopyIngroup = 0;
    bool singleCopyOutgroup = 0;
    bool relaxedSingleCopyOutgroup = 0;
    int64_t chainLengthForBigFlower = 1000000;
    int64_t longChain = 2;
    int64_t minLengthForChromosome = 1000000;
    float proportionOfUnalignedBasesForNewChromosome = 0.8;
    bool breakChainsAtReverseTandems = 1;
    int64_t maximumMedianSequenceLengthBetweenLinkedEnds = INT64_MAX;
    bool realign = 0;
    char *realignArguments = "";
    bool removeRecoverableChains = false;
    bool (*recoverableChainsFilter)(stCactusEdgeEnd *, Flower *) = NULL;
    int64_t maxRecoverableChainsIterations = 1;
    int64_t maxRecoverableChainLength = INT64_MAX;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "alignments", required_argument, 0, 'b' }, {
                "cactusDisk", required_argument, 0, 'c' }, { "lastzArguments", required_argument, 0, 'd' },
                { "help", no_argument, 0, 'h' }, { "annealingRounds", required_argument, 0, 'i' }, { "trim", required_argument, 0, 'k' }, {
                        "trimChange", required_argument, 0, 'l', }, { "minimumTreeCoverage", required_argument, 0, 'm' }, { "blockTrim",
                        required_argument, 0, 'n' }, { "deannealingRounds", required_argument, 0, 'o' }, { "minimumDegree",
                        required_argument, 0, 'p' }, { "minimumIngroupDegree", required_argument, 0, 'q' }, {
                        "minimumOutgroupDegree", required_argument, 0, 'r' }, {
                        "singleCopyIngroup", no_argument, 0, 's' }, { "singleCopyOutgroup", required_argument, 0, 't' }, {
                        "minimumSequenceLengthForBlast", required_argument, 0, 'v' }, { "maxAdjacencyComponentSizeRatio",
                        required_argument, 0, 'w' }, { "constraints", required_argument, 0, 'x' }, { "minLengthForChromosome",
                        required_argument, 0, 'y' }, { "proportionOfUnalignedBasesForNewChromosome", required_argument, 0, 'z' },
                        { "maximumMedianSequenceLengthBetweenLinkedEnds", required_argument, 0, 'A' },
                        { "realign", no_argument, 0, 'B' }, { "realignArguments", required_argument, 0, 'C' },
                        { "removeRecoverableChains", required_argument, 0, 'D' },
                        { "minimumNumberOfSpecies", required_argument, 0, 'E' },
                        { "maxRecoverableChainsIterations", required_argument, 0, 'F' },
                        { "maxRecoverableChainLength", required_argument, 0, 'G' },
                        { 0, 0, 0, 0 } };

        int option_index = 0;

        key = getopt_long(argc, argv, "a:b:c:hi:k:m:n:o:p:q:r:stv:w:x:y:z:A:BC:D:E:F:G:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                st_setLogLevelFromString(logLevelString);
                break;
            case 'b':
                alignmentsFile = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'd':
                lastzArguments = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case 'i':
                annealingRounds = getInts(optarg, &annealingRoundsLength);
                break;
            case 'o':
                meltingRounds = getInts(optarg, &meltingRoundsLength);
                break;
            case 'k':
                alignmentTrims = getInts(optarg, &alignmentTrimLength);
                break;
            case 'm':
                k = sscanf(optarg, "%f", &minimumTreeCoverage);
                assert(k == 1);
                break;
            case 'n':
                k = sscanf(optarg, "%" PRIi64 "", &blockTrim);
                assert(k == 1);
                break;
            case 'p':
                k = sscanf(optarg, "%" PRIi64 "", &minimumDegree);
                assert(k == 1);
                break;
            case 'q':
                k = sscanf(optarg, "%" PRIi64 "", &minimumIngroupDegree);
                assert(k == 1);
                break;
            case 'r':
                k = sscanf(optarg, "%" PRIi64 "", &minimumOutgroupDegree);
                assert(k == 1);
                break;
            case 's':
                singleCopyIngroup = 1;
                break;
            case 't':
                if (strcmp(optarg, "1") == 0) {
                    singleCopyOutgroup = true;
                } else if (strcmp(optarg, "relaxed") == 0) {
                    singleCopyOutgroup = false;
                    relaxedSingleCopyOutgroup = true;
                } else if (strcmp(optarg, "0") == 0) {
                    singleCopyOutgroup = false;
                    relaxedSingleCopyOutgroup = false;
                } else {
                    st_errAbort("Could not recognize singleCopyOutgroup option %s", optarg);
                }
                break;
            case 'v':
                k = sscanf(optarg, "%" PRIi64 "", &minimumSequenceLengthForBlast);
                assert(k == 1);
                break;
            case 'w':
                k = sscanf(optarg, "%f", &maximumAdjacencyComponentSizeRatio);
                assert(k == 1);
                break;
            case 'x':
                constraintsFile = stString_copy(optarg);
                break;
            case 'y':
                k = sscanf(optarg, "%" PRIi64 "", &minLengthForChromosome);
                assert(k == 1);
                break;
            case 'z':
                k = sscanf(optarg, "%f", &proportionOfUnalignedBasesForNewChromosome);
                assert(k == 1);
                break;
            case 'A':
                k = sscanf(optarg, "%" PRIi64 "", &maximumMedianSequenceLengthBetweenLinkedEnds);
                assert(k == 1);
                break;
            case 'B':
                realign = 1;
                break;
            case 'C':
                realignArguments = stString_copy(optarg);
                break;
            case 'D':
                if (strcmp(optarg, "1") == 0) {
                    removeRecoverableChains = true;
                    recoverableChainsFilter = NULL;
                } else if (strcmp(optarg, "unequalNumberOfIngroupCopies") == 0) {
                    removeRecoverableChains = true;
                    recoverableChainsFilter = stCaf_chainHasUnequalNumberOfIngroupCopies;
                } else if (strcmp(optarg, "unequalNumberOfIngroupCopiesOrNoOutgroup") == 0) {
                    removeRecoverableChains = true;
                    recoverableChainsFilter = stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup;
                } else if (strcmp(optarg, "0") == 0) {
                    removeRecoverableChains = false;
                } else {
                    st_errAbort("Could not parse removeRecoverableChains argument");
                }
                break;
            case 'E':
                k = sscanf(optarg, "%" PRIi64, &minimumNumberOfSpecies);
                if (k != 1) {
                    st_errAbort("Error parsing the minimumNumberOfSpecies argument");
                }
                break;
            case 'F':
                k = sscanf(optarg, "%" PRIi64, &maxRecoverableChainsIterations);
                if (k != 1) {
                    st_errAbort("Error parsing the maxRecoverableChainsIterations argument");
                }
                break;
            case 'G':
                k = sscanf(optarg, "%" PRIi64, &maxRecoverableChainLength);
                if (k != 1) {
                    st_errAbort("Error parsing the maxRecoverableChainLength argument");
                }
                break;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(cactusDiskDatabaseString != NULL);
    assert(minimumTreeCoverage >= 0.0);
    assert(minimumTreeCoverage <= 1.0);
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
    assert(minimumOutgroupDegree >= 0);
    assert(minimumIngroupDegree >= 0);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower disk name : %s\n", cactusDiskDatabaseString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Sort the constraints
    ///////////////////////////////////////////////////////////////////////////

    stPinchIterator *pinchIteratorForConstraints = NULL;
    if (constraintsFile != NULL) {
        pinchIteratorForConstraints = stPinchIterator_constructFromFile(constraintsFile);
        st_logInfo("Created an iterator for the alignment constaints from file: %s\n", constraintsFile);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Do the alignment
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);

    stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    if (alignmentsFile == NULL) {
        cactusDisk_preCacheStrings(cactusDisk, flowers);
    }
    char *tempFile1 = NULL;
    for (int64_t i = 0; i < stList_length(flowers); i++) {
        flower = stList_get(flowers, i);
        if (!flower_builtBlocks(flower)) { // Do nothing if the flower already has defined blocks
            st_logDebug("Processing flower: %lli\n", flower_getName(flower));

            stCaf_setFlowerForAlignmentFiltering(flower);

            //Set up the graph and add the initial alignments
            stPinchThreadSet *threadSet = stCaf_setup(flower);

            bool sortAlignments = 0;
            if (singleCopyIngroup) {
                sortAlignments = 1;
                filterFn = stCaf_filterByRepeatSpecies;
            }
            else if (singleCopyOutgroup || relaxedSingleCopyOutgroup) {
                //Here is where we get the set of outgroup threads.
                outgroupThreads = stSet_construct();
                stPinchThreadSetIt it = stPinchThreadSet_getIt(threadSet);
                stPinchThread *pinchThread;
                while ((pinchThread = stPinchThreadSetIt_getNext(&it)) != NULL) {
                    if (event_isOutgroup(stCaf_getEvent(stPinchThread_getFirst(pinchThread), flower))) {
                        stSet_insert(outgroupThreads, pinchThread);
                    }
                }

                if (stSet_size(outgroupThreads) == 0) {
                    filterFn = NULL;
                    sortAlignments = 0;
                } else if (singleCopyOutgroup) {
                    filterFn = stCaf_filterByOutgroup;
                    sortAlignments = 1;
                } else if (relaxedSingleCopyOutgroup) {
                    filterFn = stCaf_relaxedFilterByOutgroup;
                    sortAlignments = 1;
                }
            }

            //Setup the alignments
            stPinchIterator *pinchIterator;
            stList *alignmentsList = NULL;
            if (alignmentsFile != NULL) {
                assert(i == 0);
                assert(stList_length(flowers) == 1);
                if (sortAlignments) {
                    tempFile1 = getTempFile();
                    stCaf_sortCigarsFileByScoreInDescendingOrder(alignmentsFile, tempFile1);
                    pinchIterator = stPinchIterator_constructFromFile(tempFile1);
                } else {
                    pinchIterator = stPinchIterator_constructFromFile(alignmentsFile);
                }
            } else {
                if (tempFile1 == NULL) {
                    tempFile1 = getTempFile();
                }
                alignmentsList = stCaf_selfAlignFlower(flower, minimumSequenceLengthForBlast, lastzArguments, realign, realignArguments, tempFile1);
                if (sortAlignments) {
                    stCaf_sortCigarsByScoreInDescendingOrder(alignmentsList);
                }
                st_logDebug("Ran lastz and have %" PRIi64 " alignments\n", stList_length(alignmentsList));
                pinchIterator = stPinchIterator_constructFromList(alignmentsList);
            }

            for (int64_t annealingRound = 0; annealingRound < annealingRoundsLength; annealingRound++) {
                int64_t minimumChainLength = annealingRounds[annealingRound];
                int64_t alignmentTrim = annealingRound < alignmentTrimLength ? alignmentTrims[annealingRound] : 0;
                st_logDebug("Starting annealing round with a minimum chain length of %" PRIi64 " and an alignment trim of %" PRIi64 "\n", minimumChainLength, alignmentTrim);
                stPinchIterator_setTrim(pinchIterator, alignmentTrim);

                //Add back in the constraints
                if (pinchIteratorForConstraints != NULL) {
                    stCaf_anneal(threadSet, pinchIteratorForConstraints, filterFn);
                }

                //Do the annealing
                if (annealingRound == 0) {
                    stCaf_anneal(threadSet, pinchIterator, filterFn);
                } else {
                    stCaf_annealBetweenAdjacencyComponents(threadSet, pinchIterator, filterFn);
                }

                //Do the melting rounds
                for (int64_t meltingRound = 0; meltingRound < meltingRoundsLength; meltingRound++) {
                    int64_t minimumChainLengthForMeltingRound = meltingRounds[meltingRound];
                    st_logDebug("Starting melting round with a minimum chain length of %" PRIi64 " \n", minimumChainLengthForMeltingRound);
                    if (minimumChainLengthForMeltingRound >= minimumChainLength) {
                        break;
                    }
                    stCaf_melt(flower, threadSet, NULL, 0, minimumChainLengthForMeltingRound, 0, INT64_MAX);
                } st_logDebug("Last melting round of cycle with a minimum chain length of %" PRIi64 " \n", minimumChainLength);
                stCaf_melt(flower, threadSet, NULL, 0, minimumChainLength, breakChainsAtReverseTandems, maximumMedianSequenceLengthBetweenLinkedEnds);
                //This does the filtering of blocks that do not have the required species/tree-coverage/degree.
                stCaf_melt(flower, threadSet, blockFilterFn, blockTrim, 0, 0, INT64_MAX);
            }

            if (removeRecoverableChains) {
                stCaf_meltRecoverableChains(flower, threadSet, breakChainsAtReverseTandems, maximumMedianSequenceLengthBetweenLinkedEnds, recoverableChainsFilter, maxRecoverableChainsIterations, maxRecoverableChainLength);
            }

            //Sort out case when we allow blocks of degree 1
            if (minimumDegree < 2) {
                st_logDebug("Creating degree 1 blocks\n");
                stCaf_makeDegreeOneBlocks(threadSet);
                stCaf_melt(flower, threadSet, blockFilterFn, blockTrim, 0, 0, INT64_MAX);
            } else if (maximumAdjacencyComponentSizeRatio < INT64_MAX) { //Deal with giant components
                st_logDebug("Breaking up components greedily\n");
                stCaf_breakupComponentsGreedily(threadSet, maximumAdjacencyComponentSizeRatio);
            }

            //Finish up
            stCaf_finish(flower, threadSet, chainLengthForBigFlower, longChain, minLengthForChromosome,
                    proportionOfUnalignedBasesForNewChromosome); //Flower is then destroyed at this point.
            st_logInfo("Ran the cactus core script\n");

            //Cleanup
            stPinchThreadSet_destruct(threadSet);
            stPinchIterator_destruct(pinchIterator);

            if (singleCopyOutgroup && !singleCopyIngroup) {
                stSet_destruct(outgroupThreads);
            }

            if (alignmentsList != NULL) {
                stList_destruct(alignmentsList);
            }
        } else {
            st_logInfo("We've already built blocks / alignments for this flower\n");
        }
    }
    stList_destruct(flowers);
    if (tempFile1 != NULL) {
        st_system("rm %s", tempFile1);
    }

    if (constraintsFile != NULL) {
        stPinchIterator_destruct(pinchIteratorForConstraints);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Write the flower to disk.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk and %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    return 0; //Exit without clean up is quicker, enable cleanup when doing memory leak detection.

    //Destruct stuff
    startTime = time(NULL);
    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    free(cactusDiskDatabaseString);
    if (alignmentsFile != NULL) {
        free(alignmentsFile);
    }
    if (logLevelString != NULL) {
        free(logLevelString);
    }
    if (lastzArguments != NULL) {
        free(lastzArguments);
    }
    st_logInfo("Cleaned stuff up and am finished in: %" PRIi64 " seconds\n", time(NULL) - startTime);

    //while(1);
    return 0;
}
