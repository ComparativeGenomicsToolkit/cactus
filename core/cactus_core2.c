/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "pinchGraph.h"
#include "cactusGraph.h"
#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "cactus.h"
#include "pairwiseAlignment.h"
#include "cactusFlowerFunctions.h"
#include "cactus_core.h"
#include "lastzAlignments.h"
#include "sonLib.h"
#include "pairwiseAlignmentIterator.h"

void usage() {
    fprintf(stderr, "cactus_core, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --alignments : The input alignments file\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-d --lastzArguments : Lastz arguments\n");
    fprintf(stderr, "-e --writeDebugFiles : Write the debug files\n");
    fprintf(stderr, "-h --help : Print this help screen\n");

    fprintf(
            stderr,
            "-i --annealingRounds (array of ints, each greater than or equal to 1) : The rounds of annealing\n");
    fprintf(
            stderr,
            "-o --deannealingRounds (array of ints, each greater than or equal to 1 and each greater than the last) : The rounds of deannealing\n");

    fprintf(
            stderr,
            "-j --alignRepeatsAtRound (int  [0, alignUndoLoops) ) : Allow bases marked as repeats to be aligned at loop (else alignments to these bases to be excluded)\n");

    fprintf(
            stderr,
            "-k --trim (array of integers, each greater or equal to zero) : An array giving the trim for each annealing round. If the array is shorter than the annealing rounds then a trim value of 0 is assumed for annealing rounds greater than the length of the trim array\n");

    fprintf(
            stderr,
            "-m --minimumTreeCoverage : (float [0.0, 1.0]) Minimum tree coverage proportion of a block to be included in the graph\n");

    fprintf(
            stderr,
            "-n --blockTrim : (int >= 0) The number of bases to trim from the ends of each block in a chain before accepting, this filtering is done after choosing the length of chains\n");

    fprintf(
            stderr,
            "-p --minimumDegree : (int >= 0) Minimum number of sequences in a block to be included in the output graph\n");

    fprintf(
            stderr,
            "-q --requiredIngroupFraction : Fraction of ingroup events required in a block.\n");

    fprintf(
            stderr,
            "-r --requiredOutgroupFraction : Fraction of outgroup events required in a block.\n");

    fprintf(stderr,
            "-u --requiredAllFraction : Fraction of all events required in a block.\n");

    fprintf(
            stderr,
            "-s --singleCopyIngroup : Require that in-group sequences have only single coverage\n");

    fprintf(
            stderr,
            "-t --singleCopyOutgroup : Require that out-group sequences have only single coverage\n");

    fprintf(
            stderr,
            "-v --minimumSequenceLengthForBlast : The minimum length of a sequence to include when blasting\n");

    fprintf(
            stderr,
            "-w --maxAdjacencyComponentSizeRatio : The components equal or less than log(n) * of this size will be allowed in the cactus. Used to fight giant components.\n");

    fprintf(stderr,
            "-x --constraints : A file of alignments that will enforced upon the cactus\n");
}

stSortedSet *getStringSet(const char *string) {
    stList *list = stString_split(string);
    stSortedSet *strings = stSortedSet_construct3(
            (int(*)(const void *, const void *)) strcmp, free);
    for (int32_t i = 0; i < stList_length(list); i++) {
        stSortedSet_insert(strings, stString_copy(stList_get(list, i)));
    }
    stList_destruct(list);
    return strings;
}

int32_t *getInts(const char *string, int32_t *arrayLength) {
    int32_t *iA = st_malloc(sizeof(int32_t) * strlen(string));
    char *cA = stString_copy(string);
    char *cA2 = cA;
    char *cA3;
    *arrayLength = 0;
    while ((cA3 = stString_getNextWord(&cA)) != NULL) {
        int32_t i = sscanf(cA3, "%i", &iA[(*arrayLength)++]);
        assert(i == 1);
        free(cA3);
    }
    free(cA2);
    return iA;
}

int main(int argc, char *argv[]) {
    /*
     * Script for adding alignments to cactus tree.
     */
    int32_t startTime;
    stKVDatabaseConf *kvDatabaseConf;
    CactusDisk *cactusDisk;
    int key, k;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * alignmentsFile = NULL;
    char * constraintsFile = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * lastzArguments = NULL;
    int32_t minimumSequenceLengthForBlast = 1;
    CactusCoreInputParameters *cCIP = constructCactusCoreInputParameters();

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "alignments", required_argument,
                0, 'b' }, { "cactusDisk", required_argument, 0, 'c' }, {
                "lastzArguments", required_argument, 0, 'd' }, {
                "writeDebugFiles", no_argument, 0, 'e' }, { "help",
                no_argument, 0, 'h' }, { "annealingRounds", required_argument,
                0, 'i' }, { "alignRepeatsAtRound", required_argument, 0, 'j' },
                { "trim", required_argument, 0, 'k' }, { "trimChange",
                        required_argument, 0, 'l', }, { "minimumTreeCoverage",
                        required_argument, 0, 'm' }, { "blockTrim",
                        required_argument, 0, 'n' }, { "deannealingRounds",
                        required_argument, 0, 'o' }, { "minimumDegree",
                        required_argument, 0, 'p' }, {
                        "requiredIngroupFraction", required_argument, 0, 'q' },
                { "requiredOutgroupFraction", required_argument, 0, 'r' }, {
                        "requiredAllFraction", required_argument, 0, 'u' }, {
                        "singleCopyIngroup", no_argument, 0, 's' }, {
                        "singleCopyOutgroup", no_argument, 0, 't' }, {
                        "minimumSequenceLengthForBlast", required_argument, 0,
                        'v' }, { "maxAdjacencyComponentSizeRatio",
                        required_argument, 0, 'w' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        key = getopt_long(argc, argv, "a:b:c:ehi:j:k:m:n:o:p:q:r:stu:v:w:",
                long_options, &option_index);

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
            case 'e':
                cCIP->writeDebugFiles = !cCIP->writeDebugFiles;
                break;
            case 'h':
                usage();
                return 0;
            case 'i':
                free(cCIP->annealingRounds);
                cCIP->annealingRounds = getInts(optarg,
                        &cCIP->annealingRoundsLength);
                break;
            case 'o':
                free(cCIP->deannealingRounds);
                cCIP->deannealingRounds = getInts(optarg,
                        &cCIP->deannealingRoundsLength);
                break;
            case 'j':
                k = sscanf(optarg, "%i", &cCIP->alignRepeatsAtRound);
                assert(k == 1);
                break;
            case 'k':
                free(cCIP->trim);
                cCIP->trim = getInts(optarg, &cCIP->trimLength);
                break;
            case 'm':
                k = sscanf(optarg, "%f", &cCIP->minimumTreeCoverage);
                assert(k == 1);
                break;
            case 'n':
                k = sscanf(optarg, "%i", &cCIP->blockTrim);
                assert(k == 1);
                break;
            case 'p':
                k = sscanf(optarg, "%i", &cCIP->minimumDegree);
                assert(k == 1);
                break;
            case 'q':
                k = sscanf(optarg, "%f", &cCIP->requiredIngroupFraction);
                assert(k == 1);
                break;
            case 'r':
                k = sscanf(optarg, "%f", &cCIP->requiredOutgroupFraction);
                assert(k == 1);
                break;
            case 'u':
                k = sscanf(optarg, "%f", &cCIP->requiredAllFraction);
                assert(k == 1);
                break;
            case 's':
                cCIP->singleCopyIngroup = 1;
                break;
            case 't':
                cCIP->singleCopyOutgroup = 1;
                break;
            case 'v':
                k = sscanf(optarg, "%i", &minimumSequenceLengthForBlast);
                assert(k == 1);
                break;
            case 'w':
                k = sscanf(optarg, "%f", &cCIP->maxAdjacencyComponentSizeRatio);
                assert(k == 1);
                break;
            case 'x':
                constraintsFile = stString_copy(optarg);
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
    assert(cCIP->minimumTreeCoverage >= 0.0);
    assert(cCIP->minimumTreeCoverage <= 1.0);
    assert(cCIP->blockTrim >= 0);
    assert(cCIP->annealingRoundsLength >= 0);
    for (int32_t i = 0; i < cCIP->annealingRoundsLength; i++) {
        assert(cCIP->annealingRounds[i] >= 0);
    }
    assert(cCIP->deannealingRoundsLength >= 0);
    for (int32_t i = 1; i < cCIP->deannealingRoundsLength; i++) {
        assert(cCIP->deannealingRounds[i - 1] < cCIP->deannealingRounds[i]);
        assert(cCIP->deannealingRounds[i - 1] >= 1);
    }
    assert(cCIP->trimLength >= 0);
    for (int32_t i = 0; i < cCIP->trimLength; i++) {
        assert(cCIP->trim[i] >= 0);
    }
    assert(cCIP->alignRepeatsAtRound >= 0);
    assert(cCIP->requiredAllFraction >= 0);
    assert(cCIP->requiredOutgroupFraction >= 0);
    assert(cCIP->requiredIngroupFraction >= 0);

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

    kvDatabaseConf = stKVDatabaseConf_constructFromString(
            cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Sort the constraints
    ///////////////////////////////////////////////////////////////////////////

    PairwiseAlignmentIterator *pairwiseAlignmentIteratorForConstraints = NULL;
    if (pairwiseAlignmentIteratorForConstraints != NULL) {
        pairwiseAlignmentIteratorForConstraints = pairwiseAlignmentIterator_constructFromFile(
                constraintsFile);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Do the alignment
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);

    stList *flowers = parseFlowersFromStdin(cactusDisk);
    char *tempFile1 = NULL;
    for (int32_t i = 0; i < stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        if (!flower_builtBlocks(flower)) { // Do nothing if the flower already has defined blocks
            st_logDebug("Processing flower: %lli\n", flower_getName(flower));
            PairwiseAlignmentIterator *pairwiseAlignmentIterator;
            stList *alignmentsList = NULL;
            if (alignmentsFile != NULL) {
                assert(i == 0);
                assert(stList_length(flowers) == 1);
                pairwiseAlignmentIterator
                        = pairwiseAlignmentIterator_constructFromFile(
                                alignmentsFile);
            } else {
                if (tempFile1 == NULL) {
                    tempFile1 = getTempFile();
                }
                alignmentsList = selfAlignFlower(flower,
                        minimumSequenceLengthForBlast, lastzArguments,
                        tempFile1);
                st_logDebug("Ran lastz and have %i alignments\n",
                        stList_length(alignmentsList));
                pairwiseAlignmentIterator
                        = pairwiseAlignmentIterator_constructFromList(
                                alignmentsList);
            }
            exitOnFailure(
                    cactusCorePipeline(flower, cCIP, pairwiseAlignmentIterator, pairwiseAlignmentIteratorForConstraints),
                    "Failed to run the cactus core pipeline from a file of alignments\n");
            pairwiseAlignmentIterator_destruct(pairwiseAlignmentIterator);
            assert(!flower_isParentLoaded(flower));
            if(alignmentsList != NULL) {
                stList_destruct(alignmentsList);
            }
        } else {
            st_logInfo(
                    "We've already built blocks / alignments for this flower\n");
        }
    }
    stList_destruct(flowers);
    if (tempFile1 != NULL) {
        st_system("rm %s", tempFile1);
    }

    if (constraintsFile != NULL) {
        pairwiseAlignmentIterator_destruct(pairwiseAlignmentIteratorForConstraints);
    }

    ///////////////////////////////////////////////////////////////////////////
    // (9) Write the flower to disk.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk\n");

    ///////////////////////////////////////////////////////////////////////////
    //(10) Clean up.
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
    destructCactusCoreInputParameters(cCIP);

    st_logInfo("Cleaned stuff up and am finished in: %i seconds\n",
            time(NULL) - startTime);

    //while(1);

    return 0;
}
