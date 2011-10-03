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
#include "sonLib.h"

char *startAlignmentStack_fileString;
static FILE *getNextAlignment_FileHandle = NULL;

static struct PairwiseAlignment *getNextAlignment() {
    return cigarRead(getNextAlignment_FileHandle);
}

static void startAlignmentStack() {
    if (getNextAlignment_FileHandle != NULL) {
        fclose(getNextAlignment_FileHandle);
    }
    getNextAlignment_FileHandle = fopen(startAlignmentStack_fileString, "r");
}

void usage() {
    fprintf(stderr, "cactus_core, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --alignments : The input alignments file\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-d --flowerName : The name of the flower (the key in the database)\n");
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
            "-m --minimumTreeCoverage : (float [0.0, 1.0]) Minimum tree coverage proportion of an block to be included in the graph\n");

    fprintf(
            stderr,
            "-n --blockTrim : (int >= 0) The number of bases to trim from the ends of each block in a chain before accepting, this filtering is done after choosing the length of chains\n");

    fprintf(
            stderr,
            "-p --minimumDegree : (int >= 0) Minimum number of sequences in a block to be included in the output graph\n");

    fprintf(
            stderr,
            "-q --requiredSpecies : (array of event names) At least one of these events must be in a block for it to be included in the final graph\n");

    fprintf(
            stderr,
            "-t --singleCopySpecies : (array of event names) If any of the following events contribute multiple instances to a block it will be broken up\n");

    fprintf(
            stderr,
            "-u --minimumChainLength : (int >= 0) The minimum acceptable length of a chain\n");

    fprintf(stderr,
            "-v --maximumGroupSize : (int >= 0) The maximum acceptable size for a group\n");

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
    Flower *flower;
    int key, k, tempValue; //Bad hack

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * alignmentsFile = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = NULL;
    CactusCoreInputParameters *cCIP = constructCactusCoreInputParameters();
    cCIP->minimumChainLength = -1;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "alignments", required_argument,
                0, 'b' }, { "cactusDisk", required_argument, 0, 'c' }, {
                "flowerName", required_argument, 0, 'd' }, { "writeDebugFiles",
                no_argument, 0, 'e' }, { "help", no_argument, 0, 'h' }, {
                "annealingRounds", required_argument, 0, 'i' }, {
                "alignRepeatsAtRound", required_argument, 0, 'j' }, { "trim",
                required_argument, 0, 'k' }, { "trimChange", required_argument,
                0, 'l', },
                { "minimumTreeCoverage", required_argument, 0, 'm' }, {
                        "blockTrim", required_argument, 0, 'n' }, {
                        "deannealingRounds", required_argument, 0, 'o' }, {
                        "minimumDegree", required_argument, 0, 'p' }, {
                        "requiredSpecies", required_argument, 0, 'q' }, {
                        "singleCopySpecies", required_argument, 0, 't' }, {
                        "minimumChainLength", required_argument, 0, 'u' }, {
                        "maximumGroupSize", required_argument, 0, 'v' }, { 0,
                        0, 0, 0 } };

        int option_index = 0;

        key = getopt_long(argc, argv, "a:b:c:d:ehi:j:k:m:n:o:p:q:t:",
                long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                alignmentsFile = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'd':
                flowerName = stString_copy(optarg);
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
                parseRequiredSpeciesTree(optarg, cCIP);
                break;
            case 't':
                cCIP->singleCopySpecies = getStringSet(optarg);
                break;
            case 'u':
                k = sscanf(optarg, "%i", &cCIP->minimumChainLength);
                assert(k == 1);
                assert(cCIP->minimumChainLength >= 0);
                break;
            case 'v':
                k = sscanf(optarg, "%i", &tempValue);
                assert(k == 1);
                cCIP->maximumAdjacencyComponentSize = tempValue;
                assert(cCIP->maximumAdjacencyComponentSize >= 0);
                break;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(alignmentsFile != NULL);
    assert(cactusDiskDatabaseString != NULL);
    assert(flowerName != NULL);
    assert(cCIP->minimumTreeCoverage >= 0.0);
    assert(cCIP->minimumTreeCoverage <= 1.0);
    assert(cCIP->blockTrim >= 0);
    assert(cCIP->annealingRoundsLength >= 1);
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
    assert(cCIP->maximumAdjacencyComponentSize >= 0);

    if (cCIP->minimumChainLength == -1) {
        cCIP->minimumChainLength
                = cCIP->annealingRounds[cCIP->annealingRoundsLength - 1];
        st_logDebug(
                "Setting the maximum chain length to %i from the last value in the annealing rounds array\n",
                cCIP->minimumChainLength);
    }
    assert(cCIP->minimumChainLength >= 0);
    assert(
            cCIP->minimumChainLength
                    <= cCIP->annealingRounds[cCIP->annealingRoundsLength - 1]);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Pairwise alignments file : %s\n", alignmentsFile);
    st_logInfo("Flower disk name : %s\n", cactusDiskDatabaseString);
    st_logInfo("Flower name : %s\n", flowerName);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    kvDatabaseConf = stKVDatabaseConf_constructFromString(
            cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    flower = cactusDisk_getFlower(cactusDisk,
            cactusMisc_stringToName(flowerName));
    st_logInfo("Parsed the flower to be refined\n");

    if (!flower_builtBlocks(flower)) { // Do nothing if the flower already has defined blocks
        startTime = time(NULL);

        ///////////////////////////////////////////////////////////////////////////
        // Call the core program.
        ///////////////////////////////////////////////////////////////////////////

        startAlignmentStack_fileString = alignmentsFile;
        exitOnFailure(
                cactusCorePipeline(flower, cCIP, getNextAlignment,
                        startAlignmentStack),
                "Failed to run the cactus core pipeline\n");
        fclose(getNextAlignment_FileHandle);

        ///////////////////////////////////////////////////////////////////////////
        // (9) Write the flower to disk.
        ///////////////////////////////////////////////////////////////////////////

        assert(!flower_isParentLoaded(flower));
        flower_unloadParent(flower); //The parent should not have changed.
        cactusDisk_write(cactusDisk);
        st_logInfo("Updated the flower on disk\n");
    } else {
        st_logInfo("We've already built blocks / alignments for this flower\n");
    }

    ///////////////////////////////////////////////////////////////////////////
    //(10) Clean up.
    ///////////////////////////////////////////////////////////////////////////

    //Destruct stuff
    startTime = time(NULL);
    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    free(cactusDiskDatabaseString);
    free(alignmentsFile);
    free(flowerName);
    if (logLevelString != NULL) {
        free(logLevelString);
    }
    destructCactusCoreInputParameters(cCIP);

    st_logInfo("Cleaned stuff up and am finished in: %i seconds\n",
            time(NULL) - startTime);

    //while(1);

    return 0;
}
