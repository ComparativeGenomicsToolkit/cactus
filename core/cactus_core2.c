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
    fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr, "-e --writeDebugFiles : Write the debug files\n");
    fprintf(stderr, "-h --help : Print this help screen\n");

    fprintf(stderr, "-i --annealingRounds (int >= 1) : The number of rounds of alignment, undoing of over-aligned edges and recursion into adjacency connected components (groups)\n");
    fprintf(stderr, "-j --alignRepeatsAtRound (int  [0, alignUndoLoops) ) : Allow bases marked as repeats to be aligned at loop (else alignments to these bases to be excluded)\n");

    fprintf(stderr, "-k --trim (float >= 0) : The length of bases to remove from the end of each alignment\n");
    fprintf(stderr, "-l --trimChange : (float) Trim reduction, the amount to reduce the trim after each align/undo loop (to a minimum of zero)\n");

    fprintf(stderr, "-m --minimumTreeCoverage : (float [0.0, 1.0]) Minimum tree coverage proportion of an block to be included in the graph\n");

    fprintf(stderr, "-n --minimumBlockLength : (int >= 0) The minimum length of a block required to be included in the problem\n");
    fprintf(stderr, "-o --minimumBlockLengthChange : (float) The minimum-block-length increase after each align/undo loop\n");

    fprintf(stderr, "-p --minimumChainLength : (int >= 0) The minimum chain length required to be included in the problem\n");
    fprintf(stderr, "-q --minimumChainLengthChange : (float) The minimum-chain-length increase after each align/undo loop\n");

    fprintf(stderr, "-r --deannealingRounds: (float >= 1.0) The amount to increase the minimum chain length after each cactus tree undo loop\n");

    fprintf(stderr, "-s --adjacencyComponentOverlap : (int >= 0) When adding alignmnts back into the graph, the distance between the adjacency components to consider\n");
}

int main(int argc, char *argv[]) {
    /*
     * Script for adding alignments to cactus tree.
     */
    int32_t startTime;
    stKVDatabaseConf *kvDatabaseConf;
    CactusDisk *cactusDisk;
    Flower *flower;
    int key;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * alignmentsFile = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = NULL;
    CactusCoreInputParameters *cCIP = constructCactusCoreInputParameters();

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option
                long_options[] = {
                        { "logLevel", required_argument, 0, 'a' },
                        { "alignments", required_argument, 0, 'b' },
                        { "cactusDisk", required_argument, 0, 'c' },
                        { "flowerName", required_argument, 0, 'd' },
                        { "writeDebugFiles", no_argument, 0, 'e' },
                        { "help", no_argument, 0, 'h' },
                        { "annealingRounds", required_argument, 0, 'i' },
                        { "alignRepeatsAtRound", required_argument, 0, 'j' },
                        { "trim", required_argument, 0, 'k' },
                        { "trimChange", required_argument, 0, 'l', },
                        { "minimumTreeCoverage", required_argument, 0, 'm' },
                        { "minimumBlockLength", required_argument, 0, 'n' },
                        { "minimumBlockLengthChange", required_argument, 0, 'o' },
                        { "minimumChainLength", required_argument, 0, 'p' },
                        { "minimumChainLengthChange", required_argument, 0, 'q', },
                        { "deannealingRounds", required_argument, 0, 'r', },
                        { "adjacencyComponentOverlap", required_argument, 0, 's', },
                        { 0, 0, 0, 0 } };

        int option_index = 0;

        key = getopt_long(argc, argv,
                "a:b:c:d:ehi:j:k:l:m:n:o:p:q:r:s:", long_options,
                &option_index);

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
            assert(sscanf(optarg, "%i", &cCIP->annealingRounds) == 1);
            break;
        case 'j':
            assert(sscanf(optarg, "%i", &cCIP->alignRepeatsAtRound) == 1);
            break;
        case 'k':
            assert(sscanf(optarg, "%i", &cCIP->trim) == 1);
            break;
        case 'l':
            assert(sscanf(optarg, "%f", &cCIP->trimChange) == 1);
            break;
        case 'm':
            assert(sscanf(optarg, "%f", &cCIP->minimumTreeCoverage) == 1);
            break;
        case 'n':
            assert(sscanf(optarg, "%i", &cCIP->minimumBlockLength) == 1);
            break;
        case 'o':
            assert(sscanf(optarg, "%f", &cCIP->minimumBlockLengthChange) == 1);
            break;
        case 'p':
            assert(sscanf(optarg, "%i", &cCIP->minimumChainLength) == 1);
            break;
        case 'q':
            assert(sscanf(optarg, "%f", &cCIP->minimumChainLengthChange) == 1);
            break;
        case 'r':
             assert(sscanf(optarg, "%i", &cCIP->deannealingRounds) == 1);
             break;
        case 's':
             assert(sscanf(optarg, "%i", &cCIP->adjacencyComponentOverlap) == 1);
             break;

        default:
            usage();
            return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(logLevelString == NULL || strcmp(logLevelString, "CRITICAL") == 0 || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
    assert(alignmentsFile != NULL);
    assert(cactusDiskDatabaseString != NULL);
    assert(flowerName != NULL);
    assert(cCIP->minimumTreeCoverage >= 0.0);
    assert(cCIP->minimumTreeCoverage <= 1.0);
    assert(cCIP->minimumBlockLength >= 0);
    assert(cCIP->minimumChainLength >= 0);
    assert(cCIP->trim >= 0);
    assert(cCIP->annealingRounds >= 0);
    assert(cCIP->deannealingRounds >= 1);
    assert(cCIP->alignRepeatsAtRound >= 0);
    assert(cCIP->adjacencyComponentOverlap >= 0);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if (logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if (logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Pairwise alignments file : %s\n", alignmentsFile);
    st_logInfo("Flower disk name : %s\n", cactusDiskDatabaseString);
    st_logInfo("Flower name : %s\n", flowerName);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(flowerName));
    st_logInfo("Parsed the flower to be refined\n");

    if (!flower_builtBlocks(flower)) { // Do nothing if the flower already has defined blocks
        startTime = time(NULL);

        ///////////////////////////////////////////////////////////////////////////
        // Call the core program.
        ///////////////////////////////////////////////////////////////////////////

        startAlignmentStack_fileString = alignmentsFile;
        exitOnFailure(cactusCorePipeline(flower, cCIP, getNextAlignment,
                startAlignmentStack, 0),
                "Failed to run the cactus core pipeline\n");
        fclose(getNextAlignment_FileHandle);

        ///////////////////////////////////////////////////////////////////////////
        // (9) Write the flower to disk.
        ///////////////////////////////////////////////////////////////////////////

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

    st_logInfo("Cleaned stuff up and am finished in: %i seconds\n", time(NULL)
            - startTime);

    return 0;
}
