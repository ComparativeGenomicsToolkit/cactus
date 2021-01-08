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

static void usage() {
    fprintf(stderr, "cactus_caf, version 0.2\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-p --params : [Required] The cactus config file\n");
    fprintf(stderr, "-d --cactusDisk : [Required] The database config string\n");
    fprintf(stderr, "-a --alignments : [Required] The alignments file\n");
    fprintf(stderr, "-S --secondaryAlignments : The secondary alignments file\n");
    fprintf(stderr, "-c --constraintAlignments : The constraint alignments file\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    /*
     * Script for adding alignments to cactus tree.
     */

    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * alignmentsFile = NULL;
    char * secondaryAlignmentsFile = NULL;
    char * constraintAlignmentsFile = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * paramsFile = NULL;

    // To be deleted:
    //const char *referenceEventHeader = NULL;
    //char * constraintsFile = NULL;
    //char * lastzArguments = "";
    //const char *debugFileName = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "params", required_argument, 0, 'p' },
                                                { "cactusDisk", required_argument, 0, 'd' },
                                                { "alignments", required_argument, 0, 'a' },
                                                { "secondaryAlignments", required_argument, 0, 'S' },
                                                { "constraintAlignments", required_argument, 0, 'c' },
                                                { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int64_t key = getopt_long(argc, argv, "l:p:d:a:S:c:h", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'p':
                paramsFile = optarg;
                break;
            case 'd':
                cactusDiskDatabaseString = optarg;
                break;
            case 'a':
                alignmentsFile = optarg;
                break;
            case 'S':
                secondaryAlignmentsFile = optarg;
                break;
            case 'c':
                constraintAlignmentsFile = optarg;
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_logInfo("Params file: %s\n", paramsFile);
    st_logInfo("Cactus disk database string : %s\n", cactusDiskDatabaseString);
    st_logInfo("Alignments file: %s\n", alignmentsFile);
    st_logInfo("Secondary alignments file: %s\n", secondaryAlignmentsFile);
    st_logInfo("Constraint alignments file: %s\n", constraintAlignmentsFile);

    //////////////////////////////////////////////
    //Parse stuff
    //////////////////////////////////////////////

    // Load the params file
    CactusParams *params = cactusParams_load(paramsFile);
    st_logInfo("Loaded the parameters files\n");

    // Load the cactus disk
    stKVDatabaseConf *kvDatabaseConf = kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, true, true);
    st_logInfo("Set up the cactus disk\n");

    //////////////////////////////////////////////
    //Run CAF
    //////////////////////////////////////////////

    Flower *flower = NULL; // Get the root flower
    caf(flower, params, alignmentsFile, secondaryAlignmentsFile, constraintAlignmentsFile);
    st_logInfo("Ran cactus caf\n");

    ///////////////////////////////////////////////////////////////////////////
    // Write the flower to disk.
    ///////////////////////////////////////////////////////////////////////////

    st_logDebug("Writing the flowers to disk\n");
    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk and %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
}
