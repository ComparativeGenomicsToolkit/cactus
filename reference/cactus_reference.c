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

#include "cactus.h"
#include "sonLib.h"
#include "cactusReference.h"

void usage() {
    fprintf(stderr, "cactus_reference [flower names], version 0.1\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-p --params : [Required] The cactus config file\n");
    fprintf(stderr, "-d --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    /*
     * Script for adding a reference genome to a flower.
     */

    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * paramsFile = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    if(argc <= 1) {
        usage();
        return 0;
    }

    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'l' },
                { "params", required_argument, 0, 'p' },
                { "cactusDisk", required_argument, 0, 'd' },
                { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "l:d:h", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'p':
                paramsFile = optarg;
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
    //Parse stuff
    //////////////////////////////////////////////

    // Load the params file
    CactusParams *params = cactusParams_load(paramsFile);
    st_logInfo("Loaded the parameters files, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    // Load the cactus disk
    assert(cactusDiskDatabaseString != NULL);
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, false, true); //We precache the sequences
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Run the cactus reference algorithm
    ///////////////////////////////////////////////////////////////////////////

    stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    cactus_make_reference(flowers, cactusDisk, params);

    ///////////////////////////////////////////////////////////////////////////
    // Write the flower(s) back to disk.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk\n");

    ///////////////////////////////////////////////////////////////////////////
    //Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);

    return 0; //Exit without clean up is quicker, enable cleanup when doing memory leak detection.

    stKVDatabaseConf_destruct(kvDatabaseConf);
    free(cactusDiskDatabaseString);
    if (logLevelString != NULL) {
        free(logLevelString);
    }

    st_logInfo("Cleaned stuff up and am finished\n");

    //while(1);

    return 0;
}
