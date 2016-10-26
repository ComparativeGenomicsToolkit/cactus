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
#include "normal.h"

/*
 * This script is run bottom up on a cactus tree and ensures the tree is 'normalised' as we expect.
 */

void usage() {
    fprintf(stderr, "cactus_normalisation [flower names], version 0.1\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(
            stderr,
            "-d --maxNumberOfChains : The maximum number of individual chains to promote into a flower.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    /*
     * Script for adding a reference genome to a flower.
     */

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    int64_t j;
    int64_t maxNumberOfChains = 0;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'c' }, { "maxNumberOfChains", required_argument, 0, 'd' }, {
                "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:h", long_options,
                &option_index);

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
            case 'd':
                j = sscanf(optarg, "%" PRIi64 "", &maxNumberOfChains);
                assert(j == 1);
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    //assert(logLevelString == NULL || strcmp(logLevelString, "CRITICAL") == 0 || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
    assert(cactusDiskDatabaseString != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(
            cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Loop on the flowers, doing the reference genome (this process must be run bottom up)
    ///////////////////////////////////////////////////////////////////////////

    stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    preCacheNestedFlowers(cactusDisk, flowers);
    for(j = 0; j < stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);
        st_logInfo("Processing a flower\n");
        normalise(flower, maxNumberOfChains);
    }

    st_logInfo("Finished normalising the flowers\n");

    ///////////////////////////////////////////////////////////////////////////
    // Unload the parent flowers
    ///////////////////////////////////////////////////////////////////////////

    for(j = 0; j < stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);
        assert(flower != NULL);
        //flower_check(flower);
        flower_unloadParent(flower); //We have this line just in case we are loading the parent..
    }
    stList_destruct(flowers);

    st_logInfo("Unloaded the parent flowers\n");

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

    //Destruct stuff
    stKVDatabaseConf_destruct(kvDatabaseConf);
    if(logLevelString != NULL) {
        free(logLevelString);
    }
    free(cactusDiskDatabaseString);

    st_logInfo("Cleaned stuff up and am finished\n");

    //while(1);

    return 0;
}
