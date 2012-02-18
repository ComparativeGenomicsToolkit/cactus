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

void makeMaf(Flower *flower, const char *referenceEventString,
        const char *childDirectory, const char *parentDirectory,
        bool showOnlySubstitutionsWithRespectToReference,
        const char *outputFile);

void usage() {
    fprintf(stderr, "cactus_mafGenerator [flower names], version 0.1\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-g --referenceEventString : String identifying the reference event.\n");
    fprintf(stderr, "-i --childDirectory : Directory identifying child mafs.\n");
    fprintf(stderr,
            "-j --parentDirectory : Directory identifying parent mafs.\n");
    fprintf(stderr, "-k --outputFile : File to maf in.\n");
    fprintf(
            stderr,
            "-l --showOnlySubstitutionsWithRespectToReference : Put stars in place of characters that are identical to the reference.\n");
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
    char *referenceEventString =
            (char *) cactusMisc_getDefaultReferenceEventHeader();
    char *childDirectory = NULL;
    char *parentDirectory = NULL;
    char *outputFile = NULL;
    bool showOnlySubstitutionsWithRespectToReference = 0;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'c' },
                { "referenceEventString", required_argument, 0, 'g' }, {
                        "help", no_argument, 0, 'h' }, { "childDirectory",
                        required_argument, 0, 'i' }, { "parentDirectory",
                        required_argument, 0, 'j' }, { "outputFile", required_argument, 0,
                        'k' }, { "showOnlySubstitutionsWithRespectToReference",
                        no_argument, 0, 'l' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:e:g:hi:j:k:l", long_options,
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
            case 'g':
                referenceEventString = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case 'i':
                childDirectory = stString_copy(optarg);
                break;
            case 'j':
                parentDirectory = stString_copy(optarg);
                break;
            case 'k':
                outputFile = stString_copy(optarg);
                break;
            case 'l':
                showOnlySubstitutionsWithRespectToReference = 1;
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
    // Get the set of flowers to manipulate
    ///////////////////////////////////////////////////////////////////////////

    stList *flowers = parseFlowersFromStdin(cactusDisk);

    ///////////////////////////////////////////////////////////////////////////
    // Now process each flower in turn.
    ///////////////////////////////////////////////////////////////////////////

    if (outputFile != NULL && stList_length(flowers) != 1) {
        stThrowNew("RUNTIME_ERROR",
                "Output file specified, but there is not only one flower required\n");
    }
    if (childDirectory == NULL) {
        stThrowNew("RUNTIME_ERROR", "No child directory specified\n");
    }
    if (parentDirectory == NULL && outputFile == NULL) {
        stThrowNew("RUNTIME_ERROR",
                "No parent directory or output file specified\n");
    }
    for (int32_t j = 0; j < stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);
        st_logInfo("Processing a flower\n");
        makeMaf(flower, referenceEventString, childDirectory, parentDirectory,
                showOnlySubstitutionsWithRespectToReference, outputFile);
    }

    ///////////////////////////////////////////////////////////////////////////
    //Clean up.
    ///////////////////////////////////////////////////////////////////////////

    //return 0; //Exit without clean up is quicker, enable cleanup when doing memory leak detection.

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    stList_destruct(flowers);
    free(cactusDiskDatabaseString);
    free(referenceEventString);
    free(logLevelString);
    free(childDirectory);
    free(parentDirectory);

    st_logInfo("Cleaned stuff up and am finished\n");
    //while(1);
    return 0;
}
