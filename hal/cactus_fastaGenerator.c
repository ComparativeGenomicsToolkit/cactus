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
#include "hal.h"

void usage() {
    fprintf(stderr, "cactus_fastaGenerator [flower names], version 0.1\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-g --referenceEventString : String identifying the reference event.\n");
    fprintf(stderr,
                "-d --flowerName : Name of flower to print string for.\n");
    fprintf(stderr, "-k --outputFile : File to put final output in.\n");
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
    char *outputFile = NULL;
    Name flowerName = NULL_NAME;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'c' },  { "flowerName", required_argument,
                        0, 'd' },
                { "referenceEventString", required_argument, 0, 'g' }, {
                        "help", no_argument, 0, 'h' }, { "outputFile",
                        required_argument, 0, 'k' },
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:g:hk:", long_options,
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
                flowerName = cactusMisc_stringToName(optarg);
                break;
            case 'g':
                referenceEventString = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case 'k':
                outputFile = stString_copy(optarg);
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
    stKVDatabaseConf_destruct(kvDatabaseConf);
    st_logInfo("Set up the flower disk\n");


    ///////////////////////////////////////////////////////////////////////////
    // Get the set of flowers to manipulate
    ///////////////////////////////////////////////////////////////////////////

    Flower *flower = cactusDisk_getFlower(cactusDisk, flowerName);

    ///////////////////////////////////////////////////////////////////////////
    // Get the reference event name
    ///////////////////////////////////////////////////////////////////////////

    Event *referenceEvent = eventTree_getEventByHeader(
            flower_getEventTree(flower), referenceEventString);
    assert(referenceEvent != NULL);
    Name referenceEventName = event_getName(referenceEvent);

    ///////////////////////////////////////////////////////////////////////////
    // Now process each flower in turn.
    ///////////////////////////////////////////////////////////////////////////
    
    if(outputFile == NULL) {
        st_errAbort("No output file specified\n");
    }
    FILE *fileHandle = fopen(outputFile, "w");
    printFastaSequences(flower, fileHandle, referenceEventName);
    if(fileHandle != NULL) {
        fclose(fileHandle);
    }

    ///////////////////////////////////////////////////////////////////////////
    //Clean up memory
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);

    //return 0; //Exit without clean up is quicker, enable cleanup when doing memory leak detection.

    free(cactusDiskDatabaseString);
    free(referenceEventString);
    free(logLevelString);

    st_logInfo("Cleaned stuff up and am finished\n");
    //while(1);
    return 0;
}
