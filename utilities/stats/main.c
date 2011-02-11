/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * main.c
 *
 *  Created on: 15-Apr-2010
 *      Author: benedictpaten
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "treeStats.h"

/*
 * Script gathers a whole gamut of statistics about the cactus/avg datastructure and reports them in an XML formatted document.
 */

void usage() {
    fprintf(stderr, "cactus_treeStats, version 0.1\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr,
            "-e --outputFile : The file to write the stats in, XML formatted.\n");
    fprintf(stderr,
            "-f --noPerColumnStats : Do not print out per column stats.\n");
    fprintf(stderr, "-g --includeSpecies : Species to require included\n");
    fprintf(stderr, "-i --excludedSpecies : Species to require excluded\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

static void buildSet(const char *string, stSortedSet *set) {
    /*
     * Breaks a string by white space into a bunch of tokens which are then added to the given set.
     */
    char *cA = stString_copy(string);
    char **cA2 = &cA;
    char *cA3;
    while((cA3 = stString_getNextWord(cA2)) != NULL) {
        stSortedSet_insert(set, cA3);
    }
    free(cA);
}

int main(int argc, char *argv[]) {
    /*
     * The script builds a cactus tree representation of the chains and flowers.
     * The format of the output graph is dot format.
     */

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = "0";
    char * outputFile = NULL;
    bool perColumnStats = 1;
    stSortedSet *includeSpecies = stSortedSet_construct3((int (*)(const void *, const void *))strcmp, free);
    stSortedSet *excludeSpecies = stSortedSet_construct3((int (*)(const void *, const void *))strcmp, free);

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'c' }, { "flowerName", required_argument, 0, 'd' }, {
                "outputFile", required_argument, 0, 'e' }, {
                "noPerColumnStats", no_argument, 0, 'f' },
                { "includeSpecies", required_argument, 0, 'g' },
                { "help", no_argument, 0, 'h' },
                { "excludeSpecies", required_argument, 0, 'i' },{ 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:e:fg:hi:", long_options,
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
                flowerName = stString_copy(optarg);
                break;
            case 'e':
                outputFile = stString_copy(optarg);
                break;
            case 'f':
                perColumnStats = 0;
                break;
            case 'h':
                usage();
                return 0;
            case 'g':
                buildSet(optarg, includeSpecies);
                break;
            case 'i':
                buildSet(optarg, excludeSpecies);
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
    assert(flowerName != NULL);
    assert(outputFile != NULL);

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

    st_logInfo("Flower name : %s\n", flowerName);
    st_logInfo("Output graph file : %s\n", outputFile);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(
            cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    Flower *flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(
            flowerName));
    assert(flower != NULL);
    st_logInfo("Parsed the top level flower of the cactus tree to build\n");

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    reportCactusDiskStats("EMPTY", flower, fileHandle, perColumnStats, includeSpecies, excludeSpecies);
    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
