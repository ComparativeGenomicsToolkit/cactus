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
#include "cactusMatchingAlgorithms.h"
#include "adjacencyProblem.h"

void usage() {
    fprintf(stderr, "cactus_reference [flower names], version 0.1\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(
            stderr,
            "-e --matchingAlgorithm : Name of matching algorithm, either 'greedy', 'maxWeight', 'maxCardinality', 'blossom5'\n");
    fprintf(stderr,
            "-g --referenceEventString : String identifying the reference event.\n");
    fprintf(stderr,
            "-i --permutations : Number of permutations of gibss sampling, integer >= 0\n");
    fprintf(stderr, "-j --useSimulatedAnnealing : Use a cooling schedule\n");
    fprintf(stderr, "-k --theta : The value of theta\n");
    fprintf(
            stderr,
            "-l --maxNumberOfChainsBeforeSwitchingToFast : The max number of chains before switching to fast reference building mode (default INT32_MAX)\n");

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
    int32_t j;
    stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber) =
            chooseMatching_greedy;
    char *referenceEventString =
            (char *) cactusMisc_getDefaultReferenceEventHeader();
    int32_t permutations = 10;
    double theta = 0.00001;
    bool useSimulatedAnnealing = 0;
    int32_t maxNumberOfChainsBeforeSwitchingToFast = INT32_MAX;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'c' }, { "matchingAlgorithm", required_argument, 0, 'e' }, {
                "referenceEventString", required_argument, 0, 'g' }, {
                "permutations", required_argument, 0, 'i' }, {
                "useSimulatedAnnealing", no_argument, 0, 'j' }, { "theta",
                required_argument, 0, 'k' }, { "maxNumberOfChainsBeforeSwitchingToFast", required_argument, 0,
                'l' }, { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:e:g:i:jk:hl:", long_options,
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
            case 'e':
                if (strcmp("greedy", optarg) == 0) {
                    matchingAlgorithm = chooseMatching_greedy;
                } else if (strcmp("maxCardinality", optarg) == 0) {
                    matchingAlgorithm
                            = chooseMatching_maximumCardinalityMatching;
                } else if (strcmp("maxWeight", optarg) == 0) {
                    matchingAlgorithm = chooseMatching_maximumWeightMatching;
                } else if (strcmp("blossom5", optarg) == 0) {
                    matchingAlgorithm = chooseMatching_blossom5;
                } else {
                    stThrowNew(REFERENCE_BUILDING_EXCEPTION,
                            "Input error: unrecognized matching algorithm: %s",
                            optarg);
                }
                break;
            case 'g':
                referenceEventString = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case 'i':
                j = sscanf(optarg, "%i", &permutations);
                assert(j == 1);
                if (permutations < 0) {
                    stThrowNew(REFERENCE_BUILDING_EXCEPTION,
                            "Permutations is not valid %i", permutations);
                }
                break;
            case 'j':
                useSimulatedAnnealing = 1;
                break;
            case 'k':
                j = sscanf(optarg, "%lf", &theta);
                assert(j == 1);
                if (theta < 0 || theta > 1.0) {
                    stThrowNew(REFERENCE_BUILDING_EXCEPTION,
                            "The theta parameter is not valid %f", theta);
                }
                break;
            case 'l':
                j = sscanf(optarg, "%i", &maxNumberOfChainsBeforeSwitchingToFast);
                assert(j == 1);
                break;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    st_logInfo("The theta parameter has been set to %lf\n", theta);
    st_logInfo("The number of permutations is %i\n", permutations);
    st_logInfo("Simulated annealing is %i\n", useSimulatedAnnealing);

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(cactusDiskDatabaseString != NULL);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(
            cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Build the reference
    ///////////////////////////////////////////////////////////////////////////

    double (*temperatureFn)(double) =
            useSimulatedAnnealing ? exponentiallyDecreasingTemperatureFn
                    : constantTemperatureFn;

    for (j = optind; j < argc; j++) {
        const char *flowerName = argv[j];
        st_logInfo("Processing the flower named: %s\n", flowerName);
        Flower *flower = cactusDisk_getFlower(cactusDisk,
                cactusMisc_stringToName(flowerName));
        assert(flower != NULL);
        st_logInfo("Parsed the flower in which to build a reference\n");
        if (!flower_hasParentGroup(flower)) {
            buildReferenceTopDown(flower, referenceEventString, permutations,
                    matchingAlgorithm, temperatureFn, theta, maxNumberOfChainsBeforeSwitchingToFast);
        }
        Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIt)) != NULL) {
            if (group_getNestedFlower(group) != NULL) {
                buildReferenceTopDown(group_getNestedFlower(group),
                        referenceEventString, permutations, matchingAlgorithm,
                        temperatureFn, theta, maxNumberOfChainsBeforeSwitchingToFast);
            }
        }
        flower_destructGroupIterator(groupIt);
        flower_unloadParent(flower); //We have this line just in case we are loading the parent..
    }

    ///////////////////////////////////////////////////////////////////////////
    // Write the flower(s) back to disk.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk\n");

    ///////////////////////////////////////////////////////////////////////////
    //Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    free(cactusDiskDatabaseString);
    if (logLevelString != NULL) {
        free(logLevelString);
    }

    st_logInfo("Cleaned stuff up and am finished\n");

    //while(1);

    return 0;
}
