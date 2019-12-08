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
#include "stMatchingAlgorithms.h"
#include "stReferenceProblem2.h"

void usage() {
    fprintf(stderr, "cactus_reference [flower names], version 0.1\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(
    stderr,
            "-e --matchingAlgorithm : Name of matching algorithm, either 'greedy', 'maxWeight', 'maxCardinality', 'blossom5'\n");
    fprintf(stderr, "-g --referenceEventString : String identifying the reference event.\n");
    fprintf(stderr, "-i --permutations : Number of permutations of gibss sampling, integer >= 0\n");
    fprintf(stderr, "-j --useSimulatedAnnealing : Use a cooling schedule\n");
    fprintf(stderr, "-k --theta : The value of theta, higher values are more tolerant of rearrangement.\n");
    fprintf(stderr,
            "-s --phi : The value of phi, value is proportional to exponential decline in support by phylogenetic distance \n");
    fprintf(stderr, "-l --maxWalkForCalculatingZ : The max number segments along a thread before stopping calculating z-scores\n");
    fprintf(
    stderr,
            "-m --ignoredUnalignedGaps : Don't consider unaligned sequence (gaps) when calculating the score function.\n");
    fprintf(
    stderr,
            "-n --wiggle : As the chains are being added to the reference greedily, wiggle is the proportion of next best possible insertion score the actual insertion has.\n");

    fprintf(
    stderr, "-o --numberOfNs : The number of N characters to place in a scaffold gap. Default=10. Must be >=1\n");

    fprintf(
    stderr, "-p --minNumberOfSequencesToSupportAdjacency : Min number of sequences to support an ancestral adjacency call. Default=1. Must be >=0\n");

    fprintf(
    stderr, "-q --makeScaffolds : Scaffold across regions of adjacency uncertainty.\n");

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
    stList *(*matchingAlgorithm)(stList *edges, int64_t nodeNumber) =
    chooseMatching_greedy;
    char *referenceEventString = (char *) cactusMisc_getDefaultReferenceEventHeader();
    int64_t permutations = 10;
    double theta = 0.001;
    double phi = 1.0;
    bool useSimulatedAnnealing = 0;
    int64_t maxWalkForCalculatingZ = 10000;
    bool ignoreUnalignedGaps = 0;
    double wiggle = 0.95;
    int64_t numberOfNsForScaffoldGap = 10;
    int64_t minNumberOfSequencesToSupportAdjacency = 1;
    bool makeScaffolds = 0;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
        required_argument, 0, 'a' }, { "cactusDisk", required_argument, 0, 'c' }, { "matchingAlgorithm",
        required_argument, 0, 'e' }, { "referenceEventString", required_argument, 0, 'g' }, { "permutations",
        required_argument, 0, 'i' }, { "useSimulatedAnnealing", no_argument, 0, 'j' }, { "theta",
        required_argument, 0, 'k' }, { "phi",
        required_argument, 0, 's' }, { "maxWalkForCalculatingZ", required_argument, 0, 'l' }, { "ignoreUnalignedGaps",
        no_argument, 0, 'm' }, { "wiggle", required_argument, 0, 'n' }, { "numberOfNs", required_argument, 0, 'o' }, {
                "minNumberOfSequencesToSupportAdjacency", required_argument, 0, 'p' }, { "makeScaffolds", no_argument,
                0, 'q' }, { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:e:g:i:jk:hl:mn:o:p:qs:", long_options, &option_index);

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
                matchingAlgorithm = chooseMatching_maximumCardinalityMatching;
            } else if (strcmp("maxWeight", optarg) == 0) {
                matchingAlgorithm = chooseMatching_maximumWeightMatching;
            } else if (strcmp("blossom5", optarg) == 0) {
                matchingAlgorithm = chooseMatching_blossom5;
            } else {
                stThrowNew(REFERENCE_BUILDING_EXCEPTION, "Input error: unrecognized matching algorithm: %s", optarg);
            }
            break;
        case 'g':
            referenceEventString = stString_copy(optarg);
            break;
        case 'h':
            usage();
            return 0;
        case 'i':
            j = sscanf(optarg, "%" PRIi64 "", &permutations);
            assert(j == 1);
            if (permutations < 0) {
                stThrowNew(REFERENCE_BUILDING_EXCEPTION, "Permutations is not valid %" PRIi64 "", permutations);
            }
            break;
        case 'j':
            useSimulatedAnnealing = 1;
            break;
        case 'k':
            j = sscanf(optarg, "%lf", &theta);
            assert(j == 1);
            if (theta < 0 || theta > 1.0) {
                stThrowNew(REFERENCE_BUILDING_EXCEPTION, "The theta parameter is not valid %f", theta);
            }
            break;
        case 's':
            j = sscanf(optarg, "%lf", &phi);
            assert(j == 1);
            if (phi <= 0.0) {
                stThrowNew(REFERENCE_BUILDING_EXCEPTION, "The phi parameter is not valid %f", phi);
            }
            break;
        case 'l':
            j = sscanf(optarg, "%" PRIi64 "", &maxWalkForCalculatingZ);
            assert(j == 1);
            if (maxWalkForCalculatingZ <= 0) {
                stThrowNew(REFERENCE_BUILDING_EXCEPTION, "The maxWalkForCalculatingZ is not valid %" PRIi64 "\n",
                        maxWalkForCalculatingZ);
            }
            break;
        case 'm':
            ignoreUnalignedGaps = 1;
            break;
        case 'n':
            j = sscanf(optarg, "%lf", &wiggle);
            assert(j == 1);
            if (wiggle < 0 || wiggle > 1.0) {
                stThrowNew(REFERENCE_BUILDING_EXCEPTION, "The wiggle parameter is not valid %f", wiggle);
            }
            break;
        case 'o':
            j = sscanf(optarg, "%" PRIi64 "", &numberOfNsForScaffoldGap);
            assert(j == 1);
            if (numberOfNsForScaffoldGap < 1) {
                stThrowNew(REFERENCE_BUILDING_EXCEPTION,
                        "Number of Ns for scaffold gap is not valid (must be >= 1): %" PRIi64 "",
                        numberOfNsForScaffoldGap);
            }
            break;
        case 'p':
            j = sscanf(optarg, "%" PRIi64 "", &minNumberOfSequencesToSupportAdjacency);
            assert(j == 1);
            if (minNumberOfSequencesToSupportAdjacency < 0) {
                stThrowNew(REFERENCE_BUILDING_EXCEPTION,
                        "minNumberOfSequencesToSupportAdjacency is not valid (must be >= 1): %" PRIi64 "",
                        minNumberOfSequencesToSupportAdjacency);
            }
            break;
        case 'q':
            makeScaffolds = 1;
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
    st_logInfo("The ignore unaligned gaps parameter is %i\n", ignoreUnalignedGaps);
    st_logInfo("The number of permutations is %" PRIi64 "\n", permutations);
    st_logInfo("Simulated annealing is %" PRIi64 "\n", useSimulatedAnnealing);
    st_logInfo("Max number of segments in thread to calculate z-score between is %" PRIi64 "\n",
            maxWalkForCalculatingZ);
    st_logInfo("Wiggle is %f\n", wiggle);
    st_logInfo("Max number of Ns for a scaffold gap is: %" PRIi64 "\n", numberOfNsForScaffoldGap);
    st_logInfo("Min number of sequences to required to support an adjacency is: %" PRIi64 "\n",
            minNumberOfSequencesToSupportAdjacency);
    st_logInfo("Make scaffolds is: %i\n", makeScaffolds);

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(cactusDiskDatabaseString != NULL);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, false, true);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Build the reference
    ///////////////////////////////////////////////////////////////////////////

    double (*temperatureFn)(double) =
    useSimulatedAnnealing ? exponentiallyDecreasingTemperatureFn
    : constantTemperatureFn;

    FlowerStream *flowerStream = flowerWriter_getFlowerStream(cactusDisk, stdin);
    Flower *flower;
    while ((flower = flowerStream_getNext(flowerStream)) != NULL) {
        st_logInfo("Processing flower %" PRIi64 "\n", flower_getName(flower));
        // Bulk-load all the nested flowers, since we will be reading them in anyway.
        // Currently we can only cache the nested flowers of a list of flowers, so we
        // create a singleton list.
        stList *flowers = stList_construct();
        stList_append(flowers, flower);
        preCacheNestedFlowers(cactusDisk, flowers);
        stList_destruct(flowers);

        if (!flower_hasParentGroup(flower)) {
            buildReferenceTopDown(flower, referenceEventString, permutations, matchingAlgorithm, temperatureFn, theta,
                    phi, maxWalkForCalculatingZ, ignoreUnalignedGaps, wiggle, numberOfNsForScaffoldGap,
                    minNumberOfSequencesToSupportAdjacency, makeScaffolds);
            cactusDisk_addUpdateRequest(cactusDisk, flower);
        }
        Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIt)) != NULL) {
            Flower *subFlower = group_getNestedFlower(group);
            if (subFlower != NULL) {
                buildReferenceTopDown(subFlower, referenceEventString, permutations,
                        matchingAlgorithm, temperatureFn, theta, phi, maxWalkForCalculatingZ, ignoreUnalignedGaps,
                        wiggle, numberOfNsForScaffoldGap, minNumberOfSequencesToSupportAdjacency, makeScaffolds);
                cactusDisk_addUpdateRequest(cactusDisk, subFlower);
                flower_unload(subFlower);
            }
        }
        flower_destructGroupIterator(groupIt);
        assert(!flower_isParentLoaded(flower));
        cactusDisk_clearCache(cactusDisk);
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
