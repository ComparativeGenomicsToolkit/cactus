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
#include "addReferenceCoordinates.h"
#include "blockMLString.h"

void usage() {
    fprintf(stderr, "cactus_addReferenceCoordinates [flower names], version 0.1\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-c --secondaryDisk : The location of secondary disk\n");
    fprintf(stderr, "-g --referenceEventString : String identifying the reference event.\n");
    fprintf(stderr, "-j --bottomUpPhase : Do bottom up stage instead of top down.\n");
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
    char * cactusSequencesPath = NULL;
    char * secondaryDatabaseString = NULL;
    char *referenceEventString = (char *) cactusMisc_getDefaultReferenceEventHeader();
    bool bottomUpPhase = 0;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "cactusDisk", required_argument, 0, 'b' }, {"cactusSequencesPath", required_argument, 0, 'c'}, {
                "secondaryDisk", required_argument, 0, 'd' }, { "referenceEventString", required_argument, 0, 'g' }, { "help", no_argument,
                0, 'h' }, { "bottomUpPhase", no_argument, 0, 'j' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:e:g:hi:j", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'c':
                cactusSequencesPath = stString_copy(optarg);
                break;
            case 'd':
                secondaryDatabaseString = stString_copy(optarg);
                break;
            case 'g':
                referenceEventString = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case 'j':
                bottomUpPhase = 1;
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
    assert(cactusSequencesPath != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    st_logInfo("cactusSequencesPath = %s\n", cactusSequencesPath);
    st_logInfo("referenceEventString = %s", referenceEventString);
    st_logInfo("bottomUpPhase = %i", bottomUpPhase);

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct2(kvDatabaseConf, false, cactusSequencesPath, false);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    st_logInfo("Set up the flower disk\n");

    stKVDatabase *sequenceDatabase = NULL;
    if (secondaryDatabaseString != NULL) {
        kvDatabaseConf = stKVDatabaseConf_constructFromString(secondaryDatabaseString);
        sequenceDatabase = stKVDatabase_construct(kvDatabaseConf, 0);
        stKVDatabaseConf_destruct(kvDatabaseConf);
    }

    FlowerStream *flowerStream = flowerWriter_getFlowerStream(cactusDisk, stdin);
    Flower *flower;
    while ((flower = flowerStream_getNext(flowerStream)) != NULL) {
        st_logDebug("Processing flower %" PRIi64 "\n", flower_getName(flower));

        ///////////////////////////////////////////////////////////////////////////
        // Get the appropriate event names
        ///////////////////////////////////////////////////////////////////////////

        Event *referenceEvent = eventTree_getEventByHeader(flower_getEventTree(flower), referenceEventString);
        assert(referenceEvent != NULL);
        Name referenceEventName = event_getName(referenceEvent);

        ///////////////////////////////////////////////////////////////////////////
        // Now do bottom up or top down, depending
        ///////////////////////////////////////////////////////////////////////////
        stList *flowers = stList_construct();
        stList_append(flowers, flower);
        preCacheNestedFlowers(cactusDisk, flowers);
        if (bottomUpPhase) {
            assert(sequenceDatabase != NULL);

            cactusDisk_preCacheSegmentStrings(cactusDisk, flowers);
            bottomUp(flowers, sequenceDatabase, referenceEventName, !flower_hasParentGroup(flower), generateJukesCantorMatrix);

            // Unload the nested flowers to save memory. They haven't
            // been changed, so we don't write them to the cactus
            // disk.
            Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
            Group *group;
            while ((group = flower_getNextGroup(groupIt)) != NULL) {
                if (!group_isLeaf(group)) {
                    flower_unload(group_getNestedFlower(group));
                }
            }
            flower_destructGroupIterator(groupIt);
            assert(!flower_isParentLoaded(flower));

            // Write this flower to disk.
            cactusDisk_addUpdateRequest(cactusDisk, flower);

            // Unload the sequences from memory.
            cactusDisk_clearStringCache(cactusDisk);
        } else {
            topDown(flower, referenceEventName);

            // We've changed the nested flowers, but not this
            // flower. We write the nested flowers to disk, then
            // unload them to save memory. This flower will be
            // unloaded by the flower-stream code.
            Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
            Group *group;
            while ((group = flower_getNextGroup(groupIt)) != NULL) {
                if (!group_isLeaf(group)) {
                    cactusDisk_addUpdateRequest(cactusDisk, group_getNestedFlower(group));
                    flower_unload(group_getNestedFlower(group));
                }
            }
            flower_destructGroupIterator(groupIt);
        }
        stList_destruct(flowers);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Write the flower(s) back to disk.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk\n");

    ///////////////////////////////////////////////////////////////////////////
    //Clean up.
    ///////////////////////////////////////////////////////////////////////////

    if (sequenceDatabase != NULL) {
        stKVDatabase_destruct(sequenceDatabase);
    }

    cactusDisk_destruct(cactusDisk);

    return 0; //Exit without clean up is quicker, enable cleanup when doing memory leak detection.

    free(cactusDiskDatabaseString);
    free(referenceEventString);
    free(logLevelString);

    st_logInfo("Cleaned stuff up and am finished\n");

    return 0;
}
