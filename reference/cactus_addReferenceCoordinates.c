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
    fprintf(stderr, "-i --outgroupEventString : String identifying the reference event.\n");
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
    char * secondaryDatabaseString = NULL;
    char *referenceEventString = (char *) cactusMisc_getDefaultReferenceEventHeader();
    char *outgroupEventString = NULL;
    bool bottomUpPhase = 0;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "cactusDisk", required_argument, 0, 'c' }, {
                "secondaryDisk", required_argument, 0, 'd' }, { "referenceEventString", required_argument, 0, 'g' }, { "help", no_argument,
                0, 'h' }, { "outgroupEventString", required_argument, 0, 'i' }, { "bottomUpPhase", no_argument, 0, 'j' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:e:g:hi:j", long_options, &option_index);

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
                secondaryDatabaseString = stString_copy(optarg);
                break;
            case 'g':
                referenceEventString = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case 'i':
                outgroupEventString = stString_copy(optarg);
                break;
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

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Get the set of flowers to manipulate
    ///////////////////////////////////////////////////////////////////////////

    stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    preCacheNestedFlowers(cactusDisk, flowers);

    ///////////////////////////////////////////////////////////////////////////
    // Get the appropriate event names
    ///////////////////////////////////////////////////////////////////////////

    Flower *flower = stList_peek(flowers);
    Event *referenceEvent = eventTree_getEventByHeader(flower_getEventTree(flower), referenceEventString);
    assert(referenceEvent != NULL);
    Name referenceEventName = event_getName(referenceEvent);

    Name outgroupEventName = NULL_NAME;
    if (outgroupEventString != NULL) {
        Event *outgroupEvent = eventTree_getEventByHeader(flower_getEventTree(flower), outgroupEventString);
        assert(outgroupEvent != NULL);
        assert(event_isOutgroup(outgroupEvent));
        outgroupEventName = event_getName(outgroupEvent);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Now do bottom up or top down, depending
    ///////////////////////////////////////////////////////////////////////////

    if (bottomUpPhase) {
        cactusDisk_preCacheSegmentStrings(cactusDisk, flowers);
        assert(secondaryDatabaseString != NULL);
        kvDatabaseConf = stKVDatabaseConf_constructFromString(secondaryDatabaseString);
        stKVDatabase *sequenceDatabase = stKVDatabase_construct(kvDatabaseConf, 0);
        stKVDatabaseConf_destruct(kvDatabaseConf);
        Flower *flower = stList_get(flowers, 0);
        bottomUp(flowers, sequenceDatabase, referenceEventName, outgroupEventName, !flower_hasParentGroup(flower), generateJukesCantorMatrix);
        //Now unload the nested flowers.
        for (int64_t i = 0; i < stList_length(flowers); i++) {
            Flower *flower = stList_get(flowers, i);
            //Now ensure that the nested flowers are not loaded, as this will avoid writing them to disk
            Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
            Group *group;
            while ((group = flower_getNextGroup(groupIt)) != NULL) {
                if (!group_isLeaf(group)) {
                    flower_unload(group_getNestedFlower(group));
                }
            }
            flower_destructGroupIterator(groupIt);
            assert(!flower_isParentLoaded(flower));
        }
        stKVDatabase_destruct(sequenceDatabase);
    } else {
        topDown(flowers, referenceEventName);
        //Now ensure the non-nested flowers are not loaded
        for (int64_t i = 0; i < stList_length(flowers); i++) {
           flower_unload(stList_get(flowers, i)); //We haven't changed the parents
        }
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

    stList_destruct(flowers);
    free(cactusDiskDatabaseString);
    free(referenceEventString);
    free(logLevelString);

    st_logInfo("Cleaned stuff up and am finished\n");
    //while(1);
    return 0;
}
