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
#include "commonC.h"

/*
 * The script checks the flowers are structured as we expect, essentially by
 * calling flower_check for each flower in the tree. We do a couple of other tests also for the normalisation phase.
 */

/*
 * Normalisation phase. We haven't made these checks part of the api checks yet, as I don't know if we'll keep them
 * forever.
 */

static void checkTreeIsTerminalNormalised(Flower *flower) {
    /*
     * A cactus tree is terminally normalised if all leaf flowers are terminal.
     */
    if (flower_isLeaf(flower)) {
        assert(flower_getBlockNumber(flower) == 0);
        assert(flower_isTerminal(flower));
        if(flower_getGroupNumber(flower) == 0) {
            assert(!flower_hasParentGroup(flower));
        }
        else {
            assert(flower_getGroupNumber(flower) == 1);
        }
        //The following are defensive checks (in that they are implied by being a leaf)
        Group *group;
        Flower_GroupIterator *iterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(iterator)) != NULL) {
            assert(group_isLeaf(group));
        }
        flower_destructGroupIterator(iterator);
    } else {
        Group *group;
        Flower_GroupIterator *iterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(iterator)) != NULL) {
            assert(!group_isLeaf(group));
        }
        flower_destructGroupIterator(iterator);
    }
}

static void checkChainsAreMaximal(Flower *flower) {
    /*
     * Checks that each chain is maximal and consistently named..
     */
    Group *parentGroup = flower_getParentGroup(flower);
    if (parentGroup != NULL) {
        End *end;
        Flower_EndIterator *endIterator = flower_getEndIterator(flower);
        while ((end = flower_getNextEnd(endIterator)) != NULL) {
            assert(end_getOrientation(end));
            if (end_isStubEnd(end) && end_isAttached(end)) { //is an attached stub end (inherited from a higher level)
                Link *link = group_getLink(end_getGroup(end));
                if (link != NULL) { //then the flower must be a terminal flower and the link is a copy of the one in the parent..
                    Chain *chain = link_getChain(link);
                    assert(flower_isLeaf(flower));
                    assert(chain_getLength(chain) == 1);
                }
            }
        }
        flower_destructEndIterator(endIterator);
    }
}

static void checkBlocksAreMaximal(Flower *flower) {
    /*
     * Checks each block is maximal.
     */
    assert(flower != NULL);
    End *end;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        assert(end_getOrientation(end));
        if (end_isBlockEnd(end)) { //is an attached stub end (inherited from a higher level)
            Link *link = group_getLink(end_getGroup(end));
            if (link != NULL) { //then the flower must be a terminal flower and the link is a copy of the one in the parent..
                assert(!link_isTrivial(link));
            }
        }
    }
    flower_destructEndIterator(endIterator);
}

static void checkFlowerIsNotRedundant(Flower *flower) {
    /*
     * Checks that if the flower is not a leaf or the root that it contains blocks.
     */
    assert(flower_builtBlocks(flower));
    if (flower_hasParentGroup(flower) && !flower_isLeaf(flower)) {
        assert(flower_getBlockNumber(flower) > 0);
    }
}

/*
 * Random other checks.
 */

static void checkBasesAccountedFor(Flower *flower) {
    /*
     * Checks all the bases in a flower end up in child flower or a nested flower.
     */
    int64_t totalBases = flower_getTotalBaseLength(flower);
    int64_t blockBases = 0.0;
    int64_t childBases = 0.0;
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    Block_InstanceIterator *segmentIterator;
    Segment *segment;
    while ((block = flower_getNextBlock(blockIterator)) != NULL) {
        segmentIterator = block_getInstanceIterator(block);
        while ((segment = block_getNext(segmentIterator)) != NULL) {
            if (segment_getSequence(segment) != NULL) {
                blockBases += segment_getLength(segment);
            }
        }
        block_destructInstanceIterator(segmentIterator);
    }
    flower_destructBlockIterator(blockIterator);
    Flower_GroupIterator *iterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(iterator)) != NULL) {
        int64_t size = (int64_t) group_getTotalBaseLength(group);
        if (group_getNestedFlower(group) != NULL) {
            assert(!group_isLeaf(group));
            assert(flower_getTotalBaseLength(group_getNestedFlower(group)) == size);
        } else {
            assert(group_isLeaf(group));
        }
        assert(size >= 0);
        childBases += size;
    }
    flower_destructGroupIterator(iterator);
    if (blockBases + childBases != totalBases) {
        fprintf(stderr,
                "Got %i block bases, %i childBases and %i total bases\n",
                (int) blockBases, (int) childBases, (int) totalBases);
    }
    assert(blockBases + childBases == totalBases);
}

static void checkFlower(Flower *flower) {
    flower_check(flower);
    flower_checkNotEmpty(flower, 0);
    checkBasesAccountedFor(flower);
    //Normalisation checks..
    checkTreeIsTerminalNormalised(flower);
    checkChainsAreMaximal(flower);
    checkBlocksAreMaximal(flower);
    checkFlowerIsNotRedundant(flower);
}

static void checkFlowers(Flower *flower, int32_t recursive) {
    if(!flower_hasParentGroup(flower)) { //Only check if it has no parent
        checkFlower(flower);
    }

    Group *group;
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    while ((group = flower_getNextGroup(groupIt)) != NULL) { //We only check the children, to avoid constructing
        //the parent of the flower.
        if (!group_isLeaf(group)) {
            checkFlower(group_getNestedFlower(group));
            if (recursive) {
                checkFlowers(group_getNestedFlower(group), 1);
            }
        }
    }
    flower_destructGroupIterator(groupIt);
}

void usage() {
    fprintf(stderr, "cactus_tree, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-e --recursive : Check all flowers recursively\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    CactusDisk *cactusDisk;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    int32_t recursive = 0;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'c' }, { "recursive", no_argument, 0, 'e' }, { "help",
                no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key =
                getopt_long(argc, argv, "a:c:eh", long_options, &option_index);

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
                recursive = 1;
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

    assert(cactusDiskDatabaseString != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    stList *flowers = parseFlowers(argv + optind, argc - optind, cactusDisk);
    for(int32_t j = 0; j < stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);
        st_logInfo("Processing a flower\n");

        ///////////////////////////////////////////////////////////////////////////
        // Recursive check the flowers.
        ///////////////////////////////////////////////////////////////////////////

        checkFlowers(flower, recursive);
        st_logInfo("Checked the flower/\n");
    }

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    stList_destruct(flowers);
    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
