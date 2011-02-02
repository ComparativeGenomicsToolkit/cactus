/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include "cactus.h"

/*
 * Iterates through a tree an extends the cactus by creating new flowers for
 * non-terminal groups. Used during the core stage.
 */

static void addFlower(FILE *fileHandle, Flower *flower) {
    assert(!flower_builtBlocks(flower));
    assert(flower_isLeaf(flower));
    assert(flower_getBlockNumber(flower) == 0);
    assert(flower_getGroupNumber(flower) == 1);
    assert(flower_isTerminal(flower));
    int64_t size = flower_getTotalBaseLength(flower);
    assert(size >= 0);
    fprintf(fileHandle, "%s %" PRIi64 "\n", cactusMisc_nameToStringStatic(flower_getName(flower)), size);
}

static void extendFlowers(Flower *flower, FILE *fileHandle, int32_t minSizeToExtend, int32_t *flowersMade) {
    Flower_GroupIterator *groupIterator;
    Group *group;
    if (flower_builtBlocks(flower)) {
        groupIterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if (group_isLeaf(group)) {
                int64_t size = group_getTotalBaseLength(group);
                assert(size >= 0);
                if (size >= minSizeToExtend) {
                    group_makeNestedFlower(group);
                    addFlower(fileHandle, group_getNestedFlower(group));
                    (*flowersMade)++;
                }
            }
        }
        flower_destructGroupIterator(groupIterator);
    } else { //something went wrong last time, and the flower hasn't been filled in.. so we'll return it
        //again.
        addFlower(fileHandle, flower);
    }
}

int main(int argc, char *argv[]) {
    /*
     * This code iterates through the terminal groups and returns
     * a list of the new flowers.
     */
    assert(argc >= 5);

    if (strcmp(argv[1], "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    } else if (strcmp(argv[1], "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }
    st_logInfo("Set up logging\n");

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[2]);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    FILE *fileHandle = fopen(argv[3], "w");
    int32_t flowersMade = 0;

    int32_t minSizeToExtend;
    int32_t i = sscanf(argv[4], "%i", &minSizeToExtend);
    assert(i == 1); //We should parse one in correctly.
    st_logInfo("Min size to extend %i\n", minSizeToExtend);

    for (i = 5; i < argc; i++) {
        Flower *flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(argv[i]));
        assert(flower != NULL);
        st_logInfo("Parsed the flower: %s\n", argv[i]);
        extendFlowers(flower, fileHandle, minSizeToExtend, &flowersMade);
        assert(!flower_isParentLoaded(flower)); //The parent should not be loaded.
    }
    fclose(fileHandle);
    st_logInfo("We extended %i flowers\n", flowersMade);

    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flowerdisk\n");

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    st_logInfo("Am finished\n");
    return 0;
}
