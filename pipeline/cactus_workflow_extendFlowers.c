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

static void extendFlowers(Flower *flower, FILE *fileHandle, int32_t minSizeToExtend, int32_t *flowersMade) {
    Flower_GroupIterator *groupIterator;
    Group *group;
    if (flower_builtBlocks(flower)) {
        assert(flower != NULL);
        groupIterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if (!group_isLeaf(group)) {
                extendFlowers(group_getNestedFlower(group), fileHandle, minSizeToExtend, flowersMade);
            } else {
                int64_t size = group_getTotalBaseLength(group);
                assert(size >= 0);
                if (size >= minSizeToExtend) {
                    group_makeNestedFlower(group);
                    (*flowersMade)++;
                    fprintf(fileHandle, "%s %" PRIi64 "\n", cactusMisc_nameToStringStatic(group_getName(group)), size);
                }
            }
        }
        flower_destructGroupIterator(groupIterator);
    }
    else { //something went wrong last time, and the flower hasn't been filled in.. so we'll return it
        //again.
        assert(flower_getBlockNumber(flower) == 0);
        assert(flower_getGroupNumber(flower) == 1);
        int64_t size = flower_getTotalBaseLength(flower);
        assert(size >= 0);
        (*flowersMade)++;
        fprintf(fileHandle, "%s %" PRIi64 "\n", cactusMisc_nameToStringStatic(flower_getName(flower)), size);
    }
}

int main(int argc, char *argv[]) {
    /*
     * This code iterates through the terminal groups and returns
     * a list of the new flowers.
     */
    assert(argc == 6);

    if (strcmp(argv[5], "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    else if (strcmp(argv[5], "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }
    st_logInfo("Set up logging\n");

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[1]);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    Flower *flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(argv[2]));
    assert(flower != NULL);
    st_logInfo("Parsed the flower\n");

    int32_t minSizeToExtend;
    assert(sscanf(argv[4], "%i", &minSizeToExtend) == 1);

    FILE *fileHandle = fopen(argv[3], "w");
    int32_t flowersMade = 0;
    extendFlowers(flower, fileHandle, minSizeToExtend, &flowersMade);
    fclose(fileHandle);
    st_logInfo("We extended %i flowers\n", flowersMade);

    assert(!flower_isParentLoaded(flower));
    //flower_unloadParent(flower); //The parent should not be loaded.

    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flowerdisk\n");

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    st_logInfo("Am finished\n");
    return 0;
}
