/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus_workflow_shared.h"

/*
 * Iterates through a tree an extends the cactus by creating new flowers for
 * non-terminal groups. Used during the core stage.
 */

static void extendFlowers(Flower *flower, bool createRedundantFlowerLinks) {
    Flower_GroupIterator *groupIterator;
    Group *group;
    if (flower_builtBlocks(flower)) {
        groupIterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if (group_isLeaf(group)) {
                int64_t size = group_getTotalBaseLength(group);
                assert(size >= 0);
                if (size >= minFlowerSize) {
                    Flower *nestedFlower = group_makeNestedFlower(group);
                    if (createRedundantFlowerLinks) {
                        flower_setBuiltBlocks(nestedFlower, 1);
                        Group *nestedGroup = flower_getFirstGroup(nestedFlower);
                        nestedFlower = group_makeNestedFlower(nestedGroup);
                        flowerWriter_add(flowerWriter, flower_getName(nestedFlower), size);
                    } else {
                        flowerWriter_add(flowerWriter, flower_getName(nestedFlower), size);
                    }
                }
            }
        }
        flower_destructGroupIterator(groupIterator);
    } else { //We are at the top of the hierarchy.
        assert(flower_getName(flower) == 0);
        flowerWriter_add(flowerWriter, flower_getName(flower), flower_getTotalBaseLength(flower));
    }
}

static bool flowerIsLarge(Flower *flower) {
    return flower_getEndNumber(flower) > 1000;
}

int main(int argc, char *argv[]) {
    /*
     * This code iterates through the terminal groups and returns
     * a list of the new flowers.
     */
    parseArgs(argc, argv);

    st_logDebug("Starting extending flowers\n");
    stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    for (int64_t i = 0; i < stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        extendFlowers(flower, flowerIsLarge(flower));
        assert(!flower_isParentLoaded(flower)); //The parent should not be loaded.
    }
    stList_destruct(flowers);
    st_logDebug("Finish extending flowers\n");

    flowerWriter_destruct(flowerWriter);

    cactusDisk_write(cactusDisk);
    st_logDebug("Updated the cactus disk\n");
    cactusDisk_destruct(cactusDisk);

    st_logDebug("Am finished\n");
    return 0;
}
