/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus_workflow_shared.h"

int main(int argc, char *argv[]) {
    parseArgs(argc, argv);
    stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    //stList *flowers = parseFlowers(argv + 6, argc - 6, cactusDisk);
    for (int64_t i = 0; i < stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        if(!flower_isLeaf(flower)) {
            assert(flower_builtBlocks(flower)); //This recursion depends on the block structure having been properly defined for all nodes.
            Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
            Group *group;
            while ((group = flower_getNextGroup(groupIterator)) != NULL) {
                if (!group_isLeaf(group)) {
                    int64_t flowerSize = group_getTotalBaseLength(group);
                    if(flowerSize >= minFlowerSize) {
                        flowerWriter_add(flowerWriter, group_getName(group), flowerSize);
                    }
                }
            }
            flower_destructGroupIterator(groupIterator);
        }
    }
    flowerWriter_destruct(flowerWriter);
    cactusDisk_destruct(cactusDisk);

    return 0; //Avoid cleanup
    stList_destruct(flowers);
    st_logDebug("Am finished\n");
    return 0;
}
