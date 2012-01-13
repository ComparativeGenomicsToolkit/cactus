/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus_workflow_shared.h"

int main(int argc, char *argv[]) {
    parseArgs(argc, argv);
    int32_t includeTerminalFlowers;
    int32_t i = sscanf(argv[6], "%i", &includeTerminalFlowers);
    assert(i == 1);
    st_logDebug("Keeping the terminal flowers: %i\n", includeTerminalFlowers);
    stList *flowers = parseFlowers(argv + 7, argc - 7, cactusDisk);
    for (int32_t i = 0; i < stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        if(!flower_isLeaf(flower)) {
            assert(flower_builtBlocks(flower)); //This recursion depends on the block structure having been properly defined for all nodes.
            Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
            Group *group;
            while ((group = flower_getNextGroup(groupIterator)) != NULL) {
                if (!group_isLeaf(group)) {
                    Flower *nestedFlower = group_getNestedFlower(group);
                    assert(nestedFlower != NULL);
                    assert(flower_builtBlocks(nestedFlower)); //This recursion depends on the block structure having been properly defined for all nodes.
                    int64_t flowerSize = flower_getTotalBaseLength(nestedFlower);
                    if((includeTerminalFlowers || !flower_isTerminal(nestedFlower)) && flowerSize <= maxFlowerSize && flowerSize >= minFlowerSize) {
                        flowerWriter_add(flowerWriter, nestedFlower);
                    }
                }
            }
            flower_destructGroupIterator(groupIterator);
        }
    }
    stList_destruct(flowers);
    flowerWriter_destruct(flowerWriter);
    cactusDisk_destruct(cactusDisk);
    st_logDebug("Am finished\n");
    return 0;
}
