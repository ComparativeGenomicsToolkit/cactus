/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus_workflow_shared.h"

int main(int argc, char *argv[]) {
    parseArgs(argc, argv);
    FlowerStream *flowerStream = flowerWriter_getFlowerStream(cactusDisk, stdin);
    Flower *flower;
    while ((flower = flowerStream_getNext(flowerStream)) != NULL) {
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
    flowerStream_destruct(flowerStream);
    flowerWriter_destruct(flowerWriter);
    cactusDisk_destruct(cactusDisk);

    st_logDebug("Am finished\n");
    return 0;
}
