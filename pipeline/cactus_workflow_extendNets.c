#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include "cactus.h"

/*
 * Iterates through a tree an extends the cactus by creating new nets for
 * non-terminal groups. Used during the core stage.
 */

static void extendNets(Flower *net, FILE *fileHandle, int32_t minSizeToExtend) {
    Flower_GroupIterator *groupIterator;
    Group *group;
    if(flower_builtBlocks(net)) {
        assert(net != NULL);
        groupIterator = flower_getGroupIterator(net);
        while((group = flower_getNextGroup(groupIterator)) != NULL) {
            if(!group_isLeaf(group)) {
                extendNets(group_getNestedFlower(group), fileHandle, minSizeToExtend);
            }
            else {
                int64_t size = group_getTotalBaseLength(group);
                if(size >= minSizeToExtend) {
                    group_makeNestedFlower(group);
                    fprintf(fileHandle, "%s %" PRIi64 "\n", cactusMisc_nameToStringStatic(group_getName(group)), size);
                }
            }
        }
        flower_destructGroupIterator(groupIterator);
    }
    else { //something went wrong last time, and the net hasn't been filled in.. so we'll return it
        //again.
        assert(flower_getBlockNumber(net) == 0);
        fprintf(fileHandle, "%s %" PRIi64 "\n", cactusMisc_nameToStringStatic(flower_getName(net)), flower_getTotalBaseLength(net));
    }
}

int main(int argc, char *argv[]) {
    /*
     * This code iterates through the terminal groups and returns
     * a list of the new nets.
     */
    CactusDisk *netDisk;
    Flower *net;

    assert(argc == 5);
    netDisk = cactusDisk_construct(argv[1]);
    st_logInfo("Set up the net disk\n");

    net = cactusDisk_getFlower(netDisk, cactusMisc_stringToName(argv[2]));
    st_logInfo("Parsed the net\n");

    int32_t minSizeToExtend;
    assert(sscanf(argv[4], "%i", &minSizeToExtend) == 1);

    FILE *fileHandle = fopen(argv[3], "w");
    extendNets(net, fileHandle, minSizeToExtend);
    fclose(fileHandle);

    cactusDisk_write(netDisk);
    st_logInfo("Updated the netdisk\n");

    cactusDisk_destruct(netDisk);

    st_logInfo("Am finished\n");
    return 0;
}
