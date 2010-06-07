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

static void extendNets(Net *net, FILE *fileHandle, int32_t minSizeToExtend) {
    Net_GroupIterator *groupIterator;
    Group *group;
    if(net_builtBlocks(net)) {
        assert(net != NULL);
        groupIterator = net_getGroupIterator(net);
        while((group = net_getNextGroup(groupIterator)) != NULL) {
            if(!group_isLeaf(group)) {
                extendNets(group_getNestedNet(group), fileHandle, minSizeToExtend);
            }
            else {
                int64_t size = group_getTotalBaseLength(group);
                if(size >= minSizeToExtend) {
                    group_makeNestedNet(group);
                    fprintf(fileHandle, "%s %" PRIi64 "\n", netMisc_nameToStringStatic(group_getName(group)), size);
                }
            }
        }
        net_destructGroupIterator(groupIterator);
    }
    else { //something went wrong last time, and the net hasn't been filled in.. so we'll return it
        //again.
        assert(net_getBlockNumber(net) == 0);
        fprintf(fileHandle, "%s %" PRIi64 "\n", netMisc_nameToStringStatic(net_getName(net)), net_getTotalBaseLength(net));
    }
}

int main(int argc, char *argv[]) {
    /*
     * This code iterates through the terminal groups and returns
     * a list of the new nets.
     */
    NetDisk *netDisk;
    Net *net;

    assert(argc == 5);
    netDisk = netDisk_construct(argv[1]);
    st_logInfo("Set up the net disk\n");

    net = netDisk_getNet(netDisk, netMisc_stringToName(argv[2]));
    st_logInfo("Parsed the net\n");

    int32_t minSizeToExtend;
    assert(sscanf(argv[4], "%i", &minSizeToExtend) == 1);

    FILE *fileHandle = fopen(argv[3], "w");
    extendNets(net, fileHandle, minSizeToExtend);
    fclose(fileHandle);

    netDisk_write(netDisk);
    st_logInfo("Updated the netdisk\n");

    netDisk_destruct(netDisk);

    st_logInfo("Am finished\n");
    return 0;
}
