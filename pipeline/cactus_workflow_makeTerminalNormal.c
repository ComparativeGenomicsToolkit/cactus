#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cactus.h"

/*
 * Checks if a net contains terminal groups but is non-terminal itself, it
 * then makes these groups non-terminal,
 */

int main(int argc, char *argv[]) {
    /*
     * This code iterates through the leaf groups.
     */
    CactusDisk *netDisk;
    Net *net;

    assert(argc >= 2);
    netDisk = cactusDisk_construct(argv[1]);
    st_logInfo("Set up the net disk\n");
    int32_t i;
    for (i = 2; i < argc; i++) {
        net = cactusDisk_getNet(netDisk, cactusMisc_stringToName(argv[i]));
        assert(net != NULL);
        st_logInfo("Parsed net %s\n", argv[i]);
        if (!net_isTerminal(net)) {
            Net_GroupIterator *groupIterator;
            Group *group;
            groupIterator = net_getGroupIterator(net);
            while ((group = net_getNextGroup(groupIterator)) != NULL) {
                if (group_isLeaf(group)) {
                    //assert(group_getTotalBaseLength(group) == 0);
                    group_makeNestedNet(group);
                    net_setBuiltBlocks(group_getNestedNet(group), 1);
                }
            }
            net_destructGroupIterator(groupIterator);
        }
    }

    cactusDisk_write(netDisk);
    st_logInfo("Updated the netdisk\n");

    cactusDisk_destruct(netDisk);

    st_logInfo("Am finished\n");
    return 0;
}
