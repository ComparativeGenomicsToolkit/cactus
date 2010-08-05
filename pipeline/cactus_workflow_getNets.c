#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include "cactus.h"

int main(int argc, char *argv[]) {
    CactusDisk *netDisk;
    Flower *net;

    assert(argc >= 3);
    netDisk = cactusDisk_construct(argv[1]);
    st_logInfo("Set up the net disk\n");

    FILE *fileHandle = fopen(argv[2], "w");

    int32_t i;
    for (i = 3; i < argc; i++) {
        net = cactusDisk_getNet(netDisk, cactusMisc_stringToName(argv[i]));
        assert(net != NULL);
        assert(flower_builtBlocks(net)); //This recursion depends on the block structure having been properly defined for all nodes.
        st_logInfo("Parsed the net %s\n", argv[i]);
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(net);
        Group *group;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if (!group_isLeaf(group)) {
                Flower *nestedNet = group_getNestedNet(group);
                assert(nestedNet != NULL);
                assert(flower_builtBlocks(nestedNet)); //This recursion depends on the block structure having been properly defined for all nodes.
                fprintf(fileHandle, "%s %" PRIi64 " \n", cactusMisc_nameToStringStatic(flower_getName(nestedNet)), flower_getTotalBaseLength(nestedNet));
            }
        }
        flower_destructGroupIterator(groupIterator);
    }
    fclose(fileHandle);
    cactusDisk_destruct(netDisk);
    st_logInfo("Am finished\n");
    return 0;
}
