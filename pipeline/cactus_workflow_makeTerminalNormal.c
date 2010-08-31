#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cactus.h"

/*
 * Checks if a flower contains terminal groups but is non-terminal itself, it
 * then makes these groups non-terminal,
 */

int main(int argc, char *argv[]) {
    /*
     * This code iterates through the leaf groups.
     */
    assert(argc >= 2);
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[1]);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");
    int32_t i;
    for (i = 2; i < argc; i++) {
        Flower *flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(argv[i]));
        assert(flower != NULL);
        st_logInfo("Parsed flower %s\n", argv[i]);
        if (!flower_isTerminal(flower)) {
            Flower_GroupIterator *groupIterator;
            Group *group;
            groupIterator = flower_getGroupIterator(flower);
            while ((group = flower_getNextGroup(groupIterator)) != NULL) {
                if (group_isLeaf(group)) {
                    //assert(group_getTotalBaseLength(group) == 0);
                    group_makeNestedFlower(group);
                    flower_setBuiltBlocks(group_getNestedFlower(group), 1);
                }
            }
            flower_destructGroupIterator(groupIterator);
        }
    }

    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the cactus disk\n");

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    st_logInfo("Am finished\n");
    return 0;
}
