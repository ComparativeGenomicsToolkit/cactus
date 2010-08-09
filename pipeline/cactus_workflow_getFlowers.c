#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include "cactus.h"

int main(int argc, char *argv[]) {
    CactusDisk *cactusDisk;
    Flower *flower;

    assert(argc >= 3);
    cactusDisk = cactusDisk_construct(argv[1]);
    st_logInfo("Set up the flower disk\n");

    FILE *fileHandle = fopen(argv[2], "w");

    int32_t i;
    for (i = 3; i < argc; i++) {
        flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(argv[i]));
        assert(flower != NULL);
        assert(flower_builtBlocks(flower)); //This recursion depends on the block structure having been properly defined for all nodes.
        st_logInfo("Parsed the flower %s\n", argv[i]);
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if (!group_isLeaf(group)) {
                Flower *nestedFlower = group_getNestedFlower(group);
                assert(nestedFlower != NULL);
                assert(flower_builtBlocks(nestedFlower)); //This recursion depends on the block structure having been properly defined for all nodes.
                fprintf(fileHandle, "%s %" PRIi64 " \n", cactusMisc_nameToStringStatic(flower_getName(nestedFlower)), flower_getTotalBaseLength(nestedFlower));
            }
        }
        flower_destructGroupIterator(groupIterator);
    }
    fclose(fileHandle);
    cactusDisk_destruct(cactusDisk);
    st_logInfo("Am finished\n");
    return 0;
}
