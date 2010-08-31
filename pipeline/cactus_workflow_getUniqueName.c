#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cactus.h"

int main(int argc, char *argv[]) {
    /*
     * This code simply creates a unique name and puts it in a file.
     */
    assert(argc == 3);
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[1]);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    Name name = cactusDisk_getUniqueID(cactusDisk);
    FILE *fileHandle = fopen(argv[2], "w");
    fprintf(fileHandle, "%s\n", cactusMisc_nameToStringStatic(name));
    fclose(fileHandle);
    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    return 0;
}
