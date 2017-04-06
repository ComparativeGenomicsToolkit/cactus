#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include "cactus.h"

int64_t minFlowerSize;
FlowerWriter *flowerWriter;
CactusDisk *cactusDisk;

void parseArgs(int argc, char *argv[]) {
    /*
     * This code iterates through the terminal groups and returns
     * a list of the new flowers.
     */
    assert(argc >= 6);

    st_setLogLevelFromString(argv[1]);
    st_logDebug("Set up logging\n");

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[2]);
    // We disable caching for these utilities--there isn't a lot to
    // gain, and the memory advantage of "streaming" through input
    // flowers is nullified if we keep the information around in the
    // cactusDisk cache.
    cactusDisk = cactusDisk_construct2(kvDatabaseConf, false, argv[3], false);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    st_logDebug("Set up the flower disk\n");

    int64_t i = sscanf(argv[4], "%" PRId64 "", &minFlowerSize);
    assert(i == 1); //We should parse one in correctly.
    st_logDebug("Min flower size %lli\n", minFlowerSize);
    assert(minFlowerSize >= 0);

    int64_t maxSequenceSizeOfFlowerGrouping;
    i = sscanf(argv[5], "%" PRId64 "", &maxSequenceSizeOfFlowerGrouping);
    assert(i == 1); //We should parse one in correctly.
    if(maxSequenceSizeOfFlowerGrouping == -1) {
        maxSequenceSizeOfFlowerGrouping = INT64_MAX;
    }
    st_logDebug("Max size of flower grouping %lli\n", maxSequenceSizeOfFlowerGrouping);
    assert(maxSequenceSizeOfFlowerGrouping >= 0);

    int64_t maxSequenceSizeOfSecondaryFlowerGrouping;
    i = sscanf(argv[6], "%" PRId64 "", &maxSequenceSizeOfSecondaryFlowerGrouping);
    assert(i == 1); //We should parse one in correctly.
    if(maxSequenceSizeOfSecondaryFlowerGrouping == -1) {
        maxSequenceSizeOfSecondaryFlowerGrouping = INT64_MAX;
    }
    st_logDebug("Max size of flower secondary grouping %lli\n", maxSequenceSizeOfSecondaryFlowerGrouping);
    assert(maxSequenceSizeOfSecondaryFlowerGrouping >= 0);

    flowerWriter = flowerWriter_construct(stdout, maxSequenceSizeOfFlowerGrouping, maxSequenceSizeOfSecondaryFlowerGrouping);
}
