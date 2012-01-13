int64_t minFlowerSize;
int64_t maxFlowerSize;
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
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    st_logDebug("Set up the flower disk\n");

    int32_t i = sscanf(argv[3], "%" PRId64 "", &minFlowerSize);
    assert(i == 1); //We should parse one in correctly.
    st_logDebug("Min flower size %lli\n", minFlowerSize);
    assert(minFlowerSize >= 0);

    i = sscanf(argv[4], "%" PRId64 "", &maxFlowerSize);
    assert(i == 1); //We should parse one in correctly.
    if(maxFlowerSize == -1) {
        maxFlowerSize = INT64_MAX;
    }
    st_logDebug("Max flower size %lli\n", maxFlowerSize);
    assert(maxFlowerSize >= 0);

    int64_t maxSequenceSizeOfFlowerGrouping;
    i = sscanf(argv[5], "%" PRId64 "", &maxSequenceSizeOfFlowerGrouping);
    assert(i == 1); //We should parse one in correctly.
    if(maxSequenceSizeOfFlowerGrouping == -1) {
        maxSequenceSizeOfFlowerGrouping = INT64_MAX;
    }
    st_logDebug("Max size of flower grouping %lli\n", maxSequenceSizeOfFlowerGrouping);
    assert(maxSequenceSizeOfFlowerGrouping >= 0);

    flowerWriter = flowerWriter_construct(stdout, maxSequenceSizeOfFlowerGrouping);
}
