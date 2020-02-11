/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions shared by the test code.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

char *testCommon_getTmpTestDir(const char *testName) {
    return stFile_pathJoin("test-output/tmp", testName);
}


static char* getTestDatabasePath(const char *testName) {
    return stFile_pathJoin(testCommon_getTmpTestDir(testName), "cactusDisk");
}

stKVDatabaseConf *testCommon_getTemporaryKVDatabaseConf(const char *testName) {
    testCommon_deleteTemporaryKVDatabase(testName);
    char *dbPath = getTestDatabasePath(testName);
    stFile_mkdirp(dbPath);
    stKVDatabaseConf *conf = stKVDatabaseConf_constructTokyoCabinet(dbPath);
    return conf;
}

void testCommon_deleteTemporaryKVDatabase(const char *testName) {
    stFile_rmtree(getTestDatabasePath(testName));
}

CactusDisk *testCommon_getTemporaryCactusDisk(const char *testName) {
    testCommon_deleteTemporaryKVDatabase(testName);
    char *dbPath = getTestDatabasePath(testName);
    stFile_mkdirp(dbPath);
    stKVDatabaseConf *conf = stKVDatabaseConf_constructTokyoCabinet(dbPath);
    CactusDisk *cactusDisk;
    cactusDisk = cactusDisk_construct(conf, true, true);
    stKVDatabaseConf_destruct(conf);
    return cactusDisk;
}

void testCommon_deleteTemporaryCactusDisk(const char *testName,
                                          CactusDisk *cactusDisk) {
    cactusDisk_destruct(cactusDisk);
    testCommon_deleteTemporaryKVDatabase(testName);
}

Name testCommon_addThreadToFlower(Flower *flower, char *header, int64_t length) {
    char *dna = stRandom_getRandomDNAString(length, true, true, true);
    EventTree *eventTree = flower_getEventTree(flower);
    assert(eventTree != NULL);
    MetaSequence *metaSequence = metaSequence_construct(2, length, dna, header, event_getName(eventTree_getRootEvent(eventTree)), flower_getCactusDisk(flower));
    Sequence *sequence = sequence_construct(metaSequence, flower);

    End *end1 = end_construct2(0, 0, flower);
    End *end2 = end_construct2(1, 0, flower);
    Cap *cap1 = cap_construct2(end1, 1, 1, sequence);
    Cap *cap2 = cap_construct2(end2, length + 2, 1, sequence);
    cap_makeAdjacent(cap1, cap2);

    free(dna);
    return cap_getName(cap1);
}
