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

stKVDatabaseConf *testCommon_getTemporaryKVDatabaseConf() {
    testCommon_deleteTemporaryKVDatabase();
    stKVDatabaseConf *conf = stKVDatabaseConf_constructTokyoCabinet(
            "temporaryCactusDisk");
    return conf;
}

void testCommon_deleteTemporaryKVDatabase() {
    int64_t i = system("rm -rf temporaryCactusDisk");
    exitOnFailure(i, "Tried to delete the temporary KV database\n");
}

CactusDisk *testCommon_getTemporaryCactusDisk2(bool sequencesOnDisk) {
    testCommon_deleteTemporaryKVDatabase();
    stKVDatabaseConf *conf = stKVDatabaseConf_constructTokyoCabinet(
            "temporaryCactusDisk");
    CactusDisk *cactusDisk;
    if (sequencesOnDisk) {
        cactusDisk = cactusDisk_construct2(conf,
                        "cactusSequences");
    } else {
        cactusDisk = cactusDisk_construct(conf, 1);
    }
    stKVDatabaseConf_destruct(conf);
    return cactusDisk;
}

CactusDisk *testCommon_getTemporaryCactusDisk() {
    return testCommon_getTemporaryCactusDisk2(st_random() > 0.5);
}

void testCommon_deleteTemporaryCactusDisk(CactusDisk *cactusDisk) {
    cactusDisk_destruct(cactusDisk);
    testCommon_deleteTemporaryKVDatabase();
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
