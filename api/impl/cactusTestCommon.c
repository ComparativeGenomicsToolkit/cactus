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

static CactusDisk *globalCactusDisk;

stKVDatabaseConf *testCommon_getTemporaryKVDatabaseConf() {
    testCommon_deleteTemporaryKVDatabase();
    stKVDatabaseConf *conf = stKVDatabaseConf_constructTokyoCabinet(
            "temporaryCactusDisk");
    return conf;
}

void testCommon_deleteTemporaryKVDatabase() {
    if (globalCactusDisk != NULL) {
        cactusDisk_destruct(globalCactusDisk);
    }
    int64_t i = system("rm -rf temporaryCactusDisk");
    exitOnFailure(i, "Tried to delete the temporary KV database\n");
    globalCactusDisk = NULL;
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
    globalCactusDisk = cactusDisk;
    return cactusDisk;
}

CactusDisk *testCommon_getTemporaryCactusDisk() {
    return testCommon_getTemporaryCactusDisk2(st_random() > 0.5);
}

void testCommon_deleteTemporaryCactusDisk(CactusDisk *cactusDisk) {
    cactusDisk_destruct(cactusDisk);
    if (cactusDisk == globalCactusDisk) {
        globalCactusDisk = NULL;
    }
    testCommon_deleteTemporaryKVDatabase();
}
