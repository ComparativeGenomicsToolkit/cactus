/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
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
    stKVDatabaseConf *conf = stKVDatabaseConf_constructTokyoCabinet("temporaryCactusDisk");
    return conf;
}

void testCommon_deleteTemporaryKVDatabase() {
	int32_t i = system("rm -rf temporaryCactusDisk");
	exitOnFailure(i, "Tried to delete the temporary KV database\n");
}

CactusDisk *testCommon_getTemporaryCactusDisk() {
    stKVDatabaseConf *conf = stKVDatabaseConf_constructTokyoCabinet("temporaryCactusDisk");
    CactusDisk *cactusDisk = cactusDisk_construct(conf, 1);
    return cactusDisk;
}

void testCommon_deleteTemporaryCactusDisk(CactusDisk *cactusDisk) {
   cactusDisk_destruct(cactusDisk);
   testCommon_deleteTemporaryKVDatabase();
}
