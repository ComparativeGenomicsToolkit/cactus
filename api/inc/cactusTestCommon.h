/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_TEST_COMMON_H_
#define CACTUS_TEST_COMMON_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions shared by the test code.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Gets a temporary KV database conf file suitable for building a database in.
 */
stKVDatabaseConf *testCommon_getTemporaryKVDatabaseConf();

/*
 * Removes from disk and temporary KV database created in the location pointed
 * to by the testCommon_getTemporaryKVDatabaseConf file.
 */
void testCommon_deleteTemporaryKVDatabase();

/*
 * Gets a temporary cactus disk, must call testCommon_deleteTemporaryCactusDisk
 * to destroy it.
 */
CactusDisk *testCommon_getTemporaryCactusDisk();

CactusDisk *testCommon_getTemporaryCactusDisk2(bool sequencesOnDisk);

/*
 * Destroys the object given by testCommon_getTemporaryCactusDisk and removes
 * all files associated from disk.
 */
void testCommon_deleteTemporaryCactusDisk(CactusDisk *cactusDisk);

// Adds a thread with random nucleotides to the flower, and return its corresponding name in the pinch graph.
Name testCommon_addThreadToFlower(Flower *flower, char *header, int64_t length);

#endif
