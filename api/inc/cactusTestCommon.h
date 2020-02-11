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
 * Get a temporary directory for a test.
 */
char *testCommon_getTmpTestDir(const char *testName);

/*
 * Gets a temporary KV database conf file suitable for building a database in.
 */
stKVDatabaseConf *testCommon_getTemporaryKVDatabaseConf(const char *testName);

/*
 * Removes from disk and temporary KV database created in the location pointed
 * to by the testCommon_getTemporaryKVDatabaseConf file.
 */
void testCommon_deleteTemporaryKVDatabase(const char *testName);

/*
 * Gets a temporary cactus disk, must call testCommon_deleteTemporaryCactusDisk
 * to destroy it.
 */
CactusDisk *testCommon_getTemporaryCactusDisk(const char *testName);

/*
 * Destroys the object given by testCommon_getTemporaryCactusDisk and removes
 * all files associated from disk.
 */
void testCommon_deleteTemporaryCactusDisk(const char *testName,
                                          CactusDisk *cactusDisk);

// Adds a thread with random nucleotides to the flower, and return its corresponding name in the pinch graph.
Name testCommon_addThreadToFlower(Flower *flower, char *header, int64_t length);

#endif
