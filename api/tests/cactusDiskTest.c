/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk = NULL;
static stKVDatabaseConf *conf = NULL;

static void cactusDiskTestTeardown(CuTest* testCase) {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
        stKVDatabaseConf_destruct(conf);
        cactusDisk = NULL;
    }
}

static void cactusDiskTestSetup(CuTest* testCase) {
    cactusDiskTestTeardown(testCase);
    conf = testCommon_getTemporaryKVDatabaseConf(testCase->name);
    cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
}

void testCactusDisk_constructAndDestruct(CuTest* testCase) {
    cactusDiskTestSetup(testCase);
    CuAssertTrue(testCase, cactusDisk != NULL); //check the flower is actually constructed.
    cactusDiskTestTeardown(testCase);
}

void testCactusDisk_write(CuTest* testCase) {
    assert(testCase != NULL);
    cactusDiskTestSetup(testCase);
    flower_construct(cactusDisk);
    cactusDisk_write(cactusDisk);
    cactusDiskTestTeardown(testCase);
}

void testCactusDisk_getFlower(CuTest* testCase) {
    cactusDiskTestSetup(testCase);
    Flower *flower = flower_construct(cactusDisk);
    Flower *flower2 = flower_construct(cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getFlower(cactusDisk, flower_getName(flower)) == flower);
    CuAssertTrue(testCase, cactusDisk_getFlower(cactusDisk, flower_getName(flower2)) == flower2);
    //now try closing the disk, then reloading it, to see if we get the same result.
    Name name1 = flower_getName(flower);
    Name name2 = flower_getName(flower2);
    cactusDisk_write(cactusDisk);
    cactusDisk_destruct(cactusDisk);
    cactusDisk = cactusDisk_construct(conf, false, true);
    flower = cactusDisk_getFlower(cactusDisk, name1);
    flower2 = cactusDisk_getFlower(cactusDisk, name2);
    CuAssertTrue(testCase, flower != NULL);
    CuAssertTrue(testCase, flower2 != NULL);
    CuAssertTrue(testCase, flower_getName(flower) == name1);
    CuAssertTrue(testCase, flower_getName(flower2) == name2);
    cactusDiskTestTeardown(testCase);
}

void testCactusDisk_getMetaSequence(CuTest* testCase) {
    cactusDiskTestSetup(testCase);
    MetaSequence *metaSequence = metaSequence_construct(1, 10, "ACTGACTGAG",
            "FOO", 10, cactusDisk);
    MetaSequence *metaSequence2 = metaSequence_construct(2, 10, "CCCCCCCCCC",
            "BAR", 10, cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getMetaSequence(cactusDisk, metaSequence_getName(metaSequence)) == metaSequence);
    CuAssertTrue(testCase, cactusDisk_getMetaSequence(cactusDisk, metaSequence_getName(metaSequence2)) == metaSequence2);
    //now try closing the disk, then reloading it, to see if we get the same result.
    Name name1 = metaSequence_getName(metaSequence);
    Name name2 = metaSequence_getName(metaSequence2);
    cactusDisk_write(cactusDisk);
    cactusDisk_destruct(cactusDisk);
    cactusDisk = cactusDisk_construct(conf, false, true);
    metaSequence = cactusDisk_getMetaSequence(cactusDisk, name1);
    metaSequence2 = cactusDisk_getMetaSequence(cactusDisk, name2);
    CuAssertTrue(testCase, metaSequence != NULL);
    CuAssertTrue(testCase, metaSequence2 != NULL);
    CuAssertTrue(testCase, metaSequence_getName(metaSequence) == name1);
    CuAssertTrue(testCase, metaSequence_getName(metaSequence2) == name2);
    cactusDiskTestTeardown(testCase);
}

void testCactusDisk_getUniqueID(CuTest* testCase) {
    cactusDiskTestSetup(testCase);
    for (int64_t i = 0; i < 1000000; i++) { //Gets a billion ids, checks we are good.
        Name uniqueName = cactusDisk_getUniqueID(cactusDisk);
        CuAssertTrue(testCase, uniqueName > 0);
        CuAssertTrue(testCase, uniqueName < INT64_MAX);
        CuAssertTrue(testCase, uniqueName != NULL_NAME);
    }
    cactusDiskTestTeardown(testCase);
}

int testCactusDisk_getUniqueID_UniqueP(const void *a, const void *b) {
    return cactusMisc_nameCompare(cactusMisc_stringToName(a), cactusMisc_stringToName(b));
}

void testCactusDisk_getUniqueID_Unique(CuTest* testCase) {
    cactusDiskTestSetup(testCase);
    stSortedSet *uniqueNames = stSortedSet_construct3(testCactusDisk_getUniqueID_UniqueP, free);
    for (int64_t i = 0; i < 100000; i++) { //Gets a billion ids, checks we are good.
        Name uniqueName = cactusDisk_getUniqueID(cactusDisk);
        CuAssertTrue(testCase, uniqueName > 0);
        CuAssertTrue(testCase, uniqueName < INT64_MAX);
        CuAssertTrue(testCase, uniqueName != NULL_NAME);
        char *cA = cactusMisc_nameToString(uniqueName);
        CuAssertTrue(testCase, stSortedSet_search(uniqueNames, cA) == NULL);
        CuAssertTrue(testCase, cactusMisc_stringToName(cA) == uniqueName);
        stSortedSet_insert(uniqueNames, cA);
    }
    stSortedSet_destruct(uniqueNames);
    cactusDiskTestTeardown(testCase);
}

void testCactusDisk_getUniqueID_UniqueIntervals(CuTest* testCase) {
    cactusDiskTestSetup(testCase);
    stSortedSet *uniqueNames = stSortedSet_construct3(testCactusDisk_getUniqueID_UniqueP, free);
    for (int64_t i = 0; i < 10; i++) { //Gets a billion ids, checks we are good.
        int64_t intervalSize = st_randomInt(0, 100000);
        Name uniqueName = cactusDisk_getUniqueIDInterval(cactusDisk, intervalSize);
        for(int64_t j=0; j<intervalSize; j++) {
            CuAssertTrue(testCase, uniqueName > 0);
            CuAssertTrue(testCase, uniqueName < INT64_MAX);
            CuAssertTrue(testCase, uniqueName != NULL_NAME);
            char *cA = cactusMisc_nameToString(uniqueName);
            CuAssertTrue(testCase, stSortedSet_search(uniqueNames, cA) == NULL);
            CuAssertTrue(testCase, cactusMisc_stringToName(cA) == uniqueName);
            stSortedSet_insert(uniqueNames, cA);
            uniqueName++;
        }
    }
    stSortedSet_destruct(uniqueNames);
    cactusDiskTestTeardown(testCase);
}

CuSuite* cactusDiskTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testCactusDisk_write);
    SUITE_ADD_TEST(suite, testCactusDisk_getFlower);
    SUITE_ADD_TEST(suite, testCactusDisk_getMetaSequence);
    SUITE_ADD_TEST(suite, testCactusDisk_getUniqueID);
    SUITE_ADD_TEST(suite, testCactusDisk_getUniqueID_Unique);
    SUITE_ADD_TEST(suite, testCactusDisk_getUniqueID_UniqueIntervals);
    SUITE_ADD_TEST(suite, testCactusDisk_constructAndDestruct);
    return suite;
}
