/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

void testCactusDisk_constructAndDestruct(CuTest* testCase) {
    CactusDisk *cactusDisk = cactusDisk_construct()
    CuAssertTrue(testCase, cactusDisk != NULL); //check the flower is actually constructed.
    cactusDisk_destruct(cactusDisk)
}

void testCactusDisk_getFlower(CuTest* testCase) {
    CactusDisk *cactusDisk = cactusDisk_construct()
    Flower *flower = flower_construct(cactusDisk);
    Flower *flower2 = flower_construct(cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getFlower(cactusDisk, flower_getName(flower)) == flower);
    CuAssertTrue(testCase, cactusDisk_getFlower(cactusDisk, flower_getName(flower2)) == flower2);
    cactusDisk_destruct(cactusDisk)
}

void testCactusDisk_getSequence(CuTest* testCase) {
    CactusDisk *cactusDisk = cactusDisk_construct()
    Sequence *sequence = sequence_construct(1, 10, "ACTGACTGAG",
            "FOO", NULL, cactusDisk);
    Sequence *sequence2 = sequence_construct(2, 10, "CCCCCCCCCC",
            "BAR", NULL, cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getSequence(cactusDisk, sequence_getName(sequence)) == sequence);
    CuAssertTrue(testCase, cactusDisk_getSequence(cactusDisk, sequence_getName(sequence2)) == sequence2);
    //now try closing the disk, then reloading it, to see if we get the same result.
    Name name1 = sequence_getName(sequence);
    Name name2 = sequence_getName(sequence2);
    cactusDisk_write(cactusDisk);
    cactusDisk_destruct(cactusDisk);
}

void testCactusDisk_getUniqueID(CuTest* testCase) {
    CactusDisk *cactusDisk = cactusDisk_construct()
    for (int64_t i = 0; i < 1000000; i++) { //Gets a billion ids, checks we are good.
        Name uniqueName = cactusDisk_getUniqueID(cactusDisk);
        CuAssertTrue(testCase, uniqueName > 0);
        CuAssertTrue(testCase, uniqueName < INT64_MAX);
        CuAssertTrue(testCase, uniqueName != NULL_NAME);
    }
    cactusDisk_destruct(cactusDisk)
}

int testCactusDisk_getUniqueID_UniqueP(const void *a, const void *b) {
    return cactusMisc_nameCompare(cactusMisc_stringToName(a), cactusMisc_stringToName(b));
}

void testCactusDisk_getUniqueID_Unique(CuTest* testCase) {
    CactusDisk *cactusDisk = cactusDisk_construct()
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
    cactusDisk_destruct(cactusDisk)
}

void testCactusDisk_getUniqueID_UniqueIntervals(CuTest* testCase) {
    CactusDisk *cactusDisk = cactusDisk_construct()
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
    cactusDisk_destruct(cactusDisk)
}

CuSuite* cactusDiskTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testCactusDisk_getFlower);
    SUITE_ADD_TEST(suite, testCactusDisk_getSequence);
    SUITE_ADD_TEST(suite, testCactusDisk_getUniqueID);
    SUITE_ADD_TEST(suite, testCactusDisk_getUniqueID_Unique);
    SUITE_ADD_TEST(suite, testCactusDisk_getUniqueID_UniqueIntervals);
    SUITE_ADD_TEST(suite, testCactusDisk_constructAndDestruct);
    return suite;
}
