#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk = NULL;
static const char *cactusDiskFile;

static void cactusDiskTestTeardown() {
    if (cactusDisk != NULL) {
        cactusDisk_destruct(cactusDisk);
        testCommon_deleteTemporaryCactusDisk(cactusDiskFile);
        cactusDisk = NULL;
    }
}

static void cactusDiskTestSetup() {
    cactusDiskTestTeardown();
    st_setLogLevel(ST_LOGGING_DEBUG);
    cactusDiskFile = testCommon_getTemporaryCactusDisk();
    cactusDisk = cactusDisk_construct(cactusDiskFile);
}

void testCactusDisk_constructAndDestruct(CuTest* testCase) {
    cactusDiskTestSetup();
    CuAssertTrue(testCase, cactusDisk != NULL); //check the flower is actually constructed.
    cactusDiskTestTeardown();
}

void testCactusDisk_write(CuTest* testCase) {
    assert(testCase != NULL);
    cactusDiskTestSetup();
    flower_construct(cactusDisk);
    cactusDisk_write(cactusDisk);
    cactusDiskTestTeardown();
}

void testCactusDisk_getFlower(CuTest* testCase) {
    cactusDiskTestSetup();
    Flower *flower = flower_construct(cactusDisk);
    Flower *flower2 = flower_construct(cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getFlower(cactusDisk, flower_getName(flower)) == flower);
    CuAssertTrue(testCase, cactusDisk_getFlower(cactusDisk, flower_getName(flower2)) == flower2);
    //now try closing the disk, then reloading it, to see if we get the same result.
    Name name1 = flower_getName(flower);
    Name name2 = flower_getName(flower2);
    cactusDisk_write(cactusDisk);
    cactusDisk_destruct(cactusDisk);
    cactusDisk = cactusDisk_construct(cactusDiskFile);
    flower = cactusDisk_getFlower(cactusDisk, name1);
    flower2 = cactusDisk_getFlower(cactusDisk, name2);
    CuAssertTrue(testCase, flower != NULL);
    CuAssertTrue(testCase, flower2 != NULL);
    CuAssertTrue(testCase, flower_getName(flower) == name1);
    CuAssertTrue(testCase, flower_getName(flower2) == name2);
    cactusDiskTestTeardown();
}

void testCactusDisk_getMetaSequence(CuTest* testCase) {
    cactusDiskTestSetup();
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
    cactusDisk = cactusDisk_construct(cactusDiskFile);
    metaSequence = cactusDisk_getMetaSequence(cactusDisk, name1);
    metaSequence2 = cactusDisk_getMetaSequence(cactusDisk, name2);
    CuAssertTrue(testCase, metaSequence != NULL);
    CuAssertTrue(testCase, metaSequence2 != NULL);
    CuAssertTrue(testCase, metaSequence_getName(metaSequence) == name1);
    CuAssertTrue(testCase, metaSequence_getName(metaSequence2) == name2);
    cactusDiskTestTeardown();
}

void testCactusDisk_getUniqueID(CuTest* testCase) {
    cactusDiskTestSetup();
    for (int32_t i = 0; i < 1000000000; i++) { //Gets a billion ids, checks we are good.
        Name uniqueName = cactusDisk_getUniqueID(cactusDisk);
        CuAssertTrue(testCase, uniqueName > 0);
        CuAssertTrue(testCase, uniqueName < INT64_MAX);
        CuAssertTrue(testCase, uniqueName != NULL_NAME);
    }
    cactusDiskTestTeardown();
}

int testCactusDisk_getUniqueID_UniqueP(const void *a, const void *b) {
    return cactusMisc_nameCompare(cactusMisc_stringToName(a), cactusMisc_stringToName(b));
}

void testCactusDisk_getUniqueID_Unique(CuTest* testCase) {
    cactusDiskTestSetup();
    stSortedSet *uniqueNames = stSortedSet_construct3(testCactusDisk_getUniqueID_UniqueP, free);
    for (int32_t i = 0; i < 1000000; i++) { //Gets a billion ids, checks we are good.
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
    cactusDiskTestTeardown();
}

CuSuite* cactusDiskTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testCactusDisk_constructAndDestruct);
    SUITE_ADD_TEST(suite, testCactusDisk_write);
    SUITE_ADD_TEST(suite, testCactusDisk_getFlower);
    SUITE_ADD_TEST(suite, testCactusDisk_getMetaSequence);
    SUITE_ADD_TEST(suite, testCactusDisk_getUniqueID);
    SUITE_ADD_TEST(suite, testCactusDisk_getUniqueID_Unique);
    return suite;
}
