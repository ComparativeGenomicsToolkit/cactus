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

void testCactusDisk_getFlowerNumberOnDisk(CuTest* testCase) {
    cactusDiskTestSetup();
    CuAssertIntEquals(testCase, 0, cactusDisk_getFlowerNumberOnDisk(cactusDisk));
    flower_construct(cactusDisk);
    CuAssertIntEquals(testCase, 0, cactusDisk_getFlowerNumberOnDisk(cactusDisk));
    flower_construct(cactusDisk);
    CuAssertIntEquals(testCase, 0, cactusDisk_getFlowerNumberOnDisk(cactusDisk));
    cactusDisk_write(cactusDisk);
    CuAssertIntEquals(testCase, 2, cactusDisk_getFlowerNumberOnDisk(cactusDisk));
    cactusDisk_destruct(cactusDisk);
    cactusDisk = cactusDisk_construct(cactusDiskFile);
    CuAssertIntEquals(testCase, 2, cactusDisk_getFlowerNumberOnDisk(cactusDisk));
    cactusDiskTestTeardown();
}

void testCactusDisk_flowerNamesOnDiskIterator(CuTest* testCase) {
    cactusDiskTestSetup();
    CactusDisk_FlowerNameIterator *iterator = cactusDisk_getFlowerNamesOnDiskIterator(
            cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getNextFlowerName(iterator) == NULL_NAME);
    CuAssertTrue(testCase, cactusDisk_getNextFlowerName(iterator) == NULL_NAME);
    cactusDisk_destructFlowerNamesOnDiskIterator(iterator);
    cactusDiskTestTeardown();
}

void testCactusDisk_getNextFlowerName(CuTest* testCase) {
    cactusDiskTestSetup();
    Name name1 = flower_getName(flower_construct(cactusDisk));
    Name name2 = flower_getName(flower_construct(cactusDisk));
    cactusDisk_write(cactusDisk);
    cactusDisk_destruct(cactusDisk);
    cactusDisk = cactusDisk_construct(cactusDiskFile);
    Name name3 = flower_getName(flower_construct(cactusDisk));
    Name name4 = flower_getName(flower_construct(cactusDisk));
    cactusDisk_write(cactusDisk);
    CactusDisk_FlowerNameIterator *iterator = cactusDisk_getFlowerNamesOnDiskIterator(
            cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getNextFlowerName(iterator) == name1);
    CuAssertTrue(testCase, cactusDisk_getNextFlowerName(iterator) == name2);
    CuAssertTrue(testCase, cactusDisk_getNextFlowerName(iterator) == name3);
    CuAssertTrue(testCase, cactusDisk_getNextFlowerName(iterator) == name4);
    CuAssertTrue(testCase, cactusDisk_getNextFlowerName(iterator) == NULL_NAME);
    cactusDisk_destructFlowerNamesOnDiskIterator(iterator);
    cactusDiskTestTeardown();
}

void testCactusDisk_getFlowerNumberInMemory(CuTest* testCase) {
    cactusDiskTestSetup();
    CuAssertIntEquals(testCase, 0, cactusDisk_getFlowerNumberInMemory(cactusDisk));
    flower_construct(cactusDisk);
    CuAssertIntEquals(testCase, 1, cactusDisk_getFlowerNumberInMemory(cactusDisk));
    flower_construct(cactusDisk);
    CuAssertIntEquals(testCase, 2, cactusDisk_getFlowerNumberInMemory(cactusDisk));
    cactusDisk_write(cactusDisk);
    CuAssertIntEquals(testCase, 2, cactusDisk_getFlowerNumberInMemory(cactusDisk));
    cactusDisk_destruct(cactusDisk);
    cactusDisk = cactusDisk_construct(cactusDiskFile);
    CuAssertIntEquals(testCase, 0, cactusDisk_getFlowerNumberInMemory(cactusDisk));
    cactusDiskTestTeardown();
}

void testCactusDisk_flowersInMemoryIterator(CuTest* testCase) {
    cactusDiskTestSetup();
    CactusDisk_FlowerIterator *iterator = cactusDisk_getFlowersInMemoryIterator(cactusDisk);
    CuAssertTrue(testCase, iterator != NULL);
    CuAssertTrue(testCase, cactusDisk_getNextFlower(iterator) == NULL);
    cactusDisk_destructFlowersInMemoryIterator(iterator);
    cactusDiskTestTeardown();
}

void testCactusDisk_getNextAndPreviousFlower(CuTest* testCase) {
    cactusDiskTestSetup();
    Flower *flower = flower_construct(cactusDisk);
    Flower *flower2 = flower_construct(cactusDisk);
    CactusDisk_FlowerIterator *iterator = cactusDisk_getFlowersInMemoryIterator(cactusDisk);
    CuAssertTrue(testCase, iterator != NULL);
    CuAssertTrue(testCase, cactusDisk_getNextFlower(iterator) == flower);
    CuAssertTrue(testCase, cactusDisk_getNextFlower(iterator) == flower2);
    CuAssertTrue(testCase, cactusDisk_getNextFlower(iterator) == NULL);
    CuAssertTrue(testCase, cactusDisk_getPreviousFlower(iterator) == flower2);
    CuAssertTrue(testCase, cactusDisk_getPreviousFlower(iterator) == flower);
    CuAssertTrue(testCase, cactusDisk_getPreviousFlower(iterator) == NULL);
    cactusDisk_destructFlowersInMemoryIterator(iterator);
    cactusDiskTestTeardown();
}

void testCactusDisk_copyFlowerIterator(CuTest* testCase) {
    cactusDiskTestSetup();
    Flower *flower = flower_construct(cactusDisk);
    Flower *flower2 = flower_construct(cactusDisk);
    CactusDisk_FlowerIterator *iterator = cactusDisk_getFlowersInMemoryIterator(cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getNextFlower(iterator) == flower);
    CactusDisk_FlowerIterator *iterator2 = cactusDisk_copyFlowerIterator(iterator);
    CuAssertTrue(testCase, cactusDisk_getNextFlower(iterator) == flower2);
    CuAssertTrue(testCase, cactusDisk_getNextFlower(iterator) == NULL);
    CuAssertTrue(testCase, cactusDisk_getNextFlower(iterator2) == flower2);
    CuAssertTrue(testCase, cactusDisk_getNextFlower(iterator2) == NULL);
    cactusDisk_destructFlowersInMemoryIterator(iterator);
    cactusDisk_destructFlowersInMemoryIterator(iterator2);
    cactusDiskTestTeardown();
}

CuSuite* cactusDiskTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testCactusDisk_constructAndDestruct);
    SUITE_ADD_TEST(suite, testCactusDisk_write);
    SUITE_ADD_TEST(suite, testCactusDisk_getFlower);
    SUITE_ADD_TEST(suite, testCactusDisk_getFlowerNumberOnDisk);
    SUITE_ADD_TEST(suite, testCactusDisk_flowerNamesOnDiskIterator);
    SUITE_ADD_TEST(suite, testCactusDisk_getNextFlowerName);
    SUITE_ADD_TEST(suite, testCactusDisk_flowersInMemoryIterator);
    SUITE_ADD_TEST(suite, testCactusDisk_getNextAndPreviousFlower);
    SUITE_ADD_TEST(suite, testCactusDisk_copyFlowerIterator);
    return suite;
}
