#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk = NULL;
static const char *cactusDiskFile;

static void cactusDiskTestTeardown() {
    if (cactusDisk != NULL) {
        cactusDisk_destruct(cactusDisk);
        testCommon_deleteTemporaryNetDisk(cactusDiskFile);
        cactusDisk = NULL;
    }
}

static void cactusDiskTestSetup() {
    cactusDiskTestTeardown();
    st_setLogLevel(ST_LOGGING_DEBUG);
    cactusDiskFile = testCommon_getTemporaryNetDisk();
    cactusDisk = cactusDisk_construct(cactusDiskFile);
}

void testCactusDisk_constructAndDestruct(CuTest* testCase) {
    cactusDiskTestSetup();
    CuAssertTrue(testCase, cactusDisk != NULL); //check the net is actually constructed.
    cactusDiskTestTeardown();
}

void testCactusDisk_write(CuTest* testCase) {
    assert(testCase != NULL);
    cactusDiskTestSetup();
    net_construct(cactusDisk);
    cactusDisk_write(cactusDisk);
    cactusDiskTestTeardown();
}

void testCactusDisk_getNet(CuTest* testCase) {
    cactusDiskTestSetup();
    Net *net = net_construct(cactusDisk);
    Net *net2 = net_construct(cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getNet(cactusDisk, net_getName(net)) == net);
    CuAssertTrue(testCase, cactusDisk_getNet(cactusDisk, net_getName(net2)) == net2);
    //now try closing the disk, then reloading it, to see if we get the same result.
    Name name1 = net_getName(net);
    Name name2 = net_getName(net2);
    cactusDisk_write(cactusDisk);
    cactusDisk_destruct(cactusDisk);
    cactusDisk = cactusDisk_construct(cactusDiskFile);
    net = cactusDisk_getNet(cactusDisk, name1);
    net2 = cactusDisk_getNet(cactusDisk, name2);
    CuAssertTrue(testCase, net != NULL);
    CuAssertTrue(testCase, net2 != NULL);
    CuAssertTrue(testCase, net_getName(net) == name1);
    CuAssertTrue(testCase, net_getName(net2) == name2);
    cactusDiskTestTeardown();
}

void testCactusDisk_getNetNumberOnDisk(CuTest* testCase) {
    cactusDiskTestSetup();
    CuAssertIntEquals(testCase, 0, cactusDisk_getNetNumberOnDisk(cactusDisk));
    net_construct(cactusDisk);
    CuAssertIntEquals(testCase, 0, cactusDisk_getNetNumberOnDisk(cactusDisk));
    net_construct(cactusDisk);
    CuAssertIntEquals(testCase, 0, cactusDisk_getNetNumberOnDisk(cactusDisk));
    cactusDisk_write(cactusDisk);
    CuAssertIntEquals(testCase, 2, cactusDisk_getNetNumberOnDisk(cactusDisk));
    cactusDisk_destruct(cactusDisk);
    cactusDisk = cactusDisk_construct(cactusDiskFile);
    CuAssertIntEquals(testCase, 2, cactusDisk_getNetNumberOnDisk(cactusDisk));
    cactusDiskTestTeardown();
}

void testCactusDisk_netNamesOnDiskIterator(CuTest* testCase) {
    cactusDiskTestSetup();
    CactusDisk_NetNameIterator *iterator = cactusDisk_getNetNamesOnDiskIterator(
            cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getNextNetName(iterator) == NULL_NAME);
    CuAssertTrue(testCase, cactusDisk_getNextNetName(iterator) == NULL_NAME);
    cactusDisk_destructNetNamesOnDiskIterator(iterator);
    cactusDiskTestTeardown();
}

void testCactusDisk_getNextNetName(CuTest* testCase) {
    cactusDiskTestSetup();
    Name name1 = net_getName(net_construct(cactusDisk));
    Name name2 = net_getName(net_construct(cactusDisk));
    cactusDisk_write(cactusDisk);
    cactusDisk_destruct(cactusDisk);
    cactusDisk = cactusDisk_construct(cactusDiskFile);
    Name name3 = net_getName(net_construct(cactusDisk));
    Name name4 = net_getName(net_construct(cactusDisk));
    cactusDisk_write(cactusDisk);
    CactusDisk_NetNameIterator *iterator = cactusDisk_getNetNamesOnDiskIterator(
            cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getNextNetName(iterator) == name1);
    CuAssertTrue(testCase, cactusDisk_getNextNetName(iterator) == name2);
    CuAssertTrue(testCase, cactusDisk_getNextNetName(iterator) == name3);
    CuAssertTrue(testCase, cactusDisk_getNextNetName(iterator) == name4);
    CuAssertTrue(testCase, cactusDisk_getNextNetName(iterator) == NULL_NAME);
    cactusDisk_destructNetNamesOnDiskIterator(iterator);
    cactusDiskTestTeardown();
}

void testCactusDisk_getNetNumberInMemory(CuTest* testCase) {
    cactusDiskTestSetup();
    CuAssertIntEquals(testCase, 0, cactusDisk_getNetNumberInMemory(cactusDisk));
    net_construct(cactusDisk);
    CuAssertIntEquals(testCase, 1, cactusDisk_getNetNumberInMemory(cactusDisk));
    net_construct(cactusDisk);
    CuAssertIntEquals(testCase, 2, cactusDisk_getNetNumberInMemory(cactusDisk));
    cactusDisk_write(cactusDisk);
    CuAssertIntEquals(testCase, 2, cactusDisk_getNetNumberInMemory(cactusDisk));
    cactusDisk_destruct(cactusDisk);
    cactusDisk = cactusDisk_construct(cactusDiskFile);
    CuAssertIntEquals(testCase, 0, cactusDisk_getNetNumberInMemory(cactusDisk));
    cactusDiskTestTeardown();
}

void testCactusDisk_netsInMemoryIterator(CuTest* testCase) {
    cactusDiskTestSetup();
    CactusDisk_NetIterator *iterator = cactusDisk_getNetsInMemoryIterator(cactusDisk);
    CuAssertTrue(testCase, iterator != NULL);
    CuAssertTrue(testCase, cactusDisk_getNextNet(iterator) == NULL);
    cactusDisk_destructNetsInMemoryIterator(iterator);
    cactusDiskTestTeardown();
}

void testCactusDisk_getNextAndPreviousNet(CuTest* testCase) {
    cactusDiskTestSetup();
    Net *net = net_construct(cactusDisk);
    Net *net2 = net_construct(cactusDisk);
    CactusDisk_NetIterator *iterator = cactusDisk_getNetsInMemoryIterator(cactusDisk);
    CuAssertTrue(testCase, iterator != NULL);
    CuAssertTrue(testCase, cactusDisk_getNextNet(iterator) == net);
    CuAssertTrue(testCase, cactusDisk_getNextNet(iterator) == net2);
    CuAssertTrue(testCase, cactusDisk_getNextNet(iterator) == NULL);
    CuAssertTrue(testCase, cactusDisk_getPreviousNet(iterator) == net2);
    CuAssertTrue(testCase, cactusDisk_getPreviousNet(iterator) == net);
    CuAssertTrue(testCase, cactusDisk_getPreviousNet(iterator) == NULL);
    cactusDisk_destructNetsInMemoryIterator(iterator);
    cactusDiskTestTeardown();
}

void testCactusDisk_copyNetIterator(CuTest* testCase) {
    cactusDiskTestSetup();
    Net *net = net_construct(cactusDisk);
    Net *net2 = net_construct(cactusDisk);
    CactusDisk_NetIterator *iterator = cactusDisk_getNetsInMemoryIterator(cactusDisk);
    CuAssertTrue(testCase, cactusDisk_getNextNet(iterator) == net);
    CactusDisk_NetIterator *iterator2 = cactusDisk_copyNetIterator(iterator);
    CuAssertTrue(testCase, cactusDisk_getNextNet(iterator) == net2);
    CuAssertTrue(testCase, cactusDisk_getNextNet(iterator) == NULL);
    CuAssertTrue(testCase, cactusDisk_getNextNet(iterator2) == net2);
    CuAssertTrue(testCase, cactusDisk_getNextNet(iterator2) == NULL);
    cactusDisk_destructNetsInMemoryIterator(iterator);
    cactusDisk_destructNetsInMemoryIterator(iterator2);
    cactusDiskTestTeardown();
}

CuSuite* cactusDiskTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testCactusDisk_constructAndDestruct);
    SUITE_ADD_TEST(suite, testCactusDisk_write);
    SUITE_ADD_TEST(suite, testCactusDisk_getNet);
    SUITE_ADD_TEST(suite, testCactusDisk_getNetNumberOnDisk);
    SUITE_ADD_TEST(suite, testCactusDisk_netNamesOnDiskIterator);
    SUITE_ADD_TEST(suite, testCactusDisk_getNextNetName);
    SUITE_ADD_TEST(suite, testCactusDisk_netsInMemoryIterator);
    SUITE_ADD_TEST(suite, testCactusDisk_getNextAndPreviousNet);
    SUITE_ADD_TEST(suite, testCactusDisk_copyNetIterator);
    return suite;
}
