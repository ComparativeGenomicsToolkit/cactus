#include "cactusGlobalsPrivate.h"

static NetDisk *netDisk = NULL;
static const char *netDiskFile;

void cactusNetDiskTestTeardown() {
	if(netDisk != NULL) {
		netDisk_destruct(netDisk);
		testCommon_deleteTemporaryNetDisk(netDiskFile);
		netDisk = NULL;
	}
}

void cactusNetDiskTestSetup() {
	cactusNetDiskTestTeardown();
	st_setLogLevel(ST_LOGGING_DEBUG);
	netDiskFile = testCommon_getTemporaryNetDisk();
	netDisk = netDisk_construct(netDiskFile);
}

void testNetDisk_constructAndDestruct(CuTest* testCase) {
	cactusNetDiskTestSetup();
	CuAssertTrue(testCase, netDisk != NULL); //check the net is actually constructed.
	cactusNetDiskTestTeardown();
}

void testNetDisk_write(CuTest* testCase) {
	assert(testCase != NULL);
	cactusNetDiskTestSetup();
	net_construct(netDisk);
	netDisk_write(netDisk);
	cactusNetDiskTestTeardown();
}

void testNetDisk_getNet(CuTest* testCase) {
	cactusNetDiskTestSetup();
	Net *net = net_construct(netDisk);
	Net *net2 = net_construct(netDisk);
	CuAssertTrue(testCase, netDisk_getNet(netDisk, net_getName(net)) == net);
	CuAssertTrue(testCase, netDisk_getNet(netDisk, net_getName(net2)) == net2);
	//now try closing the disk, then reloading it, to see if we get the same result.
	Name name1 = net_getName(net);
	Name name2 =  net_getName(net2);
	netDisk_write(netDisk);
	netDisk_destruct(netDisk);
	netDisk = netDisk_construct(netDiskFile);
	net = netDisk_getNet(netDisk, name1);
	net2 = netDisk_getNet(netDisk, name2);
	CuAssertTrue(testCase, net != NULL);
	CuAssertTrue(testCase, net2 != NULL);
	CuAssertTrue(testCase, net_getName(net) == name1);
	CuAssertTrue(testCase, net_getName(net2) == name2);
	cactusNetDiskTestTeardown();
}

void testNetDisk_getNetNumberOnDisk(CuTest* testCase) {
	cactusNetDiskTestSetup();
	CuAssertIntEquals(testCase, 0, netDisk_getNetNumberOnDisk(netDisk));
	net_construct(netDisk);
	CuAssertIntEquals(testCase, 0, netDisk_getNetNumberOnDisk(netDisk));
	net_construct(netDisk);
	CuAssertIntEquals(testCase, 0, netDisk_getNetNumberOnDisk(netDisk));
	netDisk_write(netDisk);
	CuAssertIntEquals(testCase, 2, netDisk_getNetNumberOnDisk(netDisk));
	netDisk_destruct(netDisk);
	netDisk = netDisk_construct(netDiskFile);
	CuAssertIntEquals(testCase, 2, netDisk_getNetNumberOnDisk(netDisk));
	cactusNetDiskTestTeardown();
}

void testNetDisk_netNamesOnDiskIterator(CuTest* testCase) {
	cactusNetDiskTestSetup();
	NetDisk_NetNameIterator *iterator = netDisk_getNetNamesOnDiskIterator(netDisk);
	CuAssertTrue(testCase, netDisk_getNextNetName(iterator) == NULL_NAME);
	CuAssertTrue(testCase, netDisk_getNextNetName(iterator) == NULL_NAME);
	netDisk_destructNetNamesOnDiskIterator(iterator);
	cactusNetDiskTestTeardown();
}

void testNetDisk_getNextNetName(CuTest* testCase) {
	cactusNetDiskTestSetup();
	Name name1 = net_getName(net_construct(netDisk));
	Name name2 = net_getName(net_construct(netDisk));
	netDisk_write(netDisk);
	netDisk_destruct(netDisk);
	netDisk = netDisk_construct(netDiskFile);
	Name name3 = net_getName(net_construct(netDisk));
	Name name4 = net_getName(net_construct(netDisk));
	netDisk_write(netDisk);
	NetDisk_NetNameIterator *iterator = netDisk_getNetNamesOnDiskIterator(netDisk);
	CuAssertTrue(testCase, netDisk_getNextNetName(iterator) == name1);
	CuAssertTrue(testCase, netDisk_getNextNetName(iterator) == name2);
	CuAssertTrue(testCase, netDisk_getNextNetName(iterator) == name3);
	CuAssertTrue(testCase, netDisk_getNextNetName(iterator) == name4);
	CuAssertTrue(testCase, netDisk_getNextNetName(iterator) == NULL_NAME);
	netDisk_destructNetNamesOnDiskIterator(iterator);
	cactusNetDiskTestTeardown();
}

void testNetDisk_getNetNumberInMemory(CuTest* testCase) {
	cactusNetDiskTestSetup();
	CuAssertIntEquals(testCase, 0, netDisk_getNetNumberInMemory(netDisk));
	net_construct(netDisk);
	CuAssertIntEquals(testCase, 1, netDisk_getNetNumberInMemory(netDisk));
	net_construct(netDisk);
	CuAssertIntEquals(testCase, 2, netDisk_getNetNumberInMemory(netDisk));
	netDisk_write(netDisk);
	CuAssertIntEquals(testCase, 2, netDisk_getNetNumberInMemory(netDisk));
	netDisk_destruct(netDisk);
	netDisk = netDisk_construct(netDiskFile);
	CuAssertIntEquals(testCase, 0, netDisk_getNetNumberInMemory(netDisk));
	cactusNetDiskTestTeardown();
}

void testNetDisk_netsInMemoryIterator(CuTest* testCase) {
	cactusNetDiskTestSetup();
	NetDisk_NetIterator *iterator = netDisk_getNetsInMemoryIterator(netDisk);
	CuAssertTrue(testCase, iterator != NULL);
	CuAssertTrue(testCase, netDisk_getNextNet(iterator) == NULL);
	netDisk_destructNetsInMemoryIterator(iterator);
	cactusNetDiskTestTeardown();
}

void testNetDisk_getNextAndPreviousNet(CuTest* testCase) {
	cactusNetDiskTestSetup();
	Net *net = net_construct(netDisk);
	Net *net2 = net_construct(netDisk);
	NetDisk_NetIterator *iterator = netDisk_getNetsInMemoryIterator(netDisk);
	CuAssertTrue(testCase, iterator != NULL);
	CuAssertTrue(testCase, netDisk_getNextNet(iterator) == net);
	CuAssertTrue(testCase, netDisk_getNextNet(iterator) == net2);
	CuAssertTrue(testCase, netDisk_getNextNet(iterator) == NULL);
	CuAssertTrue(testCase, netDisk_getPreviousNet(iterator) == net2);
	CuAssertTrue(testCase, netDisk_getPreviousNet(iterator) == net);
	CuAssertTrue(testCase, netDisk_getPreviousNet(iterator) == NULL);
	netDisk_destructNetsInMemoryIterator(iterator);
	cactusNetDiskTestTeardown();
}

void testNetDisk_copyNetIterator(CuTest* testCase) {
	cactusNetDiskTestSetup();
	Net *net = net_construct(netDisk);
	Net *net2 = net_construct(netDisk);
	NetDisk_NetIterator *iterator = netDisk_getNetsInMemoryIterator(netDisk);
	CuAssertTrue(testCase, netDisk_getNextNet(iterator) == net);
	NetDisk_NetIterator *iterator2 = netDisk_copyNetIterator(iterator);
	CuAssertTrue(testCase, netDisk_getNextNet(iterator) == net2);
	CuAssertTrue(testCase, netDisk_getNextNet(iterator) == NULL);
	CuAssertTrue(testCase, netDisk_getNextNet(iterator2) == net2);
	CuAssertTrue(testCase, netDisk_getNextNet(iterator2) == NULL);
	netDisk_destructNetsInMemoryIterator(iterator);
	netDisk_destructNetsInMemoryIterator(iterator2);
	cactusNetDiskTestTeardown();
}

CuSuite* cactusNetDiskTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testNetDisk_constructAndDestruct);
	SUITE_ADD_TEST(suite, testNetDisk_write);
	SUITE_ADD_TEST(suite, testNetDisk_getNet);
	SUITE_ADD_TEST(suite, testNetDisk_getNetNumberOnDisk);
	SUITE_ADD_TEST(suite, testNetDisk_netNamesOnDiskIterator);
	SUITE_ADD_TEST(suite, testNetDisk_getNextNetName);
	SUITE_ADD_TEST(suite, testNetDisk_netsInMemoryIterator);
	SUITE_ADD_TEST(suite, testNetDisk_getNextAndPreviousNet);
	SUITE_ADD_TEST(suite, testNetDisk_copyNetIterator);
	return suite;
}
