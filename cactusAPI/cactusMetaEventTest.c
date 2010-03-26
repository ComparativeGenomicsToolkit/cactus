#include "cactusGlobalsPrivate.h"


static NetDisk *netDisk;
static MetaEvent *metaEvent;

static void cactusMetaEventTestTeardown() {
	if(netDisk != NULL) {
		netDisk_destruct(netDisk);
		testCommon_deleteTemporaryNetDisk();
		netDisk = NULL;
		metaEvent = NULL;
	}
}

static void cactusMetaEventTestSetup() {
	cactusMetaEventTestTeardown();
	netDisk = netDisk_construct(testCommon_getTemporaryNetDisk());
	metaEvent = metaEvent_construct("ROOT", netDisk);
}

void testMetaEvent_getName(CuTest* testCase) {
	cactusMetaEventTestSetup();
	CuAssertTrue(testCase, metaEvent_getName(metaEvent) != NULL_NAME);
	CuAssertTrue(testCase, netDisk_getMetaEvent(netDisk, metaEvent_getName(metaEvent)) == metaEvent);
	cactusMetaEventTestTeardown();
}

void testMetaEvent_getHeader(CuTest* testCase) {
	cactusMetaEventTestSetup();
	CuAssertStrEquals(testCase, "ROOT", metaEvent_getHeader(metaEvent));
	cactusMetaEventTestTeardown();
}

void testMetaEvent_serialisation(CuTest* testCase) {
	cactusMetaEventTestSetup();
	int32_t i;
	Name name = metaEvent_getName(metaEvent);
	CuAssertTrue(testCase, netDisk_getMetaEvent(netDisk, name) == metaEvent);
	void *vA = binaryRepresentation_makeBinaryRepresentation(metaEvent,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))metaEvent_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	metaEvent_destruct(metaEvent);
	CuAssertTrue(testCase, netDisk_getMetaEvent(netDisk, name) == NULL);
	void *vA2 = vA;
	metaEvent = metaEvent_loadFromBinaryRepresentation(&vA2, netDisk);
	CuAssertTrue(testCase, name == metaEvent_getName(metaEvent));
	CuAssertStrEquals(testCase, "ROOT", metaEvent_getHeader(metaEvent));
	CuAssertTrue(testCase, netDisk_getMetaEvent(netDisk, name) == metaEvent);
	netDisk_write(netDisk);
	metaEvent_destruct(metaEvent);
	CuAssertTrue(testCase, netDisk_getMetaEvent(netDisk, name) != NULL);
	cactusMetaEventTestTeardown();
}

CuSuite* cactusMetaEventTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testMetaEvent_getName);
	SUITE_ADD_TEST(suite, testMetaEvent_serialisation);
	SUITE_ADD_TEST(suite, testMetaEvent_getHeader);
	return suite;
}
