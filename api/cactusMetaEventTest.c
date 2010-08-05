#include "cactusGlobalsPrivate.h"


static CactusDisk *cactusDisk;
static MetaEvent *metaEvent;

static void cactusMetaEventTestTeardown() {
	if(cactusDisk != NULL) {
		cactusDisk_destruct(cactusDisk);
		testCommon_deleteTemporaryCactusDisk();
		cactusDisk = NULL;
		metaEvent = NULL;
	}
}

static void cactusMetaEventTestSetup() {
	cactusMetaEventTestTeardown();
	cactusDisk = cactusDisk_construct(testCommon_getTemporaryCactusDisk());
	metaEvent = metaEvent_construct("ROOT", cactusDisk);
}

void testMetaEvent_getName(CuTest* testCase) {
	cactusMetaEventTestSetup();
	CuAssertTrue(testCase, metaEvent_getName(metaEvent) != NULL_NAME);
	CuAssertTrue(testCase, cactusDisk_getMetaEvent(cactusDisk, metaEvent_getName(metaEvent)) == metaEvent);
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
	CuAssertTrue(testCase, cactusDisk_getMetaEvent(cactusDisk, name) == metaEvent);
	void *vA = binaryRepresentation_makeBinaryRepresentation(metaEvent,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))metaEvent_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	metaEvent_destruct(metaEvent);
	CuAssertTrue(testCase, cactusDisk_getMetaEvent(cactusDisk, name) == NULL);
	void *vA2 = vA;
	metaEvent = metaEvent_loadFromBinaryRepresentation(&vA2, cactusDisk);
	CuAssertTrue(testCase, name == metaEvent_getName(metaEvent));
	CuAssertStrEquals(testCase, "ROOT", metaEvent_getHeader(metaEvent));
	CuAssertTrue(testCase, cactusDisk_getMetaEvent(cactusDisk, name) == metaEvent);
	cactusDisk_write(cactusDisk);
	metaEvent_destruct(metaEvent);
	CuAssertTrue(testCase, cactusDisk_getMetaEvent(cactusDisk, name) != NULL);
	cactusMetaEventTestTeardown();
}

CuSuite* cactusMetaEventTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testMetaEvent_getName);
	SUITE_ADD_TEST(suite, testMetaEvent_serialisation);
	SUITE_ADD_TEST(suite, testMetaEvent_getHeader);
	return suite;
}
