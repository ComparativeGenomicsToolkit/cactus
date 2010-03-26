#include "cactusGlobalsPrivate.h"

static NetDisk *netDisk = NULL;

void cactusNetMiscTestTeardown() {
	if(netDisk != NULL) {
		netDisk_destruct(netDisk);
		testCommon_deleteTemporaryNetDisk();
		netDisk = NULL;
	}
}

static void cactusNetMiscTestSetup() {
	cactusNetMiscTestTeardown();
	netDisk = netDisk_construct(testCommon_getTemporaryNetDisk());
}

void testNetMisc_reverseComplementChar(CuTest* testCase) {
	CuAssertTrue(testCase, netMisc_reverseComplementChar('A') == 'T');
	CuAssertTrue(testCase, netMisc_reverseComplementChar('T') == 'A');
	CuAssertTrue(testCase, netMisc_reverseComplementChar('G') == 'C');
	CuAssertTrue(testCase, netMisc_reverseComplementChar('C') == 'G');
	CuAssertTrue(testCase, netMisc_reverseComplementChar('a') == 't');
	CuAssertTrue(testCase, netMisc_reverseComplementChar('t') == 'a');
	CuAssertTrue(testCase, netMisc_reverseComplementChar('g') == 'c');
	CuAssertTrue(testCase, netMisc_reverseComplementChar('c') == 'g');
	CuAssertTrue(testCase, netMisc_reverseComplementChar('N') == 'N');
	CuAssertTrue(testCase, netMisc_reverseComplementChar('X') == 'X');
}

void testNetMisc_reverseComplementString(CuTest* testCase) {
	char *cA = netMisc_reverseComplementString("ACTG");
	CuAssertStrEquals(testCase, cA, "CAGT");
	free(cA);
	cA = netMisc_reverseComplementString("");
	CuAssertStrEquals(testCase, cA, "");
	free(cA);
}

void testNetMisc_nameCompare(CuTest* testCase) {
	cactusNetMiscTestSetup();
	Name name = netDisk_getUniqueID(netDisk);
	Name name2 = netDisk_getUniqueID(netDisk);
	CuAssertTrue(testCase, netMisc_nameCompare(name, name2) == -1);
	CuAssertTrue(testCase, netMisc_nameCompare(name2, name) == 1);
	CuAssertTrue(testCase, netMisc_nameCompare(name, name) == 0);
	cactusNetMiscTestTeardown();
}

void testNetMisc_stringNameFns(CuTest* testCase) {
	cactusNetMiscTestSetup();
	int32_t i;
	for(i=0; i<100000; i++) {
		Name name = netDisk_getUniqueID(netDisk);
		CuAssertTrue(testCase, netMisc_nameCompare(netMisc_stringToName(netMisc_nameToStringStatic(name)), name) == 0);
		char *cA = netMisc_nameToString(name);
		CuAssertStrEquals(testCase, cA, netMisc_nameToStringStatic(name));
		free(cA);
	}
	cactusNetMiscTestTeardown();
}

void testNetMisc_stringNameWithOrientationFns(CuTest* testCase) {
	cactusNetMiscTestSetup();
	int32_t i;
	for(i=0; i<100000; i++) {
		Name name = netDisk_getUniqueID(netDisk);
		char *cA = netMisc_nameToStringWithOrientation(name, 0);
		char *cA2 = netMisc_nameToStringWithOrientation(name, 1);
		CuAssertStrEquals(testCase, cA, netMisc_nameToStringStaticWithOrientiation(name, 0));
		CuAssertStrEquals(testCase, cA2, netMisc_nameToStringStaticWithOrientiation(name, 1));
		CuAssertStrEquals(testCase, netMisc_nameToString(name), netMisc_nameToStringStaticWithOrientiation(name, 1));
		CuAssertStrEquals(testCase, netMisc_nameToString(name), &(netMisc_nameToStringStaticWithOrientiation(name, 0)[1]));
		CuAssertTrue(testCase, '-' == cA[0]);
		free(cA);
		free(cA2);
	}
	cactusNetMiscTestTeardown();
}

CuSuite* cactusNetMiscTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testNetMisc_reverseComplementChar);
	SUITE_ADD_TEST(suite, testNetMisc_reverseComplementString);
	SUITE_ADD_TEST(suite, testNetMisc_nameCompare);
	SUITE_ADD_TEST(suite, testNetMisc_stringNameFns);
	SUITE_ADD_TEST(suite, testNetMisc_stringNameWithOrientationFns);
	return suite;
}
