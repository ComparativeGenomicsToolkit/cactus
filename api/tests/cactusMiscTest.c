#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk = NULL;

static void cactusMiscTestTeardown() {
	if(cactusDisk != NULL) {
		cactusDisk_destruct(cactusDisk);
		testCommon_deleteTemporaryCactusDisk();
		cactusDisk = NULL;
	}
}

static void cactusMiscTestSetup() {
	cactusMiscTestTeardown();
	cactusDisk = cactusDisk_construct(testCommon_getTemporaryCactusDisk());
}

void testCactusMisc_reverseComplementChar(CuTest* testCase) {
	CuAssertTrue(testCase, cactusMisc_reverseComplementChar('A') == 'T');
	CuAssertTrue(testCase, cactusMisc_reverseComplementChar('T') == 'A');
	CuAssertTrue(testCase, cactusMisc_reverseComplementChar('G') == 'C');
	CuAssertTrue(testCase, cactusMisc_reverseComplementChar('C') == 'G');
	CuAssertTrue(testCase, cactusMisc_reverseComplementChar('a') == 't');
	CuAssertTrue(testCase, cactusMisc_reverseComplementChar('t') == 'a');
	CuAssertTrue(testCase, cactusMisc_reverseComplementChar('g') == 'c');
	CuAssertTrue(testCase, cactusMisc_reverseComplementChar('c') == 'g');
	CuAssertTrue(testCase, cactusMisc_reverseComplementChar('N') == 'N');
	CuAssertTrue(testCase, cactusMisc_reverseComplementChar('X') == 'X');
}

void testCactusMisc_reverseComplementString(CuTest* testCase) {
	char *cA = cactusMisc_reverseComplementString("ACTG");
	CuAssertStrEquals(testCase, cA, "CAGT");
	free(cA);
	cA = cactusMisc_reverseComplementString("");
	CuAssertStrEquals(testCase, cA, "");
	free(cA);
}

void testCactusMisc_nameCompare(CuTest* testCase) {
	cactusMiscTestSetup();
	Name name = cactusDisk_getUniqueID(cactusDisk);
	Name name2 = cactusDisk_getUniqueID(cactusDisk);
	CuAssertTrue(testCase, cactusMisc_nameCompare(name, name2) == -1);
	CuAssertTrue(testCase, cactusMisc_nameCompare(name2, name) == 1);
	CuAssertTrue(testCase, cactusMisc_nameCompare(name, name) == 0);
	cactusMiscTestTeardown();
}

void testCactusMisc_stringNameFns(CuTest* testCase) {
	cactusMiscTestSetup();
	int32_t i;
	for(i=0; i<100000; i++) {
		Name name = cactusDisk_getUniqueID(cactusDisk);
		CuAssertTrue(testCase, cactusMisc_nameCompare(cactusMisc_stringToName(cactusMisc_nameToStringStatic(name)), name) == 0);
		char *cA = cactusMisc_nameToString(name);
		CuAssertStrEquals(testCase, cA, cactusMisc_nameToStringStatic(name));
		free(cA);
	}
	cactusMiscTestTeardown();
}

void testCactusMisc_stringNameWithOrientationFns(CuTest* testCase) {
	cactusMiscTestSetup();
	int32_t i;
	for(i=0; i<100000; i++) {
		Name name = cactusDisk_getUniqueID(cactusDisk);
		char *cA = cactusMisc_nameToStringWithOrientation(name, 0);
		char *cA2 = cactusMisc_nameToStringWithOrientation(name, 1);
		CuAssertStrEquals(testCase, cA, cactusMisc_nameToStringStaticWithOrientiation(name, 0));
		CuAssertStrEquals(testCase, cA2, cactusMisc_nameToStringStaticWithOrientiation(name, 1));
		CuAssertStrEquals(testCase, cactusMisc_nameToString(name), cactusMisc_nameToStringStaticWithOrientiation(name, 1));
		CuAssertStrEquals(testCase, cactusMisc_nameToString(name), &(cactusMisc_nameToStringStaticWithOrientiation(name, 0)[1]));
		CuAssertTrue(testCase, '-' == cA[0]);
		free(cA);
		free(cA2);
	}
	cactusMiscTestTeardown();
}

CuSuite* cactusMiscTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testCactusMisc_reverseComplementChar);
	SUITE_ADD_TEST(suite, testCactusMisc_reverseComplementString);
	SUITE_ADD_TEST(suite, testCactusMisc_nameCompare);
	SUITE_ADD_TEST(suite, testCactusMisc_stringNameFns);
	SUITE_ADD_TEST(suite, testCactusMisc_stringNameWithOrientationFns);
	return suite;
}
