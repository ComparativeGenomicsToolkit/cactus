/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk = NULL;

static void cactusMiscTestTeardown(CuTest* testCase) {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
        cactusDisk = NULL;
    }
}

static void cactusMiscTestSetup(CuTest* testCase) {
    cactusMiscTestTeardown(testCase);
    cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
}

void testCactusMisc_nameCompare(CuTest* testCase) {
    cactusMiscTestSetup(testCase);
    Name name = cactusDisk_getUniqueID(cactusDisk);
    Name name2 = cactusDisk_getUniqueID(cactusDisk);
    CuAssertTrue(testCase, cactusMisc_nameCompare(name, name2) == -1);
    CuAssertTrue(testCase, cactusMisc_nameCompare(name2, name) == 1);
    CuAssertTrue(testCase, cactusMisc_nameCompare(name, name) == 0);
    cactusMiscTestTeardown(testCase);
}

void testCactusMisc_stringNameFns(CuTest* testCase) {
    cactusMiscTestSetup(testCase);
    int64_t i;
    for (i = 0; i < 1000000; i++) {
        Name name = cactusDisk_getUniqueID(cactusDisk);
        CuAssertTrue(
                testCase,
                cactusMisc_nameCompare(
                        cactusMisc_stringToName(
                                cactusMisc_nameToStringStatic(name)), name)
                        == 0);
        char *cA = cactusMisc_nameToString(name);
        CuAssertStrEquals(testCase, cA, cactusMisc_nameToStringStatic(name));
        free(cA);
    }
    cactusMiscTestTeardown(testCase);
}

static void testCactusCheck(CuTest* testCase) {
    return; //While we have an assert that fails in that function to provide a stack trace.
    cactusCheck(1);
    stTry {
        cactusCheck(0);
        CuAssertTrue(testCase, 0);
    } stCatch(except) {
        st_logInfo("This is the message %s\n", stExcept_getMsg(except));
        stExcept_free(except);
    } stTryEnd

    cactusCheck2(1, "This shouldn't throw an exception: %s", "blah");
    stTry {
        cactusCheck2(0, "This should throw an exception: %s", "blah");
        CuAssertTrue(testCase, 0);
    } stCatch(except) {
        st_logInfo("This is the message: %s\n", stExcept_getMsg(except));
        stExcept_free(except);
    } stTryEnd
}

CuSuite* cactusMiscTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testCactusMisc_nameCompare);
    SUITE_ADD_TEST(suite, testCactusMisc_stringNameFns);
    SUITE_ADD_TEST(suite, testCactusCheck);
    return suite;
}
