/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk = NULL;

static void cactusMiscTestTeardown() {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(cactusDisk);
        cactusDisk = NULL;
    }
}

static void cactusMiscTestSetup() {
    cactusMiscTestTeardown();
    cactusDisk = testCommon_getTemporaryCactusDisk();
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

    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('w') == 'w');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('W') == 'W');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('s') == 's');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('S') == 'S');

    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('m') == 'k');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('M') == 'K');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('k') == 'm');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('K') == 'M');

    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('r') == 'y');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('R') == 'Y');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('y') == 'r');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('Y') == 'R');

    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('b') == 'v');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('B') == 'V');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('v') == 'b');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('V') == 'B');

    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('d') == 'h');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('D') == 'H');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('h') == 'd');
    CuAssertTrue(testCase, cactusMisc_reverseComplementChar('H') == 'D');

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
    cactusMiscTestTeardown();
}

static void testCactusCheck(CuTest* testCase) {
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
    SUITE_ADD_TEST(suite, testCactusMisc_reverseComplementChar);
    SUITE_ADD_TEST(suite, testCactusMisc_reverseComplementString);
    SUITE_ADD_TEST(suite, testCactusMisc_nameCompare);
    SUITE_ADD_TEST(suite, testCactusMisc_stringNameFns);
    SUITE_ADD_TEST(suite, testCactusCheck);
    return suite;
}
