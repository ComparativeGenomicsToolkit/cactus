/*
 * stPinchGraphsTest.c
 *
 *  Created on: 11 Apr 2012
 *      Author: benedictpaten
 *
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"

static void teardown() {
}

static void setup() {
    teardown();
}

static void testStThreadSet(CuTest *testCase) {
    setup();
    teardown();
}

static void testStThread(CuTest *testCase) {
    setup();
    teardown();
}

static void testStSegment(CuTest *testCase) {
    setup();
    teardown();
}

static void testStBlock(CuTest *testCase) {
    setup();
    teardown();
}

CuSuite* stPinchGraphsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testStThreadSet);
    SUITE_ADD_TEST(suite, testStThread);
    SUITE_ADD_TEST(suite, testStSegment);
    SUITE_ADD_TEST(suite, testStBlock);
    return suite;
}
