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
#include "stCactusGraphs.h"

static void testStCactusGraph_simple(CuTest *testCase) {

}

CuSuite* stCactusGraphsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testStCactusGraph_simple);
    return suite;
}
