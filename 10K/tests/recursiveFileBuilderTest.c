/*
 * recursiveFileBuilderTest.c
 *
 *  Created on: 19 Mar 2012
 *      Author: benedictpaten
 */

#include <stdlib.h>

#include "sonLib.h"
#include "cactus.h"
#include "CuTest.h"

static void recursiveFileBuilder_test(CuTest *testCase) {

}

CuSuite* recursiveFileBuilderTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, recursiveFileBuilder_test);
    return suite;
}
