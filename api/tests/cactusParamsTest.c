/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "cactus_params_parser.h"

static void testCactusParams(CuTest *testCase) {
    char *params_file = "/Users/benedictpaten/CLionProjects/cactus/src/cactus/cactus_progressive_config.xml";

    CactusParams *p = cactusParams_load(params_file);

    const char *c = cactusParams_get_string(p, 3, "blast", "divergence", "argName");
    CuAssertStrEquals(testCase, "lastzArguments", c);

    int64_t i = cactusParams_get_int(p, 3, "bar", "pecan", "spanningTrees");
    CuAssertIntEquals(testCase, 5, i);

    double d = cactusParams_get_float(p, 3, "reference", "CactusReferenceRecursion", "maxFlowerWrapperGroupSize");
    CuAssertDblEquals(testCase, 2000000, d, 0.0);

    int64_t length;
    int64_t *l = cactusParams_get_ints(p, &length, 2, "caf", "deannealingRounds");
    CuAssertIntEquals(testCase, 2, length);
    CuAssertIntEquals(testCase, l[0], 2);
    CuAssertIntEquals(testCase, l[1], 8);

    // Test moving the root of the search
    cactusParams_set_root(p, 1, "caf");

    i = cactusParams_get_int(p, 1, "chunkSize");
    CuAssertIntEquals(testCase, 25000000, i);

    // Check we can set it back
    cactusParams_set_root(p, 0);
    i = cactusParams_get_int(p, 3, "bar", "pecan", "spanningTrees");
    CuAssertIntEquals(testCase, 5, i);

    cactusParams_destruct(p); // Cleanup
    free(l);
}

CuSuite* cactusParamsTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testCactusParams);
    return suite;
}
