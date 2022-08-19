/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "cactus_params_parser.h"

static char *params_file = "./src/cactus/cactus_progressive_config.xml";

static void testCactusParams(CuTest *testCase) {
    CactusParams *p = cactusParams_load(params_file);

    const char *c = cactusParams_get_string(p, 3, "blast", "divergence", "argName");
    CuAssertStrEquals(testCase, "lastzArguments", c);

    int64_t i = cactusParams_get_int(p, 3, "bar", "pecan", "spanningTrees");
    CuAssertIntEquals(testCase, 5, i);

    double d = cactusParams_get_float(p, 3, "bar", "poa", "partialOrderAlignmentBandFraction");
    CuAssertDblEquals(testCase, 0.05, d, 0.000001);

    int64_t length;
    int64_t *l = cactusParams_get_ints(p, &length, 2, "caf", "deannealingRounds");
    CuAssertIntEquals(testCase, 3, length);
    CuAssertIntEquals(testCase, l[0], 2);
    CuAssertIntEquals(testCase, l[1], 4);
    CuAssertIntEquals(testCase, l[2], 8);

    // Test moving the root of the search
    cactusParams_set_root(p, 1, "caf");

    i = cactusParams_get_int(p, 1, "trim");
    CuAssertIntEquals(testCase, 3, i);

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
