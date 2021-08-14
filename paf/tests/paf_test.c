#include "paf.h"
#include "CuTest.h"
#include "sonLib.h"

/*
 * Test paf library functions
 */

static void test_paf(CuTest *testCase) {

}

CuSuite* addPafTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_paf);

    return suite;
}