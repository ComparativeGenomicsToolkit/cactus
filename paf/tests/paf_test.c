#include "paf.h"
#include "CuTest.h"
#include "sonLib.h"

/*
 * Test paf library functions
 */

static char *test_paf_file = "./tests/human_chimp.paf";

static void test_paf(CuTest *testCase) {
    // Read the pafs from the test file
    FILE *fh = fopen(test_paf_file, "r");
    stList *pafs = read_pafs(fh);
    fclose(fh);

    // Check we have the right number of records / that they look okay
    CuAssertIntEquals(testCase, stList_length(pafs), 207);

    // Write the pafs to a different file
    char *test_paf_2 = "./tests/human_chimp_copy.paf";
    fh = fopen(test_paf_2, "w");
    write_pafs(fh, pafs);
    fclose(fh);

    // Read them back from the copied file
    fh = fopen(test_paf_2, "r");
    stList *pafs2 = read_pafs(fh);
    fclose(fh);
    st_system("rm -f %s", test_paf_2); // Remove the copied file

    // Check the paf records are the same
    CuAssertTrue(testCase, stList_length(pafs) == stList_length(pafs2));
    for(int64_t i=0; i<stList_length(pafs); i++) {
        Paf *paf1 = stList_get(pafs, i);
        Paf *paf2 = stList_get(pafs2, i);

        paf_check(paf1);

        st_logDebug("Checking paf1: %s\n", paf_print(paf1));
        st_logDebug("Checking paf2: %s\n", paf_print(paf2));

        // Check they are equals
        CuAssertTrue(testCase, strcmp(paf_print(paf1), paf_print(paf2)) == 0);
    }
}

CuSuite* addPafTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_paf);
    return suite;
}
