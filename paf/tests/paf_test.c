#include "paf.h"
#include "CuTest.h"
#include "sonLib.h"

/*
 * Test paf library functions
 */

static char *test_paf_file = "./tests/human_chimp.paf";

/*
 * Read all the pafs from a file.
 */
stList *read_pafs(char *paf_file) {
    FILE *fh = fopen(paf_file, "r");
    stList *pafs = stList_construct3(0, (void (*)(void *))paf_destruct);
    Paf *paf;
    while((paf = paf_read(fh)) != NULL) {
        paf_check(paf);
        stList_append(pafs, paf);
    }
    fclose(fh);
    return pafs;
}

/*
 * Write a list of pafs to a file.
 */
void write_pafs(char *paf_file, stList *pafs) {
    FILE *fh = fopen(paf_file, "w");
    for(int64_t i=0; i<stList_length(pafs); i++) {
        paf_write(stList_get(pafs, i), fh);
    }
    fclose(fh);
}

static void test_paf(CuTest *testCase) {
    // Read the pafs from the test file
    stList *pafs = read_pafs(test_paf_file);

    // Check we have the right number of records / that they look okay
    CuAssertIntEquals(testCase, stList_length(pafs), 207);

    // Write the pafs to a different file
    char *test_paf_2 = "./tests/human_chimp_copy.paf";
    write_pafs(test_paf_2, pafs);

    // Read them back from the copied file
    stList *pafs2 = read_pafs(test_paf_2);
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
