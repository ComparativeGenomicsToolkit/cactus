/*
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"

CuSuite* addFastaExtractTestSuite(void);
CuSuite* addFastaChunkAndMergeTestSuite(void);

int cactusFastaRunAllTests(void) {
    CuString *output = CuStringNew();
    CuSuite* suite = CuSuiteNew();
    CuSuiteAddSuite(suite, addFastaExtractTestSuite());
    CuSuiteAddSuite(suite, addFastaChunkAndMergeTestSuite());
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    return suite->failCount > 0;
}

int main(int argc, char *argv[]) {
    if(argc == 2) {
        st_setLogLevelFromString(argv[1]);
    }
    int i = cactusFastaRunAllTests();
    //while(1);
    return i;
}
