#include "CuTest.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

CuSuite* threeEdgeTestSuite(void);


int threeEdgeRunAllTests(void) {
    CuString *output = CuStringNew();
    CuSuite* suite = CuSuiteNew();
    CuSuiteAddSuite(suite, threeEdgeTestSuite());

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    return suite->failCount > 0;
}

int main(void) {
    return threeEdgeRunAllTests();
}
