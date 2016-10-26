/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"

CuSuite *buildReferenceTestSuite(void);
CuSuite* addReferenceCoordinatesTestSuite(void);
CuSuite* recursiveThreadBuilderTestSuite(void);

int referenceRunAllTests(void) {
    CuString *output = CuStringNew();
    CuSuite* suite = CuSuiteNew();
    CuSuiteAddSuite(suite, buildReferenceTestSuite());
    CuSuiteAddSuite(suite, addReferenceCoordinatesTestSuite());
    CuSuiteAddSuite(suite, recursiveThreadBuilderTestSuite());

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
    int i = referenceRunAllTests();
    //while(1);
    return i;
}
