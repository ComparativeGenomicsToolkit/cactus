/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sonLib.h"

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

int main(int argc, char *argv[]) {
    if(argc == 2) {
        st_setLogLevelFromString(argv[1]);
    }
    return threeEdgeRunAllTests();
}
