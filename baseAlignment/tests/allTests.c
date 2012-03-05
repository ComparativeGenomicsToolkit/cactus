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

CuSuite* stPosetAlignmentTestSuite(void);
CuSuite* pairwiseAlignmentTestSuite(void);
CuSuite* multipleAlignerTestSuite(void);
CuSuite* adjacencySequenceTestSuite(void);
CuSuite* endAlignerTestSuite(void);
CuSuite* flowerAlignerTestSuite(void);
CuSuite* pairwiseAlignmentLongTestSuite(void);
CuSuite* pairwiseAlignmentTestSuite(void);

int stBaseAlignerRunAllTests(void) {
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();
	CuSuiteAddSuite(suite, pairwiseAlignmentTestSuite());
	/*CuSuiteAddSuite(suite, stPosetAlignmentTestSuite());
	CuSuiteAddSuite(suite, multipleAlignerTestSuite());
	CuSuiteAddSuite(suite, adjacencySequenceTestSuite());
	CuSuiteAddSuite(suite, endAlignerTestSuite());
	CuSuiteAddSuite(suite, flowerAlignerTestSuite());
	CuSuiteAddSuite(suite, pairwiseAlignmentLongTestSuite());*/
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
	return stBaseAlignerRunAllTests();
}
