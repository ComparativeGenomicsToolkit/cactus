#include "cactusGlobalsPrivate.h"

CuSuite* stPosetAlignmentTestSuite(void);
CuSuite* pairwiseAlignmentTestSuite(void);
CuSuite* multipleAlignerTestSuite(void);
CuSuite* adjacencySequenceTestSuite(void);
CuSuite* endAlignerTestSuite(void);
CuSuite* netAlignerTestSuite(void);

int stBaseAlignerRunAllTests(void) {
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();
	CuSuiteAddSuite(suite, stPosetAlignmentTestSuite());
	CuSuiteAddSuite(suite, pairwiseAlignmentTestSuite());
	CuSuiteAddSuite(suite, multipleAlignerTestSuite());
	CuSuiteAddSuite(suite, adjacencySequenceTestSuite());
	CuSuiteAddSuite(suite, endAlignerTestSuite());
	CuSuiteAddSuite(suite, netAlignerTestSuite());
	CuSuiteRun(suite);
	CuSuiteSummary(suite, output);
	CuSuiteDetails(suite, output);
	printf("%s\n", output->buffer);
	return suite->failCount > 0;
}

int main(void) {
	return stBaseAlignerRunAllTests();
}
