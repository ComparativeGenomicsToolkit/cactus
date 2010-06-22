#include "cactusGlobalsPrivate.h"

CuSuite* stPosetAlignmentTestSuite(void);

int stBaseAlignerRunAllTests(void) {
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();
	CuSuiteAddSuite(suite, stPosetAlignmentTestSuite());
	CuSuiteRun(suite);
	CuSuiteSummary(suite, output);
	CuSuiteDetails(suite, output);
	printf("%s\n", output->buffer);
	return suite->failCount > 0;
}

int main(void) {
	return stBaseAlignerRunAllTests();
}
