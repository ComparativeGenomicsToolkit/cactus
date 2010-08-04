#include "cactusGlobalsPrivate.h"

CuSuite *cactusGroupTestSuite();
CuSuite *cactusSegmentTestSuite();
CuSuite *cactusBlockTestSuite();
CuSuite *cactusChainTestSuite();
CuSuite *cactusDatabaseTestSuite();
CuSuite *cactusCapTestSuite();
CuSuite *cactusEndTestSuite();
CuSuite *cactusEventTestSuite();
CuSuite *cactusEventTreeTestSuite();
CuSuite *cactusLinkTestSuite();
CuSuite *cactusMetaEventTestSuite();
CuSuite *cactusMetaSequenceTestSuite();
CuSuite *cactusDiskTestSuite();
CuSuite *cactusMiscTestSuite();
CuSuite *cactusNetTestSuite();
CuSuite *cactusFaceTestSuite();
CuSuite *cactusFaceEndTestSuite();
CuSuite *cactusSequenceTestSuite();
CuSuite *cactusSerialisationTestSuite();
CuSuite *cactusPseudoAdjacencyTestSuite();
CuSuite *cactusPseudoChromosomeTestSuite();
CuSuite *cactusReferenceTestSuite();


int cactusAPIRunAllTests(void) {
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();
	CuSuiteAddSuite(suite, cactusEventTestSuite());
	CuSuiteAddSuite(suite, cactusGroupTestSuite());
	CuSuiteAddSuite(suite, cactusSegmentTestSuite());
	CuSuiteAddSuite(suite, cactusBlockTestSuite());
	CuSuiteAddSuite(suite, cactusChainTestSuite());
	CuSuiteAddSuite(suite, cactusDatabaseTestSuite());
	CuSuiteAddSuite(suite, cactusCapTestSuite());
	CuSuiteAddSuite(suite, cactusEndTestSuite());
	CuSuiteAddSuite(suite, cactusEventTreeTestSuite());
	CuSuiteAddSuite(suite, cactusLinkTestSuite());
	CuSuiteAddSuite(suite, cactusMetaEventTestSuite());
	CuSuiteAddSuite(suite, cactusMetaSequenceTestSuite());
	CuSuiteAddSuite(suite, cactusDiskTestSuite());
	CuSuiteAddSuite(suite, cactusMiscTestSuite());
	CuSuiteAddSuite(suite, cactusNetTestSuite());
	CuSuiteAddSuite(suite, cactusFaceTestSuite());
	CuSuiteAddSuite(suite, cactusFaceEndTestSuite());
	CuSuiteAddSuite(suite, cactusSequenceTestSuite());
	CuSuiteAddSuite(suite, cactusSerialisationTestSuite());
	CuSuiteAddSuite(suite, cactusPseudoAdjacencyTestSuite());
	CuSuiteAddSuite(suite, cactusPseudoChromosomeTestSuite());
	CuSuiteAddSuite(suite, cactusReferenceTestSuite());
	CuSuiteRun(suite);
	CuSuiteSummary(suite, output);
	CuSuiteDetails(suite, output);
	printf("%s\n", output->buffer);
	return suite->failCount > 0;
}

int main(void) {
	return cactusAPIRunAllTests();
}
