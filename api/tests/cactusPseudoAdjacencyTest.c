#include "cactusReferenceTestShared.h"

static void testTeardown() {
	if(!nestedTest) {
		cactusReferenceTestSharedTeardown();
	}
}

static void testSetup() {
	if(!nestedTest) {
		cactusReferenceTestSharedSetup();
	}
}

void testPseudoAdjacency_construct(CuTest* testCase) {
	testSetup();
	assert(testCase != NULL);
	//already tested by shared code.
	testTeardown();
}

void testPseudoAdjacency_getName(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, pseudoAdjacency_getName(pseudoAdjacency1) != NULL_NAME);
	CuAssertTrue(testCase, pseudoAdjacency_getName(pseudoAdjacency2) != NULL_NAME);
	CuAssertTrue(testCase, pseudoAdjacency_getName(pseudoAdjacency1) != pseudoAdjacency_getName(pseudoAdjacency2));
	testTeardown();
}

void testPseudoAdjacency_get5End(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, pseudoAdjacency_get5End(pseudoAdjacency1) == end1);
	CuAssertTrue(testCase, pseudoAdjacency_get5End(pseudoAdjacency2) == end3);
	CuAssertTrue(testCase, pseudoAdjacency_get5End(pseudoAdjacency3) == end5);
	CuAssertTrue(testCase, pseudoAdjacency_get5End(pseudoAdjacency4) == end7);
	testTeardown();
}

void testPseudoAdjacency_get3End(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, pseudoAdjacency_get3End(pseudoAdjacency1) == end2);
	CuAssertTrue(testCase, pseudoAdjacency_get3End(pseudoAdjacency2) == end4);
	CuAssertTrue(testCase, pseudoAdjacency_get3End(pseudoAdjacency3) == end6);
	CuAssertTrue(testCase, pseudoAdjacency_get3End(pseudoAdjacency4) == end8);
	testTeardown();
}

void testPseudoAdjacency_getPseudoChromosome(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, pseudoAdjacency_getPseudoChromosome(pseudoAdjacency1) == pseudoChromosome1);
	CuAssertTrue(testCase, pseudoAdjacency_getPseudoChromosome(pseudoAdjacency2) == pseudoChromosome1);
	CuAssertTrue(testCase, pseudoAdjacency_getPseudoChromosome(pseudoAdjacency3) == pseudoChromosome1);
	CuAssertTrue(testCase, pseudoAdjacency_getPseudoChromosome(pseudoAdjacency4) == pseudoChromosome2);
	testTeardown();
}

void testPseudoAdjacency_serialisation(CuTest* testCase) {
	testSetup();
	int32_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(pseudoAdjacency1,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))pseudoAdjacency_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	pseudoAdjacency_destruct(pseudoAdjacency1);
	void *vA2 = vA;
	pseudoAdjacency1 = pseudoAdjacency_loadFromBinaryRepresentation(&vA2, pseudoChromosome1);
	free(vA);
	nestedTest = 1;

	testPseudoAdjacency_getName(testCase);
	testPseudoAdjacency_get5End(testCase);
	testPseudoAdjacency_get3End(testCase);
	testPseudoAdjacency_getPseudoChromosome(testCase);
	testPseudoAdjacency_construct(testCase);

	nestedTest = 0;
	testTeardown();
}

CuSuite* cactusPseudoAdjacencyTestSuite(void) {
	CuSuite* suite = CuSuiteNew();

	SUITE_ADD_TEST(suite, testPseudoAdjacency_getName);
	SUITE_ADD_TEST(suite, testPseudoAdjacency_get5End);
	SUITE_ADD_TEST(suite, testPseudoAdjacency_get3End);
	SUITE_ADD_TEST(suite, testPseudoAdjacency_getPseudoChromosome);
	SUITE_ADD_TEST(suite, testPseudoAdjacency_serialisation);
	SUITE_ADD_TEST(suite, testPseudoAdjacency_construct);

	return suite;
}
