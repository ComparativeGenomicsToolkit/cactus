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

void testReference_construct(CuTest* testCase) {
	testSetup();
	//tested the shared code
	testTeardown();
}

void testReference_getName(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, reference_getName(reference) != NULL_NAME);
	testTeardown();
}

void testReference_getNet(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, reference_getNet(reference) == net);
	testTeardown();
}

void testReference_getPseudoChromosomeNumber(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, reference_getPseudoChromosomeNumber(reference) == 2);
	testTeardown();
}

void testReference_getPseudoChromosome(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, reference_getPseudoChromosome(reference, pseudoChromosome_getName(pseudoChromosome1)) == pseudoChromosome1);
	CuAssertTrue(testCase, reference_getPseudoChromosome(reference, pseudoChromosome_getName(pseudoChromosome2)) == pseudoChromosome2);
	CuAssertTrue(testCase, reference_getPseudoChromosome(reference, NULL_NAME) == NULL);
	testTeardown();
}

void testReference_getFirst(CuTest* testCase) {
	testSetup();
	CuAssertTrue(testCase, reference_getFirst(reference) == pseudoChromosome1);
	testTeardown();
}

void testReference_pseudoChromosomeIterator(CuTest* testCase) {
	testSetup();
	Reference_PseudoChromosomeIterator *iterator;
	iterator = reference_getPseudoChromosomeIterator(reference);
	CuAssertTrue(testCase, reference_getNextPseudoChromosome(iterator) == pseudoChromosome1);
	Reference_PseudoChromosomeIterator *iterator2;
	iterator2 = reference_copyPseudoChromosomeIterator(iterator);

	CuAssertTrue(testCase, reference_getNextPseudoChromosome(iterator) == pseudoChromosome2);
	CuAssertTrue(testCase, reference_getNextPseudoChromosome(iterator2) == pseudoChromosome2);
	CuAssertTrue(testCase, reference_getNextPseudoChromosome(iterator) == NULL);
	CuAssertTrue(testCase, reference_getNextPseudoChromosome(iterator2) == NULL);

	CuAssertTrue(testCase, reference_getPreviousPseudoChromosome(iterator) == pseudoChromosome2);
	CuAssertTrue(testCase, reference_getPreviousPseudoChromosome(iterator2) == pseudoChromosome2);
	CuAssertTrue(testCase, reference_getPreviousPseudoChromosome(iterator) == pseudoChromosome1);
	CuAssertTrue(testCase, reference_getPreviousPseudoChromosome(iterator2) == pseudoChromosome1);
	CuAssertTrue(testCase, reference_getPreviousPseudoChromosome(iterator) == NULL);
	CuAssertTrue(testCase, reference_getPreviousPseudoChromosome(iterator2) == NULL);

	reference_destructPseudoChromosomeIterator(iterator);
	testTeardown();
}

void testReference_serialisation(CuTest* testCase) {
	testSetup();
	int32_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(reference,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))reference_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	Name name1 = pseudoChromosome_getName(pseudoChromosome1);
	Name name2 = pseudoChromosome_getName(pseudoChromosome2);
	reference_destruct(reference);
	void *vA2 = vA;
	reference = reference_loadFromBinaryRepresentation(&vA2, net);
	free(vA);
	//uglyf("Printing the following names %s %s\n", netMisc_nameToString(name1), netMisc_nameToString(name2));
	pseudoChromosome1 = reference_getPseudoChromosome(reference, name1);
	assert(pseudoChromosome1 != NULL);
	pseudoChromosome2 = reference_getPseudoChromosome(reference, name2);
	assert(pseudoChromosome2 != NULL);
	//adjacency references are now void.

	nestedTest = 1;

	testReference_getName(testCase);
	testReference_getNet(testCase);
	testReference_getPseudoChromosomeNumber(testCase);
	testReference_getPseudoChromosome(testCase);
	testReference_getFirst(testCase);
	testReference_pseudoChromosomeIterator(testCase);
	testReference_construct(testCase);

	nestedTest = 0;
	testTeardown();
}

CuSuite* cactusReferenceTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testReference_getName);
	SUITE_ADD_TEST(suite, testReference_getNet);
	SUITE_ADD_TEST(suite, testReference_getPseudoChromosomeNumber);
	SUITE_ADD_TEST(suite, testReference_getPseudoChromosome);
	SUITE_ADD_TEST(suite, testReference_getFirst);
	SUITE_ADD_TEST(suite, testReference_pseudoChromosomeIterator);
	SUITE_ADD_TEST(suite, testReference_serialisation);
	SUITE_ADD_TEST(suite, testReference_construct);

	return suite;
}
