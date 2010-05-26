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

void testReference_getEndToPseudoAdjacencyHash(CuTest* testCase) {
	testSetup();

	st_Hash *hash = reference_getEndToPseudoAdjacencyHash(reference);
	CuAssertTrue(testCase, st_hash_search(hash, end1) == pseudoAdjacency1);
	CuAssertTrue(testCase, st_hash_search(hash, end2) == pseudoAdjacency1);

	CuAssertTrue(testCase, st_hash_search(hash, end3) == pseudoAdjacency2);
	CuAssertTrue(testCase, st_hash_search(hash, end4) == pseudoAdjacency2);

	CuAssertTrue(testCase, st_hash_search(hash, end5) == pseudoAdjacency3);
	CuAssertTrue(testCase, st_hash_search(hash, end6) == pseudoAdjacency3);

	CuAssertTrue(testCase, st_hash_search(hash, end7) == pseudoAdjacency4);
	CuAssertTrue(testCase, st_hash_search(hash, end8) == pseudoAdjacency4);

	//Try static wrappers.
	CuAssertTrue(testCase, st_hash_search(hash, end_getStaticNameWrapper(end_getName(end1))) == pseudoAdjacency1);
	CuAssertTrue(testCase, st_hash_search(hash, end_getStaticNameWrapper(end_getName(end2))) == pseudoAdjacency1);

	st_hash_destruct(hash);

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

	testReference_getNet(testCase);
	testReference_getPseudoChromosomeNumber(testCase);
	testReference_getPseudoChromosome(testCase);
	testReference_getFirst(testCase);
	testReference_pseudoChromosomeIterator(testCase);
	testReference_construct(testCase);

	nestedTest = 0;
	testTeardown();
}

void testPrintCanonicalReference(CuTest *testCase) {
	testSetup();
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator = reference_getPseudoChromosomeIterator(reference);
	PseudoChromosome *pC;
	while((pC = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator = pseudoChromosome_getPseudoAdjacencyIterator(pC);
		PseudoAdjacency *pA;
		uglyf("Printing pseudo chromosome %s\n", netMisc_nameToStringStatic(pseudoChromosome_getName(pC)));
		uglyf("The 5 end of the pseudo-chromosome %s\n", netMisc_nameToStringStatic(end_getName(pseudoChromosome_get5End(pC))));
		while((pA = pseudoChromosome_getNextPseudoAdjacency(pseudoAdjacencyIterator)) != NULL) {
			uglyf("Printing pseudo adjacency %s\n", netMisc_nameToStringStatic(pseudoAdjacency_getName(pA)));
			uglyf("The pseudo adjacency containing 5 end: %s", netMisc_nameToStringStatic(end_getName(pseudoAdjacency_get5End(pA))));
			uglyf("and 3 end %s\n", netMisc_nameToStringStatic(end_getName(pseudoAdjacency_get3End(pA))));
		}
		uglyf("The 3 end of the pseudo-chromosome %s\n", netMisc_nameToStringStatic(end_getName(pseudoChromosome_get3End(pC))));
		pseudoChromosome_destructPseudoAdjacencyIterator(pseudoAdjacencyIterator);
	}
	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);

	testTeardown();
}

CuSuite* cactusReferenceTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testReference_getNet);
	SUITE_ADD_TEST(suite, testReference_getPseudoChromosomeNumber);
	SUITE_ADD_TEST(suite, testReference_getPseudoChromosome);
	SUITE_ADD_TEST(suite, testReference_getFirst);
	SUITE_ADD_TEST(suite, testReference_pseudoChromosomeIterator);
	SUITE_ADD_TEST(suite, testReference_serialisation);
	SUITE_ADD_TEST(suite, testPrintCanonicalReference);
	SUITE_ADD_TEST(suite, testReference_getEndToPseudoAdjacencyHash);
	SUITE_ADD_TEST(suite, testReference_construct);

	return suite;
}
