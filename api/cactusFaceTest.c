#include "cactusFacesTestShared.h"

static bool nestedTest = 0;

void cactusFaceTestSetup() {
	if (!nestedTest) {
		cactusFacesTestSharedSetup();
	}
}

void cactusFaceTestTeardown() {
	if (!nestedTest) {
		cactusFacesTestSharedTeardown();
	}
}

void testFace_construct(CuTest* testCase) {
	nestedTest = 0;
	cactusFaceTestSetup();
	CuAssertTrue(testCase, face != NULL);
	cactusFaceTestTeardown();
}

void testFace_getCardinal(CuTest* testCase) {
	cactusFaceTestSetup();
	CuAssertTrue(testCase, face_getCardinal(face) == 4);
	cactusFaceTestTeardown();
}

void testFace_getTopNode(CuTest* testCase) {
	cactusFaceTestSetup();
	CuAssertTrue(testCase, face_getTopNode(face, 0) == cap_getPositiveOrientation(topCap1));
	CuAssertTrue(testCase, face_getTopNode(face, 1) == cap_getPositiveOrientation(topCap2));
	CuAssertTrue(testCase, face_getTopNode(face, 2) == cap_getPositiveOrientation(topCap3));
	CuAssertTrue(testCase, face_getTopNode(face, 3) == cap_getPositiveOrientation(topCap4));
	cactusFaceTestTeardown();
}

void testFace_getDerivedDestination(CuTest* testCase) {
	cactusFaceTestSetup();
	CuAssertTrue(testCase, face_getDerivedDestinationAtIndex(face, 0, 0) == cap_getPositiveOrientation(topCap4));
	CuAssertTrue(testCase, face_getDerivedDestinationAtIndex(face, 1, 0) == cap_getPositiveOrientation(topCap3));
	CuAssertTrue(testCase, face_getDerivedDestinationAtIndex(face, 2, 0) == cap_getPositiveOrientation(topCap2));
	CuAssertTrue(testCase, face_getDerivedDestinationAtIndex(face, 3, 0) == cap_getPositiveOrientation(topCap1));
	CuAssertTrue(testCase, face_getDerivedDestination(face, 0) == cap_getPositiveOrientation(topCap4));
	CuAssertTrue(testCase, face_getDerivedDestination(face, 1) == cap_getPositiveOrientation(topCap3));
	CuAssertTrue(testCase, face_getDerivedDestination(face, 2) == cap_getPositiveOrientation(topCap2));
	CuAssertTrue(testCase, face_getDerivedDestination(face, 3) == cap_getPositiveOrientation(topCap1));
	cactusFaceTestTeardown();
}

void testFace_getBottomNodeNumber(CuTest* testCase) {
	cactusFaceTestSetup();
	CuAssertTrue(testCase, face_getBottomNodeNumber(face, 0) == 1);
	CuAssertTrue(testCase, face_getBottomNodeNumber(face, 1) == 1);
	CuAssertTrue(testCase, face_getBottomNodeNumber(face, 2) == 1);
	CuAssertTrue(testCase, face_getBottomNodeNumber(face, 3) == 1);
	cactusFaceTestTeardown();
}

void testFace_getBottomNode(CuTest* testCase) {
	cactusFaceTestSetup();
	CuAssertTrue(testCase, face_getBottomNode(face, 0, 0) == cap_getPositiveOrientation(bottomCap1));
	CuAssertTrue(testCase, face_getBottomNode(face, 1, 0) == cap_getPositiveOrientation(bottomCap2));
	CuAssertTrue(testCase, face_getBottomNode(face, 2, 0) == cap_getPositiveOrientation(bottomCap3));
	CuAssertTrue(testCase, face_getBottomNode(face, 3, 0) == cap_getPositiveOrientation(bottomCap4));
	cactusFaceTestTeardown();
}

void testFace_getName(CuTest * testCase) {
	cactusFaceTestSetup();
	CuAssertTrue(testCase, face_getName(face) == face->name);
	cactusFaceTestTeardown();
}

void testFace_serialisation(CuTest* testCase) {
	cactusFaceTestSetup();
	int32_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(face,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))face_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	face_destruct(face);
	void *vA2 = vA;
	face = face_loadFromBinaryRepresentation(&vA2, net);
	free(vA);
	nestedTest = 1;
	testFace_getCardinal(testCase);
	testFace_getTopNode(testCase);
	testFace_getDerivedDestination(testCase);
	testFace_getBottomNodeNumber(testCase);
	testFace_getBottomNode(testCase);
	testFace_getName(testCase);
	nestedTest = 0;
	cactusFaceTestTeardown();
}

CuSuite *cactusFaceTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testFace_getTopNode);
	SUITE_ADD_TEST(suite, testFace_getCardinal);
	SUITE_ADD_TEST(suite, testFace_getDerivedDestination);
	SUITE_ADD_TEST(suite, testFace_getBottomNodeNumber);
	SUITE_ADD_TEST(suite, testFace_getBottomNode);
	SUITE_ADD_TEST(suite, testFace_serialisation);
	SUITE_ADD_TEST(suite, testFace_construct);
	SUITE_ADD_TEST(suite, testFace_getName);
	return suite;
}
