#include "cactusFacesTestShared.h"
/*
 * cactusFaceEndTest.c
 *
 *  Created on: 10-May-2010
 *      Author: benedictpaten
 */

static void testTearDown() {
	cactusFacesTestSharedTeardown();
}

static void testSetup() {
	cactusFacesTestSharedSetup();
}

void testFaceEnd_getTopNode(CuTest* testCase) {
	testSetup();

	testTearDown();
}

void testFaceEnd_getFace(CuTest* testCase) {
	testSetup();
	testTearDown();
}

void testFaceEnd_getNumberOfBottomNodes(CuTest* testCase) {
	testSetup();

	testTearDown();
}

void testFaceEnd_bottomNodeIterator(CuTest* testCase) {
	testSetup();
	testTearDown();
}

CuSuite* cactusFaceEndTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testFaceEnd_getTopNode);
	SUITE_ADD_TEST(suite, testFaceEnd_getFace);
	SUITE_ADD_TEST(suite, testFaceEnd_getNumberOfBottomNodes);
	SUITE_ADD_TEST(suite, testFaceEnd_bottomNodeIterator);
	return suite;
}
