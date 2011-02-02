/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static char vA[1000];
static char *vA3;

static void cactusSerialisationTestSetup() {
	vA3 = vA;
}

static void cactusSerialisationTestTeardown() {
}

static void writeFn(const void * ptr, size_t size, size_t count) {
	memcpy(vA3, ptr, size*count);
	vA3 += size*count;
}

void testBinaryRepresentation_elementType(CuTest* testCase) {
	cactusSerialisationTestSetup();
	void *vA2 = vA;
	binaryRepresentation_writeElementType(CODE_ADJACENCY, writeFn);
	binaryRepresentation_writeElementType(CODE_LINK, writeFn);
	CuAssertTrue(testCase, binaryRepresentation_peekNextElementType(vA2) == CODE_ADJACENCY);
	CuAssertTrue(testCase, binaryRepresentation_popNextElementType(&vA2) == CODE_ADJACENCY);
	CuAssertTrue(testCase, binaryRepresentation_peekNextElementType(vA2) == CODE_LINK);
	CuAssertTrue(testCase, binaryRepresentation_popNextElementType(&vA2) == CODE_LINK);
	cactusSerialisationTestTeardown();
}

void testBinaryRepresentation_string(CuTest* testCase) {
	cactusSerialisationTestSetup();
	void *vA2 = vA;
	binaryRepresentation_writeString("HELLO I AM A STRING", writeFn);
	binaryRepresentation_writeString("GOOD_BYE", writeFn);
	CuAssertStrEquals(testCase, "HELLO I AM A STRING", binaryRepresentation_getString(&vA2));
	CuAssertStrEquals(testCase, "GOOD_BYE", binaryRepresentation_getStringStatic(&vA2));
	cactusSerialisationTestTeardown();
}

void testBinaryRepresentation_integer(CuTest* testCase) {
	cactusSerialisationTestSetup();
	void *vA2 = vA;
	binaryRepresentation_writeInteger(537869, writeFn);
	binaryRepresentation_writeInteger(720032, writeFn);
	CuAssertIntEquals(testCase, 537869, binaryRepresentation_getInteger(&vA2));
	CuAssertIntEquals(testCase, 720032, binaryRepresentation_getInteger(&vA2));
	cactusSerialisationTestTeardown();
}

void testBinaryRepresentation_64BitInteger(CuTest* testCase) {
	cactusSerialisationTestSetup();
	void *vA2 = vA;
	int64_t i = 543829676894821452;
	int64_t j = 123456789876543234;
	binaryRepresentation_write64BitInteger(i, writeFn);
	binaryRepresentation_write64BitInteger(j, writeFn);
	CuAssertTrue(testCase, i == binaryRepresentation_get64BitInteger(&vA2));
	CuAssertTrue(testCase, j == binaryRepresentation_get64BitInteger(&vA2));
	cactusSerialisationTestTeardown();
}

void testBinaryRepresentation_name(CuTest* testCase) {
	cactusSerialisationTestSetup();
	void *vA2 = vA;
	Name name1 = 543829676894821452;
	Name name2 = 123456789876543234;
	binaryRepresentation_writeName(name1, writeFn);
	binaryRepresentation_writeName(name2, writeFn);
	CuAssertTrue(testCase, name1 == binaryRepresentation_getName(&vA2));
	CuAssertTrue(testCase, name2 == binaryRepresentation_getName(&vA2));
	cactusSerialisationTestTeardown();
}

void testBinaryRepresentation_float(CuTest* testCase) {
	cactusSerialisationTestSetup();
	void *vA2 = vA;
	float i = 3.145678;
	float j = 2.714342;
	binaryRepresentation_writeFloat(i, writeFn);
	binaryRepresentation_writeFloat(j, writeFn);
	CuAssertTrue(testCase, i == binaryRepresentation_getFloat(&vA2));
	CuAssertTrue(testCase, j == binaryRepresentation_getFloat(&vA2));
	cactusSerialisationTestTeardown();
}

void testBinaryRepresentation_bool(CuTest* testCase) {
	cactusSerialisationTestSetup();
	void *vA2 = vA;
	bool i = 0;
	bool j = 1;
	binaryRepresentation_writeBool(i, writeFn);
	binaryRepresentation_writeBool(j, writeFn);
	CuAssertTrue(testCase, i == binaryRepresentation_getBool(&vA2));
	CuAssertTrue(testCase, j == binaryRepresentation_getBool(&vA2));
	cactusSerialisationTestTeardown();
}

static void testBinaryRepresentation_fn(void *object, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeInteger(*(int32_t *)object, writeFn);
}

void testBinaryRepresentation_makeBinaryRepresentation(CuTest* testCase) {
	cactusSerialisationTestSetup();
	int32_t i, j;
	i = 14314;
	void *vA = binaryRepresentation_makeBinaryRepresentation(&i, testBinaryRepresentation_fn, &j);
	void *vA2 = vA;
	CuAssertTrue(testCase, j == sizeof(int32_t));
	CuAssertIntEquals(testCase, binaryRepresentation_getInteger(&vA2), i);
	free(vA);
	cactusSerialisationTestTeardown();
}

CuSuite* cactusSerialisationTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testBinaryRepresentation_elementType);
	SUITE_ADD_TEST(suite, testBinaryRepresentation_string);
	SUITE_ADD_TEST(suite, testBinaryRepresentation_integer);
	SUITE_ADD_TEST(suite, testBinaryRepresentation_64BitInteger);
	SUITE_ADD_TEST(suite, testBinaryRepresentation_name);
	SUITE_ADD_TEST(suite, testBinaryRepresentation_float);
	SUITE_ADD_TEST(suite, testBinaryRepresentation_bool);
	SUITE_ADD_TEST(suite, testBinaryRepresentation_makeBinaryRepresentation);
	return suite;
}
