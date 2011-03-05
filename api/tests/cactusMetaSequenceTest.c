/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk;
Name eventName = 10;
static MetaSequence *metaSequence;
static const char *sequenceString = "ACTGGCACTG";
static const char *headerString = ">one";

static bool nestedTest = 0;

static void cactusMetaSequenceTestTeardown() {
	if(!nestedTest && cactusDisk != NULL) {
		testCommon_deleteTemporaryCactusDisk(cactusDisk);
		cactusDisk = NULL;
		metaSequence = NULL;
	}
}

static void cactusMetaSequenceTestSetup() {
	if(!nestedTest) {
		cactusMetaSequenceTestTeardown();
		cactusDisk = testCommon_getTemporaryCactusDisk();
		metaSequence = metaSequence_construct(1, 10, sequenceString,
					   headerString, eventName, cactusDisk);
	}
}

void testMetaSequence_getName(CuTest* testCase) {
	cactusMetaSequenceTestSetup();
	CuAssertTrue(testCase, metaSequence_getName(metaSequence) != NULL_NAME);
	CuAssertTrue(testCase, cactusDisk_getMetaSequence(cactusDisk, metaSequence_getName(metaSequence)) == metaSequence);
	cactusMetaSequenceTestTeardown();
}

void testMetaSequence_getStart(CuTest* testCase) {
	cactusMetaSequenceTestSetup();
	CuAssertIntEquals(testCase, 1, metaSequence_getStart(metaSequence));
	cactusMetaSequenceTestTeardown();
}

void testMetaSequence_getLength(CuTest* testCase) {
	cactusMetaSequenceTestSetup();
	CuAssertIntEquals(testCase, 10, metaSequence_getLength(metaSequence));
	cactusMetaSequenceTestTeardown();
}

void testMetaSequence_getEventName(CuTest* testCase) {
	cactusMetaSequenceTestSetup();
	CuAssertTrue(testCase, metaSequence_getEventName(metaSequence) == eventName);
	cactusMetaSequenceTestTeardown();
}

void testMetaSequence_getString(CuTest* testCase) {
    for(int32_t i=0; i<10; i++) {
        cactusMetaSequenceTestSetup();
        //String is ACTGGCACTG
        CuAssertStrEquals(testCase, sequenceString, metaSequence_getString(metaSequence, 1, 10, 1)); //complete sequence
        CuAssertStrEquals(testCase, "TGGC", metaSequence_getString(metaSequence, 3, 4, 1)); //sub range
        CuAssertStrEquals(testCase, "", metaSequence_getString(metaSequence, 3, 0, 1)); //zero length sub range
        CuAssertStrEquals(testCase, "CAGTGCCAGT", metaSequence_getString(metaSequence, 1, 10, 0)); //reverse complement
        CuAssertStrEquals(testCase, "GCCA", metaSequence_getString(metaSequence, 3, 4, 0)); //sub range, reverse complement
        CuAssertStrEquals(testCase, "", metaSequence_getString(metaSequence, 3, 0, 0)); //zero length sub range on reverse strand
        cactusMetaSequenceTestTeardown();
    }
}

void testMetaSequence_getHeader(CuTest* testCase) {
	cactusMetaSequenceTestSetup();
	CuAssertStrEquals(testCase, headerString, metaSequence_getHeader(metaSequence));
	cactusMetaSequenceTestTeardown();
}

void testMetaSequence_serialisation(CuTest* testCase) {
	cactusMetaSequenceTestSetup();
	int32_t i;
	Name name = metaSequence_getName(metaSequence);
	CuAssertTrue(testCase, cactusDisk_getMetaSequence(cactusDisk, name) == metaSequence);
	void *vA = binaryRepresentation_makeBinaryRepresentation(metaSequence,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))metaSequence_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	metaSequence_destruct(metaSequence);
	CuAssertTrue(testCase, cactusDisk_getMetaSequence(cactusDisk, name) == NULL);
	void *vA2 = vA;
	metaSequence = metaSequence_loadFromBinaryRepresentation(&vA2, cactusDisk);
	CuAssertTrue(testCase, name == metaSequence_getName(metaSequence));
	CuAssertStrEquals(testCase, headerString, metaSequence_getHeader(metaSequence));
	CuAssertTrue(testCase, cactusDisk_getMetaSequence(cactusDisk, name) == metaSequence);
	cactusDisk_write(cactusDisk);
	metaSequence_destruct(metaSequence);
	CuAssertTrue(testCase, cactusDisk_getMetaSequence(cactusDisk, name) != NULL);
	metaSequence = cactusDisk_getMetaSequence(cactusDisk, name);
	nestedTest = 1;
	testMetaSequence_getName(testCase);
	testMetaSequence_getStart(testCase);
	testMetaSequence_getLength(testCase);
	testMetaSequence_getEventName(testCase);
	testMetaSequence_getString(testCase);
	testMetaSequence_getHeader(testCase);
	nestedTest = 0;
	cactusMetaSequenceTestTeardown();
}

CuSuite* cactusMetaSequenceTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testMetaSequence_getName);
	SUITE_ADD_TEST(suite, testMetaSequence_getStart);
	SUITE_ADD_TEST(suite, testMetaSequence_getLength);
	SUITE_ADD_TEST(suite, testMetaSequence_getEventName);
	SUITE_ADD_TEST(suite, testMetaSequence_getString);
	SUITE_ADD_TEST(suite, testMetaSequence_serialisation);
	SUITE_ADD_TEST(suite, testMetaSequence_getHeader);
	return suite;
}
