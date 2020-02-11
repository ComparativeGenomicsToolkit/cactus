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

static void cactusMetaSequenceTestTeardown(CuTest* testCase) {
    if(!nestedTest && cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
        cactusDisk = NULL;
        metaSequence = NULL;
    }
}

static void cactusMetaSequenceTestSetup(CuTest* testCase) {
    if(!nestedTest) {
        cactusMetaSequenceTestTeardown(testCase);
        cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
        metaSequence = metaSequence_construct(1, 10, sequenceString,
                       headerString, eventName, cactusDisk);
    }
}

void testMetaSequence_getName(CuTest* testCase) {
    cactusMetaSequenceTestSetup(testCase);
    CuAssertTrue(testCase, metaSequence_getName(metaSequence) != NULL_NAME);
    CuAssertTrue(testCase, cactusDisk_getMetaSequence(cactusDisk, metaSequence_getName(metaSequence)) == metaSequence);
    cactusMetaSequenceTestTeardown(testCase);
}

void testMetaSequence_getStart(CuTest* testCase) {
    cactusMetaSequenceTestSetup(testCase);
    CuAssertIntEquals(testCase, 1, metaSequence_getStart(metaSequence));
    cactusMetaSequenceTestTeardown(testCase);
}

void testMetaSequence_getLength(CuTest* testCase) {
    cactusMetaSequenceTestSetup(testCase);
    CuAssertIntEquals(testCase, 10, metaSequence_getLength(metaSequence));
    cactusMetaSequenceTestTeardown(testCase);
}

void testMetaSequence_getEventName(CuTest* testCase) {
    cactusMetaSequenceTestSetup(testCase);
    CuAssertTrue(testCase, metaSequence_getEventName(metaSequence) == eventName);
    cactusMetaSequenceTestTeardown(testCase);
}

void testMetaSequence_getString(CuTest* testCase) {
    for(int64_t i=0; i<10; i++) {
        cactusMetaSequenceTestSetup(testCase);
        //String is ACTGGCACTG
        CuAssertStrEquals(testCase, sequenceString, metaSequence_getString(metaSequence, 1, 10, 1)); //complete sequence
        CuAssertStrEquals(testCase, "TGGC", metaSequence_getString(metaSequence, 3, 4, 1)); //sub range
        CuAssertStrEquals(testCase, "", metaSequence_getString(metaSequence, 3, 0, 1)); //zero length sub range
        CuAssertStrEquals(testCase, "CAGTGCCAGT", metaSequence_getString(metaSequence, 1, 10, 0)); //reverse complement
        CuAssertStrEquals(testCase, "GCCA", metaSequence_getString(metaSequence, 3, 4, 0)); //sub range, reverse complement
        CuAssertStrEquals(testCase, "", metaSequence_getString(metaSequence, 3, 0, 0)); //zero length sub range on reverse strand
        cactusMetaSequenceTestTeardown(testCase);
    }
}

void testMetaSequence_getHeader(CuTest* testCase) {
    cactusMetaSequenceTestSetup(testCase);
    CuAssertStrEquals(testCase, headerString, metaSequence_getHeader(metaSequence));
    cactusMetaSequenceTestTeardown(testCase);
}

void testMetaSequence_isTrivialSequence(CuTest* testCase) {
    cactusMetaSequenceTestSetup(testCase);
    CuAssertTrue(testCase, !metaSequence_isTrivialSequence(metaSequence));
    MetaSequence *metaSequence2 = metaSequence_construct3(1, 10, sequenceString,
                           headerString, eventName, 1, cactusDisk);
    CuAssertTrue(testCase, metaSequence_isTrivialSequence(metaSequence2));
    metaSequence_destruct(metaSequence2);
    cactusMetaSequenceTestTeardown(testCase);
}

void testMetaSequence_serialisation(CuTest* testCase) {
    cactusMetaSequenceTestSetup(testCase);
    int64_t i;
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
    testMetaSequence_isTrivialSequence(testCase);
    nestedTest = 0;
    cactusMetaSequenceTestTeardown(testCase);
}

CuSuite* cactusMetaSequenceTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testMetaSequence_getName);
    SUITE_ADD_TEST(suite, testMetaSequence_getStart);
    SUITE_ADD_TEST(suite, testMetaSequence_getLength);
    SUITE_ADD_TEST(suite, testMetaSequence_getEventName);
    SUITE_ADD_TEST(suite, testMetaSequence_getString);
    SUITE_ADD_TEST(suite, testMetaSequence_isTrivialSequence);
    SUITE_ADD_TEST(suite, testMetaSequence_serialisation);
    SUITE_ADD_TEST(suite, testMetaSequence_getHeader);
    return suite;
}
