/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk;
static Event *event = NULL;
static Sequence *sequence;
static const char *sequenceString = "ACTGGCACTG";
static const char *headerString = ">one";

static bool nestedTest = 0;

static void cactusSequenceTestTeardown(CuTest* testCase) {
    if(!nestedTest && cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
        cactusDisk = NULL;
        sequence = NULL;
    }
}

static void cactusSequenceTestSetup(CuTest* testCase) {
    if(!nestedTest) {
        cactusSequenceTestTeardown(testCase);
        cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
        sequence = sequence_construct(1, 10, sequenceString,
                       headerString, event, cactusDisk);
    }
}

void testSequence_getName(CuTest* testCase) {
    cactusSequenceTestSetup(testCase);
    CuAssertTrue(testCase, sequence_getName(sequence) != NULL_NAME);
    CuAssertTrue(testCase, cactusDisk_getSequence(cactusDisk, sequence_getName(sequence)) == sequence);
    cactusSequenceTestTeardown(testCase);
}

void testSequence_getStart(CuTest* testCase) {
    cactusSequenceTestSetup(testCase);
    CuAssertIntEquals(testCase, 1, sequence_getStart(sequence));
    cactusSequenceTestTeardown(testCase);
}

void testSequence_getLength(CuTest* testCase) {
    cactusSequenceTestSetup(testCase);
    CuAssertIntEquals(testCase, 10, sequence_getLength(sequence));
    cactusSequenceTestTeardown(testCase);
}

void testSequence_getEvent(CuTest* testCase) {
    cactusSequenceTestSetup(testCase);
    CuAssertTrue(testCase, sequence_getEvent(sequence) == event);
    cactusSequenceTestTeardown(testCase);
}

void testSequence_getString(CuTest* testCase) {
    for(int64_t i=0; i<10; i++) {
        cactusSequenceTestSetup(testCase);
        //String is ACTGGCACTG
        CuAssertStrEquals(testCase, sequenceString, sequence_getString(sequence, 1, 10, 1)); //complete sequence
        CuAssertStrEquals(testCase, "TGGC", sequence_getString(sequence, 3, 4, 1)); //sub range
        CuAssertStrEquals(testCase, "", sequence_getString(sequence, 3, 0, 1)); //zero length sub range
        CuAssertStrEquals(testCase, "CAGTGCCAGT", sequence_getString(sequence, 1, 10, 0)); //reverse complement
        CuAssertStrEquals(testCase, "GCCA", sequence_getString(sequence, 3, 4, 0)); //sub range, reverse complement
        CuAssertStrEquals(testCase, "", sequence_getString(sequence, 3, 0, 0)); //zero length sub range on reverse strand
        cactusSequenceTestTeardown(testCase);
    }
}

void testSequence_getHeader(CuTest* testCase) {
    cactusSequenceTestSetup(testCase);
    CuAssertStrEquals(testCase, headerString, sequence_getHeader(sequence));
    cactusSequenceTestTeardown(testCase);
}

void testSequence_isTrivialSequence(CuTest* testCase) {
    cactusSequenceTestSetup(testCase);
    CuAssertTrue(testCase, !sequence_isTrivialSequence(sequence));
    Sequence *sequence2 = sequence_construct3(1, 10, sequenceString,
                           headerString, event, 1, cactusDisk);
    CuAssertTrue(testCase, sequence_isTrivialSequence(sequence2));
    sequence_destruct(sequence2);
    cactusSequenceTestTeardown(testCase);
}

void testSequence_serialisation(CuTest* testCase) {
    cactusSequenceTestSetup(testCase);
    int64_t i;
    Name name = sequence_getName(sequence);
    CuAssertTrue(testCase, cactusDisk_getSequence(cactusDisk, name) == sequence);
    void *vA = binaryRepresentation_makeBinaryRepresentation(sequence,
            (void (*)(void *, void (*)(const void *, size_t, size_t)))sequence_writeBinaryRepresentation, &i);
    CuAssertTrue(testCase, i > 0);
    sequence_destruct(sequence);
    CuAssertTrue(testCase, cactusDisk_getSequence(cactusDisk, name) == NULL);
    void *vA2 = vA;
    sequence = sequence_loadFromBinaryRepresentation(&vA2, cactusDisk);
    CuAssertTrue(testCase, name == sequence_getName(sequence));
    CuAssertStrEquals(testCase, headerString, sequence_getHeader(sequence));
    CuAssertTrue(testCase, cactusDisk_getSequence(cactusDisk, name) == sequence);
    cactusDisk_write(cactusDisk);
    sequence_destruct(sequence);
    CuAssertTrue(testCase, cactusDisk_getSequence(cactusDisk, name) != NULL);
    sequence = cactusDisk_getSequence(cactusDisk, name);
    nestedTest = 1;
    testSequence_getName(testCase);
    testSequence_getStart(testCase);
    testSequence_getLength(testCase);
    testSequence_getEvent(testCase);
    testSequence_getString(testCase);
    testSequence_getHeader(testCase);
    testSequence_isTrivialSequence(testCase);
    nestedTest = 0;
    cactusSequenceTestTeardown(testCase);
}

CuSuite* cactusSequenceTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testSequence_getName);
    SUITE_ADD_TEST(suite, testSequence_getStart);
    SUITE_ADD_TEST(suite, testSequence_getLength);
    SUITE_ADD_TEST(suite, testSequence_getEvent);
    SUITE_ADD_TEST(suite, testSequence_getString);
    SUITE_ADD_TEST(suite, testSequence_isTrivialSequence);
    SUITE_ADD_TEST(suite, testSequence_serialisation);
    SUITE_ADD_TEST(suite, testSequence_getHeader);
    return suite;
}
