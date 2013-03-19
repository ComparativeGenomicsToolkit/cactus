/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "multipleAligner.h"
#include "stPosetAlignment.h"
#include "pairwiseAligner.h"
#include "randomSequences.h"

#include <stdlib.h>
#include <string.h>

/*
 * Test the multiple alignment code.
 */

static const char *seq1 = "AGTTC";
static const char *seq2 = "AGTTC";
static const char *seq3 = "AC";
static const char *seq4 = "";
static stList *littleSequences = NULL;
static PairwiseAlignmentParameters *pabp = NULL;

static void teardown() {
    if (littleSequences != NULL) {
        stList_destruct(littleSequences);
        littleSequences = NULL;
        pairwiseAlignmentBandingParameters_destruct(pabp);
        pabp = NULL;
    }
}

static void setup() {
    teardown();
    littleSequences = stList_construct();
    pabp = pairwiseAlignmentBandingParameters_construct();
    stList_append(littleSequences, (char *)seq1);
    stList_append(littleSequences, (char *)seq2);
    stList_append(littleSequences, (char *)seq3);
    stList_append(littleSequences, (char *)seq4);
}

stSet *makeColumns(stList *sequences);

static void test_makeColumns(CuTest *testCase) {
    setup();
    stSet *columns = makeColumns(littleSequences);
    CuAssertIntEquals(testCase, 12, stSet_size(columns));
    stSet_destruct(columns);
    teardown();
}

static void checkAlignment(CuTest *testCase, stList *sequences, stList *multipleAlignedPairs) {
    stPosetAlignment *posetAlignment = stPosetAlignment_construct(stList_length(sequences));
    //Check the aligned pairs.
    stListIterator *iterator = stList_getIterator(multipleAlignedPairs);
    stIntTuple *multipleAlignedPair;
    while ((multipleAlignedPair = stList_getNext(iterator)) != NULL) {
        CuAssertTrue(testCase, stIntTuple_length(multipleAlignedPair) == 5);
        int32_t score = stIntTuple_getPosition(multipleAlignedPair, 0);
        int32_t seqX = stIntTuple_getPosition(multipleAlignedPair, 1);
        int32_t x = stIntTuple_getPosition(multipleAlignedPair, 2);
        int32_t seqY = stIntTuple_getPosition(multipleAlignedPair, 3);
        int32_t y = stIntTuple_getPosition(multipleAlignedPair, 4);
        st_logInfo("Got aligned pair, score: %i x seq: %i x pos: %i x seq: %i y pos: %i\n", score, seqX, x, seqY, y);
        CuAssertTrue(testCase, score > 0);
        CuAssertTrue(testCase, score <= PAIR_ALIGNMENT_PROB_1);
        CuAssertTrue(testCase, seqX >= 0);
        CuAssertTrue(testCase, seqX < stList_length(sequences));
        CuAssertTrue(testCase, x >= 0);
        CuAssertTrue(testCase, x < strlen(stList_get(sequences, seqX)));
        CuAssertTrue(testCase, seqY >= 0);
        CuAssertTrue(testCase, seqY < stList_length(sequences));
        CuAssertTrue(testCase, y >= 0);
        CuAssertTrue(testCase, y < strlen(stList_get(sequences, seqY)));
        //Check we can form an alignment
        CuAssertTrue(testCase, stPosetAlignment_add(posetAlignment, seqX, x, seqY, y));
    }
    stList_destructIterator(iterator);
    stPosetAlignment_destruct(posetAlignment);
}

static void test_makeAlignmentUsingAllPairs(CuTest *testCase) {
    setup();
    stList *multipleAlignedPairs = makeAlignmentUsingAllPairs(littleSequences, 0.0, pabp);
    checkAlignment(testCase, littleSequences, multipleAlignedPairs);
    CuAssertIntEquals(testCase, 9, stList_length(multipleAlignedPairs));
    stList_destruct(multipleAlignedPairs);
    teardown();
}

stList *getRandomSequences(int32_t sequenceNumber, int32_t approxLength) {
    /*
     * Generate a random set of sequences.
     */
    stList *sequences = stList_construct3(0, free);
    char *firstSequence = getRandomSequence(approxLength);
    for (int32_t i = 0; i < sequenceNumber; i++) {
        stList_append(sequences, evolveSequence(firstSequence));
    }
    return sequences;
}

static void test_multipleAlignerAllPairsRandom(CuTest *testCase) {
    for (int32_t test = 0; test < 100; test++) {
        setup();
        stList *randomSequences = getRandomSequences(st_randomInt(0, 10), st_randomInt(0, 100));
        for (int32_t i = 0; i < stList_length(randomSequences); i++) {
            st_logInfo("Sequence to align: %s\n", stList_get(randomSequences, i));
        }
        stList *multipleAlignedPairs = makeAlignmentUsingAllPairs(randomSequences, 0.5, pabp);
        checkAlignment(testCase, randomSequences, multipleAlignedPairs);
        stList_destruct(randomSequences);
        teardown();
    }
}

static void test_multipleAlignerRandom(CuTest *testCase) {
    return;
    for (int32_t test = 0; test < 100; test++) {
        setup();
        stList *randomSequences = getRandomSequences(st_randomInt(0, 10), st_randomInt(0, 100));
        int32_t spanningTrees = st_randomInt(0, 5);
        for (int32_t i = 0; i < stList_length(randomSequences); i++) {
            st_logInfo("Sequence to align: %s\n", stList_get(randomSequences, i));
        }
        stList *multipleAlignedPairs = makeAlignment(randomSequences, spanningTrees, 10000000, 0.5, pabp);
        checkAlignment(testCase, randomSequences, multipleAlignedPairs);
        stList_destruct(randomSequences);
        teardown();
    }
}

CuSuite* multipleAlignerTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_makeColumns);
    SUITE_ADD_TEST(suite, test_makeAlignmentUsingAllPairs);
    SUITE_ADD_TEST(suite, test_multipleAlignerAllPairsRandom);
    SUITE_ADD_TEST(suite, test_multipleAlignerRandom);

    return suite;
}
