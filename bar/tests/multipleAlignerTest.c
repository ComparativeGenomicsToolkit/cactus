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
 * Test the multiple alignment code with multiple examples.
 */

//Declarations from pairwiseAlignerTest

stList *getRandomSequences(int32_t sequenceNumber, int32_t approxLength) {
    stList *sequences = stList_construct3(0, free);
    char *firstSequence = getRandomSequence(approxLength);
    for(int32_t i=0; i<sequenceNumber; i++) {
        stList_append(sequences, evolveSequence(firstSequence));
    }
    return sequences;
}

void test_multipleAlignerRandom(CuTest *testCase) {
    for(int32_t test=0; test<100; test++) {
        stList *randomSequences = getRandomSequences(st_randomInt(0, 10), st_randomInt(0, 100));
        int32_t spanningTrees = st_randomInt(0, 5);
        stPosetAlignment *posetAlignment = stPosetAlignment_construct(stList_length(randomSequences));

        for(int32_t i=0; i<stList_length(randomSequences); i++) {
            st_logInfo("Sequence to align: %s\n", stList_get(randomSequences, i));
        }

        PairwiseAlignmentParameters *pairwiseParameters = pairwiseAlignmentBandingParameters_construct();
        stList *alignedPairs = makeAlignment(randomSequences, spanningTrees, 10000000, 0.5, pairwiseParameters);
        pairwiseAlignmentBandingParameters_destruct(pairwiseParameters);
        //Check the aligned pairs.
        stListIterator *iterator = stList_getIterator(alignedPairs);
        stIntTuple *alignedPair;
        while((alignedPair = stList_getNext(iterator)) != NULL) {
            CuAssertTrue(testCase, stIntTuple_length(alignedPair) == 5);
            int32_t score = stIntTuple_getPosition(alignedPair, 0);
            int32_t seqX = stIntTuple_getPosition(alignedPair, 1);
            int32_t x = stIntTuple_getPosition(alignedPair, 2);
            int32_t seqY = stIntTuple_getPosition(alignedPair, 3);
            int32_t y = stIntTuple_getPosition(alignedPair, 4);
            st_logInfo("Got aligned pair, score: %i x seq: %i x pos: %i x seq: %i y pos: %i\n", score, seqX, x, seqY, y);
            CuAssertTrue(testCase, score > 0);
            CuAssertTrue(testCase, score <= PAIR_ALIGNMENT_PROB_1);
            CuAssertTrue(testCase, seqX >= 0);
            CuAssertTrue(testCase, seqX < stList_length(randomSequences));
            CuAssertTrue(testCase, x >= 0);
            CuAssertTrue(testCase, x < strlen(stList_get(randomSequences, seqX)));
            CuAssertTrue(testCase, seqY >= 0);
            CuAssertTrue(testCase, seqY < stList_length(randomSequences));
            CuAssertTrue(testCase, y >= 0);
            CuAssertTrue(testCase, y < strlen(stList_get(randomSequences, seqY)));
            //Check we can form an alignment
            CuAssertTrue(testCase, stPosetAlignment_add(posetAlignment, seqX, x, seqY, y));
        }
        stList_destructIterator(iterator);
        stList_destruct(randomSequences);
        stPosetAlignment_destruct(posetAlignment);
    }
}

CuSuite* multipleAlignerTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_multipleAlignerRandom);

    return suite;
}
