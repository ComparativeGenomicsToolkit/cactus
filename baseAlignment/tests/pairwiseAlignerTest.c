/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void test_diagonal(CuTest *testCase) {

}

static void test_bandIterator(CuTest *testCase) {

}

/*
 * Test the math functions.
 */

double logAdd(double x, double y);

static void test_logAdd(CuTest *testCase) {
    for (int32_t test = 0; test < 100000; test++) {
        double i = st_random();
        double j = st_random();
        double k = i + j;
        double l = exp(logAdd(log(i), log(j)));
        //st_logInfo("I got %f %f\n", k, l);
        CuAssertTrue(testCase, l < k + 0.001);
        CuAssertTrue(testCase, l > k - 0.001);
    }
}

static void test_cell(CuTest *testCase) {

}

static void test_dpDiagonal(CuTest *testCase) {

}

static void test_dpMatrix(CuTest *testCase) {

}

static void test_diagonalCalculations(CuTest *testCase) {

}

/*
 * Test the sequence functions.
 */

char *convertSequence(const char *s, int32_t sL);

static void test_convertSequence(CuTest *testCase) {
    char cA[10] = { 0, 1, 2, 3, 4, 3, 4, 1, 2, '\0' };
    CuAssertStrEquals(testCase, cA, convertSequence("AcGTntNCG", 9));
    CuAssertStrEquals(testCase, cA, convertSequence("aCGTntNcg", 9));
}

char getRandomChar() {
    char *positions = "AaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtAaCcGgTtN";
    return positions[st_randomInt(0, strlen(positions))];
}

/*
 * Creates a random DNA sequence of the given length.
 */
char *getRandomSequence(int32_t length) {
    char *seq = st_malloc((length + 1) * sizeof(char));
    for (int32_t i = 0; i < length; i++) {
        seq[i] = getRandomChar();
    }
    seq[length] = '\0';
    return seq;
}

/*
 * Transfroms the given sequence into a different sequence.
 */
char *evolveSequence(const char *startSequence) {
    //Copy sequence
    char *seq = stString_copy(startSequence);

    //Do substitutions
    for (int32_t i = 0; i < strlen(seq); i++) {
        if (st_random() > 0.8) {
            seq[i] = getRandomChar();
        }
    }

    //Do indels
    while (st_random() > 0.2) {
        char *toReplace = getRandomSequence(st_randomInt(2, 4));
        char *replacement = getRandomSequence(st_randomInt(0, 10));
        char *seq2 = stString_replace(seq, toReplace, replacement);
        free(seq);
        free(toReplace);
        free(replacement);
        seq = seq2;
    }

    return seq;
}

static void test_pairwiseAlignerRandom(CuTest *testCase) {
    for (int32_t test = 0; test < 100; test++) {
        //Make a pair of sequences
        char *seqX = getRandomSequence(st_randomInt(0, 100));
        char *seqY = evolveSequence(seqX); //stString_copy(seqX);
        int32_t seqXLength = strlen(seqX);
        int32_t seqYLength = strlen(seqY);
        st_logInfo("Sequence X to align: %s END\n", seqX);
        st_logInfo("Sequence Y to align: %s END\n", seqY);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->alignAmbiguityCharacters = st_random() > 0.5; //Do this stochastically.
        stList *alignedPairs = getAlignedPairs(seqX, seqY, p);
        //Check the aligned pairs.
        stListIterator *iterator = stList_getIterator(alignedPairs);
        stIntTuple *alignedPair;
        while ((alignedPair = stList_getNext(iterator)) != NULL) {
            CuAssertTrue(testCase, stIntTuple_length(alignedPair) == 3);
            int32_t score = stIntTuple_getPosition(alignedPair, 0);
            int32_t x = stIntTuple_getPosition(alignedPair, 1);
            int32_t y = stIntTuple_getPosition(alignedPair, 2);
            //st_logInfo("Got aligned pair, score: %i x pos: %i y pos: %i\n", score, x, y);
            CuAssertTrue(testCase, score > 0);
            CuAssertTrue(testCase, score <= PAIR_ALIGNMENT_PROB_1);
            CuAssertTrue(testCase, x >= 0);
            CuAssertTrue(testCase, x < seqXLength);
            CuAssertTrue(testCase, y >= 0);
            CuAssertTrue(testCase, y < seqYLength);
        }
        stList_destructIterator(iterator);

        //Cleanup
        free(seqX);
        free(seqY);
        stList_destruct(alignedPairs);
    }
}

/*
 * Test the blast heuristic to get the different pairs.
 */

stList *getBlastPairs(const char *sX, const char *sY, int32_t lX, int32_t lY, int32_t trim, bool repeatMask);


static void test_getBlastPairs(CuTest *testCase) {
    for (int32_t test = 0; test < 10; test++) {
        //Make a pair of sequences
        char *seqX = getRandomSequence(st_randomInt(0, 10000));
        char *seqY = evolveSequence(seqX); //stString_copy(seqX);
        int32_t seqXLength = strlen(seqX);
        int32_t seqYLength = strlen(seqY);
        st_logInfo("Sequence X to align: %s END, seq length %i\n", seqX, seqXLength);
        st_logInfo("Sequence Y to align: %s END, seq length %i\n", seqY, seqYLength);

        int32_t trim = st_randomInt(0, 5);
        bool recursive = st_random() > 0.5;
        st_logInfo("Using random trim %i, recursive %i \n", trim, recursive);

        stList *blastPairs = getBlastPairs(seqX, seqY, seqXLength, seqYLength, trim, recursive);

        st_logInfo("I got %i blast pairs\n", stList_length(blastPairs));
        int32_t pX = -1;
        int32_t pY = -1;
        for(int32_t i=0; i<stList_length(blastPairs); i++) {
            stIntTuple  *j = stList_get(blastPairs, i);
            CuAssertTrue(testCase, stIntTuple_length(j) == 2);
            int32_t x = stIntTuple_getPosition(j, 0);
            int32_t y = stIntTuple_getPosition(j, 1);
            CuAssertTrue(testCase, x >= 0);
            CuAssertTrue(testCase, y >= 0);
            CuAssertTrue(testCase, x < seqXLength);
            CuAssertTrue(testCase, y < seqYLength);
            CuAssertTrue(testCase, x > pX);
            CuAssertTrue(testCase, y > pY);
            pX = x;
            pY = y;
        }

        stList_destruct(blastPairs);
    }
}


CuSuite* pairwiseAlignmentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_diagonal);
    SUITE_ADD_TEST(suite, test_bandIterator);
    SUITE_ADD_TEST(suite, test_logAdd);
    SUITE_ADD_TEST(suite, test_cell);
    SUITE_ADD_TEST(suite, test_dpDiagonal);
    SUITE_ADD_TEST(suite, test_dpMatrix);
    SUITE_ADD_TEST(suite, test_diagonalCalculations);
    SUITE_ADD_TEST(suite, test_convertSequence);
    SUITE_ADD_TEST(suite, test_getBlastPairs);
    SUITE_ADD_TEST(suite, test_pairwiseAlignerRandom);
    return suite;
}
