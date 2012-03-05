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
#include "randomSequences.h"

static void test_diagonal(CuTest *testCase) {
    //Construct an example diagonal.
    int32_t xL = 10, yL = 20, xU = 30, yU = 0; //Coordinates of the upper and lower
    //pairs in x,y coordinates
    Diagonal d = diagonal_construct(xL + yL, xL - yL, xU - yU);
    CuAssertIntEquals(testCase, diagonal_getXay(d), xL + yL);
    CuAssertIntEquals(testCase, diagonal_getMinXmy(d), xL - yL);
    CuAssertIntEquals(testCase, diagonal_getMaxXmy(d), xU - yU);
    CuAssertIntEquals(testCase, diagonal_getWidth(d), xU - yU - (xL - yL) + 1);
    CuAssertIntEquals(testCase, diagonal_getXCoordinate(xL + yL, xL - yL), xL);
    CuAssertIntEquals(testCase, diagonal_getYCoordinate(xL + yL, xL - yL), yL);
    CuAssertTrue(testCase, diagonal_equals(d, d));
    CuAssertTrue(testCase, !diagonal_equals(d, diagonal_construct(0, 0, 0)));
    //A bogus diagonal is one such that |xay + xmy| % 2 != 0 or such that xmyR < xmyL.
    //Try constructing bogus diagonals, should throw exceptions.
    stTry {
            diagonal_construct(10, 5, 5);
            CuAssertTrue(testCase, 0);
        }stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
            }stTryEnd
    stTry {
            diagonal_construct(10, 5, 5);
            CuAssertTrue(testCase, 0);
        }stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
            }stTryEnd
    stTry {
            diagonal_construct(10, 6, 4);
            CuAssertTrue(testCase, 0);
        }stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
            }stTryEnd
}

static bool testDiagonalsEqual(Diagonal d1, Diagonal d2) {
    bool b = diagonal_equals(d1, d2);
    if(!b) {
        st_logCritical("Diagonals not equal: d1: %s, d2: %s \n", diagonal_getString(d1), diagonal_getString(d2));
    }
    return b;
}

static void test_bandIterator(CuTest *testCase) {
    stList *anchorPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    stList_append(anchorPairs, stIntTuple_construct(2, 1, 0));
    stList_append(anchorPairs, stIntTuple_construct(2, 2, 1));
    stList_append(anchorPairs, stIntTuple_construct(2, 3, 3));
    //Start the traversal
    int32_t lX = 6, lY = 5;
    BandIterator *bandIt = bandIterator_construct(anchorPairs, lX, lY, 2);

    //Forward pass
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(0, 0, 0)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(1, -1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(2, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(3, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(4, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(5, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(6, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(7, -3, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(10, 0, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(11, 1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(11, 1, 1)));

    //Go backward
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(11, 1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(10, 0, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(7, -3, 3)));
    //Now walk forward a bit
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(7, -3, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(10, 0, 2)));
    //Now carry on back again
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(10, 0, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(7, -3, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(6, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(5, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(4, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(3, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(2, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(1, -1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(0, 0, 0)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(0, 0, 0)));

    bandIterator_destruct(bandIt);

    //Try creating an illegal expansion factor
    stTry {
            bandIterator_construct(anchorPairs, lX, lY, 1);
            CuAssertTrue(testCase, 0);
        }stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
            }stTryEnd

    //Cleanup
    bandIterator_destruct(bandIt);
    stList_destruct(anchorPairs);
}

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

static void test_symbol(CuTest *testCase) {
    return;
    Symbol cA[9] = { a, c, g, t, n, t, n, c, g };
    Symbol *cA2 = symbol_convertStringToSymbols("AcGTntNCG", 9);
    for (int32_t i = 0; i < 9; i++) {
        CuAssertTrue(testCase, cA[i] == cA2[i]);
    }
    free(cA2);
}

static void test_cell(CuTest *testCase) {

}

static void test_dpDiagonal(CuTest *testCase) {

}

static void test_dpMatrix(CuTest *testCase) {

}

static void test_diagonalDPCalculations(CuTest *testCase) {

}

static void test_getAlignedPairsWithBanding(CuTest *testCase) {
    return;
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

static void test_getBlastPairs(CuTest *testCase) {
    /*
     * Test the blast heuristic to get the different pairs.
     */
    return;
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
        for (int32_t i = 0; i < stList_length(blastPairs); i++) {
            stIntTuple *j = stList_get(blastPairs, i);
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

static void test_getSplitPoints(CuTest *testCase) {

}

CuSuite* pairwiseAlignmentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_diagonal);
    SUITE_ADD_TEST(suite, test_bandIterator);
    SUITE_ADD_TEST(suite, test_logAdd);
    SUITE_ADD_TEST(suite, test_symbol);
    SUITE_ADD_TEST(suite, test_cell);
    SUITE_ADD_TEST(suite, test_dpDiagonal);
    SUITE_ADD_TEST(suite, test_dpMatrix);
    SUITE_ADD_TEST(suite, test_diagonalDPCalculations);
    SUITE_ADD_TEST(suite, test_getAlignedPairsWithBanding);
    SUITE_ADD_TEST(suite, test_getBlastPairs);
    SUITE_ADD_TEST(suite, test_getSplitPoints);

    return suite;
}
