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
    CuAssertIntEquals(testCase, diagonal_getWidth(d), (xU - yU - (xL - yL))/2 + 1);
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
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
    stTry {
            diagonal_construct(10, 5, 5);
            CuAssertTrue(testCase, 0);
        }stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
    stTry {
            diagonal_construct(10, 6, 4);
            CuAssertTrue(testCase, 0);
        }stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
}

static bool testDiagonalsEqual(Diagonal d1, Diagonal d2) {
    bool b = diagonal_equals(d1, d2);
    if (!b) {
        st_logCritical("Diagonals not equal: d1: %s, d2: %s \n",
                diagonal_getString(d1), diagonal_getString(d2));
    }
    return b;
}

static void test_bands(CuTest *testCase) {
    stList *anchorPairs = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    ///stList_append(anchorPairs, stIntTuple_construct(2, 0, 0));
    stList_append(anchorPairs, stIntTuple_construct(2, 1, 0));
    stList_append(anchorPairs, stIntTuple_construct(2, 2, 1));
    stList_append(anchorPairs, stIntTuple_construct(2, 3, 3));
    /////stList_append(anchorPairs, stIntTuple_construct(2, 5, 4));
    //Start the traversal
    int32_t lX = 6, lY = 5;
    Band *band = band_construct(anchorPairs, lX, lY, 2);
    BandIterator *bandIt = bandIterator_construct(band);

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

    //Now forward again
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(0, 0, 0)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(1, -1, 1)));
    //Now back again
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(1, -1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(0, 0, 0)));

    //Try creating an illegal expansion factor
    stTry {
            band_construct(anchorPairs, lX, lY, 1);
            CuAssertTrue(testCase, 0);
        }stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd

    //Try creating an illegal band in which pairs not all increasing on x+y.
    stTry {
            stList_append(anchorPairs, stIntTuple_construct(2, 4, 3));
            band_construct(anchorPairs, lX, lY, 2);
            CuAssertTrue(testCase, 0);
        }stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd

    //Cleanup
    bandIterator_destruct(bandIt);
    band_destruct(band);
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
    Symbol cA[9] = { a, c, g, t, n, t, n, c, g };
    Symbol *cA2 = symbol_convertStringToSymbols("AcGTntNCG", 9);
    for (int32_t i = 0; i < 9; i++) {
        CuAssertTrue(testCase, cA[i] == cA2[i]);
    }
    free(cA2);
}

static void test_cell(CuTest *testCase) {
    double lowerF[STATE_NUMBER], middleF[STATE_NUMBER], upperF[STATE_NUMBER],
            currentF[STATE_NUMBER];
    double lowerB[STATE_NUMBER], middleB[STATE_NUMBER], upperB[STATE_NUMBER],
            currentB[STATE_NUMBER];
    for (int32_t i = 0; i < STATE_NUMBER; i++) {
        middleF[i] = state_startStateProb(i);
        middleB[i] = LOG_ZERO;
        lowerF[i] = LOG_ZERO;
        lowerB[i] = LOG_ZERO;
        upperF[i] = LOG_ZERO;
        upperB[i] = LOG_ZERO;
        currentF[i] = LOG_ZERO;
        currentB[i] = state_endStateProb(i);
    }
    Symbol cX = a, cY = t;
    //Do forward
    cell_calculateForward(lowerF, NULL, NULL, middleF, cX, cY);
    cell_calculateForward(upperF, middleF, NULL, NULL, cX, cY);
    cell_calculateForward(currentF, lowerF, middleF, upperF, cX, cY);
    //Do backward
    cell_calculateBackward(currentB, lowerB, middleB, upperB, cX, cY);
    cell_calculateBackward(upperB, middleB, NULL, NULL, cX, cY);
    cell_calculateBackward(lowerB, NULL, NULL, middleB, cX, cY);
    double totalProbForward = cell_dotProduct2(currentF, state_endStateProb);
    double totalProbBackward = cell_dotProduct2(middleB, state_startStateProb);

    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.001); //Check the forward and back probabilities are about equal
}

static void test_dpDiagonal(CuTest *testCase) {
    Diagonal diagonal = diagonal_construct(3, -1, 1);
    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal);

    //Get cell
    double *c1 = dpDiagonal_getCell(dpDiagonal, -1);
    CuAssertTrue(testCase, c1 != NULL);
    double *c2 = dpDiagonal_getCell(dpDiagonal, 1);
    CuAssertTrue(testCase, c2 != NULL);
    CuAssertTrue(testCase,dpDiagonal_getCell(dpDiagonal, 3) == NULL);
    CuAssertTrue(testCase,dpDiagonal_getCell(dpDiagonal, -3) == NULL);

    dpDiagonal_initialiseValues(dpDiagonal, state_startStateProb); //Test initialise values
    for (int32_t i = 0; i < STATE_NUMBER; i++) {
        CuAssertDblEquals(testCase, c1[i], state_startStateProb(i), 0.0);
        CuAssertDblEquals(testCase, c2[i], state_startStateProb(i), 0.0);
    }

    DpDiagonal *dpDiagonal2 = dpDiagonal_clone(dpDiagonal);
    CuAssertTrue(testCase, dpDiagonal_equals(dpDiagonal, dpDiagonal2));

    dpDiagonal_destruct(dpDiagonal);
}

static void test_dpMatrix(CuTest *testCase) {
    int32_t lX = 3, lY = 2;
    DpMatrix *dpMatrix = dpMatrix_construct(lX, lY);

    CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), 0);

    for(int32_t i=-1; i<= lX + lY + 10; i++) {
        CuAssertTrue(testCase, dpMatrix_getDiagonal(dpMatrix, i) == NULL);
    }

    for(int32_t i=0; i<= lX + lY; i++) {
        DpDiagonal *dpDiagonal = dpMatrix_createDiagonal(dpMatrix, diagonal_construct(i, -i, i));
        CuAssertTrue(testCase, dpDiagonal == dpMatrix_getDiagonal(dpMatrix, i));
        CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), i+1);
    }

    for(int32_t i=lX + lY; i>= 0; i--) {
        dpMatrix_deleteDiagonal(dpMatrix, i);
        CuAssertTrue(testCase, dpMatrix_getDiagonal(dpMatrix, i) == NULL);
        CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), i);
    }

    CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), 0);

    dpMatrix_destruct(dpMatrix);
}

static void test_diagonalDPCalculations(CuTest *testCase) {
    //Sets up a complete matrix for the following example and checks the total marginal
    //probability and the posterior probabilities of the matches

    const char *sX = "ACTGTT";
    const char *sY = "ATGAATT";
    int32_t lX = strlen(sX);
    int32_t lY = strlen(sY);
    Symbol *sX2 = symbol_convertStringToSymbols(sX, lX);
    Symbol *sY2 = symbol_convertStringToSymbols(sY, lY);
    DpMatrix *dpMatrixForward = dpMatrix_construct(lX, lY);
    DpMatrix *dpMatrixBackward = dpMatrix_construct(lX, lY);
    stList *anchorPairs = stList_construct();
    Band *band = band_construct(anchorPairs, lX, lY, 2);
    BandIterator *bandIt = bandIterator_construct(band);

    //Initialise matrices
    for(int32_t i=0; i<=lX+lY; i++) {
        Diagonal d = bandIterator_getNext(bandIt);
        //intitialisation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixBackward, d));
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixForward, d));
    }
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixForward, 0), state_startStateProb);
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixBackward, lX + lY), state_endStateProb);

    //Forward algorithm
    for(int32_t i=1; i<=lX+lY; i++) {
        //Do the forward calculation
        diagonalCalculationForward(i, dpMatrixForward, sX2, sY2, lX, lY);
    }

    //Backward algorithm
    for(int32_t i=lX+lY; i>0; i--) {
        //Do the backward calculation
        diagonalCalculationBackward(i, dpMatrixBackward, sX2, sY2, lX, lY);
    }

    //Calculate total probabilities
    double totalProbForward = cell_dotProduct2(dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixForward, lX+lY), lX-lY), state_endStateProb);
    double totalProbBackward = cell_dotProduct2(dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixBackward, 0), 0), state_startStateProb);

    st_uglyf("The total probs %f %f\n", (float)totalProbForward, (float)totalProbBackward);
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.01); //Check the forward and back probabilities are about equal

    //Now do the posterior probabilities
    stList *alignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int32_t i=1; i<=lX+lY; i++) {
        diagonalCalculationPosterior(i, dpMatrixForward, dpMatrixBackward, sX2, sY2, lX, lY, 0.01, totalProbForward, alignedPairs);
    }
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *pair = stList_get(alignedPairs, i);
        st_uglyf("Pair %f %i %i\n", (float)stIntTuple_getPosition(pair, 0)/PAIR_ALIGNMENT_PROB_1, stIntTuple_getPosition(pair, 1), stIntTuple_getPosition(pair, 2));
    }
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
        PairwiseAlignmentParameters *p =
                pairwiseAlignmentBandingParameters_construct();
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
        st_logInfo("Sequence X to align: %s END, seq length %i\n", seqX,
                seqXLength);
        st_logInfo("Sequence Y to align: %s END, seq length %i\n", seqY,
                seqYLength);

        int32_t trim = st_randomInt(0, 5);
        bool recursive = st_random() > 0.5;
        st_logInfo("Using random trim %i, recursive %i \n", trim, recursive);

        stList *blastPairs = getBlastPairs(seqX, seqY, seqXLength, seqYLength,
                trim, recursive);

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
    SUITE_ADD_TEST(suite, test_bands);
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
