/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "multipleAligner.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "randomSequences.h"

static void test_diagonal(CuTest *testCase) {
    //Construct an example diagonal.
    int64_t xL = 10, yL = 20, xU = 30, yU = 0; //Coordinates of the upper and lower
    //pairs in x,y coordinates
    Diagonal d = diagonal_construct(xL + yL, xL - yL, xU - yU);
    CuAssertIntEquals(testCase, diagonal_getXay(d), xL + yL);
    CuAssertIntEquals(testCase, diagonal_getMinXmy(d), xL - yL);
    CuAssertIntEquals(testCase, diagonal_getMaxXmy(d), xU - yU);
    CuAssertIntEquals(testCase, diagonal_getWidth(d), (xU - yU - (xL - yL)) / 2 + 1);
    CuAssertIntEquals(testCase, diagonal_getXCoordinate(xL + yL, xL - yL), xL);
    CuAssertIntEquals(testCase, diagonal_getYCoordinate(xL + yL, xL - yL), yL);
    CuAssertTrue(testCase, diagonal_equals(d, d));
    CuAssertTrue(testCase, !diagonal_equals(d, diagonal_construct(0, 0, 0)));
    //A bogus diagonal is one such that |xay + xmy| % 2 != 0 or such that xmyR < xmyL.
    //Try constructing bogus diagonals, should throw exceptions.
    stTry
        {
            diagonal_construct(10, 5, 5);
            CuAssertTrue(testCase, 0);
        }
        stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
    stTry
        {
            diagonal_construct(10, 5, 5);
            CuAssertTrue(testCase, 0);
        }
        stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
    stTry
        {
            diagonal_construct(10, 6, 4);
            CuAssertTrue(testCase, 0);
        }
        stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
}

static bool testDiagonalsEqual(Diagonal d1, Diagonal d2) {
    bool b = diagonal_equals(d1, d2);
    if (!b) {
        st_logCritical("Diagonals not equal: d1: %s, d2: %s \n", diagonal_getString(d1), diagonal_getString(d2));
    }
    return b;
}

static void test_bands(CuTest *testCase) {
    stList *anchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    ///stList_append(anchorPairs, stIntTuple_construct2( 0, 0));
    stList_append(anchorPairs, stIntTuple_construct2(1, 0));
    stList_append(anchorPairs, stIntTuple_construct2(2, 1));
    stList_append(anchorPairs, stIntTuple_construct2(3, 3));
    /////stList_append(anchorPairs, stIntTuple_construct2( 5, 4));
    //Start the traversal
    int64_t lX = 6, lY = 5;
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

    //Cleanup
    bandIterator_destruct(bandIt);
    band_destruct(band);
    stList_destruct(anchorPairs);
}

static void test_logAdd(CuTest *testCase) {
    for (int64_t test = 0; test < 100000; test++) {
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
    for (int64_t i = 0; i < 9; i++) {
        CuAssertTrue(testCase, cA[i] == cA2[i]);
    }
    free(cA2);
}

static void test_cell(CuTest *testCase) {
    StateMachine *sM = stateMachine5_construct(fiveState);
    double lowerF[sM->stateNumber], middleF[sM->stateNumber], upperF[sM->stateNumber], currentF[sM->stateNumber];
    double lowerB[sM->stateNumber], middleB[sM->stateNumber], upperB[sM->stateNumber], currentB[sM->stateNumber];
    for (int64_t i = 0; i < sM->stateNumber; i++) {
        middleF[i] = sM->startStateProb(sM, i);
        middleB[i] = LOG_ZERO;
        lowerF[i] = LOG_ZERO;
        lowerB[i] = LOG_ZERO;
        upperF[i] = LOG_ZERO;
        upperB[i] = LOG_ZERO;
        currentF[i] = LOG_ZERO;
        currentB[i] = sM->endStateProb(sM, i);
    }
    Symbol cX = a, cY = t;
    //Do forward
    cell_calculateForward(sM, lowerF, NULL, NULL, middleF, cX, cY, NULL);
    cell_calculateForward(sM, upperF, middleF, NULL, NULL, cX, cY, NULL);
    cell_calculateForward(sM, currentF, lowerF, middleF, upperF, cX, cY, NULL);
    //Do backward
    cell_calculateBackward(sM, currentB, lowerB, middleB, upperB, cX, cY, NULL);
    cell_calculateBackward(sM, upperB, middleB, NULL, NULL, cX, cY, NULL);
    cell_calculateBackward(sM, lowerB, NULL, NULL, middleB, cX, cY, NULL);
    double totalProbForward = cell_dotProduct2(currentF, sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(middleB, sM, sM->startStateProb);
    st_logInfo("Total probability for cell test, forward %f and backward %f\n", totalProbForward, totalProbBackward);
    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.00001); //Check the forward and back probabilities are about equal
}

static void test_dpDiagonal(CuTest *testCase) {
    StateMachine *sM = stateMachine5_construct(fiveState);
    Diagonal diagonal = diagonal_construct(3, -1, 1);
    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal, sM->stateNumber);

    //Get cell
    double *c1 = dpDiagonal_getCell(dpDiagonal, -1);
    CuAssertTrue(testCase, c1 != NULL);
    double *c2 = dpDiagonal_getCell(dpDiagonal, 1);
    CuAssertTrue(testCase, c2 != NULL);
    CuAssertTrue(testCase, dpDiagonal_getCell(dpDiagonal, 3) == NULL);
    CuAssertTrue(testCase, dpDiagonal_getCell(dpDiagonal, -3) == NULL);

    dpDiagonal_initialiseValues(dpDiagonal, sM, sM->endStateProb); //Test initialise values
    double totalProb = LOG_ZERO;
    for (int64_t i = 0; i < sM->stateNumber; i++) {
        CuAssertDblEquals(testCase, c1[i], sM->endStateProb(sM, i), 0.0);
        CuAssertDblEquals(testCase, c2[i], sM->endStateProb(sM, i), 0.0);
        totalProb = logAdd(totalProb, 2 * c1[i]);
        totalProb = logAdd(totalProb, 2 * c2[i]);
    }

    DpDiagonal *dpDiagonal2 = dpDiagonal_clone(dpDiagonal);
    CuAssertTrue(testCase, dpDiagonal_equals(dpDiagonal, dpDiagonal2));

    CuAssertDblEquals(testCase, totalProb, dpDiagonal_dotProduct(dpDiagonal, dpDiagonal2), 0.001); //Check it runs

    dpDiagonal_destruct(dpDiagonal);
    dpDiagonal_destruct(dpDiagonal2);
}

static void test_dpMatrix(CuTest *testCase) {
    int64_t lX = 3, lY = 2;
    DpMatrix *dpMatrix = dpMatrix_construct(lX + lY, 5);

    CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), 0);

    for (int64_t i = -1; i <= lX + lY + 10; i++) {
        CuAssertTrue(testCase, dpMatrix_getDiagonal(dpMatrix, i) == NULL);
    }

    for (int64_t i = 0; i <= lX + lY; i++) {
        DpDiagonal *dpDiagonal = dpMatrix_createDiagonal(dpMatrix, diagonal_construct(i, -i, i));
        CuAssertTrue(testCase, dpDiagonal == dpMatrix_getDiagonal(dpMatrix, i));
        CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), i + 1);
    }

    for (int64_t i = lX + lY; i >= 0; i--) {
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

    const char *sX = "AGCG";
    const char *sY = "AGTTCG";
    int64_t lX = strlen(sX);
    int64_t lY = strlen(sY);
    SymbolString sX2 = symbolString_construct(sX, lX);
    SymbolString sY2 = symbolString_construct(sY, lY);
    StateMachine *sM = stateMachine5_construct(fiveState);
    DpMatrix *dpMatrixForward = dpMatrix_construct(lX + lY, sM->stateNumber);
    DpMatrix *dpMatrixBackward = dpMatrix_construct(lX + lY, sM->stateNumber);
    stList *anchorPairs = stList_construct();
    Band *band = band_construct(anchorPairs, lX, lY, 2);
    BandIterator *bandIt = bandIterator_construct(band);

    //Initialise matrices
    for (int64_t i = 0; i <= lX + lY; i++) {
        Diagonal d = bandIterator_getNext(bandIt);
        //initialisation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixBackward, d));
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(dpMatrixForward, d));
    }
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixForward, 0), sM, sM->startStateProb);
    dpDiagonal_initialiseValues(dpMatrix_getDiagonal(dpMatrixBackward, lX + lY), sM, sM->endStateProb);

    //Forward algorithm
    for (int64_t i = 1; i <= lX + lY; i++) {
        //Do the forward calculation
        diagonalCalculationForward(sM, i, dpMatrixForward, sX2, sY2);
    }

    //Backward algorithm
    for (int64_t i = lX + lY; i > 0; i--) {
        //Do the backward calculation
        diagonalCalculationBackward(sM, i, dpMatrixBackward, sX2, sY2);
    }

    //Calculate total probabilities
    double totalProbForward = cell_dotProduct2(
            dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixForward, lX + lY), lX - lY), sM, sM->endStateProb);
    double totalProbBackward = cell_dotProduct2(dpDiagonal_getCell(dpMatrix_getDiagonal(dpMatrixBackward, 0), 0), sM,
            sM->startStateProb);
    st_logInfo("Total forward and backward prob %f %f\n", (float) totalProbForward, (float) totalProbBackward);

    CuAssertDblEquals(testCase, totalProbForward, totalProbBackward, 0.001); //Check the forward and back probabilities are about equal

    //Test calculating the posterior probabilities along the diagonals of the matrix.
    for (int64_t i = 0; i <= lX + lY; i++) {
        //Calculate the total probs
        double totalDiagonalProb = diagonalCalculationTotalProbability(sM, i, dpMatrixForward, dpMatrixBackward, sX2,
                sY2);
        CuAssertDblEquals(testCase, totalProbForward, totalDiagonalProb, 0.01); //Check the forward and back probabilities are about equal
    }

    //Now do the posterior probabilities
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    void *extraArgs[1] = { alignedPairs };
    for (int64_t i = 1; i <= lX + lY; i++) {
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->threshold = 0.2;
        diagonalCalculationPosteriorMatchProbs(sM, i, dpMatrixForward, dpMatrixBackward, sX2, sY2, totalProbForward, p,
                extraArgs);
        pairwiseAlignmentBandingParameters_destruct(p);
    }

    stSortedSet *alignedPairsSet = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
            (void (*)(void *)) stIntTuple_destruct);
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(0, 0));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(1, 1));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(2, 4));
    stSortedSet_insert(alignedPairsSet, stIntTuple_construct2(3, 5));

    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *pair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(pair, 1), y = stIntTuple_get(pair, 2);
        st_logInfo("Pair %f %" PRIi64 " %" PRIi64 "\n", (float) stIntTuple_get(pair, 0) / PAIR_ALIGNMENT_PROB_1, x, y);
        CuAssertTrue(testCase, stSortedSet_search(alignedPairsSet, stIntTuple_construct2( x, y)) != NULL);
    }
    CuAssertIntEquals(testCase, stList_length(alignedPairs), 4);

}

stList *getRandomAnchorPairs(int64_t lX, int64_t lY) {
    stList *anchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    int64_t x = -1;
    int64_t y = -1;
    while (1) {
        x += st_randomInt(1, 20);
        y += st_randomInt(1, 20);
        if (x >= lX || y >= lY) {
            break;
        }
        assert(x >= 0 && x < lX);
        assert(y >= 0 && y < lY);
        stList_append(anchorPairs, stIntTuple_construct2(x, y));
    }
    return anchorPairs;
}

static void checkAlignedPairs(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY) {
    st_logInfo("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    stSortedSet *pairs = stSortedSet_construct3((int (*)(const void *, const void *)) stIntTuple_cmpFn,
            (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(blastPairs); i++) {
        stIntTuple *j = stList_get(blastPairs, i);
        CuAssertTrue(testCase, stIntTuple_length(j) == 3);

        int64_t x = stIntTuple_get(j, 1);
        int64_t y = stIntTuple_get(j, 2);
        int64_t score = stIntTuple_get(j, 0);
        CuAssertTrue(testCase, score > 0);
        CuAssertTrue(testCase, score <= PAIR_ALIGNMENT_PROB_1);

        CuAssertTrue(testCase, x >= 0);
        CuAssertTrue(testCase, y >= 0);
        CuAssertTrue(testCase, x < lX);
        CuAssertTrue(testCase, y < lY);

        //Check is unique
        stIntTuple *pair = stIntTuple_construct2(x, y);
        CuAssertTrue(testCase, stSortedSet_search(pairs, pair) == NULL);
        stSortedSet_insert(pairs, pair);
    }
    stSortedSet_destruct(pairs);
}

static void test_getAlignedPairsWithBanding(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(0, 100));
        char *sY = evolveSequence(sX); //stString_copy(seqX);
        int64_t lX = strlen(sX);
        int64_t lY = strlen(sY);
        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);
        SymbolString sX2 = symbolString_construct(sX, lX);
        SymbolString sY2 = symbolString_construct(sY, lY);
        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        p->traceBackDiagonals = st_randomInt(1, 10);
        p->minDiagsBetweenTraceBack = p->traceBackDiagonals + st_randomInt(2, 10);
        p->diagonalExpansion = st_randomInt(0, 10) * 2;
        StateMachine *sM = stateMachine5_construct(fiveState);
        stList *anchorPairs = getRandomAnchorPairs(lX, lY);

        stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
        void *extraArgs[1] = { alignedPairs };
        getPosteriorProbsWithBanding(sM, anchorPairs, sX2, sY2, p, 0, 0, diagonalCalculationPosteriorMatchProbs,
                extraArgs);
        //Check the aligned pairs.
        //Check the aligned pairs.
        checkAlignedPairs(testCase, alignedPairs, lX, lY);

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        free(sX2.sequence);
        free(sY2.sequence);
        stList_destruct(alignedPairs);
    }
}

static void checkBlastPairs(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY, bool checkNonOverlapping) {
    st_logInfo("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    int64_t pX = -1;
    int64_t pY = -1;
    for (int64_t i = 0; i < stList_length(blastPairs); i++) {
        stIntTuple *j = stList_get(blastPairs, i);
        CuAssertTrue(testCase, stIntTuple_length(j) == 2);

        int64_t x = stIntTuple_get(j, 0);
        int64_t y = stIntTuple_get(j, 1);

        CuAssertTrue(testCase, x >= 0);
        CuAssertTrue(testCase, y >= 0);
        CuAssertTrue(testCase, x < lX);
        CuAssertTrue(testCase, y < lY);
        if (checkNonOverlapping) {
            CuAssertTrue(testCase, x > pX);
            CuAssertTrue(testCase, y > pY);
        }
        pX = x;
        pY = y;
    }
}

static void test_getBlastPairs(CuTest *testCase) {
    /*
     * Test the blast heuristic to get the different pairs.
     */
    for (int64_t test = 0; test < 10; test++) {
        //Make a pair of sequences
        char *seqX = getRandomSequence(st_randomInt(0, 10000));
        char *seqY = evolveSequence(seqX); //stString_copy(seqX);
        int64_t lX = strlen(seqX), lY = strlen(seqY);
        st_logInfo("Sequence X to align: %s END, seq length %" PRIi64 "\n", seqX, lX);
        st_logInfo("Sequence Y to align: %s END, seq length %" PRIi64 "\n", seqY, lY);

        int64_t trim = st_randomInt(0, 5);
        bool repeatMask = st_random() > 0.5;
        st_logInfo("Using random trim %" PRIi64 ", recursive %" PRIi64 " \n", trim, repeatMask);

        stList *blastPairs = getBlastPairs(seqX, seqY, lX, lY, trim, repeatMask);

        checkBlastPairs(testCase, blastPairs, lX, lY, 0);
        stList_destruct(blastPairs);
        free(seqX);
        free(seqY);
    }
}

static void test_filterToRemoveOverlap(CuTest *testCase) {
    for (int64_t i = 0; i < 100; i++) {
        //Make random pairs
        int64_t lX = st_randomInt(0, 100);
        int64_t lY = st_randomInt(0, 100);
        stList *pairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
        double acceptProb = st_random();
        for (int64_t x = 0; x < lX; x++) {
            for (int64_t y = 0; y < lY; y++) {
                if (st_random() > acceptProb) {
                    stList_append(pairs, stIntTuple_construct2(x, y));
                }
            }
        }
        //Now run filter pairs
        stList *nonoverlappingPairs = filterToRemoveOverlap(pairs);

        //Check non overlapping
        checkBlastPairs(testCase, nonoverlappingPairs, lX, lY, 1);

        //Now check maximal
        stList *nonoverlappingPairs2 = stList_construct();
        for (int64_t i = 0; i < stList_length(pairs); i++) {
            stIntTuple *pair = stList_get(pairs, i);
            int64_t x = stIntTuple_get(pair, 0);
            int64_t y = stIntTuple_get(pair, 1);
            bool nonOverlapping = 1;
            for (int64_t j = 0; j < stList_length(pairs); j++) {
                stIntTuple *pair2 = stList_get(pairs, j);
                int64_t x2 = stIntTuple_get(pair2, 0);
                int64_t y2 = stIntTuple_get(pair2, 1);
                if ((x2 <= x && y2 >= y) || (x2 >= x && y2 <= y)) {
                    nonOverlapping = 0;
                    break;
                }
            }
            if (nonOverlapping) {
                stList_append(nonoverlappingPairs2, pair);
            }
        }
        stSortedSet *nonOverlappingPairsSet = stList_getSortedSet(nonoverlappingPairs,
                (int (*)(const void *, const void *)) stIntTuple_cmpFn);
        stSortedSet *nonOverlappingPairsSet2 = stList_getSortedSet(nonoverlappingPairs2,
                (int (*)(const void *, const void *)) stIntTuple_cmpFn);
        st_logDebug("The non-overlapping set sizes are %" PRIi64 " %" PRIi64 "\n",
                stSortedSet_size(nonOverlappingPairsSet), stSortedSet_size(nonOverlappingPairsSet2));
        CuAssertTrue(testCase, stSortedSet_equals(nonOverlappingPairsSet, nonOverlappingPairsSet2));

        //Cleanup
        stSortedSet_destruct(nonOverlappingPairsSet);
        stSortedSet_destruct(nonOverlappingPairsSet2);
        stList_destruct(nonoverlappingPairs2);
        stList_destruct(pairs);
        stList_destruct(nonoverlappingPairs);

    }
}

static void test_getBlastPairsWithRecursion(CuTest *testCase) {
    /*
     * Test the blast heuristic to get the different pairs.
     */
    for (int64_t test = 0; test < 10; test++) {
        //Make a pair of sequences
        char *seqX = getRandomSequence(st_randomInt(0, 10000));
        char *seqY = evolveSequence(seqX); //stString_copy(seqX);
        int64_t lX = strlen(seqX), lY = strlen(seqY);
        st_logInfo("Sequence X to align: %s END, seq length %" PRIi64 "\n", seqX, lX);
        st_logInfo("Sequence Y to align: %s END, seq length %" PRIi64 "\n", seqY, lY);

        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

        stList *blastPairs = getBlastPairsForPairwiseAlignmentParameters(seqX, seqY, lX, lY, p);

        checkBlastPairs(testCase, blastPairs, lX, lY, 1);
        stList_destruct(blastPairs);
        free(seqX);
        free(seqY);
    }
}

static void test_getSplitPoints(CuTest *testCase) {
    int64_t matrixSize = 2000 * 2000;

    stList *anchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    //Test a small region, which produces no splits
    int64_t lX = 3000;
    int64_t lY = 1000;
    stList *splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 0);
    CuAssertIntEquals(testCase, 1, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, lX, lY)));
    stList_destruct(splitPoints);

    //Test with one really big matrix with no anchors
    lX = 20000;
    lY = 25000;
    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 1, 1);
    CuAssertIntEquals(testCase, 0, stList_length(splitPoints));
    stList_destruct(splitPoints);

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 1, 0);
    CuAssertIntEquals(testCase, 1, stList_length(splitPoints));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(18000, 23000, lX, lY)));
    stList_destruct(splitPoints);

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 1);
    CuAssertIntEquals(testCase, 1, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, 2000, 2000)));
    stList_destruct(splitPoints);

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 0);
    CuAssertIntEquals(testCase, 2, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, 2000, 2000)));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 1), stIntTuple_construct4(18000, 23000, lX, lY)));
    stList_destruct(splitPoints);

    //Now test with some more points
    stList_append(anchorPairs, stIntTuple_construct2(2000, 2000)); //This should not create a split
    stList_append(anchorPairs, stIntTuple_construct2(4002, 4001)); //This should cause a split
    stList_append(anchorPairs, stIntTuple_construct2(5000, 5000)); //This should not cause a split
    stList_append(anchorPairs, stIntTuple_construct2(8000, 6000)); //Neither should this (it is maximum sized)
    stList_append(anchorPairs, stIntTuple_construct2(9000, 9000)); //Or this (it is maximum sized)
    stList_append(anchorPairs, stIntTuple_construct2(10000, 14000)); //This should create a split
    stList_append(anchorPairs, stIntTuple_construct2(15000, 15000)); //This should also create a split
    stList_append(anchorPairs, stIntTuple_construct2(16000, 16000)); //This should not, but there will be a split with the end.

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 0);

    for (int64_t i = 0; i < stList_length(splitPoints); i++) {
        stIntTuple *j = stList_get(splitPoints, i);
        st_logInfo("I got split point: x1: %" PRIi64 " y1: %" PRIi64 " x2: %" PRIi64 " y2: %" PRIi64 "\n",
                stIntTuple_get(j, 0), stIntTuple_get(j, 1), stIntTuple_get(j, 2), stIntTuple_get(j, 3));
    }

    CuAssertIntEquals(testCase, 5, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, 3001, 3001)));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 1), stIntTuple_construct4(3002, 3001, 9500, 11001)));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 2), stIntTuple_construct4(9501, 12000, 12001, 14500)));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 3), stIntTuple_construct4(13000, 14501, 18000, 18001)));
    CuAssertTrue(testCase,
            stIntTuple_equalsFn(stList_get(splitPoints, 4), stIntTuple_construct4(18001, 23000, 20000, 25000)));

    stList_destruct(splitPoints);
    stList_destruct(anchorPairs);
}

static void test_getAlignedPairs(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(0, 100));
        char *sY = evolveSequence(sX); //stString_copy(seqX);
        int64_t lX = strlen(sX);
        int64_t lY = strlen(sY);
        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        StateMachine *sM = stateMachine5_construct(fiveState);

        stList *alignedPairs = getAlignedPairs(sM, sX, sY, p, 0, 0);

        //Check the aligned pairs.
        checkAlignedPairs(testCase, alignedPairs, lX, lY);

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        stList_destruct(alignedPairs);
    }
}

static void test_getAlignedPairsWithRaggedEnds(CuTest *testCase) {
    for (int64_t test = 0; test < 1000; test++) {
        //Make a pair of sequences
        int64_t coreLength = 100, randomPortionLength = 100;
        char *sX = getRandomSequence(coreLength);
        char *sY = stString_print("%s%s%s", getRandomSequence(randomPortionLength), sX,
                getRandomSequence(randomPortionLength)); //x with an extra bit at the end.

        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();
        StateMachine *sM = stateMachine5_construct(fiveState);
        stList *alignedPairs = getAlignedPairs(sM, sX, sY, p, 1, 1);
        alignedPairs = filterPairwiseAlignmentToMakePairsOrdered(alignedPairs, sX, sY, 0.2);

        //Check the aligned pairs.
        checkAlignedPairs(testCase, alignedPairs, strlen(sX), strlen(sY));
        CuAssertIntEquals(testCase, stList_length(alignedPairs), coreLength);
        for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
            stIntTuple *j = stList_get(alignedPairs, i);
            CuAssertTrue(testCase, stIntTuple_length(j) == 3);

            int64_t x = stIntTuple_get(j, 1);
            int64_t y = stIntTuple_get(j, 2);
            CuAssertIntEquals(testCase, x + randomPortionLength, y);
        }

        //Cleanup
        stateMachine_destruct(sM);
        free(sX);
        free(sY);
        stList_destruct(alignedPairs);
    }
}

/*
 * EM training tests.
 */

static void test_hmm(CuTest *testCase, StateMachineType stateMachineType) {
    //Expectation object
    Hmm *hmm = hmm_constructEmpty(0.0, stateMachineType);

    //Add some transition expectations
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            hmm_addToTransitionExpectation(hmm, from, to, from * hmm->stateNumber + to);
        }
    }

    //Add some emission expectations
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < SYMBOL_NUMBER_NO_N; x++) {
            for (int64_t y = 0; y < SYMBOL_NUMBER_NO_N; y++) {
                hmm_addToEmissionsExpectation(hmm, state, x, y,
                        state * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N + x * SYMBOL_NUMBER_NO_N + y);
            }
        }
    }

    //Write to a file
    char *tempFile = stString_print("./temp%" PRIi64 ".hmm", st_randomInt(0, INT64_MAX));
    CuAssertTrue(testCase, !stFile_exists(tempFile)); //Quick check that we don't write over anything.
    FILE *fH = fopen(tempFile, "w");
    hmm_write(hmm, fH);
    fclose(fH);
    hmm_destruct(hmm);

    //Load from a file
    hmm = hmm_loadFromFile(tempFile);
    stFile_rmrf(tempFile);

    //Check the transition expectations
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            CuAssertTrue(testCase, hmm_getTransition(hmm, from, to) == from * hmm->stateNumber + to);
        }
    }

    //Check the emission expectations
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < SYMBOL_NUMBER_NO_N; x++) {
            for (int64_t y = 0; y < SYMBOL_NUMBER_NO_N; y++) {
                CuAssertTrue(testCase,
                        hmm_getEmissionsExpectation(hmm, state, x, y) == state * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N + x * SYMBOL_NUMBER_NO_N + y);
            }
        }
    }

    //Normalise
    hmm_normalise(hmm);

    //Recheck transitions
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            double z = from * hmm->stateNumber * hmm->stateNumber + (hmm->stateNumber * (hmm->stateNumber - 1)) / 2;
            CuAssertDblEquals(testCase, (from * hmm->stateNumber + to) / z, hmm_getTransition(hmm, from, to), 0.0);
        }
    }

    //Recheck the emissions
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < SYMBOL_NUMBER_NO_N; x++) {
            for (int64_t y = 0; y < SYMBOL_NUMBER_NO_N; y++) {
                double z = SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N * state
                        + ((SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N) * ((SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N) - 1))
                                / 2;
                CuAssertTrue(testCase,
                        hmm_getEmissionsExpectation(hmm, state, x, y) == (state * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N + x * SYMBOL_NUMBER_NO_N + y)/z);
            }
        }
    }

    //Clean up
    hmm_destruct(hmm);
}

static void test_hmm_5State(CuTest *testCase) {
    test_hmm(testCase, fiveState);
}

static void test_hmm_5StateAsymmetric(CuTest *testCase) {
    test_hmm(testCase, fiveStateAsymmetric);
}

static void test_hmm_3State(CuTest *testCase) {
    test_hmm(testCase, threeState);
}

static void test_hmm_3StateAsymmetric(CuTest *testCase) {
    test_hmm(testCase, threeStateAsymmetric);
}

static void test_em(CuTest *testCase, StateMachineType stateMachineType) {
    for (int64_t test = 0; test < 100; test++) {
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(10, 100));
        char *sY = evolveSequence(sX); //stString_copy(seqX);
        st_logInfo("Sequence X to align: %s END\n", sX);
        st_logInfo("Sequence Y to align: %s END\n", sY);

        //Now do alignment
        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

        //Currently starts from random model and iterates.
        double pLikelihood = -INFINITY;
        Hmm *hmm = hmm_constructEmpty(0.0, stateMachineType);
        hmm_randomise(hmm);
        StateMachine *sM = hmm_getStateMachine(hmm);
        hmm_destruct(hmm);

        for (int64_t iteration = 0; iteration < 10; iteration++) {
            hmm = hmm_constructEmpty(0.000000000001, stateMachineType); //The tiny pseudo count prevents overflow
            getExpectations(sM, hmm, sX, sY, p, 0, 0);
            hmm_normalise(hmm);
            //Log stuff
            for (int64_t from = 0; from < sM->stateNumber; from++) {
                for (int64_t to = 0; to < sM->stateNumber; to++) {
                    st_logInfo("Transition from %" PRIi64 " to %" PRIi64 " has expectation %f\n", from, to,
                            hmm_getTransition(hmm, from, to));
                }
            }
            for (int64_t x = 0; x < SYMBOL_NUMBER_NO_N; x++) {
                for (int64_t y = 0; y < SYMBOL_NUMBER_NO_N; y++) {
                    st_logInfo("Emission x %" PRIi64 " y %" PRIi64 " has expectation %f\n", x, y,
                            hmm_getEmissionsExpectation(hmm, sM->matchState, x, y));
                }
            }

            st_logInfo("->->-> Got expected likelihood %f for trial %" PRIi64 " and  iteration %" PRIi64 "\n",
                    hmm->likelihood, test, iteration);
            assert(pLikelihood <= hmm->likelihood * 0.95);
            CuAssertTrue(testCase, pLikelihood <= hmm->likelihood * 0.95);
            pLikelihood = hmm->likelihood;
            stateMachine_destruct(sM);
            sM = hmm_getStateMachine(hmm);
            hmm_destruct(hmm);
        }

        //Cleanup
        pairwiseAlignmentBandingParameters_destruct(p);
        free(sX);
        free(sY);
    }
}

static void test_em_5State(CuTest *testCase) {
    test_em(testCase, fiveState);
}

static void test_em_3StateAsymmetric(CuTest *testCase) {
    test_em(testCase, threeStateAsymmetric);
}

static void test_em_3State(CuTest *testCase) {
    test_em(testCase, threeState);
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
    SUITE_ADD_TEST(suite, test_getBlastPairsWithRecursion);
    SUITE_ADD_TEST(suite, test_filterToRemoveOverlap);
    SUITE_ADD_TEST(suite, test_getSplitPoints);
    SUITE_ADD_TEST(suite, test_getAlignedPairs);
    SUITE_ADD_TEST(suite, test_getAlignedPairsWithRaggedEnds);
    SUITE_ADD_TEST(suite, test_hmm_5State);
    SUITE_ADD_TEST(suite, test_hmm_5StateAsymmetric);
    SUITE_ADD_TEST(suite, test_hmm_3State);
    SUITE_ADD_TEST(suite, test_hmm_3StateAsymmetric);
    SUITE_ADD_TEST(suite, test_em_3State);
    SUITE_ADD_TEST(suite, test_em_3StateAsymmetric);
    SUITE_ADD_TEST(suite, test_em_5State);

    return suite;
}
