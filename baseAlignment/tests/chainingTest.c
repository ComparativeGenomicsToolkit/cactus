/*
 * chainingTest.c
 *
 *  Created on: 8 Mar 2012
 *      Author: benedictpaten
 */

#include "CuTest.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sonLib.h"
#include "chaining.h"

static APairArray *aPairArray = NULL;
static int32_t lX, lY;
static stList *pairs;

static void teardown() {
    if (aPairArray != NULL) {
        aPairArray_destruct(aPairArray);
        aPairArray = NULL;
        stList_destruct(pairs);
    }
}

static void setup() {
    teardown();
    lX = st_randomInt(0, 100);
    lY = st_randomInt(0, 100);
    pairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    double acceptProb = st_random();
    for (int32_t x = 0; x < lX; x++) {
        for (int32_t y = 0; y < lY; y++) {
            if (st_random() > acceptProb) {
                stList_append(pairs, stIntTuple_construct(2, x, y));
            }
        }
    }
}

int32_t max(int32_t i, int32_t j, int32_t k) {
    if (i < j) {
        i = j;
    }
    return i < k ? k : i;
}

static int32_t *initialiseMatrix() {
    int32_t *matrix = st_calloc(lX * lY, sizeof(int32_t));
    for (int32_t i = 0; i < stList_length(pairs); i++) {
        stIntTuple *pair = stList_get(pairs, i);
        int32_t x = stIntTuple_getPosition(pair, 0);
        int32_t y = stIntTuple_getPosition(pair, y);
        matrix[x * lY + y] = 1;
    }
    return matrix;
}

static int32_t *getScoreMatrixForward() {
    int32_t *matrix = initialiseMatrix();
    for (int32_t x = 1; x < lX; x++) {
        for (int32_t y = 1; y < lY; y++) {
            matrix[x * lY + y] += max(matrix[(x - 1) * lY + y - 1],
                    matrix[x * lY + y - 1], matrix[(x - 1) * lY + y]);
        }
    }
    return matrix;
}

static int32_t *getScoreMatrixBackward() {
    int32_t *matrix = initialiseMatrix();
    for (int32_t x = lX - 2; x >= 0; x--) {
        for (int32_t y = lY - 2; y >= 0; y--) {
            matrix[x * lY + y] += max(matrix[(x + 1) * lY + y + 1],
                    matrix[x * lY + y + 1], matrix[(x + 1) * lY + y]);
        }
    }
    return matrix;
}

static void test_aPairArray_calculateForwardScores(CuTest *testCase) {
    for (int32_t i = 0; i < 100; i++) {
        setup();
        int32_t *matrix = getScoreMatrixForward();
        aPairArray_calculateForwardScores(aPairArray);
        for (int32_t i = 0; i < aPairArray->length; i++) {
            APair pair = aPairArray->aPairs[i];
            CuAssertIntEquals(testCase, matrix[pair.x * lY + pair.y], pair.fScore);
        }
        free(matrix);
        teardown();
    }
}

static void test_aPairArray_calculateBackwardScores(CuTest *testCase) {
    for (int32_t i = 0; i < 100; i++) {
        setup();
        int32_t *matrix = getScoreMatrixBackward();
        aPairArray_calculateBackwardScores(aPairArray);
        for (int32_t i = 0; i < aPairArray->length; i++) {
            APair pair = aPairArray->aPairs[i];
            CuAssertIntEquals(testCase, matrix[pair.x * lY + pair.y], pair.bScore);
        }
        free(matrix);
        teardown();
    }
}

static void checkNonOverlapping(CuTest *testCase, stList *anchorPairs) {
    int32_t pX = -1;
    int32_t pY = -1;
    for (int32_t i = 0; i < stList_length(anchorPairs); i++) {
        stIntTuple *pair = stList_get(pairs, i);
        int32_t x = stIntTuple_getPosition(pair, 0);
        int32_t y = stIntTuple_getPosition(pair, 1);
        CuAssertTrue(testCase, x > pX);
        CuAssertTrue(testCase, y > pY);
        pX = x;
        pY = y;
    }
}

static void test_filterToRemoveOverlap(CuTest *testCase) {

}

static void test_getAnchorChain(CuTest *testCase) {
    for (int32_t i = 0; i < 100; i++) {
        setup();
        stList *anchorPairs = getAnchorChain(pairs, 0.8);
        checkNonOverlapping(testCase, anchorPairs);
        //Also check no pair can be slotted in the middle..
        stList_destruct(anchorPairs);
        teardown();
    }
}

CuSuite* chainingTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_aPairArray_calculateForwardScores);
    SUITE_ADD_TEST(suite, test_aPairArray_calculateBackwardScores);
    SUITE_ADD_TEST(suite, test_filterToRemoveOverlap);
    SUITE_ADD_TEST(suite, test_getAnchorChain);
    return suite;
}
