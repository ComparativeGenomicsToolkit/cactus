/*
 * diagonalIterator.c
 *
 *  Created on: 1 Mar 2012
 *      Author: benedictpaten
 */

#include "sonLib.h"

typedef struct diagonal {
    int32_t xay;
    int32_t xmyL;
    int32_t xmyR;
} Diagonal;

Diagonal diagonal_construct(int32_t xay, int32_t xmyL, int32_t xmyR) {
    Diagonal diagonal;
    diagonal.xay = xay;
    diagonal.xmyL = xmyL;
    diagonal.xmyR = xmyR;
    assert(xmyL <= xmyR);
    assert(xay >= 0);
    return diagonal;
}

inline int32_t diagonal_getXay(Diagonal diagonal) {
    return diagonal.xay;
}

inline int32_t diagonal_getMinXmy(Diagonal diagonal) {
    return diagonal.xmyL;
}

inline int32_t diagonal_getMaxXmy(Diagonal diagonal) {
    return diagonal.xmyR;
}

inline int32_t diagonal_getWidth(Diagonal diagonal) {
    return digonal.xmlR - diagonal.xmyL;
}

typedef struct _dpDiagonal {
    Diagonal diagonal;
    double *cells;
    int32_t stateNumber;
} DpDiagonal;

DpDiagonal *dpDiagonal_construct(Diagonal diagonal, int32_t stateNumber) {
    DpDiagonal *dpDiagonal = st_malloc(sizeof(DpDiagonal));
    dpDiagonal->diagonal = diagonal;
    dpDiagonal->stateNumber = stateNumber;
    dpDiagonal->cells = st_malloc(sizeof(double) * stateNumber * diagonal_getWidth(diagonal));
    return dpDiagonal;
}

double *dpDiagonal_getCell(DpDiagonal *dpDiagonal, int32_t xmy) {
    if(dpDiagonal->diagonal.xmyL < xmy || dpDiagonal->diagonal.xmyR > xmy) {
        return NULL;
    }
    return &dpDiagonal->cells[(xmy - dpDiagonal->diagonal.xmyL) * dpDiagonal->stateNumber];
}

void dpDiagonal_initialiseValues(DpDiagonal *dpDiagonal, double *stateValues) {
    int32_t width = diagonal_getWidth(dpDiagonal->diagonal);
    for(int32_t i=0; i<width; i++) {
        memcpy(dpDiagonal->cells + (width * dpDiagonal->stateNumber), stateValues, sizeof(double)*dpDiagonal->stateNumber);
    }
}

typedef struct _dpMatrix {
    DPDiagonal *diagonals;
    int32_t diagonalNumber;
    int32_t stateNumber;
    int32_t activeDiagonals;
    int64_t activeCells;
} DpMatrix;

DpMatrix *dpMatrix_construct(int32_t lX, int32_t lY, int32_t stateNumber) {
    DpMatrix *dpMatrix = st_malloc(sizeof(DpMatrix));
    dpMatrix->diagonalNumber = lX + lY;
    dpMatrix->stateNumber = stateNumber;
    dpMatrix->diagonals = st_calloc(dpMatrix->diagonalNumber, sizeof(DpDiagonal *));
    dpMatrix->activeDiagonals = 0;
    dpMatrix->activeCells = 0;
    return dpMatrix;
}

DpDiagonal dpMatrix_getDiagonal(DpMatrix *dpMatrix, int32_t xay) {
    if(xay < 0 || xay >= dpMatrix->diagonalNumber) {
        return NULL;
    }
    return dpMatrix->diagonals[xay];
}

int32_t dpMatrix_getActiveDiagonalNumber(DpMatrix *dpMatrix) {
    return dpMatrix->activeDiagonals;
}

void dpMatrix_createDiagonal(DpMatrix *dpMatrix, Diagonal diagonal) {
    assert(diagonal.xay >= 0);
    assert(diagonal.xay < dpMatrix->diagonalNumber);
    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal);
    dpMatrix->activeDiagonals[diagonal_getXay(diagonal)] = dpDiagonal;
    dpMatrix->activeDiagonals++;
    dpMatrix->activeCells += diagonal_getWidth(dpDiagonal->diagonal);
}

void dpMatrix_deleteDiagonal(DpMatrix *dpMatrix, int32_t xay) {
    assert(diagonal.xay >= 0);
    assert(diagonal.xay < dpMatrix->diagonalNumber);
    if(dpMatrix->activeDiagonals[xay] != NULL) {
        dpMatrix->activeDiagonals--;
        assert(dpMatrix->activeDiagonals >= 0);
        DpDiagonal *dpDiagonal = dpMatrix->activeDiagonals[xay];
        dpMatrix->activeCells -= diagonal_getWidth(dpDiagonal->diagonal);
        assert(dpMatrix->activeCells >= 0);
        dpDiagonal_destruct(dpMatrix->activeDiagonals[xay]);
        dpMatrix->activeDiagonals[xay] = NULL;
    }
}

int64_t dpMatrix_getTotalActiveCells(DpMatrix *dpMatrix) {
    return dpMatrix->activeCells;
}

typedef struct _bandIterator {
    stList *anchorPairs;
    stListIterator *anchorPairsIt;
    int32_t xL, xR, yL, yR, cDiag;
    int64_t totalActiveCells;
} BandIterator;

BandIterator *bandIterator_construct(stList *anchorPairs, int32_t lX, int32_t lY) {

}

BandIterator *bandIterator_clone(BandIterator *bandIterator) {

}

void bandIterator_destruct(BandIterator *bandIterator) {

}

Diagonal bandIterator_getNext(BandIterator *bandIterator);

Diagonal bandIterator_getPrevious(BandIterator *bandIterator);

void calculateForwardCell(double *current, double *lower, double *middle, double *upper) {

}

void doForwardCalculation(Diagonal diagonal, DpMatrix *dpMatrix) {
    DpDiagonal *dpDiagonal = dpMatrix_getDiagonal(dpMatrix, diagonal_getXay(diagonal));
    assert(dpDiagonal != NULL);
    DpDiagonal *dPDiagonalM1 = dpMatrix_getDiagonal(dpMatrix, diagonal_getXay(diagonal) - 1);
    assert(dpDiagonalM1 != NULL);
    DpDiagonal *dPDiagonalM2 = dpMatrix_getDiagonal(dpMatrix, diagonal_getXay(diagonal) - 2);

    int32_t xmy = diagonal_getMinXmy(diagonal) - 1;
    while (++xmy <= diagonal_getMaxXmy(diagonal)) {
        double *current = dpDiagonal_getCell(dPDiagonal, xmy);
        double *lower = dpDiagonal_getCell(dPDiagonalM1, xmy - 1);
        double *middle = dPDiagonalM2 == NULL ? NULL : dpDiagonal_getCell(ppdPDiagonalM2, xmy);
        double *upper = dpDiagonal_getCell(pdPDiagonalM1, xmy + 1);
        calculateForwardCell(current, lower, middle, upper, extraArg);
    }
}

void calculateBackwardCell(double *current, double *lower, double *middle, double *upper) {

}

void doBackwardCalculation(Diagonal diagonal, DpMatrix *dpMatrix) {
    DpDiagonal *dpDiagonal = dpMatrix_getDiagonal(dpMatrix, diagonal_getXay(diagonal));
    assert(dpDiagonal != NULL);
    DpDiagonal *dPDiagonalM1 = dpMatrix_getDiagonal(dpMatrix, diagonal_getXay(diagonal) + 1);
    assert(dpDiagonalM1 != NULL);
    DpDiagonal *dPDiagonalM2 = dpMatrix_getDiagonal(dpMatrix, diagonal_getXay(diagonal) + 2);

    int32_t xmy = diagonal_getMinXmy(diagonal) - 1;
    while (++xmy <= diagonal_getMaxXmy(diagonal)) {
        double *current = dpDiagonal_getCell(dPDiagonal, xmy);
        double *lower = dpDiagonal_getCell(dPDiagonalM1, xmy - 1);
        double *middle = dPDiagonalM2 == NULL ? NULL : dpDiagonal_getCell(ppdPDiagonalM2, xmy);
        double *upper = dpDiagonal_getCell(pdPDiagonalM1, xmy + 1);
        calculateBackwardCell(current, lower, middle, upper, extraArg);
    }
}

double calculateTotalProbability(double *current, double *lower, double *middle, double *upper) {
    current = copyCell(current);
    calculateForwardCell(current, lower, middle, upper);
    double totalProbability = sumCell(current);
    free(current);
    return totalProbability;
}

void calculatePosteriorCell(double *current, double *lower, double *middle, double *upper, double totalProbability,
        stList *alignedPairs) {
    current = copyCell(current);
    calculateForwardCell(current, lower, middle, upper);
    double posteriorProbability = current[0] - totalProbability;
}

void doPosteriorCalculation(Diagonal *diagonal, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
        double threshold, stList *alignedPairs) {
    DpDiagonal *dpDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, diagonal_getXay(diagonal));
    assert(dpDiagonal != NULL);
    DpDiagonal *dPDiagonalM1 = dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal) - 1);
    assert(dpDiagonalM1 != NULL);
    DpDiagonal *dPDiagonalM2 = dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal) - 2);

    double totalProbability = LOG_ZERO;
    int32_t xmy = diagonal_getMinXmy(diagonal) - 1;
    while (++xmy <= diagonal_getMaxXmy(diagonal)) {
        double *current = dpDiagonal_getCell(dPDiagonal, xmy);
        double *lower = dpDiagonal_getCell(dPDiagonalM1, xmy - 1);
        double *middle = dPDiagonalM2 == NULL ? NULL : dpDiagonal_getCell(ppdPDiagonalM2, xmy);
        double *upper = dpDiagonal_getCell(pdPDiagonalM1, xmy + 1);
        totalProbability = logAdd(totalProbability, calculateTotalProbability(current, lower, middle, upper));
    }

    xmy = diagonal_getMinXmy(diagonal) - 1;
    while (++xmy <= diagonal_getMaxXmy(diagonal)) {
        double *current = dpDiagonal_getCell(dPDiagonal, xmy);
        double *lower = dpDiagonal_getCell(dPDiagonalM1, xmy - 1);
        double *middle = dPDiagonalM2 == NULL ? NULL : dpDiagonal_getCell(ppdPDiagonalM2, xmy);
        double *upper = dpDiagonal_getCell(pdPDiagonalM1, xmy + 1);
        totalProbability = logAdd(totalProbability, calculateTotalProbability(current, lower, middle, upper));
    }
}

stList *getAlignedPairs(stList *anchorPairs, const char *sX, const char *sY, const int32_t lX, const int32_t lY,
        int32_t stateNumber, const double *startStates, const double *endStates, int64_t minCellsBeforeBackPass,
        int32_t skinnyDiagonalWidth, int32_t traceBackDiagonals, double threshold) {
    //Prerequisites
    assert(traceBackDiagonals >= 2);
    assert(threshold >= 0.0);
    assert(threshold <= 1.0);
    assert(skinnyDiagonalWidth >= 1);
    assert(
            minCellsBeforeBackPass > (traceBackDiagonals + skinnyDiagonalWidth) * (traceBackDiagonals
                    + skinnyDiagonalWidth));

    //This list of pairs to be returned. Not in any order, but points must be unique
    stList *alignedPairs = stList_construct3(0, (int32_t(*)(void *)) stIntTuple_destruct);

    //Primitives for the forward matrix recursion
    BandIterator *forwardBandIterator = bandIterator_construct(anchorPairs);
    DpMatrix *forwardDpMatrix = dpMatrix_construct(lX, lY, stateNumber);
    dpDiagonal_initialiseValues(&dpMatrix_getDiagonal(0), startStates); //Initialise forward matrix.

    //Backward matrix.
    DpMatrix *backwardDPMatrix = dpMatrix_construct(lX, lY, stateNumber);

    int32_t tracedBackTo = 0;
    while (1) { //Loop that moves through the matrix forward
        Diagonal diagonal = bandIterator_getNext(bandIterator);
        if (diagonal_getXay(&diagonal) > lX + lY) { //Termination
            break;
        }

        //Forward calculation
        dpMatrix_setDiagonal(forwardDpMatrix, diagonal);
        doForwardCalculation(diagonal, forwardDpMatrix);

        bool atEnd = diagonal_getXay(&diagonal) == lX + lY; //Condition true at the end of the matrix
        bool tracebackPoint = dpMatrix_getTotalActiveCells(dpMatrix) > minCellsBeforeBackPass && diagonal_getWidth(
                diagonal) <= skinnyDiagonalWidth; //Condition true when we want to do an intermediate traceback.

        //Traceback
        if (atEnd || tracebackPoint) {
            //Initialise the last row (until now) of the backward matrix to represent an end point
            dpMatrix_setDiagonal(backwardDPMatrix, &diagonal);
            dpDiagonal_initialiseValues(dpMatrix_getDiagonal(backwardDPMatrix, diagonal_getXay(&diagonal)), endStates);
            BandIterator *backwardBandIterator = bandIterator_clone(forwardBandIterator);

            if (!atEnd) { //Leave a gap in the trace back
                //bandIterator_getPrevious(bandIterator2); //Ignore the last, as that's the one we've set values for
                bandIterator_getPrevious(backwardBandIterator);
                for (int32_t i = 0; i < traceBackDiagonals; i++) {
                    Diagonal diagonal2 = bandIterator_getPrevious(backwardBandIterator);
                    dpMatrix_setDiagonal(backwardDPMatrix, &diagonal2);
                    doBackwardCalculation(&diagonal2, backwardDPMatrix);
                    dpMatrix_deleteDiagonal(backwardDPMatrix, diagonal_getXay(&diagonal2) + 2); //Delete backward diagonal after last access in backward calculation
                }
            }

            //Do walk back
            Diagonal diagonal2 = bandIterator_getPrevious(backwardBandIterator);
            int32_t tracedBackFrom = diagonal_getXay(&diagonal2);
            doPosteriorCalculation(diagonal2, forwardDpMatrix, backwardDPMatrix, threshold, alignedPairs);
            diagonal2 = bandIterator_getPrevious(backwardBandIterator);
            while (diagonal_getXay(diagonal2) > tracedBackTo) {
                dpMatrix_setDiagonal(backwardDPMatrix, &diagonal2);
                doBackwardCalculation(&diagonal2, backwardDPMatrix);
                dpMatrix_deleteDiagonal(backwardDPMatrix, diagonal_getXay(&diagonal2) + 2); //Delete backward diagonal after last access in backward calculation
                doPosteriorCalculation(&diagonal2, forwardDpMatrix, backwardDPMatrix, threshold, alignedPairs);
                dpMatrix_deleteDiagonal(forwardDpMatrix, diagonal_getXay(&diagonal2) - 1); //Delete forward diagonal after last access in posterior calculation
            }
            tracedBackTo = tracedBackFrom;
            //Delete the two last back diagonals
            dpMatrix_deleteDiagonal(backwardDPMatrix, tracedBackTo);
            dpMatrix_deleteDiagonal(backwardDPMatrix, tracedBackTo + 1);
            bandIterator_destruct(backwardBandIterator);

            //Check memory state.
            assert(dpMatrix_getActiveDiagonalNumber(backwardDPMatrix) == 0);
            if (!atEnd) {
                assert(dpMatrix_getActiveDiagonalNumber(forwardDPMatrix) == traceBackDiagonals + 1);
            }
        }
    }
    assert(tracedBackTo == lX + lY);
    //Check memory
    assert(dpMatrix_getActiveDiagonalNumber(forwardDPMatrix) == 2);
    dpMatrix_deleteDiagonal(forwardDpMatrix, 0);
    dpMatrix_deleteDiagonal(forwardDpMatrix, lX + lY);
    assert(dpMatrix_getActiveDiagonalNumber(forwardDPMatrix) == 0);
    dpMatrix_destruct(forwardDpMatrix);
    dpMatrix_destruct(backwardDPMatrix);
    bandIterator_destruct(forwardBandIterator);
}

