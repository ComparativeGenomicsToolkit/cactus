/*
 * diagonalIterator.c
 *
 *  Created on: 1 Mar 2012
 *      Author: benedictpaten
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>

#include "sonLib.h"
#include "pairwiseAligner.h"
#include "pairwiseAlignment.h"

///////////////////////////////////
///////////////////////////////////
//Diagonal
//
//Structure for working with x-y diagonal of dp matrix
///////////////////////////////////
///////////////////////////////////

const char *PAIRWISE_ALIGNMENT_EXCEPTION_ID = "PAIRWISE_ALIGNMENT_EXCEPTION";

Diagonal diagonal_construct(int32_t xay, int32_t xmyL, int32_t xmyR) {
    if ((xay + xmyL) % 2 != 0 || (xay + xmyR) % 2 != 0 || xmyL > xmyR) {
        stThrowNew(PAIRWISE_ALIGNMENT_EXCEPTION_ID,
                "Attempt to create diagonal with invalid coordinates: xay %i xmyL %i xmyR %i", xay, xmyL, xmyR);
    }
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
    return (diagonal.xmyR - diagonal.xmyL) / 2 + 1;
}

inline int32_t diagonal_getXCoordinate(int32_t xay, int32_t xmy) {
    assert((xay + xmy) % 2 == 0);
    return (xay + xmy) / 2;
}

inline int32_t diagonal_equals(Diagonal diagonal1, Diagonal diagonal2) {
    return diagonal1.xay == diagonal2.xay && diagonal1.xmyL == diagonal2.xmyL && diagonal1.xmyR == diagonal2.xmyR;
}

inline int32_t diagonal_getYCoordinate(int32_t xay, int32_t xmy) {
    assert((xay - xmy) % 2 == 0);
    return (xay - xmy) / 2;
}

inline char *diagonal_getString(Diagonal diagonal) {
    return stString_print("Diagonal, xay: %i xmyL %i, xmyR: %i", diagonal_getXay(diagonal),
            diagonal_getMinXmy(diagonal), diagonal_getMaxXmy(diagonal));
}

///////////////////////////////////
///////////////////////////////////
//Band Iterator
//
//Iterator for walking along x+y diagonals in banded fashion
//(using a set of anchor constraints)
///////////////////////////////////
///////////////////////////////////

struct _band {
    Diagonal *diagonals;
    int32_t lXalY;
};

static int32_t band_avoidOffByOne(int32_t xay, int32_t xmy) {
    return (xay + xmy) % 2 == 0 ? xmy : xmy + 1;
}

static void band_setCurrentDiagonalP(int32_t *xmy, int32_t i, int32_t j, int32_t k) {
    if (i < j) {
        *xmy += 2 * (j - i) * k;
    }
}

static Diagonal band_setCurrentDiagonal(int32_t xay, int32_t xL, int32_t yL, int32_t xU, int32_t yU) {
    Diagonal diagonal;
    diagonal.xay = xay;
    diagonal.xmyL = xL - yL;
    diagonal.xmyR = xU - yU;

    //Avoid in-between undefined x,y coordinate positions when intersecting xay and xmy.
    diagonal.xmyL = band_avoidOffByOne(xay, diagonal.xmyL);
    diagonal.xmyR = band_avoidOffByOne(xay, diagonal.xmyR);

    //Bound the xmy coordinates by the xL, yL and xU, yU band boundaries
    band_setCurrentDiagonalP(&diagonal.xmyL, diagonal_getXCoordinate(diagonal.xay, diagonal.xmyL), xL, 1);
    band_setCurrentDiagonalP(&diagonal.xmyL, yL, diagonal_getYCoordinate(diagonal.xay, diagonal.xmyL), 1);
    band_setCurrentDiagonalP(&diagonal.xmyR, xU, diagonal_getXCoordinate(diagonal.xay, diagonal.xmyR), -1);
    band_setCurrentDiagonalP(&diagonal.xmyR, diagonal_getYCoordinate(diagonal.xay, diagonal.xmyR), yU, -1);

    return diagonal;
}

static int32_t band_boundCoordinate(int32_t z, int32_t lZ) {
    return z < 0 ? 0 : (z > lZ ? lZ : z);
}

Band *band_construct(stList *anchorPairs, int32_t lX, int32_t lY, int32_t expansion) {
    assert(lX >= 0);
    assert(lY >= 0);

    if (expansion % 2 != 0) {
        stThrowNew(PAIRWISE_ALIGNMENT_EXCEPTION_ID, "Expansion around anchors is not even, expansion: %i\n", expansion);
    }

    Band *band = st_malloc(sizeof(Band));
    band->diagonals = st_malloc(sizeof(Diagonal) * (lX + lY + 1));
    band->lXalY = lX + lY;

    //Now initialise the diagonals
    int32_t anchorPairIndex = 0;
    int32_t xay = 0;
    int32_t pxay = 0, pxmy = 0;
    int32_t nxay = 0, nxmy = 0;
    int32_t xL = 0, yL = 0, xU = 0, yU = 0;

    while (xay <= band->lXalY) {
        band->diagonals[xay] = band_setCurrentDiagonal(xay, xL, yL, xU, yU);
        if (nxay == xay++) {
            //The previous diagonals become the next
            pxay = nxay;
            pxmy = nxmy;

            int32_t x = lX, y = lY;
            if (anchorPairIndex < stList_length(anchorPairs)) {
                stIntTuple *anchorPair = stList_get(anchorPairs, anchorPairIndex++);
                x = stIntTuple_getPosition(anchorPair, 0) + 1; //Plus ones, because matrix coordinates are +1 the sequence ones
                y = stIntTuple_getPosition(anchorPair, 1) + 1;
                if (x <= diagonal_getXCoordinate(pxay, pxmy) || y <= diagonal_getYCoordinate(pxay, pxmy)) {
                    stThrowNew(PAIRWISE_ALIGNMENT_EXCEPTION_ID,
                            "New anchor does not proceed previous one: x %i y %i pX %i pY %i\n", x, y,
                            diagonal_getXCoordinate(pxay, pxmy), diagonal_getYCoordinate(pxay, pxmy));
                }
            }

            nxay = x + y;
            nxmy = x - y;

            //Now call to set the lower and upper x,y coordinates
            xL = band_boundCoordinate(diagonal_getXCoordinate(pxay, pxmy - expansion), lX);
            yL = band_boundCoordinate(diagonal_getYCoordinate(nxay, nxmy - expansion), lY);
            xU = band_boundCoordinate(diagonal_getXCoordinate(nxay, nxmy + expansion), lX);
            yU = band_boundCoordinate(diagonal_getYCoordinate(pxay, pxmy + expansion), lY);
        }
    }

    return band;
}

void band_destruct(Band *band) {
    free(band->diagonals);
    free(band);
}

struct _bandIterator {
    Band *band;
    int32_t index;
};

BandIterator *bandIterator_construct(Band *band) {
    BandIterator *bandIterator = st_malloc(sizeof(BandIterator));
    bandIterator->band = band;
    bandIterator->index = 0;
    return bandIterator;
}

BandIterator *bandIterator_clone(BandIterator *bandIterator) {
    BandIterator *bandIterator2 = st_malloc(sizeof(BandIterator));
    memcpy(bandIterator2, bandIterator, sizeof(BandIterator));
    return bandIterator2;
}

void bandIterator_destruct(BandIterator *bandIterator) {
    free(bandIterator);
}

Diagonal bandIterator_getNext(BandIterator *bandIterator) {
    Diagonal diagonal =
            bandIterator->band->diagonals[bandIterator->index > bandIterator->band->lXalY ? bandIterator->band->lXalY
                    : bandIterator->index];
    if (bandIterator->index <= bandIterator->band->lXalY) {
        bandIterator->index++;
    }
    return diagonal;
}

Diagonal bandIterator_getPrevious(BandIterator *bandIterator) {
    if (bandIterator->index > 0) {
        bandIterator->index--;
    }
    return bandIterator->band->diagonals[bandIterator->index];
}

///////////////////////////////////
///////////////////////////////////
//Log Add functions
//
//Interpolation function for doing log add
///////////////////////////////////
///////////////////////////////////

#define logUnderflowThreshold 7.5
#define posteriorMatchThreshold 0.01

static inline double lookup(double x) {
    //return log (exp (x) + 1);
#ifdef BEN_DEBUG
    assert (x >= 0.00f);
    assert (x <= logUnderflowThreshold);
#endif
    if (x <= 1.00f)
        return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
    if (x <= 2.50f)
        return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
    if (x <= 4.50f)
        return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
    return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
}

double logAdd(double x, double y) {
    if (x < y)
        return (x == LOG_ZERO || y - x >= logUnderflowThreshold) ? y : lookup(y - x) + x;
    return (y == LOG_ZERO || x - y >= logUnderflowThreshold) ? x : lookup(x - y) + y;
}

///////////////////////////////////
///////////////////////////////////
//Symbols
//
//Emissions probs/functions to convert to symbol sequence
///////////////////////////////////
///////////////////////////////////

Symbol symbol_convertCharToSymbol(char i) {
    switch (i) {
        case 'A':
        case 'a':
            return a;
        case 'C':
        case 'c':
            return c;
        case 'G':
        case 'g':
            return g;
        case 'T':
        case 't':
            return t;
        default:
            return n;
    }
}

Symbol *symbol_convertStringToSymbols(const char *s, int32_t sL) {
#ifdef BEN_DEBUG
    assert(sL >= 0);
    assert(strlen(s) == sL);
#endif
    Symbol *cS = st_malloc(sL * sizeof(Symbol));
    for (int32_t i = 0; i < sL; i++) {
        cS[i] = symbol_convertCharToSymbol(s[i]);
    }
    return cS;
}

SymbolString symbolString_construct(const char *sequence, int32_t length) {
    SymbolString symbolString;
    symbolString.sequence = symbol_convertStringToSymbols(sequence, length);
    symbolString.length = length;
    return symbolString;
}

void symbolString_destruct(SymbolString s) {
    free(s.sequence);
}

static void symbol_check(Symbol c) {
#ifdef BEN_DEBUG
    assert(c >= 0 && c < 5);
#endif
}

#define EMISSION_MATCH -2.1149196655034745 //log(0.12064298095701059);
#define EMISSION_TRANSVERSION -4.5691014376830479 //log(0.010367271172731285);
#define EMISSION_TRANSITION -3.9833860032220842 //log(0.01862247669752685);
#define EMISSION_MATCH_N -3.2188758248682006 //log(0.04);
double symbol_matchProb(Symbol cX, Symbol cY) {
    symbol_check(cX);
    symbol_check(cY);
    //Symmetric matrix of transition probabilities.
    static const double matchM[25] = { EMISSION_MATCH, EMISSION_TRANSVERSION, EMISSION_TRANSITION,
            EMISSION_TRANSVERSION, EMISSION_MATCH_N, EMISSION_TRANSVERSION, EMISSION_MATCH, EMISSION_TRANSVERSION,
            EMISSION_TRANSITION, EMISSION_MATCH_N, EMISSION_TRANSITION, EMISSION_TRANSVERSION, EMISSION_MATCH,
            EMISSION_TRANSVERSION, EMISSION_MATCH_N, EMISSION_TRANSVERSION, EMISSION_TRANSITION, EMISSION_TRANSVERSION,
            EMISSION_MATCH, EMISSION_MATCH_N, EMISSION_MATCH_N, EMISSION_MATCH_N, EMISSION_MATCH_N, EMISSION_MATCH_N,
            EMISSION_MATCH_N };
    return matchM[cX * 5 + cY];
}

double symbol_gapProb(Symbol cZ) {
    symbol_check(cZ);
    return -1.6094379124341003; //log(0.2) = -1.6094379124341003
}

///////////////////////////////////
///////////////////////////////////
//State
//
//A cell is a set of states associated with an x, y coordinate.
//These functions do the state transitions for the pairwise
//alignment model.
///////////////////////////////////
///////////////////////////////////

static void state_check(State s) {
#ifdef BEN_DEBUG
    assert(c >= 0 && c < STATE_NUMBER);
#endif
}

double state_startStateProb(State state) {
    //Match state is like going to a match.
    state_check(state);
    return state == match ? 0 : LOG_ZERO;
}

double state_endStateProb(State state) {
    //End state is like to going to a match
    state_check(state);
    static const double endStates[5] = { TRANSITION_MATCH_CONTINUE, TRANSITION_MATCH_FROM_SHORT_GAP,
            TRANSITION_MATCH_FROM_SHORT_GAP, TRANSITION_MATCH_FROM_LONG_GAP, TRANSITION_MATCH_FROM_LONG_GAP };
    return endStates[state];
}

///////////////////////////////////
///////////////////////////////////
//Cell calculations
//
//A cell is a set of states associated with an x, y coordinate.
//These functions do the forward/backward calculations for the pairwise
//alignment model.
///////////////////////////////////
///////////////////////////////////

static inline void cell_calculate(double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY,
        void(*doTransition)(double *, double *, int32_t, int32_t, double, double)) {
    if (lower != NULL) {
        double eP = symbol_gapProb(cX);
        doTransition(lower, current, match, shortGapX, eP, TRANSITION_GAP_SHORT_OPEN);
        doTransition(lower, current, shortGapX, shortGapX, eP, TRANSITION_GAP_SHORT_EXTEND);
        doTransition(lower, current, shortGapY, shortGapX, eP, TRANSITION_GAP_SHORT_SWITCH);
        doTransition(lower, current, match, longGapX, eP, TRANSITION_GAP_LONG_OPEN);
        doTransition(lower, current, longGapX, longGapX, eP, TRANSITION_GAP_LONG_EXTEND);
    }
    if (middle != NULL) {
        double eP = symbol_matchProb(cX, cY);
        doTransition(middle, current, match, match, eP, TRANSITION_MATCH_CONTINUE);
        doTransition(middle, current, shortGapX, match, eP, TRANSITION_MATCH_FROM_SHORT_GAP);
        doTransition(middle, current, shortGapY, match, eP, TRANSITION_MATCH_FROM_SHORT_GAP);
        doTransition(middle, current, longGapX, match, eP, TRANSITION_MATCH_FROM_LONG_GAP);
        doTransition(middle, current, longGapY, match, eP, TRANSITION_MATCH_FROM_LONG_GAP);
    }
    if (upper != NULL) {
        double eP = symbol_gapProb(cY);
        doTransition(upper, current, match, shortGapY, eP, TRANSITION_GAP_SHORT_OPEN);
        doTransition(upper, current, shortGapY, shortGapY, eP, TRANSITION_GAP_SHORT_EXTEND);
        doTransition(upper, current, shortGapX, shortGapY, eP, TRANSITION_GAP_SHORT_SWITCH);
        doTransition(upper, current, match, longGapY, eP, TRANSITION_GAP_LONG_OPEN);
        doTransition(upper, current, longGapY, longGapY, eP, TRANSITION_GAP_LONG_EXTEND);
    }
}

static inline void doTransitionForward(double *fromCells, double *toCells, int32_t from, int32_t to, double eP,
        double tP) {
    toCells[to] = logAdd(toCells[to], fromCells[from] + (eP + tP));
}

void cell_calculateForward(double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY) {
    cell_calculate(current, lower, middle, upper, cX, cY, doTransitionForward);
}

static inline void doTransitionBackward(double *fromCells, double *toCells, int32_t from, int32_t to, double eP,
        double tP) {
    fromCells[from] = logAdd(fromCells[from], toCells[to] + (eP + tP));
}

void cell_calculateBackward(double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY) {
    cell_calculate(current, lower, middle, upper, cX, cY, doTransitionBackward);
}

double cell_dotProduct(double *cell1, double *cell2) {
    double totalProb = cell1[0] + cell2[0];
    for (int32_t i = 1; i < STATE_NUMBER; i++) {
        totalProb = logAdd(totalProb, cell1[i] + cell2[i]);
    }
    return totalProb;
}

double cell_dotProduct2(double *cell, double(*getStateValue)(State)) {
    double totalProb = cell[0] + getStateValue(0);
    for (int32_t i = 1; i < STATE_NUMBER; i++) {
        totalProb = logAdd(totalProb, cell[i] + getStateValue(i));
    }
    return totalProb;
}

///////////////////////////////////
///////////////////////////////////
//DpDiagonal
//
//Structure for storing a x-y diagonal of the dp matrix
///////////////////////////////////
///////////////////////////////////

struct _dpDiagonal {
    Diagonal diagonal;
    double *cells;
};

DpDiagonal *dpDiagonal_construct(Diagonal diagonal) {
    DpDiagonal *dpDiagonal = st_malloc(sizeof(DpDiagonal));
    dpDiagonal->diagonal = diagonal;
    dpDiagonal->cells = st_malloc(sizeof(double) * STATE_NUMBER * diagonal_getWidth(diagonal));
    return dpDiagonal;
}

DpDiagonal *dpDiagonal_clone(DpDiagonal *diagonal) {
    DpDiagonal *diagonal2 = dpDiagonal_construct(diagonal->diagonal);
    memcpy(diagonal2->cells, diagonal->cells, sizeof(double) * diagonal_getWidth(diagonal->diagonal) * STATE_NUMBER);
    return diagonal2;
}

bool dpDiagonal_equals(DpDiagonal *diagonal1, DpDiagonal *diagonal2) {
    if (!diagonal_equals(diagonal1->diagonal, diagonal2->diagonal)) {
        return 0;
    }
    for (int32_t i = 0; i < diagonal_getWidth(diagonal1->diagonal) * STATE_NUMBER; i++) {
        if (diagonal1->cells[i] != diagonal2->cells[i]) {
            return 0;
        }
    }
    return 1;
}

void dpDiagonal_destruct(DpDiagonal *dpDiagonal) {
    free(dpDiagonal->cells);
    free(dpDiagonal);
}

double *dpDiagonal_getCell(DpDiagonal *dpDiagonal, int32_t xmy) {
    if (xmy < dpDiagonal->diagonal.xmyL || xmy > dpDiagonal->diagonal.xmyR) {
        return NULL;
    }
    assert((diagonal_getXay(dpDiagonal->diagonal) + xmy) % 2 == 0);
    return &dpDiagonal->cells[((xmy - dpDiagonal->diagonal.xmyL) / 2) * STATE_NUMBER];
}

void dpDiagonal_zeroValues(DpDiagonal *diagonal) {
    for (int32_t i = 0; i < diagonal_getWidth(diagonal->diagonal) * STATE_NUMBER; i++) {
        diagonal->cells[i] = LOG_ZERO;
    }
}

void dpDiagonal_initialiseValues(DpDiagonal *diagonal, double(*getStateValue)(State)) {
    for (int32_t i = diagonal_getMinXmy(diagonal->diagonal); i <= diagonal_getMaxXmy(diagonal->diagonal); i += 2) {
        double *cell = dpDiagonal_getCell(diagonal, i);
        assert(cell != NULL);
        for (int32_t j = 0; j < STATE_NUMBER; j++) {
            cell[j] = getStateValue(j);
        }
    }
}

double dpDiagonal_dotProduct(DpDiagonal *diagonal1, DpDiagonal *diagonal2) {
    double totalProbability = LOG_ZERO;
    Diagonal diagonal = diagonal1->diagonal;
    int32_t xmy = diagonal_getMinXmy(diagonal);
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        totalProbability = logAdd(totalProbability,
                cell_dotProduct(dpDiagonal_getCell(diagonal1, xmy), dpDiagonal_getCell(diagonal2, xmy)));
        xmy += 2;
    }
    return totalProbability;
}

///////////////////////////////////
///////////////////////////////////
//DpMatrix
//
//Structure for storing dp-matrix
///////////////////////////////////
///////////////////////////////////

struct _dpMatrix {
    DpDiagonal **diagonals;
    int32_t diagonalNumber;
    int32_t activeDiagonals;
};

DpMatrix *dpMatrix_construct(int32_t diagonalNumber) {
    assert(diagonalNumber >= 0);
    DpMatrix *dpMatrix = st_malloc(sizeof(DpMatrix));
    dpMatrix->diagonalNumber = diagonalNumber;
    dpMatrix->diagonals = st_calloc(dpMatrix->diagonalNumber + 1, sizeof(DpDiagonal *));
    dpMatrix->activeDiagonals = 0;
    return dpMatrix;
}

void dpMatrix_destruct(DpMatrix *dpMatrix) {
    assert(dpMatrix->activeDiagonals == 0);
    free(dpMatrix->diagonals);
    free(dpMatrix);
}

DpDiagonal *dpMatrix_getDiagonal(DpMatrix *dpMatrix, int32_t xay) {
    if (xay < 0 || xay > dpMatrix->diagonalNumber) {
        return NULL;
    }
    return dpMatrix->diagonals[xay];
}

int32_t dpMatrix_getActiveDiagonalNumber(DpMatrix *dpMatrix) {
    return dpMatrix->activeDiagonals;
}

DpDiagonal *dpMatrix_createDiagonal(DpMatrix *dpMatrix, Diagonal diagonal) {
    assert(diagonal.xay >= 0);
    assert(diagonal.xay <= dpMatrix->diagonalNumber);
    assert(dpMatrix_getDiagonal(dpMatrix, diagonal.xay) == NULL);
    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal);
    dpMatrix->diagonals[diagonal_getXay(diagonal)] = dpDiagonal;
    dpMatrix->activeDiagonals++;
    return dpDiagonal;
}

void dpMatrix_deleteDiagonal(DpMatrix *dpMatrix, int32_t xay) {
    assert(xay >= 0);
    assert(xay <= dpMatrix->diagonalNumber);
    if (dpMatrix->diagonals[xay] != NULL) {
        dpMatrix->activeDiagonals--;
        assert(dpMatrix->activeDiagonals >= 0);
        dpDiagonal_destruct(dpMatrix->diagonals[xay]);
        dpMatrix->diagonals[xay] = NULL;
    }
}

///////////////////////////////////
///////////////////////////////////
//Diagonal DP Calculations
//
//Functions which do forward/backward/posterior calculations
//between diagonal rows of a dp-matrix
///////////////////////////////////
///////////////////////////////////

static Symbol getXCharacter(const SymbolString sX, int32_t xay, int32_t xmy) {
    int32_t x = diagonal_getXCoordinate(xay, xmy);
    assert(x >= 0 && x <= sX.length);
    return x > 0 ? sX.sequence[x - 1] : n;
}

static Symbol getYCharacter(const SymbolString sY, int32_t xay, int32_t xmy) {
    int32_t y = diagonal_getYCoordinate(xay, xmy);
    assert(y >= 0 && y <= sY.length);
    return y > 0 ? sY.sequence[y - 1] : n;
}

static void diagonalCalculation(DpDiagonal *dpDiagonal, DpDiagonal *dpDiagonalM1, DpDiagonal *dpDiagonalM2,
        const SymbolString sX, const SymbolString sY,
        void(*cellCalculation)(double *, double *, double *, double *, Symbol, Symbol)) {
    Diagonal diagonal = dpDiagonal->diagonal;
    int32_t xmy = diagonal_getMinXmy(diagonal);
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        Symbol x = getXCharacter(sX, diagonal_getXay(diagonal), xmy);
        Symbol y = getYCharacter(sY, diagonal_getXay(diagonal), xmy);
        double *current = dpDiagonal_getCell(dpDiagonal, xmy);
        double *lower = dpDiagonalM1 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM1, xmy - 1);
        double *middle = dpDiagonalM2 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM2, xmy);
        double *upper = dpDiagonalM1 == NULL ? NULL : dpDiagonal_getCell(dpDiagonalM1, xmy + 1);
        cellCalculation(current, lower, middle, upper, x, y);
        xmy += 2;
    }
}

void diagonalCalculationForward(int32_t xay, DpMatrix *dpMatrix,
        const SymbolString sX, const SymbolString sY) {
    diagonalCalculation(dpMatrix_getDiagonal(dpMatrix, xay), dpMatrix_getDiagonal(dpMatrix, xay - 1),
            dpMatrix_getDiagonal(dpMatrix, xay - 2), sX, sY, cell_calculateForward);
}

void diagonalCalculationBackward(int32_t xay, DpMatrix *dpMatrix, const SymbolString sX, const SymbolString sY) {
    diagonalCalculation(dpMatrix_getDiagonal(dpMatrix, xay), dpMatrix_getDiagonal(dpMatrix, xay - 1),
            dpMatrix_getDiagonal(dpMatrix, xay - 2), sX, sY, cell_calculateBackward);
}

double diagonalCalculationTotalProbability(int32_t xay, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
        const SymbolString sX, const SymbolString sY) {
    //Get the forward and backward diagonals
    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    double totalProbability = dpDiagonal_dotProduct(forwardDiagonal, backDiagonal);
    //Now calculate the contribution of matches through xay.
    forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay - 1);
    backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay + 1);
    if (backDiagonal != NULL && forwardDiagonal != NULL) {
        DpDiagonal *matchDiagonal = dpDiagonal_clone(backDiagonal);
        dpDiagonal_zeroValues(matchDiagonal);
        diagonalCalculation(matchDiagonal, NULL, forwardDiagonal, sX, sY, cell_calculateForward);
        totalProbability = logAdd(totalProbability, dpDiagonal_dotProduct(matchDiagonal, backDiagonal));
        dpDiagonal_destruct(matchDiagonal);
    }
    return totalProbability;
}

void diagonalCalculationPosteriorMatchProbs(int32_t xay,
        DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
       double threshold, double totalProbability,
        stList *alignedPairs) {
    //Get the forward and backward diagonals
    DpDiagonal *forwardDiagonal = dpMatrix_getDiagonal(forwardDpMatrix, xay);
    DpDiagonal *backDiagonal = dpMatrix_getDiagonal(backwardDpMatrix, xay);
    Diagonal diagonal = forwardDiagonal->diagonal;
    int32_t xmy = diagonal_getMinXmy(diagonal);
    //Walk over the cells computing the posteriors
    while (xmy <= diagonal_getMaxXmy(diagonal)) {
        int32_t x = diagonal_getXCoordinate(diagonal_getXay(diagonal), xmy);
        int32_t y = diagonal_getYCoordinate(diagonal_getXay(diagonal), xmy);
        if(x > 0 && y > 0) {
            double *cellForward = dpDiagonal_getCell(forwardDiagonal, xmy);
            double *cellBackward = dpDiagonal_getCell(backDiagonal, xmy);
            double posteriorProbability = exp((cellForward[match] + cellBackward[match]) - totalProbability);
#ifdef BEN_DEBUG
            assert(posteriorProbability > -0.01);
            assert(posteriorProbability <  1.01);
#endif
            if (posteriorProbability >= threshold) {
#ifdef BEN_DEBUG
                assert(x > 0 && x <= lX);
                assert(y > 0 && y <= lY);
#endif
                if (posteriorProbability > 1.0) {
                    posteriorProbability = 1.0;
                }
                posteriorProbability = floor(posteriorProbability * PAIR_ALIGNMENT_PROB_1);

                stList_append(alignedPairs, stIntTuple_construct(3, (int32_t) posteriorProbability, x - 1, y - 1));
            }
        }
        xmy += 2;
    }
}

///////////////////////////////////
///////////////////////////////////
//Banded alignment routine to calculate posterior match probs
//
//
///////////////////////////////////
///////////////////////////////////

stList *getAlignedPairsWithBanding(stList *anchorPairs, const SymbolString sX, const SymbolString sY,
        PairwiseAlignmentParameters *p) {
    //Prerequisites
    assert(p->traceBackDiagonals >= 1);
    assert(p->threshold >= 0.0);
    assert(p->threshold <= 1.0);
    assert(p->diagonalExpansion >= 0);
    assert(p->diagonalExpansion % 2 ==  0);
    assert(p->minDiagsBetweenTraceBack >= 2);
    assert(p->traceBackDiagonals + 1 < p->minDiagsBetweenTraceBack);

    //This list of pairs to be returned. Not in any order, but points must be unique
    stList *alignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);

    int32_t diagonalNumber = sX.length + sY.length;
    if(diagonalNumber == 0) { //Deal with trivial case
        return alignedPairs;
    }

    //Primitives for the forward matrix recursion
    Band *band = band_construct(anchorPairs, sX.length, sY.length, p->diagonalExpansion);
    BandIterator *forwardBandIterator = bandIterator_construct(band);
    DpMatrix *forwardDpMatrix = dpMatrix_construct(diagonalNumber);
    dpDiagonal_initialiseValues(dpMatrix_createDiagonal(forwardDpMatrix, bandIterator_getNext(forwardBandIterator)),
            state_startStateProb); //Initialise forward matrix.

    //Backward matrix.
    DpMatrix *backwardDpMatrix = dpMatrix_construct(diagonalNumber);

    int32_t tracedBackTo = 0;
    int32_t totalPosteriorCalculations = 0;
    while (1) { //Loop that moves through the matrix forward
        Diagonal diagonal = bandIterator_getNext(forwardBandIterator);

        //Forward calculation
        dpDiagonal_zeroValues(dpMatrix_createDiagonal(forwardDpMatrix, diagonal));
        diagonalCalculationForward(diagonal_getXay(diagonal), forwardDpMatrix, sX, sY);

        bool atEnd = diagonal_getXay(diagonal) == diagonalNumber; //Condition true at the end of the matrix
        bool tracebackPoint = diagonal_getXay(diagonal) >= tracedBackTo + p->minDiagsBetweenTraceBack && diagonal_getWidth(diagonal)
                <= p->diagonalExpansion * 2 + 1; //Condition true when we want to do an intermediate traceback.

        //Traceback
        if (atEnd || tracebackPoint) {
            //Initialise the last row (until now) of the backward matrix to represent an end point
            dpDiagonal_initialiseValues(dpMatrix_createDiagonal(backwardDpMatrix, diagonal), state_endStateProb);
            if (diagonal_getXay(diagonal) > tracedBackTo + 1) { //This is a diagonal between the place we trace back to and where we trace back from
                DpDiagonal *j = dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal) - 1);
                assert(j != NULL);
                dpDiagonal_zeroValues(dpMatrix_createDiagonal(backwardDpMatrix, j->diagonal));
            }

            //Do walk back
            BandIterator *backwardBandIterator = bandIterator_clone(forwardBandIterator);
            Diagonal diagonal2 = bandIterator_getPrevious(backwardBandIterator);
            assert(diagonal_getXay(diagonal2) == diagonal_getXay(diagonal));
            int32_t tracedBackFrom = diagonal_getXay(diagonal) - (atEnd ? 0 : p->traceBackDiagonals + 1);
            double totalProbability = LOG_ZERO;
            int32_t totalPosteriorCalculationsThisTraceback = 0;
            while (diagonal_getXay(diagonal2) > tracedBackTo) {
                //Create the earlier diagonal
                if (diagonal_getXay(diagonal2) > tracedBackTo + 2) {
                    DpDiagonal *j = dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2) - 2);
                    assert(j != NULL);
                    dpDiagonal_zeroValues(dpMatrix_createDiagonal(backwardDpMatrix, j->diagonal));
                }
                if (diagonal_getXay(diagonal2) > tracedBackTo + 1) {
                    diagonalCalculationBackward(diagonal_getXay(diagonal2), backwardDpMatrix, sX, sY);
                }
                if (diagonal_getXay(diagonal2) <= tracedBackFrom) {
#ifdef BEN_DEBUG
                    assert(dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2)) != NULL);
                    assert(dpMatrix_getDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2)-1) != NULL);
                    assert(dpMatrix_getDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2)) != NULL);
                    if(diagonal_getXay(diagonal2) != diagonalNumber) {
                        assert(dpMatrix_getDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2)+1) != NULL);
                    }
#endif
                    if(totalPosteriorCalculationsThisTraceback++ % 10 == 0) {
                        double newTotalProbability = diagonalCalculationTotalProbability(diagonal_getXay(diagonal2), forwardDpMatrix, backwardDpMatrix, sX, sY);
                        if(totalPosteriorCalculationsThisTraceback != 1) {
                            assert(totalProbability + 0.1 > newTotalProbability);
                            assert(newTotalProbability + 0.1 > newTotalProbability);
                        }
                        totalProbability = newTotalProbability;
                    }

                    diagonalCalculationPosteriorMatchProbs(diagonal_getXay(diagonal2), forwardDpMatrix,
                                                    backwardDpMatrix, p->threshold, totalProbability, alignedPairs);

                    if(diagonal_getXay(diagonal2) < tracedBackFrom || atEnd) {
                        dpMatrix_deleteDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2)); //Delete forward diagonal after last access in posterior calculation
                    }
                }
                if(diagonal_getXay(diagonal2)+1 <= diagonalNumber) {
                    dpMatrix_deleteDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2)+1); //Delete backward diagonal after last access in backward calculation
                }
                diagonal2 = bandIterator_getPrevious(backwardBandIterator);
            }
            tracedBackTo = tracedBackFrom;
            bandIterator_destruct(backwardBandIterator);
            dpMatrix_deleteDiagonal(backwardDpMatrix, diagonal_getXay(diagonal2)+1);
            dpMatrix_deleteDiagonal(forwardDpMatrix, diagonal_getXay(diagonal2));
            //Check memory state.
            assert(dpMatrix_getActiveDiagonalNumber(backwardDpMatrix) == 0);
            totalPosteriorCalculations += totalPosteriorCalculationsThisTraceback;
            if (!atEnd) {
                assert(dpMatrix_getActiveDiagonalNumber(forwardDpMatrix) == p->traceBackDiagonals + 2);
            }
        }

        if(atEnd) {
            break;
        }
    }
#ifdef BEN_DEBUG
    assert(totalPosteriorCalculations == diagonalNumber);
    assert(tracedBackTo == diagonalNumber);
    assert(dpMatrix_getActiveDiagonalNumber(backwardDpMatrix) == 0);
    assert(dpMatrix_getActiveDiagonalNumber(forwardDpMatrix) == 0);
#endif
    //Cleanup
    dpMatrix_destruct(forwardDpMatrix);
    dpMatrix_destruct(backwardDpMatrix);
    bandIterator_destruct(forwardBandIterator);
    band_destruct(band);

    return alignedPairs;
}

///////////////////////////////////
///////////////////////////////////
//Blast anchoring functions
//
//Use lastz to get sets of anchors
///////////////////////////////////
///////////////////////////////////

static int sortByXPlusYCoordinate(const void *i, const void *j) {
    int64_t k = stIntTuple_getPosition((stIntTuple *) i, 0) + stIntTuple_getPosition((stIntTuple *) i, 1);
    int64_t l = stIntTuple_getPosition((stIntTuple *) j, 0) + stIntTuple_getPosition((stIntTuple *) j, 1);
    return k > l ? 1 : (k < l ? -1 : 0);
}

static char *makeUpperCase(const char *s, int32_t l) {
    char *s2 = stString_copy(s);
    for (int32_t i = 0; i < l; i++) {
        s2[i] = toupper(s[i]);
    }
    return s2;
}

static void writeSequenceToFile(char *file, const char *name, const char *sequence) {
    FILE *fileHandle = fopen(file, "w");
    fprintf(fileHandle, ">%s\n%s\n", name, sequence);
    fclose(fileHandle);
}

stList *getBlastPairs(const char *sX, const char *sY, int32_t lX, int32_t lY, int32_t trim, bool repeatMask) {
    /*
     * Uses lastz to compute a bunch of monotonically increasing pairs such that for any pair of consecutive pairs in the list
     * (x1, y1) (x2, y2) in the set of aligned pairs x1 appears before x2 in X and y1 appears before y2 in Y.
     */
    stList *alignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct); //the list to put the output in

    if (lX == 0 || lY == 0) {
        return alignedPairs;
    }

    if (!repeatMask) {
        sX = makeUpperCase(sX, lX);
        sY = makeUpperCase(sY, lY);
    }

    //Write one sequence to file..
    char *tempFile1 = getTempFile();
    char *tempFile2 = NULL;

    writeSequenceToFile(tempFile1, "a", sX);

    char *command;

    if (lY > 10000) {
        tempFile2 = getTempFile();
        writeSequenceToFile(tempFile2, "b", sY);
        command = stString_print(
                "lastz --hspthresh=800 --chain --strand=plus --gapped --format=cigar --ambiguous=iupac %s %s",
                tempFile1, tempFile2);
    } else {
        command
                = stString_print(
                        "echo '>b\n%s\n' | lastz --hspthresh=800 --chain --strand=plus --gapped --format=cigar --ambiguous=iupac %s",
                        sY, tempFile1);
    }
    FILE *fileHandle = popen(command, "r");
    free(command);
    if (fileHandle == NULL) {
        st_errAbort("Problems with lastz pipe");
    }
    //Read from stream
    int32_t pxay = -1;
    struct PairwiseAlignment *pA;
    while ((pA = cigarRead(fileHandle)) != NULL) {
        int32_t j = pA->start1;
        int32_t k = pA->start2;

        assert(strcmp(pA->contig1, "a") == 0);
        assert(strcmp(pA->contig2, "b") == 0);
        assert(pA->strand1);
        assert(pA->strand2);

        for (int32_t i = 0; i < pA->operationList->length; i++) {
            struct AlignmentOperation *op = pA->operationList->list[i];
            if (op->opType == PAIRWISE_MATCH) {
                for (int32_t l = trim; l < op->length - trim; l++) {
                    int32_t x = j + l;
                    int32_t y = k + l;
                    if (x + y > pxay) {
                        stList_append(alignedPairs, stIntTuple_construct(2, j + l, k + l));
                        pxay = x + y;
                    }
                }
            }
            if (op->opType != PAIRWISE_INDEL_Y) {
                j += op->length;
            }
            if (op->opType != PAIRWISE_INDEL_X) {
                k += op->length;
            }
        }

        assert(j == pA->end1);
        assert(k == pA->end2);
        destructPairwiseAlignment(pA);
    }
    int32_t status = pclose(fileHandle);
    if (status != 0) {
        st_errnoAbort("pclose failed when getting rid of lastz pipe with value %i", status);
    }

    stList_sort(alignedPairs, sortByXPlusYCoordinate); //Ensure the coordinates are increasing

    //Remove old files
    st_system("rm %s", tempFile1);
    free(tempFile1);
    if (tempFile2 != NULL) {
        st_system("rm %s", tempFile2);
        free(tempFile2);
    }

    if (!repeatMask) {
        free((char *) sX);
        free((char *) sY);
    }

    return alignedPairs;
}

static void convertPairs(stList *alignedPairs2, int32_t offsetX, int32_t offsetY) {
    /*
     * Convert the coordinates of the computed pairs.
     */
    for (int32_t k = 0; k < stList_length(alignedPairs2); k++) {
        stIntTuple *i = stList_get(alignedPairs2, k);
        assert(stIntTuple_length(i) == 2);
        stList_set(
                alignedPairs2,
                k,
                stIntTuple_construct(2,  stIntTuple_getPosition(i, 0) + offsetX,
                        stIntTuple_getPosition(i, 1) + offsetY));
        stIntTuple_destruct(i);
    }
}

static void getBlastPairsForPairwiseAlignmentParametersP(const char *sX, const char *sY, int32_t pX, int32_t pY,
        int32_t x, int32_t y, PairwiseAlignmentParameters *p, stList *combinedAnchorPairs) {
    int32_t lX2 = x - pX;
    assert(lX2 >= 0);
    int32_t lY2 = y - pY;
    assert(lY2 >= 0);
    int64_t matrixSize = (int64_t) lX2 * lY2;
    if (matrixSize > p->repeatMaskMatrixBiggerThanThis) {
        char *sX2 = stString_getSubString(sX, pX, lX2);
        char *sY2 = stString_getSubString(sY, pY, lY2);
        stList *bottomLevelAnchorPairs = getBlastPairs(sX2, sY2, lX2, lY2, p->constraintDiagonalTrim, 0);
        convertPairs(bottomLevelAnchorPairs, pX, pY);
        free(sX2);
        free(sY2);
        stList_appendAll(combinedAnchorPairs, bottomLevelAnchorPairs);
        stList_setDestructor(bottomLevelAnchorPairs, NULL);
        stList_destruct(bottomLevelAnchorPairs);
    }
}

stList *getBlastPairsForPairwiseAlignmentParameters(const char *sX, const char *sY, const int32_t lX, const int32_t lY,
        PairwiseAlignmentParameters *p) {
    if ((int64_t) lX * lY <= p->anchorMatrixBiggerThanThis) {
        return stList_construct();
    }
    stList *topLevelAnchorPairs = getBlastPairs(sX, sY, lX, lY, p->constraintDiagonalTrim, 1);
    int32_t pX = 0;
    int32_t pY = 0;
    stList *combinedAnchorPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(topLevelAnchorPairs); i++) {
        stIntTuple *anchorPair = stList_get(topLevelAnchorPairs, i);
        int32_t x = stIntTuple_getPosition(anchorPair, 0);
        int32_t y = stIntTuple_getPosition(anchorPair, 1);
        getBlastPairsForPairwiseAlignmentParametersP(sX, sY, pX, pY, x, y, p, combinedAnchorPairs);
        stList_append(combinedAnchorPairs, anchorPair);
        pX = x + 1;
        pY = y + 1;
    }
    getBlastPairsForPairwiseAlignmentParametersP(sX, sY, pX, pY, lX, lY, p, combinedAnchorPairs);
    stList_setDestructor(topLevelAnchorPairs, NULL);
    stList_destruct(topLevelAnchorPairs);
    return combinedAnchorPairs;
}

///////////////////////////////////
///////////////////////////////////
//Split large gap functions
//
//Functions to split up alignment around gaps in the anchors that are too large.
///////////////////////////////////
///////////////////////////////////

static void getSplitPointsP(int32_t pX, int32_t pY, int32_t x, int32_t y, stList *splitPoints,
        PairwiseAlignmentParameters *p) {
    int32_t lX2 = x - pX;
    int32_t lY2 = y = pY;
    int64_t matrixSize = (int64_t) lX2 * lY2;
    if (matrixSize > p->splitMatrixBiggerThanThis) {
        int32_t maxSequenceLength = sqrt(p->splitMatrixBiggerThanThis);
        int32_t hX = lX2 / 2 > maxSequenceLength ? maxSequenceLength : lX2 / 2;
        int32_t hY = lY2 / 2 > maxSequenceLength ? maxSequenceLength : lY2 / 2;
        stList_append(splitPoints, stIntTuple_construct(2, pX + hX + 1, pY + hY + 1));
        stList_append(splitPoints, stIntTuple_construct(2, x - hX, y - hY));
    }
}

stList *getSplitPoints(stList *anchorPairs, int32_t lX, int32_t lY, PairwiseAlignmentParameters *p) {
    int32_t pX = 0;
    int32_t pY = 0;
    stList *splitPoints = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    stList_append(splitPoints, stIntTuple_construct(2, 0, 0));
    for (int32_t i = 0; i < stList_length(anchorPairs); i++) {
        stIntTuple *anchorPair = stList_get(anchorPairs, i);
        int32_t x = stIntTuple_getPosition(anchorPair, 0);
        int32_t y = stIntTuple_getPosition(anchorPair, 1);
        getSplitPointsP(pX, pY, x, y, splitPoints, p);
        pX = x + 1;
        pY = y + 1;
    }
    getSplitPointsP(pX, pY, lX, lY, splitPoints, p);
    stList_append(splitPoints, stIntTuple_construct(2, lX, lY));
    return splitPoints;
}

stList *splitAlignmentsByLargeGaps(stList *anchorPairs, const char *sX, const char *sY, int32_t lX, int32_t lY,
        PairwiseAlignmentParameters *p) {
    stList *splitPoints = getSplitPoints(anchorPairs, lX, lY, p);
    stListIterator *anchorPairIterator = stList_getIterator(splitPoints);
    stList *alignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    //Now to the actual alignments
    for (int32_t i = 0; i < stList_length(splitPoints); i += 2) {
        stIntTuple *from = stList_get(splitPoints, i);
        stIntTuple *to = stList_get(splitPoints, i + 1);
        int32_t x = stIntTuple_getPosition(from, 0);
        int32_t y = stIntTuple_getPosition(from, 1);
        int32_t lX2 = stIntTuple_getPosition(to, 0) - x;
        int32_t lY2 = stIntTuple_getPosition(to, 1) - y;
        char *sX2 = stString_getSubString(sX, x, lX2);
        char *sY2 = stString_getSubString(sY, y, lY2);
        stList *subListOfAnchorPoints = stList_construct();
        stIntTuple *anchorPair;
        while ((anchorPair = stList_getNext(anchorPairIterator)) != NULL) {
            int32_t xay = stIntTuple_getPosition(anchorPair, 0) + stIntTuple_getPosition(anchorPair, 1);
            assert(xay >= x + y);
            if (xay > x + y + lX2 + lY2) {
                stList_getPrevious(anchorPairIterator);
                break;
            }
            stList_append(subListOfAnchorPoints, anchorPair);
        }
        SymbolString sX3 = symbolString_construct(sX2, lX2);
        SymbolString sY3 = symbolString_construct(sY2, lY2);
        stList *subListOfAlignedPairs = getAlignedPairsWithBanding(subListOfAnchorPoints, sX3, sY3, p);
        stList_appendAll(alignedPairs, subListOfAlignedPairs);
        stList_setDestructor(subListOfAlignedPairs, NULL);
        stList_destruct(subListOfAlignedPairs);
        stList_destruct(subListOfAnchorPoints);
        free(sX2);
        free(sY2);
        symbolString_destruct(sX3);
        symbolString_destruct(sY3);
    }
    stList_destructIterator(anchorPairIterator);
    stList_destruct(splitPoints);
    return alignedPairs;
}

///////////////////////////////////
///////////////////////////////////
//Core public functions
///////////////////////////////////
///////////////////////////////////

PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters_construct() {
    PairwiseAlignmentParameters *p = st_malloc(sizeof(PairwiseAlignmentParameters));
    p->threshold = 0.01;
    p->minDiagsBetweenTraceBack = 1000;
    p->traceBackDiagonals = 20;
    p->diagonalExpansion = 10;
    p->constraintDiagonalTrim = 5;
    p->anchorMatrixBiggerThanThis = 500 * 500;
    p->repeatMaskMatrixBiggerThanThis = 500 * 500;
    p->splitMatrixBiggerThanThis = 3000 * 3000;
    p->alignAmbiguityCharacters = 1;
    return p;
}

void pairwiseAlignmentBandingParameters_destruct(PairwiseAlignmentParameters *p) {
    free(p);
}

stList *getAlignedPairs(const char *sX, const char *sY, PairwiseAlignmentParameters *p) {
    const int32_t lX = strlen(sX);
    const int32_t lY = strlen(sY);
    stList *anchorPairs = getBlastPairsForPairwiseAlignmentParameters(sX, sY, lX, lY, p);
    stList *alignedPairs = splitAlignmentsByLargeGaps(anchorPairs, sX, sY, lX, lY, p);
    stList_destruct(anchorPairs);
    return alignedPairs;
}

