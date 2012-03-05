/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * pairwiseAligner.h
 *
 *  Created on: 5 Jul 2010
 *      Author: benedictpaten
 */

#ifndef PAIRWISEALIGNER_H_
#define PAIRWISEALIGNER_H_

//The exception string
extern const char *PAIRWISE_ALIGNMENT_EXCEPTION_ID;

//Constant that gives the integer value equal to probability 1. Integer probability zero is always 0.
#define PAIR_ALIGNMENT_PROB_1 10000000

typedef struct _pairwiseAlignmentBandingParameters {
    int32_t threshold; //Minimum posterior probability of a match to be added to the output
    int32_t minDiagsBetweenTraceBack; //Minimum x+y diagonals to leave between doing traceback.
    int32_t traceBackDiagonals; //Number of diagonals to leave between trace back diagonal
    int32_t diagonalExpansion; //The number of x-y diagonals to expand around an anchor point
    int32_t constraintDiagonalTrim; //Amount to remove from a diagonal to be considered for a banding constraint
    int64_t anchorMatrixBiggerThanThis; //Search for anchors on any matrix bigger than this
    int64_t repeatMaskMatrixBiggerThanThis; //Any matrix in the anchors bigger than this is searched for anchors using non-repeat masked sequences.
    int64_t splitMatrixBiggerThanThis; //Any matrix in the anchors bigger than this is split into two.
    bool alignAmbiguityCharacters;
} PairwiseAlignmentParameters;

PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters_construct();

void pairwiseAlignmentBandingParameters_destruct(PairwiseAlignmentParameters *p);

/*
 * Gets the set of posterior match probabilities under a simple HMM model of alignment for two DNA sequences.
 */
stList *getAlignedPairs(const char *string1, const char *string2, PairwiseAlignmentParameters *p);

/*
 * Methods tested and possibly useful elsewhere
 */

////Diagonal

typedef struct _diagonal {
    int32_t xay; //x + y coordinate
    int32_t xmyL; //smallest x - y coordinate
    int32_t xmyR; //largest x - y coordinate
} Diagonal;

Diagonal diagonal_construct(int32_t xay, int32_t xmyL, int32_t xmyR);

int32_t diagonal_getXay(Diagonal diagonal);

int32_t diagonal_getMinXmy(Diagonal diagonal);

int32_t diagonal_getMaxXmy(Diagonal diagonal);

int32_t diagonal_getWidth(Diagonal diagonal);

int32_t diagonal_getXCoordinate(int32_t xay, int32_t xmy);

int32_t diagonal_getYCoordinate(int32_t xay, int32_t xmy);

int32_t diagonal_equals(Diagonal diagonal1, Diagonal diagonal2);

char *diagonal_getString(Diagonal diagonal);

////Band iterator.

typedef struct _bandIterator BandIterator;

BandIterator *bandIterator_construct(stList *anchorPairs, int32_t lX, int32_t lY, int32_t expansion);

void bandIterator_destruct(BandIterator *bandIterator);

BandIterator *bandIterator_clone(BandIterator *bandIterator);

Diagonal bandIterator_getNext(BandIterator *bandIterator);

Diagonal bandIterator_getPrevious(BandIterator *bandIterator);

//Log add

#define LOG_ZERO -INFINITY

double logAdd(double x, double y);

//Symbols and emissions

typedef enum {
    a=0,
    c=1,
    g=2,
    t=3,
    n=4
} Symbol;

Symbol symbol_convertCharToSymbol(char i);

Symbol *symbol_convertStringToSymbols(const char *s, int32_t sL);

double symbol_matchProb(Symbol cX, Symbol cY);

double symbol_gapProb(Symbol c);

//States and transitions

#define STATE_NUMBER 5

typedef enum {
    match=0,
    shortGapX=1,
    shortGapY=2,
    longGapX=3,
    longGapY=4
} State;

#define TRANSITION_MATCH_CONTINUE -0.030064059121770816 //0.9703833696510062f
#define TRANSITION_MATCH_FROM_SHORT_GAP -1.272871422049609 //1.0 - gapExtend - gapSwitch = 0.280026392297485
#define TRANSITION_MATCH_FROM_LONG_GAP -5.673280173170473 //1.0 - gapExtend = 0.00343657420938
#define TRANSITION_GAP_SHORT_OPEN -4.34381910900448 //0.0129868352330243
#define TRANSITION_GAP_SHORT_EXTEND -0.3388262689231553 //0.7126062401851738f;
#define TRANSITION_GAP_SHORT_SWITCH -4.910694825551255 //0.0073673675173412815f;
#define TRANSITION_GAP_LONG_OPEN -6.30810595366929 //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
#define TRANSITION_GAP_LONG_EXTEND -0.003442492794189331 //0.99656342579062f;

double state_startStateProb(State state);

double state_endStateProb(State state);

//Cells (states at a given coordinate(

void cell_calculateForward(double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY);

void cell_calculateBackward(double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY);

//DpDiagonal

typedef struct _dpDiagonal DpDiagonal;

DpDiagonal *dpDiagonal_construct(Diagonal diagonal);

DpDiagonal *dpDiagonal_clone(DpDiagonal *diagonal);

void dpDiagonal_destruct(DpDiagonal *dpDiagonal);

double *dpDiagonal_getCell(DpDiagonal *dpDiagonal, int32_t xmy);

double dpDiagonal_sum(DpDiagonal *diagonal);

void dpDiagonal_zeroValues(DpDiagonal *diagonal);

void dpDiagonal_initialiseValues(DpDiagonal *diagonal, double(*getStateValue)(State));

//DpMatrix

typedef struct _dpMatrix DpMatrix;

DpMatrix *dpMatrix_construct(int32_t lX, int32_t lY);

void dpMatrix_destruct(DpMatrix *dpMatrix);

DpDiagonal *dpMatrix_getDiagonal(DpMatrix *dpMatrix, int32_t xay);

int32_t dpMatrix_getActiveDiagonalNumber(DpMatrix *dpMatrix);

DpDiagonal *dpMatrix_createDiagonal(DpMatrix *dpMatrix, Diagonal diagonal);

void dpMatrix_deleteDiagonal(DpMatrix *dpMatrix, int32_t xay);

//Diagonal calculations

void diagonalCalculationForward(Diagonal diagonal, DpMatrix *dpMatrix, const Symbol *sX, const Symbol *sY, int32_t lX,
        int32_t lY);

void diagonalCalculationBackward(Diagonal diagonal, DpMatrix *dpMatrix, const Symbol *sX, const Symbol *sY, int32_t lX,
        int32_t lY);

void diagonalCalculationPosterior(Diagonal diagonal, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
        const Symbol *sX, const Symbol *sY, int32_t lX, int32_t lY, double threshold, stList *alignedPairs);

//Banded matrix calculation

stList *getAlignedPairsWithBanding(stList *anchorPairs, const Symbol *sX, const Symbol *sY,
        const int32_t lX, const int32_t lY,
        PairwiseAlignmentParameters *p);

//Blast pairs

stList *getBlastPairs(const char *sX, const char *sY, int32_t lX, int32_t lY, int32_t trim, bool repeatMask);

stList *getBlastPairsForPairwiseAlignmentParameters(const char *sX, const char *sY, const int32_t lX, const int32_t lY,
        PairwiseAlignmentParameters *p);

//Split over large gaps

stList *getSplitPoints(stList *anchorPairs, int32_t lX, int32_t lY,
        PairwiseAlignmentParameters *p);

stList *splitAlignmentsByLargeGaps(stList *anchorPairs, const char *sX, const char *sY, int32_t lX, int32_t lY,
        PairwiseAlignmentParameters *p);

#endif /* PAIRWISEALIGNER_H_ */
