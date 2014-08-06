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

#include "bioioC.h"
#include "sonLib.h"
#include "pairwiseAlignment.h"

//The exception string
extern const char *PAIRWISE_ALIGNMENT_EXCEPTION_ID;

//Constant that gives the integer value equal to probability 1. Integer probability zero is always 0.
#define PAIR_ALIGNMENT_PROB_1 10000000

typedef struct _pairwiseAlignmentBandingParameters {
    double threshold; //Minimum posterior probability of a match to be added to the output
    int64_t minDiagsBetweenTraceBack; //Minimum x+y diagonals to leave between doing traceback.
    int64_t traceBackDiagonals; //Number of diagonals to leave between trace back diagonal
    int64_t diagonalExpansion; //The number of x-y diagonals to expand around an anchor point
    int64_t constraintDiagonalTrim; //Amount to remove from a diagonal to be considered for a banding constraint
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
stList *getAlignedPairs(const char *string1, const char *string2, PairwiseAlignmentParameters *p,  bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd);

stList *convertPairwiseForwardStrandAlignmentToAnchorPairs(struct PairwiseAlignment *pA, int64_t trim);

stList *getAlignedPairsUsingAnchors(const char *sX, const char *sY, stList *anchorPairs, PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd);

/*
 * EM training stuff.
 */

typedef struct _hmm {
    double *transitions;
    double likelihood;
} Hmm;

Hmm *hmm_constructEmpty(double pseudoExpectation);

void hmm_randomise(Hmm *hmm); //Creates normalised HMM with parameters set to small random values.

void hmm_destruct(Hmm *hmmExpectations);

void hmm_write(Hmm *hmmExpectations, FILE *fileHandle);

void hmm_addToTransitionExpectation(Hmm *hmmExpectations, int64_t from, int64_t to, double p);

double hmm_getTransition(Hmm *hmmExpectations, int64_t from, int64_t to);

void hmm_setTransition(Hmm *hmm, int64_t from, int64_t to, double p);

Hmm *hmm_loadFromFile(const char *fileName);

void hmm_normalise(Hmm *hmm);

void loadTheGlobalHmm(Hmm *hmm);

void resetTheGlobalHmm();

/*
 * Expectation calculation functions for EM algorithms.
 */

void getExpectationsUsingAnchors(Hmm *hmmExpectations, const char *sX, const char *sY, stList *anchorPairs,
        PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd);

void getExpectations(Hmm *hmmExpectations, const char *sX, const char *sY, PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd);

/*
 * Methods tested and possibly useful elsewhere
 */

////Diagonal

typedef struct _diagonal {
    int64_t xay; //x + y coordinate
    int64_t xmyL; //smallest x - y coordinate
    int64_t xmyR; //largest x - y coordinate
} Diagonal;

Diagonal diagonal_construct(int64_t xay, int64_t xmyL, int64_t xmyR);

int64_t diagonal_getXay(Diagonal diagonal);

int64_t diagonal_getMinXmy(Diagonal diagonal);

int64_t diagonal_getMaxXmy(Diagonal diagonal);

int64_t diagonal_getWidth(Diagonal diagonal);

int64_t diagonal_getXCoordinate(int64_t xay, int64_t xmy);

int64_t diagonal_getYCoordinate(int64_t xay, int64_t xmy);

int64_t diagonal_equals(Diagonal diagonal1, Diagonal diagonal2);

char *diagonal_getString(Diagonal diagonal);

//Band

typedef struct _band Band;

Band *band_construct(stList *anchorPairs, int64_t lX, int64_t lY,
        int64_t expansion);

void band_destruct(Band *band);

////Band iterator.

typedef struct _bandIterator BandIterator;

BandIterator *bandIterator_construct(Band *band);

void bandIterator_destruct(BandIterator *bandIterator);

BandIterator *bandIterator_clone(BandIterator *bandIterator);

Diagonal bandIterator_getNext(BandIterator *bandIterator);

Diagonal bandIterator_getPrevious(BandIterator *bandIterator);

//Log add

#define LOG_ZERO -INFINITY

double logAdd(double x, double y);

//Symbols

typedef enum {
    a=0,
    c=1,
    g=2,
    t=3,
    n=4
} Symbol;

Symbol symbol_convertCharToSymbol(char i);

Symbol *symbol_convertStringToSymbols(const char *s, int64_t sL);

typedef struct _symbolString {
        Symbol *sequence;
                int64_t length;
} SymbolString;

SymbolString symbolString_construct(const char *sequence, int64_t length);

//Emissions

double symbol_matchProb(Symbol cX, Symbol cY);

double symbol_gapProb(Symbol c);

//States and transitions

#define STATE_NUMBER 5

#define MATCH_STATE 0

double state_startStateProb(int64_t state);

double state_endStateProb(int64_t state);

double state_raggedEndStateProb(int64_t state);

double state_raggedStartStateProb(int64_t state);

//Cells (states at a given coordinate(

void cell_calculate(double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY,
        void(*doTransition)(double *, double *, int64_t, int64_t, double, double, void *), void *extraArgs);

void cell_calculateForward(double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY, void *extraArgs);

void cell_calculateBackward(double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY, void *extraArgs);

double cell_dotProduct(double *cell1, double *cell2);

double cell_dotProduct2(double *cell1, double (*getStateValue)(int64_t));

//DpDiagonal

typedef struct _dpDiagonal DpDiagonal;

DpDiagonal *dpDiagonal_construct(Diagonal diagonal);

DpDiagonal *dpDiagonal_clone(DpDiagonal *diagonal);

bool dpDiagonal_equals(DpDiagonal *diagonal1, DpDiagonal *diagonal2);

void dpDiagonal_destruct(DpDiagonal *dpDiagonal);

double *dpDiagonal_getCell(DpDiagonal *dpDiagonal, int64_t xmy);

double dpDiagonal_dotProduct(DpDiagonal *diagonal1, DpDiagonal *diagonal2);

void dpDiagonal_zeroValues(DpDiagonal *diagonal);

void dpDiagonal_initialiseValues(DpDiagonal *diagonal, double(*getStateValue)(int64_t));

//DpMatrix

typedef struct _dpMatrix DpMatrix;

DpMatrix *dpMatrix_construct(int64_t diagonalNumber);

void dpMatrix_destruct(DpMatrix *dpMatrix);

DpDiagonal *dpMatrix_getDiagonal(DpMatrix *dpMatrix, int64_t xay);

int64_t dpMatrix_getActiveDiagonalNumber(DpMatrix *dpMatrix);

DpDiagonal *dpMatrix_createDiagonal(DpMatrix *dpMatrix, Diagonal diagonal);

void dpMatrix_deleteDiagonal(DpMatrix *dpMatrix, int64_t xay);

//Diagonal calculations

void diagonalCalculationForward(int64_t xay, DpMatrix *dpMatrix, const SymbolString sX, const SymbolString sY);

void diagonalCalculationBackward(int64_t xay, DpMatrix *dpMatrix, const SymbolString sX, const SymbolString sY);

double diagonalCalculationTotalProbability(int64_t xay, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
        const SymbolString sX, const SymbolString sY);

void diagonalCalculationPosteriorMatchProbs(int64_t xay, DpMatrix *forwardDpMatrix, DpMatrix *backwardDpMatrix,
        const SymbolString sX, const SymbolString sY,
        double totalProbability, PairwiseAlignmentParameters *p, void *extraArgs);

//Banded matrix calculation of posterior probs

void getPosteriorProbsWithBanding(stList *anchorPairs, const SymbolString sX, const SymbolString sY,
        PairwiseAlignmentParameters *p, bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd,
        void (*diagonalPosteriorProbFn)(int64_t, DpMatrix *, DpMatrix *, const SymbolString, const SymbolString,
              double, PairwiseAlignmentParameters *, void *), void *extraArgs);

//Blast pairs

stList *getBlastPairs(const char *sX, const char *sY, int64_t lX, int64_t lY, int64_t trim, bool repeatMask);

stList *getBlastPairsForPairwiseAlignmentParameters(const char *sX, const char *sY, const int64_t lX, const int64_t lY,
        PairwiseAlignmentParameters *p);

stList *filterToRemoveOverlap(stList *overlappingPairs);

//Split over large gaps

stList *getSplitPoints(stList *anchorPairs, int64_t lX, int64_t lY,
        int64_t maxMatrixSize);

void getPosteriorProbsWithBandingSplittingAlignmentsByLargeGaps(stList *anchorPairs, const char *sX, const char *sY, int64_t lX, int64_t lY,
        PairwiseAlignmentParameters *p,  bool alignmentHasRaggedLeftEnd, bool alignmentHasRaggedRightEnd,
        void (*diagonalPosteriorProbFn)(int64_t, DpMatrix *, DpMatrix *, const SymbolString, const SymbolString,
                double, PairwiseAlignmentParameters *, void *),
        void (*coordinateCorrectionFn)(), void *extraArgs);

stList *filterPairwiseAlignmentToMakePairsOrdered(stList *alignedPairs, float gapGamma);


#endif /* PAIRWISEALIGNER_H_ */
