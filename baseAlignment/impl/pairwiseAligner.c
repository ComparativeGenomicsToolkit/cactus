#include "sonLib.h"
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

///
//Sequence stuff
//

static char convertChar(char i) {
    switch(i) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return 4;
    }
}

char *convertSequence(const char *s, int32_t sL) {
    assert(sL >= 0);
    assert(strlen(s) == sL);
    char *cS = st_malloc((sL+1) * sizeof(char));
    for(int32_t i=0; i<sL; i++) {
        cS[i] = convertChar(s[i]);
    }
    cS[sL] = '\0';
    return cS;
}

///
//Basic math
///

static const double logZero = -INFINITY;
static const double logUnderflowThreshold = 7.5;

static inline double lookup (double x){
  //return log (exp (x) + 1);
  assert (x >= 0.00f);
  assert (x <= logUnderflowThreshold);
  if (x <= 1.00f) return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
  if (x <= 2.50f) return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
  if (x <= 4.50f) return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
  assert (x <= logUnderflowThreshold);
  return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
}

double logAdd(double x, double y){
  if (x < y) return (x == logZero || y - x >= logUnderflowThreshold) ? y : lookup(y-x) + x;
  return (y == logZero || x - y >= logUnderflowThreshold) ? x : lookup(x-y) + y;
}

/*
 * Calculates the addition of two log numbers then assigns them to the pointer x.
 */
static inline void logAddAndAssign(double *x, double y) {
    *x = logAdd(*x, y);
}

////
//State stuff.
////

static const int32_t cellNo = 3;

static void checkState(int32_t state) {
    assert(state >= 0);
    assert(state < cellNo);
}

static void checkPosition(int32_t z, int32_t zL) {
    assert(z >= 1);
    assert(z < zL);
}

static const double posteriorMatchThreshold = 0.01;

#define matchContinueTransition -0.027602076970648346 //0.972775379521401f
#define gapExtendTransition -0.025886909285416447 //0.974445284091146f;
#define gapSwitchTransition -7.2203887919613203 //0.0007315179552849f;
#define matchFromGapTransition -3.6959766616728253 //1.0 - gapExtend - gapSwitch = 0.024823197953569152
#define gapOpenTransition -4.2967807310002062 //(1.0 - match)/2 = 0.013612310239299985

static inline double transitionProb(int32_t from, int32_t to) {
    checkState(from);
    checkState(to);
    static const double transitions[9] = { /*Match */ matchContinueTransition, matchFromGapTransition, matchFromGapTransition,
                /*To gapX */ gapOpenTransition, gapExtendTransition, gapSwitchTransition,
                /*To gapY */ gapOpenTransition, gapSwitchTransition, gapExtendTransition };
    return transitions[to * cellNo + from];
}

static inline int32_t getTransitionOffSetX(int32_t state) {
    checkState(state);
    static int32_t offsets[] = { 1, 1, 0 };
    return offsets[state];
}

static inline int32_t getTransitionOffSetY(int32_t state) {
    checkState(state);
    static int32_t offsets[] = { 1, 0, 1 };
    return offsets[state];
}

#define gapEmission -1.6094379124341003 //log(0.2) = -1.6094379124341003
#define matchEmission -2.1149196655034745 //log(0.12064298095701059);
#define transversionEmission -4.5691014376830479 //log(0.010367271172731285);
#define transitionEmission -3.9833860032220842 //log(0.01862247669752685);
#define matchNEmission -3.2188758248682006 //log(0.04);

static inline double emissionProb(int32_t x, int32_t y, int32_t lX, int32_t lY,
        const char *sX, const char *sY, int32_t state) {
    checkState(state);
    static const double gapM[5] = { gapEmission, gapEmission, gapEmission, gapEmission, gapEmission };
    static const double matchM[25] = { matchEmission, transversionEmission, transitionEmission, transversionEmission, matchNEmission,
            transversionEmission, matchEmission, transversionEmission, transitionEmission, matchNEmission,
            transitionEmission, transversionEmission, matchEmission, transversionEmission, matchNEmission,
            transversionEmission, transitionEmission, transversionEmission, matchEmission, matchNEmission,
            matchNEmission, matchNEmission, matchNEmission, matchNEmission, matchNEmission };
    switch(state) {
        case 0:
            checkPosition(x, lX);
            checkPosition(y, lY);
            return matchM[sX[x-1] * 5 + sY[y-1]];
        case 1:
            checkPosition(x, lX);
            return gapM[(int32_t)sX[x-1]];
        case 2:
            checkPosition(y, lY);
            return gapM[(int32_t)sY[y-1]];
        default:
            assert(0);
    }
}

static inline double startStateProbs(int32_t state) {
    checkState(state);
    //static const double startProb = -1.0986122886681098; //math.log(1.0/3.0) = -1.0986122886681098
    static const double startStates[3] = { matchContinueTransition, gapOpenTransition, gapOpenTransition };
    return startStates[state];
}

static inline double endStateProbs(int32_t state) {
    checkState(state);
    static const double endProb = -1.0986122886681098; //math.log(1.0/3.0) = -1.0986122886681098
    double endStates[3] = { endProb, endProb, endProb };
    return endStates[state];
}

/////
//Forward matrix
/////

static inline double *getCell(double *m, int32_t x, int32_t y, int32_t lX) {
    if(x >= 0 && y >= 0) {
        return &(m[(y * lX + x) * cellNo]);
    }
    return NULL;
}

static inline double *getEmptyMatrix(int32_t lX, int32_t lY) {
    int32_t j = lX * lY * cellNo;
    double *m = st_malloc(j * sizeof(double));
    for(int32_t i=0; i<j; i++) {
        m[i] = logZero;
    }
    return m;
}

static double *initialiseForwardMatrix(int32_t lX, int32_t lY) {
    double *fM = getEmptyMatrix(lX, lY);
    double *cell = getCell(fM, 0, 0, lX);
    assert(cell != NULL);
    for(int32_t i=0; i<cellNo; i++) {
        cell[i] = startStateProbs(i);
    }
    return fM;
}

static inline void forwardCell(double *fM, int32_t x, int32_t y, int32_t lX,
        int32_t lY, const char *sX, const char *sY) {
    double *cell = getCell(fM, x, y, lX);
    for (int32_t to = 0; to < cellNo; to++) {
        double *pCell = getCell(fM, x - getTransitionOffSetX(to), y
                - getTransitionOffSetY(to), lX);
        if (pCell != NULL) {
            double eP = emissionProb(x, y, lX, lY, sX, sY, to);
            for (int32_t from = 0; from < cellNo; from++) {
                logAddAndAssign(&cell[to], pCell[from] + transitionProb(from, to) + eP);
            }
        }
    }
}

double *forwardMatrix(int32_t lX, int32_t lY, const char *sX, const char *sY) {
    double *fM = initialiseForwardMatrix(lX, lY);

    for (int32_t x = 0; x < lX; x++) {
        for (int32_t y = 0; y < lY; y++) {
            forwardCell(fM, x, y, lX, lY, sX, sY);
        }
    }
    return fM;
}

double totalForwardProb(double *fM, int32_t lX, int32_t lY) {
    double *cell = getCell(fM, lX-1, lY-1, lX);
    assert(cell != NULL);
    double totalProb = endStateProbs(0) + cell[0];
    for(int32_t i=1; i<cellNo; i++) {
        logAddAndAssign(&totalProb, endStateProbs(i) + cell[i]);
    }
    return totalProb;
}

/////
//Backward matrix
/////

static double *initialiseBackwardMatrix(int32_t lX, int32_t lY) {
    double *bM = getEmptyMatrix(lX, lY);
    double *cell = getCell(bM, lX-1, lY-1, lX);
    assert(cell != NULL);
    for(int32_t i=0; i<cellNo; i++) {
        cell[i] = endStateProbs(i);
    }
    return bM;
}

static inline void backwardCell(double *bM, int32_t x, int32_t y, int32_t lX,
        int32_t lY, const char *sX, const char *sY) {
    double *cell = getCell(bM, x, y, lX);
    for (int32_t to = 0; to < cellNo; to++) {
        double *pCell = getCell(bM, x - getTransitionOffSetX(to), y
                - getTransitionOffSetY(to), lX);
        if (pCell != NULL) {
            double eP = emissionProb(x, y, lX, lY, sX, sY, to);
            for (int32_t from = 0; from < cellNo; from++) {
                logAddAndAssign(&pCell[from], cell[to] + transitionProb(from, to) + eP);
            }
        }
    }
}

double *backwardMatrix(int32_t lX, int32_t lY, const char *sX, const char *sY) {
    double *bM = initialiseBackwardMatrix(lX, lY);
    for (int32_t x = lX - 1; x >= 0; x--) {
        for (int32_t y = lY - 1; y >= 0; y--) {
            backwardCell(bM, x, y, lX, lY, sX, sY);
        }
    }
    return bM;
}

double totalBackwardProb(double *fM, int32_t lX) {
    double *cell = getCell(fM, 0, 0, lX);
    double totalProb = startStateProbs(0) + cell[0];
    for(int32_t i=1; i<cellNo; i++) {
        logAddAndAssign(&totalProb, startStateProbs(i) + cell[i]);
    }
    return totalProb;
}

/////
//Posterior probabilities
/////

static inline double posteriorMatchProb(double *fM, double *bM, int32_t x, int32_t y,
        int32_t lX, int32_t lY, const char *sX, const char *sY, double totalProb) {
    int32_t to = 0;
    double *pCell = getCell(fM, x - 1, y - 1, lX);
    assert(pCell != NULL);
    double *cell = getCell(bM, x, y, lX);
    double eP = emissionProb(x, y, lX, lY, sX, sY, to);
    int32_t from = 0;
    double f = pCell[from] + transitionProb(from, to) + eP + cell[to];
    for (from = 1; from < cellNo; from++) {
        logAddAndAssign(&f, pCell[from] + transitionProb(from, to) + eP + cell[to]);
    }
    double p = exp(f - totalProb);
    //assert(p >= -0.01 && p < 1.01);
    return p;
}

static void getPosteriorProbs(double *fM, double *bM, int32_t lX, int32_t lY,
        const char *sX, const char *sY, stList *alignedPairs,
        double totalProb) {
    for (int32_t x = 1; x < lX; x++) {
        for (int32_t y = 1; y < lY; y++) {
            double f = posteriorMatchProb(fM, bM, x, y, lX, lY, sX, sY, totalProb);
            if (f > posteriorMatchThreshold) {
                if(f > 1.0) {
                    f = 1.0;
                }
                stIntTuple *alignedPair = stIntTuple_construct(3, (int32_t)(f * 1000), x - 1, y - 1);
                stList_append(alignedPairs, alignedPair);
            }
        }
    }
}

///////
//Maximal expected accuracy alignment
///////

stList *getAlignedPairs(const char *sX, const char *sY, void *parameters) {
    assert(parameters != NULL);
    //Allocate the matrices.
    int32_t lX = strlen(sX)+1;
    int32_t lY = strlen(sY)+1;

    char *cSX = convertSequence(sX, lX-1);
    char *cSY = convertSequence(sY, lY-1);

    double *fM = forwardMatrix(lX, lY, cSX, cSY);
    double *bM = backwardMatrix(lX, lY, cSX, cSY);

    double totalFProb = totalForwardProb(fM, lX, lY);
    double totalBProb = totalBackwardProb(bM, lX);
    double totalProb = (totalFProb + totalBProb)/2;
    double diff = (totalFProb - totalBProb)  / totalProb;
    if(diff < 0.0) {
        diff = -diff;
    }
    //Check they are about the same.
    //st_uglyf("This is what I aligned %f %f %f %s %s\n", totalFProb, totalBProb, diff, sX, sY);
    assert(diff < 0.001);

    //Get posterior probabilities above 0.01 threshold.
    stList *alignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    getPosteriorProbs(fM, bM, lX, lY, cSX, cSY, alignedPairs, totalProb);

    //Cleanup
    free(fM);
    free(bM);
    free(cSX);
    free(cSY);

    return alignedPairs;
}
