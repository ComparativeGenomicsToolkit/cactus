#include "sonLib.h"
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include "pairwiseAligner.h"

/*
 * Basic Pecan HMM code used by Cactus base aligner.
 */

///
//Sequence stuff
//

static char convertChar(char i) {
    switch (i) {
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
    char *cS = st_malloc((sL + 1) * sizeof(char));
    for (int32_t i = 0; i < sL; i++) {
        cS[i] = convertChar(s[i]);
    }
    cS[sL] = '\0';
    return cS;
}

///
//Basic math
///

#define logZero -INFINITY
#define logUnderflowThreshold 7.5
#define posteriorMatchThreshold 0.01

static inline double lookup(double x) {
    //return log (exp (x) + 1);
#ifdef BEN_DEBUG
    assert (x >= 0.00f);
    assert (x <= logUnderflowThreshold);
#endif
    if (x <= 1.00f)
        return ((-0.009350833524763f * x + 0.130659527668286f) * x
                + 0.498799810682272f) * x + 0.693203116424741f;
    if (x <= 2.50f)
        return ((-0.014532321752540f * x + 0.139942324101744f) * x
                + 0.495635523139337f) * x + 0.692140569840976f;
    if (x <= 4.50f)
        return ((-0.004605031767994f * x + 0.063427417320019f) * x
                + 0.695956496475118f) * x + 0.514272634594009f;
#ifdef BEN_DEBUG
    assert (x <= logUnderflowThreshold);
#endif
    return ((-0.000458661602210f * x + 0.009695946122598f) * x
            + 0.930734667215156f) * x + 0.168037164329057f;
}

double logAdd(double x, double y) {
    if (x < y)
        return (x == logZero || y - x >= logUnderflowThreshold) ? y : lookup(y
                - x) + x;
    return (y == logZero || x - y >= logUnderflowThreshold) ? x : lookup(x - y)
            + y;
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

static const int32_t cellNo = 5;

static void checkState(int32_t state) {
#ifdef BEN_DEBUG
    assert(state >= 0);
    assert(state < cellNo);
#endif
}

static void checkPosition(int32_t z, int32_t zL) {
#ifdef BEN_DEBUG
    assert(z >= 1);
    assert(z < zL);
#endif
}

#define matchContinueTransition -0.030064059121770816 //0.9703833696510062f
#define gapOpenShortTransition -4.34381910900448 //0.0129868352330243
#define gapOpenLongTransition -6.30810595366929 //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
#define gapShortExtendTransition -0.3388262689231553 //0.7126062401851738f;
#define gapShortSwitchTransition -4.910694825551255 //0.0073673675173412815f;
#define matchFromShortGapTransition -1.272871422049609 //1.0 - gapExtend - gapSwitch = 0.280026392297485
#define gapLongExtendTransition -0.003442492794189331 //0.99656342579062f;
#define matchFromLongGapTransition -5.673280173170473 //1.0 - gapExtend - gapSwitch = 0.00343657420938
static inline double transitionProb(int32_t from, int32_t to) {
    checkState(from);
    checkState(to);
    static const double transitions[25] = { /*Match */matchContinueTransition,
            matchFromShortGapTransition, matchFromShortGapTransition,
            matchFromLongGapTransition, matchFromLongGapTransition,
            /*To shortGapX */gapOpenShortTransition, gapShortExtendTransition,
            gapShortSwitchTransition, logZero, logZero,
            /*To shortGapY */gapOpenShortTransition, gapShortSwitchTransition,
            gapShortExtendTransition, logZero, logZero,
            /*To longGapX */gapOpenLongTransition, logZero, logZero,
            gapLongExtendTransition, logZero,
            /*To longGapY */gapOpenLongTransition, logZero, logZero, logZero,
            gapLongExtendTransition };
    return transitions[to * cellNo + from];
}

static inline int32_t getTransitionOffSetX(int32_t state) {
    checkState(state);
    static int32_t offsets[] = { 1, 1, 0, 1, 0 };
    return offsets[state];
}

static inline int32_t getTransitionOffSetY(int32_t state) {
    checkState(state);
    static int32_t offsets[] = { 1, 0, 1, 0, 1 };
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
    static const double gapM[5] = { gapEmission, gapEmission, gapEmission,
            gapEmission, gapEmission };
    static const double matchM[25] = { matchEmission, transversionEmission,
            transitionEmission, transversionEmission, matchNEmission,
            transversionEmission, matchEmission, transversionEmission,
            transitionEmission, matchNEmission, transitionEmission,
            transversionEmission, matchEmission, transversionEmission,
            matchNEmission, transversionEmission, transitionEmission,
            transversionEmission, matchEmission, matchNEmission,
            matchNEmission, matchNEmission, matchNEmission, matchNEmission,
            matchNEmission };
    switch (state) {
        case 0:
            checkPosition(x, lX);
            checkPosition(y, lY);
            return matchM[sX[x - 1] * 5 + sY[y - 1]];
        case 1:
        case 3:
            checkPosition(x, lX);
            return gapM[(int32_t) sX[x - 1]];
        case 2:
        case 4:
            checkPosition(y, lY);
            return gapM[(int32_t) sY[y - 1]];
        default:
            assert(0);
    }
}

static inline double startStateProbs(int32_t state) {
    checkState(state);
    //static const double startProb = -1.0986122886681098; //math.log(1.0/3.0) = -1.0986122886681098
    static const double startStates[5] = { matchContinueTransition,
            gapOpenShortTransition, gapOpenShortTransition,
            gapOpenLongTransition, gapOpenLongTransition };
    return startStates[state];
}

static inline double endStateProbs(int32_t state) {
    checkState(state);
    static const double endProb = -1.6094379124341; //math.log(1.0/5.0) = -1.6094379124341
    double endStates[5] = { endProb, endProb, endProb, endProb, endProb };
    return endStates[state];
}

/////
//Forward matrix
/////

static inline double *getCell(double *m, int32_t x, int32_t y, int32_t lX) {
    if (x >= 0 && y >= 0) {
        return &(m[(y * lX + x) * cellNo]);
    }
    return NULL;
}

static inline double *getEmptyMatrix(int32_t lX, int32_t lY) {
    int32_t j = lX * lY * cellNo;
    double *m = st_malloc(j * sizeof(double));
    for (int32_t i = 0; i < j; i++) {
        m[i] = logZero;
    }
    return m;
}

static double *initialiseForwardMatrix(int32_t lX, int32_t lY) {
    double *fM = getEmptyMatrix(lX, lY);
    double *cell = getCell(fM, 0, 0, lX);
    assert(cell != NULL);
    for (int32_t i = 0; i < cellNo; i++) {
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
                logAddAndAssign(&cell[to], pCell[from] + transitionProb(from,
                        to) + eP);
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
    double *cell = getCell(fM, lX - 1, lY - 1, lX);
    assert(cell != NULL);
    double totalProb = endStateProbs(0) + cell[0];
    for (int32_t i = 1; i < cellNo; i++) {
        logAddAndAssign(&totalProb, endStateProbs(i) + cell[i]);
    }
    return totalProb;
}

/////
//Backward matrix
/////

static double *initialiseBackwardMatrix(int32_t lX, int32_t lY) {
    double *bM = getEmptyMatrix(lX, lY);
    double *cell = getCell(bM, lX - 1, lY - 1, lX);
    assert(cell != NULL);
    for (int32_t i = 0; i < cellNo; i++) {
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
                logAddAndAssign(&pCell[from], cell[to] + transitionProb(from,
                        to) + eP);
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
    for (int32_t i = 1; i < cellNo; i++) {
        logAddAndAssign(&totalProb, startStateProbs(i) + cell[i]);
    }
    return totalProb;
}

/////
//Posterior probabilities
/////

static inline double posteriorMatchProb(double *fM, double *bM, int32_t x,
        int32_t y, int32_t lX, int32_t lY, const char *sX, const char *sY,
        double totalProb) {
    int32_t to = 0;
    double *pCell = getCell(fM, x - 1, y - 1, lX);
    assert(pCell != NULL);
    double *cell = getCell(bM, x, y, lX);
    double eP = emissionProb(x, y, lX, lY, sX, sY, to);
    int32_t from = 0;
    double f = pCell[from] + transitionProb(from, to) + eP + cell[to];
    for (from = 1; from < cellNo; from++) {
        logAddAndAssign(&f, pCell[from] + transitionProb(from, to) + eP
                + cell[to]);
    }
    double p = exp(f - totalProb);
    //assert(p >= -0.01 && p < 1.01);
    return p;
}

static void getPosteriorProbs(double *fM, double *bM, int32_t lX, int32_t lY,
        const char *sX, const char *sY, stList *alignedPairs, double totalProb) {
    for (int32_t x = 1; x < lX; x++) {
        for (int32_t y = 1; y < lY; y++) {
            double f = posteriorMatchProb(fM, bM, x, y, lX, lY, sX, sY,
                    totalProb);
            if (f > posteriorMatchThreshold) {
                if (f > 1.0) {
                    f = 1.0;
                }
                stIntTuple *alignedPair = stIntTuple_construct(3, (int32_t)(f
                        * PAIR_ALIGNMENT_PROB_1), x - 1, y - 1);
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
    int32_t lX = strlen(sX) + 1;
    int32_t lY = strlen(sY) + 1;

    char *cSX = convertSequence(sX, lX - 1);
    char *cSY = convertSequence(sY, lY - 1);

    double *fM = forwardMatrix(lX, lY, cSX, cSY);
    double *bM = backwardMatrix(lX, lY, cSX, cSY);

    double totalFProb = totalForwardProb(fM, lX, lY);
    double totalBProb = totalBackwardProb(bM, lX);
    double totalProb = (totalFProb + totalBProb) / 2;
    double diff = (totalFProb - totalBProb) / totalProb;
    if (diff < 0.0) {
        diff = -diff;
    }
    //Check they are about the same.
    //st_uglyf("This is what I aligned %f %f %f %s %s\n", totalFProb, totalBProb, diff, sX, sY);
#ifdef BEN_DEBUG
    assert(diff < 0.001);
#endif

    //Get posterior probabilities above 0.01 threshold.
    stList *alignedPairs = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    getPosteriorProbs(fM, bM, lX, lY, cSX, cSY, alignedPairs, totalProb);

    //Cleanup
    free(fM);
    free(bM);
    free(cSX);
    free(cSY);

    return alignedPairs;
}

stIntTuple *convertTuple(stIntTuple *i, int32_t offsetX, int32_t offsetY) {
    stIntTuple *j = stIntTuple_construct(3, stIntTuple_getPosition(i, 0), stIntTuple_getPosition(i, 1) + offsetX, stIntTuple_getPosition(i, 2) + offsetY);
    stIntTuple_destruct(i);
    return j;
}

stList *getAlignedPairs2(const char *sX, const char *sY, void *parameters) {
    assert(parameters != NULL);

    int32_t lX = strlen(sX) + 1;
    int32_t lY = strlen(sY) + 1;
    int32_t offsetX = 0;
    int32_t offsetY = 0;
    int32_t maxLength = 200;

    stList *alignedPairs = stList_construct3(0,
                (void(*)(void *)) stIntTuple_destruct);

    while (1) {
        //Get the appropriate prefix.
        char *sX2 = sX;
        int32_t lX2 = lX;
        if (lX2 >= maxLength) {
            sX2 = getPrefixString(sX, maxLength);
            lX2 = maxLength;
        }
        //Second prefix
        char *sY2 = sY;
        int32_t lY2 = lY;
        if (lY2 >= maxLength) {
            sY2 = getPrefixString(sY, maxLength);
            lY2 = maxLength;
        }

        stList *alignedPairs2 = getAlignedPairs(sX2, sY2, parameters);

        //Cleanup the temporary sequences
        free(sX2);
        free(sY2);

        if(lX2 < lX || lY2 < lY) { //We still have work to do on another round
            while(stList_length(alignedPairs2) > 0) {
                stIntTuple *i = convertTuple(stList_pop(alignedPairs2), offsetX, offsetY);
                int32_t diagC = stIntTuple_getPosition(i, 1) + stIntTuple_getPosition(i, 2);
                if() {

                }
                else {
                    stList_append(alignedPairs, i);
                }
            }
            stList_destruct(alignedPairs2);
        }
        else {
            while(stList_length(alignedPairs2) > 0) {
                stList_append(alignedPairs, convertTuple(stList_pop(alignedPairs2), offsetX, offsetY));
            }
            stList_destruct(alignedPairs2);
            break;
        }
    }

    return alignedPairs;
}
