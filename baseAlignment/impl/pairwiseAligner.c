/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
//This is being included to make popen work!
#define _XOPEN_SOURCE 500

#include "sonLib.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <ctype.h>
#include "pairwiseAligner.h"
#include "sonLib.h"
#include "pairwiseAlignment.h"

#ifdef _OPENMP
static const int32_t ThreadChunkSize = 150;
static inline int32_t getNumThreads(int32_t arraySize)
{
    int32_t numThreads = omp_get_max_threads();
    int32_t maxThreads = arraySize / ThreadChunkSize;
    if (maxThreads < numThreads && maxThreads > 0)
    return maxThreads;
    return numThreads;
}
#endif
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
    return ((-0.000458661602210f * x + 0.009695946122598f) * x
            + 0.930734667215156f) * x + 0.168037164329057f;
}

double logAdd(double x, double y) {
    if (x < y)
        return (x == logZero || y - x >= logUnderflowThreshold) ? y : lookup(
                y - x) + x;
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
#define matchFromLongGapTransition -5.673280173170473 //1.0 - gapExtend = 0.00343657420938
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
            return 0.0;
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
        double *pCell = getCell(fM, x - getTransitionOffSetX(to),
                y - getTransitionOffSetY(to), lX);
        if (pCell != NULL) {
            double eP = emissionProb(x, y, lX, lY, sX, sY, to);
            for (int32_t from = 0; from < cellNo; from++) {
                double tP = transitionProb(from, to);
                if (tP != logZero) {
                    logAddAndAssign(&cell[to], pCell[from] + tP + eP);
                }
            }
        }
    }
}

double *forwardMatrix(int32_t lX, int32_t lY, const char *sX, const char *sY) {
    double *fM = initialiseForwardMatrix(lX, lY);

    // all '/' diagonals intersectinog x = 0 axis
    for (int32_t x = 0; x < lX; ++x) {
        int32_t diagLen = (lY < x + 1) ? lY : x + 1;
#ifdef _OPENMP
#pragma omp parallel default(shared) if(diagLen >= ThreadChunkSize) num_threads(getNumThreads(diagLen))
#endif
        {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int32_t y = 0; y < diagLen; ++y) {
                forwardCell(fM, x - y, y, lX, lY, sX, sY);
            }
        }
    }

    // all '/' diagonals intersecting y = lY - 1 axis
    for (int32_t y = 1; y < lY; ++y) {
        int32_t diagLen = (lX < lY - y) ? lX : lY - y;
#ifdef _OPENMP
#pragma omp parallel default(shared) if(diagLen >= ThreadChunkSize) num_threads(getNumThreads(diagLen))
#endif
        {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int32_t i = 0; i < diagLen; ++i) {
                forwardCell(fM, lX - 1 - i, y + i, lX, lY, sX, sY);
            }
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
        int32_t x2 = x + getTransitionOffSetX(to);
        int32_t y2 = y + getTransitionOffSetY(to);
        if (x2 < lX && y2 < lY) {
            double *pCell = getCell(bM, x2, y2, lX);
            double eP = emissionProb(x2, y2, lX, lY, sX, sY, to);
            for (int32_t from = 0; from < cellNo; from++) {
                double tP = transitionProb(from, to);
                if (tP != LOG_ZERO) {
                    logAddAndAssign(&cell[from], pCell[to] + tP + eP);
                }
            }
        }
    }
}

double *backwardMatrix(int32_t lX, int32_t lY, const char *sX, const char *sY) {
    double *bM = initialiseBackwardMatrix(lX, lY);

    // all '/' diagonals intersecting y = lY - 1 axis
    for (int32_t y = lY - 1; y >= 1; --y) {
        int32_t diagLen = (lX < lY - y) ? lX : lY - y;
#ifdef _OPENMP
#pragma omp parallel default(shared) if(diagLen >= ThreadChunkSize) num_threads(getNumThreads(diagLen))
#endif
        {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int32_t i = 0; i < diagLen; ++i) {
                backwardCell(bM, lX - 1 - i, y + i, lX, lY, sX, sY);
            }
        }
    }

    // all '/' diagonals intersecting x = 0 axis
    for (int32_t x = lX - 1; x >= 0; --x) {
        int32_t diagLen = (lY < x + 1) ? lY : x + 1;
#ifdef _OPENMP
#pragma omp parallel default(shared) if(diagLen >= ThreadChunkSize) num_threads(getNumThreads(diagLen))
#endif
        {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int32_t y = 0; y < diagLen; ++y) {
                backwardCell(bM, x - y, y, lX, lY, sX, sY);
            }
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
        logAddAndAssign(&f,
                pCell[from] + transitionProb(from, to) + eP + cell[to]);
    }
    double p = exp(f - totalProb);
#ifdef BEN_DEBUG
    /*if(p < -0.01 || p > 1.01) {
     st_uglyf("I got a bad position, %i %i %f\n", x, y, p);
     }*/
    //assert(p >= -0.01 && p < 1.01);
#endif
    return p;
}

static void getPosteriorProbs(double *fM, double *bM, int32_t lX, int32_t lY,
        const char *sX, const char *sY, stList *alignedPairs, double totalProb,
        PairwiseAlignmentParameters *p) {
    stList** alignedPairMatrix =
            (stList**) st_malloc(lX * lY * sizeof(stList*));
#ifdef _OPENMP
#pragma omp parallel default(shared) if(lX >= ThreadChunkSize) num_threads(getNumThreads(lX))
#endif
    {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int32_t x = 1; x < lX; x++) {
            alignedPairMatrix[x] = stList_construct();
            for (int32_t y = 1; y < lY; y++) {
                double f = posteriorMatchProb(fM, bM, x, y, lX, lY, sX, sY,
                        totalProb);
                if (f >= posteriorMatchThreshold) {
                    if (f > 1.0) {
                        f = 1.0;
                    }
#ifdef BEN_DEBUG
                    assert(sX[x-1] >= 0 && sX[x-1] <= 4);
                    assert(sY[y-1] >= 0 && sY[y-1] <= 4);
#endif
                    if (p->alignAmbiguityCharacters || (sX[x - 1] < 4 && sY[y
                            - 1] < 4)) {
                        stIntTuple *alignedPair = stIntTuple_construct(3,
                                (int32_t) floor(f * PAIR_ALIGNMENT_PROB_1),
                                x - 1, y - 1);
                        stList_append(alignedPairMatrix[x], alignedPair);
                    }
                }
            }
        }
    }
    for (int32_t x = 1; x < lX; x++) {
        stList_appendAll(alignedPairs, alignedPairMatrix[x]);
        stList_destruct(alignedPairMatrix[x]);
    }
    free(alignedPairMatrix);
}

///////
//Maximal expected accuracy alignment
///////

stList *getAlignedPairs(const char *sX, const char *sY,
        PairwiseAlignmentParameters *p) {
    //Allocate the matrices.
    int32_t lX = strlen(sX) + 1;
    int32_t lY = strlen(sY) + 1;

    char *cSX = convertSequence(sX, lX - 1);
    char *cSY = convertSequence(sY, lY - 1);

    double *fM;
    double *bM;

    //#pragma omp parallel default(shared)
    {
        //#pragma omp sections
        {
            //#pragma omp section
            {
                bM = backwardMatrix(lX, lY, cSX, cSY);
            }
            //#pragma omp section
            {
                fM = forwardMatrix(lX, lY, cSX, cSY);
            }
        }
    }

    double totalFProb = totalForwardProb(fM, lX, lY);
    double totalBProb = totalBackwardProb(bM, lX);
    double totalProb = (totalFProb + totalBProb) / 2;
    double diff = (totalFProb - totalBProb) / totalProb;
    if (diff < 0.0) {
        diff = -diff;
    }
    //Check they are about the same.
#ifdef BEN_DEBUG
    //assert(diff < 0.0001);
#endif

    //Get posterior probabilities above 0.01 threshold.
    stList *alignedPairs = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    getPosteriorProbs(fM, bM, lX, lY, cSX, cSY, alignedPairs, totalProb, p);

    //Cleanup
    free(fM);
    free(bM);
    free(cSX);
    free(cSY);

    return alignedPairs;
}

char *getSubString(const char *cA, int32_t start, int32_t length) {
    char *cA2 = memcpy(st_malloc(sizeof(char) * (length + 1)), cA + start,
            length);
    cA2[length] = '\0';

#ifdef BEN_DEBUG
    for(int32_t i=0; i<length; i++) {
        assert(cA2[i] == cA[i + start]);
        assert(cA2[i] != '\0');
    }
    assert(cA2[length] == '\0');
#endif
    return cA2;
}

static int sortByXPlusYCoordinate(const void *i, const void *j) {
    int64_t k = stIntTuple_getPosition((stIntTuple *) i, 0)
            + stIntTuple_getPosition((stIntTuple *) i, 1);
    int64_t l = stIntTuple_getPosition((stIntTuple *) j, 0)
            + stIntTuple_getPosition((stIntTuple *) j, 1);
    return k > l ? 1 : (k < l ? -1 : 0);
}

static char *makeUpperCase(const char *s, int32_t l) {
    char *s2 = stString_copy(s);
    for(int32_t i=0; i<l; i++) {
        s2[i] = toupper(s[i]);
    }
    return s2;
}

stList *getBlastPairs(const char *sX, const char *sY, int32_t lX, int32_t lY,
        int32_t trim, bool repeatMask) {
    /*
     * Uses lastz to compute a bunch of monotonically increasing pairs such that for any pair of consecutive pairs in the list
     * (x1, y1) (x2, y2) in the set of aligned pairs x1 appears before x2 in X and y1 appears before y2 in Y.
     */
    stList *alignedPairs = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct); //the list to put the output in

    if (lX == 0 || lY == 0) {
        return alignedPairs;
    }

    if (!repeatMask) {
        sX = makeUpperCase(sX, lX);
        sY = makeUpperCase(sY, lY);
    }

    //Write one sequence to file..
    char *tempFile1 = getTempFile();

    FILE *fileHandle = fopen(tempFile1, "w");
    fprintf(fileHandle, ">a\n%s\n", sX);
    fclose(fileHandle);

    char
            *command =
                    stString_print(
                            "echo '>b\n%s\n' | lastz --hspthresh=800 --chain --strand=plus --gapped --format=cigar --ambiguous=iupac %s",
                            sY, tempFile1);
    fileHandle = popen(command, "r");
    free(command);
    if (fileHandle == NULL) {
        st_errAbort("Problems with lastz pipe");
    }
    //Read from stream
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
                    stList_append(alignedPairs,
                            stIntTuple_construct(2, j + l, k + l));
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
        st_errnoAbort(
                "pclose failed when getting rid of lastz pipe with value %i",
                status);
    }

    stList_sort(alignedPairs, sortByXPlusYCoordinate); //Ensure the coordinates are increasing

    //Remove old files
    st_system("rm %s", tempFile1);
    free(tempFile1);

    if (!repeatMask) {
        free((char *) sX);
        free((char *) sY);
    }

    return alignedPairs;
}

stList *getAnchorPoints(stList *alignedPairs, int32_t minRectangleSize,
        int32_t lX, int32_t lY) {
    /*
     * Filters the blast pairs so that the are spaced with dp matrices of at least minRectangle size between them.
     */
    int32_t x = -1;
    int32_t y = -1;
    stList *filteredPairs = stList_construct();
    for (int32_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int32_t x2 = stIntTuple_getPosition(alignedPair, 0);
        int32_t y2 = stIntTuple_getPosition(alignedPair, 1);
        if (x2 > x && y2 > y) { //This should never occur but deals with an error in lastz, I think....
            int64_t lRectangleSize = ((int64_t) x2 - x + 1) * (y2 - y + 1);
            int64_t rRectangeSize = ((int64_t) lX - x2) * (lY - y2);
            if (lRectangleSize >= minRectangleSize && rRectangeSize
                    >= minRectangleSize) {
                stList_append(filteredPairs, alignedPair);
                x = x2;
                y = y2;
            }
        }
    }
    assert(x < lX);
    assert(y < lY);
    return filteredPairs;
}

static void convertPairs(stList *alignedPairs2, int32_t offsetX,
        int32_t offsetY) {
    /*
     * Convert the coordinates of the computed pairs.
     */
    for (int32_t k = 0; k < stList_length(alignedPairs2); k++) {
        stIntTuple *i = stList_get(alignedPairs2, k);
        assert(stIntTuple_length(i) == 3);
        stList_set(
                alignedPairs2,
                k,
                stIntTuple_construct(3, stIntTuple_getPosition(i, 0),
                        stIntTuple_getPosition(i, 1) + offsetX,
                        stIntTuple_getPosition(i, 2) + offsetY));
        stIntTuple_destruct(i);
    }
}

static int32_t splitSequence(char *s, int32_t l, int32_t bandSize, char **sL,
        char **sR) {
    if (l > bandSize) {
        *sL = getSubString(s, 0, bandSize);
        *sR = getSubString(s, l - bandSize, bandSize);
        return l - bandSize;
    }
    *sL = s;
    *sR = s;
    return 0;
}

static stList *getAlignedPairs_Split(char *sX, char *sY, int32_t lX,
        int32_t lY, int32_t bandSize, PairwiseAlignmentParameters *p) {
    /*
     * Aligns the sequences, but if the product of there sequence lengths is greater than bandSize squared * 2
     * then a dp matrix of bandsize squared if computed in the top left part of the entire dp matrix, and a
     * corresponding square in the bottom right part of the matrix.
     */
    if ((int64_t) lX * lY <= (int64_t) bandSize * bandSize * 2) { //products can be > 2^31
        return getAlignedPairs(sX, sY, p);
    }
    st_logDebug("We found an overlarge matrix to compute: %i %i \n", lX, lY);
    char *sXL, *sXR;
    int32_t offsetX = splitSequence(sX, lX, bandSize, &sXL, &sXR);
    char *sYL, *sYR;
    int32_t offsetY = splitSequence(sY, lY, bandSize, &sYL, &sYR);
    stList *alignedPairsL = getAlignedPairs(sXL, sYL, p);
    stList *alignedPairsR = getAlignedPairs(sXR, sYR, p);
    convertPairs(alignedPairsR, offsetX, offsetY);
    assert(offsetX != 0 || offsetY != 0);
    if (offsetX != 0) {
        free(sXL);
        free(sXR);
    }
    if (offsetY != 0) {
        free(sYL);
        free(sYR);
    }
    /*
     * Get a unique set of pairs.
     */
    stSortedSet *pairs = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(alignedPairsL); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairsL, i);
        stIntTuple *coordinatePair = stIntTuple_construct(2,
                stIntTuple_getPosition(alignedPair, 1),
                stIntTuple_getPosition(alignedPair, 2));
        assert(stSortedSet_search(pairs, coordinatePair) == NULL);
        stSortedSet_insert(pairs, coordinatePair);
    }
    while (stList_length(alignedPairsR) > 0) { //empty and destroy the second list.
        stIntTuple *alignedPair = stList_pop(alignedPairsR);
        stIntTuple *coordinatePair = stIntTuple_construct(2,
                stIntTuple_getPosition(alignedPair, 1),
                stIntTuple_getPosition(alignedPair, 2));
        if (stSortedSet_search(pairs, coordinatePair) == NULL) {
            stSortedSet_insert(pairs, coordinatePair);
            stList_append(alignedPairsL, alignedPair);
        } else {
            stIntTuple_destruct(alignedPair);
            stIntTuple_destruct(coordinatePair);
        }
    }
    stSortedSet_destruct(pairs);
    stList_appendAll(alignedPairsL, alignedPairsR);
    stList_destruct(alignedPairsR);

    return alignedPairsL;
}

PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters_construct() {
    PairwiseAlignmentParameters *p = st_malloc(
            sizeof(PairwiseAlignmentParameters));
    p->maxBandingSize = 3000;
    p->minBandingSize = 500;
    p->minBandingConstraintDistance = 100;
    p->minTraceBackDiag = 50;
    p->constraintDiagonalTrim = 4;
    p->alignAmbiguityCharacters = 0;
    return p;
}

void pairwiseAlignmentBandingParameters_destruct(PairwiseAlignmentParameters *p) {
    free(p);
}

static int32_t boundCoordinateTransform(int32_t i, int32_t j, int32_t max) {
    i += j;
    return i >= max ? max - 1 : (i < 0 ? 0 : i);
}

stList *getAlignedPairs_FastP(const char *sX, const char *sY,
        PairwiseAlignmentParameters *p, bool recursive) {
    /*
     * Aligns the pair of sequences using banding constraints.
     */
    int32_t lX = strlen(sX);
    int32_t lY = strlen(sY);
    stList *alignedPairs = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    if (lX == 0 || lY == 0) {
        return alignedPairs;
    }

    int64_t minBandingSquare = (int64_t) p->minBandingSize * p->minBandingSize;

    stList *blastPairs;
    if ((int64_t) lX * lY > minBandingSquare) {
        blastPairs = getBlastPairs(sX, sY, lX, lY, p->constraintDiagonalTrim,
                recursive);
    } else {
        blastPairs = stList_construct(); //We don't bother getting anchors if the sequences are sufficiently small.
    }

    stList *bandPairs = getAnchorPoints(blastPairs,
            p->minBandingConstraintDistance * p->minBandingConstraintDistance,
            lX, lY);
    st_logDebug("We got %i aligned pairs and %i filtered pairs\n",
            stList_length(blastPairs), stList_length(bandPairs));
    stIntTuple *finalPair = stIntTuple_construct(2, lX - 1, lY - 1);
    stList_append(bandPairs, finalPair);
    stListIterator *bandIt = stList_getIterator(bandPairs);
    stIntTuple *bandPair;

    int32_t x1 = 0, x2 = -1;
    int32_t y1 = 0, y2 = -1;
    int32_t diagOffset = p->minTraceBackDiag / 2;
    while ((bandPair = stList_getNext(bandIt)) != NULL) {
        int32_t x4 = stIntTuple_getPosition(bandPair, 0);
        int32_t y4 = stIntTuple_getPosition(bandPair, 1);
        int32_t x3 = boundCoordinateTransform(x4, -diagOffset, lX);
        int32_t y3 = boundCoordinateTransform(y4, -diagOffset, lY);
        int32_t x5 = boundCoordinateTransform(x4, diagOffset, lX);
        int32_t y5 = boundCoordinateTransform(y4, diagOffset, lY);
        int32_t lX2 = x5 - x1 + 1;
        int32_t lY2 = y5 - y1 + 1;

        //st_logDebug(
        //        "The next blast square, x1: %i, x2: %i, x3: %i, x4: %i, x5: %i, y1: %i, y2: %i, y3: %i, y4: %i, y5: %i size x: %i, size y: %i\n",
        //        x1, x2, x3, x4, x5, y1, y2, y3, y4, y5, lX2, lY2);

        assert(x1 <= x5);
        assert(y1 <= y5);
        assert(x3 <= x4 && x4 <= x5);
        assert(y3 <= y4 && y4 <= y5);

        //Get the appropriate x substring
        char *sX2 = getSubString(sX, x1, lX2);

        //Get the appropriate y substring
        char *sY2 = getSubString(sY, y1, lY2);

        //Do the actual alignment.
        stList *alignedPairs2;

        if ((int64_t) lX2 * lY2 <= minBandingSquare) { //products can be > 2^31
            alignedPairs2 = getAlignedPairs(sX2, sY2, p);
        }
        else if (recursive) {
            alignedPairs2 = getAlignedPairs_FastP(sX2, sY2, p, 0);
        } else {
            alignedPairs2 = getAlignedPairs_Split(sX2, sY2, lX2, lY2,
                    p->maxBandingSize, p);
        }

        //Cleanup the temporary sequences
        free(sX2);
        free(sY2);

        //Convert the coordinates
        convertPairs(alignedPairs2, x1, y1);

        int32_t minDiag = x2 + y2;
        int32_t maxDiag = x4 + y4;
        assert(minDiag < maxDiag);

        //Add the pairs to the alignment (merging together any duplicate pairs)
        //And skip any pairs within minTraceGapDiags.
        while (stList_length(alignedPairs2) > 0) {
            stIntTuple *alignedPair = stList_pop(alignedPairs2);
            assert(stIntTuple_length(alignedPair) == 3);
            int32_t xC = stIntTuple_getPosition(alignedPair, 1);
            int32_t yC = stIntTuple_getPosition(alignedPair, 2);
            assert(xC >= 0 && xC < lX);
            assert(yC >= 0 && yC < lY);
            int32_t x_y = xC + yC;
            if (x_y > minDiag && x_y <= maxDiag) {
                stList_append(alignedPairs, alignedPair);
            } else {
                stIntTuple_destruct(alignedPair);
            }
        }
        stList_destruct(alignedPairs2);

        x1 = x3;
        x2 = x4;
        y1 = y3;
        y2 = y4;
    }
    stList_destructIterator(bandIt);
    assert(x2 == lX - 1);
    assert(y2 == lY - 1);
    stList_destruct(bandPairs);
    stList_destruct(blastPairs);
    stIntTuple_destruct(finalPair);

#ifdef BEN_DEBUG
    //Checks no coordinate pairs are unique in list.
    stSortedSet *pairs = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, (void (*)(void *))stIntTuple_destruct);
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        stIntTuple *coordinatePair = stIntTuple_construct(2, stIntTuple_getPosition(alignedPair, 1), stIntTuple_getPosition(alignedPair, 2));
        assert(stSortedSet_search(pairs, coordinatePair) == NULL);
        stSortedSet_insert(pairs, coordinatePair);
    }
    stSortedSet_destruct(pairs);
#endif

    st_logDebug("Finished creating aligned pairs, got %i pairs\n", stList_length(alignedPairs));
    return alignedPairs;
}

stList *getAlignedPairs_Fast(const char *sX, const char *sY,
        PairwiseAlignmentParameters *p) {
    return getAlignedPairs_FastP(sX, sY, p, 1);
}
