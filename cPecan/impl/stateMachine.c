/*
 * stateMachine.c
 *
 *  Created on: 1 Aug 2014
 *      Author: benedictpaten
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>

#include "bioioC.h"
#include "sonLib.h"
#include "pairwiseAligner.h"

///////////////////////////////////
///////////////////////////////////
//Em training objects.
///////////////////////////////////
///////////////////////////////////

Hmm *hmm_constructEmpty(double pseudoExpectation, StateMachineType type) {
    Hmm *hmm = st_malloc(sizeof(Hmm));
    hmm->type = type;
    switch (type) {
    case fiveState:
    case fiveStateAsymmetric:
        hmm->stateNumber = 5;
        break;
    case threeState:
    case threeStateAsymmetric:
        hmm->stateNumber = 3;
        break;
    default:
        st_errAbort("Unrecognised state type: %i\n", type);
    }
    hmm->transitions = st_malloc(hmm->stateNumber * hmm->stateNumber * sizeof(double));
    for (int64_t i = 0; i < hmm->stateNumber * hmm->stateNumber; i++) {
        hmm->transitions[i] = pseudoExpectation;
    }
    hmm->emissions = st_malloc(hmm->stateNumber * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N * sizeof(double));
    for (int64_t i = 0; i < hmm->stateNumber * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N; i++) {
        hmm->emissions[i] = pseudoExpectation;
    }
    hmm->likelihood = 0.0;
    return hmm;
}

void hmm_destruct(Hmm *hmm) {
    free(hmm->transitions);
    free(hmm->emissions);
    free(hmm);
}

static inline double *hmm_getTransition2(Hmm *hmm, int64_t from, int64_t to) {
    return &(hmm->transitions[from * hmm->stateNumber + to]);
}

double hmm_getTransition(Hmm *hmm, int64_t from, int64_t to) {
    return *hmm_getTransition2(hmm, from, to);
}

void hmm_addToTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    *hmm_getTransition2(hmm, from, to) += p;
}

void hmm_setTransition(Hmm *hmm, int64_t from, int64_t to, double p) {
    *hmm_getTransition2(hmm, from, to) = p;
}

static inline double *hmm_getEmissionsExpectation2(Hmm *hmm, int64_t state, Symbol x, Symbol y) {
    return &(hmm->emissions[state * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N + x * SYMBOL_NUMBER_NO_N + y]);
}

double hmm_getEmissionsExpectation(Hmm *hmm, int64_t state, Symbol x, Symbol y) {
    return *hmm_getEmissionsExpectation2(hmm, state, x, y);
}

void hmm_addToEmissionsExpectation(Hmm *hmm, int64_t state, Symbol x, Symbol y, double p) {
    *hmm_getEmissionsExpectation2(hmm, state, x, y) += p;
}

void hmm_setEmissionsExpectation(Hmm *hmm, int64_t state, Symbol x, Symbol y, double p) {
    *hmm_getEmissionsExpectation2(hmm, state, x, y) = p;
}

void hmm_normalise(Hmm *hmm) {
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        double total = 0.0;
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            total += hmm_getTransition(hmm, from, to);
        }
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            hmm_setTransition(hmm, from, to, hmm_getTransition(hmm, from, to) / total);
        }
    }
    //Normalise the emissions
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        double total = 0.0;
        for (int64_t x = 0; x < SYMBOL_NUMBER_NO_N; x++) {
            for (int64_t y = 0; y < SYMBOL_NUMBER_NO_N; y++) {
                total += hmm_getEmissionsExpectation(hmm, state, x, y);
            }
        }
        for (int64_t x = 0; x < SYMBOL_NUMBER_NO_N; x++) {
            for (int64_t y = 0; y < SYMBOL_NUMBER_NO_N; y++) {
                hmm_setEmissionsExpectation(hmm, state, x, y, hmm_getEmissionsExpectation(hmm, state, x, y) / total);
            }
        }
    }
}

void hmm_randomise(Hmm *hmm) {
    //Transitions
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            hmm_setTransition(hmm, from, to, st_random());
        }
    }
    //Emissions
    for (int64_t state = 0; state < hmm->stateNumber; state++) {
        for (int64_t x = 0; x < SYMBOL_NUMBER_NO_N; x++) {
            for (int64_t y = 0; y < SYMBOL_NUMBER_NO_N; y++) {
                hmm_setEmissionsExpectation(hmm, state, x, y, st_random());
            }
        }
    }

    hmm_normalise(hmm);
}

void hmm_write(Hmm *hmm, FILE *fileHandle) {
    fprintf(fileHandle, "%i\t", hmm->type);
    for (int64_t i = 0; i < hmm->stateNumber * hmm->stateNumber; i++) {
        fprintf(fileHandle, "%f\t", hmm->transitions[i]);
    }
    fprintf(fileHandle, "%f\n", hmm->likelihood);
    for (int64_t i = 0; i < hmm->stateNumber * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N; i++) {
        fprintf(fileHandle, "%f\t", hmm->emissions[i]);
    }
    fprintf(fileHandle, "\n");
}

Hmm *hmm_loadFromFile(const char *fileName) {
    FILE *fH = fopen(fileName, "r");
    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);
    if (stList_length(tokens) < 2) {
        st_errAbort("Got an empty line in the input state machine file %s\n", fileName);
    }
    int type;
    int64_t j = sscanf(stList_get(tokens, 0), "%i", &type);
    if (j != 1) {
        st_errAbort("Failed to parse state number (int) from string: %s\n", string);
    }
    Hmm *hmm = hmm_constructEmpty(0.0, type);
    if (stList_length(tokens) != hmm->stateNumber * hmm->stateNumber + 2) {
        st_errAbort(
                "Got the wrong number of transitions in the input state machine file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), hmm->stateNumber * hmm->stateNumber + 2);
    }

    for (int64_t i = 0; i < hmm->stateNumber * hmm->stateNumber; i++) {
        j = sscanf(stList_get(tokens, i + 1), "%lf", &(hmm->transitions[i]));
        if (j != 1) {
            st_errAbort("Failed to parse transition prob (float) from string: %s\n", string);
        }
    }
    j = sscanf(stList_get(tokens, stList_length(tokens) - 1), "%lf", &(hmm->likelihood));
    if (j != 1) {
        st_errAbort("Failed to parse likelihood (float) from string: %s\n", string);
    }

    //Cleanup transitions line
    free(string);
    stList_destruct(tokens);

    //Now parse the emissions line
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);

    if (stList_length(tokens) != hmm->stateNumber * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N) {
        st_errAbort(
                "Got the wrong number of emissions in the input state machine file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), hmm->stateNumber * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N);
    }

    for (int64_t i = 0; i < hmm->stateNumber * SYMBOL_NUMBER_NO_N * SYMBOL_NUMBER_NO_N; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(hmm->emissions[i]));
        if (j != 1) {
            st_errAbort("Failed to parse emission prob (float) from string: %s\n", string);
        }
    }

    //Final cleanup
    free(string);
    stList_destruct(tokens);
    fclose(fH);

    return hmm;
}

///////////////////////////////////
///////////////////////////////////
//Emissions
///////////////////////////////////
///////////////////////////////////

typedef enum {
    match = 0, shortGapX = 1, shortGapY = 2, longGapX = 3, longGapY = 4
} State;

static void state_check(StateMachine *sM, State s) {
    assert(s >= 0 && s < sM->stateNumber);
}

static void emissions_setMatchProbsToDefaults(double *emissionMatchProbs) {
    /*
     * This is used to set the emissions to reasonable values.
     */
    const double EMISSION_MATCH=-2.1149196655034745; //log(0.12064298095701059);
    const double EMISSION_TRANSVERSION=-4.5691014376830479; //log(0.010367271172731285);
    const double EMISSION_TRANSITION=-3.9833860032220842; //log(0.01862247669752685);
    //Symmetric matrix of transition probabilities.
    const double i[SYMBOL_NUMBER_NO_N*SYMBOL_NUMBER_NO_N] = {
            EMISSION_MATCH, EMISSION_TRANSVERSION, EMISSION_TRANSITION, EMISSION_TRANSVERSION,
            EMISSION_TRANSVERSION, EMISSION_MATCH, EMISSION_TRANSVERSION, EMISSION_TRANSITION,
            EMISSION_TRANSITION, EMISSION_TRANSVERSION, EMISSION_MATCH, EMISSION_TRANSVERSION,
            EMISSION_TRANSVERSION, EMISSION_TRANSITION, EMISSION_TRANSVERSION, EMISSION_MATCH };
    memcpy(emissionMatchProbs, i, sizeof(double)*SYMBOL_NUMBER_NO_N*SYMBOL_NUMBER_NO_N);
}

static void emissions_setGapProbsToDefaults(double *emissionGapProbs) {
    /*
     * This is used to set the emissions to reasonable values.
     */
    const double EMISSION_GAP = -1.6094379124341003; //log(0.2)
    const double i[4] = { EMISSION_GAP, EMISSION_GAP, EMISSION_GAP, EMISSION_GAP };
    memcpy(emissionGapProbs, i, sizeof(double)*SYMBOL_NUMBER_NO_N);
}

static void symbol_check(Symbol c) {
    assert(c >= 0 && c < SYMBOL_NUMBER);
}

static void emissions_loadMatchProbs(double *emissionMatchProbs, Hmm *hmm, int64_t matchState) {
    //Load the matches
    for(int64_t x=0; x<SYMBOL_NUMBER_NO_N; x++) {
        for(int64_t y=0; y<SYMBOL_NUMBER_NO_N; y++) {
            emissionMatchProbs[x * SYMBOL_NUMBER_NO_N + y] = log(hmm_getEmissionsExpectation(hmm, matchState, x, y));
        }
    }
}

static void emissions_loadMatchProbsSymmetrically(double *emissionMatchProbs, Hmm *hmm, int64_t matchState) {
    //Load the matches
    for(int64_t x=0; x<SYMBOL_NUMBER_NO_N; x++) {
        emissionMatchProbs[x * SYMBOL_NUMBER_NO_N + x] = log(hmm_getEmissionsExpectation(hmm, matchState, x, x));
        for(int64_t y=x+1; y<SYMBOL_NUMBER_NO_N; y++) {
            double d = log((hmm_getEmissionsExpectation(hmm, matchState, x, y) + hmm_getEmissionsExpectation(hmm, matchState, y, x))/2.0);
            emissionMatchProbs[x * SYMBOL_NUMBER_NO_N + y] = d;
            emissionMatchProbs[y * SYMBOL_NUMBER_NO_N + x] = d;
        }
    }
}

static void collapseMatrixEmissions(Hmm *hmm, int64_t state, double *gapEmissions, bool collapseToX) {
    for(int64_t x=0; x<SYMBOL_NUMBER_NO_N; x++) {
        for(int64_t y=0; y<SYMBOL_NUMBER_NO_N; y++) {
            gapEmissions[collapseToX ? x : y] += hmm_getEmissionsExpectation(hmm, state, x, y);
        }
    }
}

static void emissions_loadGapProbs(double *emissionGapProbs, Hmm *hmm,
        int64_t *xGapStates, int64_t xGapStateNo,
        int64_t *yGapStates, int64_t yGapStateNo) {
    //Initialise to 0.0
    for(int64_t i=0; i<SYMBOL_NUMBER_NO_N; i++) {
        emissionGapProbs[i] = 0.0;
    }
    //Load the probs taking the average over all the gap states
    for(int64_t i=0; i<xGapStateNo; i++) {
        collapseMatrixEmissions(hmm, xGapStates[i], emissionGapProbs, 1);
    }
    for(int64_t i=0; i<yGapStateNo; i++) {
        collapseMatrixEmissions(hmm, yGapStates[i], emissionGapProbs, 0);
    }
    //Now normalise
    double total = 0.0;
    for(int64_t i=0; i<SYMBOL_NUMBER_NO_N; i++) {
        total += emissionGapProbs[i];
    }
    for(int64_t i=0; i<SYMBOL_NUMBER_NO_N; i++) {
        emissionGapProbs[i] = log(emissionGapProbs[i]/total);
    }
}

static inline double emission_getGapProb(const double *emissionGapProbs, Symbol i) {
    symbol_check(i);
    if(i == n) {
        return -1.386294361; //log(0.25)
    }
    return emissionGapProbs[i];
}

static inline double emission_getMatchProb(const double *emissionMatchProbs, Symbol x, Symbol y) {
    symbol_check(x);
    symbol_check(y);
    if(x == n || y == n) {
        return -2.772588722; //log(0.25**2)
    }
    return emissionMatchProbs[x * SYMBOL_NUMBER_NO_N + y];
}

///////////////////////////////////
///////////////////////////////////
//Five state state-machine
///////////////////////////////////
///////////////////////////////////

//Transitions
typedef struct _StateMachine5 StateMachine5;

struct _StateMachine5 {
    StateMachine model;
    double TRANSITION_MATCH_CONTINUE; //0.9703833696510062f
    double TRANSITION_MATCH_FROM_SHORT_GAP_X; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_LONG_GAP_X; //1.0 - gapExtend = 0.00343657420938
    double TRANSITION_GAP_SHORT_OPEN_X; //0.0129868352330243
    double TRANSITION_GAP_SHORT_EXTEND_X; //0.7126062401851738f;
    double TRANSITION_GAP_SHORT_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_GAP_LONG_OPEN_X; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    double TRANSITION_GAP_LONG_EXTEND_X; //0.99656342579062f;
    double TRANSITION_GAP_LONG_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_MATCH_FROM_SHORT_GAP_Y; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_LONG_GAP_Y; //1.0 - gapExtend = 0.00343657420938
    double TRANSITION_GAP_SHORT_OPEN_Y; //0.0129868352330243
    double TRANSITION_GAP_SHORT_EXTEND_Y; //0.7126062401851738f;
    double TRANSITION_GAP_SHORT_SWITCH_TO_Y; //0.0073673675173412815f;
    double TRANSITION_GAP_LONG_OPEN_Y; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    double TRANSITION_GAP_LONG_EXTEND_Y; //0.99656342579062f;
    double TRANSITION_GAP_LONG_SWITCH_TO_Y; //0.0073673675173412815f;
    double EMISSION_MATCH_PROBS[SYMBOL_NUMBER_NO_N*SYMBOL_NUMBER_NO_N]; //Match emission probs
    double EMISSION_GAP_X_PROBS[SYMBOL_NUMBER_NO_N]; //Gap emission probs
    double EMISSION_GAP_Y_PROBS[SYMBOL_NUMBER_NO_N]; //Gap emission probs
};

static double stateMachine5_startStateProb(StateMachine *sM, int64_t state) {
    //Match state is like going to a match.
    state_check(sM, state);
    return state == match ? 0 : LOG_ZERO;
}

static double stateMachine5_raggedStartStateProb(StateMachine *sM, int64_t state) {
    state_check(sM, state);
    return (state == longGapX || state == longGapY) ? 0 : LOG_ZERO;
}

static double stateMachine5_endStateProb(StateMachine *sM, int64_t state) {
    //End state is like to going to a match
    StateMachine5 *sM5 = (StateMachine5 *) sM;
    state_check(sM, state);
    switch (state) {
    case match:
        return sM5->TRANSITION_MATCH_CONTINUE;
    case shortGapX:
        return sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X;
    case shortGapY:
        return sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y;
    case longGapX:
        return sM5->TRANSITION_MATCH_FROM_LONG_GAP_X;
    case longGapY:
        return sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y;
    }
    return 0.0;
}

static double stateMachine5_raggedEndStateProb(StateMachine *sM, int64_t state) {
    //End state is like to going to a match
    StateMachine5 *sM5 = (StateMachine5 *) sM;
    state_check(sM, state);
    switch (state) {
    case match:
        return sM5->TRANSITION_GAP_LONG_OPEN_X;
    case shortGapX:
        return sM5->TRANSITION_GAP_LONG_OPEN_X;
    case shortGapY:
        return sM5->TRANSITION_GAP_LONG_OPEN_Y;
    case longGapX:
        return sM5->TRANSITION_GAP_LONG_EXTEND_X;
    case longGapY:
        return sM5->TRANSITION_GAP_LONG_EXTEND_Y;
    }
    return 0.0;
}

static void stateMachine5_cellCalculate(StateMachine *sM, double *current, double *lower, double *middle, double *upper,
        Symbol cX, Symbol cY, void (*doTransition)(double *, double *, int64_t, int64_t, double, double, void *),
        void *extraArgs) {
    StateMachine5 *sM5 = (StateMachine5 *) sM;
    if (lower != NULL) {
        double eP = emission_getGapProb(sM5->EMISSION_GAP_X_PROBS, cX);
        doTransition(lower, current, match, shortGapX, eP, sM5->TRANSITION_GAP_SHORT_OPEN_X, extraArgs);
        doTransition(lower, current, shortGapX, shortGapX, eP, sM5->TRANSITION_GAP_SHORT_EXTEND_X, extraArgs);
        //doTransition(lower, current, shortGapY, shortGapX, eP, sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X, extraArgs);
        doTransition(lower, current, match, longGapX, eP, sM5->TRANSITION_GAP_LONG_OPEN_X, extraArgs);
        doTransition(lower, current, longGapX, longGapX, eP, sM5->TRANSITION_GAP_LONG_EXTEND_X, extraArgs);
        //doTransition(lower, current, longGapY, longGapX, eP, sM5->TRANSITION_GAP_LONG_SWITCH_TO_X, extraArgs);
    }
    if (middle != NULL) {
        double eP = emission_getMatchProb(sM5->EMISSION_MATCH_PROBS, cX, cY); //symbol_matchProb(cX, cY);
        doTransition(middle, current, match, match, eP, sM5->TRANSITION_MATCH_CONTINUE, extraArgs);
        doTransition(middle, current, shortGapX, match, eP, sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X, extraArgs);
        doTransition(middle, current, shortGapY, match, eP, sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y, extraArgs);
        doTransition(middle, current, longGapX, match, eP, sM5->TRANSITION_MATCH_FROM_LONG_GAP_X, extraArgs);
        doTransition(middle, current, longGapY, match, eP, sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y, extraArgs);
    }
    if (upper != NULL) {
        double eP = emission_getGapProb(sM5->EMISSION_GAP_Y_PROBS, cY);
        doTransition(upper, current, match, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_OPEN_Y, extraArgs);
        doTransition(upper, current, shortGapY, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_EXTEND_Y, extraArgs);
        //doTransition(upper, current, shortGapX, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y, extraArgs);
        doTransition(upper, current, match, longGapY, eP, sM5->TRANSITION_GAP_LONG_OPEN_Y, extraArgs);
        doTransition(upper, current, longGapY, longGapY, eP, sM5->TRANSITION_GAP_LONG_EXTEND_Y, extraArgs);
        //doTransition(upper, current, longGapX, longGapY, eP, sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y, extraArgs);
    }
}

StateMachine *stateMachine5_construct(StateMachineType type) {
    StateMachine5 *sM5 = st_malloc(sizeof(StateMachine5));
    sM5->TRANSITION_MATCH_CONTINUE = -0.030064059121770816; //0.9703833696510062f
    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X = -1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_X = -5.673280173170473; //1.0 - gapExtend = 0.00343657420938
    sM5->TRANSITION_GAP_SHORT_OPEN_X = -4.34381910900448; //0.0129868352330243
    sM5->TRANSITION_GAP_SHORT_EXTEND_X = -0.3388262689231553; //0.7126062401851738f;
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X = -4.910694825551255; //0.0073673675173412815f;
    sM5->TRANSITION_GAP_LONG_OPEN_X = -6.30810595366929; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    sM5->TRANSITION_GAP_LONG_EXTEND_X = -0.003442492794189331; //0.99656342579062f;
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_X = -6.30810595366929; //0.99656342579062f;

    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y = sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X;
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y = sM5->TRANSITION_MATCH_FROM_LONG_GAP_X;
    sM5->TRANSITION_GAP_SHORT_OPEN_Y = sM5->TRANSITION_GAP_SHORT_OPEN_X;
    sM5->TRANSITION_GAP_SHORT_EXTEND_Y = sM5->TRANSITION_GAP_SHORT_EXTEND_X;
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y = sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X;
    sM5->TRANSITION_GAP_LONG_OPEN_Y = sM5->TRANSITION_GAP_LONG_OPEN_X;
    sM5->TRANSITION_GAP_LONG_EXTEND_Y = sM5->TRANSITION_GAP_LONG_EXTEND_X;
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y = sM5->TRANSITION_GAP_LONG_SWITCH_TO_X;

    emissions_setMatchProbsToDefaults(sM5->EMISSION_MATCH_PROBS);
    emissions_setGapProbsToDefaults(sM5->EMISSION_GAP_X_PROBS);
    emissions_setGapProbsToDefaults(sM5->EMISSION_GAP_Y_PROBS);
    if(type != fiveState && type != fiveStateAsymmetric) {
        st_errAbort("Wrong type for five state %i", type);
    }
    sM5->model.type = type;
    sM5->model.stateNumber = 5;
    sM5->model.matchState = match;
    sM5->model.startStateProb = stateMachine5_startStateProb;
    sM5->model.endStateProb = stateMachine5_endStateProb;
    sM5->model.raggedStartStateProb = stateMachine5_raggedStartStateProb;
    sM5->model.raggedEndStateProb = stateMachine5_raggedEndStateProb;
    sM5->model.cellCalculate = stateMachine5_cellCalculate;

    return (StateMachine *) sM5;
}

static void switchDoubles(double *a, double *b) {
    double c = *a;
    *a = *b;
    *b = c;
}

static void stateMachine5_loadAsymmetric(StateMachine5 *sM5, Hmm *hmm) {
    if (hmm->type != fiveStateAsymmetric) {
        st_errAbort("Wrong hmm type");
    }
    sM5->TRANSITION_MATCH_CONTINUE = log(hmm_getTransition(hmm, match, match)); //0.9703833696510062f

    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X = log(hmm_getTransition(hmm, shortGapX, match));
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_X = log(hmm_getTransition(hmm, longGapX, match));
    sM5->TRANSITION_GAP_SHORT_OPEN_X = log(hmm_getTransition(hmm, match, shortGapX));
    sM5->TRANSITION_GAP_SHORT_EXTEND_X = log(hmm_getTransition(hmm, shortGapX, shortGapX));
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X = log(hmm_getTransition(hmm, shortGapY, shortGapX));
    sM5->TRANSITION_GAP_LONG_OPEN_X = log(hmm_getTransition(hmm, match, longGapX));
    sM5->TRANSITION_GAP_LONG_EXTEND_X = log(hmm_getTransition(hmm, longGapX, longGapX));
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_X = log(hmm_getTransition(hmm, longGapY, longGapX));

    if(sM5->TRANSITION_GAP_SHORT_EXTEND_X > sM5->TRANSITION_GAP_LONG_EXTEND_X) {
        //Switch the long and short gap parameters if one the "long states" have a smaller extend probability than the "short states", as can randomly happen during EM training.
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_EXTEND_X), &(sM5->TRANSITION_GAP_LONG_EXTEND_X));
        switchDoubles(&(sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X), &(sM5->TRANSITION_MATCH_FROM_LONG_GAP_X));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_OPEN_X), &(sM5->TRANSITION_GAP_LONG_OPEN_X));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X), &(sM5->TRANSITION_GAP_LONG_SWITCH_TO_X));
    }

    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y = log(hmm_getTransition(hmm, shortGapY, match));
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y = log(hmm_getTransition(hmm, longGapY, match));
    sM5->TRANSITION_GAP_SHORT_OPEN_Y = log(hmm_getTransition(hmm, match, shortGapY));
    sM5->TRANSITION_GAP_SHORT_EXTEND_Y = log(hmm_getTransition(hmm, shortGapY, shortGapY));
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y = log(hmm_getTransition(hmm, shortGapX, shortGapY));
    sM5->TRANSITION_GAP_LONG_OPEN_Y = log(hmm_getTransition(hmm, match, longGapY));
    sM5->TRANSITION_GAP_LONG_EXTEND_Y = log(hmm_getTransition(hmm, longGapY, longGapY));
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y = log(hmm_getTransition(hmm, longGapX, longGapY));

    if(sM5->TRANSITION_GAP_SHORT_EXTEND_Y > sM5->TRANSITION_GAP_LONG_EXTEND_Y) {
        //Switch the long and short gap parameters if one the "long states" have a smaller extend probability than the "short states", as can randomly happen during EM training.
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_EXTEND_Y), &(sM5->TRANSITION_GAP_LONG_EXTEND_Y));
        switchDoubles(&(sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y), &(sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_OPEN_Y), &(sM5->TRANSITION_GAP_LONG_OPEN_Y));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y), &(sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y));
    }

    emissions_loadMatchProbs(sM5->EMISSION_MATCH_PROBS, hmm, match);
    int64_t xGapStates[2] = { shortGapX, longGapX };
    int64_t yGapStates[2] = { shortGapY, longGapY };
    emissions_loadGapProbs(sM5->EMISSION_GAP_X_PROBS, hmm, xGapStates, 2, NULL, 0);
    emissions_loadGapProbs(sM5->EMISSION_GAP_Y_PROBS, hmm, NULL, 0, yGapStates, 2);
}

static void stateMachine5_loadSymmetric(StateMachine5 *sM5, Hmm *hmm) {
    if (hmm->type != fiveState) {
        st_errAbort("Wrong hmm type");
    }
    sM5->TRANSITION_MATCH_CONTINUE = log(hmm_getTransition(hmm, match, match)); //0.9703833696510062f
    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X = log(
            (hmm_getTransition(hmm, shortGapX, match) + hmm_getTransition(hmm, shortGapY, match)) / 2); //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_X = log(
            (hmm_getTransition(hmm, longGapX, match) + hmm_getTransition(hmm, longGapY, match)) / 2); //1.0 - gapExtend = 0.00343657420938
    sM5->TRANSITION_GAP_SHORT_OPEN_X = log(
            (hmm_getTransition(hmm, match, shortGapX) + hmm_getTransition(hmm, match, shortGapY)) / 2); //0.0129868352330243
    sM5->TRANSITION_GAP_SHORT_EXTEND_X = log(
            (hmm_getTransition(hmm, shortGapX, shortGapX) + hmm_getTransition(hmm, shortGapY, shortGapY)) / 2); //0.7126062401851738f;
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X = log(
            (hmm_getTransition(hmm, shortGapX, shortGapY) + hmm_getTransition(hmm, shortGapY, shortGapX)) / 2); //0.0073673675173412815f;
    sM5->TRANSITION_GAP_LONG_OPEN_X = log(
            (hmm_getTransition(hmm, match, longGapX) + hmm_getTransition(hmm, match, longGapY)) / 2); //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    sM5->TRANSITION_GAP_LONG_EXTEND_X = log(
            (hmm_getTransition(hmm, longGapX, longGapX) + hmm_getTransition(hmm, longGapY, longGapY)) / 2);
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_X = log(
                (hmm_getTransition(hmm, longGapX, longGapY) + hmm_getTransition(hmm, longGapY, longGapX)) / 2); //0.0073673675173412815f;

    if(sM5->TRANSITION_GAP_SHORT_EXTEND_X > sM5->TRANSITION_GAP_LONG_EXTEND_X) {
        //Switch the long and short gap parameters if one the "long states" have a smaller extend probability than the "short states", as can randomly happen during EM training.
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_EXTEND_X), &(sM5->TRANSITION_GAP_LONG_EXTEND_X));
        switchDoubles(&(sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X), &(sM5->TRANSITION_MATCH_FROM_LONG_GAP_X));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_OPEN_X), &(sM5->TRANSITION_GAP_LONG_OPEN_X));
        switchDoubles(&(sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X), &(sM5->TRANSITION_GAP_LONG_SWITCH_TO_X));
    }

    sM5->TRANSITION_MATCH_FROM_SHORT_GAP_Y = sM5->TRANSITION_MATCH_FROM_SHORT_GAP_X;
    sM5->TRANSITION_MATCH_FROM_LONG_GAP_Y = sM5->TRANSITION_MATCH_FROM_LONG_GAP_X;
    sM5->TRANSITION_GAP_SHORT_OPEN_Y = sM5->TRANSITION_GAP_SHORT_OPEN_X;
    sM5->TRANSITION_GAP_SHORT_EXTEND_Y = sM5->TRANSITION_GAP_SHORT_EXTEND_X;
    sM5->TRANSITION_GAP_SHORT_SWITCH_TO_Y = sM5->TRANSITION_GAP_SHORT_SWITCH_TO_X;
    sM5->TRANSITION_GAP_LONG_OPEN_Y = sM5->TRANSITION_GAP_LONG_OPEN_X;
    sM5->TRANSITION_GAP_LONG_EXTEND_Y = sM5->TRANSITION_GAP_LONG_EXTEND_X;
    sM5->TRANSITION_GAP_LONG_SWITCH_TO_Y = sM5->TRANSITION_GAP_LONG_SWITCH_TO_X;

    emissions_loadMatchProbsSymmetrically(sM5->EMISSION_MATCH_PROBS, hmm, match);
    int64_t xGapStates[2] = { shortGapX, longGapX };
    int64_t yGapStates[2] = { shortGapY, longGapY };
    emissions_loadGapProbs(sM5->EMISSION_GAP_X_PROBS, hmm, xGapStates, 2, yGapStates, 2);
    emissions_loadGapProbs(sM5->EMISSION_GAP_Y_PROBS, hmm, xGapStates, 2, yGapStates, 2);
}

///////////////////////////////////
///////////////////////////////////
//Three state state-machine
///////////////////////////////////
///////////////////////////////////

//Transitions
typedef struct _StateMachine3 StateMachine3;

struct _StateMachine3 {
    //3 state state machine, allowing for symmetry in x and y.
    StateMachine model;
    double TRANSITION_MATCH_CONTINUE; //0.9703833696510062f
    double TRANSITION_MATCH_FROM_GAP_X; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_GAP_Y; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_GAP_OPEN_X; //0.0129868352330243
    double TRANSITION_GAP_OPEN_Y; //0.0129868352330243
    double TRANSITION_GAP_EXTEND_X; //0.7126062401851738f;
    double TRANSITION_GAP_EXTEND_Y; //0.7126062401851738f;
    double TRANSITION_GAP_SWITCH_TO_X; //0.0073673675173412815f;
    double TRANSITION_GAP_SWITCH_TO_Y; //0.0073673675173412815f;
    double EMISSION_MATCH_PROBS[SYMBOL_NUMBER_NO_N*SYMBOL_NUMBER_NO_N]; //Match emission probs
    double EMISSION_GAP_X_PROBS[SYMBOL_NUMBER_NO_N]; //Gap X emission probs
    double EMISSION_GAP_Y_PROBS[SYMBOL_NUMBER_NO_N]; //Gap Y emission probs
};

static double stateMachine3_startStateProb(StateMachine *sM, int64_t state) {
    //Match state is like going to a match.
    state_check(sM, state);
    return state == match ? 0 : LOG_ZERO;
}

static double stateMachine3_raggedStartStateProb(StateMachine *sM, int64_t state) {
    state_check(sM, state);
    return (state == shortGapX || state == shortGapY) ? 0 : LOG_ZERO;
}

static double stateMachine3_endStateProb(StateMachine *sM, int64_t state) {
    //End state is like to going to a match
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    state_check(sM, state);
    switch (state) {
    case match:
        return sM3->TRANSITION_MATCH_CONTINUE;
    case shortGapX:
        return sM3->TRANSITION_MATCH_FROM_GAP_X;
    case shortGapY:
        return sM3->TRANSITION_MATCH_FROM_GAP_Y;
    }
    return 0.0;
}

static double stateMachine3_raggedEndStateProb(StateMachine *sM, int64_t state) {
    //End state is like to going to a match
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    state_check(sM, state);
    switch (state) {
    case match:
        return (sM3->TRANSITION_GAP_OPEN_X + sM3->TRANSITION_GAP_OPEN_Y) / 2.0;
    case shortGapX:
        return sM3->TRANSITION_GAP_EXTEND_X;
    case shortGapY:
        return sM3->TRANSITION_GAP_EXTEND_Y;
    }
    return 0.0;
}

static void stateMachine3_cellCalculate(StateMachine *sM, double *current, double *lower, double *middle, double *upper,
        Symbol cX, Symbol cY, void (*doTransition)(double *, double *, int64_t, int64_t, double, double, void *),
        void *extraArgs) {
    symbol_check(cX);
    symbol_check(cY);
    StateMachine3 *sM3 = (StateMachine3 *) sM;
    if (lower != NULL) {
        double eP = emission_getGapProb(sM3->EMISSION_GAP_X_PROBS, cX);
        doTransition(lower, current, match, shortGapX, eP, sM3->TRANSITION_GAP_OPEN_X, extraArgs);
        doTransition(lower, current, shortGapX, shortGapX, eP, sM3->TRANSITION_GAP_EXTEND_X, extraArgs);
        doTransition(lower, current, shortGapY, shortGapX, eP, sM3->TRANSITION_GAP_SWITCH_TO_X, extraArgs);
    }
    if (middle != NULL) {
        double eP = emission_getMatchProb(sM3->EMISSION_MATCH_PROBS, cX, cY);  //symbol_matchProb(cX, cY);
        doTransition(middle, current, match, match, eP, sM3->TRANSITION_MATCH_CONTINUE, extraArgs);
        doTransition(middle, current, shortGapX, match, eP, sM3->TRANSITION_MATCH_FROM_GAP_X, extraArgs);
        doTransition(middle, current, shortGapY, match, eP, sM3->TRANSITION_MATCH_FROM_GAP_Y, extraArgs);

    }
    if (upper != NULL) {
        double eP = emission_getGapProb(sM3->EMISSION_GAP_Y_PROBS, cY);
        doTransition(upper, current, match, shortGapY, eP, sM3->TRANSITION_GAP_OPEN_Y, extraArgs);
        doTransition(upper, current, shortGapY, shortGapY, eP, sM3->TRANSITION_GAP_EXTEND_Y, extraArgs);
        doTransition(upper, current, shortGapX, shortGapY, eP, sM3->TRANSITION_GAP_SWITCH_TO_Y, extraArgs);
    }
}

StateMachine *stateMachine3_construct(StateMachineType type) {
    StateMachine3 *sM3 = st_malloc(sizeof(StateMachine3));
    sM3->TRANSITION_MATCH_CONTINUE = -0.030064059121770816; //0.9703833696510062f
    sM3->TRANSITION_MATCH_FROM_GAP_X = -1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM3->TRANSITION_MATCH_FROM_GAP_Y = -1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM3->TRANSITION_GAP_OPEN_X = -4.21256642; //0.0129868352330243
    sM3->TRANSITION_GAP_OPEN_Y = -4.21256642; //0.0129868352330243
    sM3->TRANSITION_GAP_EXTEND_X = -0.3388262689231553; //0.7126062401851738f;
    sM3->TRANSITION_GAP_EXTEND_Y = -0.3388262689231553; //0.7126062401851738f;
    sM3->TRANSITION_GAP_SWITCH_TO_X = -4.910694825551255; //0.0073673675173412815f;
    sM3->TRANSITION_GAP_SWITCH_TO_Y = -4.910694825551255; //0.0073673675173412815f;
    emissions_setMatchProbsToDefaults(sM3->EMISSION_MATCH_PROBS);
    emissions_setGapProbsToDefaults(sM3->EMISSION_GAP_X_PROBS);
    emissions_setGapProbsToDefaults(sM3->EMISSION_GAP_Y_PROBS);
    if (type != threeState && type != threeStateAsymmetric) {
        st_errAbort("Tried to create a three state state-machine with the wrong type");
    }
    sM3->model.type = type;
    sM3->model.stateNumber = 3;
    sM3->model.matchState = match;
    sM3->model.startStateProb = stateMachine3_startStateProb;
    sM3->model.endStateProb = stateMachine3_endStateProb;
    sM3->model.raggedStartStateProb = stateMachine3_raggedStartStateProb;
    sM3->model.raggedEndStateProb = stateMachine3_raggedEndStateProb;
    sM3->model.cellCalculate = stateMachine3_cellCalculate;

    return (StateMachine *) sM3;
}

static void stateMachine3_loadAsymmetric(StateMachine3 *sM3, Hmm *hmm) {
    if (hmm->type != threeStateAsymmetric) {
        st_errAbort("Wrong hmm type");
    }
    sM3->TRANSITION_MATCH_CONTINUE = log(hmm_getTransition(hmm, match, match));
    sM3->TRANSITION_MATCH_FROM_GAP_X = log(hmm_getTransition(hmm, shortGapX, match));
    sM3->TRANSITION_MATCH_FROM_GAP_Y = log(hmm_getTransition(hmm, shortGapY, match));
    sM3->TRANSITION_GAP_OPEN_X = log(hmm_getTransition(hmm, match, shortGapX));
    sM3->TRANSITION_GAP_OPEN_Y = log(hmm_getTransition(hmm, match, shortGapY));
    sM3->TRANSITION_GAP_EXTEND_X = log(hmm_getTransition(hmm, shortGapX, shortGapX));
    sM3->TRANSITION_GAP_EXTEND_Y = log(hmm_getTransition(hmm, shortGapY, shortGapY));
    sM3->TRANSITION_GAP_SWITCH_TO_X = log(hmm_getTransition(hmm, shortGapY, shortGapX));
    sM3->TRANSITION_GAP_SWITCH_TO_Y = log(hmm_getTransition(hmm, shortGapX, shortGapY));
    emissions_loadMatchProbs(sM3->EMISSION_MATCH_PROBS, hmm, match);
    int64_t xGapStates[1] = { shortGapX };
    int64_t yGapStates[1] = { shortGapY };
    emissions_loadGapProbs(sM3->EMISSION_GAP_X_PROBS, hmm, xGapStates, 1, NULL, 0);
    emissions_loadGapProbs(sM3->EMISSION_GAP_Y_PROBS, hmm, NULL, 0, yGapStates, 1);
}

static void stateMachine3_loadSymmetric(StateMachine3 *sM3, Hmm *hmm) {
    if (hmm->type != threeState) {
        st_errAbort("Wrong hmm type");
    }
    sM3->TRANSITION_MATCH_CONTINUE = log(hmm_getTransition(hmm, match, match));
    sM3->TRANSITION_MATCH_FROM_GAP_X = log(
            (hmm_getTransition(hmm, shortGapX, match) + hmm_getTransition(hmm, shortGapY, match)) / 2.0);
    sM3->TRANSITION_MATCH_FROM_GAP_Y = sM3->TRANSITION_MATCH_FROM_GAP_X;
    sM3->TRANSITION_GAP_OPEN_X = log(
            (hmm_getTransition(hmm, match, shortGapX) + hmm_getTransition(hmm, match, shortGapY)) / 2.0);
    sM3->TRANSITION_GAP_OPEN_Y = sM3->TRANSITION_GAP_OPEN_X;
    sM3->TRANSITION_GAP_EXTEND_X = log(
            (hmm_getTransition(hmm, shortGapX, shortGapX) + hmm_getTransition(hmm, shortGapY, shortGapY)) / 2.0);
    sM3->TRANSITION_GAP_EXTEND_Y = sM3->TRANSITION_GAP_EXTEND_X;
    sM3->TRANSITION_GAP_SWITCH_TO_X = log(
            (hmm_getTransition(hmm, shortGapY, shortGapX) + hmm_getTransition(hmm, shortGapX, shortGapY)) / 2.0);
    sM3->TRANSITION_GAP_SWITCH_TO_Y = sM3->TRANSITION_GAP_SWITCH_TO_X;
    emissions_loadMatchProbsSymmetrically(sM3->EMISSION_MATCH_PROBS, hmm, match);
    int64_t xGapStates[2] = { shortGapX };
    int64_t yGapStates[2] = { shortGapY };
    emissions_loadGapProbs(sM3->EMISSION_GAP_X_PROBS, hmm, xGapStates, 1, yGapStates, 1);
    emissions_loadGapProbs(sM3->EMISSION_GAP_Y_PROBS, hmm, xGapStates, 1, yGapStates, 1);
}

///////////////////////////////////
///////////////////////////////////
//Public functions
///////////////////////////////////
///////////////////////////////////

StateMachine *hmm_getStateMachine(Hmm *hmm) {
    if (hmm->type == fiveState) {
        StateMachine5 *sM5 = (StateMachine5 *) stateMachine5_construct(fiveState);
        stateMachine5_loadSymmetric(sM5, hmm);
        return (StateMachine *) sM5;
    }
    if (hmm->type == fiveStateAsymmetric) {
        StateMachine5 *sM5 = (StateMachine5 *) stateMachine5_construct(fiveStateAsymmetric);
        stateMachine5_loadAsymmetric(sM5, hmm);
        return (StateMachine *) sM5;
    }
    if (hmm->type == threeStateAsymmetric) {
        StateMachine3 *sM3 = (StateMachine3 *) stateMachine3_construct(hmm->type);
        stateMachine3_loadAsymmetric(sM3, hmm);
        return (StateMachine *) sM3;
    }
    if (hmm->type == threeState) {
        StateMachine3 *sM3 = (StateMachine3 *) stateMachine3_construct(hmm->type);
        stateMachine3_loadSymmetric(sM3, hmm);
        return (StateMachine *) sM3;
    }
    return NULL;
}

void stateMachine_destruct(StateMachine *stateMachine) {
    free(stateMachine);
}
