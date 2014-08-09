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
    switch(type) {
        case fiveState:
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
    hmm->likelihood = 0.0;
    return hmm;
}

void hmm_destruct(Hmm *hmm) {
    free(hmm->transitions);
    free(hmm);
}

void hmm_addToTransitionExpectation(Hmm *hmm, int64_t from, int64_t to, double p) {
    hmm->transitions[from * hmm->stateNumber + to] += p;
}

void hmm_setTransition(Hmm *hmm, int64_t from, int64_t to, double p) {
    hmm->transitions[from * hmm->stateNumber + to] = p;
}

double hmm_getTransition(Hmm *hmm, int64_t from, int64_t to) {
    return hmm->transitions[from * hmm->stateNumber + to];
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
}

void hmm_randomise(Hmm *hmm) {
    for (int64_t from = 0; from < hmm->stateNumber; from++) {
        for (int64_t to = 0; to < hmm->stateNumber; to++) {
            hmm_setTransition(hmm, from, to, st_random());
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
        st_errAbort("Got the wrong number of transitions in the input state machine file %s, got %" PRIi64 " instead of %" PRIi64 "\n",
                fileName, stList_length(tokens), hmm->stateNumber * hmm->stateNumber + 2);
    }

    for (int64_t i = 0; i < hmm->stateNumber * hmm->stateNumber; i++) {
        j = sscanf(stList_get(tokens, i+1), "%lf", &(hmm->transitions[i]));
        if (j != 1) {
            st_errAbort("Failed to parse transition prob (float) from string: %s\n", string);
        }
    }
    j = sscanf(stList_get(tokens, stList_length(tokens) - 1), "%lf", &(hmm->likelihood));
    if (j != 1) {
        st_errAbort("Failed to parse likelihood (float) from string: %s\n", string);
    }
    fclose(fH);
    return hmm;
}

///////////////////////////////////
///////////////////////////////////
//Emissions
///////////////////////////////////
///////////////////////////////////


//Emissions
#define EMISSION_MATCH -2.1149196655034745 //log(0.12064298095701059);
#define EMISSION_TRANSVERSION -4.5691014376830479 //log(0.010367271172731285);
#define EMISSION_TRANSITION -3.9833860032220842 //log(0.01862247669752685);
#define EMISSION_MATCH_N -3.2188758248682006 //log(0.04);

static void symbol_check(Symbol c) {
    assert(c >= 0 && c < 5);
}

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
//Five state state-machine
///////////////////////////////////
///////////////////////////////////

typedef enum {
    match=0,
    shortGapX=1,
    shortGapY=2,
    longGapX=3,
    longGapY=4
} State;

static void state_check(StateMachine *sM, State s) {
    assert(s >= 0 && s < sM->stateNumber);
}

//Transitions
typedef struct _StateMachine5 StateMachine5;

struct _StateMachine5 {
    StateMachine model;
    double TRANSITION_MATCH_CONTINUE; //0.9703833696510062f
    double TRANSITION_MATCH_FROM_SHORT_GAP; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    double TRANSITION_MATCH_FROM_LONG_GAP; //1.0 - gapExtend = 0.00343657420938
    double TRANSITION_GAP_SHORT_OPEN; //0.0129868352330243
    double TRANSITION_GAP_SHORT_EXTEND; //0.7126062401851738f;
    double TRANSITION_GAP_SHORT_SWITCH; //0.0073673675173412815f;
    double TRANSITION_GAP_LONG_OPEN; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    double TRANSITION_GAP_LONG_EXTEND; //0.99656342579062f;
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
    StateMachine5 *sM5 = (StateMachine5 *)sM;
    state_check(sM, state);
    switch(state) {
        case match :
            return sM5->TRANSITION_MATCH_CONTINUE;
        case shortGapX:
        case shortGapY:
            return sM5->TRANSITION_MATCH_FROM_SHORT_GAP;
        case longGapX:
        case longGapY:
            return sM5->TRANSITION_MATCH_FROM_LONG_GAP;
    }
    return 0.0;
}

static double stateMachine5_raggedEndStateProb(StateMachine *sM, int64_t state) {
    //End state is like to going to a match
    StateMachine5 *sM5 = (StateMachine5 *)sM;
    state_check(sM, state);
    switch(state) {
        case match :
            return sM5->TRANSITION_GAP_LONG_OPEN;
        case shortGapX:
        case shortGapY:
            return sM5->TRANSITION_GAP_LONG_OPEN;
        case longGapX:
        case longGapY:
            return sM5->TRANSITION_GAP_LONG_EXTEND;
    }
    return 0.0;
}

static void stateMachine5_cellCalculate(StateMachine *sM, double *current,
        double *lower, double *middle, double *upper, Symbol cX, Symbol cY,
        void(*doTransition)(double *, double *, int64_t, int64_t, double, double, void *), void *extraArgs) {
    StateMachine5 *sM5 = (StateMachine5 *)sM;
    if (lower != NULL) {
        double eP = symbol_gapProb(cX);
        doTransition(lower, current, match, shortGapX, eP, sM5->TRANSITION_GAP_SHORT_OPEN, extraArgs);
        doTransition(lower, current, shortGapX, shortGapX, eP, sM5->TRANSITION_GAP_SHORT_EXTEND, extraArgs);
        doTransition(lower, current, shortGapY, shortGapX, eP, sM5->TRANSITION_GAP_SHORT_SWITCH, extraArgs);
        doTransition(lower, current, match, longGapX, eP, sM5->TRANSITION_GAP_LONG_OPEN, extraArgs);
        doTransition(lower, current, longGapX, longGapX, eP, sM5->TRANSITION_GAP_LONG_EXTEND, extraArgs);
    }
    if (middle != NULL) {
        double eP = symbol_matchProb(cX, cY);
        doTransition(middle, current, match, match, eP, sM5->TRANSITION_MATCH_CONTINUE, extraArgs);
        doTransition(middle, current, shortGapX, match, eP, sM5->TRANSITION_MATCH_FROM_SHORT_GAP, extraArgs);
        doTransition(middle, current, shortGapY, match, eP, sM5->TRANSITION_MATCH_FROM_SHORT_GAP, extraArgs);
        doTransition(middle, current, longGapX, match, eP, sM5->TRANSITION_MATCH_FROM_LONG_GAP, extraArgs);
        doTransition(middle, current, longGapY, match, eP, sM5->TRANSITION_MATCH_FROM_LONG_GAP, extraArgs);
    }
    if (upper != NULL) {
        double eP = symbol_gapProb(cY);
        doTransition(upper, current, match, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_OPEN, extraArgs);
        doTransition(upper, current, shortGapY, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_EXTEND, extraArgs);
        doTransition(upper, current, shortGapX, shortGapY, eP, sM5->TRANSITION_GAP_SHORT_SWITCH, extraArgs);
        doTransition(upper, current, match, longGapY, eP, sM5->TRANSITION_GAP_LONG_OPEN, extraArgs);
        doTransition(upper, current, longGapY, longGapY, eP, sM5->TRANSITION_GAP_LONG_EXTEND, extraArgs);
    }
}

StateMachine *stateMachine5_construct() {
    StateMachine5 *sM5 = st_malloc(sizeof(StateMachine5));
    sM5->TRANSITION_MATCH_CONTINUE=-0.030064059121770816; //0.9703833696510062f
    sM5->TRANSITION_MATCH_FROM_SHORT_GAP=-1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM5->TRANSITION_MATCH_FROM_LONG_GAP=-5.673280173170473; //1.0 - gapExtend = 0.00343657420938
    sM5->TRANSITION_GAP_SHORT_OPEN=-4.34381910900448; //0.0129868352330243
    sM5->TRANSITION_GAP_SHORT_EXTEND=-0.3388262689231553; //0.7126062401851738f;
    sM5->TRANSITION_GAP_SHORT_SWITCH=-4.910694825551255; //0.0073673675173412815f;
    sM5->TRANSITION_GAP_LONG_OPEN=-6.30810595366929; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    sM5->TRANSITION_GAP_LONG_EXTEND=-0.003442492794189331; //0.99656342579062f;
    sM5->model.type = fiveState;
    sM5->model.stateNumber = 5;
    sM5->model.matchState = match;
    sM5->model.startStateProb = stateMachine5_startStateProb;
    sM5->model.endStateProb = stateMachine5_endStateProb;
    sM5->model.raggedStartStateProb = stateMachine5_raggedStartStateProb;
    sM5->model.raggedEndStateProb = stateMachine5_raggedEndStateProb;
    sM5->model.cellCalculate = stateMachine5_cellCalculate;

    return (StateMachine *)sM5;
}

static void stateMachine5_load(StateMachine5 *sM5, Hmm *hmm) {
    if(hmm->type != fiveState) {
        st_errAbort("Wrong hmm type");
    }
    sM5->TRANSITION_MATCH_CONTINUE=log(hmm_getTransition(hmm, match, match)); //0.9703833696510062f
    sM5->TRANSITION_MATCH_FROM_SHORT_GAP=log((hmm_getTransition(hmm, shortGapX, match) + hmm_getTransition(hmm, shortGapY, match))/2); //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM5->TRANSITION_MATCH_FROM_LONG_GAP=log((hmm_getTransition(hmm, longGapX, match) + hmm_getTransition(hmm, longGapY, match))/2); //1.0 - gapExtend = 0.00343657420938
    sM5->TRANSITION_GAP_SHORT_OPEN=log((hmm_getTransition(hmm, match, shortGapX) + hmm_getTransition(hmm, match, shortGapY))/2); //0.0129868352330243
    sM5->TRANSITION_GAP_SHORT_EXTEND=log((hmm_getTransition(hmm, shortGapX, shortGapX) + hmm_getTransition(hmm, shortGapY, shortGapY))/2); //0.7126062401851738f;
    sM5->TRANSITION_GAP_SHORT_SWITCH=log((hmm_getTransition(hmm, shortGapX, shortGapY) + hmm_getTransition(hmm, shortGapY, shortGapX))/2); //0.0073673675173412815f;
    sM5->TRANSITION_GAP_LONG_OPEN=log((hmm_getTransition(hmm, match, longGapX) + hmm_getTransition(hmm, match, longGapY))/2); //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    sM5->TRANSITION_GAP_LONG_EXTEND=log((hmm_getTransition(hmm, longGapX, longGapX) + hmm_getTransition(hmm, longGapY, longGapY))/2);
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
    StateMachine3 *sM3 = (StateMachine3 *)sM;
    state_check(sM, state);
    switch(state) {
        case match :
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
    StateMachine3 *sM3 = (StateMachine3 *)sM;
    state_check(sM, state);
    switch(state) {
        case match :
            return (sM3->TRANSITION_GAP_OPEN_X + sM3->TRANSITION_GAP_OPEN_Y)/2.0;
        case shortGapX:
            return sM3->TRANSITION_GAP_EXTEND_X;
        case shortGapY:
            return sM3->TRANSITION_GAP_EXTEND_Y;
    }
    return 0.0;
}

static void stateMachine3_cellCalculate(StateMachine *sM, double *current,
        double *lower, double *middle, double *upper, Symbol cX, Symbol cY,
        void(*doTransition)(double *, double *, int64_t, int64_t, double, double, void *), void *extraArgs) {
    StateMachine3 *sM3 = (StateMachine3 *)sM;
    if (lower != NULL) {
        double eP = symbol_gapProb(cX);
        doTransition(lower, current, match, shortGapX, eP, sM3->TRANSITION_GAP_OPEN_X, extraArgs);
        doTransition(lower, current, shortGapX, shortGapX, eP, sM3->TRANSITION_GAP_EXTEND_X, extraArgs);
        doTransition(lower, current, shortGapY, shortGapX, eP, sM3->TRANSITION_GAP_SWITCH_TO_X, extraArgs);
    }
    if (middle != NULL) {
        double eP = symbol_matchProb(cX, cY);
        doTransition(middle, current, match, match, eP, sM3->TRANSITION_MATCH_CONTINUE, extraArgs);
        doTransition(middle, current, shortGapX, match, eP, sM3->TRANSITION_MATCH_FROM_GAP_X, extraArgs);
        doTransition(middle, current, shortGapY, match, eP, sM3->TRANSITION_MATCH_FROM_GAP_Y, extraArgs);

    }
    if (upper != NULL) {
        double eP = symbol_gapProb(cY);
        doTransition(upper, current, match, shortGapY, eP, sM3->TRANSITION_GAP_OPEN_Y, extraArgs);
        doTransition(upper, current, shortGapY, shortGapY, eP, sM3->TRANSITION_GAP_EXTEND_Y, extraArgs);
        doTransition(upper, current, shortGapX, shortGapY, eP, sM3->TRANSITION_GAP_SWITCH_TO_Y, extraArgs);
    }
}

StateMachine *stateMachine3_construct(StateMachineType type) {
    StateMachine3 *sM3 = st_malloc(sizeof(StateMachine3));
    sM3->TRANSITION_MATCH_CONTINUE=-0.030064059121770816; //0.9703833696510062f
    sM3->TRANSITION_MATCH_FROM_GAP_X=-1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM3->TRANSITION_MATCH_FROM_GAP_Y=-1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    sM3->TRANSITION_GAP_OPEN_X=-4.21256642; //0.0129868352330243
    sM3->TRANSITION_GAP_OPEN_Y=-4.21256642; //0.0129868352330243
    sM3->TRANSITION_GAP_EXTEND_X=-0.3388262689231553; //0.7126062401851738f;
    sM3->TRANSITION_GAP_EXTEND_Y=-0.3388262689231553; //0.7126062401851738f;
    sM3->TRANSITION_GAP_SWITCH_TO_X=-4.910694825551255; //0.0073673675173412815f;
    sM3->TRANSITION_GAP_SWITCH_TO_Y=-4.910694825551255; //0.0073673675173412815f;

    if(type != threeState && type != threeStateAsymmetric) {
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

    return (StateMachine *)sM3;
}

static void stateMachine3_loadAsymmetric(StateMachine3 *sM3, Hmm *hmm) {
    if(hmm->type != threeStateAsymmetric) {
        st_errAbort("Wrong hmm type");
    }
    sM3->TRANSITION_MATCH_CONTINUE=log(hmm_getTransition(hmm, match, match));
    sM3->TRANSITION_MATCH_FROM_GAP_X=log(hmm_getTransition(hmm, shortGapX, match));
    sM3->TRANSITION_MATCH_FROM_GAP_Y=log(hmm_getTransition(hmm, shortGapY, match));
    sM3->TRANSITION_GAP_OPEN_X=log(hmm_getTransition(hmm, match, shortGapX));
    sM3->TRANSITION_GAP_OPEN_Y=log(hmm_getTransition(hmm, match, shortGapY));
    sM3->TRANSITION_GAP_EXTEND_X=log(hmm_getTransition(hmm, shortGapX, shortGapX));
    sM3->TRANSITION_GAP_EXTEND_Y=log(hmm_getTransition(hmm, shortGapY, shortGapY));
    sM3->TRANSITION_GAP_SWITCH_TO_X=log(hmm_getTransition(hmm, shortGapY, shortGapX));
    sM3->TRANSITION_GAP_SWITCH_TO_Y=log(hmm_getTransition(hmm, shortGapX, shortGapY));
}

static void stateMachine3_loadSymmetric(StateMachine3 *sM3, Hmm *hmm) {
    if(hmm->type != threeState) {
        st_errAbort("Wrong hmm type");
    }
    sM3->TRANSITION_MATCH_CONTINUE=log(hmm_getTransition(hmm, match, match));
    sM3->TRANSITION_MATCH_FROM_GAP_X=log((hmm_getTransition(hmm, shortGapX, match)+hmm_getTransition(hmm, shortGapY, match))/2.0);
    sM3->TRANSITION_MATCH_FROM_GAP_Y=sM3->TRANSITION_MATCH_FROM_GAP_X;
    sM3->TRANSITION_GAP_OPEN_X=log((hmm_getTransition(hmm, match, shortGapX)+hmm_getTransition(hmm, match, shortGapY))/2.0);
    sM3->TRANSITION_GAP_OPEN_Y=sM3->TRANSITION_GAP_OPEN_X;
    sM3->TRANSITION_GAP_EXTEND_X=log((hmm_getTransition(hmm, shortGapX, shortGapX)+hmm_getTransition(hmm, shortGapY, shortGapY))/2.0);
    sM3->TRANSITION_GAP_EXTEND_Y=sM3->TRANSITION_GAP_EXTEND_X;
    sM3->TRANSITION_GAP_SWITCH_TO_X=log((hmm_getTransition(hmm, shortGapY, shortGapX)+hmm_getTransition(hmm, shortGapX, shortGapY))/2.0);
    sM3->TRANSITION_GAP_SWITCH_TO_Y=sM3->TRANSITION_GAP_SWITCH_TO_X;
}

///////////////////////////////////
///////////////////////////////////
//Public functions
///////////////////////////////////
///////////////////////////////////

StateMachine *hmm_getStateMachine(Hmm *hmm) {
    if(hmm->type == fiveState) {
        StateMachine5 *sM5 = (StateMachine5 *)stateMachine5_construct();
        stateMachine5_load(sM5, hmm);
        return (StateMachine *)sM5;
    }
    if(hmm->type == threeStateAsymmetric) {
        StateMachine3 *sM3 = (StateMachine3 *)stateMachine3_construct(hmm->type);
        stateMachine3_loadAsymmetric(sM3, hmm);
        return (StateMachine *)sM3;
    }
    if(hmm->type == threeState) {
        StateMachine3 *sM3 = (StateMachine3 *)stateMachine3_construct(hmm->type);
        stateMachine3_loadSymmetric(sM3, hmm);
        return (StateMachine *)sM3;
    }
    return NULL;
}

void stateMachine_destruct(StateMachine *stateMachine) {
    free(stateMachine);
}
