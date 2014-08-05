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

/*
 * Emissions
 */

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

/*
 * Transitions
 */

//Transitions
static double TRANSITION_MATCH_CONTINUE=-0.030064059121770816; //0.9703833696510062f
static double TRANSITION_MATCH_FROM_SHORT_GAP=-1.272871422049609; //1.0 - gapExtend - gapSwitch = 0.280026392297485
static double TRANSITION_MATCH_FROM_LONG_GAP=-5.673280173170473; //1.0 - gapExtend = 0.00343657420938
static double TRANSITION_GAP_SHORT_OPEN=-4.34381910900448; //0.0129868352330243
static double TRANSITION_GAP_SHORT_EXTEND=-0.3388262689231553; //0.7126062401851738f;
static double TRANSITION_GAP_SHORT_SWITCH=-4.910694825551255; //0.0073673675173412815f;
static double TRANSITION_GAP_LONG_OPEN=-6.30810595366929; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
static double TRANSITION_GAP_LONG_EXTEND=-0.003442492794189331; //0.99656342579062f;

typedef enum {
    match=0,
    shortGapX=1,
    shortGapY=2,
    longGapX=3,
    longGapY=4
} State;

void loadTheGlobalHmm(Hmm *hmm) {
    TRANSITION_MATCH_CONTINUE=hmm_tProb(hmm, match, match); //0.9703833696510062f
    TRANSITION_MATCH_FROM_SHORT_GAP=(hmm_tProb(hmm, shortGapX, match) + hmm_tProb(hmm, shortGapY, match))/2; //1.0 - gapExtend - gapSwitch = 0.280026392297485
    TRANSITION_MATCH_FROM_LONG_GAP=(hmm_tProb(hmm, longGapX, match) + hmm_tProb(hmm, longGapY, match))/2; //1.0 - gapExtend = 0.00343657420938
    TRANSITION_GAP_SHORT_OPEN=(hmm_tProb(hmm, match, shortGapX) + hmm_tProb(hmm, match, shortGapY))/2; //0.0129868352330243
    TRANSITION_GAP_SHORT_EXTEND=(hmm_tProb(hmm, shortGapX, shortGapX) + hmm_tProb(hmm, shortGapY, shortGapY))/2; //0.7126062401851738f;
    TRANSITION_GAP_SHORT_SWITCH=(hmm_tProb(hmm, shortGapX, shortGapY) + hmm_tProb(hmm, shortGapY, shortGapX))/2; //0.0073673675173412815f;
    TRANSITION_GAP_LONG_OPEN=(hmm_tProb(hmm, match, longGapX) + hmm_tProb(hmm, match, longGapY))/2; //(1.0 - match - 2*gapOpenShort)/2 = 0.001821479941473
    TRANSITION_GAP_LONG_EXTEND=(hmm_tProb(hmm, longGapX, longGapX) + hmm_tProb(hmm, longGapY, longGapY))/2;
}

static void state_check(State s) {
    assert(s >= 0 && s < STATE_NUMBER);
}

double state_startStateProb(int64_t state) {
    //Match state is like going to a match.
    state_check(state);
    return state == match ? 0 : LOG_ZERO;
}

double state_raggedStartStateProb(int64_t state) {
    state_check(state);
    return (state == longGapX || state == longGapY) ? 0 : LOG_ZERO;
}

double state_endStateProb(int64_t state) {
    //End state is like to going to a match
    state_check(state);
    switch(state) {
        case match :
            return TRANSITION_MATCH_CONTINUE;
        case shortGapX:
        case shortGapY:
            return TRANSITION_MATCH_FROM_SHORT_GAP;
        case longGapX:
        case longGapY:
            return TRANSITION_MATCH_FROM_LONG_GAP;
    }
    return 0.0;
}

double state_raggedEndStateProb(int64_t state) {
    //End state is like to going to a match
    state_check(state);
    switch(state) {
        case match :
            return TRANSITION_GAP_LONG_OPEN;
        case shortGapX:
        case shortGapY:
            return TRANSITION_GAP_LONG_OPEN;
        case longGapX:
        case longGapY:
            return TRANSITION_GAP_LONG_EXTEND;
    }
    return 0.0;
}

void cell_calculate(double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY,
        void(*doTransition)(double *, double *, int64_t, int64_t, double, double, void *), void *extraArgs) {
    if (lower != NULL) {
        double eP = symbol_gapProb(cX);
        doTransition(lower, current, match, shortGapX, eP, TRANSITION_GAP_SHORT_OPEN, extraArgs);
        doTransition(lower, current, shortGapX, shortGapX, eP, TRANSITION_GAP_SHORT_EXTEND, extraArgs);
        doTransition(lower, current, shortGapY, shortGapX, eP, TRANSITION_GAP_SHORT_SWITCH, extraArgs);
        doTransition(lower, current, match, longGapX, eP, TRANSITION_GAP_LONG_OPEN, extraArgs);
        doTransition(lower, current, longGapX, longGapX, eP, TRANSITION_GAP_LONG_EXTEND, extraArgs);
    }
    if (middle != NULL) {
        double eP = symbol_matchProb(cX, cY);
        doTransition(middle, current, match, match, eP, TRANSITION_MATCH_CONTINUE, extraArgs);
        doTransition(middle, current, shortGapX, match, eP, TRANSITION_MATCH_FROM_SHORT_GAP, extraArgs);
        doTransition(middle, current, shortGapY, match, eP, TRANSITION_MATCH_FROM_SHORT_GAP, extraArgs);
        doTransition(middle, current, longGapX, match, eP, TRANSITION_MATCH_FROM_LONG_GAP, extraArgs);
        doTransition(middle, current, longGapY, match, eP, TRANSITION_MATCH_FROM_LONG_GAP, extraArgs);
    }
    if (upper != NULL) {
        double eP = symbol_gapProb(cY);
        doTransition(upper, current, match, shortGapY, eP, TRANSITION_GAP_SHORT_OPEN, extraArgs);
        doTransition(upper, current, shortGapY, shortGapY, eP, TRANSITION_GAP_SHORT_EXTEND, extraArgs);
        doTransition(upper, current, shortGapX, shortGapY, eP, TRANSITION_GAP_SHORT_SWITCH, extraArgs);
        doTransition(upper, current, match, longGapY, eP, TRANSITION_GAP_LONG_OPEN, extraArgs);
        doTransition(upper, current, longGapY, longGapY, eP, TRANSITION_GAP_LONG_EXTEND, extraArgs);
    }
}
