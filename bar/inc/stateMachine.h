/*
 * stateMachine.h
 *
 *  Created on: Aug 8, 2014
 *      Author: benedictpaten
 */

#ifndef STATEMACHINE_H_
#define STATEMACHINE_H_

#include "sonLib.h"

#define SYMBOL_NUMBER 5
#define SYMBOL_NUMBER_NO_N 4

typedef enum {
    a=0,
    c=1,
    g=2,
    t=3,
    n=4
} Symbol;

/*
 * The statemachine object for computing pairwise alignments with.
 */

typedef enum {
    fiveState=0,
    fiveStateAsymmetric=1,
    threeState=2,
    threeStateAsymmetric=3
} StateMachineType;

typedef struct _stateMachine StateMachine;

struct _stateMachine {
    StateMachineType type;
    int64_t stateNumber;
    int64_t matchState;

    double (*startStateProb)(StateMachine *sM, int64_t state);

    double (*endStateProb)(StateMachine *sM, int64_t state);

    double (*raggedEndStateProb)(StateMachine *sM, int64_t state);

    double (*raggedStartStateProb)(StateMachine *sM, int64_t state);

    //Cells (states at a given coordinate(
    void (*cellCalculate)(StateMachine *sM, double *current, double *lower, double *middle, double *upper, Symbol cX, Symbol cY,
            void(*doTransition)(double *, double *, int64_t, int64_t, double, double, void *), void *extraArgs);
};

/*
 * Hmm for loading/unloading HMMs and storing expectations.
 */

typedef struct _hmm {
    StateMachineType type;
    double *transitions;
    double *emissions;
    double likelihood;
    int64_t stateNumber;
} Hmm;

Hmm *hmm_constructEmpty(double pseudoExpectation, StateMachineType type);

void hmm_randomise(Hmm *hmm); //Creates normalised HMM with parameters set to small random values.

void hmm_destruct(Hmm *hmmExpectations);

void hmm_write(Hmm *hmmExpectations, FILE *fileHandle);

void hmm_addToTransitionExpectation(Hmm *hmmExpectations, int64_t from, int64_t to, double p);

double hmm_getTransition(Hmm *hmmExpectations, int64_t from, int64_t to);

void hmm_setTransition(Hmm *hmm, int64_t from, int64_t to, double p);

void hmm_addToEmissionsExpectation(Hmm *hmmExpectations, int64_t state, Symbol x, Symbol y, double p);

double hmm_getEmissionsExpectation(Hmm *hmm, int64_t state, Symbol x, Symbol y);

void hmm_setEmissionsExpectation(Hmm *hmm, int64_t state, Symbol x, Symbol y, double p);

Hmm *hmm_loadFromFile(const char *fileName);

void hmm_normalise(Hmm *hmm);

StateMachine *hmm_getStateMachine(Hmm *hmm);

StateMachine *stateMachine5_construct(StateMachineType type);

StateMachine *stateMachine3_construct(StateMachineType type); //the type is to specify symmetric/asymmetric

void stateMachine_destruct(StateMachine *stateMachine);

#endif /* STATEMACHINE_H_ */
