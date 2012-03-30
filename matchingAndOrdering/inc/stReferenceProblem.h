/*
 * adjacencyProblem.h
 *
 *  Created on: 11 Aug 2011
 *      Author: benedictpaten
 */

#ifndef REFERENCEPROBLEM_H_
#define REFERENCEPROBLEM_H_

stList *makeReferenceGreedily(stList *stubs, stList *chains,
        double *z, int32_t nodeNumber, double *totalScore, bool fast);

void gibbsSamplingWithSimulatedAnnealing(stList *reference,
        stList *chains, double *z, int32_t permutations,
        double(*temperature)(double), bool pureGreedy);

stList *convertReferenceToAdjacencyEdges(stList *reference);

double exponentiallyDecreasingTemperatureFn(double d);

double constantTemperatureFn(double d);

void logReference(stList *reference, int32_t nodeNumber, double *zMatrix, double totalScore, const char *message);

double calculateZScoreOfReference(stList *reference, int32_t nodeNumber, double *zMatrix);

void logZScoreOfReference(stList *reference, int32_t nodeNumber, double *zMatrix);

double calculateMaxZ(int32_t nodeNumber, double *zMatrix);

double calculateZScore(int32_t n, int32_t m, int32_t k, double theta);

#endif /* REFERENCEPROBLEM_H_ */
