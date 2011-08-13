/*
 * adjacencyProblem.h
 *
 *  Created on: 11 Aug 2011
 *      Author: benedictpaten
 */

#ifndef ADJACENCYPROBLEM_H_
#define ADJACENCYPROBLEM_H_

stList *makeReferenceGreedily(stList *stubs, stList *chains,
        double *z, int32_t nodeNumber, double *totalScore);

void gibbsSamplingWithSimulatedAnnealing(stList *reference,
        stList *chains, double *z, int32_t permutations,
        double(*temperature)(double), bool pureGreedy);

stList *convertReferenceToAdjacencyEdges(stList *reference);

double exponentiallyDecreasingTemperatureFn(double d);

double constantTemperatureFn(double d);

void logReference(stList *reference, int32_t nodeNumber, double *zMatrix, double totalScore, const char *message);

double calculateZScoreOfReference(stList *reference, int32_t nodeNumber, double *zMatrix);

double calculateMaxZ(int32_t nodeNumber, double *zMatrix);

double calculateZScore(int32_t n, int32_t m, int32_t k, double theta);

#endif /* ADJACENCYPROBLEM_H_ */
