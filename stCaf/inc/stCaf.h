/*
 * stCore.h
 *
 *  Created on: 28 Apr 2012
 *      Author: benedictpaten
 */

#ifndef STCAF_H_
#define STCAF_H_

#include "sonLib.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "cactus.h"

stThreadSet *stCaf_constructEmptyPinchGraph(Flower *flower);

void stCaf_addAlignmentsToPinchGraph(stThreadSet *threadSet, stPinch *(*pinchGenerator)());

stList *stCaf_constructDeadEndComponent(Flower *flower, stThreadSet *threadSet, stSortedSet *adjacencyComponents, stHash *edgeEndsToAdjacencyComponents);

void stCaf_attachUnattachedThreadComponents(Flower *flower, stThreadSet *threadSet, stList *deadEndComponent,
        stSortedSet *adjacencyComponents, stHash *edgeEndsToAdjacencyComponents, bool markEndsAttached);

stCactusGraph *stCaf_constructCactusGraph(stSortedSet *adjacencyComponents, stList *deadEndComponent,
        stHash *edgeEndsToAdjacencyComponents, stCactusNode **startCactusNode);

void stCaf_convertCactusGraphToFlowers(stThreadSet *threadSet, stCactusNode *startCactusNode, Flower *parentFlower, stList *deadEndComponent);

#endif /* STCAF_H_ */
