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

stPinchThreadSet *stCaf_constructEmptyPinchGraph(Flower *flower);

void stCaf_addAlignmentsToPinchGraph(stPinchThreadSet *threadSet, stPinch *(*pinchGenerator)());

stList *stCaf_constructDeadEndComponent(Flower *flower, stPinchThreadSet *threadSet, stSortedSet *adjacencyComponents, stHash *edgeEndsToAdjacencyComponents);

void stCaf_attachUnattachedThreadComponents(Flower *flower, stPinchThreadSet *threadSet, stList *deadEndComponent,
        stSortedSet *adjacencyComponents, stHash *edgeEndsToAdjacencyComponents, bool markEndsAttached);

stCactusGraph *stCaf_constructCactusGraph(stSortedSet *adjacencyComponents, stList *deadEndComponent,
        stHash *edgeEndsToAdjacencyComponents, stCactusNode **startCactusNode);

void stCaf_convertCactusGraphToFlowers(stPinchThreadSet *threadSet, stCactusNode *startCactusNode, Flower *parentFlower, stList *deadEndComponent);

#endif /* STCAF_H_ */
