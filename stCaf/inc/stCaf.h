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
#include "stPinchIterator.h"
#include "stCactusGraphs.h"
#include "cactus.h"

/*
 * Create a pinch graph from a flower containing no blocks. Each thread is given the name of its 5' cap,
 * and stubs are represented by length 1 blocks.
 */
stPinchThreadSet *stCaf_constructEmptyPinchGraph(Flower *flower);

/*
 * Add the a set of alignments, represented as pinches, to the graph.
 */
void stCaf_addAlignmentsToPinchGraph(stPinchThreadSet *threadSet, stPinchIterator *pinchIterator);

/*
 * Locates the ends of all the attached ends and merges together their 'dead end' components to create a single
 * 'dead end' component, as described in the JCB cactus paper.
 */
stList *stCaf_constructDeadEndComponent(Flower *flower, stPinchThreadSet *threadSet,
        stHash *pinchEndsToAdjacencyComponents);

/*
 * Locates threads components which have no dead ends part of the dead end component, and then
 * connects them, picking the longest thread to attach them.
 */
void stCaf_attachUnattachedThreadComponents(Flower *flower, stPinchThreadSet *threadSet, stList *deadEndComponent,
        stHash *pinchEndsToAdjacencyComponents, bool markEndsAttached);

/*
 * Constructs a cactus graph from a set of pinch graph components, including the dead end component. Returns a cactus
 * graph, and assigns 'startCactusNode' to the cactus node containing the dead end component.
 */
stCactusGraph *stCaf_constructCactusGraph(stList *deadEndComponent,
        stHash *pinchEndsToAdjacencyComponents, stCactusNode **startCactusNode);

/*
 * Converts a given cactus graph/pinch graph into the cactus datastructure.
 */
void stCaf_convertCactusGraphToFlowers(stPinchThreadSet *threadSet, stCactusNode *startCactusNode, Flower *parentFlower,
        stList *deadEndComponent);

/*
 * Add the adjacencies between caps of a flower.
 * All ends must be in the flower before this function is called.
 */
void stCaf_addAdjacencies(Flower *flower);

/*
 * Basic function for filling out a flower with a given set of alignments into a cactus graph.
 */
void stCaf_core(Flower *flower, stPinchIterator *pinchIterator);

#endif /* STCAF_H_ */
