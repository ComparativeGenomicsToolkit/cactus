/*
 * giantComponent.h
 *
 *  Created on: 22 Feb 2012
 *      Author: benedictpaten
 */

#ifndef ST_GIANTCOMPONENT_H_
#define ST_GIANTCOMPONENT_H_

#include "sonLib.h"
#include "stPinchGraphs.h"

/*
 * Nodes is a list of integers representing the nodes.
 * Each edge is represented as an int tuple (weight, vertex1, vertex2).
 * Returns a sublist of the edges in edges that must deleted, so that the size of the largest component in the graph
 * is smaller than maxComponentSize.
 */
stList *stCaf_breakupComponentGreedily(stList *nodes, stList *edges, int64_t maxComponentSize);

/*
 * Break up component extra large compoonents greedily.
 */
void stCaf_breakupComponentsGreedily(stPinchThreadSet *threadSet, float maximumAdjacencyComponentSizeRatio);

#endif /* ST_GIANTCOMPONENT_H_ */
