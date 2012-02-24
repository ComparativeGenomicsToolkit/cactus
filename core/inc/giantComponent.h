/*
 * giantComponent.h
 *
 *  Created on: 22 Feb 2012
 *      Author: benedictpaten
 */

#ifndef GIANTCOMPONENT_H_
#define GIANTCOMPONENT_H_

/*
 * Nodes is a list of integers representing the nodes.
 * Each edge is represented as an int tuple (weight, vertex1, vertex2).
 * Returns a sublist of the edges in edges that must deleted, so that the size of the largest component in the graph
 * is smaller than maxComponentSize.
 */
stList *breakupComponentGreedily(stList *nodes, stList *edges, int32_t maxComponentSize);

#endif /* GIANTCOMPONENT_H_ */
