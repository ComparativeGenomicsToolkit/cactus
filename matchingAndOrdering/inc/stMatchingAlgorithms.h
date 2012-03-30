/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * externalAlgorithms.h
 *
 *  Created on: 7 Sep 2010
 *      Author: benedictpaten
 */

#ifndef EXTERNALALGORITHMS_H_
#define EXTERNALALGORITHMS_H_

#include "sonLib.h"

extern const char *MATCHING_EXCEPTION;

stIntTuple *constructEdge(int32_t node1, int32_t node2);

stIntTuple *constructWeightedEdge(int32_t node1, int32_t node2, int32_t weight);

/*
 * Uses the blossom5 maximum weight perfect matching algorithm to choose a matching
 * between the edges. The returned matching is not necessarily perfect, rather extra edges are
 * allowed to make the graph all connected. Edges are stIntTuple's of
 * length 3 of the form (node1, node2, weight), where node1 and node2 are indices
 * greater than or equal to zero and less than the total node number and weight
 * is a positive integer weight.
 */
stList *chooseMatching_blossom5(stList *edges, int32_t nodeNumber);

/*
 * Finds maximal matching (maximum cardinality) matching. Same form as the blossom algorithm.
 */
stList *chooseMatching_maximumCardinalityMatching(stList *edges, int32_t nodeNumber);

/*
 * Finds maximum weight matching. Same form as the blossum algorithm.
 */
stList *chooseMatching_maximumWeightMatching(stList *edges, int32_t nodeNumber);

/*
 * Uses a greedy algorithm to choose the matching, starting with the highest weight pair
 * edges are selected in descending order of weight until no further edges can be added to the matching.
 *
 * Same form as the blossom algorithm.
 */
stList *chooseMatching_greedy(stList *edges, int32_t nodeNumber);

/*
 * Returns number of edges with weight > 0.
 */
int32_t matchingCardinality(stList *matching);

/*
 * Returns sum of weights.
 */
int32_t matchingWeight(stList *matching);

stList *getComponents(stList *edges);

#endif /* EXTERNALALGORITHMS_H_ */
