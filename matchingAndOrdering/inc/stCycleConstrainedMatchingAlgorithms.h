/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * cycleConstrainedMatchingAlgorithms.h
 *
 *  Created on: 10 Aug 2010
 *      Author: benedictpaten
 */

#ifndef CYCLE_CONSTRAINED_MATCHING_ALGORITHMS_H_
#define CYCLE_CONSTRAINED_MATCHING_ALGORITHMS_H_

#include "sonLib.h"

/*
 * There are three types of edge, adjacency edges, stub edges and chain edges.
 * Each node is incident with either a single chain or stub edge. There are therefore an even number of nodes.
 * The adjacency edges have an associated integer valued weight >= 0, and the adjacency edges form a clique.
 * The algorithm finds a perfect matching of the nodes, such that the weight of the adjacency edges in
 * the matching is maximal and, including both the adjacency edges in the matching and the chain/stub edges,
 * there are (1) no cycles including only chain edges, and, if make stub cycles disjoint is true, 1 stub
 * edge per cycle.
 *
 * The algorithm first computes an maximum weight / cardinality perfect matching (see matching algorithm parameter),
 * then the constraints are imposed.
 *
 * Node number is the number of nodes such that nodeNumber%2 == 0.
 *
 * Adjacency edges are are stIntTuple's of
 * length 3 of the form (node1, node2, weight), where node1 and node2 are indices
 * greater than or equal to zero and less than the total node number and weight
 * is a positive integer weight.
 *
 * Stub and chain edges are stIntTuple's of
 * length 2 of the form (node1, node2).
 *
 * The matching algorithm is one used to construct an initial matching.
 *
 * The return is the matching, as a list of node pairs, of the same form as the stub end adjacency edges.
 * The list of adjacency edges may be added to with new 'pseudo' adjacency edges, if necessary.
 */
stList *getMatchingWithCyclicConstraints(stSortedSet *nodes,
        stList *adjacencyEdges, stList *stubEdges, stList *chainEdges,
        bool makeStubCyclesDisjoint,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber));

stList *makeMatchingObeyCyclicConstraints(stSortedSet *nodes,
        stList *chosenEdges,
        stSortedSet *allAdjacencyEdges, stList *nonZeroWeightAdjacencyEdges,
        stList *stubEdges, stList *chainEdges,
        bool makeStubCyclesDisjoint);

#endif /* CYCLE_CONSTRAINED_MATCHING_ALGORITHMS_H_ */

