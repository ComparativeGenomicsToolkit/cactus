/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * reference.h
 *
 *  Created on: 1 Apr 2010
 *      Author: benedictpaten
 *
 * Algorithms for building references for the cactus structure.
 */

#ifndef SHARED_H_
#define SHARED_H_

#include "sonLib.h"

stIntTuple *constructEdge(int32_t node1, int32_t node2);

stIntTuple *constructWeightedEdge(int32_t node1, int32_t node2, int32_t weight);

int compareEdgesByWeight(const void *edge, const void *edge2);

int32_t getOtherPosition(stIntTuple *edge, int32_t node);

stSortedSet *getNodeSetOfEdges(stList *edges);

void addNodeToSet(stSortedSet *nodes, int32_t node);

bool nodeInSet(stSortedSet *nodes, int32_t node);

stIntTuple *getWeightedEdgeFromSet(int32_t node1, int32_t node2,
        stSortedSet *allAdjacencyEdges);

stHash *getNodesToEdgesHash(stList *edges);

stIntTuple *getEdgeForNodes(int32_t node1, int32_t node2,
        stHash *nodesToAdjacencyEdges);

void *getItemForNode(int32_t node, stHash *nodesToItems);

bool edgeInSet(stSortedSet *edges, int32_t node1, int32_t node2);

void addEdgeToSet(stSortedSet *edges, int32_t node1, int32_t node2);

#endif /* REFERENCE_H_ */
