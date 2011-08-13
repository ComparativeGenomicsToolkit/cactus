/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 *
 *  Created on: 6 July 2011
 *      Author: benedictpaten
 *
 */

#ifndef SPARSE_MATCHING_H_
#define SPARSE_MATCHING_H_

#include "sonLib.h"

/*
 * Gets a matching for the set of nodes, which are represented as integers in the set nodes, for the
 * given matching algorithm.
 */
stList *getSparseMatching(stSortedSet *nodes,
        stList *adjacencyEdges,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber));

#endif
