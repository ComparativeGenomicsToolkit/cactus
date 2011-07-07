/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 *
 *  Created on: 6 July 2011
 *      Author: benedictpaten
 *
 */

#ifndef PERFECT_MATCHING_H_
#define PERFECT_MATCHING_H_

#include "sonLib.h"

/*
 * Computes a perfect matching. The adjacency edges must be a clique.
 */
stList *getPerfectMatching(stSortedSet *nodes,
        stList *adjacencyEdges,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber));

#endif
