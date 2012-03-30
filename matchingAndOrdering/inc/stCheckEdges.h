/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 *
 *  Created on: 6 July 2011
 *      Author: benedictpaten
 *
 */

#ifndef CHECK_EDGES_H_
#define CHECK_EDGES_H_

#include "sonLib.h"

/*
 * Check the edges all refer to nodes between 0 and node number. If length three, check final argument
 * (which is a weight), is greater than or equal to zero. Check that if coversAllNodes is non-zero, all nodes are
 * covered. Checks edges from clique if isClique is non-zero.
 */
void checkEdges(stList *edges, stSortedSet * nodes, bool coversAllNodes,
        bool isClique);

#endif
