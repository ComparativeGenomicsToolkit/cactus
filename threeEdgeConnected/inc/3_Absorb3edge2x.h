/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef ABSORB_3_EDGE_2X_H_
#define ABSORB_3_EDGE_2X_H_

/*
 * Function takes an graph represented as an adjacency list. Each vertex is represented
 * by a list in the argument "vertices", it's index in the list is its identifier.
 * Each vertice's list contain 1 length stIntTuples that indicate what the vertex is
 * connected to.
 */
stList *computeThreeEdgeConnectedComponents(stList *vertices);

#endif
