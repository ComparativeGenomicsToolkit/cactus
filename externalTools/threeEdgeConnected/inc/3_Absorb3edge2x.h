/*
 * Copyright (C) 2011 by Yung H Tsin and Nima Norouzi
 *
 * The authors of this code have given their permission
 * to distribute and redistribute this code within this software product,
 * however they require that any further incorporation or copying
 * be permitted by request only, hence the BSD/MIT license which
 * applies to much of the other code in this software product does not apply to this file.
 */

#ifndef ABSORB_3_EDGE_2X_H_
#define ABSORB_3_EDGE_2X_H_

/*
 * Function takes an graph represented as an adjacency list. Each vertex is represented
 * by a list in the argument "vertices", it's index in the list is its identifier.
 * Each vertice's list contain 1 length stIntTuples that indicate what the vertex is
 * connected to. The return value is a list of lists of nodes,
 * also represented using stLists and stIntTuples.
 */
stList *computeThreeEdgeConnectedComponents(stList *vertices);

#endif
