/*
 * pickAdjacencyPairs.h
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#ifndef PICKADJACENCYPAIRS_H_
#define PICKADJACENCYPAIRS_H_

#include "sonLib.h"
#include "cactus.h"

/*
 * Creates a set of adjacency pairs such that every end in the flower
 * is a member of one adjacency pair. Adjacency pairs are only between ends
 * in the same group.
 */
stHash *pickAdjacencyPairs(Flower *flower, Reference *reference);

#endif /* PICKADJACENCYPAIRS_H_ */
