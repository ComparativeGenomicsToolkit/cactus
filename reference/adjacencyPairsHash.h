/*
 * adjacencyPairsHash.h
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#ifndef ADJACENCYPAIRSHASH_H_
#define ADJACENCYPAIRSHASH_H_

#include "sonLib.h"
#include "cactus.h"

/*
 * Adds the adjacency pair to the hash.
 */
void adjacencyHash_add(stHash *adjacenciesHash,
		AdjacencyPair *adjacencyPair);

/*
 * Removes the adjacency pair from the hash, it does not destroy it.
 */
void adjacencyHash_remove(stHash *adjacencies,
		AdjacencyPair *adjacencyPair);

/*
 * Frees the adjacencies pairs in the adjacencies hash safely.
 */
void adjacencyHash_cleanUp(stHash *adjacencies, Flower *flower);

#endif /* ADJACENCYPAIRSHASH_H_ */
