/*
 * adjacencyPairs.h
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#ifndef ADJACENCYPAIRS_H_
#define ADJACENCYPAIRS_H_

#include "sonLib.h"
#include "cactus.h"

typedef struct _adjacencyPair AdjacencyPair;

/*
 * Constructs an adjacency pair. (an un-ordered pair of adjacencies).
 */
AdjacencyPair *adjacencyPair_construct(End *end1, End *end2);


/*
 * Gets the first end in the adjacency.
 */
End *adjacencyPair_getEnd1(AdjacencyPair *adjacencyPair);

/*
 * Gets the second end in the adjacency.
 */
End *adjacencyPair_getEnd2(AdjacencyPair *adjacencyPair);

/*
 * Gets the group the adjacency pair is in.
 */
Group *adjacencyPair_getGroup(AdjacencyPair *adjacencyPair);

/*
 * Gets the other end in the adjacency pair.
 */
End *adjacencyPair_getOtherEnd(AdjacencyPair *adjacencyPair, End *end);

/*
 * Destructs an adjacency pair.
 */
void adjacencyPair_destruct(AdjacencyPair *adjacencyPair);

/*
 * Returns the strength of the adjacency between the two ends.
 */
uint32_t adjacencyPair_getStrengthOfAdjacencyPair(AdjacencyPair *adjacencyPair);

/*
 * Hash function for adjacency pairs.
 */
uint32_t adjacencyPair_hashKey(AdjacencyPair *adjacencyPair);

/*
 * Hash equals function for adjacency pairs.
 */
int adjacencyPair_hashEqual(AdjacencyPair *adjacencyPair1,
		AdjacencyPair *adjacencyPair2);

#endif /* ADJACENCYPAIRS_H_ */
