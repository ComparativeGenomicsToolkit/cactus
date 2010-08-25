/*
 * adjacencyPairComponents.h
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#ifndef ADJACENCYPAIRCOMPONENTS_H_
#define ADJACENCYPAIRCOMPONENTS_H_

#include "sonLib.h"
#include "cactus.h"

/*
 * Gets the connected component the end is part of, as a list of adjacency pairs.
 */
stList *adjacencyHash_getConnectedComponent(stHash *adjacencies, Flower *flower, End *end);

/*
 * Gets a list of connected components of ends linked by adjacencies and blocks.
 */
stList *adjacencyHash_getConnectedComponents(stHash *adjacencies, Flower *flower);

/*
 * Gets the one or two stub ends in the component and assigns them, first to end1, then
 * to end2. If there isn't an end, sets it null.
 */
void getAttachedStubEndsInComponent(stList *component, End **end1,
		End **end2);

#endif /* ADJACENCYPAIRCOMPONENTS_H_ */
