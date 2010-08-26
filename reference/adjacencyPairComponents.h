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
 * Gets the connected component linked by adjacency pairs and hyper chains the end is part of.
 */
stList *adjacencyHash_getConnectedComponent(stHash *adjacencies, stHash *hyperChains, End *end);

/*
 * Gets a list of connected components of ends linked by adjacency pairs and hyper chains.
 */
stList *adjacencyHash_getConnectedComponents(stHash *adjacencies, stHash *hyperChains);

/*
 * Gets the one or two top stub ends in the component and assigns them, first to end1, then
 * to end2. If there isn't an end, sets it null.
 */
void getTopStubEndsInComponent(stList *component, stHash *hyperChains, End **end1,
		End **end2);

#endif /* ADJACENCYPAIRCOMPONENTS_H_ */
