/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef TRAVERSE_FLOWERS_H_
#define TRAVERSE_FLOWERS_H_

#include "sonLib.h"
#include "cactus.h"

/*
 * Used in bar recursion to recursively find all alignment subproblems
 * in the hierarchy for bar to complete.
 */
void extendFlowers(Flower *flower, stList *extendedFlowers, int64_t minFlowerSize);

/*
 * Returns a list of lists of the flowers, where the flowers in layer i are the children of those
 * in layer i-1 and the parents of those in layer i+1.
 */
stList *getFlowerHierarchyInLayers(Flower *rootFlower);

#endif /* TRAVERSE_FLOWERS_H_ */

