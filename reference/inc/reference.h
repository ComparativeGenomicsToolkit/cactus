/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * reference.h
 *
 *  Created on: 1 Apr 2010
 *      Author: benedictpaten
 *
 * Algorithms for building references for the cactus structure.
 */

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include "cactus.h"
#include "matchingAlgorithms.h"

extern const char *REFERENCE_BUILDING_EXCEPTION;

/*
 * The initial top down algorithm, which is run for each flower in the
 * cactus tree breadth first from the root.
 */
void constructReference_topDownPhase(Flower *flower, MatchingAlgorithm matchingAlgorithm);

/*
 * The finishing bottom up algorithm, which is run for each flower in the
 * cactus tree bottom up from the leaves.
 */
void constructReference_bottomUpPhase(Flower *flower);

#endif /* REFERENCE_H_ */
