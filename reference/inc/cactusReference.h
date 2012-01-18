/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
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
#include "cactusMatchingAlgorithms.h"

extern const char *REFERENCE_BUILDING_EXCEPTION;

/*
 * Construct a reference for the flower, top down.
 */
void buildReferenceTopDown(Flower *flower, const char *referenceEventHeader,
        int32_t permutations,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber),
        double (*temperature)(double),
        double theta, bool maxNumberOfChainsBeforeSwitchingToFast);

double *calculateZ(Flower *flower, stHash *endsToNodes, double theta);

#endif /* REFERENCE_H_ */
