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

/*
 * Adds sequence objects and coordinates for each reference.
 */
void addReferenceSequences(Flower *flower, Name referenceEventName, Event *outgroupEvent);

/*
 * Traverses caps in order, along a sequence, starting from a 3 prime cap.
 * The 3' function is called with an ordered list, from highest to lowest, of all copies of a given 3' cap.
 * The 5' function is called with an ordered list, from highest to lowest, of all copies of a given 5' cap.
 */
void traverseCapsInSequenceOrderFrom3PrimeCap(Cap *cap, void *extraArg,
        void(*_3PrimeFn)(stList *caps, void *extraArg),
        void(*_5PrimeFn)(stList *caps, void *extraArg));

/*
 * Gets a cap for a given event string.
 */
Cap *getCapForReferenceEvent(End *end, Name referenceEventName);

double *calculateZ(Flower *flower, stHash *endsToNodes, double theta);

#endif /* REFERENCE_H_ */
