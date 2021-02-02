/*
 * addReferenceCoordinates.h
 *
 *  Created on: 17 Jan 2012
 *      Author: benedictpaten
 */

#ifndef ADDREFERENCECOORDINATES_H_
#define ADDREFERENCECOORDINATES_H_

#include "cactus.h"
#include "recursiveThreadBuilder.h"

Cap *getCapForReferenceEvent(End *end, Name referenceEventName);

void bottomUp(Flower *flower, stKVDatabase *sequenceDatabase, Name referenceEventName, bool isTop, stMatrix *(*generateSubstitutionMatrix)(double));

void bottomUpNoDb(Flower *flower, RecordHolder *rh, Name referenceEventName,
                  bool isTop, stMatrix *(*generateSubstitutionMatrix)(double));

void topDown(Flower *flower, Name referenceEventName);

#endif /* ADDREFERENCECOORDINATES_H_ */
