/*
 * addReferenceCoordinates.h
 *
 *  Created on: 17 Jan 2012
 *      Author: benedictpaten
 */

#ifndef ADDREFERENCECOORDINATES_H_
#define ADDREFERENCECOORDINATES_H_

Cap *getCapForReferenceEvent(End *end, Name referenceEventName);

void bottomUp(stList *flowers, stKVDatabase *sequenceDatabase, Name referenceEventName, Name outgroupEventName, bool isTop);

void topDown(stList *flowers, Name referenceEventName);

#endif /* ADDREFERENCECOORDINATES_H_ */
