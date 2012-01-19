/*
 * addReferenceCoordinates.h
 *
 *  Created on: 17 Jan 2012
 *      Author: benedictpaten
 */

#ifndef ADDREFERENCECOORDINATES_H_
#define ADDREFERENCECOORDINATES_H_

Cap *getCapForReferenceEvent(End *end, Name referenceEventName);

void bottomUp(Flower *flower, Name referenceEventName);

void addSequencesAndReferenceCoordinatesToTopLevelFlower(Flower *flower, Name referenceEventName, Name outgroupEventName);

void topDown(Flower *flower, Name referenceEventName, Name outgroupEventName);

#endif /* ADDREFERENCECOORDINATES_H_ */
