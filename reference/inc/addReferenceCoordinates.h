/*
 * addReferenceCoordinates.h
 *
 *  Created on: 17 Jan 2012
 *      Author: benedictpaten
 */

#ifndef ADDREFERENCECOORDINATES_H_
#define ADDREFERENCECOORDINATES_H_

Cap *getCapForReferenceEvent(End *end, Name referenceEventName);

void bottomUp(Flower *flower, const char *childSequenceDir, const char *sequenceDir,
        Name referenceEventName, Name outgroupEventName);

void addSequences(Flower *flower, const char *sequenceDir, Name referenceEventName);

void topDown(Flower *flower, Name referenceEventName);

#endif /* ADDREFERENCECOORDINATES_H_ */
