/*
 * hal.h
 *
 *  Created on: 21 Jun 2012
 *      Author: benedictpaten
 */

#ifndef HAL_H_
#define HAL_H_

#include "sonLib.h"
#include "cactus.h"

void makeMaf(Flower *flower, RecursiveFileBuilder *recursiveFileBuilder,
        Event *referenceEvent, FILE *parentFileHandle,
        bool showOnlySubstitutionsWithRespectToReference, bool hasParent);

void makeHalFormat(stList *flowers, stKVDatabase *database, Name referenceEventName,
        FILE *fileHandle);

#endif /* HAL_H_ */
