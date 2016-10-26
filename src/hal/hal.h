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

void makeHalFormat(stList *flowers, stKVDatabase *database, Name referenceEventName,
        FILE *fileHandle);

void printFastaSequences(Flower *flower, FILE *fileHandle, Name referenceEventName);

#endif /* HAL_H_ */
