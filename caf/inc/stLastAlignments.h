/*
 * stLastAlignments.h
 *
 *  Created on: 30 Jan 2012
 *      Author: benedictpaten
 */

#ifndef ST_LASTALIGNMENT_H_
#define ST_LASTALIGNMENT_H_

#include "cactus.h"
#include "sonLib.h"

stList *stCaf_selfAlignFlower(Flower *flower, int64_t minimumSequenceLength, const char *lastArgs,
        bool realign, const char *realignArgs,
        char *tempFile1);

void stCaf_sortCigarsByScoreInDescendingOrder(stList *cigars);

void stCaf_sortCigarsFileByScoreInDescendingOrder(char *cigarsFile, char *sortedFile);

#endif /* ST_LASTALIGNMENT_H_ */
