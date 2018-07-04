/*
 * lastzAlignment.h
 *
 *  Created on: 30 Jan 2012
 *      Author: benedictpaten
 */

#ifndef ST_LASTZALIGNMENT_H_
#define ST_LASTZALIGNMENT_H_

#include "cactus.h"
#include "sonLib.h"

stList *stCaf_selfAlignFlower(Flower *flower, int64_t minimumSequenceLength, const char *lastzArgs,
        bool realign, const char *realignArgs,
        char *tempFile1);

void stCaf_sortCigarsByScoreInDescendingOrder(stList *cigars);

void stCaf_sortCigarsFileByScoreInDescendingOrder(char *cigarsFile, char *sortedFile);

#endif /* ST_LASTZALIGNMENT_H_ */
