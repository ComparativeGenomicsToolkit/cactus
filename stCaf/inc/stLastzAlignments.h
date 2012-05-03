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

stList *stCaf_selfAlignFlower(Flower *flower, int32_t minimumSequenceLength, const char *lastzArgs, char *tempFile1);

#endif /* ST_LASTZALIGNMENT_H_ */
