/*
 * lastzAlignment.h
 *
 *  Created on: 30 Jan 2012
 *      Author: benedictpaten
 */

#ifndef LASTZALIGNMENT_H_
#define LASTZALIGNMENT_H_

#include "cactus.h"
#include "sonLib.h"

stList *selfAlignFlower(Flower *flower, int32_t minimumSequenceLength, const char *lastzArgs);

#endif /* LASTZALIGNMENT_H_ */
