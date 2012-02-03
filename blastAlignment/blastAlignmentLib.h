/*
 * blastAlignmentLib.h
 *
 *  Created on: 30 Jan 2012
 *      Author: benedictpaten
 */

#ifndef BLASTALIGNMENTLIB_H_
#define BLASTALIGNMENTLIB_H_

#include "bioioC.h"
#include "cactus.h"
#include "sonLib.h"
#include "pairwiseAlignment.h"

int32_t writeFlowerSequencesInFile(Flower *flower, const char *tempFile1, int32_t minimumSequenceLength);

void convertCoordinatesOfPairwiseAlignment(struct PairwiseAlignment *pairwiseAlignment);

#endif /* BLASTALIGNMENTLIB_H_ */
