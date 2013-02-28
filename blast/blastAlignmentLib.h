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

int32_t writeFlowerSequences(Flower *flower, void(*processSequence)(const char *, const char *, int32_t), int32_t minimumSequenceLength);

void convertCoordinatesOfPairwiseAlignment(struct PairwiseAlignment *pairwiseAlignment);

void convertCoordinates(char *cigarFile, FILE *fileHandle);

void setupToChunkSequences(int64_t chunkSize2, int64_t overlapSize2, const char *chunksDir2);

void processSequenceToChunk(const char *fastaHeader, const char *sequence, int32_t length);

void finishChunkingSequences();

#endif /* BLASTALIGNMENTLIB_H_ */
