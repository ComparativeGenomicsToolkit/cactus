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

int64_t writeFlowerSequencesInFile(Flower *flower, const char *tempFile1, int64_t minimumSequenceLength);

int64_t writeFlowerSequences(Flower *flower, void(*processSequence)(const char *, const char *, int64_t), int64_t minimumSequenceLength);

void convertCoordinatesOfPairwiseAlignment(struct PairwiseAlignment *pairwiseAlignment, int convertContig1, int convertContig2);

void setupToChunkSequences(int64_t chunkSize2, int64_t overlapSize2, const char *chunksDir2);

void processSequenceToChunk(const char *fastaHeader, const char *sequence, int64_t length);

void finishChunkingSequences();

#endif /* BLASTALIGNMENTLIB_H_ */
