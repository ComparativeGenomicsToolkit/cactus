/*
 * blastAlignmentLib.c
 *
 *  Created on: 30 Jan 2012
 *      Author: benedictpaten
 */

#include "bioioC.h"
#include "cactus.h"
#include "sonLib.h"
#include "pairwiseAlignment.h"

/*
 * Converting coordinates of pairwise alignments
 */

static void convertCoordinatesP(char **contig, int64_t *start, int64_t *end) {
    stList *attributes = fastaDecodeHeader(*contig);
    //Decode attributes
    int64_t startP;
    int64_t i = sscanf((const char *) stList_peek(attributes), "%" PRIi64 "", &startP);
    (void) i;
    assert(i == 1);
    free(stList_pop(attributes));
    //Now relabel attributes
    free(*contig);
    *contig = fastaEncodeHeader(attributes);
    stList_destruct(attributes);
    *start = *start + startP;
    *end = *end + startP;
}

void convertCoordinatesOfPairwiseAlignment(struct PairwiseAlignment *pairwiseAlignment,
                                           int convertContig1,
                                           int convertContig2) {
    checkPairwiseAlignment(pairwiseAlignment);
    if(convertContig1) {
        convertCoordinatesP(&pairwiseAlignment->contig1, &pairwiseAlignment->start1, &pairwiseAlignment->end1);
    }
    if(convertContig2) {
        convertCoordinatesP(&pairwiseAlignment->contig2, &pairwiseAlignment->start2, &pairwiseAlignment->end2);
    }
    checkPairwiseAlignment(pairwiseAlignment);
}

/*
 * Routine reads in chunk up a set of sequences into overlapping sequence files.
 */

static int64_t chunkRemaining;
static FILE *chunkFileHandle = NULL;
static const char *chunksDir = NULL;
static int64_t chunkNo = 0;
static char *tempChunkFile = NULL;
static int64_t chunkSize;
static int64_t chunkOverlapSize;

void finishChunkingSequences() {
    if (chunkFileHandle != NULL) {
        fclose(chunkFileHandle);
        fprintf(stdout, "%s\n", tempChunkFile);
        free(tempChunkFile);
        tempChunkFile = NULL;
        chunkFileHandle = NULL;
    }
}

static void updateChunkRemaining(int64_t seqLength) {
    //Update remaining portion of the chunk.
    assert(seqLength >= 0);
    chunkRemaining -= seqLength;
    if (chunkRemaining <= 0) {
        finishChunkingSequences();
        chunkRemaining = chunkSize;
    }
}

static int64_t processSubsequenceChunk(char *fastaHeader, int64_t start, char *sequence, int64_t seqLength, int64_t lengthOfChunkRemaining) {
    if (chunkFileHandle == NULL) {
        tempChunkFile = stString_print("%s/%" PRIi64 "", chunksDir, chunkNo++);
        chunkFileHandle = fopen(tempChunkFile, "w");
    }

    int64_t i = 0;
    fastaHeader = stString_copy(fastaHeader);
    while (fastaHeader[i] != '\0') {
        if (fastaHeader[i] == ' ' || fastaHeader[i] == '\t') {
            fastaHeader[i] = '\0';
            break;
        }
        i++;
    }
    char *chunkHeader = stString_print("%s|%" PRIi64 "\n", fastaHeader, start);
    free(fastaHeader);
    assert(lengthOfChunkRemaining <= chunkSize);
    assert(start >= 0);
    int64_t lengthOfSubsequence = lengthOfChunkRemaining;
    if (start + lengthOfChunkRemaining > seqLength) {
        lengthOfSubsequence = seqLength - start;
    }
    assert(lengthOfSubsequence > 0);
    char c = sequence[start + lengthOfSubsequence];
    sequence[start + lengthOfSubsequence] = '\0';
    fastaWrite(&sequence[start], chunkHeader, chunkFileHandle);
    //fprintf(chunkFileHandle, "%s\n", &sequence[start]);
    free(chunkHeader);
    sequence[start + lengthOfSubsequence] = c;

    updateChunkRemaining(lengthOfSubsequence);
    return lengthOfSubsequence;
}

void processSequenceToChunk(void* dest, const char *fastaHeader, const char *sequence, int64_t sequenceLength) {
    if (sequenceLength > 0) {
        int64_t lengthOfSubsequence = processSubsequenceChunk((char *) fastaHeader, 0, (char *) sequence, sequenceLength, chunkRemaining);
        while (sequenceLength - lengthOfSubsequence > 0) {
            //Make the non overlap file
            int64_t lengthOfFollowingSubsequence = processSubsequenceChunk((char *) fastaHeader, lengthOfSubsequence, (char *) sequence, sequenceLength, chunkRemaining);

            //Make the overlap file
            if(chunkOverlapSize > 0) {
                int64_t i = lengthOfSubsequence - chunkOverlapSize / 2;
                if (i < 0) {
                    i = 0;
                }
                processSubsequenceChunk((char *) fastaHeader, i, (char *) sequence, sequenceLength, chunkOverlapSize);
            }
            lengthOfSubsequence += lengthOfFollowingSubsequence;
        }
    }
}

void setupToChunkSequences(int64_t chunkSize2, int64_t overlapSize2, const char *chunksDir2) {
    chunkSize = chunkSize2;
    assert(chunkSize > 0);
    chunkOverlapSize = overlapSize2;
    assert(chunkOverlapSize >= 0);
    chunksDir = chunksDir2;
    chunkNo = 0;
    chunkRemaining = chunkSize;
    chunkFileHandle = NULL;
}

/*
 * Get the flowers in a file.
 */

int64_t writeFlowerSequences(Flower *flower, void(*processSequence)(void *destination, const char *name, const char *seq, int64_t length), void *destination, int64_t minimumSequenceLength) {
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    int64_t sequencesWritten = 0;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(instanceIterator)) != NULL) {
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            Cap *cap2 = cap_getAdjacency(cap);
            assert(cap2 != NULL);
            assert(cap_getStrand(cap2));

            if (!cap_getSide(cap)) {
                assert(cap_getSide(cap2));
                int64_t length = cap_getCoordinate(cap2) - cap_getCoordinate(cap) - 1;
                assert(length >= 0);
                if (length >= minimumSequenceLength) {
                    Sequence *sequence = cap_getSequence(cap);
                    assert(sequence != NULL);
                    char *string = sequence_getString(sequence, cap_getCoordinate(cap) + 1, length, 1);
                    char *header = stString_print("%" PRIi64 "|%" PRIi64 "", cap_getName(cap), cap_getCoordinate(cap) + 1);
                    processSequence(destination, header, string, strlen(string));
                    free(string);
                    free(header);
                    sequencesWritten++;
                }
            }
        }
        end_destructInstanceIterator(instanceIterator);
    }
    flower_destructEndIterator(endIterator);
    return sequencesWritten;
}

/* hack to remain compatible until blastAlignmentLib can be made thread safe */
int64_t writeFlowerSequencesUnsafe(Flower *flower, void(*processSequence)(void *destination, const char *name, const char *seq, int64_t length), int64_t minimumSequenceLength) {
    return writeFlowerSequences(flower, processSequence, NULL, minimumSequenceLength);
}
    


// Delay opening so it doesn't write file if no flowers.  Not sure if this behavior is
// relied on, so keeping old behavior when making code thread safe.
typedef struct {
    FILE *fh;
    const char *fname;
} SeqFileInfo;

static void writeSequenceInFile(void *destination, const char *fastaHeader, const char *sequence, int64_t length) {
    SeqFileInfo *seqFileInfo = destination;
    if (seqFileInfo->fh == NULL) {
        seqFileInfo->fh = fopen(seqFileInfo->fname, "w");
    }
    fastaWrite((char *)sequence, (char *)fastaHeader, seqFileInfo->fh);
}

int64_t writeFlowerSequencesInFile(Flower *flower, const char *tempFile, int64_t minimumSequenceLength) {
    SeqFileInfo seqFileInfo = {NULL, tempFile};
    int64_t sequencesWritten = writeFlowerSequences(flower, writeSequenceInFile, &seqFileInfo, minimumSequenceLength);
    if (seqFileInfo.fh != NULL) {
        fclose(seqFileInfo.fh);
    }
    return sequencesWritten;
}
