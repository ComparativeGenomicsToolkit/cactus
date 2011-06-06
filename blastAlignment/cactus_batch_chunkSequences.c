/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include "bioioC.h"
#include "commonC.h"

/*
 * Assistance routine for pecan2_batch.py.
 * Routine reads in fasta file progressively, writes chunks of sequence.
 */

int32_t chunkRemaining;
FILE *fileHandle = NULL;
FILE *fileHandle2 = NULL;
char *tempSeqFile = NULL;
int32_t chunkSize;
int32_t overlapSize;
int32_t useCompression;

int32_t fn(char *fastaHeader, int32_t start, char *sequence, int32_t seqLength, int32_t length) {
    if(fileHandle == NULL) {
        tempSeqFile = getTempFile();
        fileHandle = fopen(tempSeqFile, "w");
    }

    int32_t i = 0;
    fastaHeader = stString_copy(fastaHeader);
    while (fastaHeader[i] != '\0') {
        if (fastaHeader[i] == ' ' || fastaHeader[i] == '\t') {
            fastaHeader[i] = '\0';
            break;
        }
        i++;
    }
    fprintf(fileHandle, ">%s|1|%i\n", fastaHeader, start);
    free(fastaHeader);
    assert(length <= chunkSize);
    assert(start >= 0);
    if (start + length > seqLength) {
        length = seqLength - start;
    }
    assert(length > 0);
    char c = sequence[start + length];
    sequence[start + length] = '\0';
    fprintf(fileHandle, "%s\n", &sequence[start]);
    sequence[start + length] = c;

    return length;
}

void fn4() {
    if(fileHandle != NULL) {
        fclose(fileHandle);
        fprintf(fileHandle2, "%s\n", tempSeqFile);

        //Do compression if asked for
        if (useCompression) {
            int32_t i = st_system("bzip2 %s", tempSeqFile);
            if (i != 0) {
                st_logInfo("Failed to compress file\n");
                stThrowNew("CACTUS_BATCH_CHUNK_EXCEPTION", "Failed to compress file");
            }
            i = st_system("mv %s.bz2 %s", tempSeqFile, tempSeqFile);
            if (i != 0) {
                stThrowNew("CACTUS_BATCH_CHUNK_EXCEPTION", "Failed to move file");
            }
        }

        //free(tempSeqFile);
        tempSeqFile = NULL;
        fileHandle = NULL;
    }
}

int32_t fn2(int32_t i, int32_t seqLength) {
    //Update remaining portion of the chunk.
    assert(seqLength >= 0);
    i -= seqLength;
    if (i <= 0) {
        fn4();
        return chunkSize;
    }
    return i;
}

void fn3(const char *fastaHeader, const char *sequence, int32_t length) {
    int32_t j, k, l;

    if (length > 0) {
        j = fn((char *) fastaHeader, 0, (char *) sequence, length, chunkRemaining);
        chunkRemaining = fn2(chunkRemaining, j);
        while (length - j > 0) {
            //Make the non overlap file
            k = fn((char *) fastaHeader, j, (char *) sequence, length, chunkRemaining);
            chunkRemaining = fn2(chunkRemaining, k);

            //Make the overlap file
            l = j - overlapSize / 2;
            if(l < 0) {
				l = 0;
			}
            chunkRemaining = fn2(chunkRemaining, fn((char *) fastaHeader, l, (char *) sequence, length, overlapSize));
            j += k;
        }
    }
}

int main(int argc, char *argv[]) {
    assert(argc == 7);
    fileHandle2 = fopen(argv[1], "w");
    int32_t i = sscanf(argv[2], "%i", &chunkSize);
    assert(chunkSize > 0);
    assert(i == 1);
    i = sscanf(argv[3], "%i", &overlapSize);
    assert(i == 1);
    initialiseTempFileTree(argv[4], 100, 4);
    i = sscanf(argv[5], "%i", &useCompression);
    assert(i == 1);
    chunkRemaining = chunkSize;

    FILE *fileHandle3 = fopen(argv[6], "r");
    int size = 100;
    char *cA = st_calloc(size + 1, sizeof(char));
    do {
        i = benLine(&cA, &size, fileHandle3);
        if (strlen(cA) > 0) {
            FILE *fileHandle4 = fopen(cA, "r");
            fastaReadToFunction(fileHandle4, fn3);
            fclose(fileHandle4);
        }
    } while (i != -1);
    fclose(fileHandle3);
    free(cA);
    fn4();
    fclose(fileHandle2);
    return 0;
}
