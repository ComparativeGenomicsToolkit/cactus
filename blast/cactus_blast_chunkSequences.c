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

#include "blastAlignmentLib.h"

int main(int argc, char *argv[]) {
    //log-string, chunkSize, overlapSize, dirToPutChunksIn, seqFilesX n
    assert(argc >= 5);
    st_setLogLevelFromString(argv[1]);
    int64_t chunkSize, chunkOverlapSize;
    int64_t i = sscanf(argv[2], "%" PRIi64 "", &chunkSize);
    assert(i == 1);
    i = sscanf(argv[3], "%" PRIi64 "", &chunkOverlapSize);
    assert(i == 1);
    setupToChunkSequences(chunkSize, chunkOverlapSize, argv[4]);
    for (int64_t i = 5; i < argc; i++) {
        FILE *fileHandle2 = fopen(argv[i], "r");
        fastaReadToFunction(fileHandle2, processSequenceToChunk);
        fclose(fileHandle2);
    }
    finishChunkingSequences();
    return 0;
}
