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

#include "blastAlignmentLib.c"

int main(int argc, char *argv[]) {
    assert(argc >= 5);
    st_setLogLevelFromString(argv[1]);
    int32_t chunkSize, chunkOverlapSize;
    int32_t i = sscanf(argv[2], "%i", &chunkSize);
    assert(i == 1);
    i = sscanf(argv[3], "%i", &chunkOverlapSize);
    assert(i == 1);
    setupToChunkSequences(chunkSize, chunkOverlapSize, argv[4]);
    for (int32_t i = 5; i < argc; i++) {
        FILE *fileHandle2 = fopen(argv[i], "r");
        fastaReadToFunction(fileHandle2, processSequenceToChunk);
        fclose(fileHandle2);
    }
    finishChunkingSequences();
    return 0;
}
