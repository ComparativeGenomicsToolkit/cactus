/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "commonC.h"
#include "hashTableC.h"
#include "bioioC.h"
#include "cactus.h"
#include "avl.h"
#include "pairwiseAlignment.h"
#include "blastAlignmentLib.h"

int main(int argc, char *argv[]) {
    /*
     * For each cigar in file, update the coordinates and write to the second file.
     */
    assert(argc == 4);
    FILE *fileHandleIn = fopen(argv[1], "r");
    FILE *fileHandleOut = fopen(argv[2], "w");
    int32_t roundsOfConversion;
    int32_t i = sscanf(argv[3], "%i", &roundsOfConversion);
    (void)i;
    assert(i == 1);
    assert(roundsOfConversion >= 1);
    struct PairwiseAlignment *pairwiseAlignment;
    while ((pairwiseAlignment = cigarRead(fileHandleIn)) != NULL) {
        //Correct coordinates
        for(int32_t j=0; j<roundsOfConversion; j++) {
            convertCoordinatesOfPairwiseAlignment(pairwiseAlignment);
        }
        cigarWrite(fileHandleOut, pairwiseAlignment, 0);
        destructPairwiseAlignment(pairwiseAlignment);
    }
    fclose(fileHandleIn);
    fclose(fileHandleOut);
    return 0;
}
