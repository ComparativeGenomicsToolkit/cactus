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
    assert(argc == 3);
    FILE *fileHandleIn = fopen(argv[1], "r");
    FILE *fileHandleOut = fopen(argv[2], "w");
    while ((pairwiseAlignment = cigarRead(fileHandleIn)) != NULL) {
        //Correct coordinates
        convertCoordinatesOfPairwiseAlignment(pairwiseAlignment);
        cigarWrite(fileHandleOut, pairwiseAlignment, 0);
        destructPairwiseAlignment(pairwiseAlignment);
    }
    fclose(fileHandleIn);
    fclose(fileHandleOut);
    return 0;
}
