/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>

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
    struct option opts[] = { {"onlyContig1", no_argument, NULL, '1'},
                             {"onlyContig2", no_argument, NULL, '2'},
                             {0, 0, 0, 0} };
    int convertContig1 = TRUE, convertContig2 = TRUE, flag;
    while((flag = getopt_long(argc, argv, "", opts, NULL)) != -1) {
        switch(flag) {
        case '1':
            convertContig2 = FALSE;
            break;
        case '2':
            convertContig1 = FALSE;
            break;
        }
    }
    if(!(convertContig1 || convertContig2)) {
        fprintf(stderr, "--onlyContig1 and --onlyContig2 options are "
                "mutually exclusive\n");
        return 1;
    }
    assert(argc == optind + 3);
    FILE *fileHandleIn = fopen(argv[optind], "r");
    FILE *fileHandleOut = fopen(argv[optind + 1], "w");
    int64_t roundsOfConversion;
    int64_t i = sscanf(argv[optind + 2], "%" PRIi64 "", &roundsOfConversion);
    (void)i;
    assert(i == 1);
    assert(roundsOfConversion >= 1);
    struct PairwiseAlignment *pairwiseAlignment;
    while ((pairwiseAlignment = cigarRead(fileHandleIn)) != NULL) {
        //Correct coordinates
        for(int64_t j=0; j<roundsOfConversion; j++) {
            convertCoordinatesOfPairwiseAlignment(pairwiseAlignment,
                                                  convertContig1,
                                                  convertContig2);
        }
        cigarWrite(fileHandleOut, pairwiseAlignment, 0);
        destructPairwiseAlignment(pairwiseAlignment);
    }
    fclose(fileHandleIn);
    fclose(fileHandleOut);
    return 0;
}
