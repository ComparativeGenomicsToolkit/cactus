/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "pairwiseAlignment.h"
#include "splitAlignmentOverlaps.h"

int main(int argc, char *argv[]) {
    /*
     * Each alignment has a unique first sequence interval, defined by where it starts and ends on the
     * first sequence.
     * Two alignments partially overlap if their first sequence intervals overlap but are not the same.
     * This program breaks up alignments in the input file so that there are no partial overlaps between
     * alignments, outputting the non-partially-overlapping alignments to the output file.
     */

    FILE *fileHandleIn;
    FILE *fileHandleOut;

    uint64_t maxDepth;

    if (!(argc == 3 || argc == 5)) {
        st_errAbort("Usage: %s logLevel maxDepth [inFile] [outFile]", argv[0]);
    }

    st_setLogLevelFromString(argv[1]);
    if (sscanf(argv[2], "%" PRIu64, &maxDepth) != 1) {
        st_errAbort("Couldn't understand maxDepth");
    }

    if(argc == 3) {
        fileHandleIn = stdin;
        fileHandleOut = stdout;
    }
    else {
        assert(argc == 4);
        fileHandleIn = fopen(argv[3], "r");
        fileHandleOut = fopen(argv[4], "w");
    }

    splitAlignmentOverlaps_loop(fileHandleIn, fileHandleOut, maxDepth);

    if(argc == 4) {
    	fclose(fileHandleIn);
    	fclose(fileHandleOut);
    }

    return 0;
}
