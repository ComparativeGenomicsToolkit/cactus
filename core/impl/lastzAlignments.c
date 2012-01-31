/*
 * lastzAlignments.c
 *
 *  Created on: 30 Jan 2012
 *      Author: benedictpaten
 */

#define _XOPEN_SOURCE 500

#include "bioioC.h"
#include "cactus.h"
#include "sonLib.h"
#include "pairwiseAlignment.h"
#include "blastAlignmentLib.h"

stList *selfAlignFlower(Flower *flower, int32_t minimumSequenceLength, const char *lastzArgs) {
    /*
     * Get the sequences.
     */
    stList *cigars = stList_construct3(0, (void(*)(void *)) destructPairwiseAlignment);
    char *tempFile1 = getTempFile();
    if (writeFlowerSequencesInFile(flower, tempFile1, minimumSequenceLength) != 0) {
        /*
         * Run lastz.
         */
        char *command = stString_print(
                "lastz --format=cigar %s %s[multiple][nameparse=darkspace] %s[nameparse=darkspace] --notrivial",
                lastzArgs, tempFile1, tempFile1);
        st_uglyf("I am running #%s#\n", command);
        FILE *fileHandle = popen(command, "r");
        free(command);
        if (fileHandle == NULL) {
            st_errAbort("Problems with lastz pipe");
        }

        /*
         * Process the cigars, modifying their coordinates.
         */
        //Read from stream
        struct PairwiseAlignment *pairwiseAlignment;
        while ((pairwiseAlignment = cigarRead(fileHandle)) != NULL) {
            convertCoordinatesOfPairwiseAlignment(pairwiseAlignment);
            stList_append(cigars, pairwiseAlignment);
        }
    }
    st_system("rm %s", tempFile1);

    return cigars;
}
