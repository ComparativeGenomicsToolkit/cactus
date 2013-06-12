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

stList *stCaf_selfAlignFlower(Flower *flower, int64_t minimumSequenceLength, const char *lastzArgs, char *tempFile1) {
    /*
     * Get the sequences.
     */
    stList *cigars = stList_construct3(0, (void(*)(void *)) destructPairwiseAlignment);
    //char *tempFile1 = getTempFile();
    if (writeFlowerSequencesInFile(flower, tempFile1, minimumSequenceLength) > 0) {
        /*
         * Run lastz.
         */
        char *command = stString_print(
                "lastz --format=cigar %s %s[multiple][nameparse=darkspace] %s[nameparse=darkspace] --notrivial",
                lastzArgs, tempFile1, tempFile1);
        //char *command = stString_print(
        //        "lastz --format=cigar %s %s[multiple][nameparse=darkspace] --self",
        //        lastzArgs, tempFile1);
        FILE *fileHandle = popen(command, "r");
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
        int i = pclose(fileHandle);
        if(i != 0) {
            st_errAbort("Lastz failed: %s\n", command);
        }
        free(command);
    }
    //st_system("rm %s", tempFile1);

    return cigars;
}

static int compareByScore(struct PairwiseAlignment *pA, struct PairwiseAlignment *pA2) {
    return pA->score == pA2->score ? 0 : (pA->score > pA2->score ? -1 : 1);
}

void stCaf_sortCigarsByScoreInDescendingOrder(stList *cigars) {
    stList_sort(cigars, (int (*)(const void *, const void *))compareByScore);
#ifndef NDEBUG
        double score = INT64_MAX;
        for(int64_t i=0; i<stList_length(cigars); i++) {
            struct PairwiseAlignment *pA = stList_get(cigars, i);
            assert(pA->score <= score);
            score = pA->score;
        }
#endif
}

void stCaf_sortCigarsFileByScoreInDescendingOrder(char *cigarsFile, char *sortedFile) {
    int64_t i = st_system("sort -k10nr %s > %s", cigarsFile, sortedFile);
    if(i != 0) {
        st_errAbort("Encountered unix sort error when sorting cigar alignments in file: %s\n", cigarsFile);
    }
#ifndef NDEBUG
    double score = INT64_MAX;
    FILE *fileHandle = fopen(sortedFile, "r");
    struct PairwiseAlignment *pA;
    while ((pA = cigarRead(fileHandle)) != NULL) {
        st_uglyf("I got %f\n", score);
        assert(pA->score <= score);
        score = pA->score;
    }
    fclose(fileHandle);
#endif
}
