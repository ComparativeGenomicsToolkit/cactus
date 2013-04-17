/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "endAligner.h"
#include "multipleAligner.h"
#include "adjacencySequences.h"
#include "pairwiseAligner.h"

AlignedPair *alignedPair_construct(int64_t subsequenceIdentifier1, int64_t position1, bool strand1,
        int64_t subsequenceIdentifier2, int64_t position2, bool strand2, int64_t score) {
    AlignedPair *alignedPair = st_malloc(sizeof(AlignedPair));
    alignedPair->subsequenceIdentifier = subsequenceIdentifier1;
    alignedPair->position = position1;
    alignedPair->strand = strand1;

    alignedPair->reverse = st_malloc(sizeof(AlignedPair));
    alignedPair->reverse->reverse = alignedPair;

    alignedPair->reverse->subsequenceIdentifier = subsequenceIdentifier2;
    alignedPair->reverse->position = position2;
    alignedPair->reverse->strand = strand2;

    alignedPair->score = score;
    alignedPair->reverse->score = score;

    return alignedPair;
}

void alignedPair_destruct(AlignedPair *alignedPair) {
    free(alignedPair); //We assume the reverse will be free independently.
}

static int alignedPair_cmpFnP(const AlignedPair *alignedPair1, const AlignedPair *alignedPair2) {
    int i = cactusMisc_nameCompare(alignedPair1->subsequenceIdentifier, alignedPair2->subsequenceIdentifier);
    if(i == 0) {
        i = alignedPair1->position > alignedPair2->position ? 1 : (alignedPair1->position < alignedPair2->position ? -1 : 0);
        if(i == 0) {
            i = alignedPair1->strand == alignedPair2->strand ? 0 : (alignedPair1->strand ? 1 : -1);
        }
    }
    return i;
}

int alignedPair_cmpFn(const AlignedPair *alignedPair1, const AlignedPair *alignedPair2) {
    int i = alignedPair_cmpFnP(alignedPair1, alignedPair2);
    if(i == 0) {
        i = alignedPair_cmpFnP(alignedPair1->reverse, alignedPair2->reverse);
    }
    return i;
}

stSortedSet *makeEndAlignment(End *end, int64_t spanningTrees, int64_t maxSequenceLength,
        int64_t maximumNumberOfSequencesBeforeSwitchingToFast, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    //Make an alignment of the sequences in the ends

    //Get the adjacency sequences to be aligned.
    Cap *cap;
    End_InstanceIterator *it = end_getInstanceIterator(end);
    stList *sequences = stList_construct3(0, (void (*)(void *))adjacencySequence_destruct);
    stList *strings = stList_construct();
    while((cap = end_getNext(it)) != NULL) {
        if(cap_getSide(cap)) {
            cap = cap_getReverse(cap);
        }
        AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap, maxSequenceLength);
        stList_append(sequences, adjacencySequence);
        stList_append(strings, adjacencySequence->string);
    }
    end_destructInstanceIterator(it);

    //Convert the alignment pairs to an alignment of the caps..
    if(stList_length(strings) > maximumNumberOfSequencesBeforeSwitchingToFast && spanningTrees > 1) {
        spanningTrees = 1;
    }
    stList *alignment = makeAlignment(strings, spanningTrees, 100000000, gapGamma, pairwiseAlignmentBandingParameters);
    stSortedSet *sortedAlignment =
            stSortedSet_construct3((int (*)(const void *, const void *))alignedPair_cmpFn,
            (void (*)(void *))alignedPair_destruct);

    while(stList_length(alignment) > 0) {
        stIntTuple *alignedPair = stList_pop(alignment);
        assert(stIntTuple_length(alignedPair) == 5);
        AdjacencySequence *i = stList_get(sequences, stIntTuple_get(alignedPair, 1));
        AdjacencySequence *j = stList_get(sequences, stIntTuple_get(alignedPair, 3));
        assert(i != j);
        int64_t offset1 = stIntTuple_get(alignedPair, 2);
        int64_t offset2 = stIntTuple_get(alignedPair, 4);
        AlignedPair *alignedPair2 = alignedPair_construct(
                i->subsequenceIdentifier, i->start + (i->strand ? offset1 : -offset1), i->strand,
                j->subsequenceIdentifier, j->start + (j->strand ? offset2 : -offset2), j->strand,
                stIntTuple_get(alignedPair, 0));
        assert(stSortedSet_search(sortedAlignment, alignedPair2) == NULL);
        assert(stSortedSet_search(sortedAlignment, alignedPair2->reverse) == NULL);
        stSortedSet_insert(sortedAlignment, alignedPair2);
        stSortedSet_insert(sortedAlignment, alignedPair2->reverse);
        stIntTuple_destruct(alignedPair);
    }

    //Cleanup
    stList_destruct(strings);
    stList_destruct(sequences);
    stList_destruct(alignment);

    return sortedAlignment;
}

void writeEndAlignmentToDisk(End *end, stSortedSet *endAlignment, FILE *fileHandle) {
    fprintf(fileHandle, "%s\n", cactusMisc_nameToStringStatic(end_getName(end)));
    stSortedSetIterator *it = stSortedSet_getIterator(endAlignment);
    AlignedPair *aP;
    while((aP = stSortedSet_getNext(it)) != NULL) {
        fprintf(fileHandle, "%" PRIi64 " %" PRIi64 " %i ", aP->subsequenceIdentifier, aP->position, aP->strand);
        aP = aP->reverse;
        fprintf(fileHandle, "%" PRIi64 " %" PRIi64 " %i %" PRIi64 "\n", aP->subsequenceIdentifier, aP->position, aP->strand, aP->score);
    }
    stSortedSet_destructIterator(it);
}

stSortedSet *loadEndAlignmentFromDisk(Flower *flower, FILE *fileHandle, End **end) {
    stSortedSet *endAlignment =
                stSortedSet_construct3((int (*)(const void *, const void *))alignedPair_cmpFn,
                (void (*)(void *))alignedPair_destruct);
    char *line = stFile_getLineFromFile(fileHandle);
    *end = flower_getEnd(flower, cactusMisc_stringToName(line));
    if(end == NULL) {
        st_errAbort("Failed to get %s as an end from the given flower\n", line);
    }
    free(line);
    while((line = stFile_getLineFromFile(fileHandle)) != NULL) {
        int64_t sI1, sI2;
        int64_t p1, st1, p2, st2, score;
        int64_t i = sscanf(line, "%" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 "", &sI1, &p1, &st1, &sI2, &p2, &st2, &score);
        (void)i;
        assert(i == 7);
        stSortedSet_insert(endAlignment, alignedPair_construct(sI1, p1, st1, sI2, p2, st2, score));
        free(line);
    }
    return endAlignment;
}
