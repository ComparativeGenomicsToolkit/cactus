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
        int64_t score, int64_t scoreWithoutStubs, int64_t columnId) {
    AlignedPair *alignedPair = st_malloc(sizeof(AlignedPair));
    alignedPair->subsequenceIdentifier = subsequenceIdentifier1;
    alignedPair->position = position1;
    alignedPair->strand = strand1;
    alignedPair->score = score;
    alignedPair->scoreWithoutStubs = scoreWithoutStubs;
    alignedPair->columnId = columnId;

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
        return alignedPair1->columnId > alignedPair2->columnId ? 1 : (alignedPair1->columnId < alignedPair2->columnId ? -1 : 0);
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
    stList *seqFrags = stList_construct3(0, (void (*)(void *))seqFrag_destruct);
    while((cap = end_getNext(it)) != NULL) {
        if(cap_getSide(cap)) {
            cap = cap_getReverse(cap);
        }
        AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap, maxSequenceLength);
        stList_append(sequences, adjacencySequence);
        assert(cap_getAdjacency(cap) != NULL);
        stList_append(seqFrags, seqFrag_construct(adjacencySequence->string, 0, end_getName(cap_getEnd(cap_getAdjacency(cap)))));
    }
    end_destructInstanceIterator(it);

    //Convert the alignment pairs to an alignment of the caps..
    MultipleAlignment *mA = makeAlignment(seqFrags, spanningTrees, 100000000, maximumNumberOfSequencesBeforeSwitchingToFast, gapGamma, pairwiseAlignmentBandingParameters);

    //Build the sequence position to columns scores.
    stHash *sequencePositionsToColumnsHash = getSequencePositionsToColumnsHash(mA->columns);
    stHash *sequencePositionsToColumnScoresHash = getSequencePositionsToColumnScoresHash(sequencePositionsToColumnsHash, seqFrags, mA, NULL, NULL);
    stHash *sequencePositionsToColumnScoresWithoutStubsHash = getSequencePositionsToColumnScoresHash(sequencePositionsToColumnsHash, seqFrags, mA, NULL, NULL);

    //Make an index for each column
    stHash *columnsToIndices = stHash_construct2(NULL, (void (*)(void *))stIntTuple_destruct);
    stSetIterator *columnIt = stSet_getIterator(mA->columns);
    Column *c;
    int64_t columnIndex = 0;
    while((c = stSet_getNext(columnIt)) != NULL) {
        stHash_insert(columnsToIndices, c, stIntTuple_construct1(columnIndex++));
    }
    stSet_destructIterator(columnIt);

    stSortedSet *sortedAlignment =
            stSortedSet_construct3((int (*)(const void *, const void *))alignedPair_cmpFn,
            (void (*)(void *))alignedPair_destruct);

    stHashIterator *columnIt2 = stHash_getIterator(sequencePositionsToColumnsHash);
    while((c = stHash_getNext(columnIt2)) != NULL) {
        int64_t *score = stHash_search(sequencePositionsToColumnScoresHash, c);
        int64_t *scoreWithoutStubs = stHash_search(sequencePositionsToColumnScoresWithoutStubsHash, c);
        Column *c2 = stHash_search(sequencePositionsToColumnsHash, c);
        assert(c2 != NULL);
        AdjacencySequence *i = stList_get(sequences, c->seqIndex);
        AlignedPair *alignedPair = alignedPair_construct(i->subsequenceIdentifier,
                i->start + (i->strand ? c->position : -c->position), i->strand, *score,
                /*score of alignment only to non-stubs*/ *scoreWithoutStubs,
                stIntTuple_get(stHash_search(columnsToIndices, c2), 0));
        assert(stSortedSet_search(sortedAlignment, alignedPair) == NULL);
        stSortedSet_insert(sortedAlignment, alignedPair);
    }
    stHash_destructIterator(columnIt2);

    //Cleanup
    stHash_destruct(sequencePositionsToColumnsHash);
    stHash_destruct(sequencePositionsToColumnScoresHash);
    stHash_destruct(sequencePositionsToColumnScoresWithoutStubsHash);
    stList_destruct(seqFrags);
    stList_destruct(sequences);
    stHash_destruct(columnsToIndices);
    multipleAlignment_destruct(mA);

    return sortedAlignment;
}

void writeEndAlignmentToDisk(End *end, stSortedSet *endAlignment, FILE *fileHandle) {
    fprintf(fileHandle, "%s %" PRIi64 "\n", cactusMisc_nameToStringStatic(end_getName(end)), stSortedSet_size(endAlignment));
    stSortedSetIterator *it = stSortedSet_getIterator(endAlignment);
    AlignedPair *aP;
    while((aP = stSortedSet_getNext(it)) != NULL) {
        fprintf(fileHandle, "%" PRIi64 " %" PRIi64 " %i %" PRIi64 " %" PRIi64 " %" PRIi64 "\n", aP->subsequenceIdentifier, aP->position, aP->strand, aP->score, aP->scoreWithoutStubs, aP->columnId);
    }
    stSortedSet_destructIterator(it);
}

stSortedSet *loadEndAlignmentFromDisk(Flower *flower, FILE *fileHandle, End **end) {
    stSortedSet *endAlignment =
                stSortedSet_construct3((int (*)(const void *, const void *))alignedPair_cmpFn,
                (void (*)(void *))alignedPair_destruct);
    char *line = stFile_getLineFromFile(fileHandle);
    if(line == NULL) {
        *end = NULL;
        return NULL;
    }
    Name flowerName;
    int64_t lineNumber;
    int64_t i = sscanf(line, "%" PRIi64 " %" PRIi64 "", &flowerName, &lineNumber);
    if(i != 2 || lineNumber < 0) {
        st_errAbort("We encountered a mis-specified name in loading the first line of an end alignment from the disk: '%s'\n", line);
    }
    free(line);
    *end = flower_getEnd(flower, flowerName);
    if(*end == NULL) {
        st_errAbort("We encountered an end name that is not in the database: '%s'\n", line);
    }
    for(int64_t i=0; i<lineNumber; i++) {
        line = stFile_getLineFromFile(fileHandle);
        if(line == NULL) {
            st_errAbort("Got a null line when parsing an end alignment\n");
        }
        int64_t sI1, p1, st1, columnId, score, scoreWithoutStubs;
        int64_t i = sscanf(line, "%" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 "", &sI1, &p1, &st1, &columnId, &score, &scoreWithoutStubs);
        (void)i;
        if(i != 6) {
            st_errAbort("We encountered a mis-specified name in loading an end alignment from the disk: '%s'\n", line);
        }
        stSortedSet_insert(endAlignment, alignedPair_construct(sI1, p1, st1, score, scoreWithoutStubs, columnId));
        free(line);
    }
    return endAlignment;
}

