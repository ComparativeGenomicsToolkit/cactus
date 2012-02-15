/*
 * blastAlignmentLib.c
 *
 *  Created on: 30 Jan 2012
 *      Author: benedictpaten
 */

#include "bioioC.h"
#include "cactus.h"
#include "sonLib.h"
#include "pairwiseAlignment.h"

int32_t writeFlowerSequencesInFile(Flower *flower, const char *tempFile1, int32_t minimumSequenceLength) {
    FILE *fileHandle = NULL;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    int32_t sequencesWritten = 0;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(instanceIterator)) != NULL) {
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            Cap *cap2 = cap_getAdjacency(cap);
            assert(cap2 != NULL);
            assert(cap_getStrand(cap2));

            if (!cap_getSide(cap)) {
                assert(cap_getSide(cap2));
                int32_t length = cap_getCoordinate(cap2) - cap_getCoordinate(cap) - 1;
                assert(length >= 0);
                if (length >= minimumSequenceLength) {
                    Sequence *sequence = cap_getSequence(cap);
                    assert(sequence != NULL);
                    if(fileHandle == NULL) {
                        fileHandle = fopen(tempFile1, "w");
                    }
                    char *string = sequence_getString(sequence, cap_getCoordinate(cap) + 1,
                            cap_getCoordinate(cap2) - cap_getCoordinate(cap) - 1, 1);
                    fprintf(fileHandle, ">%s|1|%i\n%s\n", cactusMisc_nameToStringStatic(sequence_getName(sequence)),
                            cap_getCoordinate(cap) + 1, string);
                    free(string);
                    sequencesWritten++;
                }
            }
        }
        end_destructInstanceIterator(instanceIterator);
    }
    flower_destructEndIterator(endIterator);
    if(fileHandle != NULL) {
        fclose(fileHandle);
    }
    return sequencesWritten;
}

static void convertCoordinatesP(char *header, char **contig1, int32_t *start, int32_t *strand) {
    struct List *attributes = fastaDecodeHeader(header);
    //uglyf(" strings: %s\n", stringsJoin(" $ ", attributes->list, attributes->length));
    assert(attributes->length >= 3);
    assert(sscanf((const char *) attributes->list[attributes->length - 1], "%i", start) == 1);
    assert(sscanf((const char *) attributes->list[attributes->length - 2], "%i", strand) == 1);
    assert(*strand == 0 || *strand == 1);
    free(attributes->list[--attributes->length]);
    free(attributes->list[--attributes->length]);
    *contig1 = fastaEncodeHeader(attributes);
    destructList(attributes);
}

void convertCoordinatesOfPairwiseAlignment(struct PairwiseAlignment *pairwiseAlignment) {
    checkPairwiseAlignment(pairwiseAlignment);

    char *contig;
    int32_t start, strand;
    convertCoordinatesP(pairwiseAlignment->contig1, &contig, &start, &strand);
    free(pairwiseAlignment->contig1);
    pairwiseAlignment->contig1 = contig;
    pairwiseAlignment->strand1 = strand ? pairwiseAlignment->strand1 : !pairwiseAlignment->strand1;
    pairwiseAlignment->start1 = strand ? pairwiseAlignment->start1 + start : start - pairwiseAlignment->start1 + 1;
    pairwiseAlignment->end1 = strand ? pairwiseAlignment->end1 + start : start - pairwiseAlignment->end1 + 1;

    convertCoordinatesP(pairwiseAlignment->contig2, &contig, &start, &strand);
    free(pairwiseAlignment->contig2);
    pairwiseAlignment->contig2 = contig;
    pairwiseAlignment->strand2 = strand ? pairwiseAlignment->strand2 : !pairwiseAlignment->strand2;
    pairwiseAlignment->start2 = strand ? pairwiseAlignment->start2 + start : start - pairwiseAlignment->start2 + 1;
    pairwiseAlignment->end2 = strand ? pairwiseAlignment->end2 + start : start - pairwiseAlignment->end2 + 1;

    checkPairwiseAlignment(pairwiseAlignment);
}
