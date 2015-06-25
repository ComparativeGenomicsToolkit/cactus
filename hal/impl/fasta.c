#include <ctype.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "cactus.h"
#include "sonLib.h"
#include "bioioC.h"

static Name globalReferenceEventName;

static int compareSequences(Sequence *sequence, Sequence *sequence2) {
    Event *event = sequence_getEvent(sequence);
    Event *event2 = sequence_getEvent(sequence2);
    int i = cactusMisc_nameCompare(event_getName(event), event_getName(event2));
    if (i != 0) {
        return event_getName(event) == globalReferenceEventName ? -1 : (event_getName(event2) == globalReferenceEventName ? 1 : i);
    }
    i = cactusMisc_nameCompare(sequence_getName(sequence), sequence_getName(sequence2));
    return i;
}

static stList *getSequences(Flower *flower, Name referenceEventName) {
    globalReferenceEventName = referenceEventName;
    stList *sequences = stList_construct();
    Sequence *sequence;
    Flower_SequenceIterator *seqIt = flower_getSequenceIterator(flower);
    while ((sequence = flower_getNextSequence(seqIt)) != NULL) {
        stList_append(sequences, sequence);
    }
    flower_destructSequenceIterator(seqIt);
    stList_sort(sequences, (int (*)(const void *, const void *))compareSequences);
    return sequences;
}

void printFastaSequences(Flower *flower, FILE *fileHandle, Name referenceEventName) {
    stList *sequences = getSequences(flower, referenceEventName);
    for(int64_t i=0; i<stList_length(sequences); i++) {
        Sequence *sequence = stList_get(sequences, i);
        if(!metaSequence_isTrivialSequence(sequence_getMetaSequence(sequence))) {
            char *string = sequence_getString(sequence, sequence_getStart(sequence),
                    sequence_getLength(sequence), 1);
            const char *header = sequence_getHeader(sequence);
            fastaWrite(string, (char *)header, fileHandle);
            //fprintf(fileHandle, ">%s\n%s\n", (char *)header, string);
            free(string);
        }
    }
    stList_destruct(sequences);
}
