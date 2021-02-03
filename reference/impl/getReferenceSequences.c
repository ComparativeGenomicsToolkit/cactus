#include "cactus.h"
#include "bioioC.h"

static char *formatSequenceHeader(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if (strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) * (1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    } else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

void getReferenceSequences(FILE *fileHandle, Flower *flower, char *referenceEventString){
    //get names of all the sequences in 'flower' for event with name 'referenceEventString'
    Sequence *sequence;
    Flower_SequenceIterator * seqIterator = flower_getSequenceIterator(flower);
    while((sequence = flower_getNextSequence(seqIterator)) != NULL)
    {
        Event* event = sequence_getEvent(sequence);
        const char* eventName = event_getHeader(event);
        if (strcmp(eventName, referenceEventString) == 0 &&
            sequence_getLength(sequence) > 0 &&
            !metaSequence_isTrivialSequence(sequence_getMetaSequence(sequence))) {
            char *sequenceHeader = formatSequenceHeader(sequence);
            st_logInfo("Sequence %s\n", sequenceHeader);
            char *string = sequence_getString(sequence, sequence_getStart(sequence), sequence_getLength(sequence), 1);
            fastaWrite(string, sequenceHeader, fileHandle);
            free(string);
            free(sequenceHeader);
        }
    }
    flower_destructSequenceIterator(seqIterator);
    return;
}
