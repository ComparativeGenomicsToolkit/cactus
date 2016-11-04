#include <ctype.h>
#include <assert.h>

#include "cactus.h"
#include "sonLib.h"

#include "recursiveThreadBuilder.h"

static char *writeMafHeaderLine(Block *block) {
    /*
     * Write the header for a MAF file.
     */
    if (block_getRootInstance(block) != NULL) {
        /* Get newick tree string with internal labels and no unary events */
        char *newickTreeString = block_makeNewickString(block, 1, 0);
        assert(newickTreeString != NULL);
        char *string = stString_print("a score=%" PRIi64 " tree='%s'\n",
                block_getLength(block) * block_getInstanceNumber(block),
                newickTreeString);
        free(newickTreeString);
        return string;
    } else {
        return stString_print("a score=%" PRIi64 "\n",
                block_getLength(block) * block_getInstanceNumber(block));
    }
}

static char *formatSequenceHeader(Sequence *sequence) {
    /*
     * Format the header of a sequence.
     */
    Event *event = sequence_getEvent(sequence);
    const char *eventHeader = event_getHeader(event);
    const char *sequenceHeader = sequence_getHeader(sequence);
    if (strlen(sequenceHeader) > 0) {
        char *header = stString_print("%s.%s", eventHeader, sequenceHeader);
        return header;
    } else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

static char *writeMafSequenceLine2(Segment *segment,
        char *segmentString) {
    Sequence *sequence = segment_getSequence(segment);
    assert(sequence != NULL);
    char *sequenceHeader = formatSequenceHeader(sequence);
    int64_t start;
    if (segment_getStrand(segment)) {
        start = segment_getStart(segment) - sequence_getStart(sequence);
    } else { //start with respect to the start of the reverse complement sequence
        start
                = (sequence_getStart(sequence) + sequence_getLength(sequence)
                        - 1) - segment_getStart(segment);
    }
    int64_t length = segment_getLength(segment);
    char *strand = segment_getStrand(segment) ? "+" : "-";
    int64_t sequenceLength = sequence_getLength(sequence);
    char *string = stString_print("s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%" PRIi64 "\t%s\n", sequenceHeader, start,
            length, strand, sequenceLength, segmentString);
    free(sequenceHeader);
    return string;
}

static char *writeMafSequenceLine(Segment *segment) {
    /*
     * Write a maf sequence line for the given segment.
     */
    char *segmentString = segment_getString(segment);
    char *string = writeMafSequenceLine2(segment, segmentString);
    free(segmentString);
    return string;
}

static char *writeMafSequenceLineShowingOnlyReferenceSubstitutions(Segment *segment, const char *referenceString) {
    /*
     * Write a maf sequence line for the given segment, but showing only differences from reference sequence.
     */
    char *segmentString = segment_getString(segment);
    for (int64_t i = 0; i < segment_getLength(segment); i++) {
        if (toupper(segmentString[i]) == toupper(referenceString[i])) {
            segmentString[i] = '*';
        }
    }
    char *string = writeMafSequenceLine2(segment, segmentString);
    free(segmentString);
    return string;
}

static bool writeMafBlock_showOnlySubstitutionsWithRespectToReference = 0;

static char *writeMafBlock(Segment *referenceSegment) {
    /*
     * Writes a maf block.
     */
    Block *block = segment_getBlock(referenceSegment);
    stList *strings = stList_construct3(0, free);
    stList_append(strings, writeMafHeaderLine(block));
    stList_append(strings, writeMafSequenceLine(referenceSegment));
    char *referenceString = segment_getString(referenceSegment);
    Block_InstanceIterator *iterator = block_getInstanceIterator(block);
    Segment *segment;
    while ((segment = block_getNext(iterator)) != NULL) {
        if (segment != referenceSegment && segment_getSequence(segment) != NULL) {
            if (writeMafBlock_showOnlySubstitutionsWithRespectToReference) {
                stList_append(strings, writeMafSequenceLineShowingOnlyReferenceSubstitutions(
                        segment, referenceString));
            } else {
                stList_append(strings, writeMafSequenceLine(segment));
            }
        }
    }
    stList_append(strings, stString_print("\n"));
    block_destructInstanceIterator(iterator);
    free(referenceString);
    char *string = stString_join2("", strings);
    stList_destruct(strings);
    return string;
}

char *writeTerminalAdjacency(Cap *cap) {
    return stString_copy("");
}

void makeMAFHeader(Flower *flower, FILE *fileHandle) {
    fprintf(fileHandle, "##maf version=1 scoring=N/A\n");
    char *cA = eventTree_makeNewickString(flower_getEventTree(flower));
    fprintf(fileHandle, "# cactus %s\n\n", cA);
    free(cA);
}

stList *getCaps(stList *flowers, Name referenceEventName) {
    stList *caps = stList_construct();
    for(int64_t i=0; i<stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        End *end;
        Flower_EndIterator *endIt = flower_getEndIterator(flower);
        while ((end = flower_getNextEnd(endIt)) != NULL) {
            if (end_isStubEnd(end)) {
                Cap *cap = end_getCapForEvent(end, referenceEventName);
                if (cap != NULL) {
                    assert(cap_getSequence(cap) != NULL);
                    cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                    if (!cap_getSide(cap)) {
                        stList_append(caps, cap);
                    }
                }
            }
        }
        flower_destructEndIterator(endIt);
    }
    return caps;
}

void makeMafFormat(stList *flowers, stKVDatabase *database, Name referenceEventName,
        FILE *fileHandle, bool showOnlySubstitutionsWithRespectToReference) {
    writeMafBlock_showOnlySubstitutionsWithRespectToReference
                = showOnlySubstitutionsWithRespectToReference;
    stList *caps = getCaps(flowers, referenceEventName);
    if(fileHandle == NULL) {
        buildRecursiveThreads(database, caps, writeMafBlock, writeTerminalAdjacency);
    }
    else {
        makeMAFHeader(stList_get(flowers, 0), fileHandle);
        stList *threadStrings = buildRecursiveThreadsInList(database, caps, writeMafBlock, writeTerminalAdjacency);
        assert(stList_length(threadStrings) == stList_length(caps));
        for(int64_t i=0; i<stList_length(threadStrings); i++) {
            char *threadString = stList_get(threadStrings, i);
            fprintf(fileHandle, "%s\n", threadString);
        }
        stList_destruct(threadStrings);
    }
    stList_destruct(caps);
}
