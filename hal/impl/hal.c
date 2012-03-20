#include <ctype.h>
#include <assert.h>

#include "cactus.h"
#include "sonLib.h"

#include "recursiveFileBuilder.h"

static Event *referenceEvent;

static void writeSequenceHeader(FILE *fileHandle, Sequence *sequence) {
    //s eventName sequenceName isReference
    Event *event = sequence_getEvent(sequence);
    assert(event != NULL);
    assert(event_getHeader(event) != NULL);
    assert(sequence_getHeader(sequence) != NULL);
    fprintf(fileHandle, "s\t'%s'\t'%s'\t%i\n", event_getHeader(event), sequence_getHeader(sequence), event_getName(event) == event_getName(referenceEvent));
}

static void writeTerminalAdjacency(FILE *fileHandle, Cap *cap) {
    //a start length reference-segment block-orientation
    Cap *adjacentCap = cap_getAdjacency(cap);
    assert(adjacentCap != NULL);
    int32_t adjacencyLength = cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) - 1;
    assert(adjacencyLength >= 0);
    fprintf(fileHandle, "a\t%i\t%i\n", cap_getCoordinate(cap) + 1, adjacencyLength);
}

static void writeSegment(FILE *fileHandle, Segment *segment) {
    //a start length reference-segment block-orientation
    fprintf(fileHandle, "a\t%i\t%i", segment_getStart(segment), segment_getLength(segment));
    Block *block = segment_getBlock(segment);
    Segment *referenceSegment = block_getSegmentForEvent(block, event_getName(referenceEvent));
    assert(referenceSegment != NULL);
    if(referenceSegment != segment) {
        fprintf(fileHandle, "\t%s\t%i", cactusMisc_nameToStringStatic(segment_getName(referenceSegment)), segment_getStrand(referenceSegment));
    }
    fprintf(fileHandle, "\n");
}

static int compareCaps(Cap *cap, Cap *cap2) {
    Event *event = cap_getEvent(cap);
    Event *event2 = cap_getEvent(cap2);
    int i = cactusMisc_nameCompare(event_getName(event), event_getName(event2));
    if(i != 0) {
        return event == referenceEvent ? 1 : (event2 == referenceEvent ? -1 : i);
    }
    Sequence *sequence = cap_getSequence(cap);
    Sequence *sequence2 = cap_getSequence(cap2);
    i = cactusMisc_nameCompare(sequence_getName(sequence), sequence_getName(sequence2));
    if(i != 0) {
        i = cap_getCoordinate(cap) > cap_getCoordinate(cap2) ? 1 : (cap_getCoordinate(cap) < cap_getCoordinate(cap2) ? -1 : 0);
    }
    return i;
}

static stSortedSet *getCaps(Flower *flower, Name referenceEventName) {
    //Get the caps in order
    stSortedSet *threadsToWrite = stSortedSet_construct3((int(*)(const void *, const void *)) compareCaps, NULL);
    End *end;
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end) && end_isAttached(end)) {
            Cap *cap = end_getCapForEvent(end, referenceEventName);
            assert(cap != NULL);
            assert(cap_getSequence(cap) != NULL);
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            if (!cap_getSide(cap)) {
                stSortedSet_insert(threadsToWrite, cap);
            }
        }
    }
    flower_destructEndIterator(endIt);
    return threadsToWrite;
}

void makeHalFormat(Flower *flower, const char *referenceEventString, const char *childDirectory,
        const char *outputFile, bool hasParent) {
    //Basic objects
    referenceEvent = eventTree_getEventByHeader(flower_getEventTree(flower), referenceEventString);
    assert(referenceEvent != NULL);

    //File in which to place output
    FILE *parentFileHandle = fopen(outputFile, "w");
    //Structure for building output
    RecursiveFileBuilder *recursiveFileBuilder = recursiveFileBuilder_construct(childDirectory, parentFileHandle,
            hasParent);

    //Build list of threads in a sorted order
    stSortedSet *threadsToWrite = getCaps(flower, event_getName(referenceEvent));
    stSortedSetIterator *setIt = stSortedSet_getIterator(threadsToWrite);
    Cap *cap;
    while ((cap = stSortedSet_getNext(setIt)) != NULL) {
        //s Event-name Seq-name Seq-header
        if (!hasParent) {
            writeSequenceHeader(parentFileHandle, cap_getSequence(cap));
        }
        //For each segment:
        //a start length reference-segment block-orientation
        recursiveFileBuilder_writeThread(recursiveFileBuilder, cap, writeSegment, writeTerminalAdjacency);
        if (!hasParent) { //Terminate with new-line
            //new-line
            fprintf(parentFileHandle, "\n");
        }
    }
    stSortedSet_destructIterator(setIt);

    //Cleanup
    fclose(parentFileHandle);
    recursiveFileBuilder_destruct(recursiveFileBuilder);
}
