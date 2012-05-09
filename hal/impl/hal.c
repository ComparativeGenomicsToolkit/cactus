#include <ctype.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "cactus.h"
#include "sonLib.h"

#include "recursiveFileBuilder.h"

static Event *globalReferenceEvent;

/*
 * Hal encodes a hierarchical alignment format.
 * .hal is as follows:
 *
 * sequence :
 *      sequenceLine segmentLines
 *
 * #Specifies a sequence
 * sequenceLine :
 *      "s\t'eventName'\t'sequenceHeader'\tisBottom\n"
 *
 * #Name taken from newick tree
 * eventHeader :
 *      "string"
 *
 * #Everything after the '>' in a fasta header,
 * sequenceHeader :
 *      "string"
 *
 * #Tells you is it is a bottom sequence, is 1 if so.
 * isBottom :
 *      "0"
 *      "1"
 *
 * #List of segment lines
 * segmentLines :
 *      segmentLine
 *      segmentLine segmentLines
 *
 * #Definition of a segment
 * segmentLine :
 *      bottomSegment
 *      topSegment
 *
 * #Is the isBottom field is 1 then all segments will be bottom segments in the sequence
 * bottomSegment :
 *      "a\tsegmentName\tstart\tlength\n"
 *
 * #Conversely, is the isBottom field is 1 then all segments will be bottom segments in the sequence
 * topSegment :
 *      #If it has a parent
 *      "a\tstart\tlength\tparentSegment\talignmentOrientation\n"
 *      #If it was an insertion
 *      "a\tstart\tlength\n"
 *
 * #Start coordinate of segment
 * start :
 *      integer >= 0
 *
 * #Length of a segment
 * length :
 *      integer >= 1
 *
 * #Segment name
 * segmentName :
 *      integer >= 0
 *
 * #Name of segment to which it aligns
 * parentSegment :
 *      segmentName
 *
 * alignmentOrientation :
 *      "0"
 *      "1"
 */

static void writeSequenceHeader(FILE *fileHandle, Sequence *sequence) {
    //s eventName sequenceName isBottom
    Event *event = sequence_getEvent(sequence);
    assert(event != NULL);
    assert(event_getHeader(event) != NULL);
    assert(sequence_getHeader(sequence) != NULL);
    fprintf(fileHandle, "s\t'%s'\t'%s'\t%i\n", event_getHeader(event), sequence_getHeader(sequence),
            event_getName(event) == event_getName(globalReferenceEvent));
}

static void writeTerminalAdjacency(FILE *fileHandle, Cap *cap) {
    //a start length reference-segment block-orientation
    Cap *adjacentCap = cap_getAdjacency(cap);
    assert(adjacentCap != NULL);
    int32_t adjacencyLength = cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) - 1;
    assert(adjacencyLength >= 0);
    Sequence *sequence = cap_getSequence(cap);
    assert(sequence != NULL);
    fprintf(fileHandle, "a\t%i\t%i\n", cap_getCoordinate(cap) - sequence_getStart(sequence), adjacencyLength);
}

static void writeSegment(FILE *fileHandle, Segment *segment) {
    Block *block = segment_getBlock(segment);
    Segment *referenceSegment = block_getSegmentForEvent(block, event_getName(globalReferenceEvent));
    assert(referenceSegment != NULL);
    Sequence *sequence = segment_getSequence(segment);
    assert(sequence != NULL);
    if (referenceSegment != segment) { //Is a bottom segment
        fprintf(fileHandle, "a\t%i\t%i\t%" PRIi64 "\t%i\n", segment_getStart(segment) - sequence_getStart(sequence), segment_getLength(segment), segment_getName(referenceSegment), segment_getStrand(referenceSegment));
    }
    else { //Is a top segment
        fprintf(fileHandle, "a\t%" PRIi64 "\t%i\t%i\n", segment_getName(segment), segment_getStart(segment) - sequence_getStart(sequence), segment_getLength(segment));
    }
}

static int compareCaps(Cap *cap, Cap *cap2) {
    Event *event = cap_getEvent(cap);
    Event *event2 = cap_getEvent(cap2);
    int i = cactusMisc_nameCompare(event_getName(event), event_getName(event2));
    if (i != 0) {
        return event == globalReferenceEvent ? 1 : (event2 == globalReferenceEvent ? -1 : i);
    }
    Sequence *sequence = cap_getSequence(cap);
    Sequence *sequence2 = cap_getSequence(cap2);
    i = cactusMisc_nameCompare(sequence_getName(sequence), sequence_getName(sequence2));
    if (i != 0) {
        i = cap_getCoordinate(cap) > cap_getCoordinate(cap2) ? 1
                : (cap_getCoordinate(cap) < cap_getCoordinate(cap2) ? -1 : 0);
    }
    return i;
}

static stSortedSet *getCaps(Flower *flower, Name globalReferenceEventName) {
    //Get the caps in order
    stSortedSet *threadsToWrite = stSortedSet_construct3((int(*)(const void *, const void *)) compareCaps, NULL);
    End *end;
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end)) { // && end_isAttached(end)) {
            Cap *cap; // = end_getCapForEvent(end, globalReferenceEventName);
            End_InstanceIterator *capIt = end_getInstanceIterator(end);
            while ((cap = end_getNext(capIt)) != NULL) {
                if (cap_getSequence(cap) != NULL) {
                    cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                    if (!cap_getSide(cap)) {
                        stSortedSet_insert(threadsToWrite, cap);
                    }
                }
            }
            end_destructInstanceIterator(capIt);
        }
    }
    flower_destructEndIterator(endIt);
    return threadsToWrite;
}

void makeHalFormat(Flower *flower, RecursiveFileBuilder *recursiveFileBuilder, Event *referenceEvent,
        FILE *parentFileHandle, bool hasParent) {
    //Cheeky global
    globalReferenceEvent = referenceEvent;
    //Build list of threads in a sorted order
    stSortedSet *threadsToWrite = getCaps(flower, event_getName(globalReferenceEvent));
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
}
