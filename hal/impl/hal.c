#include <ctype.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

#include "cactus.h"
#include "sonLib.h"
#include "recursiveThreadBuilder.h"
#include "hal3/hal3_c.h"

static Name globalReferenceEventName;

/*
 * Hal encodes a hierarchical alignment format.
 * Cactus outputs a text file that can be loaded into hal.
 * .c2h is as follows.
 *
 * sequence :
 *      sequenceLine segmentLines
 *
 * #Specifies a sequence
 * sequenceLine :
 *      "s\t'eventHeader'\t'sequenceHeader'\tisBottom\n"
 *
 * #Name taken from newick tree
 * eventHeader :
 *      string
 *
 * #Everything after the '>' in a fasta header,
 * sequenceHeader :
 *      string
 *
 * #Tells you is it is a bottom sequence, is 1 if so.
 * isBottom :
 *      0
 *      1
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
 * #If the isBottom field is 1 then all segments will be bottom segments in the sequence
 * bottomSegment :
 *      "a\tsegmentName\tstart\tlength\n"
 *
 * #Conversely, if the isBottom field is 0 then all segments will be top segments in the sequence
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
 *      0
 *      1
 */

static void writeSequenceHeader(FILE *fileHandle, Sequence *sequence) {
    //s eventName sequenceName isBottom
    Event *event = sequence_getEvent(sequence);
    assert(event != NULL);
    assert(event_getHeader(event) != NULL);
    assert(sequence_getHeader(sequence) != NULL);
    fprintf(fileHandle, "s\t'%s'\t'%s'\t%i\n", event_getHeader(event), sequence_getHeader(sequence),
            event_getName(event) == globalReferenceEventName);
}

static char *writeTerminalAdjacency(Cap *cap, void *extraArg) {
    //a start length reference-segment block-orientation
    Cap *adjacentCap = cap_getAdjacency(cap);
    assert(adjacentCap != NULL);
    int64_t adjacencyLength = cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) - 1;
    assert(adjacencyLength >= 0);
    if (adjacencyLength > 0) {
        Sequence *sequence = cap_getSequence(cap);
        assert(sequence != NULL);
        assert(cap_getEvent(cap) != NULL);
        if (event_getName(cap_getEvent(cap)) == globalReferenceEventName) {
            return stString_print("a\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\n", cap_getName(cap), cap_getCoordinate(cap) + 1 - sequence_getStart(sequence), adjacencyLength);
        }
        return stString_print("a\t%" PRIi64 "\t%" PRIi64 "\n", cap_getCoordinate(cap) + 1 - sequence_getStart(sequence), adjacencyLength);
    }
    else {
        return stString_copy("");
    }
}

static char *writeSegment(Segment *segment, void *extraArg) {
    Block *block = segment_getBlock(segment);
    Segment *referenceSegment = block_getSegmentForEvent(block, globalReferenceEventName);
    if (referenceSegment == NULL) {
        Cap *cap5 = segment_get5Cap(segment);
        Cap *cap3 = segment_get3Cap(segment);
        Sequence *sequence = cap_getSequence(cap5);
        return stString_print("a\t%" PRIi64 "\t%" PRIi64 "\n", cap_getCoordinate(cap5) - sequence_getStart(sequence), cap_getCoordinate(cap3) - cap_getCoordinate(cap5) + 1);
    }
    Sequence *sequence = segment_getSequence(segment);
    assert(sequence != NULL);
    Name eventName = event_getName(segment_getEvent(segment));
    if (referenceSegment != segment && eventName != globalReferenceEventName) { //Is a top segment
        return stString_print("a\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\n", segment_getStart(segment) - sequence_getStart(sequence), segment_getLength(segment), segment_getName(referenceSegment), segment_getStrand(referenceSegment));
    } else {
        //Is a bottom segment
        return stString_print("a\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\n", segment_getName(segment), segment_getStart(segment) - sequence_getStart(sequence), segment_getLength(segment));
    }
}

static int compareCaps(Cap *cap, Cap *cap2) {
    Event *event = cap_getEvent(cap);
    Event *event2 = cap_getEvent(cap2);
    int i = cactusMisc_nameCompare(event_getName(event), event_getName(event2));
    if (i != 0) {
        return event_getName(event) == globalReferenceEventName ? -1 : (event_getName(event2) == globalReferenceEventName ? 1 : i);
    }
    Sequence *sequence = cap_getSequence(cap);
    Sequence *sequence2 = cap_getSequence(cap2);
    i = cactusMisc_nameCompare(sequence_getName(sequence), sequence_getName(sequence2));
    if (i == 0) {
        i = cap_getCoordinate(cap) > cap_getCoordinate(cap2) ? 1 : (cap_getCoordinate(cap) < cap_getCoordinate(cap2) ? -1 : 0);
    }
    return i;
}

static stList *getCaps(Flower *flower) {
    //Get the caps in order
    stList *caps = stList_construct();
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
                        stList_append(caps, cap);
                    }
                }
            }
            end_destructInstanceIterator(capIt);
        }
    }
    flower_destructEndIterator(endIt);
    stList_sort(caps, (int(*)(const void *, const void *)) compareCaps);
    return caps;
}

void makeHalFormat(Flower *flower, stKVDatabase *database, Name referenceEventName, FILE *fileHandle) {
    globalReferenceEventName = referenceEventName;
    stList *caps = getCaps(flower);
    if (fileHandle == NULL) {
        buildRecursiveThreads(database, caps, writeSegment, writeTerminalAdjacency, NULL);
    } else {
        stList *threadStrings = buildRecursiveThreadsInList(database, caps, writeSegment, writeTerminalAdjacency, NULL);
        assert(stList_length(threadStrings) == stList_length(caps));
        for (int64_t i = 0; i < stList_length(threadStrings); i++) {
            Cap *cap = stList_get(caps, i);
            if(!sequence_isTrivialSequence(cap_getSequence(cap))) {
                char *threadString = stList_get(threadStrings, i);
                writeSequenceHeader(fileHandle, cap_getSequence(cap));
                fprintf(fileHandle, "%s\n", threadString);
            }
        }
    }
    stList_destruct(caps);
}

void makeHalFormatNoDb(Flower *flower, RecordHolder *rh, Name referenceEventName, FILE *fileHandle) {
    if (fileHandle == NULL) {
        globalReferenceEventName = referenceEventName;
        stList *caps = getCaps(flower);
        buildRecursiveThreadsNoDb(rh, caps, writeSegment, writeTerminalAdjacency, NULL);
        stList_destruct(caps);
    } else {
        makeHalOutputNoDb(flower, rh, referenceEventName, fileHandle, NULL, NULL);
    }
}

static void hal3Check(int rc, char *err, const char *fragmentDir) {
    if (rc != 0) {
        st_errAbort("writing the HAL format 3 fragment %s: %s", fragmentDir, err ? err : "unknown error");
    }
}

void makeHalOutputNoDb(Flower *flower, RecordHolder *rh, Name referenceEventName, FILE *fileHandle,
                       const char *fragmentDir, const char *fragmentTree) {
    globalReferenceEventName = referenceEventName;
    stList *caps = getCaps(flower);
    stList *threadStrings = buildRecursiveThreadsInListNoDb(rh, caps, writeSegment, writeTerminalAdjacency, NULL);
    assert(stList_length(threadStrings) == stList_length(caps));
    if (fileHandle != NULL) {
        for (int64_t i = 0; i < stList_length(threadStrings); i++) {
            Cap *cap = stList_get(caps, i);
            if(!sequence_isTrivialSequence(cap_getSequence(cap))) {
                char *threadString = stList_get(threadStrings, i);
                writeSequenceHeader(fileHandle, cap_getSequence(cap));
                fprintf(fileHandle, "%s\n", threadString);
            }
        }
    }
    if (fragmentDir != NULL) {
        // the same sequences, in the same order, as the .c2h: every one's length first,
        // then each one's segment lines and DNA
        char *err = NULL;
        hal3_c2h *h = hal3_c2h_open(fragmentDir, fragmentTree, 1, &err);
        hal3Check(h == NULL ? -1 : 0, err, fragmentDir);
        for (int64_t pass = 0; pass < 2; pass++) {
            for (int64_t i = 0; i < stList_length(threadStrings); i++) {
                Sequence *sequence = cap_getSequence((Cap *)stList_get(caps, i));
                if (sequence_isTrivialSequence(sequence)) {
                    continue;
                }
                Event *event = sequence_getEvent(sequence);
                const char *genome = event_getHeader(event), *name = sequence_getHeader(sequence);
                int64_t length = sequence_getLength(sequence);
                if (pass == 0) {
                    hal3Check(hal3_c2h_declare(h, genome, name, event_getName(event) == referenceEventName, length, &err),
                              err, fragmentDir);
                } else {
                    char *threadString = stList_get(threadStrings, i);
                    char *dna = sequence_getString(sequence, sequence_getStart(sequence), length, 1);
                    hal3Check(hal3_c2h_sequence(h, genome, name, threadString, strlen(threadString), dna, length, &err),
                              err, fragmentDir);
                    free(dna);
                }
            }
        }
        hal3Check(hal3_c2h_finish(h, &err), err, fragmentDir);
    }
    stList_destruct(threadStrings);
    stList_destruct(caps);
}