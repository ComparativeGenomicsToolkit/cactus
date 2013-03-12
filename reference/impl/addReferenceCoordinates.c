/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include "cactus.h"
#include "sonLib.h"
#include "blockConsensusString.h"
#include "recursiveThreadBuilder.h"

Cap *getCapForReferenceEvent(End *end, Name referenceEventName) {
    /*
     * Get the cap for a given event.
     */
    End_InstanceIterator *it = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(it)) != NULL) {
        if (event_getName(cap_getEvent(cap)) == referenceEventName) {
            end_destructInstanceIterator(it);
            return cap;
        }
    }
    end_destructInstanceIterator(it);
    assert(0);
    return NULL;
}

static int32_t traceThreadLength(Cap *cap) {
    /*
     * Gets the length in bases of the thread in the flower, starting from a given attached stub cap.
     * The thread length includes the lengths of adjacencies that it contains.
     */
    assert(end_isStubEnd(cap_getEnd(cap)) && end_isAttached(cap_getEnd(cap)));
    int32_t threadLength = 0;
    while (1) {
        assert(cap_getCoordinate(cap) != INT32_MAX);
        int32_t adjacencyLength = cap_getCoordinate(cap);
        threadLength += adjacencyLength;
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        assert(adjacencyLength == cap_getCoordinate(adjacentCap));
        //Traverse any block..
        if (cap_getSegment(adjacentCap) != NULL) {
            threadLength += segment_getLength(cap_getSegment(adjacentCap));
            cap = cap_getOtherSegmentCap(adjacentCap);
            assert(cap != NULL);
        } else {
            assert(end_isStubEnd(cap_getEnd(adjacentCap)));
            assert(end_isAttached(cap_getEnd(adjacentCap)));
            return threadLength;
        }
    }
}

static void setAdjacencyLengths(Cap *cap) {
    /*
     * Sets the coordinates of the caps to be equal to the length of the adjacency sequence between them.
     * Used to build the reference sequence bottom up.
     */
    while(1) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        assert(cap_getCoordinate(cap) == INT32_MAX);
        assert(cap_getCoordinate(adjacentCap) == INT32_MAX);
        assert(cap_getStrand(cap) == cap_getStrand(adjacentCap));
        assert(cap_getSide(cap) != cap_getSide(adjacentCap));
        Group *group = end_getGroup(cap_getEnd(cap));
        assert(group != NULL);
        int32_t adjacencyLength = 0;
        if (!group_isLeaf(group)) { //Adjacency is not terminal, so establish its sequence.
            Flower *nestedFlower = group_getNestedFlower(group);
            Cap *nestedCap = flower_getCap(nestedFlower, cap_getName(cap));
            assert(nestedCap != NULL);
            adjacencyLength = traceThreadLength(nestedCap);
        }
        //Set the coordinates of the caps to the adjacency size
        cap_setCoordinates(cap, adjacencyLength, cap_getStrand(cap), NULL);
        cap_setCoordinates(adjacentCap, adjacencyLength, cap_getStrand(adjacentCap), NULL);
        assert(cap_getCoordinate(cap) == cap_getCoordinate(adjacentCap));
        if((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
            break;
        }
    }
}

static char *terminalAdjacencyWriteFn(Cap *cap) {
    return stString_copy("");
}

static Name segmentWriteFnOutgroupEventName;

static char *segmentWriteFn(Segment *segment) {
    return getConsensusString(segment_getBlock(segment), segmentWriteFnOutgroupEventName);
}

static MetaSequence *addMetaSequence(Flower *flower, Cap *cap, int32_t index, char *string) {
    /*
     * Adds a meta sequence representing a top level thread to the cactus disk.
     * The sequence is all 'N's at this stage.
     */
    Event *referenceEvent = cap_getEvent(cap);
    assert(referenceEvent != NULL);
    char *sequenceName = stString_print("%s.refChr%i", event_getHeader(referenceEvent), index);
    //char *sequenceName = stString_print("refChr%i", index);
    MetaSequence *metaSequence = metaSequence_construct(1, strlen(string), string,
                                sequenceName, event_getName(referenceEvent),
                                flower_getCactusDisk(flower));
    free(sequenceName);
    return metaSequence;
}

static int32_t setCoordinates(Flower *flower, MetaSequence *metaSequence, Cap *cap, int32_t coordinate) {
    /*
     * Sets the coordinates of the reference thread and sets the bases of the actual sequence
     * that of the consensus.
     */
    Sequence *sequence = flower_getSequence(flower, metaSequence_getName(metaSequence));
    if(sequence == NULL) {
        sequence = sequence_construct(metaSequence, flower);
    }
    while (1) {
        assert(cap_getStrand(cap));
        assert(!cap_getSide(cap));
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        assert(cap_getStrand(adjacentCap));
        assert(cap_getSide(adjacentCap));
        int32_t adjacencyLength = cap_getCoordinate(cap);
        assert(adjacencyLength == cap_getCoordinate(adjacentCap));
        assert(adjacencyLength != INT32_MAX);
        assert(adjacencyLength >= 0);
        cap_setCoordinates(cap, coordinate, 1, sequence);
        coordinate += adjacencyLength + 1;
        cap_setCoordinates(adjacentCap, coordinate, 1, sequence);
        //Traverse any block..
        if ((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
            break;
        }
        coordinate += segment_getLength(cap_getSegment(adjacentCap)) - 1;
    }
    return coordinate;
}

static stList *getCaps(stList *flowers, Name referenceEventName) {
    stList *caps = stList_construct();
    for(int32_t i=0; i<stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        //Get list of caps
        Flower_EndIterator *endIt = flower_getEndIterator(flower);
        End *end;
        while ((end = flower_getNextEnd(endIt)) != NULL) {
            if (end_isStubEnd(end) && end_isAttached(end)) {
                Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
                assert(cap != NULL);
                cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                if (!cap_getSide(cap)) {
                    stList_append(caps, cap);
                }
            }
        }
        flower_destructEndIterator(endIt);
    }
    return caps;
}

void bottomUp(stList *flowers, stKVDatabase *sequenceDatabase, Name referenceEventName, Name outgroupEventName, bool isTop) {
    stList *caps = getCaps(flowers, referenceEventName);
    for(int32_t i=0; i<stList_length(caps); i++) {
        setAdjacencyLengths(stList_get(caps, i));
    }
    segmentWriteFnOutgroupEventName = outgroupEventName;
    if(isTop) {
        stList *threadStrings = buildRecursiveThreadsInList(sequenceDatabase, caps, segmentWriteFn, terminalAdjacencyWriteFn);
        assert(stList_length(threadStrings) == stList_length(caps));
        for(int32_t i=0; i<stList_length(threadStrings); i++) {
            Cap *cap = stList_get(caps, i);
            assert(cap_getStrand(cap));
            assert(!cap_getSide(cap));
            Flower *flower = end_getFlower(cap_getEnd(cap));
            MetaSequence *metaSequence = addMetaSequence(flower, cap, i, stList_get(threadStrings, i));
            int32_t endCoordinate = setCoordinates(flower, metaSequence, cap, metaSequence_getStart(metaSequence)-1);
            (void)endCoordinate;
            assert(endCoordinate == metaSequence_getLength(metaSequence) + metaSequence_getStart(metaSequence));
        }
        stList_destruct(threadStrings);
    }
    else {
        buildRecursiveThreads(sequenceDatabase, caps, segmentWriteFn, terminalAdjacencyWriteFn);
    }
    stList_destruct(caps);
}

void topDown(stList *flowers, Name referenceEventName) {
    /*
     * Run on each flower, top down. Sets the coordinates of each reference cap to the correct
     * sequence, and sets the bases of the reference sequence to be consensus bases.
     */
    for(int32_t i=0; i<stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        Flower_EndIterator *endIt = flower_getEndIterator(flower);
        End *end;
        while ((end = flower_getNextEnd(endIt)) != NULL) {
            if(end_isBlockEnd(end) || end_isAttached(end)) {
                Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
                assert(cap != NULL);
                cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                if (!cap_getSide(cap)) {
                    Sequence *sequence = cap_getSequence(cap);
                    assert(sequence != NULL);
                    assert(cap_getCoordinate(cap) != INT32_MAX);
                    Group *group = end_getGroup(end);
                    if (!group_isLeaf(group)) {
                        Flower *nestedFlower = group_getNestedFlower(group);
                        Cap *nestedCap = flower_getCap(nestedFlower, cap_getName(cap));
                        assert(nestedCap != NULL);
                        nestedCap = cap_getStrand(nestedCap) ? nestedCap : cap_getReverse(nestedCap);
                        assert(cap_getStrand(nestedCap));
                        assert(!cap_getSide(nestedCap));
                        int32_t endCoordinate = setCoordinates(nestedFlower, sequence_getMetaSequence(sequence), nestedCap, cap_getCoordinate(cap));
                        (void)endCoordinate;
                        assert(endCoordinate == cap_getCoordinate(cap_getAdjacency(cap)));
                        assert(endCoordinate == cap_getCoordinate(flower_getCap(nestedFlower, cap_getName(cap_getAdjacency(cap)))));
                    }
                }
            }
        }
        flower_destructEndIterator(endIt);
    }
}
