/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include "cactus.h"
#include "sonLib.h"
#include "blockConsensusString.h"

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

static void setAdjacencyLength(Cap *cap) {
    /*
     * Sets the coordinates of the caps to be equal to the length of the adjacency sequence between them.
     * Used to build the reference sequence bottom up.
     */
    Cap *adjacentCap = cap_getAdjacency(cap);
    assert(adjacentCap != NULL);
    assert(cap_getCoordinate(adjacentCap) == INT32_MAX);
    assert(cap_getStrand(cap) == cap_getStrand(adjacentCap));
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
}

void bottomUp(Flower *flower, Name referenceEventName) {
    /*
     * Run on each flower in turn, bottom up. Applies "setAjacencyLength" to to each
     * reference event adjacency in the given flower.
     */
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if(end_isBlockEnd(end) || end_isAttached(end)) {
            Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
            assert(cap != NULL);
            if (cap_getCoordinate(cap) == INT32_MAX) {
                setAdjacencyLength(cap);
            }
        }
    }
    flower_destructEndIterator(endIt);
}

static MetaSequence *addMetaSequence(Flower *flower, Cap *cap, int32_t *index) {
    /*
     * Adds a meta sequence representing a top level thread to the cactus disk.
     * The sequence is all 'N's at this stage.
     */
    int32_t sequenceLength = traceThreadLength(cap);
    char *sequence = st_malloc(sizeof(char) *(sequenceLength+1));
    for(int32_t i=0; i<sequenceLength; i++) {
        sequence[i] = 'N';
    }
    sequence[sequenceLength] = '\0';
    Event *referenceEvent = cap_getEvent(cap);
    assert(referenceEvent != NULL);
    char *sequenceName = stString_print("%s.%i", event_getHeader(referenceEvent),
                                (*index)++);
    MetaSequence *metaSequence = metaSequence_construct(1, sequenceLength, sequence,
                                sequenceName, event_getName(referenceEvent),
                                flower_getCactusDisk(flower));
    free(sequenceName);
    free(sequence);
    return metaSequence;
}

static stList *setCoordinates2(MetaSequence *metaSequence, stList *strings, int32_t coordinate) {
    if(stList_length(strings) > 0) {
        char *singleString = stString_join2("", strings);
        //int32_t stringLength = strlen(singleString);
        stList_destruct(strings);
        //metaSequence_setString(metaSequence, coordinate - stringLength + 1, stringLength, 1, singleString);
        free(singleString);
        strings = stList_construct3(0, free);
    }
    return strings;
}

static int32_t setCoordinates(Flower *flower, MetaSequence *metaSequence, Cap *cap, int32_t coordinate,
        Name outgroupEventName) {
    /*
     * Sets the coordinates of the reference thread and sets the bases of the actual sequence
     * that of the consensus.
     */
    Sequence *sequence = flower_getSequence(flower, metaSequence_getName(metaSequence));
    if(sequence == NULL) {
        sequence = sequence_construct(metaSequence, flower);
    }
    stList *strings = stList_construct3(0, free);
    while (1) {
        assert(!cap_getSide(cap));
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        if(!cap_getSide(adjacentCap)) {
            adjacentCap = cap_getReverse(adjacentCap);
        }
        assert(cap_getSide(adjacentCap));
        int32_t adjacencyLength = cap_getCoordinate(cap);
        assert(adjacencyLength == cap_getCoordinate(adjacentCap));
        assert(adjacencyLength != INT32_MAX);
        assert(adjacencyLength >= 0);
        cap_setCoordinates(cap, coordinate, 1, sequence);
        if(adjacencyLength > 0) {
            strings = setCoordinates2(metaSequence, strings, coordinate);
        }
        coordinate += adjacencyLength + 1;
        cap_setCoordinates(adjacentCap, coordinate, 1, sequence);
        cap_makeAdjacent(cap, adjacentCap); //This corrects the adjacency so side and strand are consistent.
        //Traverse any block..
        Segment *segment;
        if ((segment = cap_getSegment(adjacentCap)) != NULL) {
            /*
             * Set the segment string
             */
            stList_append(strings, getConsensusString(segment_getBlock(segment), outgroupEventName));
            coordinate += segment_getLength(segment)-1;
            cap = cap_getOtherSegmentCap(adjacentCap);
            assert(cap != NULL);
        } else {
            break;
        }
    }
    stList_destruct(setCoordinates2(metaSequence, strings, coordinate-1));
    return coordinate;
}

void addSequencesAndReferenceCoordinatesToTopLevelFlower(Flower *flower, Name referenceEventName, Name outgroupEventName) {
    /*
     * Adds reference (meta)sequences and then sets coordinates of the top level reference threads.
     */
    assert(!flower_hasParentGroup(flower));
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    int32_t index = 0;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end) && end_isAttached(end)) {
            Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
            if (cap_getSequence(cap) == NULL) {
                cap = cap_getSide(cap) ? cap_getReverse(cap) : cap;
                MetaSequence *metaSequence = addMetaSequence(flower, cap, &index);
                int32_t endCoordinate = setCoordinates(flower, metaSequence, cap, metaSequence_getStart(metaSequence)-1, outgroupEventName);
                (void)endCoordinate;
                assert(endCoordinate == metaSequence_getLength(metaSequence) + metaSequence_getStart(metaSequence));
            }
        }
    }
    flower_destructEndIterator(endIt);
}

void topDown(Flower *flower, Name referenceEventName, Name outgroupEventName) {
    /*
     * Run on each flower, top down. Sets the coordinates of each reference cap to the correct
     * sequence, and sets the bases of the reference sequence to be consensus bases.
     */
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if(end_isBlockEnd(end) || end_isAttached(end)) {
            Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            if (!cap_getSide(cap)) {
                assert(cap != NULL);
                Sequence *sequence = cap_getSequence(cap);
                assert(sequence != NULL);
                Group *group = end_getGroup(end);
                if (!group_isLeaf(group)) {
                    Flower *nestedFlower = group_getNestedFlower(group);
                    Cap *nestedCap = flower_getCap(nestedFlower, cap_getName(cap));
                    assert(nestedCap != NULL);
                    nestedCap = cap_getSide(nestedCap) ? cap_getReverse(nestedCap) : nestedCap;
                    assert(!cap_getSide(nestedCap));
                    int32_t endCoordinate = setCoordinates(nestedFlower, sequence_getMetaSequence(sequence), nestedCap, cap_getCoordinate(cap), outgroupEventName);
                    (void)endCoordinate;
                    assert(endCoordinate == cap_getCoordinate(cap_getAdjacency(cap)));
                    assert(endCoordinate == cap_getCoordinate(flower_getCap(nestedFlower, cap_getName(cap_getAdjacency(cap)))));
                }
            }
        }
    }
    flower_destructEndIterator(endIt);
}

