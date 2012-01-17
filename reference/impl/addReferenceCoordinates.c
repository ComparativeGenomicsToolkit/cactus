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

static bool sequenceIsInverted(Cap *cap, Cap *adjacentCap) {
    assert(cactusMisc_nameCompare(cap_getName(cap), cap_getName(adjacentCap)) != 0);
    return cactusMisc_nameCompare(cap_getName(cap), cap_getName(adjacentCap)) > 0;
}

static char *getSequenceName(const char *sequenceDir, Name flowerName, Cap *cap, Cap *adjacentCap) {
    if(sequenceIsInverted(cap, adjacentCap)) {
        Cap *cap2 = adjacentCap;
        adjacentCap = cap;
        cap = cap2;
    }
    return stString_print("%s/%lli_%lli_%lli", sequenceDir, flowerName, cap_getName(cap), cap_getName(adjacentCap));
}

static char *readSequence(const char *sequenceFile) {
    return NULL;
}

static Cap *traceAdjacency(Cap *cap) {
    while(1) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        //Traverse any block..
        if(cap_getSegment(adjacentCap) != NULL) {
            cap = cap_getOtherSegmentCap(adjacentCap);
            assert(cap != NULL);
        }
        else {
            assert(end_isStubEnd(cap_getEnd(adjacentCap)));
            assert(end_isAttached(cap_getEnd(adjacentCap)));
            return adjacentCap;
        }
    }
}

static void addReferenceSequence(Cap *cap, const char *childSequenceDir,
        FILE *sequenceHandle, Name outgroupEventName) {
    while(1) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        Group *group = end_getGroup(cap_getEnd(cap));
        assert(group != NULL);
        int32_t adjacencySize = 0;
        if(!group_isLeaf(group)) { //Adjacency is not terminal, so establish its sequence.
            char *childSequenceName = getSequenceName(childSequenceDir, group_getName(group), cap, adjacentCap);
            char *childSequence = readSequence(childSequenceName);
            if(sequenceIsInverted(cap, adjacentCap)) {
                reverseComplementString(childSequence);
            }
            adjacencySize = strlen(childSequence);
            fprintf(sequenceHandle, "%s", childSequence);
            free(childSequence);
            free(childSequenceName);
        }
        //Set the coordinates of the caps to the adjacency size
        //cap_setCoordinate(cap, adjacencySize);
        //cap_setCoordinate(adjacentCap, adjacencySize);
        //Traverse any block..
        if(cap_getSegment(adjacentCap) != NULL) {
            char *blockString = getConsensusString(segment_getBlock(cap_getSegment(cap)), outgroupEventName);
            fprintf(sequenceHandle, "%s", blockString);
            free(blockString);
            cap = cap_getOtherSegmentCap(adjacentCap);
            assert(cap != NULL);
        }
        else {
            assert(end_isStubEnd(cap_getEnd(adjacentCap)));
            assert(end_isAttached(cap_getEnd(adjacentCap)));
            break;
        }
    }
}

void bottomUp(Flower *flower, const char *childSequenceDir, const char *sequenceDir,
        Name referenceEventName, Name outgroupEventName) {
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end) && end_isAttached(end)) {
            Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
            assert(cap != NULL);
            if (cap_getCoordinate(cap) == INT32_MAX) {
                Cap *adjacentCap = traceAdjacency(cap);
                char *sequenceFileName = getSequenceName(sequenceDir, flower_getName(flower), cap, adjacentCap);
                st_logDebug("Adding the coordinates for cap %s and building a sequence in %s\n",
                        cactusMisc_nameToStringStatic(cap_getName(cap)), sequenceFileName);
                FILE *sequenceHandle = fopen(sequenceFileName, "w");
                addReferenceSequence(cap, childSequenceDir, sequenceHandle, outgroupEventName);
                fclose(sequenceHandle);
                free(sequenceFileName);
            }
        }
    }
    flower_destructEndIterator(endIt);
}

static void setCoordinates(Flower *flower, MetaSequence *metaSequence, Cap *cap, int32_t coordinate) {
    Sequence *sequence = sequence_construct(metaSequence, flower);
    while(1) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        //Set the coordinates of the caps to the adjacency size
        assert(cap_getCoordinate(cap) == cap_getCoordinate(adjacentCap));
        assert(cap_getCoordinate(cap) != INT32_MAX);
        assert(cap_getCoordinate(cap) >= 1);
        cap_setCoordinates(cap, coordinate, 1, sequence);
        coordinate += cap_getCoordinate(adjacentCap);
        cap_setCoordinates(adjacentCap, coordinate, 1, sequence);
        //Traverse any block..
        Segment *segment;
        if((segment = cap_getSegment(adjacentCap)) != NULL) {
            coordinate += segment_getLength(segment);
            cap = cap_getOtherSegmentCap(adjacentCap);
            assert(cap != NULL);
        }
        else {
           break;
        }
    }
}

void addSequences(Flower *flower, const char *sequenceDir, Name referenceEventName) {
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    int32_t index = 0;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end) && end_isAttached(end)) {
            Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
            assert(cap != NULL);
            if (cap_getSequence(cap) == NULL) {
                Cap *adjacentCap = traceAdjacency(cap);
                if(sequenceIsInverted(cap, adjacentCap)) {
                    Cap *cap2 = adjacentCap;
                    adjacentCap = cap;
                    cap = cap2;
                }
                char *sequenceFileName = getSequenceName(sequenceDir, flower_getName(flower), cap, adjacentCap);
                char *sequence = readSequence(sequenceFileName);
                //Create the sequence
                Event *event = cap_getEvent(cap);
                char *sequenceName = stString_print("%s.%i", event_getHeader(event),
                            index++);
                MetaSequence *metaSequence = metaSequence_construct(1, strlen(sequence), sequence,
                            sequenceName, event_getName(event),
                            flower_getCactusDisk(end_getFlower(cap_getEnd(cap))));
                //Now do a traversal to fill in the coordinates
                setCoordinates(flower, metaSequence, cap, 1);
                free(sequence);
                free(sequenceFileName);
                free(sequenceName);
            }
        }
    }
    flower_destructEndIterator(endIt);
}

void topDown(Flower *flower, Name referenceEventName) {
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
        cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
        if(!cap_getSide(cap)) {
            assert(cap != NULL);
            Sequence *sequence = cap_getSequence(cap);
            assert(sequence != NULL);
            Group *group = end_getGroup(end);
            if(!group_isLeaf(group)) {
                Flower *nestedFlower = group_getNestedFlower(group);
                Cap *nestedCap = flower_getCap(nestedFlower, cap_getName(cap));
                assert(nestedCap != NULL);
                setCoordinates(nestedFlower, sequence_getMetaSequence(sequence), nestedCap, cap_getCoordinate(cap));
            }
        }
    }
    flower_destructEndIterator(endIt);
}

