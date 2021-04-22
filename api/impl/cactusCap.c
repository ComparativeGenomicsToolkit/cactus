/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic cap functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Bit twiddling using the "bits" char of the cap
 * 0: strand
 * 1: forward (is the forward as opposed to the reverse copy)
 * 2: part_of_segment (is part of a segment)
 * 3: is_segment (is the segment)
 * 4: left (is the left side of a segment)
 * 5: event_not_sequence (stores the event pointer, not the sequence itself)
*/

static void cap_setBit(Cap *cap, int bit, bool value) {
    cap->bits &= ~(1UL << bit); // first clear the existing value
    if(value) {
        cap->bits |= 1UL << bit;// now set the new value
    }
}

static void cap_setBitForwardAndReverse(Cap *cap, int bit, bool value, bool invertReverse) {
    cap_setBit(cap, bit, value);
    cap_setBit(cap_getReverse(cap), bit, value ^ invertReverse);
}

static bool cap_getBit(Cap *cap, int bit) {
    return (cap->bits >> bit) & 1;
}

bool cap_getStrand(Cap *cap) {
    return cap_getBit(cap, 0);
}

bool cap_forward(Cap *cap) {
    return cap_getBit(cap, 1);
}

bool cap_partOfSegment(Cap *cap) {
    return cap_getBit(cap, 2);
}

bool cap_isSegment(Cap *cap) {
    return cap_getBit(cap, 3);
}

bool cap_left(Cap *cap) {
    return cap_getBit(cap, 4);
}

bool cap_getHasEventNotSequence(Cap *cap) {
    return cap_getBit(cap, 5);
}

SegmentCapContents *cap_getSegmentContents(Cap *cap) {
    assert(cap_partOfSegment(cap));
    switch(cap->bits & 0x1E) {
        case 0x16: // binary: 010110
            return (SegmentCapContents *)(cap+6);
        case 0x14: // binary: 010100
            return (SegmentCapContents *)(cap+5);
        case 0xE: // binary: 001110
            return (SegmentCapContents *)(cap+4);
        case 0xC: // binary: 001100
            return (SegmentCapContents *)(cap+3);
        case 0x6: // binary: 000110
            return (SegmentCapContents *)(cap+2);
        case 0x4: // binary: 000100
            return (SegmentCapContents *)(cap+1);
        default:
            assert(0);
            return NULL;
    }
}

CapContents *cap_getContents(Cap *cap) {
    assert(!cap_partOfSegment(cap));
    switch(cap->bits & 0x1E) {
        case 0x2: // binary: 000010
            return (CapContents *)(cap+2);
        case 0x0: // binary: 000000
            return (CapContents *)(cap+1);
        default:
            assert(0);
            return NULL;
    }
}

CapCoreContents *cap_getCoreContents(Cap *cap) {
    if(cap_partOfSegment(cap)) {
        return &(cap_getSegmentContents(cap)->core);
    }
    return &(cap_getContents(cap)->core);
}

Cap *cap_construct(End *end, Event *event) {
    return cap_construct3(cactusDisk_getUniqueID(flower_getCactusDisk(end_getFlower(end))), event, end);
}

Cap *cap_construct2(End *end, int64_t coordinate, bool strand, Sequence *sequence) {
    Cap *cap = cap_construct3(cactusDisk_getUniqueID(flower_getCactusDisk(end_getFlower(end))),
                              sequence_getEvent(sequence), end);
    cap_setCoordinates(cap, coordinate, strand, sequence);
    return cap;
}

static Cap *cap_construct5(Name instance, Event *event, End *end, bool setFlower) {
    assert(instance != NULL_NAME);
    assert(event != NULL);
    assert(end != NULL);
    assert(!end_partOfBlock(end));

    // Create the combined forward and reverse caps
    Cap *cap = st_calloc(1, 2*sizeof(Cap) + sizeof(CapContents));

    // see above comment to decode what is set
    // Bits: strand / forward / part_of_segment / is_segment / left / event_not_sequence
    cap->bits = 0x23; // binary: 100011
    (cap+1)->bits = 0x20; // binary: 100000

    // Set the attributes of the combined caps/segments
    cap_getCoreContents(cap)->instance = instance;
    cap_getCoreContents(cap)->coordinate = INT64_MAX;
    cap_getCoreContents(cap)->eventOrSequence = event;
    cap_getContents(cap)->end = end;

    // Add the cap to the appropriate indexes
    end_addInstance(end, cap);
    if(setFlower) {
        flower_addCap(end_getFlower(end), cap);
    }

    // Do the checks
    assert(cap_getStrand(cap));
    assert(!cap_getStrand(cap_getReverse(cap)));
    assert(cap_getReverse(cap_getReverse(cap)) == cap); // check reversal works
    assert(cap_getContents(cap) == cap_getContents(cap_getReverse(cap)));
    assert(cap_getName(cap) == instance);
    assert(cap_getCoordinate(cap) == INT64_MAX);
    assert(cap_getSequence(cap) == NULL);
    assert(cap_getAdjacency(cap) == NULL);
    assert(cap_getSegment(cap) == NULL);
    assert(cap_getEvent(cap) == event);
    assert(cap_getEnd(cap) == end);
    assert(cap_getName(cap_getReverse(cap)) == instance);
    assert(cap_getCoordinate(cap_getReverse(cap)) == INT64_MAX);
    assert(cap_getSequence(cap_getReverse(cap)) == NULL);
    assert(cap_getAdjacency(cap_getReverse(cap)) == NULL);
    assert(cap_getSegment(cap_getReverse(cap)) == NULL);
    assert(cap_getEvent(cap_getReverse(cap)) == event);
    assert(cap_getEnd(cap_getReverse(cap)) == end_getReverse(end));
    if(setFlower) {
        assert(end_getInstance(end, instance) == cap);
        assert(end_getInstance(end_getReverse(end), instance) == cap_getReverse(cap));
        assert(flower_getCap(end_getFlower(end), instance) == cap_getPositiveOrientation(cap));
    }

    return cap;
}

Cap *cap_construct3(Name instance, Event *event, End *end) {
    return cap_construct5(instance, event, end, 1);
}

void cap_setCoordinates(Cap *cap, int64_t coordinate, bool strand, Sequence *sequence) {
    assert(!cap_isSegment(cap));

    // Set the strand for all caps
    cap_setBitForwardAndReverse(cap, 0, strand, 1); // 0 is the strand bit
    if(cap_partOfSegment(cap)) {
        cap_setBitForwardAndReverse(cap_getSegment(cap), 0, strand, 1);
        cap_setBitForwardAndReverse(cap_getOtherSegmentCap(cap), 0, strand, 1);
    }

    // Set the sequence for all caps
    if (sequence != NULL) {
        // Switch to having a sequence instead of an event
        cap_getCoreContents(cap)->eventOrSequence = sequence;
        // Flip the flag for all the caps
        cap_setBitForwardAndReverse(cap, 5, 0, 0); // 5 is the eventNotSequence bit
        if(cap_partOfSegment(cap)) {
            cap_setBitForwardAndReverse(cap_getSegment(cap), 5, 0, 0);
            cap_setBitForwardAndReverse(cap_getOtherSegmentCap(cap), 5, 0, 0);
        }
    }

    // Set the coordinate
    if (cap_partOfSegment(cap)) {
        int64_t i = coordinate;
        if (!cap_left(cap) && coordinate != INT64_MAX) { // Only set if at the left end
            int64_t j = segment_getLength(cap_getSegment(cap))-1;
            i += cap_forward(cap) ^ cap_getStrand(cap) ? j : -j; // Set the coordinate ensuring we set it with respect to the
            // forward copy and strand
        }
        cap_getCoreContents(cap)->coordinate = i;
    } else {
        cap_getCoreContents(cap)->coordinate = coordinate;
    }

    assert(cap_getStrand(cap) == strand);
    assert(cap_getCoordinate(cap) == coordinate);
    assert(cap_getSequence(cap) == sequence);
}

int64_t cap_getCoordinate(Cap *cap) {
    assert(!cap_isSegment(cap));
    if(cap_partOfSegment(cap)) {
        int64_t c = cap_getCoreContents(cap)->coordinate;
        if(!cap_left(cap) && c != INT64_MAX) {
            assert(segment_getLength(cap_getSegment(cap)) > 0);
            //  The plus or minus here is dependent on if the cap if forwards and the strand
            int64_t j = segment_getLength(cap_getSegment(cap))-1;
            c += cap_forward(cap) ^ cap_getStrand(cap) ? -j : j;
        }
        return c;
    }
    return cap_getCoreContents(cap)->coordinate;
}

Cap *cap_copyConstruct2(End *end, Cap *cap, bool setFlower) {
    assert(end_getName(cap_getEnd(cap)) == end_getName(end));
    assert(end_getSide(end) == cap_getSide(cap));

    Flower *flower = end_getFlower(end);
    if (cap_getSequence(cap) != NULL) { //cap_getCoordinate(cap) != INT64_MAX) {
        Name sequenceName = sequence_getName(cap_getSequence(cap));
        Sequence *sequence = flower_getSequence(flower, sequenceName);
        if (sequence == NULL) { //add sequence to the flower.
            sequence = cactusDisk_getSequence(flower_getCactusDisk(flower), sequenceName);
            flower_addSequence(flower, sequence);
            assert(sequence != NULL);
        }
        Cap *newCap = cap_construct5(cap_getName(cap), sequence_getEvent(sequence), end, setFlower);
        cap_setCoordinates(newCap, cap_getCoordinate(cap), cap_getStrand(cap), sequence);
        return newCap;
    } else {
        Event *event = eventTree_getEvent(flower_getEventTree(flower), event_getName(cap_getEvent(cap)));
        assert(event != NULL);
        Cap *cap2 = cap_construct5(cap_getName(cap), event, end, setFlower);
        cap_setCoordinates(cap2, cap_getCoordinate(cap), cap_getStrand(cap), NULL);
        return cap2;
    }
}

Cap *cap_copyConstruct(End *end, Cap *cap) {
    return cap_copyConstruct2(end, cap, 1);
}

void cap_destruct(Cap *cap) {
    //Remove from end.
    end_removeInstance(cap_getEnd(cap), cap);

    // Free only if not part of a segment
    if(!cap_partOfSegment(cap)) {
        free(cap_forward(cap) ? cap : cap_getReverse(cap));
    }
}

Name cap_getName(Cap *cap) {
    Name name = cap_getCoreContents(cap)->instance;
    if(cap_partOfSegment(cap) && !cap_left(cap)) {
        name += cap_isSegment(cap) ? 1 : 2;
    }
    return name;
}

End *cap_getEnd(Cap *cap) {
    assert(!cap_isSegment(cap));
    End *e;
    if(cap_partOfSegment(cap)) {
        Block *b = cap_getSegmentContents(cap)->block;
        e = cap_left(cap) ? block_get5End(b) : block_get3End(b);
    }
    else {
        e = cap_getContents(cap)->end;
    }
    return cap_forward(cap) ? e : end_getReverse(e);
}

bool cap_getOrientation(Cap *cap) {
    return end_getOrientation(cap_getEnd(cap));
}

Cap *cap_getPositiveOrientation(Cap *cap) {
    return cap_getOrientation(cap) ? cap : cap_getReverse(cap);
}

Cap *cap_getReverse(Cap *cap) {
    return cap_forward(cap) ? cap+1 : cap-1;
}

Event *cap_getEvent(Cap *cap) {
    void *e = cap_getCoreContents(cap)->eventOrSequence;
    return cap_getHasEventNotSequence(cap) ? e : sequence_getEvent(e);
}

Segment *cap_getSegment(Cap *cap) {
    assert(!cap_isSegment(cap));
    if(!cap_partOfSegment(cap)) {
        return NULL;
    }
    return (cap_left(cap) ? cap + 2 : cap - 2);
}

Cap *cap_getOtherSegmentCap(Cap *cap) {
    assert(!cap_isSegment(cap));
    if(!cap_partOfSegment(cap)) {
        return NULL;
    }
    Cap *otherCap = cap_left(cap) ? cap+4 : cap-4;
    // Do some sanity checks
    assert(cap != otherCap);
    //assert(cap_getStrand(cap) == cap_getStrand(otherCap)); This is not necessarily true when setting the coordinates
    assert(cap_getSegment(cap) == cap_getSegment(otherCap));
    assert(cap_getSegmentContents(cap) == cap_getSegmentContents(otherCap));

    return otherCap;
}

bool cap_getSide(Cap *cap) {
    return end_getSide(cap_getEnd(cap));
}

Sequence *cap_getSequence(Cap *cap) {
    return cap_getHasEventNotSequence(cap) ? NULL : cap_getCoreContents(cap)->eventOrSequence;
}

static Cap **cap_getAdjacencyP(Cap *cap) {
    if(cap_partOfSegment(cap)) {
        assert(!cap_isSegment(cap));
        return cap_left(cap) ? &(cap_getSegmentContents(cap)->leftAdjacency) : &(cap_getSegmentContents(cap)->rightAdjacency);
    }
    else {
        return &(cap_getContents(cap)->adjacency);
    }
}

void cap_makeAdjacent(Cap *cap, Cap *cap2) {
    //We put them both on the same strand, as the strand is important in the pairing
    cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
    cap2 = cap_getStrand(cap2) ? cap2 : cap_getReverse(cap2);
    assert(cap != cap2);
    assert(cap_getEvent(cap) == cap_getEvent(cap2));
    cap_breakAdjacency(cap);
    cap_breakAdjacency(cap2);
    *cap_getAdjacencyP(cap) = cap_forward(cap) ? cap2 : cap_getReverse(cap2);  // store cap orientation according to forward copy
    *cap_getAdjacencyP(cap2) = cap_forward(cap2) ? cap : cap_getReverse(cap);
}


Cap *cap_getAdjacency(Cap *cap) {
    Cap *connectedCap = *cap_getAdjacencyP(cap);
    return connectedCap == NULL ? NULL : cap_forward(cap) ? connectedCap : cap_getReverse(connectedCap);
}

void cap_breakAdjacency(Cap *cap) {
    Cap **cap2 = cap_getAdjacencyP(cap);
    if (*cap2 != NULL) {
        *cap_getAdjacencyP(*cap2) = NULL;
        *cap2 = NULL;
    }
}

void cap_check(Cap *cap) {
    End *end = cap_getEnd(cap);
    cactusCheck(end_getInstance(end, cap_getName(cap)) == cap);
    cactusCheck(cap_getOrientation(cap) == end_getOrientation(end));
    cactusCheck(end_getSide(end) == cap_getSide(cap)); //This is critical, it ensures
    //that we have a consistently oriented set of caps in an end.

    //If stub end checks, there is no attached segment.
    if (end_isStubEnd(end)) {
        cactusCheck(cap_getSegment(cap) == NULL);
    } else {
        Segment *segment = cap_getSegment(cap);
        if (segment != NULL) {
            cactusCheck(cap_getOrientation(cap) == segment_getOrientation(segment));
        }
    }

    //Checks adjacencies are properly linked and have consistent coordinates and the same group
    Cap *cap2 = cap_getAdjacency(cap);
    if (cap2 != NULL) {
        cactusCheck(end_getGroup(cap_getEnd(cap2)) == end_getGroup(end)); //check they have the same group.
        cactusCheck(cap_getAdjacency(cap2) == cap); //reciprocal connection
        cactusCheck(cap_getEvent(cap) == cap_getEvent(cap2)); //common event
        cactusCheck(cap_getStrand(cap) == cap_getStrand(cap2)); //common strand
        cactusCheck(cap_getSequence(cap) == cap_getSequence(cap2)); //common sequence (which may be null)

        if (cap_getCoordinate(cap) != INT64_MAX) { //if they have a coordinate
            cactusCheck(cap_getSide(cap) != cap_getSide(cap2)); //they have to represent an interval
            if (cap_getStrand(cap)) {
                if (!cap_getSide(cap)) {
                    cactusCheck(cap_getCoordinate(cap) < cap_getCoordinate(cap2));
                } else {
                    cactusCheck(cap_getCoordinate(cap) > cap_getCoordinate(cap2));
                }
            } else {
                if (cap_getSide(cap)) {
                    cactusCheck(cap_getCoordinate(cap) < cap_getCoordinate(cap2));
                } else {
                    cactusCheck(cap_getCoordinate(cap) > cap_getCoordinate(cap2));
                }
            }
        } else {
            cactusCheck(cap_getCoordinate(cap2) == INT64_MAX);
        }
    }

    //Checks the reverse
    Cap *rCap = cap_getReverse(cap);
    cactusCheck(rCap != NULL);
    cactusCheck(cap_getReverse(rCap) == cap);
    cactusCheck(cap_getOrientation(cap) == !cap_getOrientation(rCap));
    cactusCheck(cap_getEnd(cap) == end_getReverse(cap_getEnd(rCap)));
    cactusCheck(cap_getName(cap) == cap_getName(rCap));
    cactusCheck(cap_getEvent(cap) == cap_getEvent(rCap));
    cactusCheck(cap_getEnd(cap) == end_getReverse(cap_getEnd(rCap)));
    if (cap_getSegment(cap) == NULL) {
        cactusCheck(cap_getSegment(rCap) == NULL);
    } else {
        cactusCheck(cap_getSegment(cap) == segment_getReverse(cap_getSegment(rCap)));
    }
    cactusCheck(cap_getSide(cap) == !cap_getSide(rCap));
    cactusCheck(cap_getCoordinate(cap) == cap_getCoordinate(rCap));
    cactusCheck(cap_getSequence(cap) == cap_getSequence(rCap));
    cactusCheck(cap_getStrand(cap) == !cap_getStrand(rCap));
    if (cap_getAdjacency(cap) == NULL) {
        cactusCheck(cap_getAdjacency(rCap) == NULL);
    } else {
        cactusCheck(cap_getReverse(cap_getAdjacency(rCap)) == cap_getAdjacency(cap));
    }

    //it is consistent with any copy of end in the nested flower, in terms of events, connections and coordinates.
    Flower *nestedFlower = group_getNestedFlower(end_getGroup(end));
    if (nestedFlower != NULL) {
        End *childEnd = flower_getEnd(nestedFlower, end_getName(end));
        cactusCheck(childEnd != NULL); //End must be present in child
        //Cap *childCap = end_getInstance(childEnd, cap_getName(cap));
    }
}


