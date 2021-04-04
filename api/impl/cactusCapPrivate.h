/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_END_INSTANCE_PRIVATE_H_
#define CACTUS_END_INSTANCE_PRIVATE_H_

#include "cactusGlobals.h"

typedef struct _capCoreContents {
    Name instance;
    int64_t coordinate;
    void *eventOrSequence;
} CapCoreContents;

typedef struct _capContents {
    CapCoreContents core;
    Cap *adjacency;
    End *end;
    Cap *nCap; // Links together different caps in the end
} CapContents;

typedef struct _segmentCapContents {
    CapCoreContents core;
    Cap *leftAdjacency;
    Cap *rightAdjacency;
    Block *block;
    Segment *nSegment; // Links together different segments in the block
} SegmentCapContents;

struct _cap { // Note, a segment has the same structure
    char bits;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//End instance functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * These functions deal with the internal structures
 * used to represent caps and caps+segments
 */
bool cap_forward(Cap *cap);
bool cap_partOfSegment(Cap *cap);
bool cap_isSegment(Cap *cap);
bool cap_left(Cap *cap);
bool cap_getHasEventNotSequence(Cap *cap);
SegmentCapContents *cap_getSegmentContents(Cap *cap);
SegmentCapContents *segment_getContents(Segment *segment);
CapContents *cap_getContents(Cap *cap);
CapCoreContents *cap_getCoreContents(Cap *cap);


/*
 * Constructs an cap, but not its connecting objects. Instance is the suffix m of the instance name n.m.
 */
Cap *cap_construct3(Name name, Event *event, End *end);

/*
 * As default constructor, but also sets the instance's coordinates and event.
 */
Cap *cap_construct4(Name name, End *end, int64_t startCoordinate,
        bool strand, Sequence *sequence);

/*
 * As constructor 3, but don't specify name.
 */
Cap *cap_construct5(Event *event, End *end);

/*
 * Destructs the cap, but not any connecting objects.
 */
void cap_destruct(Cap *cap);

#endif
