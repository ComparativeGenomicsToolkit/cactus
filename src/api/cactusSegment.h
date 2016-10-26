/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_ATOM_INSTANCE_H_
#define CACTUS_ATOM_INSTANCE_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic segment functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an segment and its two attached caps.
 */
Segment *segment_construct(Block *block, Event *event);

/*
 * Constructs an segment and its two attached caps, with the given coordinates.
 * Coordinates are always with respect to the positive strand. Thus the instance will
 * cover startCoordinate (inclusive) and startCoordinate + length (exclusive).
 * The instance returned by the function will be on the given strand.
 * If the strand is negative then the instance returned will have a start coordinate of
 * startCoordinate + length - 1.
 */
Segment *segment_construct2(Block *block,
		int64_t startCoordinate, bool strand, Sequence *sequence);

/*
 * Gets the encompassing block.
 */
Block *segment_getBlock(Segment *segment);

/*
 * Returns instances name.
 */
Name segment_getName(Segment *segment);

/*
 * Returns a non zero integer if the instance is oriented positively with respect to the block.
 * Else, returns zero if the instance is oriented negatively with respect to the block.
 */
bool segment_getOrientation(Segment *segment);

/*
 * Gets a positively oriented segment.
 */
Segment *segment_getPositiveOrientation(Segment *segment);

/*
 * Gets the reverse segment, giving a reversed view of the segment.
 */
Segment *segment_getReverse(Segment *segment);

/*
 * Gets the event associated with the instance.
 */
Event *segment_getEvent(Segment *segment);

/*
 * Gets the start coordinate (that which is closest to the 5' end of the strand)
 *  of the segment, returns INT64_MAX if coordinate not set.
 */
int64_t segment_getStart(Segment *segment);

/*
 * Returns non zero if one the forward strand, and zero if on the minus strand.
 */
bool segment_getStrand(Segment *segment);

/*
 * Gets the length of the segment.
 */
int64_t segment_getLength(Segment *segment);

/*
 * Gets the sequence in which the instance exists, or NULL if not set.
 */
Sequence *segment_getSequence(Segment *segment);

/*
 * Gets the actual sequence string of the segment (if the coordinates set), or returns NULL;
 */
char *segment_getString(Segment *segment);

/*
 * Gets the left cap of the segment.
 */
Cap *segment_get5Cap(Segment *segment);

/*
 * Gets the right cap of the segment.
 */
Cap *segment_get3Cap(Segment *segment);

/*
 * Returns the parent instance, or NULL, if none exists.
 */
Segment *segment_getParent(Segment *segment);

/*
 * Returns the number of children the instance has.
 */
int64_t segment_getChildNumber(Segment *segment);

/*
 * Gets a child instance.
 */
Segment *segment_getChild(Segment *segment, int64_t index);

/*
 * Links together a parent and child segment.
 */
void segment_makeParentAndChild(Segment *segmentParent, Segment *segmentChild);

/*
 * Checks (amongst other things) the following:
 * Checks the two ends have caps.
 * Checks the coordinates of the caps are consistent with the segment.
 * Checks the the segment has a parent, unless the root.
 * Checks the reverse segment is mirror of the segment.
 */
void segment_check(Segment *segment);

#endif
