/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_END_INSTANCE_H_
#define CACTUS_END_INSTANCE_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic cap functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an cap, but not its connecting objects. Instance is the suffix m of the instance name n.m.
 */
Cap *cap_construct(End *end, Event *event);

/*
 * As default constructor, but also sets the instance's coordinates and event.
 */
Cap *cap_construct2(End *end, int64_t startCoordinate, bool strand, Sequence *sequence);

/*
 * Adds coordinates to a given cap.
 */
void cap_setCoordinates(Cap *cap, int64_t coordinate, bool strand, Sequence *sequence);

/*
 * Adds the given cap into the end - will fail if the given cap is not a member of
 * an end with the same name as end.
 */
Cap *cap_copyConstruct(End *end, Cap *cap);

/*
 * Returns the instance's name.
 */
Name cap_getName(Cap *cap);

/*
 * Returns a non zero integer if the cap is oriented positively with respect to the end.
 * Else, returns zero if the cap is oriented negatively with respect to the end.
 */
bool cap_getOrientation(Cap *cap);

/*
 * Gets the positively oriented cap.
 */
Cap *cap_getPositiveOrientation(Cap *cap);

/*
 * Gets the reversed cap (the equivalent on the opposite strand, with the opposite orientation).
 */
Cap *cap_getReverse(Cap *cap);

/*
 * Gets the event associated with the cap.
 */
Event *cap_getEvent(Cap *cap);

/*
 * Gets the encompassing end.
 */
End *cap_getEnd(Cap *cap);

/*
 * Gets the segment associated with the end, or NULL, if the end has no associated block end at this level.
 */
Segment *cap_getSegment(Cap *cap);

/*
 * If the cap is part of a segment, gets the other cap connected to the end of segment,
 * else it returns NULL;
 */
Cap *cap_getOtherSegmentCap(Cap *cap);

/*
 * Gets the coordinate of the position that cap is on the end of,
 * returns INT64_MAX if coordinate not set.
 */
int64_t cap_getCoordinate(Cap *cap);

/*
 * Returns a non zero integer if the coordinate of the cap (see cap_getCoordinate)
 * is on the forward strand, and zero if on the negative strand. The return value is undefined if the coordinates are not set.
 */
bool cap_getStrand(Cap *cap);

/*
 * Returns a non zero integer if on the 5' side of the position returned by cap_getCoordinate,
 * zero if on the 3' side.
 */
bool cap_getSide(Cap *cap);

/*
 * Gets the sequence in which the instance exists, or NULL if not set.
 */
Sequence *cap_getSequence(Cap *cap);

/*
 * Sets adjacent caps (this will set the adjacency reciprocally).
 * Any previous adjacency will be set to NULL for both ends.
 */
void cap_makeAdjacent(Cap *cap, Cap *cap2);

/*
 * Sets any adjacent ends for the first adjacency of the instance to NULL;
 */
void cap_breakAdjacency(Cap *cap);

/*
 * Gets the adjacent cap.
 */
Cap *cap_getAdjacency(Cap *cap);

/*
 * Returns the top node associated with the cap. Returns NULL if
 * the cap is not attached or if the cap is a root node. This function is
 * denoted A(c) for cap c in the AVG paper.
 */
Cap *cap_getTopCap(Cap *cap);

/*
 * Gets any face in which the cap is a top node.
 */
Face *cap_getTopFace(Cap *cap);

/*
 * Returns any face-end associated in which the cap is a top node, or NULL.
 */
FaceEnd *cap_getTopFaceEnd(Cap *cap);

/*
 * Returns any face-end associated in the cap is a bottom node, or NULL.
 */
FaceEnd *cap_getBottomFaceEnd(Cap *cap);

/*
 * Gets the parent cap (in the tree of the end).
 */
Cap *cap_getParent(Cap *cap);

/*
 * Returns the number of children the cap has.
 */
int64_t cap_getChildNumber(Cap *cap);

/*
 * Gets the child cap in the tree of the end.
 */
Cap *cap_getChild(Cap *cap, int64_t index);

/*
 * Links together a parent and child cap.
 */
void cap_makeParentAndChild(Cap *capParent, Cap *capChild);

/*
 * Changes the parent cap of a child cap
 */
void cap_changeParentAndChild(Cap* newCapParent, Cap* capChild);

/*
 * Returns non zero if the cap is internal (part of an internal tree).
 */
bool cap_isInternal(Cap *cap);

/*
 * Checks the following (amongst other things):
 * If end has tree:
 *  (1) checks the cap has a parent which has an ancestral event to the caps event, unless it is the root.
 *  (2) Check caps ancestor/descendant links are proper.
 * If stub end checks, there is no attached segment.
 * Checks adjacencies are properly linked and have consistent coordinates and the same group.
 * It is consistent with any copy of end in the nested flower, in terms of connections and coordinates.
 * The reverse cap is the mirror of the cap.
 */
void cap_check(Cap *cap);


#endif
