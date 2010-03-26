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
Cap *cap_construct2(End *end, int32_t startCoordinate, bool strand, bool side, Sequence *sequence);

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
 * Gets the coordinate of the position that cap is on the end of,
 * returns INT32_MAX if coordinate not set.
 */
int32_t cap_getCoordinate(Cap *cap);

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
void cap_makeAdjacent1(Cap *cap, Cap *cap2);

/*
 * Sets alternatively adjacent caps (this will set the adjacency reciprocally).
 * Any previous alternative adjacency will be set to NULL for both ends.
 */
void cap_makeAdjacent2(Cap *cap, Cap *cap2);

/*
 * Sets any adjacent ends for the first adjacency of the instance to NULL;
 */
void cap_breakAdjacency1(Cap *cap);

/*
 * Gets the adjacent cap.
 */
Cap *cap_getAdjacency(Cap *cap);

/*
 * Gets the alternative adjacency.
 */
Cap *cap_getAdjacency2(Cap *cap);

/*
 * Gets any face associated the cap.
 */
Face *cap_getFace(Cap *cap);

/*
 * Gets the parent cap (in the tree of the end).
 */
Cap *cap_getParent(Cap *cap);

/*
 * Returns the number of children the cap has.
 */
int32_t cap_getChildNumber(Cap *cap);

/*
 * Gets the child cap in the tree of the end.
 */
Cap *cap_getChild(Cap *cap, int32_t index);

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
int32_t cap_isInternal(Cap *cap);

/*
 * Returns a non zero integer if the cap is augmented (added to accommodate an adjacency, but without an attached segment).
 */
int32_t cap_isAugmented(Cap *cap);


#endif
