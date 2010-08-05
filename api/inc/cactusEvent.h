#ifndef CACTUS_EVENT_H_
#define CACTUS_EVENT_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an event, attached on a branch from the parent event, with the given branch length.
 */
Event *event_construct(MetaEvent *metaEvent, float branchLength, Event *parentEvent, EventTree *eventTree);

/*
 * Constructs an event, on the branch between the given child and parent events. The branch length is the
 * length from the parent to the new event. In this case it should be less than the length from the parent to
 * the pre-existing child event, however if it is not, then the child branch length
 * will be set to 0 and the new branch will extend the total length of the combined branch.
 */
Event *event_construct2(MetaEvent *metaEvent, float branchLength,
		Event *parentEvent, Event *childEvent, EventTree *eventTree);

/*
 * Returns the parent event, or NULL if root.
 */
Event *event_getParent(Event *event);

/*
 * Gets the name of the event.
 */
Name event_getName(Event *event);

/*
 * Gets the associated meta event.
 */
MetaEvent *event_getMetaEvent(Event *event);

/*
 * Gets the header sequence associated with the event.
 */
const char *event_getHeader(Event *event);

/*
 * Gets the branch length.
 */
float event_getBranchLength(Event *event);

/*
 * Gets the branch length of the subtree rooted at this event, excluding the branch length of the event
 * itself.
 */
float event_getSubTreeBranchLength(Event *event);

/*
 * Gets the number of events in the sub tree of the event, excluding the event itself.
 */
int32_t event_getSubTreeEventNumber(Event *event);

/*
 * Get number of children.
 */
int32_t event_getChildNumber(Event *event);

/*
 * Gets a child by its index.
 */
Event *event_getChild(Event *event, int32_t index);

/*
 * Gets the event tree the event is part of.
 */
EventTree *event_getEventTree(Event *event);

/*
 * Returns non zero if the other event is an ancestor of the event.
 */
int32_t event_isAncestor(Event *event, Event *otherEvent);

/*
 * Returns non zero if the other event is a child of the event.
 */
int32_t event_isDescendant(Event *event, Event *otherEvent);

/*
 * Returns non zero if the other event is a sibling of the other event.
 */
int32_t event_isSibling(Event *event, Event *otherEvent);

/*
 * Checks the following:
 * Event has parent, unless it is root.
 * Each child event has event as parent.
 * Ancestor-event --> event edge is consistent with any event tree that is in the parent of the containing flower.
 */
void event_check(Event *event);

#endif
