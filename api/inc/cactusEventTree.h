/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_EVENT_TREE_H_
#define CACTUS_EVENT_TREE_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event tree functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an event tree, with given root event.
 */
EventTree *eventTree_construct(CactusDisk *cactusDisk, Name rootEventName);

/*
 * Constructs an event tree, with one root event.
 */
EventTree *eventTree_construct2(CactusDisk *cactusDisk);

/*
 * Copy constructs the event tree. Only includes
 * the unary events which are true, given the unary event function.
 */
EventTree *eventTree_copyConstruct(EventTree *eventTree, int64_t (unaryEventFilterFn)(Event *event));

/*
 * Returns the root event.
 */
Event *eventTree_getRootEvent(EventTree *eventTree);

/*
 * Gets the event with the given name.
 */
Event *eventTree_getEvent(EventTree *eventTree, Name eventName);

/*
 * Finds an event in the given event tree with the given header string.
 */
Event *eventTree_getEventByHeader(EventTree *eventTree, const char *eventHeader);

/*
 * Gets the common ancestor of two events.
 */
Event *eventTree_getCommonAncestor(Event *event, Event *event2);

/*
 * Gets the total number of events in the event tree.
 */
int64_t eventTree_getEventNumber(EventTree *eventTree);

/*
 * Gets the first event in the list.
 */
Event *eventTree_getFirst(EventTree *eventTree);

/*
 * Gets an iterator over the eventTree events.
 */
EventTree_Iterator *eventTree_getIterator(EventTree *eventTree);

/*
 * Gets the next event from the iterator.
 */
Event *eventTree_getNext(EventTree_Iterator *iterator);

/*
 * Gets the previous event from the iterator.
 */
Event *eventTree_getPrevious(EventTree_Iterator *iterator);

/*
 * Duplicates the iterator.
 */
EventTree_Iterator *eventTree_copyIterator(EventTree_Iterator *iterator);

/*
 * Destructs the iterator.
 */
void eventTree_destructIterator(EventTree_Iterator *iterator);

/*
 * Makes a newick string representation of the event tree.
 */
char *eventTree_makeNewickString(EventTree *eventTree);

/*
 * Runs event_check for each event in tree, which is sufficient to check properties of event tree.
 */
void eventTree_check(EventTree *eventTree);

/*
 * Get an stTree labeled by the event Names (unlike makeNewickString,
 * which labels the newick tree by headers.)
 */
stTree *eventTree_getStTree(EventTree *eventTree);

/*
 * Get the CactusDisk that owns this event tree.
 */
CactusDisk *eventTree_getCactusDisk(EventTree *eventTree);

#endif
