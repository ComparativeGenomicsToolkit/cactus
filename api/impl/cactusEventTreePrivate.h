/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_EVENT_TREE_PRIVATE_H_
#define CACTUS_EVENT_TREE_PRIVATE_H_

#include "cactusGlobals.h"

struct _eventTree {
	Event *rootEvent;
	stSortedSet *events;
	Flower *flower;
};


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private event tree functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs the event tree and all its events.
 */
void eventTree_destruct(EventTree *eventTree);

/*
 * Adds the cap to the event tree.
 */
void eventTree_addEvent(EventTree *eventTree, Event *event);

/*
 * Removes the instance from the event tree.
 */
void eventTree_removeEvent(EventTree *eventTree, Event *event);

/*
 * Creates a binary representation of the eventTree, returned as a char string.
 */
void eventTree_writeBinaryRepresentation(EventTree *eventTree, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a eventTree into memory from a binary representation of the eventTree.
 */
EventTree *eventTree_loadFromBinaryRepresentation(void **binaryString, Flower *flower);

#endif
