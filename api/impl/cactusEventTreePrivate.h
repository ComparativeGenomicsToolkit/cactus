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
    CactusDisk *cactusDisk;
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

#endif
