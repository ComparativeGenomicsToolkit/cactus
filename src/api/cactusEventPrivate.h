/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_EVENT_PRIVATE_H_
#define CACTUS_EVENT_PRIVATE_H_

#include "cactusGlobals.h"

struct _event {
    Name name;
    char *header;
    struct List *children;
    float branchLength;
    Event *parent;
    EventTree *eventTree;
    bool isOutgroup;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs the event and also destructs any attached child events.
 */
void event_destruct(Event *event);

/*
 * Creates a binary representation of the event, returned as a char string.
 */
void event_writeBinaryRepresentation(Event *event, void(*writeFn)(
        const void * ptr, size_t size, size_t count));

/*
 * Loads a event into memory from a binary representation of the event.
 */
Event *event_loadFromBinaryRepresentation(void **binaryString,
        EventTree *eventTree);

#endif
