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

#endif
