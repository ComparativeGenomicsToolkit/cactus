/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Event *event_construct(Name name, const char *header, float branchLength, Event *parentEvent,
        EventTree *eventTree) {
    assert(eventTree_getEvent(eventTree, name) == NULL); //the event must not already exist in the tree.
    Event *event;
    event = st_malloc(sizeof(Event));
    event->name = name;
    event->parent = parentEvent;
    event->children = constructEmptyList(0, NULL);
    event->header = stString_copy(header == NULL ? "" : header);
    event->branchLength = branchLength < 0.0 ? 0.0 : branchLength;
    event->isOutgroup = 0;
    if (parentEvent != NULL) {
        listAppend(parentEvent->children, event);
    }
    event->eventTree = eventTree;
    eventTree_addEvent(eventTree, event);
    return event;
}

Event *event_construct2(Name name, const char *header, float branchLength, Event *parentEvent,
        Event *childEvent, EventTree *eventTree) {
    Event *event;
    event = event_construct(name, header, branchLength, parentEvent, eventTree);
#ifndef NDEBUG
    assert(parentEvent != NULL);
    assert(childEvent != NULL);
    assert(listContains(parentEvent->children, childEvent));
#endif
    listRemove(parentEvent->children, childEvent);
    listAppend(event->children, childEvent);
    childEvent->parent = event;
    childEvent->branchLength = childEvent->branchLength - event->branchLength;
    if (childEvent->branchLength < 0.0) {
        childEvent->branchLength = 0.0;
    }
    return event;
}

Event *event_construct3(const char *header, float branchLength, Event *parentEvent,
        EventTree *eventTree) {
    return event_construct(cactusDisk_getUniqueID(eventTree_getCactusDisk(eventTree)),
                           header, branchLength, parentEvent, eventTree);
}

Event *event_construct4(const char *header, float branchLength, Event *parentEvent,
        Event *childEvent, EventTree *eventTree) {
    return event_construct2(cactusDisk_getUniqueID(eventTree_getCactusDisk(eventTree)),
                            header, branchLength, parentEvent, childEvent, eventTree);
}

Event *event_getParent(Event *event) {
    assert(event != NULL);
    return event->parent;
}

Name event_getName(Event *event) {
    assert(event != NULL);
    return event->name;
}

float event_getBranchLength(Event *event) {
    assert(event != NULL);
    return event->branchLength;
}

const char *event_getHeader(Event *event) {
    assert(event != NULL);
    return event->header;
}

float event_getSubTreeBranchLength(Event *event) {
    assert(event != NULL);
    int64_t i;
    Event *childEvent;
    float branchLength;

    branchLength = 0.0;
    for (i = 0; i < event_getChildNumber(event); i++) {
        childEvent = event_getChild(event, i);
        branchLength += event_getSubTreeBranchLength(childEvent)
                + event_getBranchLength(childEvent);
    }
    return branchLength;
}

int64_t event_getSubTreeEventNumber(Event *event) {
    assert(event != NULL);
    int64_t i, j;
    Event *childEvent;

    j = 0.0;
    for (i = 0; i < event_getChildNumber(event); i++) {
        childEvent = event_getChild(event, i);
        j += event_getSubTreeEventNumber(childEvent) + 1;
    }
    return j;
}

int64_t event_getChildNumber(Event *event) {
    assert(event != NULL);
    return event->children->length;
}

Event *event_getChild(Event *event, int64_t index) {
#ifndef NDEBUG
    assert(event != NULL);
    assert(index >= 0);
    assert(index < event_getChildNumber(event));
#endif
    return event->children->list[index];
}

EventTree *event_getEventTree(Event *event) {
    assert(event != NULL);
    return event->eventTree;
}

int64_t event_isAncestor(Event *event, Event *otherEvent) {
    return event != otherEvent
            && eventTree_getCommonAncestor(event, otherEvent) == otherEvent;
}

int64_t event_isDescendant(Event *event, Event *otherEvent) {
    return event != otherEvent
            && eventTree_getCommonAncestor(event, otherEvent) == event;
}

bool event_isOutgroup(Event *event) {
    return event->isOutgroup;
}

void event_setOutgroupStatus(Event *event, bool isOutgroup) {
    event->isOutgroup = isOutgroup;
}

int64_t event_isSibling(Event *event, Event *otherEvent) {
    Event *ancestorEvent;
    ancestorEvent = eventTree_getCommonAncestor(event, otherEvent);
    return event != otherEvent && ancestorEvent != event && ancestorEvent
            != otherEvent;
}

void event_check(Event *event) {
    EventTree *eventTree = event_getEventTree(event);
    Event *ancestorEvent = event_getParent(event);

    //Check event and eventree properly linked
    cactusCheck(eventTree_getEvent(event_getEventTree(event), event_getName(event)) == event);

    //Event has parent, unless it is root.
    if (eventTree_getRootEvent(eventTree) == event) {
        cactusCheck(ancestorEvent == NULL);
    } else { //not root, so must have ancestor.
        cactusCheck(ancestorEvent != NULL);
    }

    //Each child event has event as parent.
    int64_t i = 0;
    for (i = 0; i < event_getChildNumber(event); i++) {
        Event *childEvent = event_getChild(event, i);
        cactusCheck(event_getParent(childEvent) == event);
    }
}

stTree *event_getStTree(Event *event) {
    stTree *ret = stTree_construct();
    stTree_setLabel(ret, stString_print("%" PRIi64, event_getName(event)));
    stTree_setBranchLength(ret, event_getBranchLength(event));
    for(int64_t i = 0; i < event_getChildNumber(event); i++) {
        Event *child = event_getChild(event, i);
        stTree *childStTree = event_getStTree(child);
        stTree_setParent(childStTree, ret);
    }
    return ret;
}

/*
 * Private functions
 */

void event_destruct(Event *event) {
    int64_t i;
    Event *childEvent;
    Event *parentEvent = event_getParent(event);
    if (parentEvent != NULL) {
        listRemove(parentEvent->children, event);
    }
    eventTree_removeEvent(event_getEventTree(event), event);
    for (i = 0; i < event->children->length; i++) {
        childEvent = event->children->list[i];
        childEvent->parent = parentEvent;
        if (parentEvent != NULL) {
            listAppend(parentEvent->children, childEvent);
        }
    }
    destructList(event->children);
    free(event->header);
    free(event);
}

