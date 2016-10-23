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
    return event_construct(cactusDisk_getUniqueID(flower_getCactusDisk(
            eventTree_getFlower(eventTree))), header, branchLength, parentEvent,
            eventTree);
}

Event *event_construct4(const char *header, float branchLength, Event *parentEvent,
        Event *childEvent, EventTree *eventTree) {
    return event_construct2(cactusDisk_getUniqueID(flower_getCactusDisk(
                eventTree_getFlower(eventTree))), header, branchLength, parentEvent, childEvent,
                eventTree);
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

    //Ancestor-event --> event edge is consistent with any event tree that is in the parent of the containing flower.
    Group *parentGroup = flower_getParentGroup(eventTree_getFlower(
            event_getEventTree(event)));
    if (parentGroup != NULL) {
        EventTree *parentEventTree = flower_getEventTree(group_getFlower(
                parentGroup));
        Event *parentEvent = eventTree_getEvent(parentEventTree, event_getName(
                event));
        if (parentEvent != NULL) {
            if (ancestorEvent == NULL) { //the case where they are both root.
                cactusCheck(eventTree_getRootEvent(parentEventTree) == parentEvent);
            } else {
                //Check edge ancestorEvent --> event is in parent event tree.
                while (1) {
                    Event *parentAncestorEvent = eventTree_getEvent(
                            parentEventTree, event_getName(ancestorEvent));
                    if (parentAncestorEvent != NULL) {
                        cactusCheck(event_isAncestor(parentEvent, parentAncestorEvent));
                        break;
                    }
                    ancestorEvent = event_getParent(ancestorEvent);
                    cactusCheck(ancestorEvent != NULL);
                }
            }
        }
    }
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

/*
 * Serialisation functions
 */

void event_writeBinaryRepresentation(Event *event, void(*writeFn)(
        const void * ptr, size_t size, size_t count)) {
    binaryRepresentation_writeElementType(CODE_EVENT, writeFn);
    binaryRepresentation_writeName(event_getName(event_getParent(event)),
            writeFn);
    binaryRepresentation_writeName(event_getName(event), writeFn);
    binaryRepresentation_writeFloat(event_getBranchLength(event), writeFn);
    binaryRepresentation_writeString(event_getHeader(event), writeFn);
    binaryRepresentation_writeBool(event_isOutgroup(event), writeFn);
}

Event *event_loadFromBinaryRepresentation(void **binaryString,
        EventTree *eventTree) {
    Event *event, *parentEvent;
    Name name;
    float branchLength;
    char *header;

    event = NULL;
    if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_EVENT) {
        binaryRepresentation_popNextElementType(binaryString);
        parentEvent = eventTree_getEvent(eventTree,
                binaryRepresentation_getName(binaryString));
        assert(parentEvent != NULL);
        name = binaryRepresentation_getName(binaryString);
        branchLength = binaryRepresentation_getFloat(binaryString);
        header = binaryRepresentation_getString(binaryString);
        event
                = event_construct(name, header, branchLength, parentEvent,
                        eventTree);
        event_setOutgroupStatus(event, binaryRepresentation_getBool(binaryString));
        free(header);
    }
    return event;
}

