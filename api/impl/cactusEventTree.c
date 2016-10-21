/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event tree functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int eventTree_constructP(const void *o1, const void *o2) {
	return cactusMisc_nameCompare(event_getName((Event *)o1), event_getName((Event *)o2));
}

EventTree *eventTree_construct2(Flower *flower) {
	return eventTree_construct(cactusDisk_getUniqueID(flower_getCactusDisk(flower)), flower);
}

EventTree *eventTree_construct(Name rootEventName, Flower *flower) {
	EventTree *eventTree;
	eventTree = st_malloc(sizeof(EventTree));
	eventTree->events = stSortedSet_construct3(eventTree_constructP, NULL);
	eventTree->flower = flower;
	eventTree->rootEvent = event_construct(rootEventName, "ROOT", INT64_MAX, NULL, eventTree); //do this last as reciprocal call made to add the event to the events.
	flower_setEventTree(flower, eventTree);
	return eventTree;
}

void eventTree_copyConstructP(EventTree *eventTree, Event *event,
		int64_t (unaryEventFilterFn)(Event *event)) {
	int64_t i;
	Event *event2;
	for(i=0; i<event_getChildNumber(event); i++) {
		event2 = event_getChild(event, i);
		while(event_getChildNumber(event2) == 1 && unaryEventFilterFn != NULL && !unaryEventFilterFn(event2)) {
			//skip the event
			event2 = event_getChild(event2, 0);
		}
		event_setOutgroupStatus(event_construct(event_getName(event2), event_getHeader(event2), event_getBranchLength(event2),
						eventTree_getEvent(eventTree, event_getName(event)), eventTree), event_isOutgroup(event2));
		eventTree_copyConstructP(eventTree, event2, unaryEventFilterFn);
	}
}


EventTree *eventTree_copyConstruct(EventTree *eventTree, Flower *newFlower,
		int64_t (unaryEventFilterFn)(Event *event)) {
	EventTree *eventTree2;
	eventTree2 = eventTree_construct(event_getName(eventTree_getRootEvent(eventTree)), newFlower);
	eventTree_copyConstructP(eventTree2, eventTree_getRootEvent(eventTree), unaryEventFilterFn);
	return eventTree2;
}

Event *eventTree_getRootEvent(EventTree *eventTree) {
	return eventTree->rootEvent;
}

Event *eventTree_getEvent(EventTree *eventTree, Name eventName) {
	Event event;
	event.name = eventName;
	return stSortedSet_search(eventTree->events, &event);
}

Event *eventTree_getCommonAncestor(Event *event, Event *event2) {
	Event *ancestorEvent;
	struct List *list;

	assert(event != NULL);
	assert(event2 != NULL);
	assert(event_getEventTree(event) == event_getEventTree(event2));

	list = constructEmptyList(0, NULL);
	ancestorEvent = event;
	while(ancestorEvent != NULL) {
		if(ancestorEvent == event2) {
			destructList(list);
			return event2;
		}
		listAppend(list, ancestorEvent);
		ancestorEvent = event_getParent(ancestorEvent);
	}

	ancestorEvent = event2;
	while((ancestorEvent = event_getParent(ancestorEvent)) != NULL) {
		if(listContains(list, ancestorEvent)) {
			destructList(list);
			return ancestorEvent;
		}
	}
	destructList(list);
	assert(FALSE);
	return NULL;
}

Flower *eventTree_getFlower(EventTree *eventTree) {
	return eventTree->flower;
}

int64_t eventTree_getEventNumber(EventTree *eventTree) {
	return event_getSubTreeEventNumber(eventTree_getRootEvent(eventTree)) + 1;
}

Event *eventTree_getFirst(EventTree *eventTree) {
	return stSortedSet_getFirst(eventTree->events);
}

EventTree_Iterator *eventTree_getIterator(EventTree *eventTree) {
	return stSortedSet_getIterator(eventTree->events);
}

Event *eventTree_getNext(EventTree_Iterator *iterator) {
	return stSortedSet_getNext(iterator);
}

Event *eventTree_getPrevious(EventTree_Iterator *iterator) {
	return stSortedSet_getPrevious(iterator);
}

EventTree_Iterator *eventTree_copyIterator(EventTree_Iterator *iterator) {
	return stSortedSet_copyIterator(iterator);
}

void eventTree_destructIterator(EventTree_Iterator *iterator) {
	stSortedSet_destructIterator(iterator);
}

static char *eventTree_makeNewickStringP(Event *event) {
	int64_t i;
	char *cA = NULL;
	char *cA3;
	if(event_getChildNumber(event) > 0) {
		for(i=0;i<event_getChildNumber(event); i++) {
			char *cA2 = eventTree_makeNewickStringP(event_getChild(event, i));
			if(i > 0) {
				cA3 = st_malloc(sizeof(char)*(strlen(cA)+strlen(cA2)+2));
				sprintf(cA3, "%s,%s", cA, cA2);
				free(cA);
				cA = cA3;
			}
			else {
				cA = st_malloc(sizeof(char)*(strlen(cA2)+2));
				sprintf(cA, "(%s", cA2);
			}
			free(cA2);
		}
		cA3 = st_malloc(sizeof(char)*(strlen(cA) + strlen(event_getHeader(event)) + 30));
		sprintf(cA3, "%s)%s:%g", cA, event_getHeader(event), event_getBranchLength(event));
		free(cA);
		cA = cA3;
	}
	else {
		cA = st_malloc(sizeof(char)*(strlen(event_getHeader(event)) + 30));
		sprintf(cA, "%s:%g", event_getHeader(event), event_getBranchLength(event));
	}
	return cA;
}

char *eventTree_makeNewickString(EventTree *eventTree) {
	Event *rootEvent = eventTree_getRootEvent(eventTree);
	char *cA = eventTree_makeNewickStringP(rootEvent);
	char *cA2 = st_malloc(sizeof(char)*(strlen(cA) + 2));
	sprintf(cA2, "%s;", cA);
	free(cA);
	return cA2;
}

static int64_t eventTree_addSiblingUnaryEventP(Event *event, Event *event2) {
	/*
	 * Event is the new event, event2 event from the event tree we're adding to.
	 */
	assert(event != event2);
	Group *group1 = flower_getParentGroup(eventTree_getFlower(event_getEventTree(event)));
	Group *group2 = flower_getParentGroup(eventTree_getFlower(event_getEventTree(event2)));
	(void)group2;
	if(group1 != NULL) { //both events have a parent, so we can perhaps ask if one is the ancestor
		//of the other in the parent event tree.
		assert(group2 != NULL);
		Flower *parentFlower = group_getFlower(group1);
		assert(parentFlower == group_getFlower(group2));
		EventTree *parentEventTree = flower_getEventTree(parentFlower);
		Event *eventP = eventTree_getEvent(parentEventTree, event_getName(event)); //get the ancestral version of the event.
		Event *event2P = eventTree_getEvent(parentEventTree, event_getName(event2));
		if(eventP != NULL && event2P != NULL) { //we can answer who is truly ancestral because both are in the ancestral tree.
			assert(eventP != event2P);
			Event *event3 = eventTree_getCommonAncestor(eventP, event2P);
			assert(event3 == eventP || event3 == event2P); //one must be the ancestor of the other
			return event3 == eventP;
		}
	}
	else {
		assert(group2 == NULL); //they both must be root flowers.
	}
	//Maybe both events are in the sibling event tree, we can refer to that tree
	//to decide who is ancestral.
	EventTree *eventTree = event_getEventTree(event);
	Event *event2P = eventTree_getEvent(eventTree, event_getName(event2));
	if(event2P != NULL) { //event2 is in the sibling event tree, so we can decide who is ancestral.
		assert(event != event2P);
		Event *event3 = eventTree_getCommonAncestor(event, event2P);
		assert(event3 == event || event3 == event2P); //one must be the ancestor of the other
		return event3 == event;
	}

	//event2 is not in the parent or the sibling, so we should schedule it after
	//event, because the comparison might be valid for one event2's parent events..
	return 1;
}

void eventTree_addSiblingUnaryEvent(EventTree *eventTree, Event *event) {
	if(eventTree_getEvent(eventTree, event_getName(event)) == NULL) { //check it isn't already in there
		Event *event2 = event;
		do {
			assert(event_getChildNumber(event2) == 1);
			event2 = event_getChild(event2, 0);
		} while(eventTree_getEvent(eventTree, event_getName(event2)) == NULL);
		event2 = eventTree_getEvent(eventTree, event_getName(event2));
		assert(event2 != NULL);
		Event *event3 = event_getParent(event2);
		while(eventTree_addSiblingUnaryEventP(event, event3)) {
			event2 = event3;
			event3 = event_getParent(event2);
		}
		event_construct2(event_getName(event), event_getHeader(event), event_getBranchLength(event), event3, event2, eventTree);
	}
}

void eventTree_check(EventTree *eventTree) {
	//Check flower and event tree properly connected.
	cactusCheck(flower_getEventTree(eventTree_getFlower(eventTree)) == eventTree);

	Event *event;
	EventTree_Iterator *eventIterator = eventTree_getIterator(eventTree);
	while((event = eventTree_getNext(eventIterator)) != NULL) {
		event_check(event);
	}
	eventTree_destructIterator(eventIterator);
}

Event *eventTree_getEventByHeader(EventTree *eventTree, const char *eventHeader) {
    EventTree_Iterator *it = eventTree_getIterator(eventTree);
    Event *event;
    while ((event = eventTree_getNext(it)) != NULL) {
        if (strcmp(event_getHeader(event), eventHeader) == 0) {
            eventTree_destructIterator(it);
            return event;
        }
    }
    eventTree_destructIterator(it);
    return NULL;
}

/*
 * Private functions.
 */

void eventTree_destruct(EventTree *eventTree) {
	Event *event;
	flower_removeEventTree(eventTree_getFlower(eventTree), eventTree);
	while((event = eventTree_getFirst(eventTree)) != NULL) {
		event_destruct(event);
	}
	stSortedSet_destruct(eventTree->events);
	free(eventTree);
}

void eventTree_addEvent(EventTree *eventTree, Event *event) {
	stSortedSet_insert(eventTree->events, event);
}

void eventTree_removeEvent(EventTree *eventTree, Event *event) {
	stSortedSet_remove(eventTree->events, event);
}

/*
 * Serialisation functions
 */

void eventTree_writeBinaryRepresentationP(Event *event, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	int64_t i;
	event_writeBinaryRepresentation(event, writeFn);
	for(i=0; i<event_getChildNumber(event); i++) {
		eventTree_writeBinaryRepresentationP(event_getChild(event, i), writeFn);
	}
}

void eventTree_writeBinaryRepresentation(EventTree *eventTree, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	int64_t i;
	Event *event;
	event = eventTree_getRootEvent(eventTree);
	binaryRepresentation_writeElementType(CODE_EVENT_TREE, writeFn);
	binaryRepresentation_writeName(event_getName(event), writeFn);
	for(i=0; i<event_getChildNumber(event); i++) {
		eventTree_writeBinaryRepresentationP(event_getChild(event, i), writeFn);
	}
	binaryRepresentation_writeElementType(CODE_EVENT_TREE, writeFn);
}

EventTree *eventTree_loadFromBinaryRepresentation(void **binaryString, Flower *flower) {
	EventTree *eventTree;
	Name name;
	eventTree = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_EVENT_TREE) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		eventTree = eventTree_construct(name, flower);
		while(event_loadFromBinaryRepresentation(binaryString, eventTree) != NULL) {
			;
		}
		assert(binaryRepresentation_peekNextElementType(*binaryString) == CODE_EVENT_TREE);
		binaryRepresentation_popNextElementType(binaryString);
	}
	return eventTree;
}

static stTree *eventTree_getStTree_R(Event *event) {
    stTree *ret = stTree_construct();
    stTree_setLabel(ret, stString_print("%" PRIi64, event_getName(event)));
    stTree_setBranchLength(ret, event_getBranchLength(event));
    for(int64_t i = 0; i < event_getChildNumber(event); i++) {
        Event *child = event_getChild(event, i);
        stTree *childStTree = eventTree_getStTree_R(child);
        stTree_setParent(childStTree, ret);
    }
    return ret;
}

// Get species tree from event tree (labeled by the event Names),
// which requires ignoring the root event.
stTree *eventTree_getStTree(EventTree *eventTree) {
    Event *rootEvent = eventTree_getRootEvent(eventTree);
    // Need to skip the root event, since it is added onto the real
    // species tree.
    assert(event_getChildNumber(rootEvent) == 1);
    Event *speciesRoot = event_getChild(rootEvent, 0);
    return eventTree_getStTree_R(speciesRoot);
}
