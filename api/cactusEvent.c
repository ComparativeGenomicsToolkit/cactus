#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Event *event_construct(MetaEvent *metaEvent, float branchLength,
		Event *parentEvent, EventTree *eventTree) {
	assert(eventTree_getEvent(eventTree, metaEvent_getName(metaEvent)) == NULL); //the event must not already exist in the tree.
	Event *event;
	event = st_malloc(sizeof(Event));
	event->metaEvent = metaEvent;
	event->parent = parentEvent;
	event->children = constructEmptyList(0, NULL);
	event->branchLength = branchLength < 0.0 ? 0.0: branchLength;
	if(parentEvent != NULL) {
		listAppend(parentEvent->children, event);
	}
	event->eventTree = eventTree;
	eventTree_addEvent(eventTree, event);
	return event;
}

Event *event_construct2(MetaEvent *metaEvent, float branchLength,
		Event *parentEvent, Event *childEvent, EventTree *eventTree) {
	Event *event;
	event = event_construct(metaEvent, branchLength, parentEvent, eventTree);
#ifdef BEN_DEBUG
	assert(parentEvent != NULL);
	assert(childEvent != NULL);
	assert(listContains(parentEvent->children, childEvent));
#endif
	listRemove(parentEvent->children, childEvent);
	listAppend(event->children, childEvent);
	childEvent->parent = event;
	childEvent->branchLength = childEvent->branchLength - event->branchLength;
	if(childEvent->branchLength < 0.0) {
		childEvent->branchLength = 0.0;
	}
	return event;
}

Event *event_getParent(Event *event) {
	assert(event != NULL);
	return event->parent;
}

Name event_getName(Event *event) {
	assert(event != NULL);
	return metaEvent_getName(event->metaEvent);
}

MetaEvent *event_getMetaEvent(Event *event) {
	assert(event != NULL);
	return event->metaEvent;
}

const char *event_getHeader(Event *event) {
	assert(event != NULL);
	return metaEvent_getHeader(event->metaEvent);
}

float event_getBranchLength(Event *event) {
	assert(event != NULL);
	return event->branchLength;
}

float event_getSubTreeBranchLength(Event *event) {
	assert(event != NULL);
	int32_t i;
	Event *childEvent;
	float branchLength;

	branchLength = 0.0;
	for(i=0; i<event_getChildNumber(event); i++) {
		childEvent = event_getChild(event, i);
		branchLength += event_getSubTreeBranchLength(childEvent) + event_getBranchLength(childEvent);
	}
	return branchLength;
}

int32_t event_getSubTreeEventNumber(Event *event) {
	assert(event != NULL);
	int32_t i, j;
	Event *childEvent;

	j = 0.0;
	for(i=0; i<event_getChildNumber(event); i++) {
		childEvent = event_getChild(event, i);
		j += event_getSubTreeEventNumber(childEvent) + 1;
	}
	return j;
}

int32_t event_getChildNumber(Event *event) {
	assert(event != NULL);
	return event->children->length;
}

Event *event_getChild(Event *event, int32_t index) {
#ifdef BEN_DEBUG
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

int32_t event_isAncestor(Event *event, Event *otherEvent) {
	return event != otherEvent && eventTree_getCommonAncestor(event, otherEvent) == otherEvent;
}

int32_t event_isDescendant(Event *event, Event *otherEvent) {
	return event != otherEvent && eventTree_getCommonAncestor(event, otherEvent) == event;
}

int32_t event_isSibling(Event *event, Event *otherEvent) {
	Event *ancestorEvent;
	ancestorEvent = eventTree_getCommonAncestor(event, otherEvent);
	return event != otherEvent && ancestorEvent != event && ancestorEvent != otherEvent;
}

void event_check(Event *event) {
	EventTree *eventTree = event_getEventTree(event);
	Event *ancestorEvent = event_getParent(event);

	//Check event and eventree properly linked
	assert(eventTree_getEvent(event_getEventTree(event), event_getName(event)) == event);

	//Event has parent, unless it is root.
	if(eventTree_getRootEvent(eventTree) == event) {
		assert(ancestorEvent == NULL);
	}
	else { //not root, so must have ancestor.
		assert(ancestorEvent != NULL);
	}

	//Each child event has event as parent.
	int32_t i=0;
	for(i=0; i<event_getChildNumber(event); i++) {
		Event *childEvent = event_getChild(event, i);
		assert(event_getParent(childEvent) == event);
	}

	//Ancestor-event --> event edge is consistent with any event tree that is in the parent of the containing net.
	Group *parentGroup = net_getParentGroup(eventTree_getNet(event_getEventTree(event)));
	if(parentGroup != NULL) {
		EventTree *parentEventTree = net_getEventTree(group_getNet(parentGroup));
		Event *parentEvent = eventTree_getEvent(parentEventTree, event_getName(event));
		if(parentEvent != NULL) {
			if(ancestorEvent == NULL) { //the case where they are both root.
				assert(eventTree_getRootEvent(parentEventTree) == parentEvent);
			}
			else {
				//Check edge ancestorEvent --> event is in parent event tree.
				while(1) {
					Event *parentAncestorEvent = eventTree_getEvent(parentEventTree, event_getName(ancestorEvent));
					if(parentAncestorEvent != NULL) {
						assert(event_isAncestor(parentEvent, parentAncestorEvent));
						break;
					}
					ancestorEvent = event_getParent(ancestorEvent);
					assert(ancestorEvent != NULL);
				}
			}
		}
	}
}

/*
 * Private functions
 */

void event_destruct(Event *event) {
	int32_t i;
	Event *childEvent;
	Event *parentEvent = event_getParent(event);
	if(parentEvent != NULL) {
		listRemove(parentEvent->children, event);
	}
	eventTree_removeEvent(event_getEventTree(event), event);
	for(i=0; i<event->children->length; i++) {
		childEvent = event->children->list[i];
		childEvent->parent = parentEvent;
		if(parentEvent != NULL) {
			listAppend(parentEvent->children, childEvent);
		}
	}
	destructList(event->children);
	free(event);
}

Event *event_getStaticNameWrapper(Name name) {
	static Event event;
	static MetaEvent metaEvent;
	metaEvent.name = name;
	event.metaEvent = &metaEvent;
	return &event;
}

/*
 * Serialisation functions
 */

void event_writeBinaryRepresentation(Event *event,
		void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_EVENT, writeFn);
	binaryRepresentation_writeName(event_getName(event_getParent(event)), writeFn);
	binaryRepresentation_writeName(event_getName(event), writeFn);
	binaryRepresentation_writeFloat(event_getBranchLength(event), writeFn);
}

Event *event_loadFromBinaryRepresentation(void **binaryString,
		EventTree *eventTree) {
	Event *event, *parentEvent;
	MetaEvent *metaEvent;
	Name name;
	float branchLength;

	event = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_EVENT) {
		binaryRepresentation_popNextElementType(binaryString);
		parentEvent = eventTree_getEvent(eventTree, binaryRepresentation_getName(binaryString));
		assert(parentEvent != NULL);
		name = binaryRepresentation_getName(binaryString);
		branchLength = binaryRepresentation_getFloat(binaryString);
		metaEvent = cactusDisk_getMetaEvent(net_getNetDisk(eventTree_getNet(eventTree)),name);
		assert(metaEvent != NULL);
		event = event_construct(metaEvent, branchLength, parentEvent, eventTree);
	}
	return event;
}

