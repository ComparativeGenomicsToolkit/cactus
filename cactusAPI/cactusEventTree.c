#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event tree functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t eventTree_constructP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(event_getName((Event *)o1), event_getName((Event *)o2));
}

EventTree *eventTree_construct(MetaEvent *rootMetaEvent, Net *net) {
	EventTree *eventTree;
	eventTree = malloc(sizeof(EventTree));
	eventTree->events = sortedSet_construct(eventTree_constructP);
	eventTree->net = net;
	eventTree->rootEvent = event_construct(rootMetaEvent, INT32_MAX, NULL, eventTree); //do this last as reciprocal call made to add the event to the events.
	net_addEventTree(net, eventTree);
	return eventTree;
}

void eventTree_copyConstructP(EventTree *eventTree, Event *event,
		int32_t (unaryEventFilterFn)(Event *event)) {
	int32_t i;
	Event *event2;
	for(i=0; i<event_getChildNumber(event); i++) {
		event2 = event_getChild(event, i);
		while(event_getChildNumber(event2) == 1 && unaryEventFilterFn != NULL && !unaryEventFilterFn(event2)) {
			//skip the event
			event2 = event_getChild(event2, 0);
		}
		event_construct(event_getMetaEvent(event2), event_getBranchLength(event2),
						eventTree_getEvent(eventTree, event_getName(event)), eventTree);
		eventTree_copyConstructP(eventTree, event2, unaryEventFilterFn);
	}
}

EventTree *eventTree_copyConstruct(EventTree *eventTree, Net *newNet,
		int32_t (unaryEventFilterFn)(Event *event)) {
	EventTree *eventTree2;
	eventTree2 = eventTree_construct(event_getMetaEvent(eventTree_getRootEvent(eventTree)), newNet);
	eventTree_copyConstructP(eventTree2, eventTree_getRootEvent(eventTree), unaryEventFilterFn);
	return eventTree2;
}

Event *eventTree_getRootEvent(EventTree *eventTree) {
	return eventTree->rootEvent;
}

Event *eventTree_getEvent(EventTree *eventTree, Name eventName) {
	Event *event = event_getStaticNameWrapper(eventName);
	return sortedSet_find(eventTree->events, event);
}

Event *eventTree_getCommonAncestor(Event *event, Event *event2) {
	Event *ancestorEvent;
	struct List *list;

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
	assert(FALSE);
	return NULL;
}

Net *eventTree_getNet(EventTree *eventTree) {
	return eventTree->net;
}

int32_t eventTree_getEventNumber(EventTree *eventTree) {
	return event_getSubTreeEventNumber(eventTree_getRootEvent(eventTree)) + 1;
}

Event *eventTree_getFirst(EventTree *eventTree) {
	return sortedSet_getFirst(eventTree->events);
}

EventTree_Iterator *eventTree_getIterator(EventTree *eventTree) {
	return iterator_construct(eventTree->events);
}

Event *eventTree_getNext(EventTree_Iterator *iterator) {
	return iterator_getNext(iterator);
}

Event *eventTree_getPrevious(EventTree_Iterator *iterator) {
	return iterator_getPrevious(iterator);
}

EventTree_Iterator *eventTree_copyIterator(EventTree_Iterator *iterator) {
	return iterator_copy(iterator);
}

void eventTree_destructIterator(EventTree_Iterator *iterator) {
	iterator_destruct(iterator);
}

static char *eventTree_makeNewickStringP(Event *event) {
	int32_t i;
	char *cA;
	char *cA2;
	char *cA3;
	if(event_getChildNumber(event) > 0) {
		for(i=0;i<event_getChildNumber(event); i++) {
			cA2 = eventTree_makeNewickStringP(event_getChild(event, i));
			if(i > 0) {
				cA3 = malloc(sizeof(char)*(strlen(cA)+strlen(cA2)+2));
				sprintf(cA3, "%s,%s", cA, cA2);
				free(cA);
				cA = cA3;
			}
			else {
				cA = malloc(sizeof(char)*(strlen(cA2)+2));
				sprintf(cA, "(%s", cA2);
			}
		}
		cA3 = malloc(sizeof(char)*(strlen(cA) + strlen(event_getHeader(event)) + 30));
		sprintf(cA3, "%s)%s:%g", cA, netMisc_nameToStringStatic(event_getName(event)), event_getBranchLength(event));
		free(cA);
		cA = cA3;
	}
	else {
		cA = malloc(sizeof(char)*(strlen(event_getHeader(event)) + 30));
		sprintf(cA, "%s:%g", netMisc_nameToStringStatic(event_getName(event)), event_getBranchLength(event));
	}
	return cA;
}

char *eventTree_makeNewickString(EventTree *eventTree) {
	Event *rootEvent = eventTree_getRootEvent(eventTree);
	char *cA = eventTree_makeNewickStringP(rootEvent);
	char *cA2 = malloc(sizeof(char)*(strlen(cA) + 2));
	sprintf(cA2, "%s;", cA);
	free(cA);
	return cA2;
}

/*
 * Private functions.
 */

void eventTree_destruct(EventTree *eventTree) {
	Event *event;
	while((event = eventTree_getFirst(eventTree)) != NULL) {
		event_destruct(event);
	}
	sortedSet_destruct(eventTree->events, NULL);
	free(eventTree);
}

void eventTree_addEvent(EventTree *eventTree, Event *event) {
	sortedSet_insert(eventTree->events, event);
}

void eventTree_removeEvent(EventTree *eventTree, Event *event) {
	sortedSet_delete(eventTree->events, event);
}

/*
 * Serialisation functions
 */

void eventTree_writeBinaryRepresentationP(Event *event, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	int32_t i;
	event_writeBinaryRepresentation(event, writeFn);
	for(i=0; i<event_getChildNumber(event); i++) {
		eventTree_writeBinaryRepresentationP(event_getChild(event, i), writeFn);
	}
}

void eventTree_writeBinaryRepresentation(EventTree *eventTree, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	int32_t i;
	Event *event;
	event = eventTree_getRootEvent(eventTree);
	binaryRepresentation_writeElementType(CODE_EVENT_TREE, writeFn);
	binaryRepresentation_writeName(event_getName(event), writeFn);
	binaryRepresentation_writeInteger(eventTree_getEventNumber(eventTree)-1, writeFn);
	for(i=0; i<event_getChildNumber(event); i++) {
		eventTree_writeBinaryRepresentationP(event_getChild(event, i), writeFn);
	}
}

EventTree *eventTree_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	EventTree *eventTree;
	MetaEvent *metaEvent;
	eventTree = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_EVENT_TREE) {
		binaryRepresentation_popNextElementType(binaryString);
		metaEvent = netDisk_getMetaEvent(net_getNetDisk(net),
				binaryRepresentation_getName(binaryString));
		assert(metaEvent != NULL);
		eventTree = eventTree_construct(metaEvent, net);
		int32_t eventNumber = binaryRepresentation_getInteger(binaryString);
		while(eventNumber-- > 0) {
			event_loadFromBinaryRepresentation(binaryString, eventTree);
		}
	}
	return eventTree;
}
