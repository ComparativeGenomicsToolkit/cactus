#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "net2.h"
#include "net2Private.h"

/*
 * Implementation of the net API
 */

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions on a sorted set and its iterator
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct avl_table *sortedSet_construct(int32_t (*compareFn)(const void *, const void *, void *)) {
	return avl_create(compareFn, NULL, NULL);
}

void sortedSet_destruct(struct avl_table *sortedSet, void (*destructElementFn)(void *, void *)) {
	avl_destroy(sortedSet, destructElementFn);
}

void sortedSet_insert(struct avl_table *sortedSet, void *object) {
	avl_insert(sortedSet, object);
}

void *sortedSet_find(struct avl_table *sortedSet, void *object) {
	return avl_find(sortedSet, object);
}

void sortedSet_delete(struct avl_table *sortedSet, void *object) {
	avl_delete(sortedSet, object);
}

int32_t sortedSet_getLength(struct avl_table *sortedSet) {
	return avl_count(sortedSet);
}

void *sortedSet_getFirst(struct avl_table *items) {
	static struct avl_traverser iterator;
	avl_t_init(&iterator, items);
	return avl_t_first(&iterator, items);
}

struct avl_traverser *iterator_construct(struct avl_table *items) {
	struct avl_traverser *iterator;
	iterator = mallocLocal(sizeof(struct avl_traverser));
	avl_t_init(iterator, items);
	return iterator;
}

void iterator_destruct(struct avl_traverser *iterator) {
	free(iterator);
}

void *iterator_getNext(struct avl_traverser *iterator) {
	return avl_t_next(iterator);
}

struct avl_traverser *iterator_copy(struct avl_traverser *iterator) {
	struct avl_traverser *copyIterator;
	copyIterator = mallocLocal(sizeof(struct avl_traverser));
	avl_t_copy(copyIterator, iterator);
	return copyIterator;
}

void *iterator_getPrevious(struct avl_traverser *iterator) {
	return avl_t_prev(iterator);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Database functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

TCBDB *database_construct(const char *name) {
	int32_t ecode;
	TCBDB *database;
	database = tcbdbnew();
	if(!tcbdbopen(database, name, BDBOWRITER | BDBOCREAT)) {
	   ecode = tcbdbecode(database);
	   fprintf(stderr, "Opening database error: %s\n", tcbdberrmsg(ecode));
	   exit(1);
	}
	return database;
}

void database_destruct(TCBDB *database) {
	int32_t ecode;
	if(!tcbdbclose(database)){
		ecode = tcbdbecode(database);
		fprintf(stderr, "Closing database error: %s\n", tcbdberrmsg(ecode));
		exit(1);
	}
	tcbdbdel(database);
}

int32_t database_getNumberOfRecords(TCBDB *database) {
	return tcbdbrnum(database);
}

char *database_getRecord(TCBDB *database, const char *key) {
	//Return value must be freed.
	return tcbdbget2(database, key);
}

int32_t database_writeRecord(TCBDB *database, const char *key, const char *value) {
	int32_t ecode = 0;
	if(!tcbdbput2(database, key, value)){
		ecode = tcbdbecode(database);
		fprintf(stderr, "Adding net to database error: %s\n", tcbdberrmsg(ecode));
	}
	return ecode;
}

int32_t database_removeRecord(TCBDB *database, const char *key) {

}

BDBCUR *databaseIterator_construct(TCBDB *database) {
	BDBCUR *iterator;
	iterator = tcbdbcurnew(database);
	tcbdbcurfirst(iterator);
	return iterator;
}

const char *databaseIterator_getNext(BDBCUR *iterator) {
	return tcbdbcurkey2(iterator);
}

void databaseIterator_destruct(BDBCUR *iterator) {
	tcbdbcurdel(iterator);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions for serialising the objects.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void binaryRepresentation_writeElementType(int32_t elementCode, void (*writeFn)(const char *, ...)) {

}

void binaryRepresentation_writeString(const char *name, void (*writeFn)(const char *, ...)) {

}

void binaryRepresentation_writeInteger(int32_t i, void (*writeFn)(const char *, ...)) {

}

void binaryRepresentation_writeFloat(float f, void (*writeFn)(const char *, ...)) {

}

int32_t binaryRepresentation_peekNextElementType(char **binaryString) {

}

int32_t binaryRepresentation_popNextElementType(char **binaryString) {

}

char *binaryRepresentation_getString(char **binaryString) {

}

const char *binaryRepresentation_getStringStatic(char **binaryString) {

}

int32_t binaryRepresentation_getInteger(char **binaryString) {

}

float binaryRepresentation_getFloat(char **binaryString) {

}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Event *event_construct(const char *name, float branchLength,
		Event *parentEvent, EventTree *eventTree) {
	Event *event;
	event = malloc(sizeof(Event));
	event->name = stringCopy(name);
	event->parent = parentEvent;
	event->branchLength = branchLength;
	if(parentEvent != NULL) {
		listAppend(parentEvent->children, event);
	}
	event->eventTree = eventTree;
	eventTree_addEvent(eventTree, event);
	return event;
}

Event *event_construct2(const char *name, float branchLength,
		Event *parentEvent, Event *childEvent, EventTree *eventTree) {
	Event *event;
	event = event_construct(name, branchLength, parentEvent, eventTree);
#ifdef BEN_DEBUG
	assert(parentEvent != NULL);
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
	return event->parent;
}

const char *event_getName(Event *event) {
	return event->name;
}

float event_getBranchLength(Event *event) {
	return event->branchLength;
}

int32_t event_getChildNumber(Event *event) {
	return event->children->length;
}

Event *event_getChild(Event *event, int32_t index) {
#ifdef BEN_DEBUG
	assert(index >= 0);
	assert(index < event_getChildNumber(event));
#endif
	return event->children->list[index];
}

EventTree *event_getEventTree(Event *event) {
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

/*
 * Private functions
 */

void event_destruct(Event *event) {
	int32_t i;
	eventTree_removeEvent(event_getEventTree(event), event);
	for(i=0; i<event_getChildNumber(event); i++) {
		event_destruct(event_getChild(event, i));
	}
	free(event->name);
	free(event);
}

void event_writeBinaryRepresentation(Event *event,
		void (*writeFn)(const char *string, ...)) {
	binaryRepresentation_writeElementType(CODE_EVENT, writeFn);
	binaryRepresentation_writeString(event_getName(event_getParent(event)), writeFn);
	binaryRepresentation_writeString(event_getName(event), writeFn);
	binaryRepresentation_writeFloat(event_getBranchLength(event), writeFn);
}

Event *event_loadFromBinaryRepresentation(char **binaryString,
		EventTree *eventTree) {
	Event *event;
	char *parentName;
	const char *name;
	float branchLength;
	
	event = NULL;
	if(binaryRepresentation_peekNextElementType(binaryString) == CODE_EVENT) {
		binaryRepresentation_popNextElementType(binaryString);
		parentName = binaryRepresentation_getString(binaryString);
		branchLength = binaryRepresentation_getFloat(binaryString);
		name = binaryRepresentation_getStringStatic(binaryString);
		event = event_construct(name, branchLength, NULL, eventTree);
		free(parentName);
	}
	return event;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event tree functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t eventTree_constructP(const void *o1, const void *o2, void *a) {
	return strcmp(event_getName((Event *)o1), event_getName((Event *)o2));
}

EventTree *eventTree_construct(const char *rootEventName, Net *net) {
	EventTree *eventTree;
	eventTree = malloc(sizeof(EventTree));
	eventTree->rootEvent = event_construct(rootEventName, INT32_MAX, NULL, eventTree);
	eventTree->events = sortedSet_construct(eventTree_constructP);
	eventTree->net = net;
	return eventTree;
}

void eventTree_copyConstructP(EventTree *eventTree, Event *event,
		int32_t (unaryEventFilterFn)(Event *event)) {
	int32_t i;
	Event *event2;
	for(i=0; i<event_getChildNumber(event); i++) {
		event2 = event_getChild(event, i);
		while(event_getChildNumber(event2) == 1 && !unaryEventFilterFn(event2)) {
			//skip the event
			event2 = event_getChild(event2, 0);
		}
		event_construct(event_getName(event2), event_getBranchLength(event2),
						eventTree_getEvent(eventTree, event_getName(event)), eventTree);
		eventTree_copyConstructP(eventTree, event2, unaryEventFilterFn);
	}
}

EventTree *eventTree_copyConstruct(EventTree *eventTree, Net *newNet,
		int32_t (unaryEventFilterFn)(Event *event)) {
	EventTree *eventTree2;
	eventTree2 = eventTree_construct(event_getName(eventTree_getRootEvent(eventTree)), newNet);
	eventTree_copyConstructP(eventTree2, eventTree_getRootEvent(eventTree), unaryEventFilterFn);
	return eventTree2;
}

Event *eventTree_getRootEvent(EventTree *eventTree) {
	return eventTree->rootEvent;
}

Event *eventTree_getEvent(EventTree *eventTree, const char *eventName) {
	static Event event;
	event.name = (char *)eventName;
	return sortedSet_find(eventTree->events, &event);
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
			return event;
		}
	}
	assert(FALSE);
	return NULL;
}

Net *eventTree_getNet(EventTree *eventTree) {
	return eventTree->net;
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

void end_removeEvent(EventTree *eventTree, EndInstance *event) {
	sortedSet_delete(eventTree->events, event);
}


void eventTree_writeBinaryRepresentationP(Event *event, void (*writeFn)(const char *string, ...)) {
	int32_t i;
	event_writeBinaryRepresentation(event, writeFn);
	for(i=0; i<event_getChildNumber(event); i++) {
		eventTree_writeBinaryRepresentationP(event, writeFn);
	}
}

void eventTree_writeBinaryRepresentation(EventTree *eventTree, void (*writeFn)(const char *string, ...)) {
	int32_t i;
	Event *event;
	event = eventTree_getRootEvent(eventTree);
	binaryRepresentation_writeElementType(CODE_EVENT_TREE, writeFn);
	binaryRepresentation_writeString(event_getName(event), writeFn);
	for(i=0; i<event_getChildNumber(event); i++) {
		eventTree_writeBinaryRepresentationP(event_getChild(event, i), writeFn);
	}
}

EventTree *eventTree_loadFromBinaryRepresentation(char **binaryString, Net *net) {
	EventTree *eventTree;
	eventTree = NULL;
	if(binaryRepresentation_peekNextElementType(binaryString) == CODE_EVENT_TREE) {
		binaryRepresentation_popNextElementType(binaryString);
		eventTree = eventTree_construct(binaryRepresentation_getStringStatic(binaryString), net);
		while(binaryRepresentation_peekNextElementType(binaryString) == CODE_EVENT) {
			event_loadFromBinaryRepresentation(binaryString, eventTree);
		}
	}
	return eventTree;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Sequence *sequence_construct(const char *name, int32_t length,
							 const char *file, const char *eventName, NetDisk *netDisk) {
	Sequence *sequence;
#ifdef BEN_DEBUG
	assert(length >= 0);
#endif
	sequence = malloc(sizeof(Sequence));
	sequence->file = stringCopy(file);
	sequence->name = stringCopy(name);
	sequence->length = length;
	sequence->eventName = stringCopy(eventName);
	sequence->netDisk = netDisk;

	netDisk_addSequence(netDisk, sequence);

	return sequence;
}

void sequence_destruct(Sequence *sequence) {
	netDisk_unloadSequence(sequence->netDisk, sequence);
	free(sequence->file);
	free(sequence->name);
	free(sequence);
}

int32_t sequence_getLength(Sequence *sequence) {
	return sequence->length;
}

const char *sequence_getName(Sequence *sequence) {
	return sequence->name;
}

const char *sequence_getEventName(Sequence *sequence) {
	return sequence->eventName;
}

const char *sequence_getFile(Sequence *sequence) {
	return sequence->file;
}

char *sequence_makeXMLRepresentation(Sequence *sequence) {

}

Sequence *sequence_loadFromXMLRepresentation(char *xmlString, NetDisk *netDisk) {

}

/*
 * Private functions
 */

void sequence_writeBinaryRepresentation(Sequence *sequence, void (*writeFn)(const char *string, ...)) {
	binaryRepresentation_writeElementType(CODE_SEQUENCE, writeFn);
	binaryRepresentation_writeString(sequence_getName(sequence), writeFn);
	binaryRepresentation_writeInteger(sequence_getLength(sequence), writeFn);
	binaryRepresentation_writeString(sequence_getEventName(sequence), writeFn);
	binaryRepresentation_writeString(sequence_getFile(sequence), writeFn);
}

char *sequence_makeBinaryRepresentation(Sequence *sequence) {

}

Sequence *sequence_loadFromBinaryRepresentation(char **binaryString, NetDisk *netDisk) {
	Sequence *sequence;
	char *name;
	int32_t length;
	char *eventName;
	const char *fileString;
	sequence = NULL;
	if(binaryRepresentation_peekNextElementType(binaryString) == CODE_SEQUENCE) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getString(binaryString);
		length = binaryRepresentation_getInteger(binaryString);
		eventName = binaryRepresentation_getString(binaryString);
		fileString = binaryRepresentation_getStringStatic(binaryString);
		sequence = sequence_construct(name, length, fileString, eventName, netDisk);
		free(name);
		free(eventName);
	}
	return sequence;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic end instance functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

EndInstance *endInstance_construct(const char *instance, End *end) {
	EndInstance *endInstance;
#ifdef BEN_DEBUG
	assert(insinstance != NULL);
#endif

	endInstance = malloc(sizeof(EndInstance));
	endInstance->instance = stringCopy(instance);
	endInstance->end = end;

	endInstance->coordinate = INT32_MAX;
	endInstance->sequence = NULL;
	endInstance->adjacency = NULL;
	endInstance->adjacency2 = NULL;
	endInstance->operation = NULL;
	endInstance->atomInstance = NULL;
	endInstance->parent = NULL;
	endInstance->children = constructEmptyList(0, NULL);

	end_addInstance(end, endInstance);
	return endInstance;
}

EndInstance *endInstance_constructWithCoordinates(const char *instance, End *end, int32_t coordinate, Sequence *sequence) {
	EndInstance *endInstance;
	endInstance = endInstance_construct(instance, end);
	endInstance->coordinate = coordinate;
	endInstance->sequence = sequence;
	endInstance->event = sequence_getEvent(sequence);
	return endInstance;
}

void endInstance_setEvent(EndInstance *endInstance, Event *event) {
#ifdef BEN_DEBUG
	assert(endInstance_getEvent(endInstance) == NULL);
#endif
	endInstance->event = event;
}

void endInstance_destruct(EndInstance *endInstance) {
	//Remove from end.
	end_removeInstance(endInstance_getEnd(endInstance), endInstance);

	destructList(endInstance->children);
	free(endInstance->instance);
	free(endInstance);
}

const char *endInstance_getInstanceName(EndInstance *endInstance) {
	return endInstance->instance;
}

const char *endInstance_getElementName(EndInstance *endInstance) {
	return end_getName(endInstance_getEnd(endInstance));
}

char *endInstance_getCompleteName(EndInstance *endInstance) {
	return netMisc_makeCompleteName(endInstance_getElementName(endInstance), endInstance_getInstanceName(endInstance));
}

Event *endInstance_getEvent(EndInstance *endInstance) {
	return endInstance->event;
}

End *endInstance_getEnd(EndInstance *endInstance) {
	return endInstance->end;
}

AtomInstance *endInstance_getAtomInstance(EndInstance *endInstance) {
	return endInstance->atomInstance;
}

int32_t endInstance_getCoordinate(EndInstance *endInstance) {
	return endInstance->coordinate;
}

Sequence *endInstance_getSequence(EndInstance *endInstance) {
	return endInstance->sequence;
}

void endInstance_makeAdjacent1(EndInstance *endInstance, EndInstance *endInstance2) {
	endInstance_breakAdjacency1(endInstance);
	endInstance_breakAdjacency1(endInstance2);
	endInstance->adjacency = endInstance2;
	endInstance2->adjacency = endInstance;
}

void endInstance_makeAdjacent2(EndInstance *endInstance, EndInstance *endInstance2) {
	endInstance_breakAdjacency2(endInstance);
	endInstance_breakAdjacency2(endInstance2);
	endInstance->adjacency2 = endInstance2;
	endInstance2->adjacency2 = endInstance;
}

EndInstance *endInstance_getAdjacency(EndInstance *endInstance) {
	return endInstance->adjacency;
}

EndInstance *endInstance_getAdjacency2(EndInstance *endInstance) {
	return endInstance->adjacency2;
}

Operation *endInstance_getOperation(EndInstance *endInstance) {
	return endInstance->operation;
}

EndInstance *endInstance_getParent(EndInstance *endInstance) {
	return endInstance->parent;
}

int32_t endInstance_getChildNumber(EndInstance *endInstance) {
	return endInstance->children->length;
}

EndInstance *endInstance_getChild(EndInstance *endInstance, int32_t index) {
#ifdef BEN_DEBUG
	assert(endInstance_getChildNumber(endInstance) > index);
	assert(index >= 0);
#endif
	return endInstance->children->list[index];
}

void endInstance_makeParentAndChild(EndInstance *endInstanceParent, EndInstance *endInstanceChild) {
	if(!listContains(endInstanceParent->children, endInstanceChild)) { //defensive, means second calls will have no effect.
		listAppend(endInstanceParent->children, endInstanceChild);
	}
	endInstanceChild->parent = endInstanceParent;
}

int32_t endInstance_isInternal(EndInstance *endInstance) {
	return endInstance_getChildNumber(endInstance) >= 0;
}

int32_t endInstance_isAugmented(EndInstance *endInstance) {
	return end_getAtom(endInstance_getEnd(endInstance)) != NULL && endInstance_getAtomInstance(endInstance) == NULL;
}

/*
 * Private functions.
 */

void endInstance_setAtomInstance(EndInstance *endInstance, AtomInstance *atomInstance) {
	endInstance->atomInstance = atomInstance;
}

void endInstance_setOperation(EndInstance *endInstance, Operation *operation) {
	endInstance->operation = operation;
}

void endInstance_breakAdjacency1(EndInstance *endInstance) {
	EndInstance *endInstance2;
	endInstance2 = endInstance_getAdjacency(endInstance);
	if(endInstance2 != NULL) {
		endInstance2->adjacency = NULL;
		endInstance->adjacency = NULL;
	}
}

void endInstance_breakAdjacency2(EndInstance *endInstance) {
	EndInstance *endInstance2;
	endInstance2 = endInstance_getAdjacency2(endInstance);
	if(endInstance2 != NULL) {
		endInstance2->adjacency2 = NULL;
		endInstance->adjacency2 = NULL;
	}
}

void endInstance_writeBinaryRepresentationP(EndInstance *endInstance, EndInstance *endInstance2, int32_t elementType, void (*writeFn)(const char *string, ...)) {
	char *cA;
	binaryRepresentation_writeElementType(elementType, writeFn);
	cA = endInstance_getCompleteName(endInstance2);
	binaryRepresentation_writeString(cA, writeFn);
	free(cA);
}

void endInstance_writeBinaryRepresentation(EndInstance *endInstance, void (*writeFn)(const char *string, ...)) {
	EndInstance *endInstance2;
	if(endInstance_getCoordinate(endInstance) == INT32_MAX) {
		if(endInstance_getEvent(endInstance) == NULL) {
			binaryRepresentation_writeElementType(CODE_END_INSTANCE, writeFn);
			binaryRepresentation_writeString(endInstance_getInstanceName(endInstance), writeFn);
		}
		else {
			binaryRepresentation_writeElementType(CODE_END_INSTANCE_WITH_EVENT, writeFn);
			binaryRepresentation_writeString(endInstance_getInstanceName(endInstance), writeFn);
			binaryRepresentation_writeString(event_getName(endInstance_getEvent(endInstance)), writeFn);
		}
	}
	else {
		binaryRepresentation_writeElementType(CODE_END_INSTANCE_WITH_COORDINATES, writeFn);
		binaryRepresentation_writeString(endInstance_getInstanceName(endInstance), writeFn);
		binaryRepresentation_writeInteger(endInstance_getCoordinate(endInstance), writeFn);
		binaryRepresentation_writeString(sequence_getName(endInstance_getSequence(endInstance)), writeFn);
	}
	if((endInstance2 = endInstance_getAdjacency(endInstance)) != NULL) {
		endInstance_writeBinaryRepresentationP(endInstance, endInstance2, CODE_ADJACENCY, writeFn);
	}
	if((endInstance2 = endInstance_getAdjacency2(endInstance)) != NULL) {
		endInstance_writeBinaryRepresentationP(endInstance, endInstance2, CODE_ADJACENCY, writeFn);
	}
	if((endInstance2 = endInstance_getParent(endInstance)) != NULL) {
		endInstance_writeBinaryRepresentationP(endInstance, endInstance2, CODE_PARENT, writeFn);
	}
}

int32_t endInstance_loadFromBinaryRepresentationP(EndInstance *endInstance, char **binaryString, void (*linkFn)(EndInstance *, EndInstance *)) {
	const char *cA;
	EndInstance *endInstance2;
	binaryRepresentation_popNextElementType(binaryString);
	cA = binaryRepresentation_getStringStatic(binaryString);
	endInstance2 = net_getEndInstance(end_getNet(endInstance_getEnd(endInstance)), cA);
	if(endInstance2 != NULL) { //if null we'll make the adjacency when the other end is parsed.
		linkFn(endInstance2, endInstance);
		return 0;
	}
	return 1;
}

EndInstance *endInstance_loadFromBinaryRepresentation(char **binaryString, End *end) {
	EndInstance *endInstance;
	char *name;
	int32_t coordinate;
	Sequence *sequence;

	endInstance = NULL;
	if(binaryRepresentation_peekNextElementType(binaryString) == CODE_END_INSTANCE) {
		binaryRepresentation_popNextElementType(binaryString);
		endInstance = endInstance_construct(binaryRepresentation_getStringStatic(binaryString), end);
	}
	else if(binaryRepresentation_peekNextElementType(binaryString) == CODE_END_INSTANCE_WITH_EVENT) {
		binaryRepresentation_popNextElementType(binaryString);
		endInstance = endInstance_construct(binaryRepresentation_getStringStatic(binaryString), end);
		endInstance_setEvent(endInstance, eventTree_getEvent(net_getEventTree(end_getNet(end)), binaryRepresentation_getStringStatic(binaryString)));
	}
	else if(binaryRepresentation_peekNextElementType(binaryString) == CODE_END_INSTANCE_WITH_COORDINATES) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getString(binaryString);
		coordinate = binaryRepresentation_getInteger(binaryString);
		sequence = net_getSequence(end_getNet(end), binaryRepresentation_getStringStatic(binaryString));
		endInstance = endInstance_constructWithCoordinates(name, end, coordinate, sequence);
		free(name);
	}
	if(binaryRepresentation_peekNextElementType(binaryString) == CODE_ADJACENCY) {
		endInstance_loadFromBinaryRepresentationP(endInstance, binaryString, endInstance_makeAdjacent1);
	}
	if(binaryRepresentation_peekNextElementType(binaryString) == CODE_ADJACENCY) {
		endInstance_loadFromBinaryRepresentationP(endInstance, binaryString, endInstance_makeAdjacent2);
	}
	if(binaryRepresentation_peekNextElementType(binaryString) == CODE_PARENT) {
		assert(endInstance_loadFromBinaryRepresentationP(endInstance, binaryString, endInstance_makeParentAndChild) == 0);
	}

	return endInstance;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic end functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t end_constructP(const void *o1, const void *o2, void *a) {
	return strcmp(endInstance_getInstanceName((EndInstance *)o1), endInstance_getInstanceName((EndInstance *)o2));
}

End *end_construct(const char *name, Net *net) {
	End *end;
	end = malloc(sizeof(End));
	end->endInstances = sortedSet_construct(end_constructP);
	end->attachedAtom = NULL;
	end->adjacencyComponent = NULL;
	end->net = net;
	net_addEnd(net, end);
	return end;
}

End *end_copyConstruct(End *end, Net *newNet) {
	End *end2 = end_construct(end_getName(end), newNet);
	End_InstanceIterator *iterator;
	EndInstance *endInstance;
	EndInstance *endInstance2;

	//Copy the instances.
	iterator = end_getInstanceIterator(end);
	while((endInstance = end_getNext(iterator)) != NULL) {
		if(endInstance_getCoordinate(endInstance) != INT32_MAX) {
			endInstance_constructWithCoordinates(endInstance_getInstanceName(endInstance), end2,
					endInstance_getCoordinate(endInstance), endInstance_getSequence(endInstance));
		}
		else {
			endInstance_construct(endInstance_getInstanceName(endInstance), end2);
		}
	}
	end_destructInstanceIterator(iterator);

	//Copy any parent child links.
	iterator = end_getInstanceIterator(end);
	while((endInstance = end_getNext(iterator)) != NULL) {
		if((endInstance2 = endInstance_getParent(endInstance)) != NULL) {
			endInstance_linkParentAndChild(end_getInstance(end2, endInstance_getInstanceName(endInstance2)),
										   end_getInstance(end2, endInstance_getInstanceName(endInstance)));
		}
	}
	end_destructInstanceIterator(iterator);
	return end2;
}

void end_destruct(End *end) {
	EndInstance *endInstance;
	//remove from net.
	net_removeEnd(end_getNet(end), end);

	//remove instances
	while((endInstance = end_getFirst(end)) != NULL) {
		endInstance_destruct(endInstance);
	}
	//now the actual instances.
	sortedSet_destruct(end->endInstances, NULL);

	free(end->name);
	free(end);
}

const char *end_getName(End *end) {
	return end->name;
}

Net *end_getNet(End *end) {
	return end->net;
}

AdjacencyComponent *end_getAdjacencyComponent(End *end) {
	return end->adjacencyComponent;
}

int32_t end_getInstanceNumber(End *end) {
	return sortedSet_getLength(end->endInstances);
}

EndInstance *end_getInstance(End *end, const char *name) {
	static EndInstance endInstance;
	endInstance.instance = (char *)name;
	return sortedSet_find(end->endInstances, &endInstance);
}

EndInstance *end_getFirst(End *end) {
	return sortedSet_getFirst(end->endInstances);
}

End_InstanceIterator *end_getInstanceIterator(End *end) {
	return iterator_construct(end->endInstances);
}

EndInstance *end_getNext(End_InstanceIterator *iterator) {
	return iterator_getNext(iterator);
}

EndInstance *end_getPrevious(End_InstanceIterator *iterator) {
	return iterator_getPrevious(iterator);
}

End_InstanceIterator *end_copyInstanceIterator(End_InstanceIterator *iterator) {
	return iterator_copy(iterator);
}

void end_destructInstanceIterator(End_InstanceIterator *iterator) {
	iterator_destruct(iterator);
}

int32_t end_isStub(End *end) {
	return end_getName(end)[0] == '$';
}

int32_t end_isCap(End *end) {
	return !end_isStub(end) && end_getAtom(end) != NULL;
}

int32_t end_isAtomEnd(End *end) {
	return end_getAtom(end) != NULL;
}

/*
 * Private functions
 */

void end_addInstance(End *end, EndInstance *endInstance) {
	sortedSet_insert(end->endInstances, endInstance);
}

void end_removeInstance(End *end, EndInstance *endInstance) {
	sortedSet_delete(end->endInstances, endInstance);
}

void end_setAdjacencyComponent(End *end, AdjacencyComponent *adjacencyComponent) {
	//argument may be NULL
	end->adjacencyComponent = adjacencyComponent;
}

void end_writeBinaryRepresentationP(EndInstance *endInstance, void (*writeFn)(const char *string, ...)) {
	int32_t i;
	endInstance_writeBinaryRepresentation(endInstance, writeFn);
	for(i=0; i<endInstance_getChildNumber(endInstance); i++) {
		end_writeBinaryRepresentationP(endInstance_getChild(endInstance, i), writeFn);
	}
}

void end_writeBinaryRepresentation(End *end, void (*writeFn)(const char *string, ...)) {
	End_InstanceIterator *iterator;
	EndInstance *endInstance;

	binaryRepresentation_writeElementType(CODE_END, writeFn);
	binaryRepresentation_writeString(end_getName(end), writeFn);
	endInstance = end_getRootEndInstance(end);

	if(endInstance == NULL) {
		iterator = end_getInstanceIterator(end);
		while((endInstance = end_getNext(iterator)) != NULL) {
			endInstance_writeBinaryRepresentation(endInstance, writeFn);
		}
		end_destructInstanceIterator(iterator);
	}
	else {
		end_writeBinaryRepresentationP(endInstance, writeFn);
	}
}

End *end_loadFromBinaryRepresentation(char **binaryString, Net *net) {
	End *end;

	end = NULL;
	if(binaryRepresentation_peekNextElementType(binaryString) == CODE_END) {
		binaryRepresentation_popNextElementType(binaryString);
		end = end_construct(binaryRepresentation_getStringStatic(binaryString), net);
		while(binaryRepresentation_peekNextElementType(binaryString) == CODE_END_INSTANCE) {
			endInstance_loadFromBinaryRepresentation(binaryString, end);
		}
	}
	return end;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic atom instance functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

AtomInstance *atomInstance_construct(const char *instance, Atom *atom,
		EndInstance *leftEndInstance, EndInstance *rightEndInstance) {
	AtomInstance *atomInstance;
	atomInstance = malloc(sizeof(AtomInstance));
	atomInstance->rInstance = malloc(sizeof(AtomInstance));
	atomInstance->rInstance->rInstance = atomInstance;
	atomInstance->atom = atom;
	atomInstance->rInstance->atom = atom_getReverse(atom);
	atomInstance->leftEndInstance = leftEndInstance;
	atomInstance->rInstance->leftEndInstance = rightEndInstance;
	endInstance_setAtomInstance(atomInstance->leftEndInstance, atomInstance);
	endInstance_setAtomInstance(atomInstance->rInstance->leftEndInstance, atomInstance->rInstance);
	atom_addInstance(atom, atomInstance);
	return atomInstance;
}

AtomInstance *atomInstance_construct2(const char *instance, Atom *atom) {
	return atomInstance_construct(instance, atom,
			endInstance_construct(instance, atom_getLeft(atom)),
			endInstance_construct(instance, atom_getRight(atom)));
}

AtomInstance *atomInstance_construct3(const char *instance, Atom *atom,
		int32_t startCoordinate, Sequence *sequence) {
#ifdef BEN_DEBUG
	assert(startCoordinate >= 0);
	assert(startCoordinate + atom_getLength(atom) <= sequence_getLength(sequence));
#endif

	return atomInstance_construct(instance, atom,
			endInstance_constructWithCoordinates(instance, atom_getLeft(atom),
					startCoordinate, sequence),
			endInstance_constructWithCoordinates(instance, atom_getRight(atom),
						-(startCoordinate + atom_getLength(atom)), sequence));
}

void atomInstance_destruct(AtomInstance *atomInstance) {
	atom_removeInstance(atomInstance_getAtom(atomInstance), atomInstance);
	free(atomInstance->rInstance);
	free(atomInstance);
}

Atom *atomInstance_getAtom(AtomInstance *atomInstance) {
	return atomInstance->atom;
}

const char *atomInstance_getInstanceName(AtomInstance *atomInstance) {
	return endInstance_getInstanceName(atomInstance_getLeft(atomInstance));
}

const char *atomInstance_getElementName(AtomInstance *atomInstance) {
	return atom_getName(atomInstance_getAtom(atomInstance));
}

char *atomInstance_getCompleteName(AtomInstance *atomInstance) {
	return netMisc_makeCompleteName(atomInstance_getInstanceName(atomInstance), atomInstance_getElementName(atomInstance));
}

Event *atomInstance_getEvent(AtomInstance *atomInstance) {
	return endInstance_getEvent(atomInstance_getLeft(atomInstance));
}

AtomInstance *atomInstance_getReverse(AtomInstance *atomInstance) {
	return atomInstance->rInstance;
}

int32_t atomInstance_getStart(AtomInstance *atomInstance) {
	return endInstance_getCoordinate(atomInstance_getLeft(atomInstance));
}

int32_t atomInstance_getLength(AtomInstance *atomInstance) {
	return atom_getLength(atomInstance_getAtom(atomInstance));
}

Sequence *atomInstance_getSequence(AtomInstance *atomInstance) {
	return endInstance_getSequence(atomInstance_getLeft(atomInstance));
}

EndInstance *atomInstance_getLeft(AtomInstance *atomInstance) {
	return atomInstance->leftEndInstance;
}

EndInstance *atomInstance_getRight(AtomInstance *atomInstance) {
	return atomInstance_getLeft(atomInstance_getReverse(atomInstance));
}

AtomInstance *atomInstance_getParent(AtomInstance *atomInstance) {
	EndInstance *endInstance;
	AtomInstance *atomInstance2;
	endInstance = atomInstance_getLeft(atomInstance);
	while((endInstance = endInstance_getParent(endInstance)) != NULL) {
		if((atomInstance2 = endInstance_getAtomInstance(endInstance)) != NULL) {
			return atomInstance2;
		}
	}
	return NULL;
}

int32_t atomInstance_getChildNumber(AtomInstance *atomInstance) {
	return 0;
}

AtomInstance *atomInstance_getChild(AtomInstance *atomInstance, int32_t index) {
	return NULL;
}

/*
 * Private functions
 */

void atomInstance_writeBinaryRepresentation(AtomInstance *atomInstance, void (*writeFn)(const char *string, ...)) {
	binaryRepresentation_writeElementType(CODE_ATOM_INSTANCE, writeFn);
	binaryRepresentation_writeString(atomInstance_getInstanceName(atomInstance), writeFn);
}

AtomInstance *atomInstance_loadFromBinaryRepresentation(char **binaryString, Atom *atom) {
	const char *name;
	AtomInstance *atomInstance;

	atomInstance = NULL;
	if(binaryRepresentation_peekNextElementType(binaryString) == CODE_ATOM_INSTANCE) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getStringStatic(binaryString);
		atomInstance = atomInstance_construct(name, atom, end_getInstance(atom_getLeft(atom), name),
		end_getInstance(atom_getRight(atom), name));
	}
	return atomInstance;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic atom functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t atomConstruct_constructP(const void *o1, const void *o2, void *a) {
	return strcmp(atomInstance_getInstanceName((AtomInstance *)o1), atomInstance_getInstanceName((AtomInstance *)o2));
}

Atom *atom_construct(const char *name, int32_t length, Net *net) {
	Atom *atom;
	atom = malloc(sizeof(Atom));
	atom->rAtom = malloc(sizeof(Atom));
	atom->rAtom->rAtom = atom;
	atom->atomContents = malloc(sizeof(struct AtomContents));
	atom->rAtom->atomContents = atom->atomContents;

	atom->atomContents->name = stringCopy(name);
	atom->atomContents->atomInstances = sortedSet_construct(atomConstruct_constructP);
	atom->atomContents->length = length;
	atom->atomContents->net = net;

	return atom;
}

void atom_destruct(Atom *atom) {
	AtomInstance *atomInstance;
	//remove from net.
	net_removeAtom(atom_getNet(atom), atom);

	//remove instances
	while((atomInstance = atom_getFirst(atom)) != NULL) {
		atomInstance_destruct(atomInstance);
	}
	//now the actual instances.
	sortedSet_destruct(atom->atomContents->atomInstances, NULL);

	free(atom->rAtom);
	free(atom->atomContents->name);
	free(atom->atomContents);
	free(atom);
}

const char *atom_getName(Atom *atom) {
	return atom->atomContents->name;
}

int32_t atom_getLength(Atom *atom) {
	return atom->atomContents->length;
}

Net *atom_getNet(Atom *atom) {
	return atom->atomContents->net;
}

End *atom_getLeft(Atom *atom) {
	return atom->leftEnd;
}

End *atom_getRight(Atom *atom) {
	return atom->rAtom->leftEnd;
}

Atom *atom_getReverse(Atom *atom) {
	return atom->rAtom;
}

int32_t atom_getInstanceNumber(Atom *atom) {
	return sortedSet_getLength(atom->atomContents->atomInstances);
}

AtomInstance *atom_getInstance(Atom *atom, const char *name) {
	static AtomInstance atomInstance;
	static EndInstance endInstance;
	endInstance.instance = (char *)name;
	atomInstance.leftEndInstance = &endInstance;
	return sortedSet_find(atom->atomContents->atomInstances, &atomInstance);
}

AtomInstance *atom_getFirst(Atom *atom) {
	return sortedSet_getFirst(atom->atomContents->atomInstances);
}

Atom_InstanceIterator *atom_getInstanceIterator(Atom *atom) {
	return iterator_construct(atom->atomContents->atomInstances);
}

AtomInstance *atom_getNext(Atom_InstanceIterator *iterator) {
	return iterator_getNext(iterator);
}

AtomInstance *atom_getPrevious(Atom_InstanceIterator *iterator) {
	return iterator_getPrevious(iterator);
}

Atom_InstanceIterator *atom_copyInstanceIterator(Atom_InstanceIterator *iterator) {
	return iterator_copy(iterator);
}

void atom_destructInstanceIterator(Atom_InstanceIterator *atomInstanceIterator) {
	iterator_destruct(atomInstanceIterator);
}

/*
 * Private functions.
 */

void atom_addInstance(Atom *atom, AtomInstance *atomInstance) {
	sortedSet_insert(atom->atomContents->atomInstances, atomInstance);
}

void atom_removeInstance(Atom *atom, AtomInstance *atomInstance) {
	sortedSet_delete(atom->atomContents->atomInstances, atomInstance);
}

void atom_writeBinaryRepresentation(Atom *atom, void (*writeFn)(const char *string, ...)) {
	Atom_InstanceIterator *iterator;
	AtomInstance *atomInstance;

	binaryRepresentation_writeElementType(CODE_ATOM, writeFn);
	binaryRepresentation_writeString(atom_getName(atom), writeFn);
	binaryRepresentation_writeInteger(atom_getLength(atom), writeFn);
	iterator = atom_getInstanceIterator(atom);
	while((atomInstance = atom_getNext(iterator)) != NULL) {
		atomInstance_writeBinaryRepresentation(atomInstance, writeFn);
	}
	atom_destructInstanceIterator(iterator);
}

Atom *atom_loadFromBinaryRepresentation(char **binaryString, Net *net) {
	Atom *atom;
	const char *name;
	int32_t length;

	atom = NULL;
	if(binaryRepresentation_peekNextElementType(binaryString) == CODE_ATOM) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getStringStatic(binaryString);
		length = binaryRepresentation_getInteger(binaryString);
		atom = atom_construct(name, length, net);
		while(binaryRepresentation_peekNextElementType(binaryString) == CODE_ATOM_INSTANCE) {
			atomInstance_loadFromBinaryRepresentation(binaryString, atom);
		}
	}
	return atom;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic adjacency component functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t adjacencyComponent_constructP(const void *o1, const void *o2, void *a) {
	return strcmp(end_getName((End *)o1), end_getName((End *)o2));
}

AdjacencyComponent *adjacencyComponent_construct(Net *net, Net *nestedNet) {
	AdjacencyComponent *adjacencyComponent;
	adjacencyComponent = malloc(sizeof(AdjacencyComponent));

	adjacencyComponent->net = net;
	adjacencyComponent->chain = NULL;
	adjacencyComponent->nestedNetName = stringCopy(net_getName(nestedNet));
	adjacencyComponent->ends = sortedSet_construct(adjacencyComponent_constructP);
	adjacencyComponent_updateContainedEnds(adjacencyComponent);

	net_addAdjacencyComponent(net, adjacencyComponent);

	return adjacencyComponent;
}

void adjacencyComponent_destructP(End *end, void *o) {
	end_setAdjacencyComponent(end, NULL);
}

void adjacencyComponent_updateContainedEnds(AdjacencyComponent *adjacencyComponent) {
	Net *net;
	Net_EndIterator *iterator;
	End *end;
	End *end2;
	//wipe the slate clean.
	sortedSet_destruct(adjacencyComponent->ends, (void (*)(void *, void *))adjacencyComponent_destructP);
	adjacencyComponent->ends = sortedSet_construct(adjacencyComponent_constructP);
	//now calculate the ends
	net = adjacencyComponent_getNet(adjacencyComponent);
	iterator = net_getEndIterator(adjacencyComponent_getNestedNet(adjacencyComponent));
	while((end = net_getNextEnd(iterator)) != NULL) {
		if((end2 = net_getEnd(net, end_getName(end))) != NULL) {
			sortedSet_insert(adjacencyComponent->ends, end2);
			end_setAdjacencyComponent(end2, adjacencyComponent);
		}
	}
	net_destructEndIterator(iterator);
}

void adjacencyComponent_destruct(AdjacencyComponent *adjacencyComponent) {
	//Detach from the parent net.
	net_removeAdjacencyComponent(adjacencyComponent_getNet(adjacencyComponent), adjacencyComponent);
	sortedSet_destruct(adjacencyComponent->ends, NULL);
	//Free the memory
	free(adjacencyComponent->nestedNetName);
	free(adjacencyComponent);
}

Net *adjacencyComponent_getNet(AdjacencyComponent *adjacencyComponent) {
	return adjacencyComponent->net;
}

const char *adjacencyComponent_getNestedNetName(AdjacencyComponent *adjacencyComponent) {
	return adjacencyComponent->nestedNetName;
}

Net *adjacencyComponent_getNestedNet(AdjacencyComponent *adjacencyComponent) {
	return netDisk_getNet(net_getNetDisk(adjacencyComponent_getNestedNet(adjacencyComponent)), adjacencyComponent->nestedNetName);
}

Chain *adjacencyComponent_getChain(AdjacencyComponent *adjacencyComponent) {
	return adjacencyComponent->chain;
}

End *adjacencyComponent_getEnd(AdjacencyComponent *adjacencyComponent, const char *name) {
	static End end;
	end.name = (char *)name;
	return sortedSet_find(adjacencyComponent->ends, &end);
}

int32_t adjacencyComponent_getEndNumber(AdjacencyComponent *adjacencyComponent) {
	return sortedSet_getLength(adjacencyComponent->ends);
}

AdjacencyComponent_EndIterator *adjacencyComponent_getEndIterator(AdjacencyComponent *adjacencyComponent) {
	return iterator_construct(adjacencyComponent->ends);
}

End *adjacencyComponent_getNextEnd(AdjacencyComponent_EndIterator *endIterator) {
	return iterator_getNext(endIterator);
}

End *adjacencyComponent_getPreviousEnd(AdjacencyComponent_EndIterator *endIterator) {
	return iterator_getPrevious(endIterator);
}

AdjacencyComponent_EndIterator *adjacencyComponent_copyEndIterator(AdjacencyComponent_EndIterator *endIterator) {
	return iterator_copy(endIterator);
}

void adjacencyComponent_destructEndIterator(AdjacencyComponent_EndIterator *endIterator) {
	iterator_destruct(endIterator);
}

/*
 * Private functions.
 */

void adjacencyComponent_setChain(AdjacencyComponent *adjacencyComponent, Chain *chain) {
	//argument may be NULL
	adjacencyComponent->chain = chain;
}

void adjacencyComponent_writeBinaryRepresentation(AdjacencyComponent *adjacencyComponent, void (*writeFn)(const char *string, ...)) {
	Chain *chain;
	End *end;
	AdjacencyComponent_EndIterator *iterator;

	binaryRepresentation_writeElementType(CODE_ADJACENCY_COMPONENT, writeFn);
	binaryRepresentation_writeString(adjacencyComponent_getNestedNetName(adjacencyComponent), writeFn);
	chain = adjacencyComponent_getChain(adjacencyComponent);
	if(chain == NULL) {
		binaryRepresentation_writeElementType(CODE_ADJACENCY_COMPONENT_CHAIN, writeFn);
		binaryRepresentation_writeInteger(chain_getIndex(chain), writeFn);
	}
	iterator = adjacencyComponent_getEndIterator(adjacencyComponent);
	while((end = adjacencyComponent_getNextEnd(iterator)) != NULL) {
		binaryRepresentation_writeElementType(CODE_ADJACENCY_COMPONENT_END, writeFn);
		binaryRepresentation_writeString(end_getName(end), writeFn);
	}
	adjacencyComponent_destructEndIterator(iterator);
}

AdjacencyComponent *adjacencyComponent_loadFromBinaryRepresentation(char **binaryString, Net *net) {

}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic link functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Link *link_construct(End *leftEnd, End *rightEnd, AdjacencyComponent *adjacencyComponent, Chain *parentChain) {
	Link *link;
	link = malloc(sizeof(Link));

	link->leftEnd = leftEnd;
	link->rightEnd = rightEnd;
	link->chain = parentChain;
	link->adjacencyComponent = adjacencyComponent;

	chain_addLink(parentChain, link); //will set the link indices.
	return link;
}

void link_destruct(Link *link) {
	Link *link2;
	link_getChain(link)->linkNumber = link->linkIndex;
	if(link->pLink == NULL) {
		link_getChain(link)->link = NULL;
	}
	else {
		link->pLink->nLink = NULL;
	}
	while(link != NULL) {
		link2 = link->nLink;
		link = link->nLink;
		free(link2);
	}
}

Link *link_getNextLink(Link *link) {
	return link->nLink;
}

Link *link_getPreviousLink(Link *link) {
	return link->pLink;
}

AdjacencyComponent *link_getAdjacencyComponent(Link *link) {
	return link->adjacencyComponent;
}

End *link_getLeft(Link *link) {
	return link->leftEnd;
}

End *link_getRight(Link *link) {
	return link->rightEnd;
}

Chain *link_getChain(Link *link) {
	return link->chain;
}

int32_t link_getIndex(Link *link) {
	return link->linkIndex;
}

/*
 * Private functions
 */

void link_writeBinaryRepresentation(Link *link, void (*writeFn)(const char *string, ...)) {

}

Link *link_loadFromBinaryRepresentation(char **binaryString, Chain *chain) {

}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic chain functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Chain *chain_construct(Net *net) {
	Chain *chain;
	chain = malloc(sizeof(Chain));
	chain->net = net;
	chain->link = NULL;
	chain->linkNumber = 0;
	net_addChain(net, chain);
	return chain;
}

void chain_destruct(Chain *chain) {
	net_removeChain(chain_getNet(chain), chain);
	if(chain->link != NULL) {
		link_destruct(chain->link);
	}
	free(chain);
}

Link *chain_getLink(Chain *chain, int32_t linkIndex) {
	int32_t i;
	Link *link;

#ifdef BEN_DEBUG
	assert(linkIndex >= 0);
	assert(linkIndex < chain->linkNumber);
#endif

	i=0;
	link = chain->link;
	while(i++ < linkIndex) {
		link = link->nLink;
	}
	return link;
}

int32_t chain_getLength(Chain *chain) {
	return chain->linkNumber;
}

int32_t chain_getIndex(Chain *chain) {
	return chain->chainIndex;
}

Net *chain_getNet(Chain *chain) {
	return chain->net;
}

/*
 * Private functions
 */

void chain_addLink(Chain *chain, Link *childLink) {
	Link *pLink;
	if(chain->linkNumber != 0) {
		pLink = chain_getLink(chain, chain->linkNumber -1);
		pLink->nLink = childLink;
		childLink->pLink = pLink;
	}
	else {
		childLink->pLink = NULL;
		chain->link = childLink;
	}
	childLink->nLink = NULL;
	childLink->linkIndex = chain->linkNumber++;
}

void chain_setIndex(Chain *chain, int32_t index) {
	chain->chainIndex = index;
}

void chain_writeBinaryRepresentation(Chain *chain, void (*writeFn)(const char *string, ...)) {

}

Chain *chain_loadFromBinaryRepresentation(char **binaryString, Net *net) {

}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic operation functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Operation *operation_construct(Net *net) {
	Operation *operation;
	operation = malloc(sizeof(Operation));

	net_addOperation(net, operation);
	return operation;
}

void operation_destruct(Operation *operation) {
	net_removeOperation(operation_getNet(operation), operation);
	free(operation);
}

Net *operation_getNet(Operation *operation) {
	return operation->net;
}

/*
 * Private functions
 */
void operation_setIndex(Operation *operation, int32_t index) {
	operation->index = index;
}

void *operation_writeBinaryRepresentation(Operation *operation, void (*writeFn)(const char *string, ...)) {

}

Operation *operation_loadFromBinaryRepresentation(char **binaryString, Net *net) {

}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t net_constructSequencesP(const void *o1, const void *o2, void *a) {
	return strcmp(sequence_getName((Sequence *)o1), sequence_getName((Sequence *)o2));
}

int32_t net_constructEndsP(const void *o1, const void *o2, void *a) {
	return strcmp(end_getName((End *)o1), end_getName((End *)o2));
}

int32_t net_constructAtomsP(const void *o1, const void *o2, void *a) {
	return strcmp(atom_getName((Atom *)o1), atom_getName((Atom *)o2));
}

int32_t net_constructAdjacencyComponentsP(const void *o1, const void *o2, void *a) {
	return strcmp(adjacencyComponent_getNestedNetName((AdjacencyComponent *)o1),
			adjacencyComponent_getNestedNetName((AdjacencyComponent *)o2));
}

int32_t net_constructChainsP(const void *o1, const void *o2, void *a) {
	return chain_getIndex((Chain *)o1) - chain_getIndex((Chain *)o2);
}

int32_t net_constructOperationsP(const void *o1, const void *o2, void *a) {
	return operation_getIndex((Operation *)o1) - operation_getIndex((Operation *)o2);
}

Net *net_construct(const char *name, NetDisk *netDisk) {
	Net *net;
	net = malloc(sizeof(Net));

	net->name = stringCopy(name);

	net->sequences = sortedSet_construct(net_constructSequencesP);
	net->ends = sortedSet_construct(net_constructEndsP);
	net->atoms = sortedSet_construct(net_constructAtomsP);
	net->adjacencyComponents = sortedSet_construct(net_constructAdjacencyComponentsP);
	net->chains = sortedSet_construct(net_constructChainsP);
	net->operations = sortedSet_construct(net_constructOperationsP);

	net->parentNetName = NULL;
	net->netDisk = netDisk;
	net->operationIndex = 0;
	net->chainIndex = 0;

	netDisk_addNet(net->netDisk, net);

	return net;
}

const char *net_getName(Net *net) {
	return net->name;
}

NetDisk *net_getNetDisk(Net *net) {
	return net->netDisk;
}

void net_destruct(Net *net, int32_t recursive) {
	Net_AdjacencyComponentIterator *iterator;
	End *end;
	Atom *atom;
	AdjacencyComponent *adjacencyComponent;
	Chain *chain;
	Operation *operation;

	if(recursive) {
		iterator = net_getAdjacencyComponentIterator(net);
		while((adjacencyComponent = net_getNextAdjacencyComponent(iterator)) != NULL) {
			net_destruct(adjacencyComponent_getNestedNet(adjacencyComponent), recursive);
		}
		net_destructAdjacencyComponentIterator(iterator);
	}

	netDisk_unloadNet(net->netDisk, net);

	sortedSet_destruct(net->sequences, NULL);

	while((end = net_getFirstEnd(net)) != NULL) {
		end_destruct(end);
	}
	sortedSet_destruct(net->ends, NULL);

	while((atom = net_getFirstAtom(net)) != NULL) {
		atom_destruct(atom);
	}
	sortedSet_destruct(net->atoms, NULL);

	while((adjacencyComponent = net_getFirstAdjacencyComponent(net)) != NULL) {
		adjacencyComponent_destruct(adjacencyComponent);
	}
	sortedSet_destruct(net->adjacencyComponents, NULL);

	while((chain = net_getFirstChain(net)) != NULL) {
		chain_destruct(chain);
	}
	sortedSet_destruct(net->chains, NULL);

	while((operation = net_getFirstOperation(net)) != NULL) {
		operation_destruct(operation);
	}
	sortedSet_destruct(net->operations, NULL);

	free(net->name);
	if(net->parentNetName != NULL) {
		free(net->parentNetName);
	}
	free(net);
}

void net_addSequence(Net *net, Sequence *sequence) {
	sortedSet_insert(net->sequences, sequence);
}

Sequence *net_getFirstSequence(Net *net) {
	return sortedSet_getFirst(net->sequences);
}

Sequence *net_getSequence(Net *net, const char *name) {
	static Sequence sequence;
	sequence.name = (char *)name;
	return sortedSet_find(net->sequences, &sequence);
}

int32_t net_getSequenceNumber(Net *net) {
	return sortedSet_getLength(net->sequences);
}

Net_SequenceIterator *net_getSequenceIterator(Net *net) {
	return iterator_construct(net->sequences);
}

Sequence *net_getNextSequence(Net_SequenceIterator *sequenceIterator) {
	return iterator_getNext(sequenceIterator);
}

Sequence *net_getPreviousSequence(Net_SequenceIterator *sequenceIterator) {
	return iterator_getPrevious(sequenceIterator);
}

Net_SequenceIterator *net_copySequenceIterator(Net_SequenceIterator *sequenceIterator) {
	return iterator_copy(sequenceIterator);
}

void net_destructSequenceIterator(Net_SequenceIterator *sequenceIterator) {
	iterator_destruct(sequenceIterator);
}

End *net_getFirstEnd(Net *net) {
	return sortedSet_getFirst(net->ends);
}

End *net_getEnd(Net *net, const char *name) {
	static End end;
	end.name = (char *)name;
	return sortedSet_find(net->ends, &end);
}

EndInstance *net_getEndInstance(Net *net, const char *completeName) {
	End *end;
	end = net_getEnd(net, netMisc_getElementNameStatic(completeName));
	if(end != NULL) { //if null we'll make the adjacency when the other end is parsed.
		return end_getInstance(end, netMisc_getInstanceNameStatic(completeName));
	}
	return NULL;
}

int32_t net_getEndNumber(Net *net) {
	return sortedSet_getLength(net->ends);
}

Net_EndIterator *net_getEndIterator(Net *net) {
	return iterator_construct(net->ends);
}

End *net_getNextEnd(Net_EndIterator *endIterator) {
	return iterator_getNext(endIterator);
}

End *net_getPreviousEnd(Net_EndIterator *endIterator) {
	return iterator_getPrevious(endIterator);
}

Net_EndIterator *net_copyEndIterator(Net_EndIterator *endIterator) {
	return iterator_copy(endIterator);
}

void net_destructEndIterator(Net_EndIterator *endIterator) {
	iterator_destruct(endIterator);
}

Atom *net_getFirstAtom(Net *net) {
	return sortedSet_getFirst(net->atoms);
}

Atom *net_getAtom(Net *net, const char *name) {
	static Atom atom;
	static struct AtomContents atomContents;
	atom.atomContents = &atomContents;
	atomContents.name = (char *)name;
	return sortedSet_find(net->atoms, &atom);
}

AtomInstance *net_getAtomInstance(Net *net, const char *completeName) {
	Atom *atom;
	atom = net_getAtom(net, netMisc_getElementNameStatic(completeName));
	if(atom != NULL) { //if null we'll make the adjacency when the other end is parsed.
		return atom_getInstance(atom, netMisc_getInstanceNameStatic(completeName));
	}
	return NULL;
}

int32_t net_getAtomNumber(Net *net) {
	return sortedSet_getLength(net->atoms);
}

Net_AtomIterator *net_getAtomIterator(Net *net) {
	return iterator_construct(net->atoms);
}

Atom *net_getNextAtom(Net_AtomIterator *atomIterator) {
	return iterator_getNext(atomIterator);
}

Atom *net_getPreviousAtom(Net_AtomIterator *atomIterator) {
	return iterator_getPrevious(atomIterator);
}

Net_AtomIterator *net_copyAtomIterator(Net_AtomIterator *atomIterator) {
	return iterator_copy(atomIterator);
}

void net_destructAtomIterator(Net_AtomIterator *atomIterator) {
	iterator_destruct(atomIterator);
}

AdjacencyComponent *net_getFirstAdjacencyComponent(Net *net) {
	return sortedSet_getFirst(net->adjacencyComponents);
}

AdjacencyComponent *net_getAdjacencyComponent(Net *net, const char *netName) {
	static AdjacencyComponent adjacencyComponent;
	adjacencyComponent.nestedNetName = (char *)netName;
	return sortedSet_find(net->adjacencyComponents, &adjacencyComponent);
}

int32_t net_getAdjacencyComponentNumber(Net *net) {
	return sortedSet_getLength(net->adjacencyComponents);
}

Net_AdjacencyComponentIterator *net_getAdjacencyComponentIterator(Net *net) {
	return iterator_construct(net->adjacencyComponents);
}

AdjacencyComponent *net_getNextAdjacencyComponent(Net_AdjacencyComponentIterator *adjacencyComponentIterator) {
	return iterator_getNext(adjacencyComponentIterator);
}

AdjacencyComponent *net_getPreviousAdjacencyComponent(Net_AdjacencyComponentIterator *adjacencyComponentIterator) {
	return iterator_getPrevious(adjacencyComponentIterator);
}

Net_AdjacencyComponentIterator *net_copyAdjacencyComponentIterator(Net_AdjacencyComponentIterator *adjacencyComponentIterator) {
	return iterator_copy(adjacencyComponentIterator);
}

void net_destructAdjacencyComponentIterator(Net_AdjacencyComponentIterator *adjacencyComponentIterator) {
	iterator_destruct(adjacencyComponentIterator);
}

AdjacencyComponent *net_getParentAdjacencyComponent(Net *net) {
	Net *net2;
	net2 = netDisk_getNet(net_getNetDisk(net), net->parentNetName);
	return net_getAdjacencyComponent(net2, net_getName(net));
}

Chain *net_getFirstChain(Net *net) {
	return sortedSet_getFirst(net->chains);
}

Chain *net_getChain(Net *net, int32_t index) {
	static Chain chain;
	chain.chainIndex = index;
	return sortedSet_find(net->chains, &chain);
}

int32_t net_getChainNumber(Net *net) {
	return sortedSet_getLength(net->chains);
}

Net_ChainIterator *net_getChainIterator(Net *net) {
	return iterator_construct(net->chains);
}

Chain *net_getNextChain(Net_ChainIterator *chainIterator) {
	return iterator_getNext(chainIterator);
}

Chain *net_getPreviousChain(Net_ChainIterator *chainIterator) {
	return iterator_getPrevious(chainIterator);
}

Net_ChainIterator *net_copyChainIterator(Net_ChainIterator *chainIterator) {
	return iterator_copy(chainIterator);
}

void net_destructChainIterator(Net_ChainIterator *chainIterator) {
	iterator_destruct(chainIterator);
}

Operation *net_getFirstOperation(Net *net) {
	return sortedSet_getFirst(net->operations);
}

Operation *net_getOperation(Net *net, int32_t index) {
	static Operation operation;
	operation.index = index;
	return sortedSet_find(net->operations, &operation);
}

int32_t net_getOperationNumber(Net *net) {
	return sortedSet_getLength(net->operations);
}

Net_OperationIterator *net_getOperationIterator(Net *net) {
	return iterator_construct(net->operations);
}

Operation *net_getNextOperation(Net_OperationIterator *operationIterator) {
	return iterator_getNext(operationIterator);
}

Operation *net_getPreviousOperation(Net_OperationIterator *operationIterator) {
	return iterator_getPrevious(operationIterator);
}

Net_OperationIterator *net_copyOperationIterator(Net_OperationIterator *operationIterator) {
	return iterator_copy(operationIterator);
}

void net_destructOperationIterator(Net_OperationIterator *operationIterator) {
	iterator_destruct(operationIterator);
}

char *net_makeXMLRepresentation(Net *net) {
}


Net *net_loadFromXMLRepresentation(char *xmlString, NetDisk *netDisk) {

}

/*
 * Private functions
 */

void net_addEventTree(Net *net, EventTree *eventTree) {
	net->eventTree = eventTree;
}

void net_addAtom(Net *net, Atom *atom) {
	sortedSet_insert(net->atoms, atom);
}

void net_removeAtom(Net *net, Atom *atom) {
	sortedSet_delete(net->atoms, atom);
}

void net_addEnd(Net *net, End *end) {
	sortedSet_insert(net->ends, end);
}

void net_removeEnd(Net *net, End *end) {
	sortedSet_delete(net->ends, end);
}

void net_addChain(Net *net, Chain *chain) {
	chain_setIndex(chain, net->chainIndex++);
	sortedSet_insert(net->chains, chain);
}

void net_removeChain(Net *net, Chain *chain) {
	sortedSet_delete(net->chains, chain);
}

void net_addAdjacencyComponent(Net *net, AdjacencyComponent *adjacencyComponent) {
	sortedSet_insert(net->adjacencyComponents, adjacencyComponent);
}

void net_removeAdjacencyComponent(Net *net, AdjacencyComponent *adjacencyComponent) {
	sortedSet_delete(net->adjacencyComponents, adjacencyComponent);
}

void net_setParentAdjacencyComponent(Net *net, AdjacencyComponent *adjacencyComponent) {
	net->parentNetName = stringCopy(net_getName(adjacencyComponent_getNet(adjacencyComponent)));
}

void net_addOperation(Net *net, Operation *operation) {
	operation_setIndex(operation, net_getOperationNumber(net));
	sortedSet_insert(net->operations, operation);
}

void net_removeOperation(Net *net, Operation *operation) {
	sortedSet_delete(net->operations, operation);
}

void net_writeBinaryRepresentation(Net *net, void (*writeFn)(const char *string, ...)) {
	Net_EndIterator *endIterator;
	Net_AtomIterator *atomIterator;
	Net_AdjacencyComponentIterator *adjacencyComponentIterator;
	Net_ChainIterator *chainIterator;
	Net_OperationIterator *operationIterator;
	End *end;
	Atom *atom;
	AdjacencyComponent *adjacencyComponent;
	Chain *chain;
	Operation *operation;

	binaryRepresentation_writeString(net_getName(net), writeFn);

	endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		end_writeBinaryRepresentation(end, writeFn);
	}
	net_destructEndIterator(endIterator);

	atomIterator = net_getAtomIterator(net);
	while((atom = net_getNextAtom(atomIterator)) != NULL) {
		atom_writeBinaryRepresentation(atom, writeFn);
	}
	net_destructAtomIterator(atomIterator);

	adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		adjacencyComponent_writeBinaryRepresentation(adjacencyComponent, writeFn);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

	chainIterator = net_getChainIterator(net);
	while((chain = net_getNextChain(chainIterator)) != NULL) {
		chain_writeBinaryRepresentation(chain, writeFn);
	}
	net_destructChainIterator(chainIterator);

	operationIterator = net_getOperationIterator(net);
	while((operation = net_getNextOperation(operationIterator)) != NULL) {
		operation_writeBinaryRepresentation(operation, writeFn);
	}
	net_destructOperationIterator(operationIterator);
}

char *net_makeBinaryRepresentation(Net *net) {

}

Net *net_loadFromBinaryRepresentation(char **binaryString, NetDisk *netDisk) {
	Net *net;
	net = net_construct(binaryRepresentation_getStringStatic(binaryString), netDisk);
	while(end_loadFromBinaryRepresentation(binaryString, net) != NULL);
	while(atom_loadFromBinaryRepresentation(binaryString, net) != NULL);
	while(adjacencyComponent_loadFromBinaryRepresentation(binaryString, net) != NULL);
	while(chain_loadFromBinaryRepresentation(binaryString, net) != NULL);
	while(operation_loadFromBinaryRepresentation(binaryString, net) != NULL);
	return net;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t netDisk_constructNetsP(const void *o1, const void *o2, void *a) {
	return strcmp(net_getName((Net *)o1), net_getName((Net *)o2));
}

int32_t netDisk_constructSequencesP(const void *o1, const void *o2, void *a) {
	return strcmp(sequence_getName((Sequence *)o1), sequence_getName((Sequence *)o2));
}

NetDisk *netDisk_construct(const char *netDiskFile) {
	NetDisk *netDisk;
	netDisk = malloc(sizeof(NetDisk));

	//construct lists of in memory objects
	netDisk->sequences = sortedSet_construct(net_constructSequencesP);
	netDisk->nets = sortedSet_construct(net_constructSequencesP);

	//open the sequences database
	netDisk->sequencesDatabase = database_construct(netDisk->sequencesDatabaseName);
	netDisk->netsDatabase = database_construct(netDisk->netsDatabaseName);

	return netDisk;
}

void netDisk_destruct(NetDisk *netDisk){
	Sequence *sequence;
	Net *net;

	while((net = netDisk_getFirstNetInMemory(netDisk)) != NULL) {
		net_destruct(net, FALSE);
	}
	sortedSet_destruct(netDisk->nets, NULL);

	while((sequence = netDisk_getFirstSequenceInMemory(netDisk)) != NULL) {
		sequence_destruct(sequence);
	}
	sortedSet_destruct(netDisk->sequences, NULL);

	//close DBs
	database_destruct(netDisk->sequencesDatabase);
	database_destruct(netDisk->netsDatabase);

	free(netDisk);
}

int32_t netDisk_write(NetDisk *netDisk){
	NetDisk_NetIterator *netIterator;
	NetDisk_SequenceIterator *sequenceIterator;
	char *cA;
	Net *net;
	Sequence *sequence;
	int32_t ecode;

	netIterator = netDisk_getNetInMemoryIterator(netDisk);
	while((net = netDisk_getNextNet(netIterator)) != NULL) {
		cA = net_makeBinaryRepresentation(net);
		if((ecode = database_writeRecord(netDisk->netsDatabase, net_getName(net), cA)) != 0) {
			return ecode;
		}
		free(cA);
	}
	netDisk_destructNetIterator(netIterator);

	sequenceIterator = netDisk_getSequenceInMemoryIterator(netDisk);
	while((sequence = netDisk_getNextSequence(sequenceIterator)) != NULL) {
		cA = sequence_makeBinaryRepresentation(sequence);
		if((ecode = database_writeRecord(netDisk->sequencesDatabase, sequence_getName(sequence), cA)) != 0) {
			return ecode;
		}
		free(cA);
	}
	netDisk_destructSequenceIterator(sequenceIterator);
	return 0;
}

Sequence *netDisk_getSequence(NetDisk *netDisk, const char *sequenceName) {
	char *cA;
	char *cA2;
	Sequence *sequence;

	//try in memory list first.
	if((sequence = netDisk_getSequenceInMemory(netDisk, sequenceName)) != NULL) {
		return sequence;
	}
	//else try the database.
	cA = database_getRecord(netDisk->sequencesDatabase, sequenceName);
	if(cA == NULL) {
		return NULL;
	}
	else {
		cA2 = cA;
		sequence = sequence_loadFromBinaryRepresentation(&cA2, netDisk);
		free(cA);
		return sequence;
	}
}

int32_t netDisk_getSequenceNumberOnDisk(NetDisk *netDisk) {
	return database_getNumberOfRecords(netDisk->sequencesDatabase);
}

NetDisk_SequenceNameIterator *netDisk_getSequenceNameIterator(NetDisk *netDisk) {
	return databaseIterator_construct(netDisk->sequencesDatabase);
}

const char *netDisk_getNextSequenceName(NetDisk_SequenceNameIterator *sequenceIterator) {
	return databaseIterator_getNext(sequenceIterator);
}

void netDisk_destructSequenceNameIterator(NetDisk_SequenceNameIterator *sequenceIterator) {
	databaseIterator_destruct(sequenceIterator);
}

Sequence *netDisk_getSequenceInMemory(NetDisk *netDisk, const char *sequenceName) {
	static Sequence sequence;
	sequence.name = (char *)sequenceName;
	return sortedSet_find(netDisk->sequences, &sequence);
}

Sequence *netDisk_getFirstSequenceInMemory(NetDisk *netDisk) {
	return sortedSet_getFirst(netDisk->sequences);
}

int32_t netDisk_getSequenceNumberInMemory(NetDisk *netDisk) {
	return sortedSet_getLength(netDisk->sequences);
}

NetDisk_SequenceIterator *netDisk_getSequenceInMemoryIterator(NetDisk *netDisk) {
	return iterator_construct(netDisk->sequences);
}

Sequence *netDisk_getNextSequence(NetDisk_SequenceIterator *sequenceIterator) {
	return iterator_getNext(sequenceIterator);
}

Sequence *netDisk_getPreviousSequence(NetDisk_SequenceIterator *sequenceIterator) {
	return iterator_getPrevious(sequenceIterator);
}

NetDisk_SequenceIterator *netDisk_copySequenceIterator(NetDisk_SequenceIterator *sequenceIterator) {
	return iterator_copy(sequenceIterator);
}

void netDisk_destructSequenceIterator(NetDisk_SequenceIterator *sequenceIterator) {
	return iterator_destruct(sequenceIterator);
}

Net *netDisk_getNet(NetDisk *netDisk, const char *netName) {
	char *cA;
	char *cA2;
	Net *net;

	//try in memory list first.
	if((net = netDisk_getNetInMemory(netDisk, netName)) != NULL) {
		return net;
	}
	//else try the database.
	cA = database_getRecord(netDisk->netsDatabase, netName);
	if(cA == NULL) {
		return NULL;
	}
	else {
		cA2 = cA;
		net = net_loadFromBinaryRepresentation(&cA2, netDisk);
		free(cA);
		return net;
	}
}

int32_t netDisk_getNetNumberOnDisk(NetDisk *netDisk) {
	return database_getNumberOfRecords(netDisk->netsDatabase);
}

NetDisk_NetNameIterator *netDisk_getNetNameIterator(NetDisk *netDisk) {
	return databaseIterator_construct(netDisk->netsDatabase);
}

const char *netDisk_getNextNetName(NetDisk_NetNameIterator *netIterator) {
	return databaseIterator_getNext(netIterator);
}

void netDisk_destructNetNameIterator(NetDisk_NetNameIterator *netIterator) {
	databaseIterator_destruct(netIterator);
}

Net *netDisk_getNetInMemory(NetDisk *netDisk, const char *netName) {
	static Net net;
	net.name = (char *)netName;
	return sortedSet_find(netDisk->nets, &net);
}

Net *netDisk_getFirstNetInMemory(NetDisk *netDisk) {
	return sortedSet_getFirst(netDisk->nets);
}

int32_t netDisk_getNetNumberInMemory(NetDisk *netDisk) {
	return sortedSet_getLength(netDisk->nets);
}

NetDisk_NetIterator *netDisk_getNetInMemoryIterator(NetDisk *netDisk) {
	return iterator_construct(netDisk->nets);
}

Net *netDisk_getNextNet(NetDisk_NetIterator *netIterator) {
	return iterator_getNext(netIterator);
}

Net *netDisk_getPreviousNet(NetDisk_NetIterator *netIterator) {
	return iterator_getPrevious(netIterator);
}

NetDisk_NetIterator *netDisk_copyNetIterator(NetDisk_NetIterator *netIterator) {
	return iterator_copy(netIterator);
}

void netDisk_destructNetIterator(NetDisk_NetIterator *netIterator) {
	iterator_destruct(netIterator);
}

/*
 * Private functions.
 */

int32_t netDisk_deleteSequenceFromDisk(NetDisk *netDisk, const char *sequenceName) {
	return database_removeRecord(netDisk->sequencesDatabase, sequenceName);
}

int32_t netDisk_deleteNetFromDisk(NetDisk *netDisk, const char *netName) {
	return database_removeRecord(netDisk->netsDatabase, netName);
}

void netDisk_addSequence(NetDisk *netDisk, Sequence *sequence) {
#ifdef BEN_DEBUG
	assert(netDisk_getSequence(netDisk) == NULL);
#endif
	sortedSet_insert(netDisk->sequences, sequence);
}

void netDisk_unloadSequence(NetDisk *netDisk, Sequence *sequence) {
#ifdef BEN_DEBUG
	assert(netDisk_getSequenceInMemroy(netDisk) != NULL);
#endif
	sortedSet_delete(netDisk->sequences, sequence);
}

void netDisk_addNet(NetDisk *netDisk, Net *net) {
#ifdef BEN_DEBUG
	assert(netDisk_getNet(netDisk) == NULL);
#endif
	sortedSet_insert(netDisk->nets, net);
}

void netDisk_unloadNet(NetDisk *netDisk, Net *net) {
#ifdef BEN_DEBUG
	assert(netDisk_getNetInMemory(netDisk) == NULL);
#endif
	sortedSet_delete(netDisk->nets, net);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Useful utility functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

char *netMisc_getInstanceName(const char *completeName) {
	char *cA;
	char *cA2;
	char *cA3;
#ifdef BEN_DEBUG
	assert(completeName != NULL);
#endif
	cA = stringCopy(completeName);
	cA2 = cA;
#ifdef BEN_DEBUG
	assert(strlen(cA2) > 0);
#endif
	while(cA2[0] != '.') {
		cA2++;
#ifdef BEN_DEBUG
		assert(strlen(cA2) > 0);
#endif
	}
#ifdef BEN_DEBUG
	assert(cA2[0] == '.');
#endif
	cA3 = stringCopy(cA2+1);
	free(cA);
	return cA3;
}

static char *netMisc_getInstanceNameStatic_instanceName = NULL;
const char *netMisc_getInstanceNameStatic(const char *completeName) {
	if(netMisc_getInstanceNameStatic_instanceName != NULL) {
		free(netMisc_getInstanceNameStatic_instanceName);
	}
	netMisc_getInstanceNameStatic_instanceName = netMisc_getInstanceName(completeName);
	return netMisc_getInstanceNameStatic_instanceName;
}

char *netMisc_getElementName(const char *completeName) {
	char *cA2;
	int32_t i;

	cA2 = (char *)mallocLocal(sizeof(char)*(strlen(completeName)+1));

	i=0;
	while(completeName[i] != '.' && completeName[i] != '\0') {
		cA2[i] = completeName[i];
		i++;
	}

#ifdef BEN_DEBUG
	assert(i <= (int32_t)strlen(completeName));
#endif

	cA2[i] = '\0';
	return cA2;
}

static char *netMisc_getElementNameStatic_elementName = NULL;
const char *netMisc_getElementNameStatic(const char *completeName) {
	if(netMisc_getElementNameStatic_elementName != NULL) {
		free(netMisc_getElementNameStatic_elementName);
	}
	netMisc_getElementNameStatic_elementName = netMisc_getElementName(completeName);
	return netMisc_getElementNameStatic_elementName;
}

char *netMisc_makeCompleteName(const char *elementName, const char *instanceName) {
	char *cA;
	cA = malloc(sizeof(char) * (strlen(elementName) + strlen(instanceName) + 2));
	sprintf(cA, "%s.%s", elementName, instanceName);
	return cA;
}

static char *netMisc_makeCompleteNameStatic_completeName = NULL;
const char *netMisc_makeCompleteNameStatic(const char *elementName, const char *instanceName) {
	if(netMisc_makeCompleteNameStatic_completeName != NULL) {
		free(netMisc_makeCompleteNameStatic_completeName);
	}
	netMisc_makeCompleteNameStatic_completeName = netMisc_makeCompleteName(elementName, instanceName);
	return netMisc_makeCompleteNameStatic_completeName;
}
