#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "net.h"
#include "netPrivate.h"

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
	return tcbdbout2(database, key);
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

float event_getSubTreeBranchLength(Event *event) {
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

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event tree functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t eventTree_constructP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
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

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

MetaSequence *metaSequence_construct2(const char *name, int32_t start,
		int32_t length, int64_t fileOffset, const char *header,
		const char *eventName, NetDisk *netDisk) {
	MetaSequence *metaSequence;

	metaSequence = malloc(sizeof(MetaSequence));
	metaSequence->name = stringCopy(name);
	assert(length >= 0);
	metaSequence->start = start;
	metaSequence->length = length;
	metaSequence->fileOffset = fileOffset;
	metaSequence->eventName = stringCopy(eventName);
	metaSequence->netDisk = netDisk;
	metaSequence->referenceCount = 0;
	metaSequence->header = stringCopy(header);

	netDisk_addMetaSequence(netDisk, metaSequence);
	return metaSequence;
}

MetaSequence *metaSequence_construct(const char *name, int32_t start, int32_t length,
		const char *string, const char *header, const char *eventName, NetDisk *netDisk) {
	int64_t fileOffset;
	fileOffset = netDisk_addString(netDisk, string, length);
	return metaSequence_construct2(name, start, length, fileOffset, header, eventName, netDisk);
}

void metaSequence_destruct(MetaSequence *metaSequence) {
	assert(metaSequence->referenceCount == 0);
	netDisk_unloadMetaSequence(metaSequence->netDisk, metaSequence);
	free(metaSequence->name);
	free(metaSequence->eventName);
	free(metaSequence->header);
	free(metaSequence);
}

const char *metaSequence_getName(MetaSequence *metaSequence) {
	return metaSequence->name;
}

int32_t metaSequence_getStart(MetaSequence *metaSequence) {
	return metaSequence->start;
}

int32_t metaSequence_getLength(MetaSequence *metaSequence) {
	return metaSequence->length;
}

const char *metaSequence_getEventName(MetaSequence *metaSequence) {
	return metaSequence->eventName;
}

char *metaSequence_getString(MetaSequence *metaSequence, int32_t start, int32_t length, int32_t strand) {
	assert(start >= metaSequence_getStart(metaSequence));
	assert(start < metaSequence_getStart(metaSequence) + metaSequence_getLength(metaSequence));
	return netDisk_getString(metaSequence->netDisk, metaSequence->fileOffset, start, length, strand);
}

const char *metaSequence_getHeader(MetaSequence *metaSequence) {
	return metaSequence->header;
}

int64_t metaSequence_getFileOffset(MetaSequence *metaSequence) {
	return metaSequence->fileOffset;
}

void metaSequence_increaseReferenceCount(MetaSequence *metaSequence) {
	metaSequence->referenceCount++;
}

void metaSequence_decreaseReferenceCountAndDestructIfZero(MetaSequence *metaSequence) {
	metaSequence->referenceCount--;
	assert(metaSequence->referenceCount >= 0);
	if(metaSequence->referenceCount == 0) {
		metaSequence_destruct(metaSequence);
	}
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Sequence *sequence_construct2(MetaSequence *metaSequence, Net *net) {
	Sequence *sequence;
	sequence = malloc(sizeof(Sequence));
	sequence->metaSequence = metaSequence;
	metaSequence_increaseReferenceCount(metaSequence);
	sequence->net = net;
	net_addSequence(net, sequence);

	return sequence;
}

Sequence *sequence_construct(const char *name, int32_t start, int32_t length,
							 const char *string, const char *header, Event *event, Net *net) {
	return sequence_construct2(metaSequence_construct(name, start, length, string, header,
			event_getName(event), net_getNetDisk(net)), net);
}

Sequence *sequence_copyConstruct(Sequence *sequence, Net *newNet) {
	return sequence_construct2(sequence->metaSequence, newNet);
}

void sequence_destruct(Sequence *sequence) {
	net_removeSequence(sequence_getNet(sequence), sequence);
	metaSequence_decreaseReferenceCountAndDestructIfZero(sequence->metaSequence);
	free(sequence);
}

int32_t sequence_getStart(Sequence *sequence) {
	return metaSequence_getStart(sequence->metaSequence);
}

int32_t sequence_getLength(Sequence *sequence) {
	return metaSequence_getLength(sequence->metaSequence);
}

const char *sequence_getName(Sequence *sequence) {
	return metaSequence_getName(sequence->metaSequence);
}

Event *sequence_getEvent(Sequence *sequence) {
	return eventTree_getEvent(net_getEventTree(sequence_getNet(sequence)), metaSequence_getEventName(sequence->metaSequence));
}

Net *sequence_getNet(Sequence *sequence) {
	return sequence->net;
}

char *sequence_getString(Sequence *sequence, int32_t start, int32_t length, int32_t strand) {
	return metaSequence_getString(sequence->metaSequence, start, length, strand);
}

const char *sequence_getHeader(Sequence *sequence) {
	return metaSequence_getHeader(sequence->metaSequence);
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

	endInstance = malloc(sizeof(EndInstance));
	endInstance->endInstanceContents = malloc(sizeof(EndInstanceContents));
	endInstance->rEndInstance = malloc(sizeof(EndInstance));
	endInstance->rEndInstance->rEndInstance = endInstance;
	endInstance->rEndInstance->endInstanceContents = endInstance->endInstanceContents;

	end->orientation = 1;
	end->rEnd->orientation = -1;

	endInstance->end = end;
	endInstance->rEndInstance->end = end_getReverse(end);

	endInstance->endInstanceContents->instance = stringCopy(instance);
	endInstance->endInstanceContents->coordinate = INT32_MAX;
	endInstance->endInstanceContents->sequence = NULL;
	endInstance->endInstanceContents->adjacency = NULL;
	endInstance->endInstanceContents->adjacency2 = NULL;
	endInstance->endInstanceContents->operation = NULL;
	endInstance->endInstanceContents->atomInstance = NULL;
	endInstance->endInstanceContents->parent = NULL;
	endInstance->endInstanceContents->children = constructEmptyList(0, NULL);

	end_addInstance(end, endInstance);
	return endInstance;
}

EndInstance *endInstance_construct2(const char *instance, End *end,
		int32_t coordinate, int32_t strand, int32_t side, Sequence *sequence) {
	EndInstance *endInstance;
	endInstance = endInstance_construct(instance, end);
	endInstance->endInstanceContents->coordinate = coordinate;
	assert(strand != 0);
	assert(side != 0);
	endInstance->endInstanceContents->strand = strand;
	endInstance->endInstanceContents->side = side;
	endInstance->endInstanceContents->sequence = sequence;
	endInstance->endInstanceContents->event = sequence_getEvent(sequence);
	return endInstance;
}

void endInstance_setEvent(EndInstance *endInstance, Event *event) {
#ifdef BEN_DEBUG
	assert(endInstance_getEvent(endInstance) == NULL);
#endif
	endInstance->endInstanceContents->event = event;
}

void endInstance_destruct(EndInstance *endInstance) {
	//Remove from end.
	end_removeInstance(endInstance_getEnd(endInstance), endInstance);

	destructList(endInstance->endInstanceContents->children);
	free(endInstance->rEndInstance);
	free(endInstance->endInstanceContents->instance);
	free(endInstance->endInstanceContents);
	free(endInstance);
}

const char *endInstance_getInstanceName(EndInstance *endInstance) {
	return endInstance->endInstanceContents->instance;
}

const char *endInstance_getElementName(EndInstance *endInstance) {
	return end_getName(endInstance_getEnd(endInstance));
}

char *endInstance_getCompleteName(EndInstance *endInstance) {
	return netMisc_makeCompleteName(endInstance_getElementName(endInstance), endInstance_getInstanceName(endInstance), 1);
}

int32_t endInstance_getOrientation(EndInstance *endInstance) {
	return end_getOrientation(endInstance_getEnd(endInstance));
}

char *endInstance_getCompleteNameWithOrientation(EndInstance *endInstance) {
	return netMisc_makeCompleteName(endInstance_getElementName(endInstance),
			endInstance_getInstanceName(endInstance), endInstance_getOrientation(endInstance));
}

EndInstance *endInstance_getReverse(EndInstance *endInstance) {
	return endInstance->rEndInstance;
}

Event *endInstance_getEvent(EndInstance *endInstance) {
	return endInstance->endInstanceContents->event;
}

End *endInstance_getEnd(EndInstance *endInstance) {
	return endInstance->end;
}

AtomInstance *endInstance_getAtomInstance(EndInstance *endInstance) {
	return endInstance_getOrientation(endInstance) ?
			endInstance->endInstanceContents->atomInstance :
	(endInstance->endInstanceContents->atomInstance != NULL ? atomInstance_getReverse(endInstance->endInstanceContents->atomInstance) : NULL);
}

int32_t endInstance_getCoordinate(EndInstance *endInstance) {
	return endInstance->endInstanceContents->coordinate;
}

int32_t endInstance_getStrand(EndInstance *endInstance) {
	return endInstance_getOrientation(endInstance) ? endInstance->endInstanceContents->strand : !endInstance->endInstanceContents->strand;
}

int32_t endInstance_getSide(EndInstance *endInstance) {
	return endInstance_getOrientation(endInstance) ? endInstance->endInstanceContents->side : !endInstance->endInstanceContents->side;
}

Sequence *endInstance_getSequence(EndInstance *endInstance) {
	return endInstance->endInstanceContents->sequence;
}

void endInstance_makeAdjacent1(EndInstance *endInstance, EndInstance *endInstance2) {
	endInstance_breakAdjacency1(endInstance);
	endInstance_breakAdjacency1(endInstance2);
	endInstance->endInstanceContents->adjacency = endInstance2;
	endInstance2->endInstanceContents->adjacency = endInstance;
}

void endInstance_makeAdjacent2(EndInstance *endInstance, EndInstance *endInstance2) {
	endInstance_breakAdjacency2(endInstance);
	endInstance_breakAdjacency2(endInstance2);
	endInstance->endInstanceContents->adjacency2 = endInstance2;
	endInstance2->endInstanceContents->adjacency2 = endInstance;
}

EndInstance *endInstance_getP(EndInstance *endInstance, EndInstance *connectedEndInstance) {
	return connectedEndInstance == NULL ? NULL :
			endInstance_getOrientation(endInstance) ? connectedEndInstance : endInstance_getReverse(connectedEndInstance);
}

EndInstance *endInstance_getAdjacency(EndInstance *endInstance) {
	return endInstance_getP(endInstance, endInstance->endInstanceContents->adjacency);
}

EndInstance *endInstance_getAdjacency2(EndInstance *endInstance) {
	return endInstance_getP(endInstance, endInstance->endInstanceContents->adjacency2);
}

Operation *endInstance_getOperation(EndInstance *endInstance) {
	return endInstance->endInstanceContents->operation;
}

EndInstance *endInstance_getParent(EndInstance *endInstance) {
	EndInstance *e = endInstance->endInstanceContents->parent;
	return e == NULL ? NULL :
		endInstance_getOrientation(endInstance) ? e : endInstance_getReverse(e);
}

int32_t endInstance_getChildNumber(EndInstance *endInstance) {
	return endInstance->endInstanceContents->children->length;
}

EndInstance *endInstance_getChild(EndInstance *endInstance, int32_t index) {
#ifdef BEN_DEBUG
	assert(endInstance_getChildNumber(endInstance) > index);
	assert(index >= 0);
#endif
	return endInstance_getP(endInstance, endInstance->endInstanceContents->children->list[index]);
}

void endInstance_makeParentAndChild(EndInstance *endInstanceParent, EndInstance *endInstanceChild) {
	assert(endInstance_getOrientation(endInstanceParent) == endInstance_getOrientation(endInstanceChild));
	if(!listContains(endInstanceParent->endInstanceContents->children, endInstanceChild)) { //defensive, means second calls will have no effect.
		listAppend(endInstanceParent->endInstanceContents->children, endInstanceChild);
	}
	endInstanceChild->endInstanceContents->parent = endInstanceParent;
}

int32_t endInstance_isInternal(EndInstance *endInstance) {
	return endInstance_getChildNumber(endInstance) > 0;
}

int32_t endInstance_isAugmented(EndInstance *endInstance) {
	return end_getAtom(endInstance_getEnd(endInstance)) != NULL && endInstance_getAtomInstance(endInstance) == NULL;
}

/*
 * Private functions.
 */

void endInstance_setAtomInstance(EndInstance *endInstance, AtomInstance *atomInstance) {
	endInstance->endInstanceContents->atomInstance = atomInstance;
}

void endInstance_setOperation(EndInstance *endInstance, Operation *operation) {
	endInstance->endInstanceContents->operation = operation;
}

void endInstance_breakAdjacency1(EndInstance *endInstance) {
	EndInstance *endInstance2;
	endInstance2 = endInstance_getAdjacency(endInstance);
	if(endInstance2 != NULL) {
		endInstance2->endInstanceContents->adjacency = NULL;
		endInstance->endInstanceContents->adjacency = NULL;
	}
}

void endInstance_breakAdjacency2(EndInstance *endInstance) {
	EndInstance *endInstance2;
	endInstance2 = endInstance_getAdjacency2(endInstance);
	if(endInstance2 != NULL) {
		endInstance2->endInstanceContents->adjacency2 = NULL;
		endInstance->endInstanceContents->adjacency2 = NULL;
	}
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic end functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t end_constructP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(endInstance_getInstanceName((EndInstance *)o1), endInstance_getInstanceName((EndInstance *)o2));
}

End *end_construct(const char *name, Net *net) {
	End *end;
	end = malloc(sizeof(End));
	end->rEnd = malloc(sizeof(End));
	end->rEnd->rEnd = end;
	end->endContents = malloc(sizeof(EndContents));
	end->rEnd->endContents = end->endContents;

	end->orientation = 1;
	end->rEnd->orientation = -1;

	end->endContents->rootInstance = NULL;
	end->endContents->name = stringCopy(name);
	end->endContents->endInstances = sortedSet_construct(end_constructP);
	end->endContents->attachedAtom = NULL;
	end->endContents->adjacencyComponent = NULL;
	end->endContents->net = net;
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
			endInstance_construct2(endInstance_getInstanceName(endInstance), end2,
					endInstance_getCoordinate(endInstance), endInstance_getStrand(endInstance),
					endInstance_getSide(endInstance), endInstance_getSequence(endInstance));
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
			endInstance_makeParentAndChild(end_getInstance(end2, endInstance_getInstanceName(endInstance2)),
										   end_getInstance(end2, endInstance_getInstanceName(endInstance)));
		}
	}
	end_destructInstanceIterator(iterator);

	//Copy root.
	end_setRootInstance(end2, end_getRootInstance(end));
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
	sortedSet_destruct(end->endContents->endInstances, NULL);

	free(end->endContents->name);
	free(end->endContents);
	free(end->rEnd);
	free(end);
}

const char *end_getName(End *end) {
	return end->endContents->name;
}

int32_t end_getOrientation(End *end) {
	return end->orientation;
}

char *end_getNameWithOrientation(End *end) {
	return netMisc_getNameWithOrientation(end_getName(end), end_getOrientation(end));
}

End *end_getReverse(End *end) {
	return end->rEnd;
}

Net *end_getNet(End *end) {
	return end->endContents->net;
}

Atom *end_getAtom(End *end) {
	Atom *a = end->endContents->attachedAtom;
	return a == NULL || end_getOrientation(end) ? a : atom_getReverse(a);
}

AdjacencyComponent *end_getAdjacencyComponent(End *end) {
	return end->endContents->adjacencyComponent;
}

int32_t end_getInstanceNumber(End *end) {
	return sortedSet_getLength(end->endContents->endInstances);
}

EndInstance *end_getInstanceP(End *end, EndInstance *connectedEndInstance) {
	return end_getOrientation(end) ? connectedEndInstance : connectedEndInstance;
}

EndInstance *end_getInstance(End *end, const char *name) {
	static EndInstance endInstance;
	static EndInstanceContents endInstanceContents;
	endInstance.endInstanceContents = &endInstanceContents;
	endInstanceContents.instance = (char *)name;
	return end_getInstanceP(end, sortedSet_find(end->endContents->endInstances, &endInstance));
}

EndInstance *end_getFirst(End *end) {
	return end_getInstanceP(end, sortedSet_getFirst(end->endContents->endInstances));
}

EndInstance *end_getRootInstance(End *end) {
	return end_getInstanceP(end, end->endContents->rootInstance);
}

void end_setRootInstance(End *end, EndInstance *endInstance) {
	end->endContents->rootInstance = endInstance_getOrientation(endInstance) ? endInstance : endInstance_getReverse(endInstance);
}

End_InstanceIterator *end_getInstanceIterator(End *end) {
	End_InstanceIterator *iterator;
	iterator = malloc(sizeof(struct _end_instanceIterator));
	iterator->end = end;
	iterator->iterator = iterator_construct(end->endContents->endInstances);
	return iterator;
}

EndInstance *end_getNext(End_InstanceIterator *iterator) {
	return end_getInstanceP(iterator->end, iterator_getNext(iterator->iterator));
}

EndInstance *end_getPrevious(End_InstanceIterator *iterator) {
	return end_getInstanceP(iterator->end, iterator_getPrevious(iterator->iterator));
}

End_InstanceIterator *end_copyInstanceIterator(End_InstanceIterator *iterator) {
	End_InstanceIterator *iterator2;
	iterator2 = malloc(sizeof(struct _end_instanceIterator));
	iterator2->end = iterator->end;
	iterator2->iterator = iterator_copy(iterator->iterator);
	return iterator2;
}

void end_destructInstanceIterator(End_InstanceIterator *iterator) {
	iterator_destruct(iterator->iterator);
	free(iterator);
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
	assert(end_getOrientation(end) == endInstance_getOrientation(endInstance));
	sortedSet_insert(end->endContents->endInstances, endInstance);
}

void end_removeInstance(End *end, EndInstance *endInstance) {
	assert(end_getOrientation(end) == endInstance_getOrientation(endInstance));
	sortedSet_delete(end->endContents->endInstances, endInstance);
}

void end_setAdjacencyComponent(End *end, AdjacencyComponent *adjacencyComponent) {
	//argument may be NULL
	end->endContents->adjacencyComponent = adjacencyComponent;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic atom instance functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

AtomInstance *atomInstance_construct(Atom *atom,
		EndInstance *leftEndInstance, EndInstance *rightEndInstance) {
	AtomInstance *atomInstance;
	atomInstance = malloc(sizeof(AtomInstance));
	atomInstance->rInstance = malloc(sizeof(AtomInstance));
	atomInstance->rInstance->rInstance = atomInstance;
	atomInstance->atom = atom;
	atomInstance->rInstance->atom = atom_getReverse(atom);
	atomInstance->_5EndInstance = leftEndInstance;
	atomInstance->rInstance->_5EndInstance = rightEndInstance;
	endInstance_setAtomInstance(atomInstance->_5EndInstance, atomInstance);
	endInstance_setAtomInstance(atomInstance->rInstance->_5EndInstance, atomInstance->rInstance);
	atom_addInstance(atom, atomInstance);
	return atomInstance;
}

AtomInstance *atomInstance_construct2(const char *instance, Atom *atom) {
	return atomInstance_construct(atom,
			endInstance_construct(instance, atom_get5End(atom)),
			endInstance_construct(instance, atom_get3End(atom)));
}

AtomInstance *atomInstance_construct3(const char *instance, Atom *atom,
		int32_t startCoordinate, int32_t strand, Sequence *sequence) {
#ifdef BEN_DEBUG
	assert(startCoordinate >= 0);
	assert(startCoordinate + atom_getLength(atom) <= sequence_getLength(sequence));
#endif

	return atomInstance_construct(atom,
			endInstance_construct2(instance, atom_get5End(atom),
					startCoordinate, strand, 1, sequence),
			endInstance_construct2(instance, atom_get3End(atom),
						startCoordinate + atom_getLength(atom) - 1, strand, -1, sequence));
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
	return endInstance_getInstanceName(atomInstance_get5End(atomInstance));
}

const char *atomInstance_getElementName(AtomInstance *atomInstance) {
	return atom_getName(atomInstance_getAtom(atomInstance));
}

char *atomInstance_getCompleteName(AtomInstance *atomInstance) {
	return netMisc_makeCompleteName(atomInstance_getInstanceName(atomInstance), atomInstance_getElementName(atomInstance), 1);
}

int32_t atomInstance_getOrientation(AtomInstance *atomInstance) {
	return atom_getOrientation(atomInstance_getAtom(atomInstance));
}

char *atomInstance_getCompleteNameWithOrientation(AtomInstance *atomInstance) {
	return netMisc_makeCompleteName(atomInstance_getInstanceName(atomInstance), atomInstance_getElementName(atomInstance), atomInstance_getOrientation(atomInstance));
}

AtomInstance *atomInstance_getReverse(AtomInstance *atomInstance) {
	return atomInstance->rInstance;
}

Event *atomInstance_getEvent(AtomInstance *atomInstance) {
	return endInstance_getEvent(atomInstance_get5End(atomInstance));
}

int32_t atomInstance_getStart(AtomInstance *atomInstance) {
	return endInstance_getCoordinate(atomInstance_get5End(atomInstance));
}

int32_t atomInstance_getStrand(AtomInstance *atomInstance) {
	return endInstance_getStrand(atomInstance_get5End(atomInstance));
}

int32_t atomInstance_getLength(AtomInstance *atomInstance) {
	return atom_getLength(atomInstance_getAtom(atomInstance));
}

Sequence *atomInstance_getSequence(AtomInstance *atomInstance) {
	return endInstance_getSequence(atomInstance_get5End(atomInstance));
}

EndInstance *atomInstance_get5End(AtomInstance *atomInstance) {
	return atomInstance->_5EndInstance;
}

EndInstance *atomInstance_get3End(AtomInstance *atomInstance) {
	return atomInstance_get5End(atomInstance_getReverse(atomInstance));
}

AtomInstance *atomInstance_getParent(AtomInstance *atomInstance) {
	EndInstance *endInstance;
	AtomInstance *atomInstance2;
	endInstance = atomInstance_get5End(atomInstance);
	while((endInstance = endInstance_getParent(endInstance)) != NULL) {
		if((atomInstance2 = endInstance_getAtomInstance(endInstance)) != NULL) {
			return atomInstance2;
		}
	}
	return NULL;
}

int32_t atomInstance_getChildNumber(AtomInstance *atomInstance) {
	return endInstance_getChildNumber(atomInstance_get5End(atomInstance));
}

AtomInstance *atomInstance_getChild(AtomInstance *atomInstance, int32_t index) {
	EndInstance *endInstance;
	AtomInstance *atomInstance2;
	endInstance = endInstance_getChild(atomInstance_get5End(atomInstance), index);
	while(endInstance_isAugmented(endInstance)) {
		assert(endInstance_getChildNumber(endInstance) == 1);
		endInstance = endInstance_getChild(endInstance, 0);
	}
	atomInstance2 = endInstance_getAtomInstance(endInstance);
	assert(atomInstance_getOrientation(atomInstance) == atomInstance_getOrientation(atomInstance2));
	return atomInstance2;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic atom functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t atomConstruct_constructP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(atomInstance_getInstanceName((AtomInstance *)o1), atomInstance_getInstanceName((AtomInstance *)o2));
}

Atom *atom_construct(const char *name, int32_t length, Net *net) {
	Atom *atom;
	atom = malloc(sizeof(Atom));
	atom->rAtom = malloc(sizeof(Atom));
	atom->rAtom->rAtom = atom;
	atom->atomContents = malloc(sizeof(AtomContents));
	atom->rAtom->atomContents = atom->atomContents;

	atom->orientation = 1;
	atom->rAtom->orientation = -1;

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

int32_t atom_getOrientation(Atom *atom) {
	return atom->orientation;
}

char *atom_getNameWithOrientation(Atom *atom) {
	return netMisc_getNameWithOrientation(atom_getName(atom), atom_getOrientation(atom));
}

Atom *atom_getReverse(Atom *atom) {
	return atom->rAtom;
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

End *atom_get5End(Atom *atom) {
	return atom->_5End;
}

End *atom_get3End(Atom *atom) {
	return atom->rAtom->_5End;
}

int32_t atom_getInstanceNumber(Atom *atom) {
	return sortedSet_getLength(atom->atomContents->atomInstances);
}

AtomInstance *atom_getInstanceP(Atom *atom, AtomInstance *connectedAtomInstance) {
	return atom_getOrientation(atom) || connectedAtomInstance == NULL ? connectedAtomInstance : atomInstance_getReverse(connectedAtomInstance);
}

AtomInstance *atom_getInstance(Atom *atom, const char *name) {
	static AtomInstance atomInstance;
	static EndInstance endInstance;
	static EndInstanceContents endInstanceContents;
	endInstanceContents.instance = (char *)name;
	endInstance.endInstanceContents = &endInstanceContents;
	atomInstance._5EndInstance = &endInstance;
	return atom_getInstanceP(atom, sortedSet_find(atom->atomContents->atomInstances, &atomInstance));
}

AtomInstance *atom_getFirst(Atom *atom) {
	return atom_getInstanceP(atom, sortedSet_getFirst(atom->atomContents->atomInstances));
}

Atom_InstanceIterator *atom_getInstanceIterator(Atom *atom) {
	Atom_InstanceIterator *iterator;
	iterator = malloc(sizeof(struct _atom_instanceIterator));
	iterator->atom = atom;
	iterator->iterator = iterator_construct(atom->atomContents->atomInstances);
	return iterator;
}

AtomInstance *atom_getNext(Atom_InstanceIterator *iterator) {
	return atom_getInstanceP(iterator->atom, iterator_getNext(iterator->iterator));
}

AtomInstance *atom_getPrevious(Atom_InstanceIterator *iterator) {
	return atom_getInstanceP(iterator->atom, iterator_getPrevious(iterator->iterator));
}

Atom_InstanceIterator *atom_copyInstanceIterator(Atom_InstanceIterator *iterator) {
	Atom_InstanceIterator *iterator2;
	iterator2 = malloc(sizeof(struct _atom_instanceIterator));
	iterator2->atom = iterator->atom;
	iterator2->iterator = iterator_copy(iterator->iterator);
	return iterator2;
}

void atom_destructInstanceIterator(Atom_InstanceIterator *iterator) {
	iterator_destruct(iterator->iterator);
	free(iterator);
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

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic adjacency component functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t adjacencyComponent_constructP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(end_getName((End *)o1), end_getName((End *)o2));
}

AdjacencyComponent *adjacencyComponent_construct(Net *net, Net *nestedNet) {
	AdjacencyComponent *adjacencyComponent;

	adjacencyComponent = adjacencyComponent_construct2(net, net_getName(nestedNet));
	adjacencyComponent_updateContainedEnds(adjacencyComponent);
	return adjacencyComponent;
}

void adjacencyComponent_destructP(End *end, void *o) {
	assert(o == NULL);
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
			adjacencyComponent_addEnd(adjacencyComponent, end);
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
	static EndContents endContents;
	end.endContents = &endContents;
	endContents.name = (char *)name;
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

AdjacencyComponent *adjacencyComponent_construct2(Net *net, const char *nestedNetName) {
	AdjacencyComponent *adjacencyComponent;
	adjacencyComponent = malloc(sizeof(AdjacencyComponent));

	adjacencyComponent->net = net;
	adjacencyComponent->chain = NULL;
	adjacencyComponent->nestedNetName = stringCopy(nestedNetName);
	adjacencyComponent->ends = sortedSet_construct(adjacencyComponent_constructP);
	net_addAdjacencyComponent(net, adjacencyComponent);

	return adjacencyComponent;
}

void adjacencyComponent_addEnd(AdjacencyComponent *adjacencyComponent, End *end) {
	sortedSet_insert(adjacencyComponent->ends, end);
	end_setAdjacencyComponent(end, adjacencyComponent);
}

void adjacencyComponent_setChain(AdjacencyComponent *adjacencyComponent, Chain *chain) {
	//argument may be NULL
	adjacencyComponent->chain = chain;
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

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic chain functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Chain *chain_construct(Net *net, const char *name) {
	Chain *chain;
	chain = malloc(sizeof(Chain));
	chain->name = stringCopy(name);
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
	free(chain->name);
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

const char *chain_getName(Chain *chain) {
	return chain->name;
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

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic operation functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Operation *operation_construct(Net *net, const char *name) {
	Operation *operation;
	operation = malloc(sizeof(Operation));
	operation->name = stringCopy(name);
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

const char *operation_getName(Operation *operation) {
	return operation->name;
}

/*
 * Private functions
 */
void operation_setIndex(Operation *operation, int32_t index) {
	operation->index = index;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t net_constructSequencesP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(sequence_getName((Sequence *)o1), sequence_getName((Sequence *)o2));
}

int32_t net_constructEndsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(end_getName((End *)o1), end_getName((End *)o2));
}

int32_t net_constructAtomsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(atom_getName((Atom *)o1), atom_getName((Atom *)o2));
}

int32_t net_constructAdjacencyComponentsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(adjacencyComponent_getNestedNetName((AdjacencyComponent *)o1),
			adjacencyComponent_getNestedNetName((AdjacencyComponent *)o2));
}

int32_t net_constructChainsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(chain_getName((Chain *)o1), chain_getName((Chain *)o2));
}

int32_t net_constructOperationsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(operation_getName((Operation *)o1), operation_getName((Operation *)o2));
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


void net_destruct(Net *net, int32_t recursive) {
	Net_AdjacencyComponentIterator *iterator;
	Sequence *sequence;
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

	eventTree_destruct(net_getEventTree(net));

	sortedSet_destruct(net->sequences, NULL);
	while((sequence = net_getFirstSequence(net)) != NULL) {
		sequence_destruct(sequence);
	}

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

const char *net_getName(Net *net) {
	return net->name;
}

NetDisk *net_getNetDisk(Net *net) {
	return net->netDisk;
}

EventTree *net_getEventTree(Net *net) {
	return net->eventTree;
}

Sequence *net_getFirstSequence(Net *net) {
	return sortedSet_getFirst(net->sequences);
}

Sequence *net_getSequence(Net *net, const char *name) {
	static Sequence sequence;
	static MetaSequence metaSequence;
	sequence.metaSequence = &metaSequence;
	metaSequence.name = (char *)name;
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
	static EndContents endContents;
	end.endContents = &endContents;
	endContents.name = (char *)name;
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
	static AtomContents atomContents;
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

Chain *net_getChain(Net *net, const char *name) {
	static Chain chain;
	chain.name = (char *)name;
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

Operation *net_getOperation(Net *net, const char *name) {
	static Operation operation;
	operation.name = (char *)name;
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

/*
 * Private functions
 */

void net_addEventTree(Net *net, EventTree *eventTree) {
	net->eventTree = eventTree;
}

void net_addSequence(Net *net, Sequence *sequence) {
	sortedSet_insert(net->sequences, sequence);
}

void net_removeSequence(Net *net, Sequence *sequence) {
	sortedSet_delete(net->sequences, sequence);
}

void net_addEnd(Net *net, End *end) {
	sortedSet_insert(net->ends, end);
}

void net_removeEnd(Net *net, End *end) {
	sortedSet_delete(net->ends, end);
}

void net_addAtom(Net *net, Atom *atom) {
	sortedSet_insert(net->atoms, atom);
}

void net_removeAtom(Net *net, Atom *atom) {
	sortedSet_delete(net->atoms, atom);
}

void net_addChain(Net *net, Chain *chain) {
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
	sortedSet_insert(net->operations, operation);
}

void net_removeOperation(Net *net, Operation *operation) {
	sortedSet_delete(net->operations, operation);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t netDisk_constructNetsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(net_getName((Net *)o1), net_getName((Net *)o2));
}

int32_t netDisk_constructMetaSequencesP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return strcmp(metaSequence_getName((MetaSequence *)o1), metaSequence_getName((MetaSequence *)o2));
}

NetDisk *netDisk_construct(const char *netDiskFile) {
	NetDisk *netDisk;
	netDisk = malloc(sizeof(NetDisk));

	//construct lists of in memory objects
	netDisk->metaSequences = sortedSet_construct(netDisk_constructMetaSequencesP);
	netDisk->nets = sortedSet_construct(netDisk_constructNetsP);

	//the files to write the databases in
	netDisk->metaSequencesDatabaseName = pathJoin(netDiskFile, "sequences");
	netDisk->netsDatabaseName = pathJoin(netDiskFile, "nets");
	netDisk->iDDatabaseName = pathJoin(netDiskFile, "uniqueIDs");

	//open the sequences database
	netDisk->metaSequencesDatabase = database_construct(netDisk->metaSequencesDatabaseName);
	netDisk->netsDatabase = database_construct(netDisk->netsDatabaseName);
	netDisk->iDDatabase = database_construct(netDisk->iDDatabaseName);

	//construct the string file
	netDisk->stringFile = pathJoin(netDiskFile, "strings");
	netDisk->stringFileLength = 0;

	//initialise the unique ids.
	netDisk_getUniqueID(netDisk);

	return netDisk;
}

void netDisk_destruct(NetDisk *netDisk){
	Net *net;

	while((net = netDisk_getFirstNetInMemory(netDisk)) != NULL) {
		net_destruct(net, FALSE);
	}
	sortedSet_destruct(netDisk->nets, NULL);

	sortedSet_destruct(netDisk->metaSequences, NULL);

	//close DBs
	database_destruct(netDisk->metaSequencesDatabase);
	database_destruct(netDisk->netsDatabase);

	//free string names
	free(netDisk->metaSequencesDatabaseName);
	free(netDisk->netsDatabaseName);

	free(netDisk);
}

int32_t netDisk_write(NetDisk *netDisk){
	NetDisk_NetIterator *netIterator;
	struct avl_traverser *metaSequenceIterator;
	char *cA;
	Net *net;
	MetaSequence *metaSequence;
	int32_t ecode;

	netIterator = netDisk_getNetInMemoryIterator(netDisk);
	while((net = netDisk_getNextNet(netIterator)) != NULL) {
		cA = netMisc_makeBinaryRepresentation(net,
				(void (*)(void *, void (*)(const char *string, ...)))net_writeBinaryRepresentation);
		if((ecode = database_writeRecord(netDisk->netsDatabase, net_getName(net), cA)) != 0) {
			return ecode;
		}
		free(cA);
	}
	netDisk_destructNetIterator(netIterator);

	metaSequenceIterator = iterator_construct(netDisk->metaSequences);
	while((metaSequence = iterator_getNext(metaSequenceIterator)) != NULL) {
		cA = netMisc_makeBinaryRepresentation(metaSequence,
				(void (*)(void *, void (*)(const char *string, ...)))metaSequence_writeBinaryRepresentation);
		if((ecode = database_writeRecord(netDisk->metaSequencesDatabase, metaSequence_getName(metaSequence), cA)) != 0) {
			return ecode;
		}
		free(cA);
	}
	iterator_destruct(metaSequenceIterator);
	return 0;
}

void *netDisk_getObject(NetDisk *netDisk, TCBDB *database, void *(*getObjectInMemory)(NetDisk *, const char *),
		void *(*loadFromBinaryRepresentation)(char **, NetDisk *), const char *objectName) {
	char *cA;
	char *cA2;
	void *object;

	//try in memory list first.
	if((object = getObjectInMemory(netDisk, objectName)) != NULL) {
		return object;
	}
	//else try the database.
	cA = database_getRecord(database, objectName);
	if(cA == NULL) {
		return NULL;
	}
	else {
		cA2 = cA;
		object = loadFromBinaryRepresentation(&cA2, netDisk);
		free(cA);
		return object;
	}
}

Net *netDisk_getNet(NetDisk *netDisk, const char *netName) {
	return netDisk_getObject(netDisk, netDisk->netsDatabase,
			(void *(*)(NetDisk *, const char *))netDisk_getNetInMemory,
			(void *(*)(char **, NetDisk *))net_loadFromBinaryRepresentation, netName);
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

void netDisk_addNet(NetDisk *netDisk, Net *net) {
	assert(netDisk_getNet(netDisk, net_getName(net)) == NULL);
	sortedSet_insert(netDisk->nets, net);
}

int32_t netDisk_deleteNetFromDisk(NetDisk *netDisk, const char *netName) {
	return database_removeRecord(netDisk->netsDatabase, netName);
}

void netDisk_unloadNet(NetDisk *netDisk, Net *net) {
	assert(netDisk_getNetInMemory(netDisk, net_getName(net)) != NULL);
	sortedSet_delete(netDisk->nets, net);
}

void netDisk_addMetaSequence(NetDisk *netDisk, MetaSequence *metaSequence) {
#ifdef BEN_DEBUG
	assert(netDisk_getMetaSequence(netDisk, metaSequence_getName(metaSequence)) == NULL);
#endif
	sortedSet_insert(netDisk->metaSequences, metaSequence);
}

int32_t netDisk_deleteMetaSequenceFromDisk(NetDisk *netDisk, const char *metaSequenceName) {
	return database_removeRecord(netDisk->metaSequencesDatabase, metaSequenceName);
}

void netDisk_unloadMetaSequence(NetDisk *netDisk, MetaSequence *metaSequence) {
	sortedSet_delete(netDisk->metaSequences, metaSequence);
}

MetaSequence *netDisk_getMetaSequenceInMemory(NetDisk *netDisk, const char *metaSequenceName) {
	static MetaSequence metaSequence;
	metaSequence.name = (void *)metaSequenceName;
	return sortedSet_find(netDisk->metaSequences, &metaSequence);
}

MetaSequence *netDisk_getMetaSequence(NetDisk *netDisk, const char *metaSequenceName) {
return netDisk_getObject(netDisk, netDisk->netsDatabase,
		(void *(*)(NetDisk *, const char *))netDisk_getMetaSequenceInMemory,
		(void *(*)(char **, NetDisk *))metaSequence_loadFromBinaryRepresentation, metaSequenceName);
}

int64_t netDisk_addString(NetDisk *netDisk, const char *string, int32_t length) {
	int64_t fileOffset;
	FILE *fileHandle;

	assert(length == (int32_t)strlen(string));
	fileOffset = netDisk->stringFileLength;
	netDisk->stringFileLength += length;
	fileHandle = fopen(netDisk->stringFile, "a");
	fprintf(fileHandle, "%s", string);
	fclose(fileHandle);
	return fileOffset;
}

char *netDisk_getString(NetDisk *netDisk, int64_t offset, int32_t start, int32_t length, int32_t strand) {
	FILE *fileHandle;
	char *cA;
	char *cA2;

	fileHandle = fopen(netDisk->stringFile, "r");
	fseek(fileHandle, offset+start, SEEK_SET);
	cA = malloc(sizeof(char)*(length+1));
	fread(cA, sizeof(char), length, fileHandle);
	cA[length] = '\0';
	fclose(fileHandle);

	if(!strand) {
		cA2 = netMisc_reverseComplementString(cA);
		free(cA);
		return cA2;
	}
	return cA;
}

void netDisk_getBlockOfUniqueIDs(NetDisk *netDisk) {
	char *cA;
	char *cA2;
	cA = database_getRecord(netDisk->iDDatabase, "uniqueID");
	if(cA == NULL) {
		database_writeRecord(netDisk->iDDatabase, "uniqueID", "0");
		cA = database_getRecord(netDisk->iDDatabase, "uniqueID");
	}
	cA2 = cA;
	database_
	netDisk->uniqueNumber = binaryRepresentation_get64BitInteger(&cA2);
	netDisk->uniqueNumber = database_getRecord(netDisk->);
	netDisk->maxUniqueNumber = netDisk->uniqueNumber + MAX_INCREMENT_NUMBER;
	free(cA);
}

int64_t netDisk_getUniqueID(NetDisk *netDisk) {
	assert(netDisk->uniqueNumber <= netDisk->maxUniqueNumber);
	if(netDisk->uniqueNumber == netDisk->maxUniqueNumber) {
		netDisk_getBlockOfUniqueIDs(netDisk);
	}
	return netDisk->uniqueNumber++;
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

char *netMisc_makeCompleteName(const char *elementName, const char *instanceName, int32_t orientation) {
	char *cA;
	if(orientation > 0) {
		cA = malloc(sizeof(char) * (strlen(elementName) + strlen(instanceName) + 2));
		sprintf(cA, "%s.%s", elementName, instanceName);
	}
	else {
		cA = malloc(sizeof(char) * (strlen(elementName) + strlen(instanceName) + 3));
		sprintf(cA, "-%s.%s", elementName, instanceName);
	}
	return cA;
}

static char *netMisc_makeCompleteNameStatic_completeName = NULL;
const char *netMisc_makeCompleteNameStatic(const char *elementName, const char *instanceName, int32_t orientation) {
	if(netMisc_makeCompleteNameStatic_completeName != NULL) {
		free(netMisc_makeCompleteNameStatic_completeName);
	}
	netMisc_makeCompleteNameStatic_completeName = netMisc_makeCompleteName(elementName, instanceName, orientation);
	return netMisc_makeCompleteNameStatic_completeName;
}

char *netMisc_getNameWithOrientation(const char *name, int32_t orientation) {
	char *cA;
	assert(strlen(name) > 0);
	assert(orientation != 0);
	if(name[0] == '-') {
		assert(orientation < 0);
		return stringCopy(name);
	}
	if(orientation > 1) {
		return stringCopy(name);
	}
	cA = malloc(sizeof(char) * (strlen(name) + 2));
	sprintf(cA, "-%s", name);
	return cA;
}

/*
 * Private utility functions.
 */

char netMisc_makeBinaryRepresentationP_cA[100000];
int32_t netMisc_makeBinaryRepresentationP_i = 0;
void netMisc_makeBinaryRepresentationP(const char *string, ...) {
	/*
	 * Records the cummulative size of the substrings written out in creating the net.
	 */
	//return;
    va_list ap;
    va_start(ap, string);
    sprintf(netMisc_makeBinaryRepresentationP_cA, string, ap);
    netMisc_makeBinaryRepresentationP_i += strlen(netMisc_makeBinaryRepresentationP_cA);
    va_end(ap);
}

char *netMisc_makeBinaryRepresentationP2_cA = NULL;
void netMisc_makeBinaryRepresentationP2(const char *string, ...) {
	/*
	 * Cummulates all the sequences into one.
	 */
	//return;
    va_list ap;
    va_start(ap, string);
    sprintf(netMisc_makeBinaryRepresentationP2_cA, string, ap);
    netMisc_makeBinaryRepresentationP2_cA += strlen(netMisc_makeBinaryRepresentationP2_cA);
    va_end(ap);
}

char *netMisc_makeBinaryRepresentation(void *object, void (*writeBinaryRepresentation)(void *, void (*writeFn)(const char *string, ...))) {
	char *cA;
	netMisc_makeBinaryRepresentationP_i = 0;
	writeBinaryRepresentation(object, netMisc_makeBinaryRepresentationP);
	cA = malloc(sizeof(char)*(netMisc_makeBinaryRepresentationP_i + 1));
	netMisc_makeBinaryRepresentationP2_cA = cA;
	net_writeBinaryRepresentation(object, netMisc_makeBinaryRepresentationP2);
	return cA;
}

char netMisc_reverseComplementChar(char c) {
	switch(c) {
		case 'a':
			return 't';
		case 'c':
			return 'g';
		case 'g':
			return 'c';
		case 't':
			return 'a';
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		default:
			return c;
	}
}

char *netMisc_reverseComplementString(const char *string) {
	int32_t i, j;

	j = strlen(string);
	char *cA;

	cA = malloc(sizeof(string) *(j+1));
	for(i=0; i<j; i++) {
		cA[i] = netMisc_reverseComplementChar(string[j-1-i]);
	}
	return cA;
}

