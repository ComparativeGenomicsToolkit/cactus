#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "net.h"
#include "netPrivate.h"
#include "bioioC.h"

/*
 * Implementation of code to interface with the C code.
 */

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions for serialising the objects.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void binaryRepresentation_writeElementType(char elementCode, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	writeFn(&elementCode, sizeof(char), 1);
}

void binaryRepresentation_writeString(const char *name, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	int32_t i;
	i = strlen(name);
	writeFn(&i, sizeof(int32_t), 1);
	writeFn(&name, sizeof(char), i);
}

void binaryRepresentation_writeInteger(int32_t i, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	writeFn(&i, sizeof(int32_t), 1);
}

void binaryRepresentation_write64BitInteger(int64_t i, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	writeFn(&i, sizeof(int64_t), 1);
}

void binaryRepresentation_writeFloat(float f, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	writeFn(&f, sizeof(float), 1);
}

char binaryRepresentation_peekNextElementType(void *binaryString) {
	return *((char *)binaryString);
}

char binaryRepresentation_popNextElementType(void **binaryString) {
	return binaryRepresentation_getChar(binaryString);
}

char binaryRepresentation_getChar(void **binaryString) {
	char *c;
	c = *binaryString;
	*binaryString = c + 1;
	return *c;
}

char *binaryRepresentation_getString(void **binaryString) {
	int32_t i, j;
	char *cA;
	i = binaryRepresentation_getInteger(binaryString);
	cA = malloc(sizeof(char)*i);
	for(j=0; j<i; j++) {
		cA[j] = binaryRepresentation_popNextElementType(binaryString);
	}
	return cA;
}

char *binaryRepresentation_getStringStatic_cA = NULL;
const char *binaryRepresentation_getStringStatic(void **binaryString) {
	if(binaryRepresentation_getStringStatic_cA != NULL) {
		free(binaryRepresentation_getStringStatic_cA);
	}
	binaryRepresentation_getStringStatic_cA = binaryRepresentation_getString(binaryString);
	return binaryRepresentation_getStringStatic_cA;
}

int32_t binaryRepresentation_getInteger(void **binaryString) {
	int32_t *i;
	i = *binaryString;
	*binaryString = i + 1;
	return *i;
}

int64_t binaryRepresentation_get64BitInteger(void **binaryString) {
	int64_t *i;
	i = *binaryString;
	*binaryString = i + 1;
	return *i;
}

float binaryRepresentation_getFloat(void **binaryString) {
	float *i;
	i = *binaryString;
	*binaryString = i + 1;
	return *i;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void metaEvent_writeBinaryRepresentation(MetaEvent *metaEvent, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_META_EVENT, writeFn);
	binaryRepresentation_writeString(metaEvent_getName(metaEvent), writeFn);
	binaryRepresentation_writeString(metaEvent_getHeader(metaEvent), writeFn);
}

MetaEvent *metaEvent_loadFromBinaryRepresentation(void **binaryString, NetDisk *netDisk) {
	MetaEvent *metaEvent;
	char *name;
	char *header;

	metaEvent = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_META_EVENT) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getString(binaryString);
		header = binaryRepresentation_getString(binaryString);
		metaEvent = metaEvent_construct(name, header, netDisk);
		free(name);
		free(header);
	}
	return metaEvent;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void event_writeBinaryRepresentation(Event *event,
		void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_EVENT, writeFn);
	binaryRepresentation_writeString(event_getName(event_getParent(event)), writeFn);
	binaryRepresentation_writeString(event_getName(event), writeFn);
	binaryRepresentation_writeFloat(event_getBranchLength(event), writeFn);
}

Event *event_loadFromBinaryRepresentation(void **binaryString,
		EventTree *eventTree) {
	Event *event;
	char *parentName;
	float branchLength;

	event = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_EVENT) {
		binaryRepresentation_popNextElementType(binaryString);
		parentName = binaryRepresentation_getString(binaryString);
		branchLength = binaryRepresentation_getFloat(binaryString);
		event = event_construct(netDisk_getMetaEvent(net_getNetDisk(eventTree_getNet(eventTree)),
				binaryRepresentation_getStringStatic(binaryString)), branchLength,
				eventTree_getEvent(eventTree, parentName), eventTree);
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

void eventTree_writeBinaryRepresentationP(Event *event, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	int32_t i;
	event_writeBinaryRepresentation(event, writeFn);
	for(i=0; i<event_getChildNumber(event); i++) {
		eventTree_writeBinaryRepresentationP(event, writeFn);
	}
}

void eventTree_writeBinaryRepresentation(EventTree *eventTree, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	int32_t i;
	Event *event;
	event = eventTree_getRootEvent(eventTree);
	binaryRepresentation_writeElementType(CODE_EVENT_TREE, writeFn);
	binaryRepresentation_writeString(event_getName(event), writeFn);
	for(i=0; i<event_getChildNumber(event); i++) {
		eventTree_writeBinaryRepresentationP(event_getChild(event, i), writeFn);
	}
}

EventTree *eventTree_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	EventTree *eventTree;
	eventTree = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_EVENT_TREE) {
		binaryRepresentation_popNextElementType(binaryString);
		eventTree = eventTree_construct(netDisk_getMetaEvent(net_getNetDisk(eventTree_getNet(eventTree)),
				binaryRepresentation_getStringStatic(binaryString)), net);
		while(event_loadFromBinaryRepresentation(binaryString, eventTree) != NULL);
	}
	return eventTree;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic meta sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void metaSequence_writeBinaryRepresentation(MetaSequence *metaSequence,
		void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_META_SEQUENCE, writeFn);
	binaryRepresentation_writeInteger(metaSequence_getStart(metaSequence), writeFn);
	binaryRepresentation_writeInteger(metaSequence_getLength(metaSequence), writeFn);
	binaryRepresentation_writeString(metaSequence_getEventName(metaSequence), writeFn);
	binaryRepresentation_write64BitInteger(metaSequence_getFileOffset(metaSequence), writeFn);
	binaryRepresentation_writeString(metaSequence_getHeader(metaSequence), writeFn);
}

MetaSequence *metaSequence_loadFromBinaryRepresentation(void **binaryString,
		NetDisk *netDisk) {
	MetaSequence *metaSequence;
	char *name;
	int32_t start;
	int32_t length;
	int64_t fileOffset;
	char *eventName;
	char *header;

	metaSequence = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_META_SEQUENCE) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getString(binaryString);
		start = binaryRepresentation_getInteger(binaryString);
		length = binaryRepresentation_getInteger(binaryString);
		eventName = binaryRepresentation_getString(binaryString);
		fileOffset = binaryRepresentation_get64BitInteger(binaryString);
		header = binaryRepresentation_getString(binaryString);
		metaSequence = metaSequence_construct2(name, start, length,
				fileOffset, header, eventName, netDisk);
		free(name);
		free(eventName);
		free(header);
	}
	return metaSequence;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void sequence_writeBinaryRepresentation(Sequence *sequence, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_SEQUENCE, writeFn);
	binaryRepresentation_writeString(sequence_getName(sequence), writeFn);
}

Sequence *sequence_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	Sequence *sequence;

	sequence = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_SEQUENCE) {
		binaryRepresentation_popNextElementType(binaryString);
		sequence = sequence_construct(netDisk_getMetaSequence(net_getNetDisk(net), binaryRepresentation_getStringStatic(binaryString)), net);
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

void endInstance_writeBinaryRepresentationP(EndInstance *endInstance2, int32_t elementType, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	char *cA;
	binaryRepresentation_writeElementType(elementType, writeFn);
	cA = endInstance_getCompleteName(endInstance2);
	binaryRepresentation_writeString(cA, writeFn);
	free(cA);
}

void endInstance_writeBinaryRepresentation(EndInstance *endInstance, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
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
		binaryRepresentation_writeInteger(endInstance_getStrand(endInstance), writeFn);
		binaryRepresentation_writeInteger(endInstance_getSide(endInstance), writeFn);
		binaryRepresentation_writeString(sequence_getName(endInstance_getSequence(endInstance)), writeFn);
	}
	if((endInstance2 = endInstance_getAdjacency(endInstance)) != NULL) {
		endInstance_writeBinaryRepresentationP(endInstance2, CODE_ADJACENCY, writeFn);
	}
	if((endInstance2 = endInstance_getAdjacency2(endInstance)) != NULL) {
		endInstance_writeBinaryRepresentationP(endInstance2, CODE_ADJACENCY, writeFn);
	}
	if((endInstance2 = endInstance_getParent(endInstance)) != NULL) {
		endInstance_writeBinaryRepresentationP(endInstance2, CODE_PARENT, writeFn);
	}
}

int32_t endInstance_loadFromBinaryRepresentationP(EndInstance *endInstance, void **binaryString, void (*linkFn)(EndInstance *, EndInstance *)) {
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

EndInstance *endInstance_loadFromBinaryRepresentation(void **binaryString, End *end) {
	EndInstance *endInstance;
	char *name;
	int32_t coordinate;
	int32_t strand;
	int32_t side;
	Sequence *sequence;

	endInstance = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_END_INSTANCE) {
		binaryRepresentation_popNextElementType(binaryString);
		endInstance = endInstance_construct(binaryRepresentation_getStringStatic(binaryString), end);
	}
	else if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_END_INSTANCE_WITH_EVENT) {
		binaryRepresentation_popNextElementType(binaryString);
		endInstance = endInstance_construct(binaryRepresentation_getStringStatic(binaryString), end);
		endInstance_setEvent(endInstance, eventTree_getEvent(net_getEventTree(end_getNet(end)), binaryRepresentation_getStringStatic(binaryString)));
	}
	else if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_END_INSTANCE_WITH_COORDINATES) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getString(binaryString);
		coordinate = binaryRepresentation_getInteger(binaryString);
		strand = binaryRepresentation_getInteger(binaryString);
		side = binaryRepresentation_getInteger(binaryString);
		sequence = net_getSequence(end_getNet(end), binaryRepresentation_getStringStatic(binaryString));
		endInstance = endInstance_construct2(name, end, coordinate, strand, side, sequence);
		free(name);
	}
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_ADJACENCY) {
		endInstance_loadFromBinaryRepresentationP(endInstance, binaryString, endInstance_makeAdjacent1);
	}
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_ADJACENCY) {
		endInstance_loadFromBinaryRepresentationP(endInstance, binaryString, endInstance_makeAdjacent2);
	}
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_PARENT) {
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

void end_writeBinaryRepresentationP(EndInstance *endInstance, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	int32_t i;
	endInstance_writeBinaryRepresentation(endInstance, writeFn);
	for(i=0; i<endInstance_getChildNumber(endInstance); i++) {
		end_writeBinaryRepresentationP(endInstance_getChild(endInstance, i), writeFn);
	}
}

void end_writeBinaryRepresentation(End *end, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	End_InstanceIterator *iterator;
	EndInstance *endInstance;

	endInstance = end_getRootInstance(end);
	binaryRepresentation_writeElementType(endInstance == NULL ?
			CODE_END_WITHOUT_PHYLOGENY : CODE_END_WITH_PHYLOGENY, writeFn);
	binaryRepresentation_writeString(end_getName(end), writeFn);

	if(endInstance == NULL) {
		iterator = end_getInstanceIterator(end);
		while((endInstance = end_getNext(iterator)) != NULL) {
			assert(endInstance_getParent(endInstance) == NULL);
			endInstance_writeBinaryRepresentation(endInstance, writeFn);
		}
		end_destructInstanceIterator(iterator);
	}
	else {
		end_writeBinaryRepresentationP(endInstance, writeFn);
	}
}

End *end_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	End *end;

	end = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_END_WITHOUT_PHYLOGENY) {
		binaryRepresentation_popNextElementType(binaryString);
		end = end_construct(binaryRepresentation_getStringStatic(binaryString), net);
		while(endInstance_loadFromBinaryRepresentation(binaryString, end) != NULL);
	}
	else {
		if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_END_WITH_PHYLOGENY) {
			binaryRepresentation_popNextElementType(binaryString);
			end = end_construct(binaryRepresentation_getStringStatic(binaryString), net);
			end_setRootInstance(end, endInstance_loadFromBinaryRepresentation(binaryString, end));
			while(endInstance_loadFromBinaryRepresentation(binaryString, end) != NULL);
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

void atomInstance_writeBinaryRepresentation(AtomInstance *atomInstance, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_ATOM_INSTANCE, writeFn);
	binaryRepresentation_writeString(atomInstance_getInstanceName(atomInstance), writeFn);
}

AtomInstance *atomInstance_loadFromBinaryRepresentation(void **binaryString, Atom *atom) {
	const char *name;
	AtomInstance *atomInstance;

	atomInstance = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_ATOM_INSTANCE) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getStringStatic(binaryString);
		atomInstance = atomInstance_construct(atom, end_getInstance(atom_get5End(atom), name),
		end_getInstance(atom_get3End(atom), name));
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

void atom_writeBinaryRepresentation(Atom *atom, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
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

Atom *atom_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	Atom *atom;
	const char *name;
	int32_t length;

	atom = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_ATOM) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getStringStatic(binaryString);
		length = binaryRepresentation_getInteger(binaryString);
		atom = atom_construct(name, length, net);
		while(atomInstance_loadFromBinaryRepresentation(binaryString, atom) != NULL);
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

void adjacencyComponent_writeBinaryRepresentation(AdjacencyComponent *adjacencyComponent, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	Chain *chain;
	End *end;
	AdjacencyComponent_EndIterator *iterator;

	binaryRepresentation_writeElementType(CODE_ADJACENCY_COMPONENT, writeFn);
	binaryRepresentation_writeString(adjacencyComponent_getNestedNetName(adjacencyComponent), writeFn);
	chain = adjacencyComponent_getChain(adjacencyComponent);
	iterator = adjacencyComponent_getEndIterator(adjacencyComponent);
	while((end = adjacencyComponent_getNextEnd(iterator)) != NULL) {
		binaryRepresentation_writeElementType(CODE_ADJACENCY_COMPONENT_END, writeFn);
		binaryRepresentation_writeString(end_getName(end), writeFn);
	}
	adjacencyComponent_destructEndIterator(iterator);
}

AdjacencyComponent *adjacencyComponent_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	AdjacencyComponent *adjacencyComponent;

	adjacencyComponent = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_ADJACENCY_COMPONENT) {
		binaryRepresentation_popNextElementType(binaryString);
		adjacencyComponent = adjacencyComponent_construct2(net, binaryRepresentation_getStringStatic(binaryString));
		while(binaryRepresentation_peekNextElementType(*binaryString) == CODE_ADJACENCY_COMPONENT_END) {
			binaryRepresentation_popNextElementType(binaryString);
			adjacencyComponent_addEnd(adjacencyComponent, net_getEnd(net, binaryRepresentation_getStringStatic(binaryString)));
		}
	}
	return adjacencyComponent;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic link functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void link_writeBinaryRepresentation(Link *link, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_LINK, writeFn);
	binaryRepresentation_writeString(adjacencyComponent_getNestedNetName(link_getAdjacencyComponent(link)), writeFn);
	binaryRepresentation_writeString(end_getName(link_getLeft(link)), writeFn);
	binaryRepresentation_writeString(end_getName(link_getRight(link)), writeFn);
}

Link *link_loadFromBinaryRepresentation(void **binaryString, Chain *chain) {
	AdjacencyComponent *adjacencyComponent;
	End *leftEnd;
	End *rightEnd;
	Link *link;

	link = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_LINK) {
		binaryRepresentation_popNextElementType(binaryString);
		adjacencyComponent = net_getAdjacencyComponent(chain_getNet(chain),
				binaryRepresentation_getStringStatic(binaryString));
		leftEnd = net_getEnd(chain_getNet(chain),
				binaryRepresentation_getStringStatic(binaryString));
		rightEnd = net_getEnd(chain_getNet(chain),
						binaryRepresentation_getStringStatic(binaryString));
		link = link_construct(leftEnd, rightEnd, adjacencyComponent, chain);
	}
	return link;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic chain functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void chain_writeBinaryRepresentation(Chain *chain, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	Link *link;
	binaryRepresentation_writeElementType(CODE_CHAIN, writeFn);
	binaryRepresentation_writeString(chain_getName(chain), writeFn);
	link = chain_getLink(chain, 0);
	while(link != NULL) {
		link_writeBinaryRepresentation(link, writeFn);
		link = link_getNextLink(link);
	}
}

Chain *chain_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	Chain *chain;

	chain = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_CHAIN) {
		binaryRepresentation_popNextElementType(binaryString);
		chain = chain_construct(net, binaryRepresentation_getStringStatic(binaryString));
		while(link_loadFromBinaryRepresentation(binaryString, chain) != NULL);
	}
	return chain;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic operation functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void operation_writeBinaryRepresentation(Operation *operation, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_OPERATION, writeFn);
	binaryRepresentation_writeString(operation_getName(operation), writeFn);
}

Operation *operation_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	Operation *operation;

	operation = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_OPERATION) {
		binaryRepresentation_popNextElementType(binaryString);
		operation = operation_construct(net, binaryRepresentation_getStringStatic(binaryString));
	}
	return operation;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void net_writeBinaryRepresentation(Net *net, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	Net_SequenceIterator *sequenceIterator;
	Net_EndIterator *endIterator;
	Net_AtomIterator *atomIterator;
	Net_AdjacencyComponentIterator *adjacencyComponentIterator;
	Net_ChainIterator *chainIterator;
	Net_OperationIterator *operationIterator;
	Sequence *sequence;
	End *end;
	Atom *atom;
	AdjacencyComponent *adjacencyComponent;
	Chain *chain;
	Operation *operation;

	binaryRepresentation_writeString(net_getName(net), writeFn);

	eventTree_writeBinaryRepresentation(net_getEventTree(net), writeFn);

	sequenceIterator = net_getSequenceIterator(net);
	while((sequence = net_getNextSequence(sequenceIterator)) != NULL) {
		sequence_writeBinaryRepresentation(sequence, writeFn);
	}
	net_destructSequenceIterator(sequenceIterator);

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

	operationIterator = net_getOperationIterator(net);
	while((operation = net_getNextOperation(operationIterator)) != NULL) {
		operation_writeBinaryRepresentation(operation, writeFn);
	}
	net_destructOperationIterator(operationIterator);

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
}

Net *net_loadFromBinaryRepresentation(void **binaryString, NetDisk *netDisk) {
	Net *net;
	net = net_construct(binaryRepresentation_getStringStatic(binaryString), netDisk);
	assert(eventTree_loadFromBinaryRepresentation(binaryString, net) != NULL);
	while(sequence_loadFromBinaryRepresentation(binaryString, net) != NULL);
	while(end_loadFromBinaryRepresentation(binaryString, net) != NULL);
	while(atom_loadFromBinaryRepresentation(binaryString, net) != NULL);
	while(operation_loadFromBinaryRepresentation(binaryString, net) != NULL);
	while(adjacencyComponent_loadFromBinaryRepresentation(binaryString, net) != NULL);
	while(chain_loadFromBinaryRepresentation(binaryString, net) != NULL);
	return net;
}

