#ifndef NET_PRIVATE_H_
#define NET_PRIVATE_H_

#include "net.h"
#include "commonC.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic data structure declarations
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct _metaEvent {
	char *name;
	char *header;
	NetDisk *netDisk;
	int32_t referenceCount;
};

struct _event {
	MetaEvent *metaEvent;
	struct List *children;
	float branchLength;
	Event *parent;
	EventTree *eventTree;
};

struct _eventTree {
	Event *rootEvent;
	struct avl_table *events;
	Net *net;
};

struct _metaSequence {
	char *name;
	int64_t fileOffset;
	int32_t start;
	int32_t length;
	char *eventName;
	NetDisk *netDisk;
	char *header;
	int32_t referenceCount;
};

struct _sequence {
	MetaSequence *metaSequence;
	Net *net;
};

typedef struct _endInstanceContents {
	char *instance;
	int32_t coordinate;
	int32_t strand;
	int32_t side;
	Event *event;
	Sequence *sequence;
	EndInstance *adjacency;
	EndInstance *adjacency2;
	Operation *operation;
	AtomInstance *atomInstance;
	EndInstance *parent;
	struct List *children;
} EndInstanceContents;

struct _endInstance {
	EndInstanceContents *endInstanceContents;
	End *end;
	EndInstance *rEndInstance;
};

typedef struct _endContents {
	EndInstance *rootInstance;
	char *name;
	Atom *attachedAtom;
	struct avl_table *endInstances;
	AdjacencyComponent *adjacencyComponent;
	Net *net;
} EndContents;

struct _end_instanceIterator {
	struct avl_traverser *iterator;
	End *end;
};

struct _end {
	EndContents *endContents;
	int32_t orientation;
	End *rEnd;
};

struct _atomInstance {
	EndInstance *_5EndInstance;
	AtomInstance *rInstance;
	Atom *atom;
};

typedef struct _atomContents {
	char *name;
	struct avl_table *atomInstances;
	int32_t length;
	Net *net;
} AtomContents;

struct _atom_instanceIterator {
	struct avl_traverser *iterator;
	Atom *atom;
};

struct _atom {
	AtomContents *atomContents;
	int32_t orientation;
	End *_5End;
	Atom *rAtom;
};

struct _adjacencyComponent {
	Net *net;
	Chain *chain;
	char *nestedNetName;
	struct avl_table *ends;
};

struct _link {
	End *leftEnd;
	End *rightEnd;
	Chain *chain;
	AdjacencyComponent *adjacencyComponent;
	//previous link in the chain.
	Link *pLink;
	//next link in the chain.
	Link *nLink;
	//index of the link in the chain.
	int32_t linkIndex;
};

struct _chain {
	char *name;
	Net *net;
	Link *link;
	int32_t linkNumber;
	int32_t chainIndex;
};

struct _operation {
	char *name;
	Net *net;
	int32_t index;
};

struct _net {
	char *name;
	EventTree *eventTree;
	struct avl_table *sequences;
	struct avl_table *ends;
	struct avl_table *atoms;
	struct avl_table *adjacencyComponents;
	struct avl_table *chains;
	struct avl_table *operations;
	char *parentNetName;
	NetDisk *netDisk;
	int32_t operationIndex;
	int32_t chainIndex;
};

struct _netDisk {
	char *stringFile;
	int64_t stringFileLength;
	char *netsDatabaseName;
	char *metaDataDatabaseName;
	char *iDDatabaseName;
	TCBDB *netsDatabase;
	TCBDB *metaDataDatabase;
	TCBDB *iDDatabase;
	struct avl_table *metaSequences;
	struct avl_table *metaEvents;
	struct avl_table *nets;
	int64_t uniqueNumber;
	int64_t maxUniqueNumber;
};

///
//Private functions.
///

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions on a sorted set and its iterator
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a sorted set, using the given comparison function.
 */
struct avl_table *sortedSet_construct(int32_t (*compareFn)(const void *, const void *, void *));

/*
 * Destructs the sorted set, applying the destruct function to each element.
 */
void sortedSet_destruct(struct avl_table *sortedSet, void (*destructElementFn)(void *, void *));

/*
 * Inserts the object into the sorted set.
 */
void sortedSet_insert(struct avl_table *sortedSet, void *object);

/*
 * Finds the objects in the sorted set, or returns null.
 */
void *sortedSet_find(struct avl_table *sortedSet, void *object);

/*
 * Deletes the object in the sorted set.
 */
void sortedSet_delete(struct avl_table *sortedSet, void *object);

/*
 * Gets the number of elements in the sorted set.
 */
int32_t sortedSet_getLength(struct avl_table *sortedSet);

/*
 * Gets the first element (with lowest value), in the sorted set.
 */
void *sortedSet_getFirst(struct avl_table *items);

/*
 * Constructs an iterator for the sorted set.
 */
struct avl_traverser *iterator_construct(struct avl_table *items);

/*
 * Destructs an iterator for the sorted set.
 */
void iterator_destruct(struct avl_traverser *iterator);

/*
 * Gets next element in the sorted set.
 */
void *iterator_getNext(struct avl_traverser *iterator);

/*
 * Gets the previous element in the sorted set.
 */
void *iterator_getPrevious(struct avl_traverser *iterator);

/*
 * Copies the iterator.
 */
struct avl_traverser *iterator_copy(struct avl_traverser *iterator);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Database functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a sorted-set database object.
 */
TCBDB *database_construct(const char *name);

/*
 * Destructs a database.
 */
void database_destruct(TCBDB *database);

/*
 * Returns number of records in database.
 */
int32_t database_getNumberOfRecords(TCBDB *database);

/*
 * Gets a record from the database, given the key. The record is in newly allocated memory, and must be freed.
 */
void *database_getRecord(TCBDB *database, const char *key);

/*
 * Writes a key value record to the database.
 */
int32_t database_writeRecord(TCBDB *database, const char *key, const void *value, int32_t sizeOfRecord);

/*
 * Removes a record from the database.
 */
int32_t database_removeRecord(TCBDB *database, const char *key);

/*
 * Constructs an iterator over the sorted database records.
 */
BDBCUR *databaseIterator_construct(TCBDB *database);

/*
 * Gets the next element from the database iterator.
 */
const char *databaseIterator_getNext(BDBCUR *iterator);

/*
 * Destructs a database iterator.
 */
void databaseIterator_destruct(BDBCUR *iterator);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions for serialising the objects.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Codes used to define which objects are being encoded/decoded from a binary database stream.
 */

#define CODE_META_EVENT 19
#define CODE_EVENT 0
#define CODE_POINTER 1
#define CODE_EVENT_TREE 2
#define CODE_META_SEQUENCE 3
#define CODE_SEQUENCE 4
#define CODE_ADJACENCY 5
#define CODE_PARENT 6
#define CODE_END_INSTANCE 7
#define CODE_END_INSTANCE_WITH_EVENT 8
#define CODE_END_INSTANCE_WITH_COORDINATES 9
#define CODE_END_WITHOUT_PHYLOGENY 10
#define CODE_END_WITH_PHYLOGENY 11
#define CODE_ATOM_INSTANCE 12
#define CODE_ATOM 13
#define CODE_ADJACENCY_COMPONENT 14
#define CODE_ADJACENCY_COMPONENT_END 15
#define CODE_LINK 16
#define CODE_CHAIN 17
#define CODE_OPERATION 18

/*
 * Writes a code for the element type.
 */
void binaryRepresentation_writeElementType(char elementCode, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Writes a string to the binary stream.
 */
void binaryRepresentation_writeString(const char *string, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Writes an integer to the binary stream
 */
void binaryRepresentation_writeInteger(int32_t i, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Writes an integer to the binary stream
 */
void binaryRepresentation_write64BitInteger(int64_t i, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Writes a float to the binary stream.
 */
void binaryRepresentation_writeFloat(float f, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Returns indicating which element is next, but does not increment the string pointer.
 */
char binaryRepresentation_peekNextElementType(void *binaryString);

/*
 * Returns indicating which element is next, while incrementing the string pointer.
 */
char binaryRepresentation_popNextElementType(void **binaryString);

/*
 * Returns a char from the binary string, moving up the pointer.
 */
char binaryRepresentation_getChar(void **binaryString);

/*
 * Parses out a string, returning it in a newly allocated string which must be freed.
 */
char *binaryRepresentation_getString(void **binaryString);

/*
 * Parses out a string, placing the memory in a buffer owned by the function. Thid buffer
 * will be overidden by the next call to the function.
 */
const char *binaryRepresentation_getStringStatic(void **binaryString);

/*
 * Parses an integer from binary string.
 */
int32_t binaryRepresentation_getInteger(void **binaryString);

/*
 * Parses an 64bit integer from binary string.
 */
int64_t binaryRepresentation_get64BitInteger(void **binaryString);

/*
 * Parses a float from the binary string.
 */
float binaryRepresentation_getFloat(void **binaryString);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs a meta event.
 */
void metaEvent_destruct(MetaEvent *metaEvent);

/*
 * Increases the number of references (held by event objects), to one.
 */
void metaEvent_increaseReferenceCount(MetaEvent *metaEvent);

/*
 * Descrease the number of references, by one. If it gets to zero then the object will
 * be destroyed.
 */
void metaEvent_decreaseReferenceCountAndDestructIfZero(MetaEvent *metaEvent);

/*
 * Creates a binary representation of the eventTree, returned as a char string.
 */
void metaEvent_writeBinaryRepresentation(MetaEvent *metaEvent, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a eventTree into memory from a binary representation of the eventTree.
 */
MetaEvent *metaEvent_loadFromBinaryRepresentation(void **binaryString, NetDisk *netDisk);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event functions.
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
void event_writeBinaryRepresentation(Event *event, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a event into memory from a binary representation of the event.
 */
Event *event_loadFromBinaryRepresentation(void **binaryString, EventTree *eventTree);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event tree functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs the event tree and all its events.
 */
void eventTree_destruct(EventTree *eventTree);

/*
 * Adds the end instance to the event tree.
 */
void eventTree_addEvent(EventTree *eventTree, Event *event);

/*
 * Removes the instance from the event tree.
 */
void eventTree_removeEvent(EventTree *eventTree, Event *event);

/*
 * Creates a binary representation of the eventTree, returned as a char string.
 */
void eventTree_writeBinaryRepresentation(EventTree *eventTree, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a eventTree into memory from a binary representation of the eventTree.
 */
EventTree *eventTree_loadFromBinaryRepresentation(void **binaryString, Net *net);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a meta sequence using an existing reference to a sequence in the sequence file.
 */
MetaSequence *metaSequence_construct2(const char *name, int32_t start, int32_t length, int64_t fileOffset, const char *header,
		const char *eventName, NetDisk *netDisk);

/*
 * Destructs a meta sequence.
 */
void metaSequence_destruct(MetaSequence *metaSequence);

/*
 * Gets the file offset location of the string backing the metasequence.
 */
int64_t metaSequence_getFileOffset(MetaSequence *metaSequence);

/*
 * Increases the number of references (held by sequence objects), to one.
 */
void metaSequence_increaseReferenceCount(MetaSequence *metaSequence);

/*
 * Descrease the number of references, by one. If it gets to zero then the object will
 * be destroyed.
 */
void metaSequence_decreaseReferenceCountAndDestructIfZero(MetaSequence *metaSequence);

/*
 * Creates a binary representation of the eventTree, returned as a char string.
 */
void metaSequence_writeBinaryRepresentation(MetaSequence *metaSequence, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a eventTree into memory from a binary representation of the eventTree.
 */
MetaSequence *metaSequence_loadFromBinaryRepresentation(void **binaryString, NetDisk *netDisk);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Write a binary representation of the sequence to the write function.
 */
void sequence_writeBinaryRepresentation(Sequence *sequence, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a sequence into memory from a binary representation of the sequence.
 */
Sequence *sequence_loadFromBinaryRepresentation(void **binaryString, Net *net);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//End instance functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs the end instance, but not any connecting objects.
 */
void endInstance_destruct(EndInstance *endInstance);

/*
 * Gets the atom instance associated with the end, or NULL, if the end has no associated atom end at this level.
 * The atom instance returned will have the end instance on its left side.
 */
void endInstance_setAtomInstance(EndInstance *endInstance, AtomInstance *atomInstance);

/*
 * Sets the operation associated with the end instance.
 */
void endInstance_setOperation(EndInstance *endInstance, Operation *operation);

/*
 * Sets any adjacent ends for the first adjacency of the instance to NULL;
 */
void endInstance_breakAdjacency1(EndInstance *endInstance);

/*
 * Sets any adjacent ends for the alternative adjacency of the instance to NULL;
 */
void endInstance_breakAdjacency2(EndInstance *endInstance);

/*
 * Write a binary representation of the endInstance to the write function.
 */
void endInstance_writeBinaryRepresentation(EndInstance *endInstance, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
EndInstance *endInstance_loadFromBinaryRepresentation(void **binaryString, End *end);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//End functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs the end and any contained end instances.
 */
void end_destruct(End *end);

/*
 * Adds the end instance to the end.
 */
void end_addInstance(End *end, EndInstance *endInstance);

/*
 * Removes the instance from the end.
 */
void end_removeInstance(End *end, EndInstance *endInstance);

/*
 * Sets the adjacency component that the end is part of.
 */
void end_setAdjacencyComponent(End *end, AdjacencyComponent *adjacencyComponent);

/*
 * Write a binary representation of the end to the write function.
 */
void end_writeBinaryRepresentation(End *end, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
End *end_loadFromBinaryRepresentation(void **binaryString, Net *net);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Atom instance functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destruct the atom instance, does not destruct ends.
 */
void atomInstance_destruct(AtomInstance *atomInstance);

/*
 * Write a binary representation of the atomInstance to the write function.
 */
void atomInstance_writeBinaryRepresentation(AtomInstance *atomInstance, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
AtomInstance *atomInstance_loadFromBinaryRepresentation(void **binaryString, Atom *atom);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Atom functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs the atom and all atom instances it contains.
 */
void atom_destruct(Atom *atom);

/*
 * Adds in the instance to the atom.
 */
void atom_addInstance(Atom *atom, AtomInstance *atomInstance);

/*
 * Removes the instance from the atom.
 */
void atom_removeInstance(Atom *atom, AtomInstance *atomInstance);

/*
 * Write a binary representation of the atom to the write function.
 */
void atom_writeBinaryRepresentation(Atom *atom, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Atom *atom_loadFromBinaryRepresentation(void **binaryString, Net *net);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Adjacency component functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a nested net without having the nested net loaded in memory.
 */
AdjacencyComponent *adjacencyComponent_construct2(Net *net, const char *nestedNetName);

/*
 * Destructs an adjacency component.
 */
void adjacencyComponent_destruct(AdjacencyComponent *adjacencyComponent);

/*
 * Adds an end to the adjacency component. This is private, because it's normally calculated
 * by a call to the adjacencyComponent_updateContainedEnds function.
 */
void adjacencyComponent_addEnd(AdjacencyComponent *adjacencyComponent, End *end);

/*
 * Sets the chain the adjacency component is part of.
 */
void adjacencyComponent_setChain(AdjacencyComponent *adjacencyComponent, Chain *chain);

/*
 * Write a binary representation of the adjacencyComponent to the write function.
 */
void adjacencyComponent_writeBinaryRepresentation(AdjacencyComponent *adjacencyComponent, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
AdjacencyComponent *adjacencyComponent_loadFromBinaryRepresentation(void **binaryString, Net *net);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Link functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs the link and all subsequent nLinks.
 */
void link_destruct(Link *link);

/*
 * Write a binary representation of the link to the write function.
 */
void link_writeBinaryRepresentation(Link *link, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Link *link_loadFromBinaryRepresentation(void **binaryString, Chain *chain);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Chain functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs the chain.
 */
void chain_destruct(Chain *chain);

/*
 * Add the link to the chain.
 */
void chain_addLink(Chain *chain, Link *childLink);

/*
 * Write a binary representation of the chain to the write function.
 */
void chain_writeBinaryRepresentation(Chain *chain, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Chain *chain_loadFromBinaryRepresentation(void **binaryString, Net *net);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Operation functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs the operation.
 */
void operation_destruct(Operation *operation);

/*
 * Sets the operation's index.
 */
void operation_setIndex(Operation *operation, int32_t index);

/*
 * Write a binary representation of the operation to the write function.
 */
void operation_writeBinaryRepresentation(Operation *operation, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Operation *operation_loadFromBinaryRepresentation(void **binaryString, Net *net);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Net functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Adds the event tree for the net to the net.
 */
void net_addEventTree(Net *net, EventTree *eventTree);

/*
 * Adds the sequence to the net.
 */
void net_addSequence(Net *net, Sequence *sequence);

/*
 * Removes the sequence from the net.
 */
void net_removeSequence(Net *net, Sequence *sequence);

/*
 * Adds the atom to the net.
 */
void net_addAtom(Net *net, Atom *atom);

/*
 * Remove the atom from the net.
 */
void net_removeAtom(Net *net, Atom *atom);

/*
 * Adds the end to the net.
 */
void net_addEnd(Net *net, End *end);

/*
 * Remove the end from the net.
 */
void net_removeEnd(Net *net, End *end);

/*
 * Adds the adjacency component to the net.
 */
void net_addAdjacencyComponent(Net *net, AdjacencyComponent *adjacencyComponent);

/*
 * Removes an empty adjacency component from the net.
 */
void net_removeAdjacencyComponent(Net *net, AdjacencyComponent *adjacencyComponent);

/*
 * Sets the parent adjacency component of the net.
 */
void net_setParentAdjacencyComponent(Net *net, AdjacencyComponent *adjacencyComponent);

/*
 * Adds the chain to the net.
 */
void net_addChain(Net *net, Chain *chain);

/*
 * Remove the chain from the net.
 */
void net_removeChain(Net *net, Chain *chain);

/*
 * Adds the operation to the net.
 */
void net_addOperation(Net *net, Operation *operation);

/*
 * Remove the operation from the net.
 */
void net_removeOperation(Net *net, Operation *operation);

/*
 * Write a binary representation of the net to the write function.
 */
void net_writeBinaryRepresentation(Net *net, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Net *net_loadFromBinaryRepresentation(void **binaryString, NetDisk *netDisk);


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Net disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Adds a newly constructed net to the memory of the netDisk.
 */
void netDisk_addNet(NetDisk *netDisk, Net *net);

/*
 * Removes the net from the disk.
 */
int32_t netDisk_deleteNetFromDisk(NetDisk *netDisk, const char *netName);

/*
 * Registers the net is being freed from memory.
 */
void netDisk_unloadNet(NetDisk *netDisk, Net *net);

/*
 * Functions on meta sequences.
 */

/*
 * Adds a newly constructed meta sequence to the memory of the netDisk.
 */
void netDisk_addMetaSequence(NetDisk *netDisk, MetaSequence *metaSequence);

/*
 * Removes the meta sequence from the disk.
 */
int32_t netDisk_deleteMetaSequenceFromDisk(NetDisk *netDisk, const char *metaSequenceName);

/*
 * Registers the meta sequence is being freed from memory.
 */
void netDisk_unloadMetaSequence(NetDisk *netDisk, MetaSequence *metaSequence);

/*
 * Gets a meta sequence from the set held in memory.
 */
MetaSequence *netDisk_getMetaSequenceInMemory(NetDisk *netDisk, const char *metaSequenceName);

/*
 * Gets the meta sequence for an object.
 */
MetaSequence *netDisk_getMetaSequence(NetDisk *netDisk, const char *metaSequenceName);

/*
 * Functions on meta events.
 */

/*
 * Adds a newly constructed meta sequence to the memory of the netDisk.
 */
void netDisk_addMetaEvent(NetDisk *netDisk, MetaEvent *metaEvent);

/*
 * Deletes the meta event from the disk.
 */
int32_t netDisk_deleteMetaEventFromDisk(NetDisk *netDisk, const char *metaEventName);

/*
 * Registers the meta event is being freed from memory.
 */
void netDisk_unloadMetaEvent(NetDisk *netDisk, MetaEvent *metaEvent);

/*
 * Gets a meta event from the pool of meta events in memory, or returns null.
 */
MetaEvent *netDisk_getMetaEventInMemory(NetDisk *netDisk, const char *metaEventName);

/*
 * Gets a meta event from either disk or memory.
 */
MetaEvent *netDisk_getMetaEvent(NetDisk *netDisk, const char *metaEventName);

/*
 * Functions on strings stored by the net disk.
 */

/*
 * Adds the sequence string to the bucket of sequence.
 *
 * This function is NOT thread safe.
 */
int64_t netDisk_addString(NetDisk *netDisk, const char *string, int32_t length);

/*
 * Retrieves a string from the bucket of sequence.
 */
char *netDisk_getString(NetDisk *netDisk, int64_t offset, int32_t start, int32_t length, int32_t strand);

/*
 * Retrieves the next unique ID.
 */
//int64_t netDisk_getUniqueID(NetDisk *netDisk);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Useful utility functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Makes a binary representation of an object, using a passed function which writes
 * out the representation of the considered object.
 */
void *netMisc_makeBinaryRepresentation(void *object, void (*writeBinaryRepresentation)(void *, void (*writeFn)(const void * ptr, size_t size, size_t count)), int32_t *recordSize);

#endif
