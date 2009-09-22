#ifndef NET_PRIVATE_H_
#define NET_PRIVATE_H_

#include "net2.h"
#include "commonC.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic data structure declarations
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct _event {
	char *name;
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

struct _sequence {
	char *name;
	int32_t length;
	char *eventName;
	char *file;
	NetDisk *netDisk;
};

struct _endInstance {
	char *instance;
	End *end;
	int32_t coordinate;
	Event *event;
	Sequence *sequence;
	EndInstance *adjacency;
	EndInstance *adjacency2;
	Operation *operation;
	AtomInstance *atomInstance;
	EndInstance *parent;
	struct List *children;
};

struct _end {
	char *name;
	struct avl_table *endInstances;
	Atom *attachedAtom;
	AdjacencyComponent *adjacencyComponent;
	Net *net;
};

struct _atomInstance {
	EndInstance *leftEndInstance;
	AtomInstance *rInstance;
	Atom *atom;
};

struct AtomContents {
	char *name;
	struct avl_table *atomInstances;
	int32_t length;
	Net *net;
};

struct _atom {
	struct AtomContents *atomContents;
	End *leftEnd;
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
	Net *net;
	Link *link;
	int32_t linkNumber;
	int32_t chainIndex;
};

struct _operation {
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
	char *netsDatabaseName;
	char *sequencesDatabaseName;
	TCBDB *netsDatabase;
	TCBDB *sequencesDatabase;
	struct avl_table *sequences;
	struct avl_table *nets;
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
char *database_getRecord(TCBDB *database, const char *key);

/*
 * Writes a key value record to the database.
 */
int32_t database_writeRecord(TCBDB *database, const char *key, const char *value);

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

#define CODE_EVENT 0
#define CODE_POINTER 1
#define CODE_EVENT_TREE 1
#define CODE_SEQUENCE 0
#define CODE_ADJACENCY 3
#define CODE_PARENT 3
#define CODE_END_INSTANCE 1
#define CODE_END_INSTANCE_WITH_EVENT 1
#define CODE_END_INSTANCE_WITH_COORDINATES 2
#define CODE_END 2
#define CODE_ATOM_INSTANCE 3
#define CODE_ATOM 3
#define CODE_ADJACENCY_COMPONENT 2
#define CODE_ADJACENCY_COMPONENT_CHAIN 2
#define CODE_ADJACENCY_COMPONENT_END 2
#define CODE_LINK 2
#define CODE_CHAIN 3
#define CODE_OPERATION 3

/*
 * Writes a code for the element type.
 */
void binaryRepresentation_writeElementType(int32_t elementCode, void (*writeFn)(const char *, ...));

/*
 * Writes a string, containing no white space, to the binary stream.
 */
void binaryRepresentation_writeString(const char *string, void (*writeFn)(const char *, ...));

/*
 * Writes an integer to the binary stream
 */
void binaryRepresentation_writeInteger(int32_t i, void (*writeFn)(const char *, ...));

/*
 * Writes a float to the binary stream.
 */
void binaryRepresentation_writeFloat(float f, void (*writeFn)(const char *, ...));

/*
 * Returns indicating which element is next, but does not increment the string pointer.
 */
int32_t binaryRepresentation_peekNextElementType(char **binaryString);

/*
 * Returns indicating which element is next, while incrementing the string pointer.
 */
int32_t binaryRepresentation_popNextElementType(char **binaryString);

/*
 * Parses out a string, returning it in a newly allocated string which must be freed.
 */
char *binaryRepresentation_getString(char **binaryString);

/*
 * Parses out a string, placing the memory in a buffer owned by the function. Thid buffer
 * will be overidden by the next call to the function.
 */
const char *binaryRepresentation_getStringStatic(char **binaryString);

/*
 * Parses an integer from binary string.
 */
int32_t binaryRepresentation_getInteger(char **binaryString);

/*
 * Parses a float from the binary string.
 */
float binaryRepresentation_getFloat(char **binaryString);

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
void event_writeBinaryRepresentation(Event *event, void (*writeFn)(const char *string, ...));

/*
 * Loads a event into memory from a binary representation of the event.
 */
Event *event_loadFromBinaryRepresentation(char **binaryString, EventTree *eventTree);

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
void eventTree_writeBinaryRepresentation(EventTree *eventTree, void (*writeFn)(const char *string, ...));

/*
 * Loads a eventTree into memory from a binary representation of the eventTree.
 */
EventTree *eventTree_loadFromBinaryRepresentation(char **binaryString, Net *net);

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
void sequence_writeBinaryRepresentation(Sequence *sequence, void (*writeFn)(const char *string, ...));

/*
 * Creates a binary representation of the sequence, returned as a char string.
 */
char *sequence_makeBinaryRepresentation(Sequence *sequence);

/*
 * Loads a sequence into memory from a binary representation of the sequence.
 */
Sequence *sequence_loadFromBinaryRepresentation(char **binaryString, NetDisk *netDisk);

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
void endInstance_writeBinaryRepresentation(EndInstance *endInstance, void (*writeFn)(const char *string, ...));

/*
 * Loads a net into memory from a binary representation of the net.
 */
EndInstance *endInstance_loadFromBinaryRepresentation(char **binaryString, End *end);

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
void end_writeBinaryRepresentation(End *end, void (*writeFn)(const char *string, ...));

/*
 * Loads a net into memory from a binary representation of the net.
 */
End *end_loadFromBinaryRepresentation(char **binaryString, Net *net);

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
void atomInstance_writeBinaryRepresentation(AtomInstance *atomInstance, void (*writeFn)(const char *string, ...));

/*
 * Loads a net into memory from a binary representation of the net.
 */
AtomInstance *atomInstance_loadFromBinaryRepresentation(char **binaryString, Atom *atom);

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
void atom_writeBinaryRepresentation(Atom *atom, void (*writeFn)(const char *string, ...));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Atom *atom_loadFromBinaryRepresentation(char **binaryString, Net *net);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Adjacency component functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs an adjacency component.
 */
void adjacencyComponent_destruct(AdjacencyComponent *adjacencyComponent);

/*
 * Sets the chain the adjacency component is part of.
 */
void adjacencyComponent_setChain(AdjacencyComponent *adjacencyComponent, Chain *chain);

/*
 * Write a binary representation of the adjacencyComponent to the write function.
 */
void adjacencyComponent_writeBinaryRepresentation(AdjacencyComponent *adjacencyComponent, void (*writeFn)(const char *string, ...));

/*
 * Loads a net into memory from a binary representation of the net.
 */
AdjacencyComponent *adjacencyComponent_loadFromBinaryRepresentation(char **binaryString, Net *net);

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
void link_writeBinaryRepresentation(Link *link, void (*writeFn)(const char *string, ...));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Link *link_loadFromBinaryRepresentation(char **binaryString, Chain *chain);

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
 * Sets the chain's index.
 */
void chain_setIndex(Chain *chain, int32_t index);

/*
 * Write a binary representation of the chain to the write function.
 */
void chain_writeBinaryRepresentation(Chain *chain, void (*writeFn)(const char *string, ...));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Chain *chain_loadFromBinaryRepresentation(char **binaryString, Net *net);

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
void *operation_writeBinaryRepresentation(Operation *operation, void (*writeFn)(const char *string, ...));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Operation *operation_loadFromBinaryRepresentation(char **binaryString, Net *net);

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
void net_writeBinaryRepresentation(Net *net, void (*writeFn)(const char *string, ...));

/*
 * Creates a binary representation of the net, returned as a char string.
 */
char *net_makeBinaryRepresentation(Net *net);

/*
 * Loads a net into memory from a binary representation of the net.
 */
Net *net_loadFromBinaryRepresentation(char **binaryString, NetDisk *netDisk);


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Net disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Removes the sequence from the disk.
 */
int32_t netDisk_deleteSequenceFromDisk(NetDisk *netDisk, const char *sequenceName);

/*
 * Removes the net from the disk.
 */
int32_t netDisk_deleteNetFromDisk(NetDisk *netDisk, const char *netName);

/*
 * Adds a newly constructed sequence to the memory of the netDisk.
 */
void netDisk_addSequence(NetDisk *netDisk, Sequence *sequence);

/*
 * Registers the sequence is being freed from memory.
 */
void netDisk_unloadSequence(NetDisk *netDisk, Sequence *sequence);

/*
 * Adds a newly constructed net to the memory of the netDisk.
 */
void netDisk_addNet(NetDisk *netDisk, Net *net);

/*
 * Registers the net is being freed from memory.
 */
void netDisk_unloadNet(NetDisk *netDisk, Net *net);

#endif
