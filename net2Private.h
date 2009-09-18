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

struct _sequence {
	char *name;
	int32_t length;
	char *file;
	NetDisk *netDisk;
};

struct _endInstance {
	char *instance;
	End *end;
	int32_t coordinate;
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
