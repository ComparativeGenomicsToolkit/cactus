#ifndef NET_PRIVATE_H_
#define NET_PRIVATE_H_

#include <inttypes.h>
#include "net2.h"
#include "avl.h"

struct _sequence {
	char *name;
	int32_t length;
	char *file;
	NetDisk *netDisk;
};

struct _endInstance {
	/*
	 * Represents an end instance.
	 */
	char *instance;
	int32_t coordinate;
	Sequence *sequence;
	EndInstance *adjacency;
	EndInstance *adjacency2;
	Operation *operation;
	AtomInstance *atomInstance;
	EndInstance *parent;
	struct List *children;
	End *end;
};

struct _end {
	/*
	 * Represents a cap/stub.
	 */
	char *name;
	struct avl_table *endInstances;
	Atom *attachedAtom;
	AdjacencyComponent *adjacencyComponent;
	Net *net;
};

struct AtomContents {
	/*
	 * Holds a single copy of the internals of each atom (to avoid duplicating in reverse).
	 */
	char *name;
	struct avl_table *atomInstances;
	int32_t length;
	Net *net;
};

struct _atom {
	/*
	 * Represents an atom.
	 */
	struct AtomContents *atomContents;
	End *leftEnd;
	Atom *rAtom;
};

struct _atomInstance {
	/*
	 * Represents an instance of an atom.
	 */
	EndInstance *leftEndInstance;
	AtomInstance *rInstance;
	Atom *atom;
};

struct _adjacencyComponent {
	Net *net;
	Chain *chain;
	char *nestedNetName;
	struct avl_table *ends;
};


struct _chain {
	/*
	 * A doubly linked list which can hold sub chains (to form the the net structure).
	 */
	End *leftEnd;
	End *rightEnd;
	AdjacencyComponent *adjacencyComponent;
	//previous link in the chain.
	Chain *pLink;
	//next link in the chain.
	Chain *nLink;
	//Pointer to contained subproblem
	Net *net;
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
	int32_t chainIndex;
};

struct _netDisk {
	char *fileName;
	Net *net;
};

///
//Private functions.
///

/*
 * Adds the instance to the atom.
 */

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
 * Adds in the instance to the atom.
 */
void atom_addInstance(Atom *atom, AtomInstance *atomInstance);

/*
 * Removes the instance from the atom.
 */
void atom_removeInstance(Atom *atom, AtomInstance *atomInstance);

/*
 * Sets the chain the adjacency component is part of.
 */
void adjacencyComponent_setChain(AdjacencyComponent *adjacencyComponent, Chain *chain);

/*
 * Sets the chain's index.
 */
void chain_setIndex(Chain *chain, int32_t index);

/*
 * Sets the operation's index.
 */
void operation_setIndex(Operation *operation, int32_t index);

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
 * Adds the operation to the net.
 */
void net_addOperation(Net *net, Operation *operation);

/*
 * Remove the operation from the net.
 */
void net_removeOperation(Net *net, Operation *operation);

/*
 * Adds a newly constructed sequence to the memory of the netDisk.
 */
void netDisk_addSequence(NetDisk *netDisk, Sequence *sequence);

/*
 * Registers the sequence is being freed from memory.
 */
void netDisk_unloadSequence(NetDisk *netDisk, Sequence *sequence);

/*
 * Adds the chain to the net.
 */
void net_addChain(Net *net, Chain *chain);

/*
 * Remove the chain from the net.
 */
void net_removeChain(Net *net, Chain *chain);

/*
 * Adds a newly constructed net to the memory of the netDisk.
 */
void netDisk_addNet(NetDisk *netDisk, Net *net);

/*
 * Registers the net is being freed from memory.
 */
void netDisk_unloadNet(NetDisk *netDisk, Net *net);

#endif
