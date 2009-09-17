#ifndef NET_H_
#define NET_H_

#include <inttypes.h>

/*
 * Includes for Tokyo Cabinet.
 */
#include <tcutil.h>
#include <tcbdb.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "avl.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic data structure declarations (contents hidden)
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

typedef struct _sequence Sequence;
typedef struct _end End;
typedef struct _endInstance EndInstance;
typedef struct _atomInstance AtomInstance;
typedef struct _atom Atom;
typedef struct _adjacencyComponent AdjacencyComponent;
typedef struct _link Link;
typedef struct _chain Chain;
typedef struct _operation Operation;
typedef struct _net Net;
typedef struct _netDisk NetDisk;

typedef struct avl_traverser End_InstanceIterator;
typedef struct avl_traverser Atom_InstanceIterator;
typedef struct avl_traverser AdjacencyComponent_EndIterator;
typedef struct avl_traverser Net_SequenceIterator;
typedef struct avl_traverser Net_EndIterator;
typedef struct avl_traverser Net_AtomIterator;
typedef struct avl_traverser Net_AdjacencyComponentIterator;
typedef struct avl_traverser Net_ChainIterator;
typedef struct avl_traverser Net_OperationIterator;
typedef BDBCUR NetDisk_SequenceNameIterator;
typedef struct avl_traverser NetDisk_SequenceIterator;
typedef BDBCUR NetDisk_NetNameIterator;
typedef struct avl_traverser NetDisk_NetIterator;

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a sequence. If the sequence name is already in netDisk, then an assert error is created.
 */
Sequence *sequence_construct(const char *name, int32_t length, const char *file, NetDisk *netDisk);

/*
 * Destructs the sequence.
 */
void sequence_destruct(Sequence *sequence);

/*
 * Gets the length of the sequence.
 */
int32_t sequence_getLength(Sequence *sequence);

/*
 * Gets the name of the sequence.
 */
const char *sequence_getName(Sequence *sequence);

/*
 * Gets the file containing the sequence.
 */
const char *sequence_getFile(Sequence *sequence);

/*
 * Creates a binary representation of the sequence, returned as a char string.
 */
char *sequence_makeBinaryRepresentation(Sequence *sequence);

/*
 * Creates an XML representation of the sequence, returned as a char string.
 */
char *sequence_makeXMLRepresentation(Sequence *sequence);

/*
 * Loads a sequence into memory from a binary representation of the sequence.
 */
Sequence *sequence_loadFromBinaryRepresentation(char *binaryString, NetDisk *netDisk);

/*
 * Loads a sequence into memory from an XML representation of the sequence.
 */
Sequence *sequence_loadFromXMLRepresentation(char *xmlString, NetDisk *netDisk);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic end instance functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an end instance, but not its connecting objects. Instance is the suffix m of the instance name n.m.
 */
EndInstance *endInstance_construct(const char *instance, End *end);

/*
 * As default constructor, but also sets the instance's coordinates.
 */
EndInstance *endInstance_constructWithCoordinates(const char *instance, End *end, int32_t coordinate, Sequence *sequence);

/*
 * Returns the m part of an instance's n.m name.
 */
const char *endInstance_getInstanceName(EndInstance *endInstance);

/*
 * Returns the n part of an instance's n.m name.
 */
const char *endInstance_getElementName(EndInstance *endInstance);

/*
 * Gets the complete name of an instance. This involves a new memory allocation, you are therefore responsible for
 * cleaning up the string's memory.
 */
char *endInstance_getCompleteName(EndInstance *endInstance);

/*
 * Gets the encompassing end.
 */
End *endInstance_getEnd(EndInstance *endInstance);

/*
 * Gets the atom instance associated with the end, or NULL, if the end has no associated atom end at this level.
 * The atom instance returned will have the end instance on its left side.
 */
AtomInstance *endInstance_getAtomInstance(EndInstance *endInstance);

/*
 * Gets the coordinate of the end instance, returns INT32_MAX if coordinate not set.
 */
int32_t endInstance_getCoordinate(EndInstance *endInstance);

/*
 * Gets the sequence in which the instance exists, or NULL if not set.
 */
Sequence *endInstance_getSequence(EndInstance *endInstance);

/*
 * Sets adjacent end instances (this will set the adjacency reciprocally).
 * Any previous adjacency will be set to NULL for both ends.
 */
void endInstance_makeAdjacent1(EndInstance *endInstance, EndInstance *endInstance2);

/*
 * Sets alternatively adjacent end instances (this will set the adjacency reciprocally).
 * Any previous alternative adjacency will be set to NULL for both ends.
 */
void endInstance_makeAdjacent2(EndInstance *endInstance, EndInstance *endInstance2);

/*
 * Gets the adjacent end instance.
 */
EndInstance *endInstance_getAdjacency(EndInstance *endInstance);

/*
 * Gets the alternative adjacency.
 */
EndInstance *endInstance_getAdjacency2(EndInstance *endInstance);

/*
 * Gets any operation associated the end instance.
 */
Operation *endInstance_getOperation(EndInstance *endInstance);

/*
 * Gets the parent end instance (in the tree of the end).
 */
EndInstance *endInstance_getParent(EndInstance *endInstance);

/*
 * Returns the number of children the end instance has.
 */
int32_t endInstance_getChildNumber(EndInstance *endInstance);

/*
 * Gets the child end instance in the tree of the end.
 */
EndInstance *endInstance_getChild(EndInstance *endInstance, int32_t index);

/*
 * Links together a parent and child end instance.
 */
void endInstance_linkParentAndChild(EndInstance *endInstanceParent, EndInstance *endInstanceChild);

/*
 * Returns non zero if the end instance is internal (part of an internal tree).
 */
int32_t endInstance_isInternal(EndInstance *endInstance);

/*
 * Returns non zero if the end instance is augmented (added to accommodate an adjacency, but without an attached atom instance).
 */
int32_t endInstance_isAugmented(EndInstance *endInstance);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic end functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs the end, but not any attached atom.
 */
End *end_construct(const char *name, Net *net);

/*
 * Copies the end, but not any attached atom. Replaces the net attached to the end with the given
 * 'newNet'.
 */
End *end_copyConstruct(End *end, Net *newNet);

/*
 *	Name of the end.
 */
const char *end_getName(End *end);

/*
 * Gets the net the end is part of.
 */
Net *end_getNet(End *end);

/*
 * Gets the atom the end is on the left of.
 */
Atom *end_getAtom(End *end);

/*
 * Gets the adjacency component that the end is part of.
 */
AdjacencyComponent *end_getAdjacencyComponent(End *end);

/*
 * Returns the number of end instances the end contains.
 */
int32_t end_getInstanceNumber(End *end);

/*
 * Gets an instance using its instance name as a key. Instance name is m of full name n.m.
 */
EndInstance *end_getInstance(End *end, const char *instanceName);

/*
 * Gets the first instance in the end, or NULL if none.
 */
EndInstance *end_getFirst(End *end);

/*
 * Gets an iterator over the end instances.
 */
End_InstanceIterator *end_getInstanceIterator(End *end);

/*
 * Gets the next end instance from the iterator.
 */
EndInstance *end_getNext(End_InstanceIterator *iterator);

/*
 * Gets the previous end instance from the iterator.
 */
EndInstance *end_getPrevious(End_InstanceIterator *iterator);

/*
 * Duplicates the iterator.
 */
End_InstanceIterator *end_copyInstanceIterator(End_InstanceIterator *iterator);

/*
 * Destructs the iterator.
 */
void end_destructInstanceIterator(End_InstanceIterator *end);

/*
 * Return non zero if the end represents a stub.
 */
int32_t end_isStub(End *end);

/*
 * Return non zero if the end represents a cap (which is an atom end or psuedo telomere from a higher level).
 */
int32_t end_isCap(End *end);

/*
 * Return non zero if the end represents the end of an atom represented in the net of the end instance (at this level).
 */
int32_t end_isAtomEnd(End *end);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic atom instance functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs atom instance with two end instances. Instance is the suffix m of the instance name n.m.
 */
AtomInstance *atomInstance_construct(const char *instance, Atom *atom);

/*
 * As default constructor, but also sets the instance's coordinates.
 */
AtomInstance *atomInstance_constructWithCoordinates(const char *instance, Atom *atom, int32_t startCoordinate, Sequence *sequenceName);

/*
 * Gets the encompassing atom.
 */
Atom *atomInstance_getAtom(AtomInstance *atomInstance);

/*
 * Returns the m part of an instance's n.m name.
 */
const char *atomInstance_getInstanceName(AtomInstance *atomInstance);

/*
 * Returns the n part of an instance's n.m name.
 */
const char *atomInstance_getElementName(AtomInstance *atomInstance);

/*
 * Gets the complete name of an instance. This involves a new memory allocation, you are therefore responsible for
 * cleaning up the string's memory.
 */
char *atomInstance_getCompleteName(AtomInstance *atomInstance);

/*
 * Gets the reverse atom instance, giving a reversed view of the atom instance.
 */
AtomInstance *atomInstance_getReverse(AtomInstance *atomInstance);

/*
 * Gets the start coordinate of the atom instance, returns INT32_MAX if coordinate not set.
 */
int32_t atomInstance_getStart(AtomInstance *atomInstance);

/*
 * Gets the length of the atom instance.
 */
int32_t atomInstance_getLength(AtomInstance *atomInstance);

/*
 * Gets the sequence in which the instance exists, or NULL if not set.
 */
Sequence *atomInstance_getSequence(AtomInstance *atomInstance);

/*
 * Gets the left end instance of the atom instance.
 */
EndInstance *atomInstance_getLeft(AtomInstance *atomInstance);

/*
 * Gets the right end instance of the atom instance.
 */
EndInstance *atomInstance_getRight(AtomInstance *atomInstance);

/*
 * Returns the parent instance, or NULL, if none exists.
 */
AtomInstance *atomInstance_getParent(AtomInstance *atomInstance);

/*
 * Returns the number of children the instance has.
 */
int32_t atomInstance_getChildNumber(AtomInstance *atomInstance);

/*
 * Gets a child instance.
 */
AtomInstance *atomInstance_getChild(AtomInstance *atomInstance, int32_t index);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic atom functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs the atom, but not its ends.
 */
Atom *atom_construct(const char *name, int32_t length, Net *net);

/*
 * Returns string name of the atom.
 */
const char *atom_getName(Atom *atom);

/*
 * Returns the length in bases of the atom.
 */
int32_t atom_getLength(Atom *atom);

/*
 * Gets the net the atom is part of.
 */
Net *atom_getNet(Atom *atom);

/*
 * Gets the left end of the atom.
 */
End *atom_getLeft(Atom *atom);

/*
 * Gets the right end of the atom.
 */
End *atom_getRight(Atom *atom);

/*
 * Returns a reversed view of the atom.
 */
Atom *atom_getReverse(Atom *atom);

/*
 * Returns the number of instances (including any internal instances), the atom contains.
 */
int32_t atom_getInstanceNumber(Atom *atom);

/*
 * Gets the atom instance using its instance name as a key. Instance name is m of full name n.m.
 */
AtomInstance *atom_getInstance(Atom *atom, const char *instanceName);

/*
 * Gets the first atom instance in the list.
 */
AtomInstance *atom_getFirst(Atom *atom);

/*
 * Gets an iterator to iterate over the atom instances.
 */
Atom_InstanceIterator *atom_getInstanceIterator(Atom *atom);

/*
 * Gets the next atom instance, or NULL if none remaining.
 */
AtomInstance *atom_getNext(Atom_InstanceIterator *iterator);

/*
 * Gets the previous atom instance, or NULL if none remaining.
 */
AtomInstance *atom_getPrevious(Atom_InstanceIterator *iterator);

/*
 * Duplicates the iterator.
 */
Atom_InstanceIterator *atom_copyInstanceIterator(Atom_InstanceIterator *iterator);

/*
 * Destructs the iterator - should always be coupled with the iterator.
 */
void atom_destructInstanceIterator(Atom_InstanceIterator *atom);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic adjacency component functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an adjacency component.
 */
AdjacencyComponent *adjacencyComponent_construct(Net *net, Net *nestedNet);

/*
 * Updates the adjacency component's set of ends to contain the intersection of ends
 * contained in both the parent net and the nested net of the adjacency component.
 */
void adjacencyComponent_updateContainedEnds(AdjacencyComponent *adjacencyComponent);

/*
 *  Gets the net the adjacency component is part of.
 */
Net *adjacencyComponent_getNet(AdjacencyComponent *adjacencyComponent);

/*
 * Gets the name of the nested net the adjacency component contains.
 */
const char *adjacencyComponent_getNestedNetName(AdjacencyComponent *adjacencyComponent);

/*
 * Gets the nested net the adjacency component contains.
 */
Net *adjacencyComponent_getNestedNet(AdjacencyComponent *adjacencyComponent);

/*
 * Gets the chain the adjacency component is part of, or NULL, if not part of a chain.
 */
Chain *adjacencyComponent_getChain(AdjacencyComponent *adjacencyComponent);

/*
 * Gets an end by name
 */
End *adjacencyComponent_getEnd(AdjacencyComponent *adjacencyComponent, const char *name);

/*
 * Returns the number of ends.
 */
int32_t adjacencyComponent_getEndNumber(AdjacencyComponent *adjacencyComponent);

/*
 * Gets an iterator to iterate through the ends in the adjacency component.
 */
AdjacencyComponent_EndIterator *adjacencyComponent_getEndIterator(AdjacencyComponent *adjacencyComponent);

/*
 * Gets the next end from the iterator.
 */
End *adjacencyComponent_getNextEnd(AdjacencyComponent_EndIterator *endIterator);

/*
 * Gets the previous end from the iterator.
 */
End *adjacencyComponent_getPreviousEnd(AdjacencyComponent_EndIterator *endIterator);

/*
 * Duplicates the iterator.
 */
AdjacencyComponent_EndIterator *adjacencyComponent_copyEndIterator(AdjacencyComponent_EndIterator *endIterator);

/*
 * Destructs the iterator.
 */
void adjacencyComponent_destructEndIterator(AdjacencyComponent_EndIterator *endIterator);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic link functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Construct a link.
 */
Link *link_construct(End *leftEnd, End *rightEnd, AdjacencyComponent *adjacencyComponent, Chain *parentChain);

/*
 * Gets the next link in the link.
 */
Link *link_getNextLink(Link *link);

/*
 * Gets the prior link in the link.
 */
Link *link_getPreviousLink(Link *link);

/*
 * Gets the nested net the link contains.
 */
AdjacencyComponent *link_getAdjacencyComponent(Link *link);

/*
 * Gets the left end of the link in the link.
 */
End *link_getLeft(Link *link);

/*
 * Gets the right end of the link in the link.
 */
End *link_getRight(Link *link);

/*
 * Gets the chain the link is part of.
 */
Chain *link_getChain(Link *link);

/*
 * Gets the index of the link in the chain..
 */
int32_t link_getIndex(Link *link);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic chain functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a chain, which in turn holds links.
 */
Chain *chain_construct(Net *net);

/*
 * Gets a link in the chain.
 */
Link *chain_getLink(Chain *chain, int32_t linkIndex);

/*
 * Returns the number of links in the chain.
 */
int32_t chain_getLength(Chain *chain);

/*
 * Gets the index of the chain in the net.
 */
int32_t chain_getIndex(Chain *chain);

/*
 * Gets the parent net of the chain.
 */
Net *chain_getNet(Chain *chain);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic operation functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an operation.
 */
Operation *operation_construct(Net *net);

/*
 * Gets the net it is part of.
 */
Net *operation_getNet(Operation *opetation);

/*
 * Get the index of the operation.
 */
int32_t operation_getIndex(Operation *operation);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs the net.
 */
Net *net_construct(const char *name, NetDisk *netDisk);

/*
 * Destructs the net, and all the elements it contains. If recursive the function will destroy all
 * loaded nested nets.
 */
void net_destruct(Net *net, int32_t recursive);

/*
 * Gets the name of the net.
 */
const char *net_getName(Net *net);

/*
 * Gets the parent net disk.
 */
NetDisk *net_getNetDisk(Net *net);

/*
 * Adds the sequence to the net.
 */
void net_addSequence(Net *net, Sequence *sequence);

/*
 * Gets the 'first' sequence.
 */
Sequence *net_getFirstSequence(Net *net);

/*
 * Gets an sequence by its name.
 */
Sequence *net_getSequence(Net *net, const char *name);

/*
 * Returns the number of sequences.
 */
int32_t net_getSequenceNumber(Net *net);

/*
 * Gets an iterator to iterate through the chains in the net, at this level.
 */
Net_SequenceIterator *net_getSequenceIterator(Net *net);

/*
 * Gets the next sequence from the iterator.
 */
Sequence *net_getNextSequence(Net_SequenceIterator *sequenceIterator);

/*
 * Gets the previous sequence from the iterator.
 */
Sequence *net_getPreviousSequence(Net_SequenceIterator *sequenceIterator);

/*
 * Duplicates the iterator.
 */
Net_SequenceIterator *net_copySequenceIterator(Net_SequenceIterator *sequenceIterator);

/*
 * Destructs the iterator.
 */
void net_destructSequenceIterator(Net_SequenceIterator *sequenceIterator);

/*
 * Gets the 'first' end.
 */
End *net_getFirstEnd(Net *net);

/*
 * Gets an end by name.
 */
End *net_getEnd(Net *net, const char *name);

/*
 * Returns the number of ends.
 */
int32_t net_getEndNumber(Net *net);

/*
 * Gets an iterator to iterate through the ends in the net, at this level.
 */
Net_EndIterator *net_getEndIterator(Net *net);

/*
 * Gets the next end from the iterator.
 */
End *net_getNextEnd(Net_EndIterator *endIterator);

/*
 * Gets the previous end from the iterator.
 */
End *net_getPreviousEnd(Net_EndIterator *endIterator);

/*
 * Duplicates the iterator.
 */
Net_EndIterator *net_copyEndIterator(Net_EndIterator *endIterator);

/*
 * Destructs the iterator.
 */
void net_destructEndIterator(Net_EndIterator *endIterator);

/*
 * Gets the 'first' atom.
 */
Atom *net_getFirstAtom(Net *net);

/*
 * Gets an atom by name.
 */
Atom *net_getAtom(Net *net, const char *name);

/*
 * Returns the number of atoms.
 */
int32_t net_getAtomNumber(Net *net);

/*
 * Gets an iterator to iterate through the atoms in the net, at this level.
 */
Net_AtomIterator *net_getAtomIterator(Net *net);

/*
 * Gets the next atom from the iterator.
 */
Atom *net_getNextAtom(Net_AtomIterator *atomIterator);

/*
 * Gets the previous atom from the iterator.
 */
Atom *net_getPreviousAtom(Net_AtomIterator *atomIterator);

/*
 * Duplicates the iterator
 */
Net_AtomIterator *net_copyAtomIterator(Net_AtomIterator *atomIterator);

/*
 * Destructs the iterator.
 */
void net_destructAtomIterator(Net_AtomIterator *atomIterator);

/*
 * Gets the 'first' adjacency component.
 */
AdjacencyComponent *net_getFirstAdjacencyComponent(Net *net);

/*
 * Gets an adjacency component by the name of the nested net it contains.
 */
AdjacencyComponent *net_getAdjacencyComponent(Net *net, const char *netName);

/*
 * Returns the number of adjacency components.
 */
int32_t net_getAdjacencyComponentNumber(Net *net);

/*
 * Gets an iterator to iterate through the adjacency components in the net, at this level.
 */
Net_AdjacencyComponentIterator *net_getAdjacencyComponentIterator(Net *net);

/*
 * Gets the next adjacency component from the iterator.
 */
AdjacencyComponent *net_getNextAdjacencyComponent(Net_AdjacencyComponentIterator *adjacencyComponentIterator);

/*
 * Gets the previous adjacency component from the iterator.
 */
AdjacencyComponent *net_getPreviousAdjacencyComponent(Net_AdjacencyComponentIterator *adjacencyComponentIterator);

/*
 * Duplicates the iterator.
 */
Net_AdjacencyComponentIterator *net_copyAdjacencyComponentIterator(Net_AdjacencyComponentIterator *adjacencyComponentIterator);

/*
 * Destructs the iterator.
 */
void net_destructAdjacencyComponentIterator(Net_AdjacencyComponentIterator *adjacencyComponentIterator);

/*
 * Gets the parent adjacency component of the net.
 */
AdjacencyComponent *net_getParentAdjacencyComponent(Net *net);

/*
 * Gets the 'first' chain.
 */
Chain *net_getFirstChain(Net *net);

/*
 * Gets a chain by its index
 */
Chain *net_getChain(Net *net, int32_t index);

/*
 * Returns the number of chains.
 */
int32_t net_getChainNumber(Net *net);

/*
 * Gets an iterator to iterate through the chains in the net, at this level.
 */
Net_ChainIterator *net_getChainIterator(Net *net);

/*
 * Gets the next chain from the iterator.
 */
Chain *net_getNextChain(Net_ChainIterator *chainIterator);

/*
 * Gets the previous chain from the iterator.
 */
Chain *net_getPreviousChain(Net_ChainIterator *chainIterator);

/*
 * Duplicates the iterator.
 */
Net_ChainIterator *net_copyChainIterator(Net_ChainIterator *chainIterator);

/*
 * Destructs the iterator.
 */
void net_destructChainIterator(Net_ChainIterator *chainIterator);

/*
 * Gets the 'first' operation.
 */
Operation *net_getFirstOperation(Net *net);

/*
 * Gets an chain by index.
 */
Operation *net_getOperation(Net *net, int32_t index);

/*
 * Returns the number of operations.
 */
int32_t net_getOperationNumber(Net *net);

/*
 * Gets an iterator to iterate through the operations in the net, at this level.
 */
Net_OperationIterator *net_getOperationIterator(Net *net);

/*
 * Gets the next operation from the iterator.
 */
Operation *net_getNextOperation(Net_OperationIterator *operationIterator);

/*
 * Gets the previous operation from the iterator.
 */
Operation *net_getPreviousOperation(Net_OperationIterator *operationIterator);

/*
 * Duplicates the iterator.
 */
Net_OperationIterator *net_copyOperationIterator(Net_OperationIterator *operationIterator);

/*
 * Destructs the iterator.
 */
void net_destructOperationIterator(Net_OperationIterator *operationIterator);


/*
 * Creates a binary representation of the net, returned as a char string.
 */
char *net_makeBinaryRepresentation(Net *net);

/*
 * Creates an XML representation of the net, returned as a char string.
 */
char *net_makeXMLRepresentation(Net *net);

/*
 * Loads a net into memory from a binary representation of the net.
 */
Net *net_loadFromBinaryRepresentation(char *binaryString, NetDisk *netDisk);

/*
 * Loads a net into memory from an XML representation of the net.
 */
Net *net_loadFromXMLRepresentation(char *xmlString, NetDisk *netDisk);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a net disk to load nets.
 */
NetDisk *netDisk_construct(const char *netDiskFile);

/*
 * Destructs the net disk, and all open nets and sequences.
 */
void netDisk_destruct(NetDisk *netDisk);

/*
 * Writes the updated state of the parts of the net disk in memory to disk.
 * Returns 0 for success, non-zero for failure.
 */
int32_t netDisk_write(NetDisk *netDisk);

/*
 * Gets a sequence the netDisk contains. If the sequence is not in memory it will be loaded. If not in memory or on disk, returns NULL.
 */
Sequence *netDisk_getSequence(NetDisk *netDisk, const char *sequenceName);

/*
 * Returns the number of sequences on disk.
 */
int32_t netDisk_getSequenceNumberOnDisk(NetDisk *netDisk);

/*
 * Gets an iterator to iterate through the sequence names currently on disk.
 */
NetDisk_SequenceNameIterator *netDisk_getSequenceNameIterator(NetDisk *netDisk);

/*
 * Gets the next sequence name from the iterator.
 */
const char *netDisk_getNextSequenceName(NetDisk_SequenceNameIterator *sequenceIterator);

/*
 * Destructs the iterator.
 */
void netDisk_destructSequenceNameIterator(NetDisk_SequenceNameIterator *sequenceIterator);

/*
 * Gets a sequence the netDisk contains that is currently in memory. Returns NULL if not in memory.
 */
Sequence *netDisk_getSequenceInMemory(NetDisk *netDisk, const char *sequenceName);

/*
 * Gets the first sequence in the list of sequences in memory, or returns NULL if the list is empty.
 */
Sequence *netDisk_getFirstSequenceInMemory(NetDisk *netDisk);

/*
 * Returns the number of sequences currently in memory.
 */
int32_t netDisk_getSequenceNumberInMemory(NetDisk *netDisk);

/*
 * Gets an iterator to iterate through the sequences currently in memory.
 */
NetDisk_SequenceIterator *netDisk_getSequenceInMemoryIterator(NetDisk *netDisk);

/*
 * Gets the next sequence from the iterator.
 */
Sequence *netDisk_getNextSequence(NetDisk_SequenceIterator *sequenceIterator);

/*
 * Gets the previous sequence from the iterator.
 */
Sequence *netDisk_getPreviousSequence(NetDisk_SequenceIterator *sequenceIterator);

/*
 * Duplicates the iterator.
 */
NetDisk_SequenceIterator *netDisk_copySequenceIterator(NetDisk_SequenceIterator *sequenceIterator);

/*
 * Destructs the iterator.
 */
void netDisk_destructSequenceIterator(NetDisk_SequenceIterator *sequenceIterator);

/*
 * Gets a net the netDisk contains. If the net is not in memory it will be loaded. If not in memory or on disk, returns NULL.
 */
Net *netDisk_getNet(NetDisk *netDisk, const char *netName);

/*
 * Returns the number of nets on disk.
 */
int32_t netDisk_getNetNumberOnDisk(NetDisk *netDisk);

/*
 * Gets an iterator to iterate through the net names currently on disk.
 */
NetDisk_NetNameIterator *netDisk_getNetNameIterator(NetDisk *netDisk);

/*
 * Gets the next net name from the iterator.
 */
const char *netDisk_getNextNetName(NetDisk_NetNameIterator *netIterator);

/*
 * Destructs the iterator.
 */
void netDisk_destructNetNameIterator(NetDisk_NetNameIterator *netIterator);

/*
 * Gets a net the netDisk contains that is currently in memory. Returns NULL if not in memory.
 */
Net *netDisk_getNetInMemory(NetDisk *netDisk, const char *netName);

/*
 * Gets the first net in the list of nets in memory, or returns NULL if the list is empty.
 */
Net *netDisk_getFirstNetInMemory(NetDisk *netDisk);

/*
 * Returns the number of nets currently in memory.
 */
int32_t netDisk_getNetNumberInMemory(NetDisk *netDisk);

/*
 * Gets an iterator to iterate through the nets currently in memory.
 */
NetDisk_NetIterator *netDisk_getNetInMemoryIterator(NetDisk *netDisk);

/*
 * Gets the next net from the iterator.
 */
Net *netDisk_getNextNet(NetDisk_NetIterator *netIterator);

/*
 * Gets the previous net from the iterator.
 */
Net *netDisk_getPreviousNet(NetDisk_NetIterator *netIterator);

/*
 * Duplicates the iterator.
 */
NetDisk_NetIterator *netDisk_copyNetIterator(NetDisk_NetIterator *netIterator);

/*
 * Destructs the iterator.
 */
void netDisk_destructNetIterator(NetDisk_NetIterator *netIterator);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Useful utility functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * The m part of a complete instance name n.m is returned. This involves memory allocation,
 * you are responsible for cleaning up the memory.
 */
char *netMisc_getInstanceName(const char *completeName);

/*
 * The m part of a complete instance name n.m is returned.
 * The memory for the string is owned by the function, so you needn't clean it up.
 * However, this memory will be overwritten with each call to the function.
 */
const char *netMisc_getInstanceNameStatic(const char *completeName);

/*
 * The n part of a complete instance name n.m is returned. This involves memory allocation,
 * you are responsible for cleaning up the memory.
 */
char *netMisc_getElementName(const char *completeName);

/*
 * The n part of a complete instance name n.m is returned.
 * The memory for the string is owned by the function, so you needn't clean it up.
 * However, this memory will be overwritten with each call to the function.
 */
const char *netMisc_getElementNameStatic(const char *completeName);

/*
 * Concatenates an element name n and instance name m to form a complete name of the form n.m .
 * This involves memory allocation, you are responsible for cleaning up the memory.
 */
char *netMisc_makeCompleteName(const char *elementName, const char *instanceName);

/*
 * Concatenates an element name n and instance name m to form a complete name of the form n.m .
 * The memory for the string is owned by the function, so you needn't clean it up.
 * However, this memory will be overwritten with each call to the function.
 */
const char *netMisc_makeCompleteNameStatic(const char *elementName, const char *instanceName);

#endif
