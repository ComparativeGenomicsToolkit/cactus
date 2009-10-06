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

typedef struct _metaEvent MetaEvent;
typedef struct _event Event;
typedef struct _eventTree EventTree;
typedef struct _metaSequence MetaSequence;
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

typedef struct avl_traverser EventTree_Iterator;
typedef struct _end_instanceIterator End_InstanceIterator;
typedef struct _atom_instanceIterator Atom_InstanceIterator;
typedef struct avl_traverser AdjacencyComponent_EndIterator;
typedef struct avl_traverser Net_SequenceIterator;
typedef struct avl_traverser Net_EndIterator;
typedef struct avl_traverser Net_AtomIterator;
typedef struct avl_traverser Net_AdjacencyComponentIterator;
typedef struct avl_traverser Net_ChainIterator;
typedef struct avl_traverser Net_OperationIterator;
typedef BDBCUR NetDisk_NetNameIterator;
typedef struct avl_traverser NetDisk_NetIterator;

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a meta event, which contains all the essential info for a event.
 */
MetaEvent *metaEvent_construct(const char *name, const char *header,
		NetDisk *netDisk);

/*
 * Gets the name of the event.
 */
const char *metaEvent_getName(MetaEvent *metaEvent);

/*
 * Gets the header line associated with the meta event.
 */
const char *metaEvent_getHeader(MetaEvent *metaEvent);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an event, attached on a branch from the parent event, with the given branch length.
 */
Event *event_construct(MetaEvent *metaEvent, float branchLength, Event *parentEvent, EventTree *eventTree);

/*
 * Constructs an event, on the branch between the given child and parent events. The branch length is the
 * length from the parent to the new event. In this case it should be less than the length from the parent to
 * the pre-existing child event.
 */
Event *event_construct2(MetaEvent *metaEvent, float branchLength,
		Event *parentEvent, Event *childEvent, EventTree *eventTree);

/*
 * Returns the parent event, or NULL if root.
 */
Event *event_getParent(Event *event);

/*
 * Gets the name of the event.
 */
const char *event_getName(Event *event);

/*
 * Gets the associated meta event.
 */
MetaEvent *event_getMetaEvent(Event *event);

/*
 * Gets the header sequence associated with the event.
 */
const char *event_getHeader(Event *event);

/*
 * Gets the branch length.
 */
float event_getBranchLength(Event *event);

/*
 * Gets the branch length of the subtree rooted at this event, excluding the branch length of the event
 * itself.
 */
float event_getSubTreeBranchLength(Event *event);

/*
 * Gets the number of events in the sub tree of the event, excluding the event itself.
 */
int32_t event_getSubTreeEventNumber(Event *event);

/*
 * Get number of children.
 */
int32_t event_getChildNumber(Event *event);

/*
 * Gets a child by its index.
 */
Event *event_getChild(Event *event, int32_t index);

/*
 * Gets the event tree the event is part of.
 */
EventTree *event_getEventTree(Event *event);

/*
 * Returns non zero if the other event is an ancestor of the event.
 */
int32_t event_isAncestor(Event *event, Event *otherEvent);

/*
 * Returns non zero if the other event is a child of the event.
 */
int32_t event_isDescendant(Event *event, Event *otherEvent);

/*
 * Returns non zero if the other event is a sibling of the other event.
 */
int32_t event_isSibling(Event *event, Event *otherEvent);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic event tree functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an event tree, with one root event.
 */
EventTree *eventTree_construct(MetaEvent *rootMetaEvent, Net *net);

/*
 * Copy constructs the event tree, replacing the existing net with the newNet. Only includes
 * the unary events which are true, given the unary event function.
 */
EventTree *eventTree_copyConstruct(EventTree *eventTree, Net *newNet, int32_t (unaryEventFilterFn)(Event *event));

/*
 * Returns the root event.
 */
Event *eventTree_getRootEvent(EventTree *eventTree);

/*
 * Gets the event with the given name.
 */
Event *eventTree_getEvent(EventTree *eventTree, const char *eventName);

/*
 * Gets the common ancestor of two events.
 */
Event *eventTree_getCommonAncestor(Event *event, Event *event2);

/*
 * Gets the parent net.
 */
Net *eventTree_getNet(EventTree *eventTree);

/*
 * Gets the total number of events in the event tree.
 */
int32_t eventTree_getEventNumber(EventTree *eventTree);

/*
 * Gets the first event in the list.
 */
Event *eventTree_getFirst(EventTree *eventTree);

/*
 * Gets an iterator over the eventTree events.
 */
EventTree_Iterator *eventTree_getIterator(EventTree *eventTree);

/*
 * Gets the next event from the iterator.
 */
Event *eventTree_getNext(EventTree_Iterator *iterator);

/*
 * Gets the previous event from the iterator.
 */
Event *eventTree_getPrevious(EventTree_Iterator *iterator);

/*
 * Duplicates the iterator.
 */
EventTree_Iterator *eventTree_copyIterator(EventTree_Iterator *iterator);

/*
 * Destructs the iterator.
 */
void eventTree_destructIterator(EventTree_Iterator *iterator);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


/*
 * Constructs a meta sequence, which contains all the essential info for a sequence.
 *
 * This function is NOT thread safe, do not try to have concurrent instances of this function!
 */
MetaSequence *metaSequence_construct(const char *name, int32_t start, int32_t length, const char *string, const char *header,
		const char *eventName, NetDisk *netDisk);

/*
 * Gets the name of the sequence.
 */
const char *metaSequence_getName(MetaSequence *metaSequence);

/*
 * Gets the start coordinate of the sequence.
 */
int32_t metaSequence_getStart(MetaSequence *metaSequence);

/*
 * Gets the length of the sequence.
 */
int32_t metaSequence_getLength(MetaSequence *metaSequence);

/*
 * Gets the associated event name.
 */
const char *metaSequence_getEventName(MetaSequence *metaSequence);

/*
 * Gets a string for representing a subsequence of the meta sequence.
 */
char *metaSequence_getString(MetaSequence *metaSequence, int32_t start, int32_t length, int32_t strand);

/*
 * Gets the header line associated with the meta sequence.
 */
const char *metaSequence_getHeader(MetaSequence *metaSequence);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Creates a sequence for a net, wrapping a meta sequence.
 */
Sequence *sequence_construct(MetaSequence *metaSequence, Net *net);

/*
 * Destructs the sequence.
 */
void sequence_destruct(Sequence *sequence);

/*
 * Gets the associated meta sequence.
 */
MetaSequence *sequence_getMetaSequence(Sequence *sequence);

/*
 * Gets the start coordinate of the sequence (inclusive).
 */
int32_t sequence_getStart(Sequence *sequence);

/*
 * Gets the length of the sequence.
 */
int32_t sequence_getLength(Sequence *sequence);

/*
 * Gets the name of the sequence.
 */
const char *sequence_getName(Sequence *sequence);

/*
 * Gets the name associated with the sequence.
 */
Event *sequence_getEvent(Sequence *sequence);

/*
 * Gets a sub string of the the sequence, indexes must be equal to or greater than the start coordinate,
 * and less than the start coordinate plus the sequences length.
 * If the strand is negative then the sequence returned will be the reverse complement sequence, traversing
 * in the opposite direction.
 *
 * The returned string must be freed.
 */
char *sequence_getString(Sequence *sequence, int32_t start, int32_t length, int32_t strand);

/*
 * Gets the header line associated with the sequence.
 */
const char *sequence_getHeader(Sequence *sequence);

/*
 * Gets the net the sequence is associated with.
 */
Net *sequence_getNet(Sequence *sequence);

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
 * As default constructor, but also sets the instance's coordinates and event.
 */
EndInstance *endInstance_construct2(const char *instance, End *end, int32_t startCoordinate, int32_t strand, int32_t side, Sequence *sequence);

/*
 * Sets the event associated with the event. Will create an error if the event is already set.
 */
void endInstance_setEvent(EndInstance *endInstance, Event *event);

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
 * Returns a non zero integer if the end instance is oriented positively with respect to the end.
 * Else, returns zero if the end instance is oriented negatively with respect to the end.
 *
 * Thus, if you want to know the sign of instance, call this function and prefix a minus sign to the name if
 * this function returns zero.
 */
int32_t endInstance_getOrientation(EndInstance *endInstance);

/*
 * Gets the complete name of an instance, as in endInstance_getCompleteName, but prefixes
 * the name a sign if the instance is the reverse orientation.
 *
 * You must clean up the returned string.
 */
char *endInstance_getCompleteNameWithOrientation(EndInstance *endInstance);

/*
 * Gets the reversed end instance (the equivalent on the opposite strand, with the opposite orientation).
 */
EndInstance *endInstance_getReverse(EndInstance *endInstance);

/*
 * Gets the event associated with the endInstance.
 */
Event *endInstance_getEvent(EndInstance *endInstance);

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
 * Gets the coordinate of the position that end instance is on the end of,
 * returns INT32_MAX if coordinate not set.
 */
int32_t endInstance_getCoordinate(EndInstance *endInstance);

/*
 * Returns a non zero integer if the coordinate of the end instance (see endInstance_getCoordinate)
 * is on the forward strand, and zero if on the negative strand.
 */
int32_t endInstance_getStrand(EndInstance *endInstance);

/*
 * Returns a non zero integer if on the 5' side of the position returned by endInstance_getCoordinate,
 * zero if on the 3' side.
 */
int32_t endInstance_getSide(EndInstance *endInstance);

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
 *
 */
void endInstance_makeParentAndChild(EndInstance *endInstanceParent, EndInstance *endInstanceChild);

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
 * Returns a positive integer if the end is oriented positively.
 * Else, returns a negative integer.
 * The orientation is arbitrary (it is not explicitly with respect to anything else), but is consistent.
 *
 * All instances returned by the iterator of the end_getInstance() will have the same orientation as the
 * parent end.
 *
 * Thus, if you want to know the sign of end, call this function and prefix a minus sign to the name if
 * this function returns a negative integer.
 */
int32_t end_getOrientation(End *end);

/*
 *	Gets the name of the end, as in end_getName(), but prefixes a minus sign to the name
 *	if in the reverse orientation.
 *	You must clean up the returned string.
 */
char *end_getNameWithOrientation(End *end);

/*
 * Returns a reverse strand view of the end (in the opposite orientation).
 */
End *end_getReverse(End *end);

/*
 * Gets the net the end is part of.
 */
Net *end_getNet(End *end);

/*
 * Gets the atom the end is on the side of.
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
 * Gets the root end instance of the end, if it is set, or returns NULL;
 */
EndInstance *end_getRootInstance(End *end);

/*
 * Sets the root end instance of the end. Will throw an error if the endInstance
 * is not part of the end, or already has a parent.
 */
void end_setRootInstance(End *end, EndInstance *endInstance);

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
 * Constructs atom instance with the two end instances, which must both have the same instance name.
 */
AtomInstance *atomInstance_construct(Atom *atom,
		EndInstance *leftEndInstance, EndInstance *rightEndInstance);

/*
 * Constructs an atom instance and its two attached end instances.
 */
AtomInstance *atomInstance_construct2(const char *instance, Atom *atom);

/*
 * Constructs an atom instance and its two attached end instances, with the given coordinates.
 */
AtomInstance *atomInstance_construct3(const char *instance, Atom *atom,
		int32_t startCoordinate, int32_t strand, Sequence *sequence);

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
 * Returns a non zero integer if the instance is oriented positively with respect to the atom.
 * Else, returns zero if the instance is oriented negatively with respect to the atom.
 */
int32_t atomInstance_getOrientation(AtomInstance *atomInstance);

/*
 * Gets the complete name of an instance, as in atomInstance_getCompleteName, but prefixes
 * the name a sign if the instance is the reverse orientation.
 */
char *atomInstance_getCompleteNameWithOrientation(AtomInstance *atomInstance);

/*
 * Gets the reverse atom instance, giving a reversed view of the atom instance.
 */
AtomInstance *atomInstance_getReverse(AtomInstance *atomInstance);

/*
 * Gets the event associated with the instance.
 */
Event *atomInstance_getEvent(AtomInstance *atomInstance);

/*
 * Gets the start coordinate (that which is closest to the 5' end of the strand)
 *  of the atom instance, returns INT32_MAX if coordinate not set.
 */
int32_t atomInstance_getStart(AtomInstance *atomInstance);

/*
 * Returns non zero if one the forward strand, and zero if on the minus strand.
 */
int32_t atomInstance_getStrand(AtomInstance *atomInstance);

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
EndInstance *atomInstance_get5End(AtomInstance *atomInstance);

/*
 * Gets the right end instance of the atom instance.
 */
EndInstance *atomInstance_get3End(AtomInstance *atomInstance);

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
 * Returns a positive integer if the atom is oriented positively.
 * Else, returns a negative integer.
 * The orientation is arbitrary (it is not explicitly with respect to anything else), but is consistent.
 *
 * All instances returned by the iterator of the atom_getInstance() will have the same orientation as the
 * parent end.
 *
 * Thus, if you want to know the sign of end, call this function and prefix a minus sign to the name if
 * this function returns a negative integer.
 */
int32_t atom_getOrientation(Atom *atom);

/*
 *	Gets the name of the atom, as in atom_getName(), but prefixes a minus sign to the name
 *	if in the reverse orientation.
 *	You must clean up the returned string.
 */
char *atom_getNameWithOrientation(Atom *atom);

/*
 * Returns a reversed view of the atom (in the opposite orientation).
 */
Atom *atom_getReverse(Atom *atom);

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
End *atom_get5End(Atom *atom);

/*
 * Gets the right end of the atom.
 */
End *atom_get3End(Atom *atom);

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
Chain *chain_construct(Net *net, const char *name);

/*
 * Gets a link in the chain.
 */
Link *chain_getLink(Chain *chain, int32_t linkIndex);

/*
 * Returns the number of links in the chain.
 */
int32_t chain_getLength(Chain *chain);

/*
 * Gets the name of the chain in the net.
 */
const char *chain_getName(Chain *chain);

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
Operation *operation_construct(Net *net, const char *name);

/*
 * Gets the net it is part of.
 */
Net *operation_getNet(Operation *opetation);

/*
 * Get the name of the operation.
 */
const char *operation_getName(Operation *operation);

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
 * Gets the net tree associated with the event tree.
 */
EventTree *net_getEventTree(Net *net);

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
 * Gets an end instance by its complete name.
 */
EndInstance *net_getEndInstance(Net *net, const char *completeName);

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
 * Gets an atom instance by its complete name.
 */
AtomInstance *net_getAtomInstance(Net *net, const char *completeName);

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
Chain *net_getChain(Net *net, const char *name);

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
Operation *net_getOperation(Net *net, const char *name);

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
char *netMisc_makeCompleteName(const char *elementName, const char *instanceName, int32_t orientation);

/*
 * Concatenates an element name n and instance name m to form a complete name of the form n.m .
 * The memory for the string is owned by the function, so you needn't clean it up.
 * However, this memory will be overwritten with each call to the function.
 */
const char *netMisc_makeCompleteNameStatic(const char *elementName, const char *instanceName, int32_t orientation);

/*
 * Adds an orientation to a name. The returned string must be manually freed.
 */
char *netMisc_getNameWithOrientation(const char *name, int32_t orientation);

/*
 * Computes the reverse complement character of a ACTGactg, returning other characters unmodified.
 */
char netMisc_reverseComplementChar(char c);

/*
 * Computes the reverse complement of the string, returning the r-c string in newly allocated memory that must be freed.
 */
char *netMisc_reverseComplementString(const char *string);
#endif
