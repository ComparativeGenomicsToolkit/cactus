#ifndef CACTUS_ADJACENCY_COMPONENT_PRIVATE_H_
#define CACTUS_ADJACENCY_COMPONENT_PRIVATE_H_

#include "cactusGlobals.h"

struct _group {
	Net *net;
	Link *link;
	Name name;
	struct avl_table *ends;
	bool terminalGroup;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Group functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a nested net without having the nested net loaded in memory.
 */
Group *group_construct3(Net *net, Name nestedNetName, bool terminalGroup);

/*
 * Destructs an group.
 */
void group_destruct(Group *group);

/*
 * Updates the group's set of ends to contain the intersection of ends
 * contained in both the parent net and the nested net of the group.
 *
 * This function will create an assertion error if the group is terminal.
 */
void group_updateContainedEnds(Group *group);

/*
 * Sets the link the group is part of.
 */
void group_setLink(Group *group, Link *link);

/*
 * Write a binary representation of the group to the write function.
 */
void group_writeBinaryRepresentation(Group *group, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Group *group_loadFromBinaryRepresentation(void **binaryString, Net *net);

/*
 * Get a static instance (from the heap) with the netName set.
 */
Group *group_getStaticNameWrapper(Name netName);

#endif
