/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_ADJACENCY_COMPONENT_PRIVATE_H_
#define CACTUS_ADJACENCY_COMPONENT_PRIVATE_H_

#include "cactusGlobals.h"

struct _group {
	Flower *flower;
	Link *link;
	Name name;
	stSortedSet *ends;
	bool leafGroup;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Group functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a nested flower without having the nested flower loaded in memory.
 */
Group *group_construct4(Flower *flower, Name nestedFlowerName, bool terminalGroup);

/*
 * Updates the group's set of ends to contain the intersection of ends
 * contained in both the parent flower and the nested flower of the group.
 *
 * This function will create an assertion error if the group is terminal.
 */
void group_updateContainedEnds(Group *group);

/*
 * Sets the link the group is part of.
 */
void group_setLink(Group *group, Link *link);

/*
 * Removes the end from the group.
 */
void group_removeEnd(Group *group, End *end);

/*
 * Write a binary representation of the group to the write function.
 */
void group_writeBinaryRepresentation(Group *group, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a flower into memory from a binary representation of the flower.
 */
Group *group_loadFromBinaryRepresentation(void **binaryString, Flower *flower);

/*
 * Sets the flower containing the group (the public function is end_setGroup).
 */
void group_setFlower(Group *group, Flower *flower);

/*
 * Adds an end to the group (the public function is end_setGroup).
 */
void group_addEnd(Group *group, End *end);

#endif
