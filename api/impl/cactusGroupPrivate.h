/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_ADJACENCY_COMPONENT_PRIVATE_H_
#define CACTUS_ADJACENCY_COMPONENT_PRIVATE_H_

#include "cactusGlobals.h"

struct _group {
    void *flowerOrChain;
	//Flower *flower;
	//Link *link; // this becomes next link in the chain
	Link *nLink;
	Name name;
	//stSortedSet *ends;
	End *firstEnd; // If a link, this becomes the 5end and the second is the 3end in the chain
	//bool leafGroup; // this becomes an array of bools  including a flag indicating if it's a link
    char bits; // 0 bit: is leaf, 1 bit: is link
};

struct _group_endIterator {
    Group *group;
    End *end;
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
//void group_setLink(Group *group, Link *link);

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

/*
 * Set as a leaf group
 */
void group_setLeaf(Group *group, bool isLeaf);

/*
 * Set the group as a link
 */
void group_setLink(Group *group, bool isLink);

#endif
