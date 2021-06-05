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
	Link *nLink;
	Name name;
	End *firstEnd; // If a link, this becomes the 5end and the second is the 3end in the chain
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
 * Removes the end from the group.
 */
void group_removeEnd(Group *group, End *end);

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
