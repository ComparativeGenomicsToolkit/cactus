/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_FLOWER_PRIVATE_H_
#define CACTUS_FLOWER_PRIVATE_H_

#include "cactusGlobals.h"

struct _flower {
    Name name;
    stList *ends;
    stList *caps;
    stList *groups;
    stList *chains;
    stList *sequences;
    Name parentFlowerName;
    CactusDisk *cactusDisk;
    bool builtBlocks;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private flower functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Adds the event tree for the flower to the flower.
 * If an previous event tree exists for the flower
 * it will call eventTree_destruct for the existing tree
 * (which should not exist without the flower).
 */
void flower_setEventTree(Flower *flower, EventTree *eventTree);

/*
 * This function is called by eventTree_destruct and cleans up the reference.
 */
void flower_removeEventTree(Flower *flower, EventTree *eventTree);

/*
 * Adds the cap to the flower.
 */
void flower_addCap(Flower *flower, Cap *cap);

/*
 * Bulk add a set of ends to the flower.
 */
void flower_bulkAddCaps(Flower *flower, stList *capsToAdd);

/*
 * Adds the end to the flower.
 */
void flower_addEnd(Flower *flower, End *end);

/*
 * Bulk add a set of ends to the flower.
 */
void flower_bulkAddEnds(Flower *flower, stList *endsToAdd);

/*
 * Remove the end from the flower.
 */
void flower_removeEnd(Flower *flower, End *end);

/*
 * Adds the group to the flower.
 */
void flower_addGroup(Flower *flower, Group *group);

/*
 * Removes an empty group from the flower.
 */
void flower_removeGroup(Flower *flower, Group *group);

/*
 * Sets the parent group of the flower.
 */
void flower_setParentGroup(Flower *flower, Group *group);

/*
 * Adds the chain to the flower.
 */
void flower_addChain(Flower *flower, Chain *chain);

/*
 * Remove the chain from the flower.
 */
void flower_removeChain(Flower *flower, Chain *chain);

#endif
