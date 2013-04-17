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
    EventTree *eventTree;
    stSortedSet *sequences;
    stSortedSet *ends;
    stSortedSet *caps;
    stSortedSet *blocks;
    stSortedSet *segments;
    stSortedSet *groups;
    stSortedSet *chains;
    stSortedSet *faces;
    Name parentFlowerName;
    CactusDisk *cactusDisk;
    int64_t faceIndex;
    int64_t chainIndex;
    bool builtBlocks;
    bool builtTrees;
    bool builtFaces;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private flower functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructs the flower, and all the elements it contains. If recursive the function will destroy all
 * loaded nested flowers.
 */
void flower_destruct(Flower *flower, int64_t recursive);

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
 * Adds the sequence to the flower.
 */
void flower_addSequence(Flower *flower, Sequence *sequence);

/*
 * Removes the sequence from the flower.
 */
void flower_removeSequence(Flower *flower, Sequence *sequence);

/*
 * Adds the segment to the flower.
 */
void flower_addSegment(Flower *flower, Segment *segment);

/*
 * Remove the segment from the flower.
 */
void flower_removeSegment(Flower *flower, Segment *segment);

/*
 * Adds the block to the flower.
 */
void flower_addBlock(Flower *flower, Block *block);

/*
 * Remove the block from the flower.
 */
void flower_removeBlock(Flower *flower, Block *block);

/*
 * Adds the cap to the flower.
 */
void flower_addCap(Flower *flower, Cap *cap);

/*
 * Remove the cap from the flower.
 */
void flower_removeCap(Flower *flower, Cap *cap);

/*
 * Adds the end to the flower.
 */
void flower_addEnd(Flower *flower, End *end);

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

/*
 * Adds the face to the flower.
 */
void flower_addFace(Flower *flower, Face *face);

/*
 * Remove the face from the flower.
 */
void flower_removeFace(Flower *flower, Face *face);

/*
 * This function constructs faces for the flower. If faces are already created then
 * they will be first delete.
 */
void flower_reconstructFaces(Flower * flower);

/*
 * Destroys all faces in the flower.
 */
void flower_destructFaces(Flower *flower);

/*
 * Write a binary representation of the flower to the write function.
 */
void flower_writeBinaryRepresentation(Flower *flower, void(*writeFn)(const void * ptr,
        size_t size, size_t count));

/*
 * Loads a flower into memory from a binary representation of the flower.
 */
Flower *flower_loadFromBinaryRepresentation(void **binaryString, CactusDisk *cactusDisk);

#endif
