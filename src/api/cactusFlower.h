/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_FLOWER_H_
#define CACTUS_FLOWER_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic flower functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs the flower.
 */
Flower *flower_construct(CactusDisk *cactusDisk);

/*
 * Constructs the flower.
 */
Flower *flower_construct2(Name name, CactusDisk *cactusDisk);

/*
 * Gets the name of the flower.
 */
Name flower_getName(Flower *flower);

/*
 * Gets the parent flower disk.
 */
CactusDisk *flower_getCactusDisk(Flower *flower);

/*
 * Gets the flower tree associated with the event tree.
 */
EventTree *flower_getEventTree(Flower *flower);

/*
 * Gets the 'first' sequence.
 */
Sequence *flower_getFirstSequence(Flower *flower);

/*
 * Gets an sequence by its name.
 */
Sequence *flower_getSequence(Flower *flower, Name name);

/*
 * Returns the number of sequences.
 */
int64_t flower_getSequenceNumber(Flower *flower);

/*
 * Gets an iterator to iterate through the chains in the flower, at this level.
 */
Flower_SequenceIterator *flower_getSequenceIterator(Flower *flower);

/*
 * Gets the next sequence from the iterator.
 */
Sequence *flower_getNextSequence(Flower_SequenceIterator *sequenceIterator);

/*
 * Gets the previous sequence from the iterator.
 */
Sequence *flower_getPreviousSequence(Flower_SequenceIterator *sequenceIterator);

/*
 * Duplicates the iterator.
 */
Flower_SequenceIterator *flower_copySequenceIterator(
        Flower_SequenceIterator *sequenceIterator);

/*
 * Destructs the iterator.
 */
void flower_destructSequenceIterator(Flower_SequenceIterator *sequenceIterator);

/*
 *  Gets the 'first' cap.
 */
Cap *flower_getFirstCap(Flower *flower);

/*
 * Gets an cap by its  name.
 */
Cap *flower_getCap(Flower *flower, Name name);

/*
 * Returns the number of caps.
 */
int64_t flower_getCapNumber(Flower *flower);

/*
 * Gets an iterator to iterate through the caps in the flower, at this level.
 */
Flower_CapIterator *flower_getCapIterator(Flower *flower);

/*
 * Gets the next cap from the iterator.
 */
Cap *flower_getNextCap(Flower_CapIterator *capIterator);

/*
 * Gets the previous cap from the iterator.
 */
Cap *flower_getPreviousCap(Flower_CapIterator *capIterator);

/*
 * Duplicates the iterator.
 */
Flower_CapIterator *flower_copyCapIterator(Flower_CapIterator *capIterator);

/*
 * Destructs the iterator.
 */
void flower_destructCapIterator(Flower_CapIterator *capIterator);

/*
 * Gets the 'first' end.
 */
End *flower_getFirstEnd(Flower *flower);

/*
 * Gets an end by name.
 */
End *flower_getEnd(Flower *flower, Name name);

/*
 * Returns the number of ends.
 */
int64_t flower_getEndNumber(Flower *flower);

/*
 * Sugar for flower_getBlockNumber(flower)*2
 */
int64_t flower_getBlockEndNumber(Flower *flower);

/*
 * Returns number of stub ends in the flower.
 */
int64_t flower_getStubEndNumber(Flower *flower);

/*
 * Returns number of free stub ends in the flower.
 */
int64_t flower_getFreeStubEndNumber(Flower *flower);

/*
 * Returns the number of attached stub ends in the flower.
 */
int64_t flower_getAttachedStubEndNumber(Flower *flower);

/*
 * Gets an iterator to iterate through the ends in the flower, at this level.
 */
Flower_EndIterator *flower_getEndIterator(Flower *flower);

/*
 * Gets the next end from the iterator.
 */
End *flower_getNextEnd(Flower_EndIterator *endIterator);

/*
 * Gets the previous end from the iterator.
 */
End *flower_getPreviousEnd(Flower_EndIterator *endIterator);

/*
 * Duplicates the iterator.
 */
Flower_EndIterator *flower_copyEndIterator(Flower_EndIterator *endIterator);

/*
 * Destructs the iterator.
 */
void flower_destructEndIterator(Flower_EndIterator *endIterator);

/*
 * Gets the 'first' block.
 */
Block *flower_getFirstBlock(Flower *flower);

/*
 * Gets an block by name.
 */
Block *flower_getBlock(Flower *flower, Name name);

/*
 *  Gets the 'first' segment.
 */
Segment *flower_getFirstSegment(Flower *flower);

/*
 * Gets an segment by its name.
 */
Segment *flower_getSegment(Flower *flower, Name name);

/*
 * Returns the number of segments.
 */
int64_t flower_getSegmentNumber(Flower *flower);

/*
 * Gets an iterator to iterate through the segments in the flower, at this level.
 */
Flower_SegmentIterator *flower_getSegmentIterator(Flower *flower);

/*
 * Gets the next segment from the iterator.
 */
Segment *flower_getNextSegment(Flower_SegmentIterator *segmentIterator);

/*
 * Gets the previous segment from the iterator.
 */
Segment *flower_getPreviousSegment(Flower_SegmentIterator *segmentIterator);

/*
 * Duplicates the iterator.
 */
Flower_SegmentIterator *flower_copySegmentIterator(
        Flower_SegmentIterator *segmentIterator);

/*
 * Destructs the iterator.
 */
void flower_destructSegmentIterator(Flower_SegmentIterator *segmentIterator);

/*
 * Returns the number of blocks.
 */
int64_t flower_getBlockNumber(Flower *flower);

/*
 * Gets an iterator to iterate through the blocks in the flower, at this level.
 */
Flower_BlockIterator *flower_getBlockIterator(Flower *flower);

/*
 * Gets the next block from the iterator.
 */
Block *flower_getNextBlock(Flower_BlockIterator *blockIterator);

/*
 * Gets the previous block from the iterator.
 */
Block *flower_getPreviousBlock(Flower_BlockIterator *blockIterator);

/*
 * Duplicates the iterator
 */
Flower_BlockIterator *flower_copyBlockIterator(
        Flower_BlockIterator *blockIterator);

/*
 * Destructs the iterator.
 */
void flower_destructBlockIterator(Flower_BlockIterator *blockIterator);

/*
 * Gets the 'first' group.
 */
Group *flower_getFirstGroup(Flower *flower);

/*
 * Gets an group by the name of the nested flower it contains.
 */
Group *flower_getGroup(Flower *flower, Name flowerName);

/*
 * Returns the number of groups.
 */
int64_t flower_getGroupNumber(Flower *flower);

/*
 * Gets an iterator to iterate through the groups in the flower, at this level.
 */
Flower_GroupIterator *flower_getGroupIterator(Flower *flower);

/*
 * Gets the next group from the iterator.
 */
Group *flower_getNextGroup(Flower_GroupIterator *groupIterator);

/*
 * Gets the previous group from the iterator.
 */
Group *flower_getPreviousGroup(Flower_GroupIterator *groupIterator);

/*
 * Duplicates the iterator.
 */
Flower_GroupIterator *flower_copyGroupIterator(
        Flower_GroupIterator *groupIterator);

/*
 * Destructs the iterator.
 */
void flower_destructGroupIterator(Flower_GroupIterator *groupIterator);

/*
 * Returns non-zero if the flower has a parent group (only the root should have no parent group once the
 * structure is finished).
 */
bool flower_hasParentGroup(Flower *flower);

/*
 * Gets the parent group of the flower.
 */
Group *flower_getParentGroup(Flower *flower);

/*
 * Gets the 'first' chain.
 */
Chain *flower_getFirstChain(Flower *flower);

/*
 * Gets a chain by its index
 */
Chain *flower_getChain(Flower *flower, Name name);

/*
 * Returns the number of chains.
 */
int64_t flower_getChainNumber(Flower *flower);

/*
 * Returns the number of non-trivial chains.
 */
int64_t flower_getTrivialChainNumber(Flower *flower);

/*
 * Gets an iterator to iterate through the chains in the flower, at this level.
 */
Flower_ChainIterator *flower_getChainIterator(Flower *flower);

/*
 * Gets the next chain from the iterator.
 */
Chain *flower_getNextChain(Flower_ChainIterator *chainIterator);

/*
 * Gets the previous chain from the iterator.
 */
Chain *flower_getPreviousChain(Flower_ChainIterator *chainIterator);

/*
 * Duplicates the iterator.
 */
Flower_ChainIterator *flower_copyChainIterator(
        Flower_ChainIterator *chainIterator);

/*
 * Destructs the iterator.
 */
void flower_destructChainIterator(Flower_ChainIterator *chainIterator);

/*
 * Gets the 'first' face.
 */
Face *flower_getFirstFace(Flower *flower);

/*
 * Returns the number of faces.
 */
int64_t flower_getFaceNumber(Flower *flower);

/*
 * Gets an iterator to iterate through the faces in the flower, at this level.
 */
Flower_FaceIterator *flower_getFaceIterator(Flower *flower);

/*
 * Gets the next face from the iterator.
 */
Face *flower_getNextFace(Flower_FaceIterator *faceIterator);

/*
 * Gets the previous face from the iterator.
 */
Face *flower_getPreviousFace(Flower_FaceIterator *faceIterator);

/*
 * Duplicates the iterator.
 */
Flower_FaceIterator *flower_copyFaceIterator(Flower_FaceIterator *faceIterator);

/*
 * Destructs the iterator.
 */
void flower_destructFaceIterator(Flower_FaceIterator *faceIterator);

/*
 * Get flower size, in terms of total bases it contains. Looks only at threads have defined sequences.
 */
int64_t flower_getTotalBaseLength(Flower *flower);

/*
 * Merges together the two flowers and there parent groups.
 *
 * Only works if both parent groups do not have links. Merging together groups that
 * are in links means breaking the chains, which it currently will not do.
 */
//Flower *flower_mergeFlowers(Flower *flower1, Flower *flower2);

/*
 * Runs check function for each type of object contained in the flower.
 * Checks that flower_builtTrees and flower_builtFaces are correctly set.
 */
void flower_check(Flower *flower);

/*
 * Checks that the flower contains at least one end, unless it is the parent problem
 * and the whole reconstruction is empty. Also checks that each group contains at least one attached/block end.
 */
void flower_checkNotEmpty(Flower *flower, bool recursive);

/*
 * Runs flower_check for the given flower and all nested flowers.
 */
void flower_checkRecursive(Flower *flower);

/*
 * Returns non-zero iff the blocks for the flower have been added (i.e. no further
 * alignment will be added to the flower).
 */
bool flower_builtBlocks(Flower *flower);

/*
 * Switches the status of flower_buildBlocks(). By default flower_builtBlocks returns
 * 0.
 */
void flower_setBuiltBlocks(Flower *flower, bool b);

/*
 * Returns non-zero iff every end tree and block tree is well defined.
 */
bool flower_builtTrees(Flower *flower);

/*
 * Switches the status of flower_buildTrees(). By default flower_builtTrees returns
 * 0.
 */
void flower_setBuiltTrees(Flower *flower, bool b);

/*
 * Returns non-zero iff faces are being built for the tree. If this is set
 * then faces are built for the AVG.
 */
bool flower_builtFaces(Flower *flower);

/*
 * Switches the status of flower_builtFaces(). By default flower_builtFaces returns
 * 0.
 */
void flower_setBuildFaces(Flower *flower, bool b);

/*
 * Returns non-zero iff the flower has no nested flowers.
 */
bool flower_isLeaf(Flower *flower);

/*
 * Returns non-zero iff the flower is a leaf, contains zero or one group and contains only stub ends.
 */
bool flower_isTerminal(Flower *flower);

/*
 * Deletes the flower if it is not the root and not a leaf and contains no blocks. Returns
 * non-zero if deleted.
 */
bool flower_removeIfRedundant(Flower *flower);

/*
 * If the flower contains no ends and is not the root flower of the tree the flower is removed, along
 * with any children. Returns non-zero if the flower is removed, else returns false.
 */
bool flower_deleteIfEmpty(Flower *flower);

/*
 * Deletes the flower and its children.
 */
void flower_delete(Flower *flower);

/*
 * Deletes the flower and its children. Can optionally specify that flower is not on disk (i.e. if you know you just created it),
 * and therefore avoid adding to the database communication.
 */
void flower_delete2(Flower *flower, bool isOnDisk);

/*
 * Ensures that all terminal groups have an attached leaf flower.
 */
void flower_makeTerminalNormal(Flower *flower);

/*
 * Returns non zero if the parent is loaded.
 */
bool flower_isParentLoaded(Flower *flower);

/*
 * If the parent of the flower is memory it is unloaded. This should be used with care, as if the parent has not yet
 * be written to disk then it will not exist after the end of session.
 */
void flower_unloadParent(Flower *flower);

/*
 * Unload flowers.
 */
void flower_unload(Flower *flower);

#endif
