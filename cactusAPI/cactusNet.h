#ifndef CACTUS_NET_H_
#define CACTUS_NET_H_

#include "cactusGlobals.h"

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
Net *net_construct(NetDisk *netDisk);

/*
 * Gets the name of the net.
 */
Name net_getName(Net *net);

/*
 * Gets the parent net disk.
 */
NetDisk *net_getNetDisk(Net *net);

/*
 * Gets the net tree associated with the event tree.
 */
EventTree *net_getEventTree(Net *net);

/*
 * Gets the 'first' sequence.
 */
Sequence *net_getFirstSequence(Net *net);

/*
 * Gets an sequence by its name.
 */
Sequence *net_getSequence(Net *net, Name name);

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
 *  Gets the 'first' cap.
 */
Cap *net_getFirstCap(Net *net);

/*
 * Gets an cap by its complete name.
 */
Cap *net_getCap(Net *net, Name name);

/*
 * Returns the number of caps.
 */
int32_t net_getCapNumber(Net *net);

/*
 * Gets an iterator to iterate through the caps in the net, at this level.
 */
Net_CapIterator *net_getCapIterator(Net *net);

/*
 * Gets the next cap from the iterator.
 */
Cap *net_getNextCap(Net_CapIterator *capIterator);

/*
 * Gets the previous cap from the iterator.
 */
Cap *net_getPreviousCap(Net_CapIterator *capIterator);

/*
 * Duplicates the iterator.
 */
Net_CapIterator *net_copyCapIterator(Net_CapIterator *capIterator);

/*
 * Destructs the iterator.
 */
void net_destructCapIterator(Net_CapIterator *capIterator);

/*
 * Gets the 'first' end.
 */
End *net_getFirstEnd(Net *net);

/*
 * Gets an end by name.
 */
End *net_getEnd(Net *net, Name name);

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
 * Gets the 'first' block.
 */
Block *net_getFirstBlock(Net *net);

/*
 * Gets an block by name.
 */
Block *net_getBlock(Net *net, Name name);

/*
 * Gets an segment by its complete name.
 */
Segment *net_getSegment(Net *net, Name completeName);

/*
 * Returns the number of blocks.
 */
int32_t net_getBlockNumber(Net *net);

/*
 * Gets an iterator to iterate through the blocks in the net, at this level.
 */
Net_BlockIterator *net_getBlockIterator(Net *net);

/*
 * Gets the next block from the iterator.
 */
Block *net_getNextBlock(Net_BlockIterator *blockIterator);

/*
 * Gets the previous block from the iterator.
 */
Block *net_getPreviousBlock(Net_BlockIterator *blockIterator);

/*
 * Duplicates the iterator
 */
Net_BlockIterator *net_copyBlockIterator(Net_BlockIterator *blockIterator);

/*
 * Destructs the iterator.
 */
void net_destructBlockIterator(Net_BlockIterator *blockIterator);

/*
 * Gets the 'first' group.
 */
Group *net_getFirstGroup(Net *net);

/*
 * Gets an group by the name of the nested net it contains.
 */
Group *net_getGroup(Net *net, Name netName);

/*
 * Returns the number of groups.
 */
int32_t net_getGroupNumber(Net *net);

/*
 * Gets an iterator to iterate through the groups in the net, at this level.
 */
Net_GroupIterator *net_getGroupIterator(Net *net);

/*
 * Gets the next group from the iterator.
 */
Group *net_getNextGroup(Net_GroupIterator *groupIterator);

/*
 * Gets the previous group from the iterator.
 */
Group *net_getPreviousGroup(Net_GroupIterator *groupIterator);

/*
 * Duplicates the iterator.
 */
Net_GroupIterator *net_copyGroupIterator(Net_GroupIterator *groupIterator);

/*
 * Destructs the iterator.
 */
void net_destructGroupIterator(Net_GroupIterator *groupIterator);

/*
 * Gets the parent group of the net.
 */
Group *net_getParentGroup(Net *net);

/*
 * Gets the 'first' chain.
 */
Chain *net_getFirstChain(Net *net);

/*
 * Gets a chain by its index
 */
Chain *net_getChain(Net *net, Name name);

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
 * Gets the 'first' face.
 */
Face *net_getFirstFace(Net *net);

/*
 * Gets an chain by index.
 */
Face *net_getFace(Net *net, Name name);

/*
 * Returns the number of faces.
 */
int32_t net_getFaceNumber(Net *net);

/*
 * Gets an iterator to iterate through the faces in the net, at this level.
 */
Net_FaceIterator *net_getFaceIterator(Net *net);

/*
 * Gets the next face from the iterator.
 */
Face *net_getNextFace(Net_FaceIterator *faceIterator);

/*
 * Gets the previous face from the iterator.
 */
Face *net_getPreviousFace(Net_FaceIterator *faceIterator);

/*
 * Duplicates the iterator.
 */
Net_FaceIterator *net_copyFaceIterator(Net_FaceIterator *faceIterator);

/*
 * Destructs the iterator.
 */
void net_destructFaceIterator(Net_FaceIterator *faceIterator);

/*
 * Get net size, in terms of total bases it contains.
 */
int64_t net_getTotalBaseLength(Net *net);

/*
 * Gets the 'first' reference.
 */
Reference *net_getFirstReference(Net *net);

/*
 * Gets an chain by index.
 */
Reference *net_getReference(Net *net, Name name);

/*
 * Returns the number of references.
 */
int32_t net_getReferenceNumber(Net *net);

/*
 * Gets an iterator to iterate through the references in the net, at this level.
 */
Net_ReferenceIterator *net_getReferenceIterator(Net *net);

/*
 * Gets the next reference from the iterator.
 */
Reference *net_getNextReference(Net_ReferenceIterator *referenceIterator);

/*
 * Gets the previous reference from the iterator.
 */
Reference *net_getPreviousReference(Net_ReferenceIterator *referenceIterator);

/*
 * Duplicates the iterator.
 */
Net_ReferenceIterator *net_copyReferenceIterator(Net_ReferenceIterator *referenceIterator);

/*
 * Destructs the iterator.
 */
void net_destructReferenceIterator(Net_ReferenceIterator *referenceIterator);


#endif
