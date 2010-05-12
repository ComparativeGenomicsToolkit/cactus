#ifndef CACTUS_NET_PRIVATE_H_
#define CACTUS_NET_PRIVATE_H_

#include "cactusGlobals.h"

struct _net {
	Name name;
	EventTree *eventTree;
	Reference *reference;
	struct avl_table *sequences;
	struct avl_table *ends;
	struct avl_table *caps;
	struct avl_table *blocks;
	struct avl_table *segments;
	struct avl_table *groups;
	struct avl_table *chains;
	struct avl_table *faces;
	Name parentNetName;
	NetDisk *netDisk;
	int32_t faceIndex;
	int32_t chainIndex;
	bool builtBlocks;
	bool builtTrees;
	bool builtFaces;
};


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private net functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs the net.
 */
Net *net_construct2(Name name, NetDisk *netDisk);

/*
 * Destructs the net, and all the elements it contains. If recursive the function will destroy all
 * loaded nested nets.
 */
void net_destruct(Net *net, int32_t recursive);

/*
 * Adds the event tree for the net to the net.
 * If an previous event tree exists for the net
 * it will call eventTree_destruct for the existing tree
 * (which should not exist without the net).
 */
void net_setEventTree(Net *net, EventTree *eventTree);

/*
 * This function is called by eventTree_destruct and cleans up the reference.
 */
void net_removeEventTree(Net *net, EventTree *eventTree);

/*
 * Adds the sequence to the net.
 */
void net_addSequence(Net *net, Sequence *sequence);

/*
 * Removes the sequence from the net.
 */
void net_removeSequence(Net *net, Sequence *sequence);

/*
 * Adds the segment to the net.
 */
void net_addSegment(Net *net, Segment *segment);

/*
 * Remove the segment from the net.
 */
void net_removeSegment(Net *net, Segment *segment);

/*
 * Adds the block to the net.
 */
void net_addBlock(Net *net, Block *block);

/*
 * Remove the block from the net.
 */
void net_removeBlock(Net *net, Block *block);

/*
 * Adds the cap to the net.
 */
void net_addCap(Net *net, Cap *cap);

/*
 * Remove the cap from the net.
 */
void net_removeCap(Net *net, Cap *cap);

/*
 * Adds the end to the net.
 */
void net_addEnd(Net *net, End *end);

/*
 * Remove the end from the net.
 */
void net_removeEnd(Net *net, End *end);

/*
 * Adds the group to the net.
 */
void net_addGroup(Net *net, Group *group);

/*
 * Removes an empty group from the net.
 */
void net_removeGroup(Net *net, Group *group);

/*
 * Sets the parent group of the net.
 */
void net_setParentGroup(Net *net, Group *group);

/*
 * Adds the chain to the net.
 */
void net_addChain(Net *net, Chain *chain);

/*
 * Remove the chain from the net.
 */
void net_removeChain(Net *net, Chain *chain);

/*
 * Adds the face to the net.
 */
void net_addFace(Net *net, Face *face);

/*
 * Remove the face from the net.
 */
void net_removeFace(Net *net, Face *face);

/*
 * Adds the event tree for the net to the net.
 * If an previous event tree exists for the net
 * it will call Reference_destruct for the existing tree
 * (which should not exist without the net).
 */
void net_setReference(Net *net, Reference *reference);

/*
 * This function is called by Reference_destruct and cleans up the reference.
 */
void net_removeReference(Net *net, Reference *reference);

/*
 * This function constructs faces for the net. If faces are already created then
 * they will be first delete.
 */
void net_reconstructFaces(Net * net);

/*
 * Destroys all faces in the net.
 */
void net_destructFaces(Net *net);

/*
 * Write a binary representation of the net to the write function.
 */
void net_writeBinaryRepresentation(Net *net, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Net *net_loadFromBinaryRepresentation(void **binaryString, NetDisk *netDisk);

#endif
