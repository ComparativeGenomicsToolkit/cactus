#ifndef CACTUS_END_PRIVATE_H_
#define CACTUS_END_PRIVATE_H_

#include "cactusGlobals.h"

typedef struct _endContents {
	Cap *rootInstance;
	bool isStub;
	bool isAttached;
	Name name;
	Block *attachedBlock;
	struct avl_table *caps;
	Group *group;
	Net *net;
} EndContents;

struct _end_instanceIterator {
	struct avl_traverser *iterator;
	End *end;
};

struct _end {
	EndContents *endContents;
	bool orientation;
	End *rEnd;
};


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//End functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs the end, but not any attached block.
 */
End *end_construct2(Name name, int32_t isStub, int32_t isAttached, Net *net);

/*
 * Destructs the end and any contained caps.
 */
void end_destruct(End *end);

/*
 * Sets the attached block.
 */
void end_setBlock(End *end, Block *block);

/*
 * Adds the cap to the end.
 */
void end_addInstance(End *end, Cap *cap);

/*
 * Removes the instance from the end.
 */
void end_removeInstance(End *end, Cap *cap);

/*
 * Sets the group that the end is part of.
 */
void end_setGroup(End *end, Group *group);

/*
 * Write a binary representation of the end to the write function.
 */
void end_writeBinaryRepresentation(End *end, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
End *end_loadFromBinaryRepresentation(void **binaryString, Net *net);

/*
 * Get a static instance (from the heap) with the name set.
 */
End *end_getStaticNameWrapper(Name name);

#endif
