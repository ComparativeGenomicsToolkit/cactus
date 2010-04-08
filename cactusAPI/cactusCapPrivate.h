#ifndef CACTUS_END_INSTANCE_PRIVATE_H_
#define CACTUS_END_INSTANCE_PRIVATE_H_

#include "cactusGlobals.h"

typedef struct _capContents {
	Name instance;
	int32_t coordinate;
	bool strand;
	bool side;
	Event *event;
	Sequence *sequence;
	Cap *adjacency;
	Cap *adjacency2;
	Face *face;
	Segment *segment;
	Cap *parent;
	struct List *children;
} CapContents;

struct _cap {
	CapContents *capContents;
	End *end;
	Cap *rCap;
};


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//End instance functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an cap, but not its connecting objects. Instance is the suffix m of the instance name n.m.
 */
Cap *cap_construct3(Name name, Event *event, End *end);

/*
 * As default constructor, but also sets the instance's coordinates and event.
 */
Cap *cap_construct4(Name name, End *end, int32_t startCoordinate, int32_t strand, int32_t side, Sequence *sequence);

/*
 * Destructs the cap, but not any connecting objects.
 */
void cap_destruct(Cap *cap);

/*
 * Gets the segment associated with the end, or NULL, if the end has no associated block end at this level.
 * The segment returned will have the cap on its left side.
 */
void cap_setSegment(Cap *cap, Segment *segment);

/*
 * Sets the face associated with the cap.
 */
void cap_setFace(Cap *cap, Face *face);

/*
 * Sets any adjacent ends for the alternative adjacency of the instance to NULL;
 */
void cap_breakAdjacency2(Cap *cap);

/*
 * Write a binary representation of the cap to the write function.
 */
void cap_writeBinaryRepresentation(Cap *cap, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a net into memory from a binary representation of the net.
 */
Cap *cap_loadFromBinaryRepresentation(void **binaryString, End *end);

/*
 * Get a static instance (from the heap) with the name set.
 */
Cap *cap_getStaticNameWrapper(Name name);

/*
 * Sets the event associated with the cap. Dangerous method, must be used carefully.
 */
void cap_setEvent(Cap *cap, Event *event);

/*
 * Sets the sequence associated with the cap. Dangerous method, must be used carefully.
 */
void cap_setSequence(Cap *cap, Sequence *sequence);


#endif
