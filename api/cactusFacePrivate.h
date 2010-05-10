#ifndef _CACTUS_FACE_PRIVATE_H_
#define _CACTUS_FACE_PRIVATE_H_

#include "cactusGlobals.h"

struct _face {
	// Net
	Net * net;

	// Serialisation name
	Name name;

	// Number of top nodes in face, length of subsequent arrays
	int32_t cardinal;

	// Array of pointers of top nodes of the face
	Cap **topNodes;

	// Array of array of pointers of bottom nodes of the face
	// One array for each top node
	Cap ***bottomNodes;

	// Number of lifted edges which are colinear to the
	// corresponding top node's adjacency edge
	int32_t *bottomNodeNumbers;

	// Destination of the isolated lifted edge
	Cap ***derivedEdgeDestinations;

	//Pointer to memory for faceEnds..
	FaceEnd **faceEnds;
};

struct _face_FaceEndIterator {
	Face *face;
	int32_t index;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private face functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Private constructor
 */
Face * face_construct2(Name name, Net * net);

/*
 * Destructs the face
 */
void face_destruct(Face * face);

/*
 * Creates a binary representation of the face, returned as a char string.
 */
void face_writeBinaryRepresentation(Face * face,
				    void (*writeFn) (const void *ptr,
						     size_t size,
						     size_t count));

/*
 * Loads a face into memory from a binary representation of the face.
 */
Face *face_loadFromBinaryRepresentation(void **binaryString, Net * net);


/*
 * Get a static instance (from the heap) with the name set.
 */
Face *face_getStaticNameWrapper(Name name);

/*
 * Sets the net associated with the face.
 */
void face_setNet(Face *face, Net *net);

/*
 * Gets the face end associated with the top node of the cap. The public
 * way todo this is cap_getTopFaceEnd(cap);
 */
FaceEnd *face_getFaceEndForTopNode(Face *face, Cap *cap);

#endif
