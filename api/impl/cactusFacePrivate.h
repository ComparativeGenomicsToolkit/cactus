/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _CACTUS_FACE_PRIVATE_H_
#define _CACTUS_FACE_PRIVATE_H_

#include "cactusGlobals.h"

struct _face {
	// Flower
	Flower * flower;

	// Number of top nodes in face, length of subsequent arrays
	int64_t cardinal;

	// Array of pointers of top nodes of the face
	Cap **topNodes;

	// Array of array of pointers of bottom nodes of the face
	// One array for each top node
	Cap ***bottomNodes;

	// Number of lifted edges which are colinear to the
	// corresponding top node's adjacency edge
	int64_t *bottomNodeNumbers;

	// Destination of the isolated lifted edge
	Cap ***derivedEdgeDestinations;

	//Pointer to memory for faceEnds..
	FaceEnd **faceEnds;
};

struct _face_FaceEndIterator {
	Face *face;
	int64_t index;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private face functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs faces in flower
 */
Face * face_construct(Flower * flower);

/*
 * Get selected derived destinations for selected top node
 */
Cap * face_getDerivedDestinationAtIndex(Face * face, int64_t topIndex, int64_t derivedIndex);

/*
 * Get non-null derived destination of selected top node (useful for simple faces)
 */
Cap * face_getDerivedDestination(Face * face, int64_t index);

/*
 * Get the number of bottom nodes for the selected top node
 */
int64_t face_getBottomNodeNumber(Face * face, int64_t topIndex);

/*
 * Get selected bottom node from selected top node in face
 */
Cap * face_getBottomNode(Face * face, int64_t topNodeIndex, int64_t bottomNodeIndex);

/*
 * Allocate arrays to allow for data
 */
void face_allocateSpace(Face * face, int64_t cardinal);

/*
 * Sets the selected top node
 */
void face_setTopNode(Face * face, int64_t topIndex, Cap * topNode);

/*
 * Set bottom node count and allocate space to store pointers
 */
void face_setBottomNodeNumber(Face * face, int64_t topIndex, int64_t number);

/*
 * Sets the derived edge destination for a given top node in face
 */
void face_setDerivedDestination(Face * face, int64_t topIndex, int64_t bottomIndex, Cap * destination);

/*
 * Adds bottom node to selected top node in face
 */
void face_addBottomNode(Face * face, int64_t topIndex, Cap * bottomNode);

/*
 * Adds a top node with a single bottom node
 */
void face_engineerArtificialNodes(Face * face, Cap * topNode, Cap * bottomNode, int64_t nonDerived);

/*
 * Gets the face end associated with the top node of the cap. The public
 * way to do this is cap_getTopFaceEnd(cap);
 */
FaceEnd *face_getFaceEndForTopNode(Face *face, Cap *cap);

#endif
