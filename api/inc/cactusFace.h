/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _CACTUS_FACE_H_
#define _CACTUS_FACE_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic face functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Destructor
 */
void face_destruct(Face * face);

/*
 * Returns the flower that contains the face.
 */
Flower *face_getFlower(Face *face);

/*
 * Returns cardinal
 */
int32_t face_getCardinal(Face * face);

/*
 * Get selected top node
 */
Cap * face_getTopNode(Face * face, int32_t index);

/*
 * Gets an iterator over the face ends in the face.
 */
Face_FaceEndIterator *face_getFaceEndIterator(Face *face);

/*
 * Gets next face end.
 */
FaceEnd *face_getNextFaceEnd(Face_FaceEndIterator *iterator);

/*
 * Gets the previous face end.
 */
FaceEnd *face_getPreviousFaceEnd(Face_FaceEndIterator *iterator);

/*
 * Duplicates the iterator.
 */
Face_FaceEndIterator *face_copyFaceEndIterator(Face_FaceEndIterator *iterator);

/*
 * Deletes the iterator.
 */
void face_destructFaceEndIterator(Face_FaceEndIterator *iterator);

/*
 * Tests if simple
 */
int32_t face_isSimple(Face * face);

/*
 * Tests if regular
 */
int32_t face_isRegular(Face * face);

/*
 * Tests if canonical
 */
int32_t face_isCanonical(Face * face);

/*
 * Checks the face has the properties that we expect. Creates an assertion error if not.
 */
void face_check(Face *face);

/*
 * Checks that the set of faces is as we expect - with a face object created
 * for each valid face structure in the DNA history graph.
 */
void face_checkFaces(Flower *flower);

#endif
