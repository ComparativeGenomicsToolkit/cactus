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
 * Returns the net that contains the face.
 */
Net *face_getNet(Face *face);

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
 * Returns face name
 */
Name face_getName(Face * face);

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
 * Checks that the set of faces is as we expect - with a face created
 * for each non-trivial face.
 */
void face_checkFaces(Net *net);

#endif
