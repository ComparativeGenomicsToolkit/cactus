/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * FaceEnd.h
 *
 *  Created on: 10-May-2010
 *      Author: benedictpaten
 */

#ifndef FACEEND_H_
#define FACEEND_H_

#include "cactusGlobals.h"

/*
 * Gets the top node associated with the face end.
 */
Cap *faceEnd_getTopNode(FaceEnd *faceEnd);

/*
 * Gets the face containing the face end.
 */
Face *faceEnd_getFace(FaceEnd *faceEnd);

/*
 * Gets the number of bottom nodes in the face.
 */
int64_t faceEnd_getNumberOfBottomNodes(FaceEnd *faceEnd);

/*
 * Gets an iterator over the bottom nodes of the face end.
 */
FaceEnd_BottomNodeIterator *faceEnd_getBottomNodeIterator(FaceEnd *faceEnd);

/*
 * Gets the next cap in the face end.
 */
Cap *faceEnd_getNextBottomNode(FaceEnd_BottomNodeIterator *iterator);

/*
 * Gets the previous cap in the face end.
 */
Cap *faceEnd_getPreviousBottomNode(FaceEnd_BottomNodeIterator *iterator);

/*
 * Copies the iterator.
 */
FaceEnd_BottomNodeIterator *faceEnd_copyBottomNodeIterator(FaceEnd_BottomNodeIterator *iterator);

/*
 * Destructs the face bottom node iterator.
 */
void faceEnd_destructBottomNodeIterator(FaceEnd_BottomNodeIterator *iterator);

/*
 * Checks the structure of the face end.
 */
void faceEnd_check(FaceEnd *faceEnd);

#endif /* FACEEND_H_ */
