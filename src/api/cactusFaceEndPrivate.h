/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * faceEndPrivate.h
 *
 *  Created on: 10-May-2010
 *      Author: benedictpaten
 */

#ifndef FACEENDPRIVATE_H_
#define FACEENDPRIVATE_H_

#include "cactusGlobals.h"

struct _faceEnd {
	Face *face;
	int64_t topNodeIndex;
};

struct _faceEndIterator {
	FaceEnd *faceEnd;
	int64_t bottomNodeIndex;
};

/*
 * Creates the face end object, currently it just sits on top of the face object.
 */
FaceEnd *faceEnd_construct(Face *face, int64_t topNodeIndex);

/*
 * Destructs the memory for the face end object.
 */
void faceEnd_destruct(FaceEnd *faceEnd);


#endif /* FACEENDPRIVATE_H_ */
