/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"
/*
 * FaceEnd.c
 *
 *  Created on: 10-May-2010
 *      Author: benedictpaten
 */

Cap *faceEnd_getTopNode(FaceEnd *faceEnd) {
	return face_getTopNode(faceEnd_getFace(faceEnd), faceEnd->topNodeIndex);
}

Face *faceEnd_getFace(FaceEnd *faceEnd) {
	return faceEnd->face;
}

int64_t faceEnd_getNumberOfBottomNodes(FaceEnd *faceEnd) {
	return face_getBottomNodeNumber(faceEnd_getFace(faceEnd), faceEnd->topNodeIndex);
}

FaceEnd_BottomNodeIterator *faceEnd_getBottomNodeIterator(FaceEnd *faceEnd) {
	FaceEnd_BottomNodeIterator *iterator = st_malloc(sizeof(FaceEnd_BottomNodeIterator));
	iterator->faceEnd = faceEnd;
	iterator->bottomNodeIndex = 0;
	return iterator;
}

Cap *faceEnd_getNextBottomNode(FaceEnd_BottomNodeIterator *iterator) {
	Face *face = faceEnd_getFace(iterator->faceEnd);
	if(iterator->bottomNodeIndex < face_getBottomNodeNumber(face, iterator->faceEnd->topNodeIndex)) {
		return face_getBottomNode(face, iterator->faceEnd->topNodeIndex, iterator->bottomNodeIndex++);
	}
	return NULL;
}

Cap *faceEnd_getPreviousBottomNode(FaceEnd_BottomNodeIterator *iterator) {
	Face *face = faceEnd_getFace(iterator->faceEnd);
	assert(iterator->bottomNodeIndex <= face_getBottomNodeNumber(face, iterator->faceEnd->topNodeIndex));
	if(iterator->bottomNodeIndex-1 >= 0) {
		return face_getBottomNode(face, iterator->faceEnd->topNodeIndex, --iterator->bottomNodeIndex);
	}
	return NULL;
}

FaceEnd_BottomNodeIterator *faceEnd_copyBottomNodeIterator(FaceEnd_BottomNodeIterator *iterator) {
	FaceEnd_BottomNodeIterator *iterator2 = faceEnd_getBottomNodeIterator(iterator->faceEnd);
	iterator2->bottomNodeIndex = iterator->bottomNodeIndex;
	return iterator2;
}

void faceEnd_destructBottomNodeIterator(FaceEnd_BottomNodeIterator *iterator) {
	free(iterator);
}

void faceEnd_check(FaceEnd *faceEnd) {
#ifndef NDEBUG
	Cap *topCap = faceEnd_getTopNode(faceEnd);
	assert(cap_getTopFaceEnd(topCap) == faceEnd);
	Cap *bottomCap;
	FaceEnd_BottomNodeIterator *iterator = faceEnd_getBottomNodeIterator(faceEnd);
	while((bottomCap = faceEnd_getNextBottomNode(iterator)) != NULL) {
		assert(cap_getBottomFaceEnd(bottomCap) == faceEnd);
	}
	faceEnd_destructBottomNodeIterator(iterator);
#endif
}

/*
 * Private functions
 */

FaceEnd *faceEnd_construct(Face *face, int64_t topNodeIndex) {
	FaceEnd *faceEnd = st_malloc(sizeof(FaceEnd));
	assert(topNodeIndex >= 0);
	assert(topNodeIndex <= face_getCardinal(face));
	faceEnd->face = face;
	faceEnd->topNodeIndex = topNodeIndex;
	return faceEnd;
}

void faceEnd_destruct(FaceEnd *faceEnd) {
	free(faceEnd);
}
