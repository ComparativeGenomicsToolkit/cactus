/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef BUILD_FACES_H_
#define BUILD_FACES_H_

#include "cactus.h"

/*
 * The works:
 * 	Build faces from current adjacencies
 *	Simplify faces,
 * 	Regularize faces,
 * 	Canonize faces.
 */
void buildFaces_buildAndProcessFaces(Flower * flower);

/*
 * Simplifies a given face
 */
void buildFaces_simplify(Face * face, Flower * flower);

/* 
 * Simplifies all the faces in the flower
 */
void buildFaces_simplifyFaces(Flower * flower);

/*
 * Isolates into a regular and trivial faces
 */
void face_isolate(Face * face, Flower * flower);

/*
 * Isolates all the faces in the flower
 */
void face_isolateFaces(Flower * flower);

/*
 * Canonizes face into a regular cycle
 */
void face_canonize(Face * face, Flower * flower);

/*
 * Canonizes all the faces in the flower
 */
void face_canonizeFaces(Flower * flower);

/////
//Misc functions
////

/*
 * Creates a new free stub end in which to place caps.
 */
End *createNewFreeStubEnd(Flower *flower);

#endif
