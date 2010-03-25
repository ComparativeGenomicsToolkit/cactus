#ifndef BUILD_FACES_H_ 
#define BUILD_FACES_H_

#include "cactus.h"

/*
 * The works: 
 * 	Build faces from current adjacencies
 * 	Regularize faces,
 * 	Canonize faces.
 */
void buildFaces_buildAndProcessFaces(Net * net);

/*
 * Isolates into a regular and trivial faces
 */
void face_isolate(Face * face, Net * net);

/*
 * Isolates all the faces in the net
 */
void face_isolateFaces(Net * net);

/*
 * Canonizes face into a regular cycle
 */
void face_canonize(Face * face, Net * net);

/*
 * Canonizes all the faces in the net
 */
void face_canonizeFaces(Net * net);

/*
 * Tests if simple
 */
void face_isSimple(Face * face);

/*
 * Tests if regular
 */
int32_t face_isRegular(Face * face);

/*
 * Tests if canonical
 */
int32_t face_isCanonical(Face * face);
#endif
