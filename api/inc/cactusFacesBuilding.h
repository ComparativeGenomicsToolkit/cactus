/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _CACTUS_FACESBUILDING_H_
#define _CACTUS_FACESBUILDING_H_

/*
 * Destroy then rebuild all the faces in a flower 
 */
void flower_reconstructFaces(Flower * flower);

/*
 * Constructs a face locally from a given Cap but without precomputed liftedEdges
 */
void buildFaces_reconstructFromCap(Cap * startingCap, Flower * flower);

/*
 * Destroy all the faces in a flower
 */
void flower_destructFaces(Flower *flower);
#endif
