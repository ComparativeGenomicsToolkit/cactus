#ifndef CACTUS_PULLDOWN_H_
#define CACTUS_PULLDOWN_H_

#include "cactus.h"
#include "cactusGlobals.h"

/*
 * Take a flower, and determine if any pulldowns are necessary. 
 * Returns the Name of either a Segment, or a Cap w/ adjacency, 
 *  that needs a pulldown if such an object exists.
 * Returns NULL_NAME if no pulldowns are necessary.
 * */
Name pulldown_getFlowerPulldownName(Flower *flower);

// Calculate the weight/badness for a given segment in some Flower.
int32_t pulldown_getSegmentWeight(Segment *segment);
int32_t pulldown_getSegmentBadness(Segment *segment);

// Calculate the weight/badness for the adjacency of a Cap in some flower.
int32_t pulldown_getCapAdjacencyWeight(Cap *cap);
int32_t pulldown_getCapAdjacencyBadness(Cap *cap);

/*
 * Take a flower, and perform pulldowns as necessary to insure junctions
 * are consistently attached.
 * */
void pulldown_doFlowerPulldowns(Flower * flower);

#endif
