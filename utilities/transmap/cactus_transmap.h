#ifndef TRANSMAP_H_ 
#define TRANSMAP_H_ 

#include "cactus.h"

/* 
 * Checks the breaks in a sequence of blocks from the present back to a 
 * specified event. 
 * Returns an array of booleans which is 1 element shorter than the array
 * of blocks (fence post principle)
 */
int32_t * transmap_syntenyWasConservedSinceEvent(Block ** blocks, int32_t blockCount, Event * event);

/* 
 * Checks the breaks in a sequence of blocks at a given time point
 * Returns an array of booleans which is 1 element shorter than the array
 * of blocks (fence post principle)
 */
int32_t * transmap_syntenyWasConservedAtEvent(Block ** blocks, int32_t blockCount, Event * event);

#endif
