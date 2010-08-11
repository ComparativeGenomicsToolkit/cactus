/*
 * balanceTangles.h
 *
 *  Created on: 10 Aug 2010
 *      Author: benedictpaten
 */

#ifndef BALANCETANGLES_H_
#define BALANCETANGLES_H_

#include "sonLib.h"
#include "pinchGraph.h"
#include "cactusGraph.h"
#include "cactus.h"

/*
 * Ensures there are an even number of attached stub ends in each group of the flower.
 * Assumes that the flower contains an even number of attached stub ends and that
 * each group contains one or more attached stub ends or block ends.
 * Call a tangle with an odd number of non-free stub ends an 'odd tangle'.
 * By the above assumptions, there must be an even number N of odd tangles.
 * The algorithm does the following:
 * (1) Create N/2 empty blocks.
 * (2) For each odd tangle:
 *    (a) add one of the empty block ends (checking the side of the block end added)
 *    (b) if the tangle is now a link create a chain (possibly extending another chain)
 *    (c) if the tangle is not a leaf copy the added end into the nested flower
 *    and randomly assign it to an odd tangle in that flower (there must be one or more,
 *    as there were previously an odd number of ends). Repeat this process recursively until
 *    you hit a leaf tangle.
 * (3) Calls the procedure recursively for each nested flower.
 */
void balanceTangles(Flower *flower);


#endif /* BALANCETANGLES_H_ */
