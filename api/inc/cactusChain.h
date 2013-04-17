/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_CHAIN_H_
#define CACTUS_CHAIN_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic chain functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a chain, which in turn holds links.
 */
Chain *chain_construct(Flower *flower);

/*
 * Destructs the chain. Does not mess with groups, should be clean.
 */
void chain_destruct(Chain *chain);

/*
 * Gets the first link in the chain.
 */
Link *chain_getFirst(Chain *chain);

/*
 * Gets the last link in the chain.
 */
Link *chain_getLast(Chain *chain);

/*
 * Returns the number of links in the chain.
 */
int64_t chain_getLength(Chain *chain);

/*
 * Get a list of the blocks, in order, in the chain. Block number is initialised with the
 * length of the array. The array must be freed.
 */
Block **chain_getBlockChain(Chain *chain, int64_t *blockNumber);

/*
 * Gets the name of the chain in the flower.
 */
Name chain_getName(Chain *chain);

/*
 * Gets the parent flower of the chain.
 */
Flower *chain_getFlower(Chain *chain);

/*
 * Calculates the average number of bases in an instance of the chain (including those bases in the links)
 */
double chain_getAverageInstanceBaseLength(Chain *chain);

/*
 * Returns non-zero if the chain is circular.
 */
bool chain_isCircular(Chain *chain);

/*
 * Checks (amongst other things) the following:
 * That each link is properly contained in the chain.
 * Links and the contained ends are properly connected.
 * That each contiguous pair of link groups are bridged by a block.
 * If a block end is at the 5 or 3 prime end of a chain the other end of the
 * block is not in a link group (otherwise the chain is not maximal).
 * That stub ends are not the ends of the links in the chain.
 */
void chain_check(Chain *chain);

/*
 * Pushes the chain into the higher level flower. Will not work if the chain to be promoted is nested
 * within a chain of the higher level flower.
 */
void chain_promote(Chain *chain);

#endif
