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
Chain *chain_construct(Net *net);

/*
 * Gets a link in the chain.
 */
Link *chain_getLink(Chain *chain, int32_t linkIndex);

/*
 * Returns the number of links in the chain.
 */
int32_t chain_getLength(Chain *chain);

/*
 * Get a list of the blocks, in order, in the chain. Block number is initialised with the
 * length of the array. The array must be freed.
 */
Block **chain_getBlockChain(Chain *chain, int32_t *blockNumber);

/*
 * Gets the name of the chain in the net.
 */
Name chain_getName(Chain *chain);

/*
 * Gets the parent net of the chain.
 */
Net *chain_getNet(Chain *chain);

/*
 * Calculates the average number of bases in an instance of the chain (including those bases in the links)
 */
double chain_getAverageInstanceBaseLength(Chain *chain);

#endif
