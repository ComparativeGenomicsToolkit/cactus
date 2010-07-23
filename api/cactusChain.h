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
 * Ensures that all chains in this net are not part of a higher level chain, by promoting
 * those chains which are.
 *
 * Let l be a link involving one attached stub-end and one block end. For any l in a
 * net n, if the attached end in l is in a link in the parent of n (if n has a parent), then this
 * function promotes the chain containing l to the higher level.
 */
void chain_promoteChainsThatExtendHigherLevelChains(Net *net);

#endif
