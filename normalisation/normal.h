/*
 * normal.h
 *
 *  Created on: 17 Jul 2010
 *      Author: benedictpaten
 */

#ifndef NORMAL_H_
#define NORMAL_H_

#include "cactus.h"

/*
 * Ensures that all terminal groups have an attached leaf flower.
 */
void makeTerminalNormal(Flower *flower);

/*
 * Ensures that all chains in this flower are not part of a higher level chain, by promoting
 * those chains which are.
 *
 * Let l be a link involving one attached stub-end and one block end. For any l in a
 * flower n, if the attached end in l is in a link in the parent of n (if n has a parent), then this
 * function promotes the chain containing l to the higher level.
 */
void promoteChainsThatExtendHigherLevelChains(Flower *flower);

/*
 * Promotes chains from the flower into any parent flower while the number of chains in the parent
 * flower is less than maxNumberOfChains.
 */
void promoteChainsToFillParents(Flower *flower, int32_t maxNumberOfChains);

#endif /* NORMAL_H_ */
