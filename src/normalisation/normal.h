/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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
 * Ensures that all chains in this flower can not be extended by chains within its nested flowers,
 * by promoting those chains in the nested flowers which extend chains in the flower.
 *
 * Let l be a link involving one attached stub-end and one block end. Let m be the
 * input flower and n be a nested flower of n. For any l in any
 * n, if the attached end in l is in a link in m, then this
 * function promotes the chain containing l into m.
 */
void promoteNestedChainsThatExtendChains(Flower *flower);

/*
 * Promotes chains from the nested flowers of the input flower into the input flower while the number of chains in the input
 * flower is less than maxNumberOfChains (this includes trivial chains including just one block).
 */
void promoteNestedChainsToFillFlower(Flower *flower, int64_t maxNumberOfChains);

/*
 * Removes all links which are trivial from the flower.
 */
void removeTrivialLinks(Flower *flower);

/*
 * Function which runs the above functions and others to ensure the children
 * of the given flower are all 'normalised'. When this function is run for each flower
 * in a cactus graph, bottom up from leaves to the root the resulting cactus tree will be normalised.
 * Several iterations are needed if we wish to push the maxNumberOfChains up into the parents of the graph.
 */
void normalise(Flower *flower, int64_t maxNumberOfChains);

#endif /* NORMAL_H_ */
