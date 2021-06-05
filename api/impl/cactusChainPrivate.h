/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_CHAIN_PRIVATE_H_
#define CACTUS_CHAIN_PRIVATE_H_

#include "cactusGlobals.h"

struct _chain {
	Name name;
	Flower *flower;
	Link *link;
	Link *endLink; // Could be removed if we build back to front
	//int64_t linkNumber; // can easily be removed
	//int64_t chainIndex;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Chain functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a chain, which in turn holds links.
 */
Chain *chain_construct2(Name name, Flower *flower);

/*
 * Add the link to the chain.
 */
void chain_addLink(Chain *chain, Link *childLink);

/*
 * Sets the flower containing the chain.
 */
void chain_setFlower(Chain *chain, Flower *flower);

/*
 * Joins two chains together where the _5Chain abuts at the 3' end with the _3Chain.
 */
void chain_join(Chain *_5Chain, Chain *_3Chain);

#endif
