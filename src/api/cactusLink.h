/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_LINK_H_
#define CACTUS_LINK_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic link functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Construct a link.
 */
Link *link_construct(End *_3End, End *_5End, Group *group, Chain *parentChain);

/*
 * Gets the next link in the link.
 */
Link *link_getNextLink(Link *link);

/*
 * Gets the prior link in the link.
 */
Link *link_getPreviousLink(Link *link);

/*
 * Gets the nested flower the link contains.
 */
Group *link_getGroup(Link *link);

/*
 * Gets the left end of the link in the link, which will
 * be positively oriented and a 3' end.
 */
End *link_get3End(Link *link);

/*
 * Gets the right end of the link in the link, which will
 * be positively oriented and a 5' end.
 */
End *link_get5End(Link *link);

/*
 * Gets the chain the link is part of.
 */
Chain *link_getChain(Link *link);

/*
 * Destroys the link and splits any containing chain into two chains, before and after,
 * unless they have zero length.
 */
void link_split(Link *link);

/*
 * Returns true if and only if the two ends of the link are block ends and they are always adjacent (no self loops), and
 * all the adjacencies are trivial (empty).
 */
bool link_isTrivial(Link *link);

/*
 * Removes the link from the chain and destroys the ends and any nested problem is the link is trivial.
 * Only works before trees have been built. Return non zero is the link is removed.
 */
bool link_mergeIfTrivial(Link *link);

#endif
