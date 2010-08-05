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
 * Gets the index of the link in the chain..
 */
int32_t link_getIndex(Link *link);

/*
 * Destroys the link and splits any containing chain into two chains, before and after,
 * unless they have zero length.
 */
void link_split(Link *link);

#endif
