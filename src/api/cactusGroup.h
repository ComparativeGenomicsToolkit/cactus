/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_ADJACENCY_COMPONENT_H_
#define CACTUS_ADJACENCY_COMPONENT_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic group functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a non-leaf group.
 */
Group *group_construct(Flower *flower, Flower *nestedFlower);

/*
 * Constructs an group without a nested flower. The group is thus
 * labelled a leaf (see group_isLeaf(group)).
 */
Group *group_construct2(Flower *flower);

/*
 * As group_construct2, but allowing the name to be specified.
 */
Group *group_construct3(Flower *flower, Name name);

/*
 * Destructs a group.
 */
void group_destruct(Group *group);

/*
 *  Gets the flower the group is part of.
 */
Flower *group_getFlower(Group *group);

/*
 * Gets the name of the group. This name will also be the name of the nested flower,
 * if the group has one.
 */
Name group_getName(Group *group);

/*
 * Equivalent to group_getNestedFlower(group) == NULL.
 */
bool group_isLeaf(Group *group);

/*
 * Converts a leaf group into a non-leaf group, copying the event
 * tree, but nothing else.
 */
Flower *group_makeEmptyNestedFlower(Group *group);

/*
 * Converts a leaf group into a non-leaf group,
 * constructing a nested flower containing the appropriate ends. The leaf adjacencies will
 * be set and all the ends will be in one new leaf group.
 *
 * Will fail if the problem is already a non-leaf group.
 * Returns the nested flower.
 */
Flower *group_makeNestedFlower(Group *group);

/*
 * Gets the nested flower the group contains, or NULL if it doesn't contain one.
 */
Flower *group_getNestedFlower(Group *group);

/*
 * Gets the link the group is part of, or NULL, if not part of a chain.
 */
Link *group_getLink(Group *group);

/*
 * group_isTangle(group) == (group_getLink(group) == NULL)
 */
bool group_isTangle(Group *group);

/*
 * group_isLink(group) != group_isTangle(group);
 */
bool group_isLink(Group *group);

/*
 * Gets the first end in the group.
 */
End *group_getFirstEnd(Group *group);

/*
 * Gets an end by name
 */
End *group_getEnd(Group *group, Name name);

/*
 * Returns the number of ends.
 */
int64_t group_getEndNumber(Group *group);

/*
 * Returns the number of stub ends in the group.
 */
int64_t group_getStubEndNumber(Group *group);

/*
 * Returns the number of attached stub ends.
 */
int64_t group_getAttachedStubEndNumber(Group *group);

/*
 * Returns the number of free stub ends.
 */
int64_t group_getFreeStubEndNumber(Group *group);

/*
 * Returns the number of block ends.
 */
int64_t group_getBlockEndNumber(Group *group);

/*
 * Gets an iterator to iterate through the ends in the group.
 */
Group_EndIterator *group_getEndIterator(Group *group);

/*
 * Gets the next end from the iterator.
 */
End *group_getNextEnd(Group_EndIterator *endIterator);

/*
 * Gets the previous end from the iterator.
 */
End *group_getPreviousEnd(Group_EndIterator *endIterator);

/*
 * Duplicates the iterator.
 */
Group_EndIterator *group_copyEndIterator(Group_EndIterator *endIterator);

/*
 * Destructs the iterator.
 */
void group_destructEndIterator(Group_EndIterator *endIterator);

/*
 * Gets the total number of bases in the group for threads that have defined sequences.
 */
int64_t group_getTotalBaseLength(Group *group);

/*
 * Merges together the two groups and there nested flowers, if they have them.
 *
 * Only works if both groups do not have links. Merging together groups that
 * are in links means breaking the chains, which it currently will not do.
 */
//Group *group_mergeGroups(Group *group1, Group *group2);

/*
 * Checks (amongst other things) the following:
 * That the ends of the groups are doubly linked to the ends (so every end is in only one link).
 * That, if terminal has no nested flower,
 * else that any nested flower contains the correct set of stub ends.
 * That if the group has only two non-free stub-ends:
 *  	it has a link group, with containing chain.
 * else:
 * 		it is not a link group, with no containing chain.
 */
void group_check(Group *group);

/*
 * If the group contains only two non-free stub ends (either block or attached stub ends),
 * then the function constructs a chain and puts the two ends in the link.
 */
void group_constructChainForLink(Group *group);

#endif
