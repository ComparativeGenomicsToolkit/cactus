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
 * Constructs a non-terminal group.
 */
Group *group_construct(Net *net, Net *nestedNet);

/*
 * Constructs an group without a nested net. The group is thus
 * labelled terminal (see group_isTerminal(group)).
 */
Group *group_construct2(Net *net);

/*
 * A non-terminal group is one with out a nested net holding a subproblem and further recursion.
 * Returns 1 if the group is terminal, else returns zero.
 */
bool group_isTerminal(Group *group);

/*
 * Converts a terminal group into a non-terminal group,
 * constructing a nested net containing the appropriate ends. The leaf adjacencies will
 * be set.
 *
 * Will fail if the problem is already non-terminal.
 */
void group_makeNonTerminal(Group *group);

/*
 * Adds an end to the group.
 */
void group_addEnd(Group *group, End *end);

/*
 *  Gets the net the group is part of.
 */
Net *group_getNet(Group *group);

/*
 * Gets the name of the group. This name will also be the name of the nested net,
 * if the group has one.
 */
Name group_getName(Group *group);

/*
 * Gets the nested net the group contains, or NULL if it doesn't contain one.
 */
Net *group_getNestedNet(Group *group);

/*
 * Gets the link the group is part of, or NULL, if not part of a chain.
 */
Link *group_getLink(Group *group);

/*
 * Gets an end by name
 */
End *group_getEnd(Group *group, Name name);

/*
 * Returns the number of ends.
 */
int32_t group_getEndNumber(Group *group);

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
 * Gets the total number of bases in the group. Only works if the adjacencies have been set for the
 * leaf instances of the ends in the group.
 */
int64_t group_getTotalBaseLength(Group *group);

/*
 * Merges together the two groups and there nested nets, if they have them.
 *
 * Only works if both groups do not have links. Merging together groups that
 * are in links means breaking the chains, which it currently will not do.
 */
void group_mergeGroups(Group *group1, Group *group2);

#endif
