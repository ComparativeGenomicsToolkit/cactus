/*
 * phylogeny.h
 *
 *  Created on: 31 Mar 2010
 *      Author: benedictpaten
 */

#ifndef PHYLOGENY_H_
#define PHYLOGENY_H_

/*
 * This takes a binary block tree with leaves labelled with valid leaf events
 * from the given event tree. The internal nodes are then greedily reconciled to
 * the most recent possible event in the event tree.
 *
 * The return value is the event which the root binary node is reconciled to.
 */
Event *reconcile(struct BinaryTree *blockTree, EventTree *eventTree);


#endif /* PHYLOGENY_H_ */
