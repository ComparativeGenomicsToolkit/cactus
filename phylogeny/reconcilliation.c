/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * reconcilliation.c
 *
 *  Created on: 31-Mar-2010
 *      Author: benedictpaten
 */
#include <assert.h>

#include "cactus.h"
#include "commonC.h"

Event *labelAndReturnEvent(struct BinaryTree *blockTree, Event *event) {
	assert(event != NULL);
	blockTree->label = cactusMisc_nameToString(event_getName(event));
	return event;
}

Event *reconcile(struct BinaryTree *blockTree, EventTree *eventTree,
		Event *(*getEventFromLeaf)(struct BinaryTree *binaryTree, void *extraArg), void *extraArg) {
	if(blockTree->internal) { //internal node
		Event *leftChild = reconcile(blockTree->left, eventTree, getEventFromLeaf, extraArg);
		Event *rightChild = reconcile(blockTree->right, eventTree, getEventFromLeaf, extraArg);
		if(leftChild != NULL) {
			if(rightChild != NULL) {
				Event *commonAncestor = eventTree_getCommonAncestor(leftChild, rightChild);
				assert(commonAncestor != eventTree_getRootEvent(eventTree));
				if(commonAncestor == leftChild || commonAncestor == rightChild) {
					Event *parentEvent = event_getParent(commonAncestor);
					if(event_getChildNumber(parentEvent) != 1 || parentEvent == eventTree_getRootEvent(eventTree)) {
						//we need to invent a new event before the parent speciation or root of tree
						Event *newEvent = event_construct4(NULL, 0.01, parentEvent, commonAncestor, eventTree);
						return labelAndReturnEvent(blockTree, newEvent);
					}
					else { //there already exists a valid unary event..
						return labelAndReturnEvent(blockTree, parentEvent);
					}
				}
				else { //speciation case
					return labelAndReturnEvent(blockTree, commonAncestor);
				}
			}
			else {
				return leftChild;
			}
		}
		else {
			return rightChild;
		}
	}
	else { //leaf case
		Event *childEvent = getEventFromLeaf(blockTree, extraArg);
		return childEvent;
	}
}
