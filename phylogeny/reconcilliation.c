/*
 * reconcilliation.c
 *
 *  Created on: 31-Mar-2010
 *      Author: benedictpaten
 */

#include "cactus.h"
#include "commonC.h"

Event *reconcile(struct BinaryTree *blockTree, EventTree *eventTree) {
	if(blockTree->internal) { //internal node
		Event *leftChild = reconcile(blockTree->left, eventTree);
		Event *rightChild = reconcile(blockTree->right, eventTree);
		Event *commonAncestor = eventTree_getCommonAncestor(leftChild, rightChild);
		if(commonAncestor == leftChild || commonAncestor == rightChild) {
			Event *parentEvent = event_getParent(commonAncestor);
			if(event_getChildNumber(parentEvent) != 1 || parentEvent == eventTree_getRootEvent(eventTree)) {
				//we need to invent a new event before the parent speciation or root of tree
				MetaEvent *metaEvent = metaEvent_construct(NULL, net_getNetDisk(eventTree_getNet(eventTree)));
				Event *newEvent = event_construct(metaEvent, 0.01, parentEvent, eventTree);
				blockTree->label = netMisc_nameToString(event_getName(newEvent));
				return newEvent;
			}
			else { //there already exists a valid unary event..
				blockTree->label = netMisc_nameToString(event_getName(parentEvent));
				return parentEvent;
			}
		}
		else { //speciation case
			blockTree->label = netMisc_nameToString(event_getName(commonAncestor));
			return commonAncestor;
		}
	}
	else { //leaf case
		return eventTree_getEvent(eventTree, netMisc_stringToName(blockTree->label));
	}
}
