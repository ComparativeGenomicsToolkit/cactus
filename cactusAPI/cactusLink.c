#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic link functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Link *link_construct(End *leftEnd, End *rightEnd, Group *group, Chain *parentChain) {
	Link *link;
	link = malloc(sizeof(Link));

	leftEnd = end_getPositiveOrientation(leftEnd);
	rightEnd = end_getPositiveOrientation(rightEnd);
	assert(leftEnd != rightEnd);

	link->leftEnd = leftEnd;
	link->rightEnd = rightEnd;
	link->chain = parentChain;
	link->group = group;

	//Checks.
	assert(group_getEnd(group, end_getName(leftEnd)) == leftEnd);
	assert(group_getEnd(group, end_getName(rightEnd)) == rightEnd);

	chain_addLink(parentChain, link); //will set the link indices.
	group_setLink(group, link);
	return link;
}

Link *link_getNextLink(Link *link) {
	return link->nLink;
}

Link *link_getPreviousLink(Link *link) {
	return link->pLink;
}

Group *link_getGroup(Link *link) {
	return link->group;
}

End *link_getLeft(Link *link) {
	return link->leftEnd;
}

End *link_getRight(Link *link) {
	return link->rightEnd;
}

Chain *link_getChain(Link *link) {
	return link->chain;
}

int32_t link_getIndex(Link *link) {
	return link->linkIndex;
}

/*
 * Private functions.
 */

void link_destruct(Link *link) {
	Link *link2;
	link_getChain(link)->linkNumber = link->linkIndex;
	if(link->pLink == NULL) {
		link_getChain(link)->link = NULL;
	}
	else {
		link->pLink->nLink = NULL;
	}
	while(link != NULL) {
		group_setLink(link_getGroup(link), NULL);
		link2 = link;
		link = link->nLink;
		free(link2);
	}
}

void link_split(Link *link) {
	Chain *chain = link_getChain(link);
	Chain *chain1 = chain_construct(chain_getNet(chain));
	Chain *chain2 = chain_construct(chain_getNet(chain));
	int32_t i = 0;
	while(i<chain_getLength(chain)) {
		Link *link2 = chain_getLink(chain, i++);
		if(link2 == link) {
			break;
		}
		link_construct(link_getLeft(link2), link_getRight(link2), end_getGroup(link_getLeft(link2)), chain1);
	}
	while(i<chain_getLength(chain)) {
		Link *link2 = chain_getLink(chain, i++);
		assert(link2 != link);
		link_construct(link_getLeft(link2), link_getRight(link2), end_getGroup(link_getLeft(link2)), chain2);
	}
	assert(chain_getLength(chain1) + chain_getLength(chain2) == chain_getLength(chain)-1);
	chain_destruct(chain);
	if(chain_getLength(chain1) == 0) {
		chain_destruct(chain1);
	}
	if(chain_getLength(chain2) == 0) {
		chain_destruct(chain2);
	}
}

/*
 * Serialisation functions.
 */

void link_writeBinaryRepresentation(Link *link, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_LINK, writeFn);
	binaryRepresentation_writeName(group_getName(link_getGroup(link)), writeFn);
	binaryRepresentation_writeName(end_getName(link_getLeft(link)), writeFn);
	binaryRepresentation_writeName(end_getName(link_getRight(link)), writeFn);
}

Link *link_loadFromBinaryRepresentation(void **binaryString, Chain *chain) {
	Group *group;
	End *leftEnd;
	End *rightEnd;
	Link *link;

	link = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_LINK) {
		binaryRepresentation_popNextElementType(binaryString);
		group = net_getGroup(chain_getNet(chain),
				binaryRepresentation_getName(binaryString));
		leftEnd = net_getEnd(chain_getNet(chain),
				binaryRepresentation_getName(binaryString));
		rightEnd = net_getEnd(chain_getNet(chain),
						binaryRepresentation_getName(binaryString));
		link = link_construct(leftEnd, rightEnd, group, chain);
	}
	return link;
}
