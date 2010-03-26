#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic chain functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Chain *chain_construct(Net *net) {
	return chain_construct2(netDisk_getUniqueID(net_getNetDisk(net)), net);
}

Chain *chain_construct2(Name name, Net *net) {
	Chain *chain;
	chain = malloc(sizeof(Chain));
	chain->name = name;
	chain->net = net;
	chain->link = NULL;
	chain->linkNumber = 0;
	net_addChain(net, chain);
	return chain;
}

void chain_destruct(Chain *chain) {
	net_removeChain(chain_getNet(chain), chain);
	if(chain->link != NULL) {
		link_destruct(chain->link);
	}
	free(chain);
}

Link *chain_getLink(Chain *chain, int32_t linkIndex) {
	int32_t i;
	Link *link;

#ifdef BEN_DEBUG
	assert(linkIndex >= 0);
	assert(linkIndex < chain->linkNumber);
#endif

	i=0;
	link = chain->link;
	while(i++ < linkIndex) {
		link = link->nLink;
	}
	return link;
}

int32_t chain_getLength(Chain *chain) {
	return chain->linkNumber;
}

Block **chain_getBlockChain(Chain *chain, int32_t *blockNumber) {
	int32_t i;
	Link *link;
	End *end;
	Block *block;
	struct List *blocks = constructEmptyList(0, NULL);
	for(i=0; i<chain_getLength(chain); i++) {
		link = chain_getLink(chain, i);
		end = link_getLeft(link);
		block = end_getBlock(end);
		if(block != NULL) {
			listAppend(blocks, block);
		}
	}
	if(chain_getLength(chain) > 0) {
		link = chain_getLink(chain, chain_getLength(chain)-1);
		end = link_getRight(link);
		block = end_getBlock(end);
		if(block != NULL) {
			listAppend(blocks, block);
		}
	}
	i = sizeof(void *)*(blocks->length+1);
	Block **blockChain = memcpy(malloc(i), blocks->list, i);
	*blockNumber = blocks->length;
	destructList(blocks);
	return blockChain;
}

Name chain_getName(Chain *chain) {
	return chain->name;
}

Net *chain_getNet(Chain *chain) {
	return chain->net;
}

double chain_getAverageInstanceBaseLength(Chain *chain) {
	Block **blocks;
	int32_t i, j, l;
	double k = 0.0;
	blocks = chain_getBlockChain(chain, &i);
	l = 0;
	for(j=0; j<i; j++) {
		k += block_getLength(blocks[j]);
	}
	free(blocks);
	for(i=0; i<chain_getLength(chain); i++) {
		Net *nestedNet = group_getNestedNet(link_getGroup(chain_getLink(chain, i)));
		if(nestedNet != NULL) {
			k += (net_getTotalBaseLength(nestedNet) / net_getSequenceNumber(nestedNet));
		}
	}
	return k;
}

/*
 * Private functions
 */

void chain_addLink(Chain *chain, Link *childLink) {
	Link *pLink;
	if(chain->linkNumber != 0) {
		pLink = chain_getLink(chain, chain->linkNumber -1);
		pLink->nLink = childLink;
		childLink->pLink = pLink;
		assert(link_getRight(pLink) != link_getLeft(childLink));
		assert(end_getBlock(link_getRight(pLink)) != NULL);
		assert(end_getBlock(link_getLeft(childLink)) != NULL);
		assert(end_getBlock(link_getRight(pLink)) == end_getBlock(link_getLeft(childLink)));
	}
	else {
		childLink->pLink = NULL;
		chain->link = childLink;
	}
	childLink->nLink = NULL;
	childLink->linkIndex = chain->linkNumber++;
}

/*
 * Serialisation functions.
 */

void chain_writeBinaryRepresentation(Chain *chain, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	Link *link;
	binaryRepresentation_writeElementType(CODE_CHAIN, writeFn);
	binaryRepresentation_writeName(chain_getName(chain), writeFn);
	link = chain_getLink(chain, 0);
	while(link != NULL) {
		link_writeBinaryRepresentation(link, writeFn);
		link = link_getNextLink(link);
	}
}

Chain *chain_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	Chain *chain;

	chain = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_CHAIN) {
		binaryRepresentation_popNextElementType(binaryString);
		chain = chain_construct2(binaryRepresentation_getName(binaryString), net);
		while(link_loadFromBinaryRepresentation(binaryString, chain) != NULL);
	}
	return chain;
}

Chain *chain_getStaticNameWrapper(Name name) {
	static Chain chain;
	chain.name = name;
	return &chain;
}
