#include "normal.h"
#include "sonLib.h"

void promoteChainsThatExtendHigherLevelChains(Flower *flower) {
	Group *parentGroup = flower_getParentGroup(flower);
	if (parentGroup != NULL) {
		//Find links which need to be promoted
		Chain *chain;
		Flower_ChainIterator *chainIt = flower_getChainIterator(flower);
		stList *chains = stList_construct();
		while ((chain = flower_getNextChain(chainIt)) != NULL) {
			assert(chain_getLength(chain) >= 1);
			assert(chain_getFlower(chain) == flower);
			stList_append(chains, chain);
		}
		flower_destructChainIterator(chainIt);
		while (stList_length(chains) > 0) {
			chain = stList_pop(chains);
			End *_3End = link_get3End(chain_getLink(chain, 0));
			End *_5End = link_get5End(chain_getLink(chain, chain_getLength(
					chain) - 1));
			if (end_isStubEnd(_3End) || end_isStubEnd(_5End)) { //Is part of higher chain..
				chain_promote(chain);
			}
		}
	}
}

static int promoteChainsFillParentsP(Chain *chain1, Chain *chain2) {
	return chain_getLength(chain1) - chain_getLength(chain2);
}

static int promoteChainsFillParentsP2(Block *block1, Block *block2) {
	return block_getLength(block1) * block_getInstanceNumber(block1)
			- block_getLength(block2) * block_getInstanceNumber(block2);
}

void promoteChainsToFillParents(Flower *flower, int32_t maxNumberOfChains) {
	Group *parentGroup = flower_getParentGroup(flower);
	if (parentGroup != NULL && group_isTangle(parentGroup)) {
		//Find links which need to be promoted
		Chain *chain;
		Flower_ChainIterator *chainIt = flower_getChainIterator(flower);
		stList *chains = stList_construct();
		while ((chain = flower_getNextChain(chainIt)) != NULL) {
			assert(chain_getLength(chain) >= 1);
			assert(chain_getFlower(chain) == flower);
			stList_append(chains, chain);
		}
		flower_destructChainIterator(chainIt);
		stList_sort(chains,
				(int(*)(const void *, const void *)) promoteChainsFillParentsP);
		Flower *parentFlower = group_getFlower(parentGroup);
		int32_t chainLength = INT32_MAX;
		while (stList_length(chains) > 0 &&
				flower_getChainNumber(parentFlower) + flower_getTrivialChainNumber(parentFlower)
				< maxNumberOfChains) {
			chain = stList_pop(chains);
			assert(chainLength >= chain_getLength(chain));
			chainLength = chain_getLength(chain);
			chain_promote(chain);
		}
		stList_destruct(chains);

		//Now promote the trivial chains.
		Block *block;
		Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
		stList *blocks = stList_construct();
		while ((block = flower_getNextBlock(blockIt)) != NULL) {
			assert(block_getFlower(block) == flower);
			if(block_isTrivialChain(block)) {
				stList_append(blocks, block);
			}
		}
		flower_destructBlockIterator(blockIt);
		stList_sort(blocks,
				(int(*)(const void *, const void *)) promoteChainsFillParentsP2);

		int32_t blockCoverage = INT32_MAX;
		while (stList_length(blocks) > 0 &&
				flower_getChainNumber(parentFlower) + flower_getTrivialChainNumber(parentFlower)
				< maxNumberOfChains) {
			block = stList_pop(blocks);
			assert(blockCoverage >= block_getLength(block) * block_getInstanceNumber(block));
			blockCoverage = block_getLength(block) * block_getInstanceNumber(block);
			block_promote(block);
		}
		stList_destruct(blocks);
	}
}
