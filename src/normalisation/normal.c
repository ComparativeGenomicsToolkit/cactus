/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "normal.h"
#include "sonLib.h"

void removeTrivialLinks(Flower *flower) {
    Chain *chain;
    Flower_ChainIterator *chainIt = flower_getChainIterator(flower);
    stList *chainsToDelete = stList_construct3(0, (void(*)(void *)) chain_destruct);
    while ((chain = flower_getNextChain(chainIt)) != NULL) {
        Link *link = chain_getLast(chain);
        while(link != NULL) {
        	Link *link2 = link;
        	link = link_getPreviousLink(link);
            link_mergeIfTrivial(link2);
        }
        if (chain_getLength(chain) == 0) {
            stList_append(chainsToDelete, chain);
        }
    }
    flower_destructChainIterator(chainIt);
    stList_destruct(chainsToDelete);
}

static void promoteChainsThatExtendHigherLevelChainsP(Flower *flower, Group *parentGroup) {
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
        End *_3End = link_get3End(chain_getFirst(chain));
        End *_5End = link_get5End(chain_getLast(chain));
        if (end_isStubEnd(_3End) || end_isStubEnd(_5End)) { //Is part of higher chain..
            chain_promote(chain);
        }
    }
    stList_destruct(chains);
}

static stList *getNestedFlowers(Flower *flower) {
    /*
     * Gets a list of nested flowers for the given flower.
     */
    Group *group;
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    stList *childFlowers = stList_construct();
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (!group_isLeaf(group)) {
            stList_append(childFlowers, group_getNestedFlower(group));
        }
    }
    flower_destructGroupIterator(groupIt);
    return childFlowers;
}

void promoteNestedChainsThatExtendChains(Flower *flower) {
    stList *childFlowers = getNestedFlowers(flower);
    while (stList_length(childFlowers) > 0) {
        Flower *childFlower = stList_pop(childFlowers);
        promoteChainsThatExtendHigherLevelChainsP(childFlower, flower_getParentGroup(childFlower));
    }
    stList_destruct(childFlowers);
}

static int promoteChainsFillParentsP(Chain *chain1, Chain *chain2) {
    return chain_getLength(chain1) - chain_getLength(chain2);
}

static int promoteChainsFillParentsP2(Block *block1, Block *block2) {
    return block_getLength(block1) * block_getInstanceNumber(block1) - block_getLength(block2)
            * block_getInstanceNumber(block2);
}

static void promoteChainsToFillParentsP(Flower *flower, Group *parentGroup, int64_t maxNumberOfChains) {
    Chain *chain;
    Flower_ChainIterator *chainIt = flower_getChainIterator(flower);
    stList *chains = stList_construct();
    while ((chain = flower_getNextChain(chainIt)) != NULL) {
        assert(chain_getLength(chain) >= 1);
        assert(chain_getFlower(chain) == flower);
        stList_append(chains, chain);
    }
    flower_destructChainIterator(chainIt);
    stList_sort(chains, (int(*)(const void *, const void *)) promoteChainsFillParentsP);
    Flower *parentFlower = group_getFlower(parentGroup);
    assert(group_isTangle(parentGroup)); //Completely redundant check for old bug
    assert(group_getLink(parentGroup) == NULL);
    int64_t chainLength = INT64_MAX;
    while (group_isTangle(parentGroup) && stList_length(chains) > 0 && (flower_getChainNumber(parentFlower)
            + flower_getTrivialChainNumber(parentFlower) < maxNumberOfChains || flower_getAttachedStubEndNumber(flower)
            == 0)) {
        chain = stList_pop(chains);
        assert(chainLength >= chain_getLength(chain));
        assert(chain_getFlower(chain) == flower);
        assert(flower_getParentGroup(flower) == parentGroup);
        assert(group_getLink(parentGroup) == NULL); //Should never become a chain as if tangle it should stay a tangle as all chains are already maximal.
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
        if (block_isTrivialChain(block)) {
            stList_append(blocks, block);
        }
    }
    flower_destructBlockIterator(blockIt);
    stList_sort(blocks, (int(*)(const void *, const void *)) promoteChainsFillParentsP2);
    int64_t blockCoverage = INT64_MAX;
    while (group_isTangle(parentGroup) && stList_length(blocks) > 0 && (flower_getChainNumber(parentFlower)
            + flower_getTrivialChainNumber(parentFlower) < maxNumberOfChains || flower_getAttachedStubEndNumber(flower)
            == 0)) {
        block = stList_pop(blocks);
        assert(block_getFlower(block) == flower);
        assert(flower_getParentGroup(flower) == parentGroup);
        assert(group_getLink(parentGroup) == NULL); //Should never become a chain as if tangle it should stay a tangle as all chains are already maximal.
        assert(blockCoverage >= block_getLength(block) * block_getInstanceNumber(block));
        blockCoverage = block_getLength(block) * block_getInstanceNumber(block);
        block_promote(block);
    }
    stList_destruct(blocks);
}

void promoteNestedChainsToFillFlower(Flower *flower, int64_t maxNumberOfChains) {
    stList *childFlowers = getNestedFlowers(flower);
    while (stList_length(childFlowers) > 0) {
        Flower *childFlower = stList_pop(childFlowers);
        Group *group = flower_getParentGroup(childFlower);
        assert(group_getFlower(group) == flower);
        if (group_isTangle(group)) {
            assert(group_getFlower(group) == flower);
            assert(group_getLink(group) == NULL);
            promoteChainsToFillParentsP(childFlower, group, maxNumberOfChains);
        }
    }
    stList_destruct(childFlowers);
}

void normalise(Flower *flower, int64_t maxNumberOfChains) {
    promoteNestedChainsThatExtendChains(flower);
    promoteNestedChainsToFillFlower(flower, maxNumberOfChains);

    //Now we normalise the children of the flower..
    stList *childFlowers = getNestedFlowers(flower);
    while (stList_length(childFlowers) > 0) {
        Flower *childFlower = stList_pop(childFlowers);
        if (!flower_deleteIfEmpty(childFlower)) { //If we delete the flower we need not run the remaining functions..
            flower_makeTerminalNormal(childFlower);
            flower_removeIfRedundant(childFlower); //This may destroy the flower, but nots its children..
        }
    }
    stList_destruct(childFlowers);
    //Remove any trivial links.
    removeTrivialLinks(flower);
    //Finally we make the flower itself terminal normal (this is to handle the case for the root).
    flower_makeTerminalNormal(flower);
}
