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
    chain = st_malloc(sizeof(Chain));
    chain->name = name;
    chain->net = net;
    chain->link = NULL;
    chain->linkNumber = 0;
    net_addChain(net, chain);
    return chain;
}

void chain_destruct(Chain *chain) {
    net_removeChain(chain_getNet(chain), chain);
    if (chain->link != NULL) {
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

    i = 0;
    link = chain->link;
    while (i++ < linkIndex) {
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
    for (i = 0; i < chain_getLength(chain); i++) {
        link = chain_getLink(chain, i);
        end = link_get5End(link);
        block = end_getBlock(end);
        if (block != NULL) {
            assert(block_getOrientation(block));
            listAppend(blocks, block);
        }
    }
    if (chain_getLength(chain) > 0) {
        link = chain_getLink(chain, chain_getLength(chain) - 1);
        end = link_get3End(link);
        block = end_getBlock(end);
        if (block != NULL) {
            assert(block_getOrientation(block));
            listAppend(blocks, block);
        }
    }
    i = sizeof(void *) * (blocks->length + 1);
    Block **blockChain = memcpy(st_malloc(i), blocks->list, i);
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
    for (j = 0; j < i; j++) {
        k += block_getLength(blocks[j]);
    }
    free(blocks);
    for (i = 0; i < chain_getLength(chain); i++) {
        Net *nestedNet = group_getNestedNet(link_getGroup(chain_getLink(chain,
                i)));
        if (nestedNet != NULL) {
            k += (net_getTotalBaseLength(nestedNet) / net_getSequenceNumber(
                    nestedNet));
        }
    }
    return k;
}

void chain_check(Chain *chain) {
    Link *link = NULL, *pLink = NULL;
    int32_t i;
    assert(chain_getLength(chain) > 0);
    for (i = 0; i < chain_getLength(chain); i++) {
        link = chain_getLink(chain, i);
        //That each link is properly contained in the chain.
        assert(chain == link_getChain(link));
        End *_5End = link_get5End(link);
        End *_3End = link_get3End(link);
        assert(_5End != NULL);
        assert(_3End != NULL);
        //Links and the contained ends are properly connected.
        assert(group_getLink(end_getGroup(_5End)) == link);
        assert(group_getLink(end_getGroup(_3End)) == link);
        //Check the orientations
        assert(end_getOrientation(_5End));
        assert(end_getOrientation(_3End));

        //Check stub ends are not free stubs.
        if (end_isStubEnd(_5End)) {
            assert(end_isAttached(_5End));
        }
        if (end_isStubEnd(_3End)) {
            assert(end_isAttached(_3End));
        }

        //That the ends are consistently oriented
        assert(end_getSide(_5End) != end_getSide(_3End));
        assert(!end_getSide(_5End));
        assert(end_getSide(_3End));

        //That each contiguous pair of link groups are bridged by a block.
        if (pLink != NULL) {
            assert(end_isBlockEnd(link_get3End(pLink)));
            assert(end_isBlockEnd(_5End));
            assert(end_getOtherBlockEnd(link_get3End(pLink)) == _5End);
        } else {
            if (end_isBlockEnd(_5End)) {
                //If a block end is at the 5 prime end of a chain the other end of the
                //block is not in a link group (otherwise the chain is not maximal).
                assert(group_getLink(end_getGroup(end_getOtherBlockEnd(_5End))) == NULL);
            }
        }
        pLink = link;
    }
    //If a block end is at the 3 prime end of a chain the other end of the
    //block is not in a link group (otherwise the chain is not maximal).
    assert(link != NULL);
    if (end_isBlockEnd(link_get3End(link))) {
        assert(group_getLink(end_getGroup(end_getOtherBlockEnd(link_get3End(link)))) == NULL);
    }
}

/*
 * Private functions
 */

void chain_addLink(Chain *chain, Link *childLink) {
    Link *pLink;
    if (chain->linkNumber != 0) {
        pLink = chain_getLink(chain, chain->linkNumber - 1);
        pLink->nLink = childLink;
        childLink->pLink = pLink;
        assert(link_get3End(pLink) != link_get5End(childLink));
        assert(end_getBlock(link_get3End(pLink)) != NULL);
        assert(end_getBlock(link_get5End(childLink)) != NULL);
        assert(end_getBlock(link_get3End(pLink)) == end_getBlock(link_get5End(childLink)));
    } else {
        childLink->pLink = NULL;
        chain->link = childLink;
    }
    childLink->nLink = NULL;
    childLink->linkIndex = chain->linkNumber++;
}

void chain_setNet(Chain *chain, Net *net) {
    net_removeChain(chain_getNet(chain), chain);
    chain->net = net;
    net_addChain(net, chain);
}

/*
 * Serialisation functions.
 */

void chain_writeBinaryRepresentation(Chain *chain, void(*writeFn)(
        const void * ptr, size_t size, size_t count)) {
    Link *link;
    binaryRepresentation_writeElementType(CODE_CHAIN, writeFn);
    binaryRepresentation_writeName(chain_getName(chain), writeFn);
    link = chain_getLink(chain, 0);
    while (link != NULL) {
        link_writeBinaryRepresentation(link, writeFn);
        link = link_getNextLink(link);
    }
}

Chain *chain_loadFromBinaryRepresentation(void **binaryString, Net *net) {
    Chain *chain;

    chain = NULL;
    if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_CHAIN) {
        binaryRepresentation_popNextElementType(binaryString);
        chain = chain_construct2(binaryRepresentation_getName(binaryString),
                net);
        while (link_loadFromBinaryRepresentation(binaryString, chain) != NULL)
            ;
    }
    return chain;
}

Chain *chain_getStaticNameWrapper(Name name) {
    static Chain chain;
    chain.name = name;
    return &chain;
}

void chain_rename(Chain *chain, Name newName) {
    assert(net_getChain(chain_getNet(chain), newName) == NULL);
    chain->name = newName;
    net_removeChain(chain_getNet(chain), chain);
    net_addChain(chain_getNet(chain), chain);
}
