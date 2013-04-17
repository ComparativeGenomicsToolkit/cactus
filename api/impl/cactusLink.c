/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic link functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Link *link_construct(End *_3End, End *_5End, Group *group, Chain *parentChain) {
    Link *link;
    link = st_malloc(sizeof(Link));

    _3End = end_getPositiveOrientation(_3End);
    _5End = end_getPositiveOrientation(_5End);
    assert(_3End != _5End);

    link->_3End = _3End;
    link->_5End = _5End;
    link->chain = parentChain;
    link->group = group;

    //Checks.
    assert(group_getEnd(group, end_getName(_3End)) == _3End);
    assert(group_getEnd(group, end_getName(_5End)) == _5End);
    //The following ensure that our chain is oriented along an (arbitrary) but consistent direction.
    assert(end_getSide(_3End) != end_getSide(_5End));
    assert(!end_getSide(_3End));
    assert(end_getSide(_5End));

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

End *link_get3End(Link *link) {
    return link->_3End;
}

End *link_get5End(Link *link) {
    return link->_5End;
}

Chain *link_getChain(Link *link) {
    return link->chain;
}

/*
 * Private functions.
 */

void link_destruct(Link *link) {
	group_setLink(link_getGroup(link), NULL);
    int64_t i = 1;
    Link *link2 = link->nLink;
    while (link2 != NULL) {
        group_setLink(link_getGroup(link2), NULL);
        Link *link3 = link2;
        link2 = link2->nLink;
        free(link3);
        i++;
    }
    Chain *chain = link_getChain(link);
    chain->linkNumber -= i;
    assert(chain->linkNumber >= 0);
    if (link->pLink == NULL) {
        chain->link = NULL;
        chain->endLink = NULL;
    } else {
        link->pLink->nLink = NULL;
        chain->endLink = link->pLink;
    }
    free(link);
}

static void link_splitP(struct List *list, Flower *flower) {
    if (list->length > 0) {
        Chain *chain = chain_construct(flower);
        int64_t i;
        assert(list->length % 2 == 0);
        for (i = 0; i < list->length; i += 2) {
            link_construct(list->list[i], list->list[i + 1], end_getGroup(
                    list->list[i]), chain);
        }
    }
    destructList(list);
}

void link_split(Link *link) {
    Chain *chain = link_getChain(link);
    struct List *list1 = constructEmptyList(0, NULL), *list2 =
            constructEmptyList(0, NULL);
    Link *link2 = chain_getFirst(chain);
    while (link2 != NULL) {
        if (link2 == link) {
            link2 = link_getNextLink(link2);
            break;
        }
        listAppend(list1, link_get3End(link2));
        listAppend(list1, link_get5End(link2));
        link2 = link_getNextLink(link2);
    }
    while (link2 != NULL) {
        assert(link2 != link);
        listAppend(list2, link_get3End(link2));
        listAppend(list2, link_get5End(link2));
        link2 = link_getNextLink(link2);
    }
    assert(list1->length + list2->length + 2 == chain_getLength(chain) * 2);
    Flower *flower = chain_getFlower(chain);
    chain_destruct(chain);
    link_splitP(list1, flower);
    link_splitP(list2, flower);
}

bool link_isTrivialP(End *_3End, End *_5End) {
    Cap *_3Cap;
    End_InstanceIterator *capIt = end_getInstanceIterator(_3End);
    while ((_3Cap = end_getNext(capIt)) != NULL) {
        assert(cap_getOrientation(_3Cap));
        Cap *_5Cap = cap_getAdjacency(_3Cap);
        if (_5Cap != NULL) {
            assert(cap_getEnd(_5Cap) != end_getReverse(_5End));
            if (cap_getEnd(_5Cap) != _5End) { //The adjacency must be a self adjacency as it contains no free stubs.
                assert(cap_getEnd(_5Cap) == end_getReverse(_3End));
                end_destructInstanceIterator(capIt);
                return 0;
            }
            if (llabs(cap_getCoordinate(_3Cap) - cap_getCoordinate(_5Cap)) != 1) { //There is a nontrivial adjacency
                end_destructInstanceIterator(capIt);
                return 0;
            }
        }
    }
    end_destructInstanceIterator(capIt);
    return 1;
}

bool link_isTrivial(Link *link) {
    End *_3End = link_get3End(link);
    End *_5End = link_get5End(link);
    assert(!end_getSide(_3End));
    assert(end_getSide(_5End));
    Group *group = link_getGroup(link);
    if (group_getEndNumber(group) == 2) {
        if (end_isBlockEnd(_3End) && end_isBlockEnd(_5End)) { //Must both be block ends
            if (end_getInstanceNumber(_3End) == end_getInstanceNumber(_5End)) { //Must each be connected to other.
                return link_isTrivialP(_3End, _5End) && link_isTrivialP(_5End,
                        _3End);
            }
        }
    }
    return 0;
}

bool link_mergeIfTrivial(Link *link) {
    Flower *flower = chain_getFlower(link_getChain(link));
    (void)flower;
    assert(flower_builtBlocks(flower));
    assert(!flower_builtTrees(flower));
    assert(!flower_builtFaces(flower));
    if (link_isTrivial(link)) {
        End *_3End = link_get3End(link);
        End *_5End = link_get5End(link);
        Group *group = link_getGroup(link);
        assert(group_getEndNumber(group) == 2); //Can not contain any stubs!
        Chain *chain = link_getChain(link);
        Block *_5Block = end_getBlock(_3End);
        Block *_3Block = end_getBlock(_5End);
        Flower *flower = group_getFlower(group);
        CactusDisk *cactusDisk = flower_getCactusDisk(flower);

        //First eliminate the link
        if (link->pLink != NULL) {
            link->pLink->nLink = link->nLink;
        }
        if (link->nLink != NULL) {
            link->nLink->pLink = link->pLink;
        }
        if (chain->link == link) {
            chain->link = link->nLink;
        }
        if(chain->endLink == link) {
            chain->endLink = link->pLink;
        }
        chain->linkNumber--;
        free(link); //We do our own destruction of the object..
        assert(chain_getLength(chain) >= 0);

        //Get rid of the block ends by making a new block..

        //First prepare what we need for the new block
        End *outer5End = block_get5End(_5Block);
        End *outer3End = block_get3End(_3Block);
        assert(end_getOrientation(outer5End));
        assert(end_getOrientation(outer3End));
        int64_t newBlockLength = block_getLength(_5Block) + block_getLength(
                _3Block);
        stHash *newSegments = stHash_construct();
        //This works out the caps for the new segments..
        End_InstanceIterator *capIt = end_getInstanceIterator(outer5End);
        Cap *_5Cap;
        while ((_5Cap = end_getNext(capIt)) != NULL) {
            assert(cap_getOrientation(_5Cap));
            Cap *_3Cap = cap_getOtherSegmentCap(cap_getAdjacency(
                    cap_getOtherSegmentCap(_5Cap)));
            assert(cap_getEnd(_3Cap) == outer3End);
            assert(cap_getOrientation(_3Cap)); //redundant
            stHash_insert(newSegments, _5Cap, _3Cap);
        }
        end_destructInstanceIterator(capIt);

        //Now get rid of the blocks and the segments and outer ends.
        block_destruct(_5Block);
        block_destruct(_3Block);
        end_destruct(_3End);
        end_destruct(_5End);

        //Construct the merged block
        Block *mergedBlock = block_construct2(
                cactusDisk_getUniqueID(cactusDisk), newBlockLength, outer5End,
                outer3End, flower);
        capIt = end_getInstanceIterator(outer5End);
        while ((_5Cap = end_getNext(capIt)) != NULL) {
            Cap *_3Cap = stHash_remove(newSegments, _5Cap);
            assert(_3Cap != NULL);
            assert(cap_getEnd(_3Cap) == outer3End);
            segment_construct3(cactusDisk_getUniqueID(cactusDisk), mergedBlock,
                    _5Cap, _3Cap);
        }
        end_destructInstanceIterator(capIt);
        assert(stHash_size(newSegments) == 0);
        stHash_destruct(newSegments);

        //Finally eliminate the group
        if (!group_isLeaf(group)) {
            Flower *nestedFlower = group_getNestedFlower(group);
            assert(flower_isTerminal(nestedFlower));
            cactusDisk_deleteFlowerFromDisk(cactusDisk, nestedFlower);
            flower_destruct(nestedFlower, 0);
        }
        assert(group_getEndNumber(group) == 0);
        group_destruct(group);

        return 1;
    }
    return 0;
}

/*
 * Serialisation functions.
 */

void link_writeBinaryRepresentation(Link *link, void(*writeFn)(
        const void * ptr, size_t size, size_t count)) {
    binaryRepresentation_writeElementType(CODE_LINK, writeFn);
    binaryRepresentation_writeName(group_getName(link_getGroup(link)), writeFn);
    binaryRepresentation_writeName(end_getName(link_get3End(link)), writeFn);
    binaryRepresentation_writeName(end_getName(link_get5End(link)), writeFn);
}

Link *link_loadFromBinaryRepresentation(void **binaryString, Chain *chain) {
    Group *group;
    End *leftEnd;
    End *rightEnd;
    Link *link;

    link = NULL;
    if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_LINK) {
        binaryRepresentation_popNextElementType(binaryString);
        group = flower_getGroup(chain_getFlower(chain),
                binaryRepresentation_getName(binaryString));
        leftEnd = flower_getEnd(chain_getFlower(chain),
                binaryRepresentation_getName(binaryString));
        rightEnd = flower_getEnd(chain_getFlower(chain),
                binaryRepresentation_getName(binaryString));
        link = link_construct(leftEnd, rightEnd, group, chain);
    }
    return link;
}
