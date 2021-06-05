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
    assert(!group_isLink(group));
    assert(chain_getFlower(parentChain) == group_getFlower(group));

    // Orient the ends in the link
    _3End = end_getPositiveOrientation(_3End);
    _5End = end_getPositiveOrientation(_5End);
    assert(_3End != _5End);

    // Get a list of the ends in the group excluding the 5' and 3' ends of the link
    stList *extraEnds = stList_construct();
    Group_EndIterator *endIt = group_getEndIterator(group);
    End *end;
    while((end = group_getNextEnd(endIt)) != NULL) {
        end = end_getPositiveOrientation(end);
        if(end != _3End && end != _5End) {
            stList_append(extraEnds, end);
        }
    }
    group_destructEndIterator(endIt);

    // Order the ends in the group so that the 5' and 3' ends are the first
    // and second, respectively, ends in the list
    group->firstEnd = NULL;
    while(stList_length(extraEnds) > 0) {
        group_addEnd(group, stList_pop(extraEnds));
    }
    group_addEnd(group, _3End);
    group_addEnd(group, _5End);
    stList_destruct(extraEnds);

    // Finally set as a link
    group_setLink(group, 1);
    group->nLink = NULL;
    assert(group_isLink(group));
    group->flowerOrChain = parentChain; // Now set the pointer to the chain
    chain_addLink(parentChain, group); //will set the link indices.

    // Doe some checks
    assert(link_get5End(group) == _5End);
    assert(link_get3End(group) == _3End);
    assert(link_getGroup(group) == group);
    assert(link_getChain(group) == parentChain);
    assert(group_getFlower(group) == chain_getFlower(parentChain));

    return group;
}

Link *link_getNextLink(Link *link) {
    return link->nLink;
}

Group *link_getGroup(Link *link) {
    return link; //link->group;
}

End *link_get3End(Link *link) {
    return *getNextEndPointer(link->firstEnd); //link->_3End;
}

End *link_get5End(Link *link) {
    return link->firstEnd; //_5End;
}

Chain *link_getChain(Link *link) {
    return link->flowerOrChain; //link->chain;
}

/*
 * Private functions.
 */

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

void link_destruct(Link *link) {
	Chain *chain = link_getChain(link);
    Link *link2 = link->nLink;
    while (link2 != NULL) {
        Link *link3 = link2;
        link2 = link2->nLink;
        group_setLink(link3, 0);
        link3->flowerOrChain = chain_getFlower(chain);
        link3->nLink = NULL;
    }
    link->nLink = NULL;
    link->flowerOrChain = chain_getFlower(chain);
    group_setLink(link, 0);
}

