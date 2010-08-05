#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic link functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Link *link_construct(End *_3End, End *_5End, Group *group,
        Chain *parentChain) {
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

int32_t link_getIndex(Link *link) {
    return link->linkIndex;
}

/*
 * Private functions.
 */

void link_destruct(Link *link) {
    Link *link2;
    link_getChain(link)->linkNumber = link->linkIndex;
    if (link->pLink == NULL) {
        link_getChain(link)->link = NULL;
    } else {
        link->pLink->nLink = NULL;
    }
    while (link != NULL) {
        group_setLink(link_getGroup(link), NULL);
        link2 = link;
        link = link->nLink;
        free(link2);
    }
}

static void link_splitP(struct List *list, Flower *flower) {
    if (list->length > 0) {
        Chain *chain = chain_construct(flower);
        int32_t i;
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
    int32_t i = 0;
    while (i < chain_getLength(chain)) {
        Link *link2 = chain_getLink(chain, i++);
        if (link2 == link) {
            break;
        }
        listAppend(list1, link_get3End(link2));
        listAppend(list1, link_get5End(link2));
    }
    while (i < chain_getLength(chain)) {
        Link *link2 = chain_getLink(chain, i++);
        assert(link2 != link);
        listAppend(list2, link_get3End(link2));
        listAppend(list2, link_get5End(link2));
    }
    assert(list1->length + list2->length + 2 == chain_getLength(chain)*2);
    Flower *flower = chain_getFlower(chain);
    chain_destruct(chain);
    link_splitP(list1, flower);
    link_splitP(list2, flower);
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
        group = flower_getGroup(chain_getFlower(chain), binaryRepresentation_getName(
                binaryString));
        leftEnd = flower_getEnd(chain_getFlower(chain), binaryRepresentation_getName(
                binaryString));
        rightEnd = flower_getEnd(chain_getFlower(chain),
                binaryRepresentation_getName(binaryString));
        link = link_construct(leftEnd, rightEnd, group, chain);
    }
    return link;
}
