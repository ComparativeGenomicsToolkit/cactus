/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic group functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static void group_setBit(Group *group, int bit, bool value) {
    group->bits &= ~(1UL << bit); // first clear the existing value
    if(value) {
        group->bits |= 1UL << bit;// now set the new value
    }
}

static bool group_getBit(Group *group, int bit) {
    return (group->bits >> bit) & 1;
}

bool group_isLeaf(Group *group) {
    return group_getBit(group, 0);
}

void group_setLeaf(Group *group, bool isLeaf) {
    group_setBit(group, 0, isLeaf);
}

bool group_isLink(Group *group) {
    return group_getBit(group, 1);
}

void group_setLink(Group *group, bool isLink) {
    group_setBit(group, 1, isLink);
}

Group *group_construct(Flower *flower, Flower *nestedFlower) {
    Group *group = group_construct4(flower, flower_getName(nestedFlower), 0);

    // Get the ends that are in the nested flower and add them to this group
    Flower_EndIterator *iterator = flower_getEndIterator(group_getNestedFlower(group));
    End *end, *end2;
    while ((end = flower_getNextEnd(iterator)) != NULL) {
        if ((end2 = flower_getEnd(flower, end_getName(end))) != NULL) {
            end_setGroup(end2, group);
        }
    }
    flower_destructEndIterator(iterator);

    flower_setParentGroup(nestedFlower, group);

    return group;
}

Group *group_construct3(Flower *flower, Name name) {
    return group_construct4(flower, name, 1);
}

Group *group_construct2(Flower *flower) {
    return group_construct3(flower, cactusDisk_getUniqueID(flower_getCactusDisk(flower)));
}

static void copyAdjacencies(Group *group, Flower *nestedFlower) {
    /*
     * Copies of the adjacencies from the caps in the ends in the group to the nested flower.
     */
    assert(flower_getParentGroup(nestedFlower) == group);
    Group_EndIterator *endIterator = group_getEndIterator(group);
    End *end;
    while ((end = group_getNextEnd(endIterator)) != NULL) {
        End *nestedEnd = flower_getEnd(nestedFlower, end_getName(end));
        assert(nestedEnd != NULL);
        Cap *cap, *adjacentCap, *nestedCap, *nestedAdjacentCap;
        End_InstanceIterator *capIterator = end_getInstanceIterator(end);
        while ((cap = end_getNext(capIterator)) != NULL) {
            adjacentCap = cap_getAdjacency(cap);
            if (adjacentCap != NULL) {
                nestedCap = end_getInstance(nestedEnd, cap_getName(cap));
                nestedAdjacentCap = flower_getCap(nestedFlower, cap_getName(adjacentCap));
                assert(nestedCap != NULL);
                assert(nestedAdjacentCap != NULL);
                nestedAdjacentCap
                        = cap_getOrientation(adjacentCap) == cap_getOrientation(nestedAdjacentCap) ? nestedAdjacentCap
                                : cap_getReverse(nestedAdjacentCap);
                assert(cap_getOrientation(cap));
                assert(cap_getOrientation(cap) == cap_getOrientation(nestedCap));
                assert(cap_getOrientation(adjacentCap) == cap_getOrientation(nestedAdjacentCap));
                assert(end_getFlower(cap_getEnd(nestedCap)) == nestedFlower);
                assert(end_getFlower(cap_getEnd(nestedAdjacentCap)) == nestedFlower);
                cap_makeAdjacent(nestedCap, nestedAdjacentCap);
            }
        }
        end_destructInstanceIterator(capIterator);
    }
    group_destructEndIterator(endIterator);
}

Flower *group_makeEmptyNestedFlower(Group *group) {
    assert(group_isLeaf(group));
    group_setLeaf(group, 0);
    Flower *nestedFlower = flower_construct2(group_getName(group), flower_getCactusDisk(group_getFlower(group)));
    flower_setParentGroup(nestedFlower, group);
    return nestedFlower;
}

Flower *group_makeNestedFlower(Group *group) {
    Flower *nestedFlower = group_makeEmptyNestedFlower(group);
    Group *nestedGroup = group_construct2(nestedFlower);
    //Add the ends to the nested flower.
    Group_EndIterator *endIterator = group_getEndIterator(group);
    End *end;
    while ((end = group_getNextEnd(endIterator)) != NULL) {
        assert(end_getOrientation(end));
        end_setGroup(end_copyConstruct(end, nestedFlower), nestedGroup);
    }
    group_destructEndIterator(endIterator);
    //Now add adjacencies between the caps, mirroring the parent adjacencies.
    copyAdjacencies(group, nestedFlower);
    //Create the trivial chain for the ends..
    group_constructChainForLink(nestedGroup);
    assert(group_getTotalBaseLength(group) == flower_getTotalBaseLength(nestedFlower));
    return nestedFlower;
}

End **getNextEndPointer(End *end) {
    if(end_partOfBlock(end)) {
        if(end_left(end)) {
            return &(end_getBlockEndContents(end)->nEndLeft);
        }
        else {
            return &(end_getBlockEndContents(end)->nEndRight);
        }
    }
    return &(end_getContents(end)->nEnd);
}

void group_addEnd(Group *group, End *end) {
    end = end_getPositiveOrientation(end);

    // Get the pointer we need
    End **end2;
    if(group_isLink(group)) {
        end2 = getNextEndPointer(*getNextEndPointer(group->firstEnd));
    }
    else {
        end2 = &(group->firstEnd);
    }

    End **nEnd = getNextEndPointer(end);
    *nEnd = *end2;
    *end2 = end;
}

void group_destruct(Group *group) {
    //Detach from the parent flower.
    flower_removeGroup(group_getFlower(group), group);
    while (!group_isEmpty(group)) {
        end_setGroup(group_getFirstEnd(group), NULL);
    }
    //Free the memory
    free(group);
}

Flower *group_getFlower(Group *group) {
    return group_isLink(group) ? chain_getFlower(group->flowerOrChain) : group->flowerOrChain;
}

Name group_getName(Group *group) {
    return group->name;
}

Flower *group_getNestedFlower(Group *group) {
    return group_isLeaf(group) ? NULL : cactusDisk_getFlower(flower_getCactusDisk(group_getFlower(group)), group->name);
}

Link *group_getLink(Group *group) {
    return group_isLink(group) ? group : NULL; //group->link;
}

bool group_isTangle(Group *group) {
    return !group_isLink(group); //group_getLink(group) == NULL;
}

End *group_getFirstEnd(Group *group) {
    return group->firstEnd; //stSortedSet_getFirst(group->ends);
}

End *group_getEnd(Group *group, Name name) {
    Flower *flower = group_getFlower(group);
    assert(flower != NULL);
    End *end = flower_getEnd(flower, name);
    if(end != NULL) {
        assert(end_getName(end) == name);
        assert(end == end_getPositiveOrientation(end));
        if(end_getGroup(end) == group) {
            return end;
        }
    }
    return NULL;
}

bool group_isEmpty(Group *group) {
    return group->firstEnd == NULL;
}

int64_t group_getEndNumber(Group *group) {
    End *end = group->firstEnd;
    int64_t totalEnds = 0;
    while(end != NULL) {
        totalEnds++; end = *getNextEndPointer(end);
    }
    return totalEnds;
}

Group_EndIterator *group_getEndIterator(Group *group) {
    Group_EndIterator *iterator = st_malloc(sizeof(Group_EndIterator));
    iterator->group = group;
    iterator->end = group->firstEnd;
    return iterator;
}

End *group_getNextEnd(Group_EndIterator *endIterator) {
    End *end = endIterator->end;
    if(end == NULL) {
        return NULL;
    }
    endIterator->end = *getNextEndPointer(end);
    assert(end_getPositiveOrientation(end) == end);
    return end;
}

void group_destructEndIterator(Group_EndIterator *endIterator) {
    free(endIterator);
    //stSortedSet_destructIterator(endIterator);
}

int64_t group_getTotalBaseLength(Group *group) {
    Group_EndIterator *endIterator = group_getEndIterator(group);
    End *end;
    int64_t totalLength = 0;
    while ((end = group_getNextEnd(endIterator)) != NULL) {
        End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(instanceIterator)) != NULL) {
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            if (!cap_getSide(cap) && cap_getSequence(cap) != NULL) {
                Cap *cap2 = cap_getAdjacency(cap);
                assert(cap2 != NULL);
                assert(cap_getStrand(cap2));
                assert(cap_getSide(cap2));
                assert(end_getGroup(cap_getEnd(cap2)) == group);
                int64_t length = cap_getCoordinate(cap2) - cap_getCoordinate(cap) - 1;
                assert(length >= 0);
                totalLength += length;
            }
        }
        end_destructInstanceIterator(instanceIterator);
    }
    group_destructEndIterator(endIterator);
    return totalLength;
}

void group_check(Group *group) {
    Flower *flower = group_getFlower(group);

    //Check flower and group properly connected.
    cactusCheck(flower_getGroup(flower, group_getName(group)) == group);

    Group_EndIterator *endIterator = group_getEndIterator(group);
    End *end;
    int64_t nonFree = 0;
    while ((end = group_getNextEnd(endIterator)) != NULL) {
        //That the ends of the groups are doubly linked to the ends (so every end is in only one link).
        cactusCheck(end_getGroup(end) == group);
        if (end_isAttached(end) || end_isBlockEnd(end)) {
            cactusCheck(end_isBlockEnd(end) || (end_isStubEnd(end) && end_isAttached(end)));
            nonFree++;
        }
    }
    group_destructEndIterator(endIterator);

    Link *link = group_getLink(group);
    if (nonFree == 2) { //We get rid of this now, as the reference can create new links, which we do not want to count as chains.
        //cactusCheck(link != NULL); // has only two non-free ends, is a link therefore
        //cactusCheck(group_isLink(group));
    } else {
        cactusCheck(group_isTangle(group));
        cactusCheck(link == NULL); // can not be a link!
    }

    if (group_isLeaf(group)) { //If terminal has no nested flower
        cactusCheck(group_getNestedFlower(group) == NULL);
    } else { //else that any nested flower contains the correct set of stub ends.
        Flower *nestedFlower = group_getNestedFlower(group);
        cactusCheck(nestedFlower != NULL);
        endIterator = group_getEndIterator(group);
        while ((end = group_getNextEnd(endIterator)) != NULL) {
            End *end2 = flower_getEnd(nestedFlower, end_getName(end));
            cactusCheck(end2 != NULL);
            cactusCheck(end_isStubEnd(end2));
            if (end_isBlockEnd(end) || end_isAttached(end)) {
                end_isAttached(end2);
            } else {
                end_isFree(end2);
            }
        }
        group_destructEndIterator(endIterator);
    }
}

void group_constructChainForLink(Group *group) {
    //The following constructs a chain, if necessary.
    if (group_getLink(group) == NULL) {
        Group_EndIterator *endIt = group_getEndIterator(group);
        End *end;
        int64_t i = 0;
        while ((end = group_getNextEnd(endIt)) != NULL) {
            if (end_isAttached(end) || end_isBlockEnd(end)) {
                i++;
            }
        }
        group_destructEndIterator(endIt);
        if (i == 2) {
            End *_3End = NULL, *_5End = NULL;
            endIt = group_getEndIterator(group);
            while ((end = group_getNextEnd(endIt)) != NULL) {
                if (end_isAttached(end) || end_isBlockEnd(end)) {
                    if (end_getSide(end)) {
                        assert(_5End == NULL);
                        _5End = end;
                    } else {
                        assert(_3End == NULL);
                        _3End = end;
                    }
                }
            }
            assert(_3End != NULL && _5End != NULL);
            assert(!end_getSide(_3End) && end_getSide(_5End));
            group_destructEndIterator(endIt);
            Link *_3Link = NULL;
            if (end_isBlockEnd(_3End)) {
                _3Link = group_getLink(end_getGroup(end_getOtherBlockEnd(_3End)));
            }
            Chain *chain;
            if (_3Link != NULL) {
                chain = link_getChain(_3Link);
                assert(_3Link == chain_getLast(chain));
            } else {
                chain = chain_construct(group_getFlower(group));
            }
            link_construct(_3End, _5End, group, chain);

            //Now see if we must join the chain
            if (end_isBlockEnd(_5End)) {
                assert(end_getSide(_5End));
                assert(!end_getSide(end_getOtherBlockEnd(_5End)));
                Link *_5Link = group_getLink(end_getGroup(end_getOtherBlockEnd(_5End)));
                if (_5Link != NULL) {
                    Chain *_5Chain = link_getChain(_5Link);
                    assert(_5Chain != NULL);
                    assert(_5Link == chain_getFirst(_5Chain));
                    if (chain != _5Chain) { //We don't want to merge a circle
                        chain_join(chain, _5Chain);
                    }
                }
            }
            assert(group_isLink(group));

        }
    }
}

int64_t group_getStubEndNumber(Group *group) {
    return group_getEndNumber(group) - group_getBlockEndNumber(group);
}

int64_t group_getAttachedStubEndNumber(Group *group) {
    Group_EndIterator *endIt = group_getEndIterator(group);
    End *end;
    int64_t i = 0;
    while ((end = group_getNextEnd(endIt)) != NULL) {
        if (end_isAttached(end)) {
            assert(end_isStubEnd(end));
            i++;
        }
    }
    group_destructEndIterator(endIt);
    return i;
}

int64_t group_getFreeStubEndNumber(Group *group) {
    int64_t i = group_getStubEndNumber(group) - group_getAttachedStubEndNumber(group);
    assert(i >= 0);
    return i;
}

int64_t group_getBlockEndNumber(Group *group) {
    Group_EndIterator *endIt = group_getEndIterator(group);
    End *end;
    int64_t i = 0;
    while ((end = group_getNextEnd(endIt)) != NULL) {
        if (end_isBlockEnd(end)) {
            i++;
        }
    }
    group_destructEndIterator(endIt);
    return i;
}

/*
 * Private functions.
 */

Group *group_construct4(Flower *flower, Name name, bool terminalGroup) {
    Group *group;
    group = st_calloc(1, sizeof(Group));
    group_setLeaf(group, terminalGroup);
    assert(group_isLeaf(group) == terminalGroup);
    assert(!group_isLink(group));
    group->flowerOrChain = flower;
    group->name = name;
    flower_addGroup(flower, group);

    assert(group_getFlower(group) == flower);

    return group;
}

void group_removeEnd(Group *group, End *end) {
    End **endP = &(group->firstEnd);
    while(*endP != NULL) {
        if(end_getName(end) == end_getName(*endP)) {
            End **nEnd = getNextEndPointer(*endP);
            (*endP) = *nEnd;  // Splice it out
            *nEnd = NULL;
            return;
        }
        endP = getNextEndPointer(*endP);
    }
}

/*
 * Stuff to remove
 */

End *group_getPreviousEnd(Group_EndIterator *endIterator) {
    assert(0);
    //return stSortedSet_getPrevious(endIterator);
}

Group_EndIterator *group_copyEndIterator(Group_EndIterator *endIterator) {
    assert(0);
    //return stSortedSet_copyIterator(endIterator);
}

void group_setFlower(Group *group, Flower *flower) {
    assert(0);
    /*flower_removeGroup(group_getFlower(group), group);
    group->flower = flower;
    flower_addGroup(flower, group);
    Flower *nestedFlower = group_getNestedFlower(group);
    if (nestedFlower != NULL) { //we re-do this link, because the parent flower has changed.
        flower_setParentGroup(nestedFlower, group);
    }*/
}

//void group_setLink(Group *group, Link *link) {
//    assert(0);
    //argument may be NULL
    /*group->link = link;
    if (link != NULL) {
        assert(group_getEnd(group, end_getName(link_get3End(link))) == link_get3End(link));
        assert(group_getEnd(group, end_getName(link_get5End(link))) == link_get5End(link));
    }*/
//}

/*
 * Serialisation functions
 */

void group_writeBinaryRepresentation(Group *group, void(*writeFn)(const void * ptr, size_t size, size_t count)) {
    End *end;
    Group_EndIterator *iterator;

    binaryRepresentation_writeElementType(CODE_GROUP, writeFn);
    binaryRepresentation_writeBool(group_isLeaf(group), writeFn);
    binaryRepresentation_writeName(group_getName(group), writeFn);
    iterator = group_getEndIterator(group);
    while ((end = group_getNextEnd(iterator)) != NULL) {
        binaryRepresentation_writeElementType(CODE_GROUP_END, writeFn);
        binaryRepresentation_writeName(end_getName(end), writeFn);
    }
    group_destructEndIterator(iterator);
    binaryRepresentation_writeElementType(CODE_GROUP, writeFn);
}

Group *group_loadFromBinaryRepresentation(void **binaryString, Flower *flower) {
    Group *group;

    group = NULL;
    if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_GROUP) {
        binaryRepresentation_popNextElementType(binaryString);
        bool terminalGroup = binaryRepresentation_getBool(binaryString);
        Name name = binaryRepresentation_getName(binaryString);
        group = group_construct4(flower, name, terminalGroup);
        while (binaryRepresentation_peekNextElementType(*binaryString) == CODE_GROUP_END) {
            binaryRepresentation_popNextElementType(binaryString);
            end_setGroup(flower_getEnd(flower, binaryRepresentation_getName(binaryString)), group);
        }
        assert(binaryRepresentation_peekNextElementType(*binaryString) == CODE_GROUP);
        binaryRepresentation_popNextElementType(binaryString);
    }
    return group;
}



