/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic end functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Bit twiddling using the "bits" char of the end/block
 * 0: orientation (is the forward as opposed to the reverse copy)
 * 1: part_of_block (is part of a block)
 * 2: is_block (is the block)
 * 3: left (is the left end of a end-block-end)
 * 4: is_attached (for non-segment blocks)
 * 5: side (is the 5' or 3' copy)
*/

static void end_setBit(End *end, int bit, bool value) {
    end->bits &= ~(1UL << bit); // first clear the existing value
    if(value) {
        end->bits |= 1UL << bit;// now set the new value
    }
}

static void end_setBitForwardAndReverse(End *end, int bit, bool value, bool invertReverse) {
    end_setBit(end, bit, value);
    end_setBit(end_getReverse(end), bit, value ^ invertReverse);
}

static bool end_getBit(End *end, int bit) {
    return (end->bits >> bit) & 1;
}

bool end_getOrientation(End *end) {
    return end_getBit(end, 0);
}

bool end_partOfBlock(End *end) {
    return end_getBit(end, 1);
}

bool end_isBlock(End *end) {
    return end_getBit(end, 2);
}

bool end_isBlockEnd(End *end) {
    return end_partOfBlock(end) && !end_isBlock(end);
}

bool end_isStubEnd(End *end) {
    return !end_isBlockEnd(end);
}

bool end_left(End *end) {
    return end_getBit(end, 3);
}

bool end_isAttached(End *end) {
    return end_getBit(end, 4);
}

bool end_isFree(End *end) {
    return !end_isAttached(end);
}

bool end_getSide(End *end) {
    return end_getBit(end, 5);
}

void end_makeAttached(End *end) {
    assert(end_isStubEnd(end));
    assert(end_isFree(end));
    assert(flower_getName(end_getFlower(end)) == 0);
    if(end_getGroup(end) != NULL) {
        assert(group_isLeaf(end_getGroup(end)));
    }
    // Set the attached bit
    end_setBitForwardAndReverse(end, 4, 1, 0);
}

BlockEndContents *end_getBlockEndContents(End *end) {
    assert(end_partOfBlock(end));
    // Bits: (0) orientation / (1) part_of_block / (2) is_block / (3) left / (4) is_attached / (5) side
    switch(end->bits) {
        case 0x2B: // binary: 101011
            return (BlockEndContents *)(end+6);
        case 0xA: // binary: 001010
            return (BlockEndContents *)(end+5);
        case 0x7: // binary: 000111
            return (BlockEndContents *)(end+4);
        case 0x6: // binary: 000110
            return (BlockEndContents *)(end+3);
        case 0x3: // binary: 000011
            return (BlockEndContents *)(end+2);
        case 0x22: // binary: 100010
            return (BlockEndContents *)(end+1);
        default:
            assert(0);
            return NULL;
    }
}

EndContents *end_getContents(End *end) {
    assert(!end_partOfBlock(end));
    // Bits: (0) orientation / (1) part_of_block / (2) is_block / (3) left / (4) is_attached / (5) side
    switch(end->bits & 0x3) {
        case 0x1: // binary: 000001
            return (EndContents *)(end+2);
        case 0x0: // binary: 000000
            return (EndContents *)(end+1);
        default:
            assert(0);
            return NULL;
    }
}

End *end_construct(bool isAttached, Flower *flower) {
    return end_construct3(cactusDisk_getUniqueID(flower_getCactusDisk(flower)),
            isAttached, 1, flower);
}

End *end_construct2(bool side, bool isAttached, Flower *flower) {
    return end_construct3(cactusDisk_getUniqueID(flower_getCactusDisk(flower)),
            isAttached, side, flower);
}

End *end_construct3(Name name, int64_t isAttached,
        int64_t side, Flower *flower) {
    End *end = st_calloc(1, 2*sizeof(End) + sizeof(EndContents));
    // see above comment to decode what is set
    // Bits: (0) orientation / (1) part_of_block / (2) is_block / (3) left / (4) is_attached / (5) side
    end->bits = 1; // binary 000001
    (end+1)->bits = 0; // binary 000000
    end_setBitForwardAndReverse(end, 4, isAttached, 0);
    end_setBitForwardAndReverse(end, 5, side, 1);

    // Set attributes of the shared contents object
    end_getContents(end)->name = name;
    end_getContents(end)->flower = flower;

    // Add to the flower
    flower_addEnd(flower, end);

    // Checks
    assert(end_getSide(end) == side);
    assert(end_getSide(end_getReverse(end)) != side);
    assert(end_getReverse(end_getReverse(end)) == end);
    assert(end_getFlower(end) == flower);
    assert(end_getFlower(end_getReverse(end)) == flower);
    assert(end_isAttached(end) == isAttached);
    assert(end_isAttached(end_getReverse(end)) == isAttached);
    assert(end_isStubEnd(end));
    assert(end_isStubEnd(end_getReverse(end)));
    assert(flower_getEnd(flower, end_getName(end)) == end);

    return end;
}

End *end_copyConstruct(End *end, Flower *newFlower) {
    end = end_getPositiveOrientation(end);
    assert(flower_getEnd(newFlower, end_getName(end)) == NULL);

    End *end2 = end_construct3(end_getName(end), end_isBlockEnd(end) ? 1
            : end_isAttached(end), end_getSide(end), newFlower);
    //Copy the instances.
    End_InstanceIterator *iterator = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(iterator)) != NULL) {
        cap_copyConstruct(end2, cap);
    }
    end_destructInstanceIterator(iterator);

    return end2;
}

void end_destruct(End *end) {
    /*
     * This method is the only way to clean up ends / end-blocks. It must be called
     * on both ends of a block
     */

    //remove from flower.
    flower_removeEnd(end_getFlower(end), end);

    //remove from group.
    end_setGroup(end, NULL);

    if(!end_partOfBlock(end)) { // If not part of a block
        //remove instances
        Cap *cap;
        while ((cap = end_getFirst(end)) != NULL) {
            cap_destruct(cap);
        }

        free(end_getOrientation(end) ? end : end_getReverse(end));
    }
    else if(end_left(end)) { // is the left end of a block
        Block *block = end_getBlock(end);

        //remove instances
        Segment *segment;
        while((segment = block_getFirst(block)) != NULL) {
            segment_destruct(segment);
        }

        free(block_getOrientation(block) ? block-2 : block-3);
    }
}

Name end_getName(End *end) {
    if(end_partOfBlock(end)) {
        return end_getBlockEndContents(end)->name + (end_left(end) ? 0 : 2);
    }
    return end_getContents(end)->name;
}

End *end_getPositiveOrientation(End *end) {
    return end_getOrientation(end) ? end : end_getReverse(end);
}

End *end_getReverse(End *end) {
    return end_getOrientation(end) ? end+1 : end-1; 
}

Flower *end_getFlower(End *end) {
    if(end_partOfBlock(end)) {
        return end_getBlockEndContents(end)->flower;
    }
    return end_getContents(end)->flower;
}

Block *end_getBlock(End *end) {
    if(end_partOfBlock(end)) {
        assert(!end_isBlock(end));
        return end_left(end) ? end+2 : end-2;
    }
    return NULL;
}

End *end_getOtherBlockEnd(End *end) {
    if (!end_partOfBlock(end)) {
        return NULL; //the end must be block end to return the other end of a block!
    }
    assert(!end_isBlock(end));
    return end_left(end) ? end+4 : end-4;
}

Group *end_getGroup(End *end) {
    if(end_partOfBlock(end)) {
        assert(!end_isBlock(end));
        BlockEndContents *i = end_getBlockEndContents(end);
        return end_left(end) ? i->leftGroup : i->rightGroup;
    }
    return end_getContents(end)->group;
}

bool end_isEmpty(End *end) {
    return end_getFirst(end) == NULL;
}

int64_t end_getInstanceNumber(End *end) {
    assert(!end_isBlock(end));
    if(end_partOfBlock(end)) {
        assert(!end_isBlock(end));
        return block_getInstanceNumber(end_getBlock(end));
    }

    Cap *cap = end_getContents(end)->firstCap;
    int64_t totalCaps = 0;
    while(cap != NULL) {
        totalCaps++; cap = cap_getContents(cap)->nCap;
    }
    return totalCaps;
}

Cap *end_getInstance(End *end, Name name) {
    Flower *flower = end_getFlower(end);
    assert(flower != NULL);
    Cap *cap = flower_getCap(flower, name);
    if(cap != NULL) {
        assert(cap_getName(cap) == name);
        assert(cap == cap_getPositiveOrientation(cap));
        cap = end_getOrientation(end) ? cap : cap_getReverse(cap);
        if(cap_getEnd(cap) == end) {
            return cap;
        }
    }
    return NULL;
}

End_InstanceIterator *end_getInstanceIterator(End *end) {
    End_InstanceIterator *iterator;
    iterator = st_malloc(sizeof(struct _end_instanceIterator));
    iterator->end = end;
    iterator->cap = end_partOfBlock(end) ? block_getFirst(end_getBlock(end)) : end_getContents(end)->firstCap;
    return iterator;
}

Cap *end_getNext(End_InstanceIterator *iterator) {
    Cap *cap = iterator->cap;
    if(cap == NULL) {
        return NULL;
    }
    // If the we're iterating through a cap segment
    if(end_partOfBlock(iterator->end)) {
        iterator->cap = segment_getContents(cap)->nSegment;
        // Get the segment in the correction orientation
        Segment *segment = block_getOrientation(end_getBlock(iterator->end)) == segment_getOrientation(cap) ? cap : segment_getReverse(cap);
        // Get the cap for the desired side of the segment
        cap = end_getSide(iterator->end) ? segment_get5Cap(segment) : segment_get3Cap(segment);
        assert(cap_getEnd(cap) == iterator->end);
        return cap;
    }

    iterator->cap = cap_getContents(cap)->nCap;
    return end_getOrientation(iterator->end) ? cap : cap_getReverse(cap);
}

Cap *end_getFirst(End *end) {
    // If is it attached to a block
    if(end_partOfBlock(end)) {
        Block *b = end_getBlock(end);
        Segment *s = block_getFirst(b);
        Cap *c = end_getSide(end) ? segment_get5Cap(s) : segment_get3Cap(s);
        assert(cap_getEnd(c) == end);
        return c;
    }

    Cap *c = end_getContents(end)->firstCap;
    return end_getOrientation(end) ? c : cap_getReverse(c);
}

void end_destructInstanceIterator(End_InstanceIterator *iterator) {
    free(iterator);
}

void end_setGroup(End *end, Group *group) {
    if (end_getGroup(end) != NULL) {
        group_removeEnd(end_getGroup(end), end);
    }
    if(end_partOfBlock(end)) {
        assert(!end_isBlock(end));
        BlockEndContents *b = end_getBlockEndContents(end);
        if(end_left(end)) {
            b->leftGroup = group;
        }
        else {
            b->rightGroup = group;
        }
    }
    else {
        end_getContents(end)->group = group;
    }
    if (group != NULL) {
        group_addEnd(group, end);
    }
}

void end_check(End *end) {
    //Check is connected to flower properly
    assert(flower_getEnd(end_getFlower(end), end_getName(end)) == end_getPositiveOrientation(end));

    //check end is part of group..
    Group *group = end_getGroup(end);
    assert(group != NULL);
    assert(group_getEnd(group, end_getName(end)) == end_getPositiveOrientation(end));

    if (end_isBlockEnd(end)) {
        assert(!end_isStubEnd(end));
        assert(end_isFree(end));
        //Check block..
        Block *block = end_getBlock(end);
        assert(block != NULL);
        assert(block_getOrientation(block) == end_getOrientation(end));
        //check not attached
        assert(end_isFree(end));
        assert(!end_isAttached(end));
        //Check sides correspond..
        if (end_getSide(end)) {
            assert(block_get5End(block) == end);
        } else {
            assert(block_get3End(block) == end);
        }
    } else {
        assert(end_isStubEnd(end)); //Is stub end:
        //there must be no attached block.
        assert(end_getBlock(end) == NULL);
        Group *parentGroup = flower_getParentGroup(end_getFlower(end));
        if (parentGroup != NULL) {
            // if attached the is inherited from a parent flower to the containing flower.
            End *parentEnd = group_getEnd(parentGroup, end_getName(end));
            assert(end_getOrientation(parentEnd));
            if (end_isAttached(end)) {
                assert(parentEnd != NULL);
            }
            if (parentEnd != NULL) {
                assert(end_getSide(parentEnd) == end_getSide(end_getPositiveOrientation(end)));
            }
        }
    }

    //Check reverse, not comprehensively, perhaps.
    End *rEnd = end_getReverse(end);
    assert(rEnd != NULL);
    assert(end_getReverse(rEnd) == end);
    assert(end_getOrientation(end) == !end_getOrientation(rEnd));
    assert(end_getSide(end) == !end_getSide(rEnd));
    assert(end_getName(end) == end_getName(rEnd));
    assert(end_getInstanceNumber(end) == end_getInstanceNumber(rEnd));
    assert(end_isAttached(end) == end_isAttached(rEnd));
    assert(end_isStubEnd(end) == end_isStubEnd(rEnd));
    if (!end_isEmpty(end)) {
        assert(end_getFirst(end) == cap_getReverse(end_getFirst(rEnd)));
    }

    //For each segment calls segment_check.
    End_InstanceIterator *iterator = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(iterator)) != NULL) {
        cap_check(cap);
    }
    end_destructInstanceIterator(iterator);
}

Cap *end_getCapForEvent(End *end, Name eventName) {
    /*
     * Get the cap for a given event.
     */
    End_InstanceIterator *it = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(it)) != NULL) {
        if (event_getName(cap_getEvent(cap)) == eventName) {
            end_destructInstanceIterator(it);
            return cap;
        }
    }
    end_destructInstanceIterator(it);
    return NULL;
}

/*
 * Private functions
 */

void end_addInstance(End *end, Cap *cap) {
    if(!end_partOfBlock(end)) {
        assert(cap_getContents(cap)->nCap == NULL);
        assert(!cap_partOfSegment(cap));
        cap_getContents(cap)->nCap = end_getContents(end)->firstCap;
        end_getContents(end)->firstCap = cap_getPositiveOrientation(cap);
    }
    else {
        assert(!end_isBlock(end));
    }
}

void end_removeInstance(End *end, Cap *cap) {
    if(!end_partOfBlock(end)) {
        assert(!cap_partOfSegment(cap));
        Cap **capP = &(end_getContents(end)->firstCap);
        while (*capP != NULL) {
            if (cap_getName(cap) == cap_getName(*capP)) {
                (*capP) = cap_getContents(*capP)->nCap; // Splice it out
                return;
            }
            capP = &(cap_getContents(*capP)->nCap);
        }
    }
    else {
        assert(!end_isBlock(end));
    }
}

void end_setFlower(End *end, Flower *flower) {
    flower_removeEnd(end_getFlower(end), end);
    if(end_partOfBlock(end)) {
        assert(0); // todo: This needs fixing if this happens
        end_getBlockEndContents(end)->flower = flower;
    }
    else {
        end_getContents(end)->flower = flower;
    }
    flower_addEnd(flower, end);
}

uint64_t end_hashKey(const void *o) {
    return end_getName((End *) o);
}

int end_hashEqualsKey(const void *o, const void *o2) {
    End *end1 = (End *) o;
    End *end2 = (End *) o2;
    return end_getName(end1) == end_getName(end2) && end_getOrientation(end1)
                                                     == end_getOrientation(end2);
}
