/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic flower functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static int flower_constructSequencesP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(sequence_getName((Sequence *) o1), sequence_getName((Sequence *) o2));
}

static int flower_constructCapsP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(cap_getName((Cap *) o1), cap_getName((Cap *) o2));
}

static int flower_constructEndsP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(end_getName((End *) o1), end_getName((End *) o2));
}

static int flower_constructGroupsP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(group_getName((Group *) o1), group_getName((Group *) o2));
}

static int flower_constructChainsP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(chain_getName((Chain *) o1), chain_getName((Chain *) o2));
}

static Flower *flower_construct3(Name name, CactusDisk *cactusDisk) {
    Flower *flower;
    flower = st_malloc(sizeof(Flower));
    flower->name = name;
    flower->sequences = stList_construct3(0, NULL);
    flower->caps = stList_construct3(0, NULL);
    flower->caps2 = NULL;
    flower->ends = stList_construct3(0, NULL);
    flower->groups = stList_construct3(0, NULL);
    flower->chains = stList_construct3(0, NULL);
    flower->parentFlowerName = NULL_NAME;
    flower->cactusDisk = cactusDisk;
    flower->builtBlocks = 0;
    cactusDisk_addFlower(flower->cactusDisk, flower);

    return flower;
}

Flower *flower_construct2(Name name, CactusDisk *cactusDisk) {
    Flower *flower = flower_construct3(name, cactusDisk);
    return flower;
}

Flower *flower_construct(CactusDisk *cactusDisk) {
    return flower_construct2(cactusDisk_getUniqueID(cactusDisk), cactusDisk);
}

void flower_destruct(Flower *flower, int64_t recursive, bool removeFromParentGroup) {
    Flower_GroupIterator *iterator;
    Sequence *sequence;
    End *end;
    Group *group;
    Chain *chain;
    Flower *nestedFlower;

    if(removeFromParentGroup) {
        Group *parentGroup = flower_getParentGroup(flower);
        if (parentGroup != NULL) {
            group_setLeaf(parentGroup, 1);
        }
    }

    if (recursive) {
        iterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(iterator)) != NULL) {
            nestedFlower = group_getNestedFlower(group);
            if (nestedFlower != NULL) {
                flower_destruct(nestedFlower, recursive, 0);
            }
        }
        flower_destructGroupIterator(iterator);
    }

    cactusDisk_removeFlower(flower->cactusDisk, flower);

    while ((sequence = flower_getFirstSequence(flower)) != NULL) {
        flower_removeSequence(flower, sequence);
    }
    stList_destruct(flower->sequences);

    while ((chain = flower_getFirstChain(flower)) != NULL) {
        assert(chain_getFlower(chain) == flower);
        chain_destruct(chain);
    }
    stList_destruct(flower->chains);

    while ((group = flower_getFirstGroup(flower)) != NULL) {
        assert(group_getFlower(group) == flower);
        group_destruct(group);
    }
    stList_destruct(flower->groups);

    // This cleans up the ends, blocks, caps and segments contained in the flower
    while ((end = flower_getFirstEnd(flower)) != NULL) {
        end_destruct(end);
    }
    stList_destruct(flower->caps);
    stList_destruct(flower->ends);

    free(flower);
}

Name flower_getName(Flower *flower) {
    return flower->name;
}

CactusDisk *flower_getCactusDisk(Flower *flower) {
    return flower->cactusDisk;
}

EventTree *flower_getEventTree(Flower *flower) {
    return cactusDisk_getEventTree(flower->cactusDisk);
}

Sequence *flower_getFirstSequence(Flower *flower) {
    return stList_length(flower->sequences) > 0 ? stList_peek(flower->sequences) : NULL; //stSortedSet_getFirst(flower->sequences);
}

Sequence *flower_getSequence(Flower *flower, Name name) {
    Sequence sequence;
    sequence.name = name;
    return stList_binarySearch(flower->sequences, &sequence, flower_constructSequencesP);
}

int64_t flower_getSequenceNumber(Flower *flower) {
    return stList_length(flower->sequences);
}

Flower_SequenceIterator *flower_getSequenceIterator(Flower *flower) {
    return stList_getIterator(flower->sequences);
}

Sequence *flower_getNextSequence(Flower_SequenceIterator *sequenceIterator) {
    return stList_getNext(sequenceIterator);
}

Sequence *flower_getPreviousSequence(Flower_SequenceIterator *sequenceIterator) {
    return stList_getPrevious(sequenceIterator);
}

Flower_SequenceIterator *flower_copySequenceIterator(Flower_SequenceIterator *sequenceIterator) {
    return stList_copyIterator(sequenceIterator);
}

void flower_destructSequenceIterator(Flower_SequenceIterator *sequenceIterator) {
    stList_destructIterator(sequenceIterator);
}

Cap *flower_getFirstCap(Flower *flower) {
    return stList_length(flower->caps) > 0 ? stList_peek(flower->caps) : NULL;
}

static int sort_caps(const void *a, const void *b) {
    return cactusMisc_nameCompare(cap_getName((Cap*)a), cap_getName((Cap*)b));
}

Cap *flower_getCap(Flower *flower, Name name) {
    CapContents capContents[2];
    Cap *cap = (Cap *)(&capContents); // Very ugly cast
    cap->bits = 2; // binary: 000010
    cap_getCoreContents(cap)->instance = name;
    assert(cap_getName(cap) == name);
    Cap *result = stList_binarySearch(flower->caps, cap, flower_constructCapsP);    
    if (result == NULL && flower->caps2 != NULL) {
        // fall back on unsorted caps2 list
        int64_t caps2_size = stList_length(flower->caps2);
        for (int64_t i = 0; i < caps2_size && result == NULL; ++i) {
            Cap *cur_cap = stList_get(flower->caps2, i);
            if (sort_caps(cap, cur_cap) == 0) {
                result = cur_cap;
            }
        }
        // some defragmentation
        if (caps2_size > MAX_FLOWER_LAZY_CAPS_SIZE) {
            stList_appendAll(flower->caps, flower->caps2);
            stList_sort(flower->caps, sort_caps);
            stList_destruct(flower->caps2);
            flower->caps2 = stList_construct3(0, NULL);            
        }
    }
    return result;
}

int64_t flower_getCapNumber(Flower *flower) {
    return stList_length(flower->caps);
}

Flower_CapIterator *flower_getCapIterator(Flower *flower) {
    return stList_getIterator(flower->caps);
}

Cap *flower_getNextCap(Flower_CapIterator *capIterator) {
    return stList_getNext(capIterator);
}

Cap *flower_getPreviousCap(Flower_CapIterator *capIterator) {
    return stList_getPrevious(capIterator);
}

Flower_CapIterator *flower_copyCapIterator(Flower_CapIterator *capIterator) {
    return stList_copyIterator(capIterator);
}

void flower_destructCapIterator(Flower_CapIterator *capIterator) {
    stList_destructIterator(capIterator);
}

End *flower_getFirstEnd(Flower *flower) {
    return stList_length(flower->ends) > 0 ? stList_peek(flower->ends) : NULL;
}

End *flower_getEnd(Flower *flower, Name name) {
    EndContents endContents[2];
    End *end = (End *)(&endContents); // Very ugly cast
    end->bits = 0x1;
    end_getContents(end)->name = name;
    assert(end_getName(end) == name);
    return stList_binarySearch(flower->ends, end, flower_constructEndsP);
}

Block *flower_getBlock(Flower *flower, Name name) {
    End *end = flower_getEnd(flower, name-1);
    if(end == NULL) {
        return NULL;
    }
    assert(end_partOfBlock(end));
    assert(end_left(end));
    return end_getBlock(end);
}

int64_t flower_getEndNumber(Flower *flower) {
    return stList_length(flower->ends);
}

int64_t flower_getBlockEndNumber(Flower *flower) {
    End *end;
    Flower_EndIterator *iterator = flower_getEndIterator(flower);
    int64_t i = 0;
    while ((end = flower_getNextEnd(iterator)) != NULL) {
        if (end_partOfBlock(end)) {
            i++;
        }
    }
    flower_destructEndIterator(iterator);
    assert(i % 2 == 0);
    return i;
}

int64_t flower_getBlockNumber(Flower *flower) {
    return flower_getBlockEndNumber(flower)/2;
}

int64_t flower_getStubEndNumber(Flower *flower) {
    return flower_getEndNumber(flower) - flower_getBlockEndNumber(flower);
}

int64_t flower_getFreeStubEndNumber(Flower *flower) {
    End *end;
    Flower_EndIterator *iterator = flower_getEndIterator(flower);
    int64_t i = 0;
    while ((end = flower_getNextEnd(iterator)) != NULL) {
        if (end_isStubEnd(end) && end_isFree(end)) {
            i++;
        }
    }
    flower_destructEndIterator(iterator);
    return i;
}

int64_t flower_getAttachedStubEndNumber(Flower *flower) {
    return flower_getStubEndNumber(flower) - flower_getFreeStubEndNumber(flower);
}

Flower_EndIterator *flower_getEndIterator(Flower *flower) {
    return stList_getIterator(flower->ends);
}

End *flower_getNextEnd(Flower_EndIterator *endIterator) {
    return stList_getNext(endIterator);
}

End *flower_getPreviousEnd(Flower_EndIterator *endIterator) {
    return stList_getPrevious(endIterator);
}

Flower_EndIterator *flower_copyEndIterator(Flower_EndIterator *endIterator) {
    return stList_copyIterator(endIterator);
}

void flower_destructEndIterator(Flower_EndIterator *endIterator) {
    stList_destructIterator(endIterator);
}

Group *flower_getFirstGroup(Flower *flower) {
    return stList_length(flower->groups) > 0 ? stList_peek(flower->groups): NULL;
}

Group *flower_getGroup(Flower *flower, Name flowerName) {
    Group group;
    group.name = flowerName;
    return stList_binarySearch(flower->groups, &group, flower_constructGroupsP);
}

int64_t flower_getGroupNumber(Flower *flower) {
    return stList_length(flower->groups);
}

Flower_GroupIterator *flower_getGroupIterator(Flower *flower) {
    return stList_getIterator(flower->groups);
}

Group *flower_getNextGroup(Flower_GroupIterator *groupIterator) {
    return stList_getNext(groupIterator);
}

Group *flower_getPreviousGroup(Flower_GroupIterator *groupIterator) {
    return stList_getPrevious(groupIterator);
}

Flower_GroupIterator *flower_copyGroupIterator(Flower_GroupIterator *groupIterator) {
    return stList_copyIterator(groupIterator);
}

void flower_destructGroupIterator(Flower_GroupIterator *groupIterator) {
    stList_destructIterator(groupIterator);
}

bool flower_hasParentGroup(Flower *flower) {
    return flower->parentFlowerName != NULL_NAME;
}

Group *flower_getParentGroup(Flower *flower) {
    if (flower->parentFlowerName == NULL_NAME) {
        return NULL;
    }
    Flower *flower2 = cactusDisk_getFlower(flower_getCactusDisk(flower), flower->parentFlowerName);
    assert(flower2 != NULL);
    return flower_getGroup(flower2, flower_getName(flower));
}

Chain *flower_getFirstChain(Flower *flower) {
    return stList_length(flower->chains) > 0 ? stList_peek(flower->chains) : NULL;
}

Chain *flower_getChain(Flower *flower, Name name) {
    Chain chain;
    chain.name = name;
    return stList_binarySearch(flower->chains, &chain, flower_constructChainsP);
}

int64_t flower_getChainNumber(Flower *flower) {
    return stList_length(flower->chains);
}

int64_t flower_getTrivialChainNumber(Flower *flower) {
    int64_t i = 0;
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if(end_partOfBlock(end) && end_left(end)) {
            if (block_isTrivialChain(end_getBlock(end))) {
                i++;
            }
        }
    }
    flower_destructEndIterator(endIt);
    return i;
}

Flower_ChainIterator *flower_getChainIterator(Flower *flower) {
    return stList_getIterator(flower->chains);
}

Chain *flower_getNextChain(Flower_ChainIterator *chainIterator) {
    return stList_getNext(chainIterator);
}

Chain *flower_getPreviousChain(Flower_ChainIterator *chainIterator) {
    return stList_getPrevious(chainIterator);
}

Flower_ChainIterator *flower_copyChainIterator(Flower_ChainIterator *chainIterator) {
    return stList_copyIterator(chainIterator);
}

void flower_destructChainIterator(Flower_ChainIterator *chainIterator) {
    stList_destructIterator(chainIterator);
}

int64_t flower_getTotalBaseLength(Flower *flower) {
    /*
     * The implementation of this function is very like that in group_getTotalBaseLength, with a few differences. Consider merging them.
     */
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    int64_t totalLength = 0;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        if (!end_isBlockEnd(end)) {
            End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
            Cap *cap;
            while ((cap = end_getNext(instanceIterator)) != NULL) {
                cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                if (!cap_getSide(cap) && cap_getSequence(cap) != NULL) {
                    Cap *cap2 = cap_getAdjacency(cap);
                    assert(cap2 != NULL);
                    while (end_isBlockEnd(cap_getEnd(cap2))) {
                        Segment *segment = cap_getSegment(cap2);
                        assert(segment != NULL);
                        assert(segment_get5Cap(segment) == cap2);
                        cap2 = cap_getAdjacency(segment_get3Cap(segment));
                        assert(cap2 != NULL);
                        assert(cap_getStrand(cap2));
                        assert(cap_getSide(cap2));
                    }
                    assert(cap_getStrand(cap2));
                    assert(cap_getSide(cap2));
                    int64_t length = cap_getCoordinate(cap2) - cap_getCoordinate(cap) - 1;
                    assert(length >= 0);
                    totalLength += length;
                }
            }
            end_destructInstanceIterator(instanceIterator);
        }
    }
    flower_destructEndIterator(endIterator);
    return totalLength;
}

void flower_checkNotEmpty(Flower *flower, bool recursive) {
    //First check the flower is not empty, unless it is the parent group.
    if (flower_hasParentGroup(flower)) {
        assert(flower_getGroupNumber(flower) > 0);
        assert(flower_getEndNumber(flower) > 0);
        assert(flower_getAttachedStubEndNumber(flower) > 0); //We must have some ends to tie us into the parent problem + flower_getBlockEndNumber(flower) > 0);
    }
    //Now Checks that each group contains at least one end and call recursive.
    Group *group;
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        assert(!group_isEmpty(group));
        assert(group_getAttachedStubEndNumber(group) + group_getBlockEndNumber(group) > 0);
        if (recursive && !group_isLeaf(group)) {
            flower_checkNotEmpty(group_getNestedFlower(group), 1);
        }
    }
    flower_destructGroupIterator(groupIt);
}

void flower_check(Flower *flower) {
    eventTree_check(flower_getEventTree(flower));

    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        group_check(group);
    }
    flower_destructGroupIterator(groupIterator);

    Flower_ChainIterator *chainIterator = flower_getChainIterator(flower);
    Chain *chain;
    while ((chain = flower_getNextChain(chainIterator)) != NULL) {
        chain_check(chain);
    }
    flower_destructChainIterator(chainIterator);

    //We check built trees in here.
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        end_check(end);
        end_check(end_getReverse(end)); //We will test everything backwards also.

        // Check any attached block
        if(end_partOfBlock(end) && end_left(end)) {
            assert(flower_builtBlocks(flower));
            block_check(end_getBlock(end));
            block_check(block_getReverse(end_getBlock(end))); //We will test everything backwards also.
        }
    }
    flower_destructEndIterator(endIterator);

    if (!flower_builtBlocks(flower)) {
        cactusCheck(flower_isLeaf(flower)); //Defensive
        cactusCheck(flower_isTerminal(flower)); //Checks that a flower without built blocks is a leaf and does not
        //contain any blocks.
    }
}

void flower_checkRecursive(Flower *flower) {
    flower_check(flower);
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (!group_isLeaf(group)) {
            flower_checkRecursive(group_getNestedFlower(group));
        }
    }
    flower_destructGroupIterator(groupIt);
}

bool flower_builtBlocks(Flower *flower) {
    return flower->builtBlocks;
}

void flower_setBuiltBlocks(Flower *flower, bool b) {
    flower->builtBlocks = b;
}

bool flower_isLeaf(Flower *flower) {
    Group *group;
    Flower_GroupIterator *iterator = flower_getGroupIterator(flower);
    while ((group = flower_getNextGroup(iterator)) != NULL) {
        if (!group_isLeaf(group)) {
            flower_destructGroupIterator(iterator);
            return 0;
        }
    }
    flower_destructGroupIterator(iterator);
    return 1;
}

bool flower_isTerminal(Flower *flower) {
    return flower_isLeaf(flower) && flower_getGroupNumber(flower) <= 1 && flower_getStubEndNumber(flower) == flower_getEndNumber(flower);
}

/*
 * Private functions
 */

static void removeFromFlower(stList *l, void *item) {
    assert(stList_length(l) > 0);
    if(stList_peek(l) == item) {
        stList_pop(l);
    }
    else {
        assert(stList_contains(l, item));
        stList_removeItem(l, item);
    }
}

static void swap(stList *l, int64_t i) {
    void *o = stList_get(l, i+1);
    stList_set(l, i+1, stList_get(l, i));
    stList_set(l, i, o);
}

void flower_addSequence(Flower *flower, Sequence *sequence) {
    stList_append(flower->sequences, sequence);
    // Now ensure we have fixed the sort
    int64_t i = stList_length(flower->sequences)-1;
    while(--i >= 0 && sequence_getName(stList_get(flower->sequences, i)) > sequence_getName(sequence)) {
        swap(flower->sequences, i);
    }
}

void flower_removeSequence(Flower *flower, Sequence *sequence) {
    removeFromFlower(flower->sequences, sequence);
}

void flower_setLazyCaps(Flower *flower, bool b) {
    if (b == true) {
        // activate the second list, and keep a copy of the first list in caps2
        assert(flower->caps2 == NULL);
        flower->caps2 = stList_construct3(0, NULL);
    } else {
        assert(flower->caps2 != NULL);
        // merge up the lists
        stList_appendAll(flower->caps, flower->caps2);
        stList_destruct(flower->caps2);
        flower->caps2 = NULL;
        stList_sort(flower->caps, sort_caps);
    }
}


void flower_bulkAddCaps(Flower *flower, stList *capsToAdd) {
    if(stList_length(capsToAdd) > 0) {
        stList_appendAll(flower->caps, capsToAdd);
        stList_sort(flower->caps, sort_caps);
    }
}

void flower_addCap(Flower *flower, Cap *cap) {
    cap = cap_getPositiveOrientation(cap);
    if (flower->caps2 != NULL) {
        // if we're in lazy mode, add to unsorted list
        stList_append(flower->caps2, cap);
    } else {        
        stList_append(flower->caps, cap);
        // Now ensure we have fixed the sort
        int64_t i = stList_length(flower->caps)-1;
        while(--i >= 0 && cap_getName(stList_get(flower->caps, i)) > cap_getName(cap)) {
            swap(flower->caps, i);
        }
    }
}

static int sort_ends(const void *a, const void *b) {
    return cactusMisc_nameCompare(end_getName((End*)a), end_getName((End*)b));
}

void flower_bulkAddEnds(Flower *flower, stList *endsToAdd) {
    if(stList_length(endsToAdd) > 0) {
        stList_appendAll(flower->ends, endsToAdd);
        stList_sort(flower->ends, sort_ends);
    }
}

void flower_addEnd(Flower *flower, End *end) {
    end = end_getPositiveOrientation(end);
    stList_append(flower->ends, end);
    // Now ensure we have fixed the sort
    int64_t i = stList_length(flower->ends)-1;
    while(--i >= 0 && end_getName(stList_get(flower->ends, i)) > end_getName(end)) {
        //fprintf(stderr, "adding a end out of order: %" PRIi64 " flower name: %" PRIi64 "\n", end_getName(end), flower_getName(flower));
        swap(flower->ends, i);
    }
}

void flower_removeEnd(Flower *flower, End *end) {
    removeFromFlower(flower->ends, end);
}

void flower_addChain(Flower *flower, Chain *chain) {
    stList_append(flower->chains, chain);
    // Now ensure we have fixed the sort
    int64_t i = stList_length(flower->chains)-1;
    while(--i >= 0 && chain_getName(stList_get(flower->chains, i)) > chain_getName(chain)) {
        swap(flower->chains, i);
    }
}

void flower_removeChain(Flower *flower, Chain *chain) {
    removeFromFlower(flower->chains, chain);
}

void flower_addGroup(Flower *flower, Group *group) {
    stList_append(flower->groups, group);
    // Now ensure we have fixed the sort
    int64_t i = stList_length(flower->groups)-1;
    while(--i >= 0 && group_getName(stList_get(flower->groups, i)) > group_getName(group)) {
        swap(flower->groups, i);
    }
}

void flower_removeGroup(Flower *flower, Group *group) {
    removeFromFlower(flower->groups, group);
}

void flower_setParentGroup(Flower *flower, Group *group) {
    flower->parentFlowerName = flower_getName(group_getFlower(group));
}

