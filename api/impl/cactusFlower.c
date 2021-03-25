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

/*static int flower_constructFacesP(const void *o1, const void *o2) {
    assert(o1 != NULL);
    assert(o2 != NULL);
    return o1 == o2 ? 0 : o1 > o2 ? 1 : -1;
}*/

static Flower *flower_construct3(Name name, CactusDisk *cactusDisk) {
    Flower *flower;
    flower = st_malloc(sizeof(Flower));

    flower->name = name;

    flower->sequences = stList_construct3(0, NULL); //stSortedSet_construct3(flower_constructSequencesP, NULL);
    flower->caps = stList_construct3(0, NULL); //stSortedSet_construct3(flower_constructCapsP, NULL);
    flower->ends = stList_construct3(0, NULL); //stSortedSet_construct3(flower_constructEndsP, NULL);
    //flower->segments = NULL; //stSortedSet_construct3(flower_constructSegmentsP, NULL);
    //flower->blocks = stSortedSet_construct3(flower_constructBlocksP, NULL);
    flower->groups = stList_construct3(0, NULL); // stSortedSet_construct3(flower_constructGroupsP, NULL);
    flower->chains = stList_construct3(0, NULL); //stSortedSet_construct3(flower_constructChainsP, NULL);
    //flower->faces = NULL; //stSortedSet_construct3(flower_constructFacesP, NULL);

    flower->parentFlowerName = NULL_NAME;
    flower->cactusDisk = cactusDisk;
    //flower->faceIndex = 0;
    //flower->chainIndex = 0;

    flower->builtBlocks = 0;
    //flower->builtFaces = 0;
    //flower->builtTrees = 0;

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

void flower_destruct(Flower *flower, int64_t recursive) {
    Flower_GroupIterator *iterator;
    Sequence *sequence;
    End *end;
    Group *group;
    Chain *chain;
    Flower *nestedFlower;

    if (recursive) {
        iterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(iterator)) != NULL) {
            nestedFlower = group_getNestedFlower(group);
            if (nestedFlower != NULL) {
                flower_destruct(nestedFlower, recursive);
            }
        }
        flower_destructGroupIterator(iterator);
    }

    cactusDisk_removeFlower(flower->cactusDisk, flower);

    while ((sequence = flower_getFirstSequence(flower)) != NULL) {
        flower_removeSequence(flower, sequence);
    }
    stList_destruct(flower->sequences);

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

    while ((chain = flower_getFirstChain(flower)) != NULL) {
        chain_destruct(chain);
    }
    stList_destruct(flower->chains);

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
    //Sequence sequence;
    Sequence sequence;
    //sequence.sequence = &sequence;
    sequence.name = name;
    //return cactusDisk_getSequence(flower->cactusDisk, name);
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
    return stList_length(flower->caps) > 0 ? stList_peek(flower->caps) : NULL; //stSortedSet_getFirst(flower->caps);
}

Cap *flower_getCap(Flower *flower, Name name) {
    //Cap cap;
    CapContents capContents[2];
    Cap *cap = (Cap *)(&capContents); // Very ugly hack
    cap->bits = 2; // binary: 000010
    cap_getCoreContents(cap)->instance = name;
    assert(cap_getName(cap) == name);
    return stList_binarySearch(flower->caps, cap, flower_constructCapsP); //stSortedSet_search(flower->caps, cap);
}

int64_t flower_getCapNumber(Flower *flower) {
    return stList_length(flower->caps); //stSortedSet_size(flower->caps);
}

Flower_CapIterator *flower_getCapIterator(Flower *flower) {
    return stList_getIterator(flower->caps); //stSortedSet_getIterator(flower->caps);
}

Cap *flower_getNextCap(Flower_CapIterator *capIterator) {
    return stList_getNext(capIterator); //stSortedSet_getNext(capIterator);
}

Cap *flower_getPreviousCap(Flower_CapIterator *capIterator) {
    return stList_getPrevious(capIterator); //stSortedSet_getPrevious(capIterator);
}

Flower_CapIterator *flower_copyCapIterator(Flower_CapIterator *capIterator) {
    return stList_copyIterator(capIterator); //stSortedSet_copyIterator(capIterator);
}

void flower_destructCapIterator(Flower_CapIterator *capIterator) {
    stList_destructIterator(capIterator); //stSortedSet_destructIterator(capIterator);
}

End *flower_getFirstEnd(Flower *flower) {
    return stList_length(flower->ends) > 0 ? stList_peek(flower->ends) : NULL;
}

End *flower_getEnd(Flower *flower, Name name) {
    //fprintf(stderr, "Starting flower get end\n");
    EndContents endContents[2];
    End *end = (End *)(&endContents); // Very ugly hack
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
    return stList_length(flower->groups) > 0 ? stList_peek(flower->groups): NULL; //stSortedSet_getFirst(flower->groups);
}

Group *flower_getGroup(Flower *flower, Name flowerName) {
    Group group;
    group.name = flowerName;
    return stList_binarySearch(flower->groups, &group, flower_constructGroupsP); //stSortedSet_search(flower->groups, &group);
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
    return stList_length(flower->chains) > 0 ? stList_peek(flower->chains) : NULL; //stSortedSet_getFirst(flower->chains);
}

Chain *flower_getChain(Flower *flower, Name name) {
    Chain chain;
    chain.name = name;
    return stList_binarySearch(flower->chains, &chain, flower_constructChainsP); //stSortedSet_search(flower->chains, &chain);
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

Face *flower_getFirstFace(Flower *flower) {
    assert(0);
    return NULL; //stSortedSet_getFirst(flower->faces);
}

int64_t flower_getFaceNumber(Flower *flower) {
    assert(0);
    return 0; //stSortedSet_size(flower->faces);
}

Flower_FaceIterator *flower_getFaceIterator(Flower *flower) {
    assert(0);
    return NULL;
    //return stSortedSet_getIterator(flower->faces);
}

Face *flower_getNextFace(Flower_FaceIterator *faceIterator) {
    assert(0);
    return NULL;
    //return stSortedSet_getNext(faceIterator);
}

Face *flower_getPreviousFace(Flower_FaceIterator *faceIterator) {
    assert(0);
    return NULL;
    //return stSortedSet_getPrevious(faceIterator);
}

Flower_FaceIterator *flower_copyFaceIterator(Flower_FaceIterator *faceIterator) {
    assert(0);
    return NULL;
    //return stSortedSet_copyIterator(faceIterator);
}

void flower_destructFaceIterator(Flower_FaceIterator *faceIterator) {
    assert(0);
    //stSortedSet_destructIterator(faceIterator);
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

    /*if (flower_builtFaces(flower)) {
        Flower_FaceIterator *faceIterator = flower_getFaceIterator(flower);
        Face *face;
        while ((face = flower_getNextFace(faceIterator)) != NULL) {
            face_check(face);
        }
        flower_destructFaceIterator(faceIterator);
        face_checkFaces(flower);
    } else {
        cactusCheck(flower_getFaceNumber(flower) == 0);
    }*/

    if (!flower_builtBlocks(flower)) {
        cactusCheck(flower_isLeaf(flower)); //Defensive
        cactusCheck(flower_isTerminal(flower)); //Checks that a flower without built blocks is a leaf and does not
        //contain any blocks.
    }

    /*Flower_SequenceIterator *sequenceIterator = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(sequenceIterator)) != NULL) {
        sequence_check(sequence);
    }
    flower_destructSequenceIterator(sequenceIterator);*/
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

bool flower_builtTrees(Flower *flower) {
    return 0; //flower->builtTrees;
}

void flower_setBuiltTrees(Flower *flower, bool b) {
    assert(0);
    //flower->builtTrees = b;
}

bool flower_builtFaces(Flower *flower) {
    return 0; //flower->builtFaces;
}

void flower_setBuildFaces(Flower *flower, bool b) {
    assert(0);
    //flower->builtFaces = b;
    //if (flower_builtFaces(flower)) {
    //    flower_reconstructFaces(flower);
    //}
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

bool flower_removeIfRedundant(Flower *flower) {
    assert(0);
    return 0;
}

void flower_delete2(Flower *flower, bool isOnDisk) {
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (!group_isLeaf(group)) {
            flower_delete2(group_getNestedFlower(group), isOnDisk);
        }
    }
    flower_destructGroupIterator(groupIt);
    Group *parentGroup = flower_getParentGroup(flower);
    if(parentGroup != NULL) {
        group_setLeaf(parentGroup, 1);
    }
    //This needs modification so that we don't do this directly..
    if(isOnDisk) {
        cactusDisk_deleteFlowerFromDisk(flower_getCactusDisk(flower), flower);
    }
    flower_destruct(flower, 0);
}

void flower_delete(Flower *flower) {
    flower_delete2(flower, 1);
}

bool flower_deleteIfEmpty(Flower *flower) {
    if (flower_getEndNumber(flower) == 0 && flower_getParentGroup(flower) != NULL) { //contains nothing useful..
        assert(flower_getChainNumber(flower) == 0);
        while (flower_getGroupNumber(flower) > 0) {
            Group *group = flower_getFirstGroup(flower);
            if (!group_isLeaf(group)) {
                bool i = flower_deleteIfEmpty(group_getNestedFlower(group));
                (void) i;
                assert(i);
            }
        }
        assert(flower_getGroupNumber(flower) == 0);
        //This needs modification so that we don't do this directly..
        cactusDisk_deleteFlowerFromDisk(flower_getCactusDisk(flower), flower);
        Group *parentGroup = flower_getParentGroup(flower);
        group_destruct(parentGroup);
        flower_destruct(flower, 0);
        return 1;
    }
    return 0;
}

void flower_makeTerminalNormal(Flower *flower) {
    if (!flower_isTerminal(flower)) {
        Flower_GroupIterator *groupIterator;
        Group *group;
        groupIterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if (group_isLeaf(group)) {
                //assert(group_getTotalBaseLength(group) == 0);
                Flower *nestedFlower = group_makeNestedFlower(group);
                flower_setBuiltBlocks(nestedFlower, flower_builtBlocks(flower));
                //flower_setBuiltTrees(nestedFlower, flower_builtTrees(flower));
                //flower_setBuildFaces(nestedFlower, flower_builtFaces(flower));
            }
        }
        flower_destructGroupIterator(groupIterator);
    }
}

void flower_unload(Flower *flower) {
    flower_destruct(flower, 0);
}

bool flower_isParentLoaded(Flower *flower) {
    Name parentName = flower->parentFlowerName;
    if (parentName != NULL_NAME) {
        CactusDisk *cactusDisk = flower_getCactusDisk(flower);
        return cactusDisk_flowerIsLoaded(cactusDisk, parentName);
    }
    return 0;
}

void flower_unloadParent(Flower *flower) {
    Name parentName = flower->parentFlowerName;
    if (parentName != NULL_NAME) {
        CactusDisk *cactusDisk = flower_getCactusDisk(flower);
        if (cactusDisk_flowerIsLoaded(cactusDisk, parentName)) {
            Flower *parentFlower = cactusDisk_getFlower(cactusDisk, parentName);
            flower_unload(parentFlower);
        }
    }
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
    //assert(stSortedSet_search(flower->sequences, sequence) == NULL);
    //stSortedSet_insert(flower->sequences, sequence);
}

void flower_removeSequence(Flower *flower, Sequence *sequence) {
    removeFromFlower(flower->sequences, sequence);
}

void flower_addCap(Flower *flower, Cap *cap) {
    cap = cap_getPositiveOrientation(cap);
    stList_append(flower->caps, cap);
    // Now ensure we have fixed the sort
    int64_t i = stList_length(flower->caps)-1;
    while(--i >= 0 && cap_getName(stList_get(flower->caps, i)) > cap_getName(cap)) {
        swap(flower->caps, i);
    }
}

void flower_addEnd(Flower *flower, End *end) {
    end = end_getPositiveOrientation(end);
    stList_append(flower->ends, end);
    // Now ensure we have fixed the sort
    int64_t i = stList_length(flower->ends)-1;
    while(--i >= 0 && end_getName(stList_get(flower->ends, i)) > end_getName(end)) {
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
    //assert(stSortedSet_search(flower->chains, chain) == NULL);
    //stSortedSet_insert(flower->chains, chain);
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
    //assert(stSortedSet_search(flower->groups, group) == NULL);
    //stSortedSet_insert(flower->groups, group);
}

void flower_removeGroup(Flower *flower, Group *group) {
    removeFromFlower(flower->groups, group);
}

void flower_setParentGroup(Flower *flower, Group *group) {
    //assert(flower->parentFlowerName == NULL_NAME); we can change this if merging the parent flowers, so this no longer applies.
    flower->parentFlowerName = flower_getName(group_getFlower(group));
}

void flower_addFace(Flower *flower, Face *face) {
    assert(0);
    //assert(stSortedSet_search(flower->faces, face) == NULL);
    //stSortedSet_insert(flower->faces, face);
}

void flower_removeFace(Flower *flower, Face *face) {
    assert(0);
    //assert(stSortedSet_search(flower->faces, face) != NULL);
    //stSortedSet_remove(flower->faces, face);
}

/*
 * Serialisation functions.
 */

void flower_writeBinaryRepresentation(Flower *flower, void(*writeFn)(const void * ptr, size_t size, size_t count)) {
    /*Flower_SequenceIterator *sequenceIterator;
    Flower_EndIterator *endIterator;
    Flower_BlockIterator *blockIterator;
    Flower_GroupIterator *groupIterator;
    Flower_ChainIterator *chainIterator;
    Sequence *sequence;
    End *end;
    Block *block;
    Group *group;
    Chain *chain;

    binaryRepresentation_writeElementType(CODE_FLOWER, writeFn);
    binaryRepresentation_writeName(flower_getName(flower), writeFn);
    binaryRepresentation_writeBool(flower_builtBlocks(flower), writeFn);
    //binaryRepresentation_writeBool(flower_builtTrees(flower), writeFn);
    //binaryRepresentation_writeBool(flower_builtFaces(flower), writeFn);
    binaryRepresentation_writeName(flower->parentFlowerName, writeFn);

    sequenceIterator = flower_getSequenceIterator(flower);
    while ((sequence = flower_getNextSequence(sequenceIterator)) != NULL) {
        sequence_writeBinaryRepresentation(sequence, writeFn);
    }
    flower_destructSequenceIterator(sequenceIterator);

    endIterator = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        end_writeBinaryRepresentation(end, writeFn);
    }
    flower_destructEndIterator(endIterator);

    blockIterator = flower_getBlockIterator(flower);
    while ((block = flower_getNextBlock(blockIterator)) != NULL) {
        block_writeBinaryRepresentation(block, writeFn);
    }
    flower_destructBlockIterator(blockIterator);

    groupIterator = flower_getGroupIterator(flower);
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        group_writeBinaryRepresentation(group, writeFn);
    }
    flower_destructGroupIterator(groupIterator);

    chainIterator = flower_getChainIterator(flower);
    while ((chain = flower_getNextChain(chainIterator)) != NULL) {
        chain_writeBinaryRepresentation(chain, writeFn);
    }
    flower_destructChainIterator(chainIterator);

    binaryRepresentation_writeElementType(CODE_FLOWER, writeFn); //this avoids interpretting things wrong.*/
}

Flower *flower_loadFromBinaryRepresentation(void **binaryString, CactusDisk *cactusDisk) {
    Flower *flower = NULL;
    //bool buildFaces;
    /*if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_FLOWER) {
        binaryRepresentation_popNextElementType(binaryString);
        flower = flower_construct3(binaryRepresentation_getName(binaryString), cactusDisk);
        flower_setBuiltBlocks(flower, binaryRepresentation_getBool(binaryString));
        //flower_setBuiltTrees(flower, binaryRepresentation_getBool(binaryString));
        //buildFaces = binaryRepresentation_getBool(binaryString);
        flower->parentFlowerName = binaryRepresentation_getName(binaryString);
        while (sequence_loadFromBinaryRepresentation(binaryString, flower) != NULL)
            ;
        while (end_loadFromBinaryRepresentation(binaryString, flower) != NULL)
            ;
        while (block_loadFromBinaryRepresentation(binaryString, flower) != NULL)
            ;
        while (group_loadFromBinaryRepresentation(binaryString, flower) != NULL)
            ;
        while (chain_loadFromBinaryRepresentation(binaryString, flower) != NULL)
            ;
        //flower_setBuildFaces(flower, buildFaces);
        assert(binaryRepresentation_popNextElementType(binaryString) == CODE_FLOWER);
    }*/
    return flower;
}
