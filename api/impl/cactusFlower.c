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

static int flower_constructSegmentsP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(segment_getName((Segment *) o1), segment_getName((Segment *) o2));
}

static int flower_constructBlocksP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(block_getName((Block *) o1), block_getName((Block *) o2));
}

static int flower_constructGroupsP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(group_getName((Group *) o1), group_getName((Group *) o2));
}

static int flower_constructChainsP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(chain_getName((Chain *) o1), chain_getName((Chain *) o2));
}

static int flower_constructFacesP(const void *o1, const void *o2) {
    assert(o1 != NULL);
    assert(o2 != NULL);
    return o1 == o2 ? 0 : o1 > o2 ? 1 : -1;
}

static Flower *flower_construct3(Name name, CactusDisk *cactusDisk) {
    Flower *flower;
    flower = st_malloc(sizeof(Flower));

    flower->name = name;

    flower->sequences = stSortedSet_construct3(flower_constructSequencesP, NULL);
    flower->caps = stSortedSet_construct3(flower_constructCapsP, NULL);
    flower->ends = stSortedSet_construct3(flower_constructEndsP, NULL);
    flower->segments = stSortedSet_construct3(flower_constructSegmentsP, NULL);
    flower->blocks = stSortedSet_construct3(flower_constructBlocksP, NULL);
    flower->groups = stSortedSet_construct3(flower_constructGroupsP, NULL);
    flower->chains = stSortedSet_construct3(flower_constructChainsP, NULL);
    flower->faces = stSortedSet_construct3(flower_constructFacesP, NULL);

    flower->parentFlowerName = NULL_NAME;
    flower->cactusDisk = cactusDisk;
    flower->faceIndex = 0;
    flower->chainIndex = 0;

    flower->builtBlocks = 0;
    flower->builtFaces = 0;
    flower->builtTrees = 0;

    cactusDisk_addFlower(flower->cactusDisk, flower);

    flower->eventTree = NULL;

    return flower;
}

Flower *flower_construct2(Name name, CactusDisk *cactusDisk) {
    Flower *flower = flower_construct3(name, cactusDisk);
    flower->eventTree = eventTree_construct2(flower);
    return flower;
}

Flower *flower_construct(CactusDisk *cactusDisk) {
    return flower_construct2(cactusDisk_getUniqueID(cactusDisk), cactusDisk);
}

void flower_destruct(Flower *flower, int64_t recursive) {
    Flower_GroupIterator *iterator;
    Sequence *sequence;
    End *end;
    Block *block;
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

    flower_destructFaces(flower);
    stSortedSet_destruct(flower->faces);

    assert(flower_getEventTree(flower) != NULL);
    eventTree_destruct(flower_getEventTree(flower));

    while ((sequence = flower_getFirstSequence(flower)) != NULL) {
        sequence_destruct(sequence);
    }
    stSortedSet_destruct(flower->sequences);

    while ((chain = flower_getFirstChain(flower)) != NULL) {
        chain_destruct(chain);
    }
    stSortedSet_destruct(flower->chains);

    while ((end = flower_getFirstEnd(flower)) != NULL) {
        end_destruct(end);
    }
    stSortedSet_destruct(flower->caps);
    stSortedSet_destruct(flower->ends);

    while ((block = flower_getFirstBlock(flower)) != NULL) {
        block_destruct(block);
    }
    stSortedSet_destruct(flower->segments);
    stSortedSet_destruct(flower->blocks);

    while ((group = flower_getFirstGroup(flower)) != NULL) {
        group_destruct(group);
    }
    stSortedSet_destruct(flower->groups);

    free(flower);
}

Name flower_getName(Flower *flower) {
    return flower->name;
}

CactusDisk *flower_getCactusDisk(Flower *flower) {
    return flower->cactusDisk;
}

EventTree *flower_getEventTree(Flower *flower) {
    return flower->eventTree;
}

Sequence *flower_getFirstSequence(Flower *flower) {
    return stSortedSet_getFirst(flower->sequences);
}

Sequence *flower_getSequence(Flower *flower, Name name) {
    Sequence sequence;
    MetaSequence metaSequence;
    sequence.metaSequence = &metaSequence;
    metaSequence.name = name;
    return stSortedSet_search(flower->sequences, &sequence);
}

int64_t flower_getSequenceNumber(Flower *flower) {
    return stSortedSet_size(flower->sequences);
}

Flower_SequenceIterator *flower_getSequenceIterator(Flower *flower) {
    return stSortedSet_getIterator(flower->sequences);
}

Sequence *flower_getNextSequence(Flower_SequenceIterator *sequenceIterator) {
    return stSortedSet_getNext(sequenceIterator);
}

Sequence *flower_getPreviousSequence(Flower_SequenceIterator *sequenceIterator) {
    return stSortedSet_getPrevious(sequenceIterator);
}

Flower_SequenceIterator *flower_copySequenceIterator(Flower_SequenceIterator *sequenceIterator) {
    return stSortedSet_copyIterator(sequenceIterator);
}

void flower_destructSequenceIterator(Flower_SequenceIterator *sequenceIterator) {
    stSortedSet_destructIterator(sequenceIterator);
}

Cap *flower_getFirstCap(Flower *flower) {
    return stSortedSet_getFirst(flower->caps);
}

Cap *flower_getCap(Flower *flower, Name name) {
    Cap cap;
    CapContents capContents;
    cap.capContents = &capContents;
    cap.capContents->instance = name;
    return stSortedSet_search(flower->caps, &cap);
}

int64_t flower_getCapNumber(Flower *flower) {
    return stSortedSet_size(flower->caps);
}

Flower_CapIterator *flower_getCapIterator(Flower *flower) {
    return stSortedSet_getIterator(flower->caps);
}

Cap *flower_getNextCap(Flower_CapIterator *capIterator) {
    return stSortedSet_getNext(capIterator);
}

Cap *flower_getPreviousCap(Flower_CapIterator *capIterator) {
    return stSortedSet_getPrevious(capIterator);
}

Flower_CapIterator *flower_copyCapIterator(Flower_CapIterator *capIterator) {
    return stSortedSet_copyIterator(capIterator);
}

void flower_destructCapIterator(Flower_CapIterator *capIterator) {
    stSortedSet_destructIterator(capIterator);
}

End *flower_getFirstEnd(Flower *flower) {
    return stSortedSet_getFirst(flower->ends);
}

End *flower_getEnd(Flower *flower, Name name) {
    End end;
    EndContents endContents;
    end.endContents = &endContents;
    endContents.name = name;
    end.orientation = 1;
    return stSortedSet_search(flower->ends, &end);
}

int64_t flower_getEndNumber(Flower *flower) {
    return stSortedSet_size(flower->ends);
}

int64_t flower_getBlockEndNumber(Flower *flower) {
    return flower_getBlockNumber(flower) * 2;
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
    return stSortedSet_getIterator(flower->ends);
}

End *flower_getNextEnd(Flower_EndIterator *endIterator) {
    return stSortedSet_getNext(endIterator);
}

End *flower_getPreviousEnd(Flower_EndIterator *endIterator) {
    return stSortedSet_getPrevious(endIterator);
}

Flower_EndIterator *flower_copyEndIterator(Flower_EndIterator *endIterator) {
    return stSortedSet_copyIterator(endIterator);
}

void flower_destructEndIterator(Flower_EndIterator *endIterator) {
    stSortedSet_destructIterator(endIterator);
}

Segment *flower_getFirstSegment(Flower *flower) {
    return stSortedSet_getFirst(flower->segments);
}

Segment *flower_getSegment(Flower *flower, Name name) {
    Segment segment;
    segment.name = name;
    return stSortedSet_search(flower->segments, &segment);
}

int64_t flower_getSegmentNumber(Flower *flower) {
    return stSortedSet_size(flower->segments);
}

Flower_SegmentIterator *flower_getSegmentIterator(Flower *flower) {
    return stSortedSet_getIterator(flower->segments);
}

Segment *flower_getNextSegment(Flower_SegmentIterator *segmentIterator) {
    return stSortedSet_getNext(segmentIterator);
}

Segment *flower_getPreviousSegment(Flower_SegmentIterator *segmentIterator) {
    return stSortedSet_getPrevious(segmentIterator);
}

Flower_SegmentIterator *flower_copySegmentIterator(Flower_SegmentIterator *segmentIterator) {
    return stSortedSet_copyIterator(segmentIterator);
}

void flower_destructSegmentIterator(Flower_SegmentIterator *segmentIterator) {
    stSortedSet_destructIterator(segmentIterator);
}

Block *flower_getFirstBlock(Flower *flower) {
    return stSortedSet_getFirst(flower->blocks);
}

Block *flower_getBlock(Flower *flower, Name name) {
    Block block;
    BlockContents blockContents;
    block.blockContents = &blockContents;
    blockContents.name = name;
    return stSortedSet_search(flower->blocks, &block);
}

int64_t flower_getBlockNumber(Flower *flower) {
    return stSortedSet_size(flower->blocks);
}

Flower_BlockIterator *flower_getBlockIterator(Flower *flower) {
    return stSortedSet_getIterator(flower->blocks);
}

Block *flower_getNextBlock(Flower_BlockIterator *blockIterator) {
    return stSortedSet_getNext(blockIterator);
}

Block *flower_getPreviousBlock(Flower_BlockIterator *blockIterator) {
    return stSortedSet_getPrevious(blockIterator);
}

Flower_BlockIterator *flower_copyBlockIterator(Flower_BlockIterator *blockIterator) {
    return stSortedSet_copyIterator(blockIterator);
}

void flower_destructBlockIterator(Flower_BlockIterator *blockIterator) {
    stSortedSet_destructIterator(blockIterator);
}

Group *flower_getFirstGroup(Flower *flower) {
    return stSortedSet_getFirst(flower->groups);
}

Group *flower_getGroup(Flower *flower, Name flowerName) {
    Group group;
    group.name = flowerName;
    return stSortedSet_search(flower->groups, &group);
}

int64_t flower_getGroupNumber(Flower *flower) {
    return stSortedSet_size(flower->groups);
}

Flower_GroupIterator *flower_getGroupIterator(Flower *flower) {
    return stSortedSet_getIterator(flower->groups);
}

Group *flower_getNextGroup(Flower_GroupIterator *groupIterator) {
    return stSortedSet_getNext(groupIterator);
}

Group *flower_getPreviousGroup(Flower_GroupIterator *groupIterator) {
    return stSortedSet_getPrevious(groupIterator);
}

Flower_GroupIterator *flower_copyGroupIterator(Flower_GroupIterator *groupIterator) {
    return stSortedSet_copyIterator(groupIterator);
}

void flower_destructGroupIterator(Flower_GroupIterator *groupIterator) {
    stSortedSet_destructIterator(groupIterator);
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
    return stSortedSet_getFirst(flower->chains);
}

Chain *flower_getChain(Flower *flower, Name name) {
    Chain chain;
    chain.name = name;
    return stSortedSet_search(flower->chains, &chain);
}

int64_t flower_getChainNumber(Flower *flower) {
    return stSortedSet_size(flower->chains);
}

int64_t flower_getTrivialChainNumber(Flower *flower) {
    int64_t i = 0;
    Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIt)) != NULL) {
        if (block_isTrivialChain(block)) {
            i++;
        }
    }
    flower_destructBlockIterator(blockIt);
    return i;
}

Flower_ChainIterator *flower_getChainIterator(Flower *flower) {
    return stSortedSet_getIterator(flower->chains);
}

Chain *flower_getNextChain(Flower_ChainIterator *chainIterator) {
    return stSortedSet_getNext(chainIterator);
}

Chain *flower_getPreviousChain(Flower_ChainIterator *chainIterator) {
    return stSortedSet_getPrevious(chainIterator);
}

Flower_ChainIterator *flower_copyChainIterator(Flower_ChainIterator *chainIterator) {
    return stSortedSet_copyIterator(chainIterator);
}

void flower_destructChainIterator(Flower_ChainIterator *chainIterator) {
    stSortedSet_destructIterator(chainIterator);
}

Face *flower_getFirstFace(Flower *flower) {
    return stSortedSet_getFirst(flower->faces);
}

int64_t flower_getFaceNumber(Flower *flower) {
    return stSortedSet_size(flower->faces);
}

Flower_FaceIterator *flower_getFaceIterator(Flower *flower) {
    return stSortedSet_getIterator(flower->faces);
}

Face *flower_getNextFace(Flower_FaceIterator *faceIterator) {
    return stSortedSet_getNext(faceIterator);
}

Face *flower_getPreviousFace(Flower_FaceIterator *faceIterator) {
    return stSortedSet_getPrevious(faceIterator);
}

Flower_FaceIterator *flower_copyFaceIterator(Flower_FaceIterator *faceIterator) {
    return stSortedSet_copyIterator(faceIterator);
}

void flower_destructFaceIterator(Flower_FaceIterator *faceIterator) {
    stSortedSet_destructIterator(faceIterator);
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
        assert(group_getEndNumber(group) > 0);
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
    flower_destructCapIterator(chainIterator);

    //We check built trees in here.
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        end_check(end);
        end_check(end_getReverse(end)); //We will test everything backwards also.
    }
    flower_destructEndIterator(endIterator);

    if (flower_builtFaces(flower)) {
        Flower_FaceIterator *faceIterator = flower_getFaceIterator(flower);
        Face *face;
        while ((face = flower_getNextFace(faceIterator)) != NULL) {
            face_check(face);
        }
        flower_destructFaceIterator(faceIterator);
        face_checkFaces(flower);
    } else {
        cactusCheck(flower_getFaceNumber(flower) == 0);
    }

    if (flower_builtBlocks(flower)) { //Note that a flower for which the blocks are not yet built must be a leaf.
        Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
        Block *block;
        while ((block = flower_getNextBlock(blockIterator)) != NULL) {
            block_check(block);
            block_check(block_getReverse(block)); //We will test everything backwards also.
        }
        flower_destructBlockIterator(blockIterator);
    } else {
        cactusCheck(flower_isLeaf(flower)); //Defensive
        cactusCheck(flower_isTerminal(flower)); //Checks that a flower without built blocks is a leaf and does not
        //contain any blocks.
    }

    Flower_SequenceIterator *sequenceIterator = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(sequenceIterator)) != NULL) {
        sequence_check(sequence);
    }
    flower_destructSequenceIterator(sequenceIterator);
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
    return flower->builtTrees;
}

void flower_setBuiltTrees(Flower *flower, bool b) {
    flower->builtTrees = b;
}

bool flower_builtFaces(Flower *flower) {
    return flower->builtFaces;
}

void flower_setBuildFaces(Flower *flower, bool b) {
    flower->builtFaces = b;
    if (flower_builtFaces(flower)) {
        flower_reconstructFaces(flower);
    }
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
    if (!flower_isLeaf(flower) && flower_getParentGroup(flower) != NULL && flower_getBlockNumber(flower) == 0) { //We will remove this flower..
        Group *parentGroup = flower_getParentGroup(flower); //This group will be destructed
        //Deal with any parent chain..
        if (group_isLink(parentGroup)) {
            link_split(group_getLink(parentGroup));
        }
        Flower *parentFlower = group_getFlower(parentGroup); //We will add the groups in the flower to the parent

        /*
         * For each group in the flower we take its nested flower and attach it to the parent.
         */
        Group *group;
        Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(groupIt)) != NULL) {
            if (!group_isLeaf(group)) {
                //Copy the group into the parent..
                Flower *nestedFlower = group_getNestedFlower(group);
                assert(nestedFlower != NULL);
                Group *newParentGroup = group_construct(parentFlower, nestedFlower);
                flower_setParentGroup(nestedFlower, newParentGroup);
                group_constructChainForLink(newParentGroup);
            } else {
                Group *newParentGroup = group_construct2(parentFlower);
                End *end;
                Group_EndIterator *endIt = group_getEndIterator(group);
                while ((end = group_getNextEnd(endIt)) != NULL) {
                    End *parentEnd = flower_getEnd(parentFlower, end_getName(end));
                    assert(parentEnd != NULL);
                    end_setGroup(parentEnd, newParentGroup);
                }
                group_destructEndIterator(endIt);
                group_constructChainForLink(newParentGroup);
            }
        }
        flower_destructGroupIterator(groupIt);

        //The group attached to the flower should now be empty
        assert(group_getEndNumber(parentGroup) == 0);
        group_destruct(parentGroup);

        //Now wipe the flower out..
        cactusDisk_deleteFlowerFromDisk(flower_getCactusDisk(flower), flower);
        flower_destruct(flower, 0);
        return 1;
    }
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
        parentGroup->leafGroup = 1;
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
        assert(flower_getBlockNumber(flower) == 0);
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
                flower_setBuiltTrees(nestedFlower, flower_builtTrees(flower));
                flower_setBuildFaces(nestedFlower, flower_builtFaces(flower));
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

void flower_setEventTree(Flower *flower, EventTree *eventTree) {
    if (flower_getEventTree(flower) != NULL) {
        eventTree_destruct(flower_getEventTree(flower));
    }
    flower->eventTree = eventTree;
}

void flower_removeEventTree(Flower *flower, EventTree *eventTree) {
    assert(flower_getEventTree(flower) == eventTree);
    flower->eventTree = NULL;
}

void flower_addSequence(Flower *flower, Sequence *sequence) {
    assert(stSortedSet_search(flower->sequences, sequence) == NULL);
    stSortedSet_insert(flower->sequences, sequence);
}

void flower_removeSequence(Flower *flower, Sequence *sequence) {
    assert(stSortedSet_search(flower->sequences, sequence) != NULL);
    stSortedSet_remove(flower->sequences, sequence);
}

void flower_addCap(Flower *flower, Cap *cap) {
    cap = cap_getPositiveOrientation(cap);
    assert(stSortedSet_search(flower->caps, cap) == NULL);
    stSortedSet_insert(flower->caps, cap);
}

void flower_removeCap(Flower *flower, Cap *cap) {
    cap = cap_getPositiveOrientation(cap);
    assert(stSortedSet_search(flower->caps, cap) != NULL);
    stSortedSet_remove(flower->caps, cap);
}

void flower_addEnd(Flower *flower, End *end) {
    end = end_getPositiveOrientation(end);
    assert(stSortedSet_search(flower->ends, end) == NULL);
    stSortedSet_insert(flower->ends, end);
}

void flower_removeEnd(Flower *flower, End *end) {
    end = end_getPositiveOrientation(end);
    assert(stSortedSet_search(flower->ends, end) != NULL);
    stSortedSet_remove(flower->ends, end);
}

void flower_addSegment(Flower *flower, Segment *segment) {
    segment = segment_getPositiveOrientation(segment);
    assert(stSortedSet_search(flower->segments, segment) == NULL);
    stSortedSet_insert(flower->segments, segment);
}

void flower_removeSegment(Flower *flower, Segment *segment) {
    segment = segment_getPositiveOrientation(segment);
    assert(stSortedSet_search(flower->segments, segment) != NULL);
    stSortedSet_remove(flower->segments, segment);
}

void flower_addBlock(Flower *flower, Block *block) {
    block = block_getPositiveOrientation(block);
    assert(stSortedSet_search(flower->blocks, block) == NULL);
    stSortedSet_insert(flower->blocks, block);
}

void flower_removeBlock(Flower *flower, Block *block) {
    block = block_getPositiveOrientation(block);
    assert(stSortedSet_search(flower->blocks, block) != NULL);
    stSortedSet_remove(flower->blocks, block);
}

void flower_addChain(Flower *flower, Chain *chain) {
    assert(stSortedSet_search(flower->chains, chain) == NULL);
    stSortedSet_insert(flower->chains, chain);
}

void flower_removeChain(Flower *flower, Chain *chain) {
    assert(stSortedSet_search(flower->chains, chain) != NULL);
    stSortedSet_remove(flower->chains, chain);
}

void flower_addGroup(Flower *flower, Group *group) {
    assert(stSortedSet_search(flower->groups, group) == NULL);
    stSortedSet_insert(flower->groups, group);
}

void flower_removeGroup(Flower *flower, Group *group) {
    assert(stSortedSet_search(flower->groups, group) != NULL);
    stSortedSet_remove(flower->groups, group);
}

void flower_setParentGroup(Flower *flower, Group *group) {
    //assert(flower->parentFlowerName == NULL_NAME); we can change this if merging the parent flowers, so this no longer applies.
    flower->parentFlowerName = flower_getName(group_getFlower(group));
}

void flower_addFace(Flower *flower, Face *face) {
    assert(stSortedSet_search(flower->faces, face) == NULL);
    stSortedSet_insert(flower->faces, face);
}

void flower_removeFace(Flower *flower, Face *face) {
    assert(stSortedSet_search(flower->faces, face) != NULL);
    stSortedSet_remove(flower->faces, face);
}

/*
 * Serialisation functions.
 */

void flower_writeBinaryRepresentation(Flower *flower, void(*writeFn)(const void * ptr, size_t size, size_t count)) {
    Flower_SequenceIterator *sequenceIterator;
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
    binaryRepresentation_writeBool(flower_builtTrees(flower), writeFn);
    binaryRepresentation_writeBool(flower_builtFaces(flower), writeFn);
    binaryRepresentation_writeName(flower->parentFlowerName, writeFn);

    if (flower_getEventTree(flower) != NULL) {
        eventTree_writeBinaryRepresentation(flower_getEventTree(flower), writeFn);
    }

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

    binaryRepresentation_writeElementType(CODE_FLOWER, writeFn); //this avoids interpretting things wrong.
}

Flower *flower_loadFromBinaryRepresentation(void **binaryString, CactusDisk *cactusDisk) {
    Flower *flower = NULL;
    bool buildFaces;
    if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_FLOWER) {
        binaryRepresentation_popNextElementType(binaryString);
        flower = flower_construct3(binaryRepresentation_getName(binaryString), cactusDisk); //Constructed without an event tree.
        flower_setBuiltBlocks(flower, binaryRepresentation_getBool(binaryString));
        flower_setBuiltTrees(flower, binaryRepresentation_getBool(binaryString));
        buildFaces = binaryRepresentation_getBool(binaryString);
        flower->parentFlowerName = binaryRepresentation_getName(binaryString);
        eventTree_loadFromBinaryRepresentation(binaryString, flower);
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
        flower_setBuildFaces(flower, buildFaces);
        assert(binaryRepresentation_popNextElementType(binaryString) == CODE_FLOWER);
    }
    return flower;
}
