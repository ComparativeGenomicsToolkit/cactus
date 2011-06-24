#include "cactus.h"
#include "sonLib.h"

Cap *getCapForEvent(End *end, Event *event) {

}

/*
 * Functions to finally set the coordinates.
 */

void setCoordinates1(stList *caps, bool side, void **extraArgs) {
    MetaSequence *metaSequence = extraArgs[0];
    int32_t *coordinate = extraArgs[1];
    for (int32_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        Flower *flower = end_getFlower(cap_getEnd(cap));
        Sequence *sequence;
        if ((sequence = flower_getSequence(flower,
                metaSequence_getName(metaSequence))) == NULL) {
            sequence = sequence_construct(metaSequence, flower, side);
        }
        cap_setCoordinates(sequence, *coordinate);
    }
}

void setCoordinates2(stList *caps, bool side, void **extraArgs) {
    MetaSequence *metaSequence = extraArgs[0];
    int32_t *coordinate = extraArgs[1];
    assert(side);
    setCoordinates1(caps, sequence, coordinate, 1);
    Cap *cap = stList_get(caps, 0);
    Segment *segment = cap_getSegment(cap);
    assert(segment != NULL);
    *coordinate += segment_getLength(segment);
}

/*
 * Function to get length of the new sequence.
 */
void getSequenceLength(stList *caps, bool side, void **extraArgs) {
    int32_t *coordinate = extraArgs[0];
    if (side) {
        Cap *cap = stList_get(caps, 0);
        Segment *segment = cap_getSegment(cap);
        assert(segment != NULL);
        *coordinate += segment_getLength(segment);
    }
}

/*
 * Function to construct the sequence
 */

char *getConsensusString(Block *block) {
    int32_t *iA = st_calloc(block_getLength(block) * 5, sizeof(int32_t));

    Segment *segment;
    Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
    upperCounts = st_calloc(block_getLength(block), sizeof(int32_t));
    while ((segment = block_getNext(instanceIt)) != NULL) {
        if (segment_getSequence(segment) != NULL) {
            char *string = segment_getString(segment);
            assert(string != NULL);
            for (int32_t i = 0; i < segment_getLength(segment); i++) {
                upperCount[i] += toupper(string[i]) == string[i] ? 1 : 0;
                switch (toupper(string[i])) {
                    case 'A':
                        iA[i * 5]++;
                        break;
                    case 'C':
                        iA[i * 5 + 1]++;
                        break;
                    case 'G':
                        iA[i * 5 + 2]++;
                        break;
                    case 'T':
                        iA[i * 5 + 3]++;
                        break;
                    default:
                        iA[i * 5 + 4]++;
                }
            }
        }
    }
    block_destructInstanceIterator(instanceIt);

    char *sequence = st_malloc(sizeof(char) * (block_getLength(block) + 1));
    for (int32_t i = 0; i < block_getLength(block); i++) {
        int32_t base = 4;
        int32_t baseCount = iA[i*5 + 4];
        for(int32_t j=0; j<4; j++) {
            int32_t baseCount2 = iA[i*5 + j];
            if(iA[i*5 + j] > baseCount || (iA[i*5 + j] >= baseCount && st_random() > 0.5)) {
                base = j;
                baseCount = iA[i*5 + j];
            }
        }
        switch(base) {
            case 0:
                base = 'a';
                break;
            case 1:
                base = 'c';
                break;
            case 2:
                base = 'g';
                break;
            case 3:
                base = 't';
                break;
            default:
                base = 'N';
        }
        sequence[i] = upperCount[i] >= block_getInstanceNumber(block)/2 ? toupper(base) : base;
    }
    free(iA);

    return sequence;
}

void getSequence(stList *caps, bool side, void **extraArgs) {
    int32_t *coordinate = extraArgs[0];
    char *sequence = extraArgs[1];
    if (side) {
        Cap *cap = stList_get(caps, 0);
        Segment *segment = cap_getSegment(cap);
        assert(segment != NULL);
        char *string = getConsensusString(segment_getBlock(segment));
        for (int32_t i = 0; i < block_getLength(block); i++) {
            sequence[(*coordinate)++] = string[i];
        }
        free(string);
    }
}

Cap *getCapUp(Cap *cap) {
    while (1) {
        if (end_isBlockEnd(cap_getEnd(cap))) {
            return cap;
        }
        Group *parentGroup = flower_getParentGroup(
                end_getFlower(cap_getEnd(cap)));
        if (parentGroup == NULL) {
            return cap;
        }
        Cap *parentCap = flower_getCap(group_getFlower(parentGroup),
                cap_getName(cap));
        cap = cap_getOrientation(cap) ? parentCap : cap_getReverse(parentCap);
    }
}

stList *getCapsDown(Cap *cap) {
    stList *caps = stList_construct();
    while (1) {
        stList_append(caps, cap);
        Flower *nestedFlower = group_getNestedFlower(
                end_getGroup(cap_getEnd(cap)));
        if (nestedFlower != NULL) {
            Cap *nestedCap = flower_getCap(nestedFlower, cap_getName(cap));
            assert(nestedCap != NULL);
            cap = cap_getOrientation(cap) ? nestedCap : cap_getReverse(
                    nestedCap);
        } else {
            break;
        }
    }
    return cap;
}

void traverseCapsInOrder(Flower *flower, Cap *cap, void **extraArgs,
        void(*fn1)(stList *caps, bool side, void **extraArgs),
        void(*fn2)(stList *caps, bool side, void **extraArgs)) {
    assert(cap_getStrand(cap));
    while (1) {
        stList *caps = getCapsDown(cap);
        fn1(caps, metaSequence, coordinate, 0);
        cap = getCapUp(cap_getAdjacency(stList_pop(caps)));
        stList_destruct(caps);
        assert(cap != NULL);
        if (cap_getSegment(cap) != NULL) {
            caps = getCapsDown(cap);
            fn2(caps, metaSequence, coordinate);
            stList_destruct(caps);
            cap = cap_getOtherSegmentCap(cap);
        } else {
            assert(end_isStubEnd(cap_getEnd(cap)));
            break;
        }
    }
}

void addSequence(Flower *flower, Cap *cap) {
    int32_t coordinate = 0;
    traverseCapsInOrder(flower, cap, &coordinate, getSequenceLength, getSequenceLength);
    //Allocate the sequence..

}


