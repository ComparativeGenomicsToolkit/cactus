#include <stdio.h>
#include <ctype.h>
#include "cactus.h"
#include "sonLib.h"

////////////////////////////////////
////////////////////////////////////
//Calculate the length of a reference sequence.
////////////////////////////////////
////////////////////////////////////

static void getSequenceLength(stList *caps, int32_t *coordinate) {
    Cap *cap = stList_get(caps, 0);
    assert(cap_getSide(cap));
    Segment *segment = cap_getSegment(cap);
    if (segment != NULL) {
        *coordinate += segment_getLength(segment);
    }
}

////////////////////////////////////
////////////////////////////////////
//Calculate the bases of a reference sequence.
////////////////////////////////////
////////////////////////////////////

char *getConsensusStringP(stList *strings, int32_t blockLength) {
    //Matrix to store the number of occurrences of each base type, for each column in the block
    int32_t *baseCounts = st_calloc(blockLength * 5, sizeof(int32_t));
    //Array storing number of bases that are upper case letters (non repetitive)..
    int32_t *upperCounts = st_calloc(blockLength, sizeof(int32_t));

    for (int32_t j = 0; j < stList_length(strings); j++) {
        char *string = stList_get(strings, j);
        for (int32_t i = 0; i < blockLength; i++) {
            upperCounts[i] += toupper(string[i]) == string[i] ? 1 : 0;
            switch (toupper(string[i])) {
                case 'A':
                    baseCounts[i * 5]++;
                    break;
                case 'C':
                    baseCounts[i * 5 + 1]++;
                    break;
                case 'G':
                    baseCounts[i * 5 + 2]++;
                    break;
                case 'T':
                    baseCounts[i * 5 + 3]++;
                    break;
                default:
                    baseCounts[i * 5 + 4]++;
            }
        }
    }

    char *string = st_malloc(sizeof(char) * (blockLength + 1));
    string[blockLength] = '\0';
    for (int32_t i = 0; i < blockLength; i++) {
        int32_t base = 4;
        int32_t baseCount = baseCounts[i * 5 + 4];
        for (int32_t j = 0; j < 4; j++) {
            int32_t k = baseCounts[i * 5 + j];
            if (k > baseCount || (base != 4 && k == baseCount && st_random()
                    > 0.5)) {
                base = j;
                baseCount = k;
            }
        }
        switch (base) {
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
            case 4:
                base = 'n';
                break;
            default:
                assert(0);

        }
        string[i]
                = upperCounts[i] >= ((double) stList_length(strings)) / 2 ? toupper(
                        base)
                        : base;
    }
    free(baseCounts);
    free(upperCounts);

    return string;
}

static char *getConsensusString(Block *block) {
    /*
     * Returns a consensus string for a block.
     */

    //Get the strings.
    Segment *segment;
    Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
    stList *strings = stList_construct3(0, free);
    while ((segment = block_getNext(instanceIt)) != NULL) {
        if (segment_getSequence(segment) != NULL) {
            stList_append(strings, segment_getString(segment));
        }
    }
    block_destructInstanceIterator(instanceIt);

    char *string = getConsensusStringP(strings, block_getLength(block));
    stList_destruct(strings);
    return string;
}

static void getString(stList *caps, void **extraArgs) {
    int32_t *coordinate = extraArgs[0];
    char *string = extraArgs[1];
    Cap *cap = stList_get(caps, 0);
    assert(cap_getSide(cap));
    Segment *segment = cap_getSegment(cap);
    if (segment != NULL) {
        char *segmentString = getConsensusString(segment_getBlock(segment));
        for (int32_t i = 0; i < segment_getLength(segment); i++) {
            string[(*coordinate)++] = segmentString[i];
        }
        free(segmentString);
    }
}

////////////////////////////////////
////////////////////////////////////
//Fill in the bases of a reference sequence
////////////////////////////////////
////////////////////////////////////

static void setCoordinates1(stList *caps, void **extraArgs) {
    int32_t *coordinate = extraArgs[0];
    MetaSequence *metaSequence = extraArgs[1];
    for (int32_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        //Add the metasequence to the flower if it isn't already present.
        Flower *flower = end_getFlower(cap_getEnd(cap));
        Sequence *sequence;
        if ((sequence = flower_getSequence(flower,
                metaSequence_getName(metaSequence))) == NULL) {
            sequence = sequence_construct(metaSequence, flower);
        }
        //Now set the coordinates
        assert(cap_getSequence(cap) == NULL);
        assert(cap_getCoordinate(cap) == INT32_MAX);
        cap_setCoordinates(cap, *coordinate, 1, sequence);
        //Check its all set okay..
        assert(cap_getCoordinate(cap) == *coordinate);
        assert(cap_getSequence(cap) == sequence);
        assert(cap_getStrand(cap));
        assert(cap_getCoordinate(cap_getReverse(cap)) == *coordinate);
        assert(cap_getSequence(cap_getReverse(cap)) == sequence);
        assert(!cap_getStrand(cap_getReverse(cap)));
    }
    *coordinate += 1;
}

static void setCoordinates2(stList *caps, void **extraArgs) {
    int32_t *coordinate = extraArgs[0];
    setCoordinates1(caps, extraArgs);
    for (int32_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        assert(cap_getAdjacency(cap) != NULL);
        assert(cap_getSequence(cap_getAdjacency(cap)) != NULL);
        cap_makeAdjacent(cap, cap_getAdjacency(cap)); //This corrects the adjacency if we have modified the strand..
    }
    Cap *cap = stList_get(caps, 0);
    assert(cap_getSide(cap));
    assert(cap_getStrand(cap));
    Segment *segment = cap_getSegment(cap);
    if (segment != NULL) {
        *coordinate += segment_getLength(segment) - 2;
    }
}

////////////////////////////////////
////////////////////////////////////
//Functions to traverse the caps of a sequence
////////////////////////////////////
////////////////////////////////////

static Cap *getCapUp(Cap *cap) {
    /*
     * Gets the highest level version of a cap.
     */
    while (1) {
        assert(cap != NULL);
        cap = cap_getSide(cap) ? cap : cap_getReverse(cap);
        if (end_isBlockEnd(cap_getEnd(cap))) {
            return cap;
        }
        Group *parentGroup = flower_getParentGroup(
                end_getFlower(cap_getEnd(cap)));
        if (parentGroup == NULL) {
            return cap;
        }
        cap = flower_getCap(group_getFlower(parentGroup), cap_getName(cap));
    }
}

static stList *getCapsDown(Cap *cap, bool side) {
    /*
     * Gets a list of a cap in progressively lower level problems, such that
     * the first member of the return list will be the given cap, and the last member
     * will be the lowest level version of the cap.
     */
    stList *caps = stList_construct();
    assert(end_isAttached(cap_getEnd(cap)) || end_isBlockEnd(cap_getEnd(cap)));
    if (end_isStubEnd(cap_getEnd(cap))) {
        assert(flower_getParentGroup(end_getFlower(cap_getEnd(cap))) == NULL);
    }
    while (1) {
        stList_append(caps,
                cap_getSide(cap) == side ? cap : cap_getReverse(cap));
        assert(end_getGroup(cap_getEnd(cap)) != NULL);
        Flower *nestedFlower = group_getNestedFlower(
                end_getGroup(cap_getEnd(cap)));
        if (nestedFlower != NULL) {
            assert(
                    flower_getEnd(nestedFlower, end_getName(cap_getEnd(cap)))
                            != NULL);
            cap = flower_getCap(nestedFlower, cap_getName(cap));
            assert(cap != NULL);
        } else {
            break;
        }
    }
    return caps;
}

void traverseCapsInSequenceOrderFrom3PrimeCap(Cap *cap, void *extraArg,
        void(*_3PrimeFn)(stList *caps, void *extraArg),
        void(*_5PrimeFn)(stList *caps, void *extraArg)) {
    assert(end_isStubEnd(cap_getEnd(cap)));
    assert(end_isAttached(cap_getEnd(cap)));
    assert(flower_getParentGroup(end_getFlower(cap_getEnd(cap))) == NULL);
    while (1) {
        //Call 3' function
        stList *caps = getCapsDown(cap, 0);
        if (_3PrimeFn != NULL) {
            _3PrimeFn(caps, extraArg);
        }
        //Get the adjacent 5 prime cap
        assert(group_isLeaf(end_getGroup(cap_getEnd(stList_peek(caps)))));
        cap = getCapUp(cap_getAdjacency(stList_peek(caps)));
        stList_destruct(caps);
        //Now call 5' function
        caps = getCapsDown(cap, 1);
        if (_5PrimeFn != NULL) {
            _5PrimeFn(caps, extraArg);
        }
        stList_destruct(caps);
        if (cap_getSegment(cap) != NULL) { //Get the opposite 3 prime cap.
            cap = cap_getOtherSegmentCap(cap);
            assert(cap != NULL);
        } else {
            assert(end_isStubEnd(cap_getEnd(cap)));
            assert(end_isAttached(cap_getEnd(cap)));
            assert(
                    flower_getParentGroup(end_getFlower(cap_getEnd(cap)))
                            == NULL);
            break;
        }
    }
}

////////////////////////////////////
////////////////////////////////////
//Function to fill in the coordinates of a reference sequence.
////////////////////////////////////
////////////////////////////////////

static void addReferenceSequence(Cap *cap, int32_t *index) {
    /*
     * Add a sequence for the given cap.
     */

    /*
     * Get the sequence length required.
     */
    int32_t length = 0;
    traverseCapsInSequenceOrderFrom3PrimeCap(cap, &length, NULL,
            (void(*)(stList *, void *)) getSequenceLength);
    st_logDebug("Calculated the length of the sequence: %i\n", length);

    /*
     * Build the string for the sequence.
     */
    char *string = st_malloc(sizeof(char) * (length + 1));
    string[length] = '\0';
    int32_t coordinate = 0;
    void *extraArgs[] = { &coordinate, string };
    traverseCapsInSequenceOrderFrom3PrimeCap(cap, extraArgs, NULL,
            (void(*)(stList *, void *)) getString);
    st_logDebug("Built the string for the sequence\n");

    /*
     * Make the meta sequence.
     */
    Event *event = cap_getEvent(cap);
    char *sequenceName = stString_print("%s_%i", event_getHeader(event), (*index)++);
    MetaSequence *metaSequence = metaSequence_construct(1, length, string,
            sequenceName, event_getName(event),
            flower_getCactusDisk(end_getFlower(cap_getEnd(cap))));
    free(sequenceName);

    /*
     * Now fill in the actual sequence;
     */
    coordinate = 0;
    void *extraArgs2[] = { &coordinate, metaSequence };
    traverseCapsInSequenceOrderFrom3PrimeCap(cap, extraArgs2,
            (void(*)(stList *, void *)) setCoordinates1,
            (void(*)(stList *, void *)) setCoordinates2);
    assert(coordinate == length + 2);
    st_logDebug("Set the coordinates for the new sequence\n");

    /*
     * Cleanup
     */
    free(string);
}

Cap *getCapForReferenceEvent(End *end, Name referenceEventName) {
    /*
     * Get the cap for a given event.
     */
    End_InstanceIterator *it = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(it)) != NULL) {
        if (event_getName(cap_getEvent(cap)) == referenceEventName) {
            end_destructInstanceIterator(it);
            return cap;
        }
    }
    end_destructInstanceIterator(it);
    assert(0);
    return NULL;
}

void addReferenceSequences(Flower *flower, Name referenceEventName) {
#ifdef BEN_DEBUG
    flower_checkRecursive(flower);
#endif

    int32_t index = 0;
    assert(flower_getParentGroup(flower) == NULL);
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end) && end_isAttached(end)) {
            Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
            assert(cap != NULL);
            if (cap_getSequence(cap) == NULL) {
                st_logDebug("Adding the coordinates for cap %s\n",
                        cactusMisc_nameToStringStatic(cap_getName(cap)));
                addReferenceSequence(cap, &index);
            }
        }
    }
    flower_destructEndIterator(endIt);
#ifdef BEN_DEBUG
    flower_checkRecursive(flower);
#endif
}
