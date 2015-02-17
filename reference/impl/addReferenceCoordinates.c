/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include "cactus.h"
#include "sonLib.h"
#include "recursiveThreadBuilder.h"
#include "blockMLString.h"

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
    //assert(0);
    return NULL;
}

static int64_t traceThreadLength(Cap *cap, Cap **terminatingCap) {
    /*
     * Gets the length in bases of the thread in the flower, starting from a given attached stub cap.
     * The thread length includes the lengths of adjacencies that it contains.
     *
     * Terminating cap is initialised with the final cap on the thread from cap.
     */
    assert(end_isStubEnd(cap_getEnd(cap)));
    int64_t threadLength = 0;
    while (1) {
        assert(cap_getCoordinate(cap) != INT64_MAX);
        int64_t adjacencyLength = cap_getCoordinate(cap);
        threadLength += adjacencyLength;
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        assert(adjacencyLength == cap_getCoordinate(adjacentCap));
        //Traverse any block..
        if (cap_getSegment(adjacentCap) != NULL) {
            threadLength += segment_getLength(cap_getSegment(adjacentCap));
            cap = cap_getOtherSegmentCap(adjacentCap);
            assert(cap != NULL);
        } else {
            assert(end_isStubEnd(cap_getEnd(adjacentCap)));
            *terminatingCap = adjacentCap;
            return threadLength;
        }
    }
    return 1;
}

static Cap *copyCapToParent(Cap *cap, stList *recoveredCaps) {
    /*
     * Get the adjacent stub end by looking at the reference adjacency in the parent.
     */
    End *end = cap_getEnd(cap);
    assert(end != NULL);
    Group *parentGroup = flower_getParentGroup(end_getFlower(end));
    assert(parentGroup != NULL);
    End *copiedEnd = end_copyConstruct(end, group_getFlower(parentGroup));
    end_setGroup(copiedEnd, parentGroup); //Set group
    Cap *copiedCap = end_getInstance(copiedEnd, cap_getName(cap));
    assert(copiedCap != NULL);
    copiedCap = cap_getStrand(copiedCap) ? copiedCap : cap_getReverse(copiedCap);
    if (!cap_getSide(copiedCap)) {
        stList_append(recoveredCaps, copiedCap);
    }
    return copiedCap;
}

static void setAdjacencyLength(Cap *cap, Cap *adjacentCap, int64_t adjacencyLength) {
    //Set the coordinates of the caps to the adjacency size
    cap_setCoordinates(cap, adjacencyLength, cap_getStrand(cap), NULL);
    cap_setCoordinates(adjacentCap, adjacencyLength, cap_getStrand(adjacentCap),
    NULL);
    assert(cap_getCoordinate(cap) == cap_getCoordinate(adjacentCap));
}

static void setAdjacencyLengthsAndRecoverNewCapsAndBrokenAdjacencies(Cap *cap, stList *recoveredCaps) {
    /*
     * Sets the coordinates of the caps to be equal to the length of the adjacency sequence between them.
     * Used to build the reference sequence bottom up.
     *
     * One complexity is that a reference thread between the two caps
     * in each flower f may be broken into two in the children of f.
     * Therefore, for each flower f first identify attached stub ends present in the children of f that are
     * not present in f and copy them into f, reattaching the reference caps as needed.
     */
    while (1) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        assert(cap_getCoordinate(cap) == INT64_MAX);
        assert(cap_getCoordinate(adjacentCap) == INT64_MAX);
        assert(cap_getStrand(cap) == cap_getStrand(adjacentCap));
        assert(cap_getSide(cap) != cap_getSide(adjacentCap));
        Group *group = end_getGroup(cap_getEnd(cap));
        assert(group != NULL);
        if (!group_isLeaf(group)) { //Adjacency is not terminal, so establish its sequence.
            Flower *nestedFlower = group_getNestedFlower(group);
            Cap *nestedCap = flower_getCap(nestedFlower, cap_getName(cap));
            assert(nestedCap != NULL);
            Cap *nestedAdjacentCap = flower_getCap(nestedFlower, cap_getName(adjacentCap));
            assert(nestedAdjacentCap != NULL);
            Cap *breakerCap;
            int64_t adjacencyLength = traceThreadLength(nestedCap, &breakerCap);
            assert(cap_getOrientation(nestedAdjacentCap));
            if (cap_getPositiveOrientation(breakerCap) != nestedAdjacentCap) { //The thread is broken at the lower level.
                //Copy cap into higher level graph.
                breakerCap = copyCapToParent(breakerCap, recoveredCaps);
                assert(cap_getSide(breakerCap));
                cap_makeAdjacent(cap, breakerCap);
                setAdjacencyLength(cap, breakerCap, adjacencyLength);
                adjacencyLength = traceThreadLength(nestedAdjacentCap, &breakerCap);
                assert(cap_getPositiveOrientation(breakerCap) != cap);
                breakerCap = copyCapToParent(breakerCap, recoveredCaps);
                assert(!cap_getSide(breakerCap));
                cap_makeAdjacent(breakerCap, adjacentCap);
                setAdjacencyLength(adjacentCap, breakerCap, adjacencyLength);
            } else { //The thread is not broken at the lower level
                setAdjacencyLength(cap, adjacentCap, adjacencyLength);
            }
        } else {
            //Set the coordinates of the caps to the adjacency size
            setAdjacencyLength(cap, adjacentCap, 0);
        }
        if ((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
            break;
        }
    }
}

static void recoverBrokenAdjacencies(Flower *flower, stList *recoveredCaps, Name referenceEventName) {
    /*
     * Find reference intervals that are book-ended by stubs created in a child flower.
     */
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while((group = flower_getNextGroup(groupIt)) != NULL) {
        Flower *nestedFlower;
        if((nestedFlower = group_getNestedFlower(group)) != NULL) {
            Flower_EndIterator *endIt = flower_getEndIterator(nestedFlower);
            End *childEnd;
            while((childEnd = flower_getNextEnd(endIt)) != NULL) {
                if(end_isStubEnd(childEnd) && flower_getEnd(flower, end_getName(childEnd)) == NULL) { //We have a thread we need to promote
                    Cap *childCap = getCapForReferenceEvent(childEnd, referenceEventName); //The cap in the reference
                    assert(childCap != NULL);
                    assert(!end_isAttached(childEnd));
                    childCap = cap_getStrand(childCap) ? childCap : cap_getReverse(childCap);
                    if (!cap_getSide(childCap)) {
                        Cap *adjacentChildCap = NULL;
                        int64_t adjacencyLength = traceThreadLength(childCap, &adjacentChildCap);
                        Cap *cap = copyCapToParent(childCap, recoveredCaps);
                        assert(adjacentChildCap != NULL);
                        assert(!end_isAttached(cap_getEnd(adjacentChildCap)));
                        assert(!cap_getSide(cap));
                        Cap *adjacentCap = copyCapToParent(adjacentChildCap, recoveredCaps);
                        cap_makeAdjacent(cap, adjacentCap);
                        setAdjacencyLength(cap, adjacentCap, adjacencyLength);
                    }
                }
            }
            flower_destructEndIterator(endIt);
        }
    }
    flower_destructGroupIterator(groupIt);
}

static char *terminalAdjacencyWriteFn(Cap *cap) {
    return stString_copy("");
}

//static Name segmentWriteFnOutgroupEventName;
static stHash *segmentWriteFn_flowerToPhylogeneticTreeHash;


static char *segmentWriteFn(Segment *segment) {
    stTree *phylogeneticTree = stHash_search(segmentWriteFn_flowerToPhylogeneticTreeHash, block_getFlower(segment_getBlock(segment)));
    assert(phylogeneticTree != NULL);
    char *segmentString = getMaximumLikelihoodString(phylogeneticTree, segment_getBlock(segment));
    //char *segmentString = getConsensusString(segment_getBlock(segment), segmentWriteFnOutgroupEventName);
    //We append a zero to a segment string if it is part of block containing only a reference segment, else we append a 1.
    //We use these boolean values to determine if a sequence contains only these trivial strings, and is therefore trivial.
    char *appendedSegmentString = stString_print("%s%c ", segmentString, block_getInstanceNumber(segment_getBlock(segment)) == 1 ? '0' : '1');
    free(segmentString);
    return appendedSegmentString;
}

/*
 * A thread is trivial if all the segments it contains come from blocks containing only a reference segment.
 * These reference only segments represent scaffold gaps. At the same time, it processes the thread string
 * to remove the boolean values use to indicate if a thread is trivial or not.
 */
static bool isTrivialString(char **threadString) {
    stList *strings = stString_split(*threadString); //Split splits into individual segments.
    bool trivialString = 1;
    for(int64_t i=0; i<stList_length(strings); i++) {
        char *segmentString = stList_get(strings, i);
        int64_t j = strlen(segmentString)-1; //Location of the boolean value within a segment.
        assert(j > 0);
        assert(segmentString[j] == '0' || segmentString[j] == '1');
        if(segmentString[j] == '1') { //Found a non-trivial segment, hence the thread is non-trivial.
            trivialString = 0;
        }
        segmentString[j] = '\0';
    }
    free(*threadString); //Free old thread string. Doing it this way is a bit more memory efficient, as we don't keep two copies of the string around.
    *threadString = stString_join2("", strings); //Concatenation makes one sequence, now without the booleans.
    stList_destruct(strings);
    return trivialString;
}

static MetaSequence *addMetaSequence(Flower *flower, Cap *cap, int64_t index, char *string, bool trivialString) {
    /*
     * Adds a meta sequence representing a top level thread to the cactus disk.
     * The sequence is all 'N's at this stage.
     */
    Event *referenceEvent = cap_getEvent(cap);
    assert(referenceEvent != NULL);
    char *sequenceName = stString_print("%srefChr%" PRIi64 "", event_getHeader(referenceEvent), index);
    //char *sequenceName = stString_print("refChr%" PRIi64 "", index);
    MetaSequence *metaSequence = metaSequence_construct3(1, strlen(string), string, sequenceName,
            event_getName(referenceEvent), trivialString, flower_getCactusDisk(flower));
    free(sequenceName);
    return metaSequence;
}

static int64_t setCoordinates(Flower *flower, MetaSequence *metaSequence, Cap *cap, int64_t coordinate) {
    /*
     * Sets the coordinates of the reference thread and sets the bases of the actual sequence
     * that of the consensus.
     */
    Sequence *sequence = flower_getSequence(flower, metaSequence_getName(metaSequence));
    if (sequence == NULL) {
        sequence = sequence_construct(metaSequence, flower);
    }
    while (1) {
        assert(cap_getStrand(cap));
        assert(!cap_getSide(cap));
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        assert(cap_getStrand(adjacentCap));
        assert(cap_getSide(adjacentCap));
        int64_t adjacencyLength = cap_getCoordinate(cap);
        assert(adjacencyLength == cap_getCoordinate(adjacentCap));
        assert(adjacencyLength != INT64_MAX);
        assert(adjacencyLength >= 0);
        cap_setCoordinates(cap, coordinate, 1, sequence);
        coordinate += adjacencyLength + 1;
        cap_setCoordinates(adjacentCap, coordinate, 1, sequence);
        //Traverse any block..
        if ((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
            break;
        }
        coordinate += segment_getLength(cap_getSegment(adjacentCap)) - 1;
    }
    return coordinate;
}

static stList *getCaps(stList *flowers, Name referenceEventName) {
    stList *caps = stList_construct();
    for (int64_t i = 0; i < stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        //Get list of caps
        Flower_EndIterator *endIt = flower_getEndIterator(flower);
        End *end;
        while ((end = flower_getNextEnd(endIt)) != NULL) {
            if (end_isStubEnd(end)) {
                Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
                if(cap != NULL) {
                    cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                    if (!cap_getSide(cap)) {
                        stList_append(caps, cap);
                    }
                }
            }
        }
        flower_destructEndIterator(endIt);
    }
    return caps;
}

void bottomUp(stList *flowers, stKVDatabase *sequenceDatabase, Name referenceEventName, Name outgroupEventName,
        bool isTop, stMatrix *(*generateSubstitutionMatrix)(double)) {
    /*
     * A reference thread between the two caps
     * in each flower f may be broken into two in the children of f.
     * Therefore, for each flower f first identify attached stub ends present in the children of f that are
     * not present in f and copy them into f, reattaching the reference caps as needed.
     */
    stList *caps = getCaps(flowers, referenceEventName);
    for (int64_t i = stList_length(caps) - 1; i >= 0; i--) { //Start from end, as we add to this list.
        setAdjacencyLengthsAndRecoverNewCapsAndBrokenAdjacencies(stList_get(caps, i), caps);
    }
    for(int64_t i=0; i<stList_length(flowers); i++) {
        recoverBrokenAdjacencies(stList_get(flowers, i), caps, referenceEventName);
    }

    //Build the phylogenetic event trees for base calling.
    segmentWriteFn_flowerToPhylogeneticTreeHash = stHash_construct2(NULL, (void (*)(void *))cleanupPhylogeneticTree);
    for(int64_t i=0; i<stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        Event *refEvent = eventTree_getEvent(flower_getEventTree(flower), referenceEventName);
        assert(refEvent != NULL);
        stHash_insert(segmentWriteFn_flowerToPhylogeneticTreeHash, flower, getPhylogeneticTreeRootedAtGivenEvent(refEvent, generateSubstitutionMatrix));
    }

    if (isTop) {
        stList *threadStrings = buildRecursiveThreadsInList(sequenceDatabase, caps, segmentWriteFn,
                terminalAdjacencyWriteFn);
        assert(stList_length(threadStrings) == stList_length(caps));

        int64_t nonTrivialSeqIndex = 0, trivialSeqIndex = stList_length(threadStrings); //These are used as indices for the names of trivial and non-trivial sequences.
        for (int64_t i = 0; i < stList_length(threadStrings); i++) {
            Cap *cap = stList_get(caps, i);
            assert(cap_getStrand(cap));
            assert(!cap_getSide(cap));
            Flower *flower = end_getFlower(cap_getEnd(cap));
            char *threadString = stList_get(threadStrings, i);
            bool trivialString = isTrivialString(&threadString); //This alters the original string
            MetaSequence *metaSequence = addMetaSequence(flower, cap, trivialString ? trivialSeqIndex++ : nonTrivialSeqIndex++,
                    threadString, trivialString);
            free(threadString);
            int64_t endCoordinate = setCoordinates(flower, metaSequence, cap, metaSequence_getStart(metaSequence) - 1);
            (void) endCoordinate;
            assert(endCoordinate == metaSequence_getLength(metaSequence) + metaSequence_getStart(metaSequence));
        }
        stList_setDestructor(threadStrings, NULL); //The strings are already cleaned up by the above loop
        stList_destruct(threadStrings);
    } else {
        buildRecursiveThreads(sequenceDatabase, caps, segmentWriteFn, terminalAdjacencyWriteFn);
    }
    stHash_destruct(segmentWriteFn_flowerToPhylogeneticTreeHash);
    stList_destruct(caps);
}

void topDown(stList *flowers, Name referenceEventName) {
    /*
     * Run on each flower, top down. Sets the coordinates of each reference cap to the correct
     * sequence, and sets the bases of the reference sequence to be consensus bases.
     */
    for (int64_t i = 0; i < stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        Flower_EndIterator *endIt = flower_getEndIterator(flower);
        End *end;
        while ((end = flower_getNextEnd(endIt)) != NULL) {
            Cap *cap = getCapForReferenceEvent(end, referenceEventName); //The cap in the reference
            if (cap != NULL) {
                cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                if (!cap_getSide(cap)) {
                    assert(cap_getCoordinate(cap) != INT64_MAX);
                    Sequence *sequence = cap_getSequence(cap);
                    assert(sequence != NULL);
                    Group *group = end_getGroup(end);
                    if (!group_isLeaf(group)) {
                        Flower *nestedFlower = group_getNestedFlower(group);
                        Cap *nestedCap = flower_getCap(nestedFlower, cap_getName(cap));
                        assert(nestedCap != NULL);
                        nestedCap = cap_getStrand(nestedCap) ? nestedCap : cap_getReverse(nestedCap);
                        assert(cap_getStrand(nestedCap));
                        assert(!cap_getSide(nestedCap));
                        int64_t endCoordinate = setCoordinates(nestedFlower, sequence_getMetaSequence(sequence),
                                nestedCap, cap_getCoordinate(cap));
                        (void) endCoordinate;
                        assert(endCoordinate == cap_getCoordinate(cap_getAdjacency(cap)));
                        assert(endCoordinate
                                        == cap_getCoordinate(
                                                flower_getCap(nestedFlower, cap_getName(cap_getAdjacency(cap)))));
                    }
                }
            }
        }
        flower_destructEndIterator(endIt);
    }
}
