/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "endAligner.h"
#include "cactus.h"
#include "sonLib.h"
#include "adjacencySequences.h"
#include "pairwiseAligner.h"

stList *getInducedAlignment(stSortedSet *endAlignment, AdjacencySequence *adjacencySequence) {
    /*
     * Gets an ordered list of pairs from the end alignment for the given adjacency sequence.
     */
    stList *inducedAlignment = stList_construct();
    if (adjacencySequence->strand) {
        AlignedPair *alignedPair = alignedPair_construct(adjacencySequence->subsequenceIdentifier,
                adjacencySequence->start - 1, adjacencySequence->strand, 0, 0, 0, 0, 0);
        AlignedPair *alignedPair2 = stSortedSet_searchGreaterThan(endAlignment, alignedPair);
        if (alignedPair2 != NULL) {
            stSortedSetIterator *it = stSortedSet_getIteratorFrom(endAlignment, alignedPair2);
            while ((alignedPair2 = stSortedSet_getNext(it)) != NULL) {
                if (alignedPair2->subsequenceIdentifier == adjacencySequence->subsequenceIdentifier) {
                    if (alignedPair2->position >= adjacencySequence->start + adjacencySequence->length) {
                        break;
                    }
                    assert(alignedPair2->position >= adjacencySequence->start);
                    if (alignedPair2->strand == adjacencySequence->strand) {
                        stList_append(inducedAlignment, alignedPair2);
                    }
                } else {
                    break;
                }
            }
            stSortedSet_destructIterator(it);
        }
        alignedPair_destruct(alignedPair->reverse);
        alignedPair_destruct(alignedPair);
    } else {
        AlignedPair *alignedPair = alignedPair_construct(adjacencySequence->subsequenceIdentifier,
                adjacencySequence->start + 1, adjacencySequence->strand, 0, 0, 0, 0, 0);
        AlignedPair *alignedPair2 = stSortedSet_searchLessThan(endAlignment, alignedPair);
        if (alignedPair2 != NULL) {
            stSortedSetIterator *it = stSortedSet_getIteratorFrom(endAlignment, alignedPair2);
            stSortedSet_getNext(it);
            stSortedSet_getNext(it); //shift it up..
            while ((alignedPair2 = stSortedSet_getPrevious(it)) != NULL) {
                if (alignedPair2->subsequenceIdentifier == adjacencySequence->subsequenceIdentifier) {
                    if (alignedPair2->position <= adjacencySequence->start - adjacencySequence->length) {
                        break;
                    }
                    assert(alignedPair2->position <= adjacencySequence->start);
                    if (alignedPair2->strand == adjacencySequence->strand) {
                        stList_append(inducedAlignment, alignedPair2);
                    }
                } else {
                    break;
                }
            }
            stSortedSet_destructIterator(it);
        }
        alignedPair_destruct(alignedPair->reverse);
        alignedPair_destruct(alignedPair);
    }
    /*
     * Check the induced alignment
     */
    for (int64_t i = 0; i < stList_length(inducedAlignment); i++) {
        AlignedPair *alignedPair = stList_get(inducedAlignment, i);
        (void) alignedPair;
        assert(alignedPair->subsequenceIdentifier == adjacencySequence->subsequenceIdentifier);
        assert(alignedPair->strand == adjacencySequence->strand);
        if (adjacencySequence->strand) {
            assert(alignedPair->position >= adjacencySequence->start);
            assert(alignedPair->position < adjacencySequence->start + adjacencySequence->length);
        } else {
            assert(alignedPair->position <= adjacencySequence->start);
            assert(alignedPair->position > adjacencySequence->start - adjacencySequence->length);
        }
    }
    return inducedAlignment;
}

/*
 * Runs along and cumulate the score of the pairs, traversing forward through the induced alignment.
 */
static int64_t *cumulateScoreForward(stList *inducedAlignment1) {
    int64_t *iA = st_malloc(sizeof(int64_t) * stList_length(inducedAlignment1));
    int64_t totalScore = 0;
    for (int64_t i = 0; i < stList_length(inducedAlignment1); i++) {
        AlignedPair *alignedPair = stList_get(inducedAlignment1, i);
        totalScore += alignedPair->score;
        iA[i] = totalScore;
    }
    return iA;
}

/*
 * Runs along and cumulate the score of the pairs, traversing backward through the induced alignment.
 */
static int64_t *cumulateScoreBackward(stList *inducedAlignment1) {
    int64_t *iA = st_malloc(sizeof(int64_t) * stList_length(inducedAlignment1));
    int64_t totalScore = 0;
    for (int64_t i = stList_length(inducedAlignment1) - 1; i >= 0; i--) {
        AlignedPair *alignedPair = stList_get(inducedAlignment1, i);
        totalScore += alignedPair->score;
        iA[i] = totalScore;
    }
    return iA;
}

/*
 * Chooses a point along the adjacency sequence at which to filter the two alignments,
 */
static int64_t getCutOff(stList *inducedAlignment1, stList *inducedAlignment2, int64_t *cutOff1, int64_t *cutOff2) {
    int64_t *cScore1 = cumulateScoreForward(inducedAlignment1);
    int64_t *cScore2 = cumulateScoreBackward(inducedAlignment2);

    //Check the score arrays for sanity..
    for (int64_t i = 1; i < stList_length(inducedAlignment1); i++) {
        assert(cScore1[i - 1] < cScore1[i]);
    }
    for (int64_t i = 1; i < stList_length(inducedAlignment2); i++) {
        assert(cScore2[i - 1] > cScore2[i]);
    }

    //Find the score threshold by iterating along each alignment.
    *cutOff1 = 0;
    *cutOff2 = 0;
    int64_t maxScore = -1;
    if (stList_length(inducedAlignment2) > 0) {
        maxScore = cScore2[0];
    }
    int64_t j = 0;
    int64_t pPos1 = INT64_MIN, pPos2 = INT64_MIN;
    for (int64_t i = 0; i < stList_length(inducedAlignment1); i++) {
        AlignedPair *alignedPair1 = stList_get(inducedAlignment1, i);
        assert(alignedPair1->strand);
        assert(pPos1 <= alignedPair1->position);
        pPos1 = alignedPair1->position;
        if (j < stList_length(inducedAlignment2)) {
            do {
                AlignedPair *alignedPair2 = stList_get(inducedAlignment2, j);
                assert(!alignedPair2->strand);
                assert(pPos2 <= alignedPair2->position);
                pPos2 = alignedPair2->position;
                if (alignedPair1->position < alignedPair2->position) {
                    if (cScore1[i] + cScore2[j] >= maxScore) {
                        maxScore = cScore1[i] + cScore2[j];
                        *cutOff1 = i + 1;
                        *cutOff2 = j;
                    }
                    break;
                } else {
                    j++;
                }
            } while (j < stList_length(inducedAlignment2));
        } else {
            if (cScore1[i] >= maxScore) {
                *cutOff1 = stList_length(inducedAlignment1);
                *cutOff2 = j;
                assert(cScore1[stList_length(inducedAlignment1) - 1] >= maxScore);
                maxScore = cScore1[stList_length(inducedAlignment1) - 1];
                break;
            }
        }
    }
    //Cleanup
    free(cScore1);
    free(cScore2);

    return maxScore;
}

void updateDeletedPairs(int64_t subsequenceIdentifier, stHash *deletedAlignedPairCounts) {
	/*
	 * Adds one to count for the given sequenceIdentifier;
	 */
    stIntTuple *i = stIntTuple_construct1(subsequenceIdentifier);
    int64_t *j = stHash_search(deletedAlignedPairCounts, i);
    if(j == NULL) {
        j = st_calloc(1, sizeof(int64_t));
        stHash_insert(deletedAlignedPairCounts, i, j);
    }
    else {
        stIntTuple_destruct(i);
    }
    (*j)++;
}

static void pruneAlignmentsP(stList *inducedAlignment, stSortedSet *endAlignment, int64_t start, int64_t end,
        stSortedSet *pairsToDelete, stHash *deletedAlignedPairCounts) {
    for (int64_t i = start; i < end; i++) {
        AlignedPair *alignedPair = stList_get(inducedAlignment, i);
        if (stSortedSet_search(endAlignment, alignedPair) != NULL) { //can be missing if we are pruning the reverse strand alignment at the same time
            assert(stSortedSet_search(endAlignment, alignedPair->reverse) != NULL);
            updateDeletedPairs(alignedPair->subsequenceIdentifier, deletedAlignedPairCounts);
            updateDeletedPairs(alignedPair->reverse->subsequenceIdentifier, deletedAlignedPairCounts);
            stSortedSet_remove(endAlignment, alignedPair);
            stSortedSet_remove(endAlignment, alignedPair->reverse);
            if (stSortedSet_search(pairsToDelete, alignedPair) == NULL) { // &&
                assert(stSortedSet_search(pairsToDelete, alignedPair->reverse) == NULL);
                stSortedSet_insert(pairsToDelete, alignedPair);
                stSortedSet_insert(pairsToDelete, alignedPair->reverse);
            } else {
                assert(stSortedSet_search(pairsToDelete, alignedPair->reverse) != NULL);
            }
        }
    }
}

static void pruneAlignments(Cap *cap, stList *inducedAlignment1, stList *inducedAlignment2, stSortedSet *endAlignment1,
        stSortedSet *endAlignment2, void *deletedAlignedPairCounts) {
    /*
     * Chooses a point along the adjacency sequence at which to filter the two alignments,
     * then filters the aligned pairs by this point.
     */
    int64_t cutOff1 = 0, cutOff2 = 0;
    getCutOff(inducedAlignment1, inducedAlignment2, &cutOff1, &cutOff2);
    stSortedSet *pairsToDelete = stSortedSet_construct2((void(*)(void *)) alignedPair_destruct);
    //Now do the actual filtering of the alignments.
    pruneAlignmentsP(inducedAlignment1, endAlignment1, cutOff1, stList_length(inducedAlignment1), pairsToDelete, deletedAlignedPairCounts);
    pruneAlignmentsP(inducedAlignment2, endAlignment2, 0, cutOff2, pairsToDelete, deletedAlignedPairCounts);
    stSortedSet_destruct(pairsToDelete);
}

void getScore(Cap *cap, stList *inducedAlignment1, stList *inducedAlignment2, stSortedSet *endAlignment1,
        stSortedSet *endAlignment2, void *capScoresFnHash) {

    int64_t i, j;
    int64_t *maxScore = st_malloc(sizeof(int64_t));
    maxScore[0] = getCutOff(inducedAlignment1, inducedAlignment2, &i, &j);
    assert(cap != NULL);
    assert(stHash_search(capScoresFnHash, cap) == NULL);
    stHash_insert(capScoresFnHash, cap, maxScore);
}

static int sortCapsFn(const void *cap1, const void *cap2, void *capScoresFnHash) {
    assert(stHash_search((stHash *)capScoresFnHash, (void *) cap1) != NULL);
    assert(stHash_search((stHash *)capScoresFnHash, (void *) cap2) != NULL);
    int64_t i = ((int64_t *) stHash_search((stHash *)capScoresFnHash, (void *) cap1))[0] - ((int64_t *) stHash_search(
            (stHash *)capScoresFnHash, (void *) cap2))[0];
    return (i > 0) ? 1 : ((i < 0) ? -1 : 0); 
}

bool isAlignedToStubSequence(AlignedPair *alignedPair, Flower *flower) {
	Cap *cap = flower_getCap(flower, alignedPair->reverse->subsequenceIdentifier);
    assert(cap != NULL);
    End *end1 = cap_getEnd(cap), *end2 = cap_getEnd(cap_getAdjacency(cap));
    assert(end1 != NULL && end2 != NULL);
    return (end_isStubEnd(end1) && end_isFree(end1)) || (end_isStubEnd(end2) && end_isFree(end2));
} 

static int64_t findFirstNonStubAlignment(Flower *flower, stList *inducedAlignment, bool reverse) {
    AlignedPair *pAlignedPair = NULL;
    int64_t j = -1;
    for (int64_t i = reverse ? stList_length(inducedAlignment) - 1 : 0; i < stList_length(inducedAlignment) && i >= 0; i
            += reverse ? -1 : 1) {
        AlignedPair *alignedPair = stList_get(inducedAlignment, i);
        assert(isAlignedToStubSequence(alignedPair->reverse, flower));
        assert(pAlignedPair == NULL || pAlignedPair->subsequenceIdentifier == alignedPair->subsequenceIdentifier);
        if (pAlignedPair == NULL || pAlignedPair->position != alignedPair->position) {
            pAlignedPair = alignedPair;
            j = i;
        }
        if(!isAlignedToStubSequence(alignedPair, flower)) {
            assert(j != -1);
            return j;
        }
    }
    return (reverse ? -1 : stList_length(inducedAlignment));
}

static void pruneStubAlignments(Cap *cap, stList *inducedAlignment1, stList *inducedAlignment2,
        stSortedSet *endAlignment1, stSortedSet *endAlignment2, void *deletedAlignedPairCounts) {
    assert(cap != NULL);
    End *end = cap_getEnd(cap);
    assert(cap_getAdjacency(cap) != NULL);
    End *adjacentEnd = cap_getEnd(cap_getAdjacency(cap));
    assert(end != NULL);
    assert(adjacentEnd != NULL);
    int64_t cutOff1 = stList_length(inducedAlignment1) - 1;
    int64_t cutOff2 = 0;
    if (end_isStubEnd(adjacentEnd) && end_isFree(adjacentEnd)) {
        cutOff1 = findFirstNonStubAlignment(end_getFlower(end), inducedAlignment1, 1);
        assert(stList_length(inducedAlignment2) == 0);
        cutOff2 = stList_length(inducedAlignment2);
    }
    if (end_isStubEnd(end) && end_isFree(end)) {
        assert(stList_length(inducedAlignment1) == 0);
        cutOff1 = -1;
        cutOff2 = findFirstNonStubAlignment(end_getFlower(end), inducedAlignment2, 0);
    }
    stSortedSet *pairsToDelete = stSortedSet_construct2((void(*)(void *)) alignedPair_destruct);
    //Now do the actual filtering of the alignments.
    pruneAlignmentsP(inducedAlignment1, endAlignment1, cutOff1 + 1, stList_length(inducedAlignment1), pairsToDelete, deletedAlignedPairCounts);
    pruneAlignmentsP(inducedAlignment2, endAlignment2, 0, cutOff2, pairsToDelete, deletedAlignedPairCounts);
    stSortedSet_destruct(pairsToDelete);
}

/*
 * Outer control functions that coordinate bar algorithm.
 */

static int makeFlowerAlignmentP(Cap *cap, stHash *endAlignments,
        void(*fn)(Cap *, stList *, stList *, stSortedSet *, stSortedSet *, void *), void *extraArg) {
    stSortedSet *endAlignment1 = stHash_search(endAlignments, end_getPositiveOrientation(cap_getEnd(cap)));
    assert(endAlignment1 != NULL);

    Cap *adjacentCap = cap_getAdjacency(cap);
    assert(adjacentCap != NULL);
    assert(cap_getSide(adjacentCap));
    assert(cap_getStrand(adjacentCap));
    adjacentCap = cap_getReverse(adjacentCap);
    stSortedSet *endAlignment2 = stHash_search(endAlignments, end_getPositiveOrientation(cap_getEnd(adjacentCap)));
    assert(endAlignment2 != NULL);

    AdjacencySequence *adjacencySequence1 = adjacencySequence_construct(cap, INT64_MAX);
    AdjacencySequence *adjacencySequence2 = adjacencySequence_construct(adjacentCap, INT64_MAX);
    assert(adjacencySequence1->length == adjacencySequence2->length);
    assert(adjacencySequence1->subsequenceIdentifier == adjacencySequence2->subsequenceIdentifier);
    assert(adjacencySequence1->strand == !adjacencySequence2->strand);
    assert(adjacencySequence2->start == adjacencySequence1->start + adjacencySequence1->length - 1);

    stList *inducedAlignment1 = getInducedAlignment(endAlignment1, adjacencySequence1);
    stList *inducedAlignment2 = getInducedAlignment(endAlignment2, adjacencySequence2);
    stList_reverse(inducedAlignment2);

    fn(cap, inducedAlignment1, inducedAlignment2, endAlignment1, endAlignment2, extraArg);

    //Cleanup.
    adjacencySequence_destruct(adjacencySequence1);
    adjacencySequence_destruct(adjacencySequence2);
    stList_destruct(inducedAlignment1);
    stList_destruct(inducedAlignment2);
    return 1;
}

static stSortedSet *makeFlowerAlignment2(Flower *flower, stHash *endAlignments, bool pruneOutStubAlignments) {
    /*
     * Makes the alignments of the ends, in "endAlignments", consistent with one another using the bar algorithm.
     */

    //Get the subsequences in the alignment that need to be pruned.
    End *end;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    stHash *capScoresFnHash = stHash_construct2(NULL, free);
    stList *caps = stList_construct();
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        Cap *cap;
        End_InstanceIterator *capIterator = end_getInstanceIterator(end);
        while ((cap = end_getNext(capIterator)) != NULL) {
            if (cap_getSide(cap)) {
                cap = cap_getReverse(cap);
            }
            if (cap_getStrand(cap)) {
                if (makeFlowerAlignmentP(cap, endAlignments, getScore, capScoresFnHash)) {
                    stList_append(caps, cap);
                }
            }
        }
        end_destructInstanceIterator(capIterator);
    }
    flower_destructEndIterator(endIterator);
    assert(stHash_size(capScoresFnHash) == stList_length(caps));
    stList_sort2(caps, sortCapsFn, capScoresFnHash); //sorts the caps in ascending order according to their cut off score.

    //Now do the actual pruning
    stHash *deletedAlignedPairCounts = stHash_construct3((uint64_t (*)(const void *))stIntTuple_hashKey,
            (int (*)(const void *, const void *))stIntTuple_equalsFn, (void (*)(void *))stIntTuple_destruct, free);
    stList *freeStubCaps = stList_construct(); //Caps that we'll use when pruning the stub only ends of alignments.
    while (stList_length(caps) > 0) {
        //Pick cap with greatest number of deleted aligned pairs.
        //This biases the bar algorithm to pick cutpoints that are consistent
        //with previously selected cutpoints.
        Cap *cap = NULL;
        int64_t deletedPairsForChosenCap = 0;
        int64_t capIndex = INT64_MAX;
        for(int64_t i=0; i<stList_length(caps); i++) {
            Cap *cap2 = stList_get(caps, i);
            assert(!cap_getSide(cap2));
            stIntTuple *sequenceIdentifier = stIntTuple_construct1(cap_getName(cap_getStrand(cap2) ? cap2 : cap_getAdjacency(cap2)));
            int64_t *deletedPairsCount = stHash_search(deletedAlignedPairCounts, sequenceIdentifier);
            stIntTuple_destruct(sequenceIdentifier);
            if(cap == NULL || deletedPairsForChosenCap == 0 || (deletedPairsCount != NULL && *deletedPairsCount >= deletedPairsForChosenCap)) {
                cap = cap2;
                capIndex = i;
                if(deletedPairsCount != NULL) {
                    deletedPairsForChosenCap = *deletedPairsCount;
                }
            }
        }
        //Having chosen the cap, remove it
        assert(cap != NULL);
        stList_remove(caps, capIndex);
        //Do the filtering.
        makeFlowerAlignmentP(cap, endAlignments, pruneAlignments, deletedAlignedPairCounts);
        assert(cap_getAdjacency(cap) != NULL);
        if ((end_isFree(cap_getEnd(cap)) && end_isStubEnd(cap_getEnd(cap))) || (end_isFree(
                cap_getEnd(cap_getAdjacency(cap))) && end_isStubEnd(cap_getEnd(cap_getAdjacency(cap))))) {
            stList_append(freeStubCaps, cap);
        } 
    }
    stList_destruct(caps);
    stHash_destruct(capScoresFnHash);

    if (pruneOutStubAlignments) { //This is used to remove matches only containing stub sequences at end of an end alignment.
    	while (stList_length(freeStubCaps) > 0) {
        	makeFlowerAlignmentP(stList_pop(freeStubCaps), endAlignments, pruneStubAlignments, deletedAlignedPairCounts);
        }
    }
    stList_destruct(freeStubCaps);

    //Now convert to set of final aligned pairs to return.
    stSortedSet *sortedAlignment = stSortedSet_construct3((int(*)(const void *, const void *)) alignedPair_cmpFn,
            (void(*)(void *)) alignedPair_destruct);
    stList *endAlignmentsList = stHash_getValues(endAlignments);
    while (stList_length(endAlignmentsList) > 0) {
        stSortedSet *endAlignment = stList_pop(endAlignmentsList);
        while (stSortedSet_size(endAlignment) > 0) {
            AlignedPair *alignedPair = stSortedSet_getFirst(endAlignment);
            stSortedSet_remove(endAlignment, alignedPair);
            assert(stSortedSet_search(endAlignment, alignedPair->reverse) != NULL);
            stSortedSet_remove(endAlignment, alignedPair->reverse);
            assert(stSortedSet_search(sortedAlignment, alignedPair) == NULL);
            assert(stSortedSet_search(sortedAlignment, alignedPair->reverse) == NULL);
            stSortedSet_insert(sortedAlignment, alignedPair);
            stSortedSet_insert(sortedAlignment, alignedPair->reverse);
        }
    }
    stList_destruct(endAlignmentsList);
    stHash_destruct(endAlignments);
    stHash_destruct(deletedAlignedPairCounts);

    return sortedAlignment;
}

/*
 * Functions to decide which ends need alignments.
 */

int64_t getMaxAdjacencyLength(Flower *flower) {
    assert(flower_getGroupNumber(flower) <= 1);
    Flower_CapIterator *capIt = flower_getCapIterator(flower);
    Cap *cap;
    int64_t maxAdjacencyLength = 0;
    while ((cap = flower_getNextCap(capIt)) != NULL) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        int64_t adjacencyLength = llabs(cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap)) - 1;
        assert(adjacencyLength >= 0);
        if (adjacencyLength > maxAdjacencyLength) {
            maxAdjacencyLength = adjacencyLength;
        }
    }
    flower_destructCapIterator(capIt);
    return maxAdjacencyLength;
}

End *getDominantEnd(Flower *flower) {
    /*
     * Returns an end, if exists, that has cap involved in every adjacency, else returns null.
     */
    //return NULL;
    assert(flower_getGroupNumber(flower) <= 1);
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    int64_t maxInstanceNumber = 0;
    End *dominantEnd = NULL;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_getInstanceNumber(end) > maxInstanceNumber) {
            maxInstanceNumber = end_getInstanceNumber(end);
            dominantEnd = end;
        }
    }
    flower_destructEndIterator(endIt);
    if (dominantEnd == NULL) {
        return NULL;
    }
    assert(end_getOrientation(dominantEnd));
    if (end_getInstanceNumber(dominantEnd) * 2 < flower_getCapNumber(flower)) {
        return NULL;
    }
    Cap *cap;
    Flower_CapIterator *capIt = flower_getCapIterator(flower);
    while ((cap = flower_getNextCap(capIt)) != NULL) {
        assert(cap_getAdjacency(cap) != NULL);
        if (end_getPositiveOrientation(cap_getEnd(cap)) != dominantEnd && end_getPositiveOrientation(
                cap_getEnd(cap_getAdjacency(cap))) != dominantEnd) {
            flower_destructCapIterator(capIt);
            return NULL;
        }
    }
    flower_destructCapIterator(capIt);
    return dominantEnd;
}

static stSortedSet *getEndsToAlign(Flower *flower, int64_t maxSequenceLength) {
    /*
     * Gets a set of the ends that we need to construct actual alignments for.
     */
    stSortedSet *endsToAlign = stSortedSet_construct();
    End *dominantEnd = getDominantEnd(flower); //an end to which all adjacencies are incident
    if (dominantEnd != NULL && getMaxAdjacencyLength(flower) <= 2 * maxSequenceLength) {
        stSortedSet_insert(endsToAlign, dominantEnd);
    } else {
        End *end;
        Flower_EndIterator *endIterator = flower_getEndIterator(flower);
        while ((end = flower_getNextEnd(endIterator)) != NULL) {
            stSortedSet_insert(endsToAlign, end);
        }
        flower_destructEndIterator(endIterator);
    }
    return endsToAlign;
}

/*
 * Functions that either create end alignments or load end alignments into memory from disk, and which
 * then call the makeFlowerAlignment2 consistency generating function.
 */

static void computeMissingEndAlignments(StateMachine *sM, Flower *flower, stHash *endAlignments, int64_t spanningTrees,
        int64_t maxSequenceLength, bool useProgressiveMerging, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters,
        int64_t poaWindow) {
    /*
     * Creates end alignments for the ends that
     * do not have an alignment in the "endAlignments" hash, only creating
     * non-trivial end alignments for those specified by "getEndsToAlign".
     */
    //Make the end alignments, representing each as an adjacency alignment.
    stSortedSet *endsToAlign = getEndsToAlign(flower, maxSequenceLength);
    End *end;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        if (stHash_search(endAlignments, end) == NULL) {
            if (stSortedSet_search(endsToAlign, end) != NULL) {
                stHash_insert(
                        endAlignments,
                        end,
                        makeEndAlignment(sM, end, spanningTrees, maxSequenceLength,
                                useProgressiveMerging, gapGamma,
                                pairwiseAlignmentBandingParameters, poaWindow));
            } else {
                stHash_insert(endAlignments, end, stSortedSet_construct());
            }
        }
    }
    flower_destructEndIterator(endIterator);
    stSortedSet_destruct(endsToAlign);
}

stSortedSet *makeFlowerAlignment(StateMachine *sM, Flower *flower, int64_t spanningTrees, int64_t maxSequenceLength,
        bool useProgressiveMerging, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters, bool pruneOutStubAlignments, int64_t poaWindow) {
    stHash *endAlignments = stHash_construct2(NULL, (void(*)(void *)) stSortedSet_destruct);
    computeMissingEndAlignments(sM, flower, endAlignments, spanningTrees, maxSequenceLength,
            useProgressiveMerging, gapGamma, pairwiseAlignmentBandingParameters, poaWindow);
    return makeFlowerAlignment2(flower, endAlignments, pruneOutStubAlignments);
}

static void loadEndAlignments(Flower *flower, stHash *endAlignments, stList *listOfEndAlignments) {
    /*
     * Load alignments from given list of files and add them to the "endAlignments" hash.
     */
    for (int64_t i = 0; i < stList_length(listOfEndAlignments); i++) {
        End *end;
        FILE *fileHandle = fopen(stList_get(listOfEndAlignments, i), "r");
        stSortedSet *alignment;
        while((alignment = loadEndAlignmentFromDisk(flower, fileHandle, &end)) != NULL) {
            assert(stHash_search(endAlignments, end) == NULL);
            stHash_insert(endAlignments, end, alignment);
        }
        fclose(fileHandle);
    }
}

stSortedSet *makeFlowerAlignment3(StateMachine *sM, Flower *flower, stList *listOfEndAlignmentFiles, int64_t spanningTrees,
        int64_t maxSequenceLength, bool useProgressiveMerging, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters, bool pruneOutStubAlignments, int64_t poaWindow) {
    stHash *endAlignments = stHash_construct2(NULL, (void(*)(void *)) stSortedSet_destruct);
    if(listOfEndAlignmentFiles != NULL) {
        loadEndAlignments(flower, endAlignments, listOfEndAlignmentFiles);
    }
    computeMissingEndAlignments(sM, flower, endAlignments, spanningTrees, maxSequenceLength,
            useProgressiveMerging, gapGamma, pairwiseAlignmentBandingParameters, poaWindow);
    return makeFlowerAlignment2(flower, endAlignments, pruneOutStubAlignments);
}

/*
 * Functions for calculating large end alignments that should be computed separately for parallelism.
 */

int64_t getTotalAdjacencyLength(End *end) {
    /*
     * Gets the total length of unaligned sequences on adjacencies
     * incident with the instances of the end.
     */
    End_InstanceIterator *capIt = end_getInstanceIterator(end);
    Cap *cap;
    int64_t totalAdjacencyLength = 0;
    while ((cap = end_getNext(capIt)) != NULL) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        totalAdjacencyLength += llabs(cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap)) - 1;
    }
    end_destructInstanceIterator(capIt);
    return totalAdjacencyLength;
}

stSortedSet *getEndsToAlignSeparately(Flower *flower, int64_t maxSequenceLength, int64_t largeEndSize) {
    /*
     * Picks a set of end alignments that contain more than "largeEndSize" bases and, if there are more
     * than 2 of them, returns them in a set.
     */
    stSortedSet *endsToAlign = getEndsToAlign(flower, maxSequenceLength);
    stSortedSetIterator *it = stSortedSet_getIterator(endsToAlign);
    End *end;
    stSortedSet *largeEndsToAlign = stSortedSet_construct();
    while ((end = stSortedSet_getNext(it)) != NULL) {
        if (getTotalAdjacencyLength(end) >= largeEndSize) {
            stSortedSet_insert(largeEndsToAlign, end);
        }
    }
    stSortedSet_destruct(endsToAlign);
    if (stSortedSet_size(largeEndsToAlign) <= 1) {
        stSortedSet_destruct(largeEndsToAlign);
        return stSortedSet_construct();
    }
    return largeEndsToAlign;
}
