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

/*
 * Following functions deal with trimming the ends of stub alignments.
 */
static stSortedSet *getFreeStubAdjacencySequences(Flower *flower) {
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    stSortedSet *freeStubAdjacencySequences = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn, (void(*)(void *)) stIntTuple_destruct);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        if (end_isFree(end) && end_isStubEnd(end)) {
            Cap *cap;
            End_InstanceIterator *capIterator = end_getInstanceIterator(end);
            while ((cap = end_getNext(capIterator)) != NULL) {
                if (cap_getAdjacency(cap) != NULL) {
                    cap = cap_getSide(cap) ? cap_getReverse(cap) : cap;
                    AdjacencySequence *aS = adjacencySequence_construct(cap, INT64_MAX);
                    if(cap_getStrand(cap)) {
                        stSortedSet_insert(freeStubAdjacencySequences, stIntTuple_construct3(aS->subsequenceIdentifier, aS->start, aS->start + aS->length));
                    }
                    else { //The start coordinate is the end of the subsequence, kind of bizarrely.
                        stSortedSet_insert(freeStubAdjacencySequences, stIntTuple_construct3(aS->subsequenceIdentifier, aS->start - aS->length + 1, aS->start + 1));
                    }
                    adjacencySequence_destruct(aS);
                } else {
                    assert(cap_getSequence(cap) == NULL);
                }
            }
            end_destructInstanceIterator(capIterator);
        }
    }
    flower_destructEndIterator(endIterator);
    return freeStubAdjacencySequences;
}

static bool isAlignedOnlyToStubSequence(AlignedPair *alignedPair, stSet *stubColumns) {
    stIntTuple *k = stIntTuple_construct1(alignedPair->columnId);
    bool b = stSet_search(stubColumns, k) != NULL;
    stIntTuple_destruct(k);
    return b;
}

stList *getInducedAlignment(stSortedSet *endAlignment, AdjacencySequence *adjacencySequence) {
    /*
     * Gets an ordered list of pairs from the end alignment for the given adjacency sequence.
     */
    stList *inducedAlignment = stList_construct();
    if (adjacencySequence->strand) {
        AlignedPair *alignedPair = alignedPair_construct(adjacencySequence->subsequenceIdentifier,
                adjacencySequence->start - 1, adjacencySequence->strand, 0, 0, 0);
        AlignedPair *alignedPair2 = stSortedSet_searchGreaterThan(endAlignment, alignedPair);
        if (alignedPair2 != NULL) {
            stSortedSetIterator *it = stSortedSet_getIteratorFrom(endAlignment, alignedPair2);
            while ((alignedPair2 = stSortedSet_getNext(it)) != NULL) {
                if (alignedPair2->subsequenceIdentifier == adjacencySequence->subsequenceIdentifier) {
                    if (alignedPair2->position >= adjacencySequence->start + adjacencySequence->length) {
                        break;
                    }
                    assert(alignedPair2->position >= adjacencySequence->start);
                    assert(alignedPair2->strand == adjacencySequence->strand);
                    stList_append(inducedAlignment, alignedPair2);
                } else {
                    break;
                }
            }
            stSortedSet_destructIterator(it);
        }
        alignedPair_destruct(alignedPair);
    } else {
        AlignedPair *alignedPair = alignedPair_construct(adjacencySequence->subsequenceIdentifier,
                adjacencySequence->start + 1, adjacencySequence->strand, 0, 0, 0);
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
                    assert(alignedPair2->strand == adjacencySequence->strand);
                    stList_append(inducedAlignment, alignedPair2);
                } else {
                    break;
                }
            }
            stSortedSet_destructIterator(it);
        }
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
static int64_t *cumulateScoreForward(stList *inducedAlignment1, bool scoreWithoutStubs) {
    int64_t *iA = st_malloc(sizeof(int64_t) * stList_length(inducedAlignment1));
    int64_t totalScore = 0;
    for (int64_t i = 0; i < stList_length(inducedAlignment1); i++) {
        AlignedPair *alignedPair = stList_get(inducedAlignment1, i);
        totalScore += scoreWithoutStubs ? alignedPair->scoreWithoutStubs : alignedPair->score;
        iA[i] = totalScore;
    }
    return iA;
}

/*
 * Runs along and cumulate the score of the pairs, traversing backward through the induced alignment.
 */
static int64_t *cumulateScoreBackward(stList *inducedAlignment1,  bool scoreWithoutStubs) {
    int64_t *iA = st_malloc(sizeof(int64_t) * stList_length(inducedAlignment1));
    int64_t totalScore = 0;
    for (int64_t i = stList_length(inducedAlignment1) - 1; i >= 0; i--) {
        AlignedPair *alignedPair = stList_get(inducedAlignment1, i);
        totalScore += scoreWithoutStubs ? alignedPair->scoreWithoutStubs : alignedPair->score;
        iA[i] = totalScore;
    }
    return iA;
}

/*
 * Chooses a point along the adjacency sequence at which to filter the two alignments,
 */
static int64_t getCutOff(stList *inducedAlignment1, stList *inducedAlignment2, int64_t *cutOff1, int64_t *cutOff2,
        bool excludeFreeStubAdjacencySequences) {
    int64_t *cScore1, *cScore1B;
    int64_t *cScore2, *cScore2B;
    if(excludeFreeStubAdjacencySequences) { //We calculate arrays with and without including alignments to stub sequences
        cScore1 = cumulateScoreForward(inducedAlignment1, 1);
        cScore2 = cumulateScoreBackward(inducedAlignment2, 1);
        cScore1B = cumulateScoreForward(inducedAlignment1, 0);
        cScore2B = cumulateScoreBackward(inducedAlignment2, 0);
    }
    else {
        cScore1 = cumulateScoreForward(inducedAlignment1, 0);
        cScore2 = cumulateScoreBackward(inducedAlignment2, 0);
        cScore1B = cScore1;
        cScore2B = cScore2;
    }

    //Check the score arrays for sanity..
    for (int64_t i = 1; i < stList_length(inducedAlignment1); i++) {
        assert(cScore1B[i - 1] < cScore1B[i]);
    }
    for (int64_t i = 1; i < stList_length(inducedAlignment2); i++) {
        assert(cScore2B[i - 1] > cScore2B[i]);
    }

    //Find the score threshold by iterating along each alignment.
    *cutOff1 = 0;
    *cutOff2 = 0;
    int64_t maxScore = -1, maxScoreB = -1;
    if (stList_length(inducedAlignment2) > 0) {
        maxScore = cScore2[0];
        maxScoreB = cScore2B[0];
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
                    if ((cScore1[i] + cScore2[j] > maxScore) || (cScore1[i] + cScore2[j] == maxScore && cScore1B[i] + cScore2B[j] >= maxScoreB)) {
                        maxScore = cScore1[i] + cScore2[j];
                        maxScoreB = cScore1B[i] + cScore2B[j];
                        *cutOff1 = i + 1;
                        *cutOff2 = j;
                    }
                    break;
                } else {
                    j++;
                }
            } while (j < stList_length(inducedAlignment2));
        } else {
            if (cScore1[i] > maxScore || (cScore1[i] == maxScore && cScore1B[i] >= maxScoreB)) {
                *cutOff1 = stList_length(inducedAlignment1);
                *cutOff2 = j;
                assert(cScore1[stList_length(inducedAlignment1) - 1] >= maxScore);
                maxScore = cScore1[stList_length(inducedAlignment1) - 1];
                maxScoreB = cScore1B[stList_length(inducedAlignment1) - 1];
                break;
            }
        }
    }
    //Cleanup
    free(cScore1);
    free(cScore2);
    if(excludeFreeStubAdjacencySequences) {
        free(cScore1B);
        free(cScore2B);
    }

    return maxScoreB;
}

static void pruneAlignmentsP(stList *inducedAlignment, stSortedSet *endAlignment, int64_t start, int64_t end) {
    for (int64_t i = start; i < end; i++) {
        AlignedPair *alignedPair = stList_get(inducedAlignment, i);
        assert(stSortedSet_search(endAlignment, alignedPair) != NULL);
        stSortedSet_remove(endAlignment, alignedPair);
        alignedPair_destruct(alignedPair);
    }
}

static void pruneAlignments(Cap *cap, stList *inducedAlignment1, stList *inducedAlignment2, stSortedSet *endAlignment1,
        stSortedSet *endAlignment2, void *extraArg) {
    /*
     * Chooses a point along the adjacency sequence at which to filter the two alignments,
     * then filters the aligned pairs by this point.
     */
    int64_t cutOff1 = 0, cutOff2 = 0;
    getCutOff(inducedAlignment1, inducedAlignment2, &cutOff1, &cutOff2, 1);
    //Now do the actual filtering of the alignments.
    pruneAlignmentsP(inducedAlignment1, endAlignment1, cutOff1, stList_length(inducedAlignment1));
    pruneAlignmentsP(inducedAlignment2, endAlignment2, 0, cutOff2);
}

static stHash *capScoresFn_Hash = NULL;
void getScore(Cap *cap, stList *inducedAlignment1, stList *inducedAlignment2, stSortedSet *endAlignment1,
        stSortedSet *endAlignment2, void *extraArg) {

    int64_t i, j;
    int64_t *maxScore = st_malloc(sizeof(int64_t));
    maxScore[0] = getCutOff(inducedAlignment1, inducedAlignment2, &i, &j, 0);
    assert(cap != NULL);
    assert(stHash_search(capScoresFn_Hash, cap) == NULL);
    stHash_insert(capScoresFn_Hash, cap, maxScore);
}

static int sortCapsFn(const void *cap1, const void *cap2) {
    assert(stHash_search(capScoresFn_Hash, (void *) cap1) != NULL);
    assert(stHash_search(capScoresFn_Hash, (void *) cap2) != NULL);
    int64_t i = ((int64_t *) stHash_search(capScoresFn_Hash, (void *) cap1))[0] - ((int64_t *) stHash_search(
            capScoresFn_Hash, (void *) cap2))[0];
    return (i > 0) ? 1 : ((i < 0) ? -1 : 0);
}

static int64_t findFirstNonStubAlignment(stList *inducedAlignment, bool reverse,
        stSet *stubColumns) {
    AlignedPair *pAlignedPair = NULL;
    int64_t j = -1;
    for (int64_t i = reverse ? stList_length(inducedAlignment) - 1 : 0; i < stList_length(inducedAlignment) && i >= 0; i
            += reverse ? -1 : 1) {
        AlignedPair *alignedPair = stList_get(inducedAlignment, i);
        assert(pAlignedPair == NULL || pAlignedPair->subsequenceIdentifier == alignedPair->subsequenceIdentifier);
        if (pAlignedPair == NULL || pAlignedPair->position != alignedPair->position) {
            pAlignedPair = alignedPair;
            j = i;
        }
        if(!isAlignedOnlyToStubSequence(alignedPair, stubColumns)) {
            assert(j != -1);
            return j;
        }
    }
    return (reverse ? -1 : stList_length(inducedAlignment));
}

static void pruneStubAlignments(Cap *cap, stList *inducedAlignment1, stList *inducedAlignment2,
        stSortedSet *endAlignment1, stSortedSet *endAlignment2, void *stubColumns) {
    assert(cap != NULL);
    End *end = cap_getEnd(cap);
    assert(cap_getAdjacency(cap) != NULL);
    End *adjacentEnd = cap_getEnd(cap_getAdjacency(cap));
    assert(end != NULL);
    assert(adjacentEnd != NULL);
    int64_t cutOff1 = stList_length(inducedAlignment1) - 1;
    int64_t cutOff2 = 0;
    if (end_isStubEnd(adjacentEnd) && end_isFree(adjacentEnd)) {
        cutOff1 = findFirstNonStubAlignment(inducedAlignment1, 1, stubColumns);
        assert(stList_length(inducedAlignment2) == 0);
        cutOff2 = stList_length(inducedAlignment2);
    }
    if (end_isStubEnd(end) && end_isFree(end)) {
        assert(stList_length(inducedAlignment1) == 0);
        cutOff1 = -1;
        cutOff2 = findFirstNonStubAlignment(inducedAlignment2, 0, stubColumns);
    }
    //Now do the actual filtering of the alignments.
    pruneAlignmentsP(inducedAlignment1, endAlignment1, cutOff1 + 1, stList_length(inducedAlignment1));
    pruneAlignmentsP(inducedAlignment2, endAlignment2, 0, cutOff2);
}

/*
 * Outer control functions that coordinate bar algorithm.
 */

static int makeFlowerAlignmentP(Cap *cap, stHash *endAlignments,
        void(*fn)(Cap *, stList *, stList *, stSortedSet *, stSortedSet *, void *),
        void *extraArg) {
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

static uint64_t alignedPair_columnFn(const void *a) {
    return ((AlignedPair *)a)->columnId;
}


static int alignedPair_columnEqualsFn(const void *a, const void *b) {
    return ((AlignedPair *)a)->columnId == ((AlignedPair *)b)->columnId;
}

static stSortedSet *makeFlowerAlignment2(Flower *flower, stHash *endAlignments, bool pruneOutStubAlignments) {
    /*
     * Makes the alignments of the ends, in "endAlignments", consistent with one another using the bar algorithm.
     */

    //Prune end alignments
    End *end;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    capScoresFn_Hash = stHash_construct2(NULL, free);
    stList *caps = stList_construct();
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        Cap *cap;
        End_InstanceIterator *capIterator = end_getInstanceIterator(end);
        while ((cap = end_getNext(capIterator)) != NULL) {
            if (cap_getSide(cap)) {
                cap = cap_getReverse(cap);
            }
            if (cap_getStrand(cap)) {
                if (makeFlowerAlignmentP(cap, endAlignments, getScore, NULL)) {
                    stList_append(caps, cap);
                }
            }
        }
        end_destructInstanceIterator(capIterator);
    }
    flower_destructEndIterator(endIterator);
    assert(stHash_size(capScoresFn_Hash) == stList_length(caps));
    stList_sort(caps, sortCapsFn);
    int64_t score = INT64_MAX;
    stList *freeStubCaps = stList_construct();
    while (stList_length(caps) > 0) {
        Cap *cap = stList_pop(caps);
        int64_t score2 = ((int64_t *) stHash_search(capScoresFn_Hash, cap))[0];
        assert(score2 <= score);
        score = score2;

        makeFlowerAlignmentP(cap, endAlignments, pruneAlignments, NULL);
        assert(cap_getAdjacency(cap) != NULL);
        if ((end_isFree(cap_getEnd(cap)) && end_isStubEnd(cap_getEnd(cap))) || (end_isFree(
                cap_getEnd(cap_getAdjacency(cap))) && end_isStubEnd(cap_getEnd(cap_getAdjacency(cap))))) {
            stList_append(freeStubCaps, cap);
        }
    }
    stList_destruct(caps);
    stHash_destruct(capScoresFn_Hash);

    if (pruneOutStubAlignments) { //This is used to remove matches only containing stub sequences at end of an end alignment.
        //Get sequences that are stubs.
        stSortedSet *freeStubAdjacencySequences = getFreeStubAdjacencySequences(flower);

        //Build hash of columns only containing alignments between stubs.
        stSet *stubColumns = stSet_construct();

        while (stList_length(freeStubCaps) > 0) {
        	makeFlowerAlignmentP(stList_pop(freeStubCaps), endAlignments, pruneStubAlignments, stubColumns);
        }

        stSortedSet_destruct(freeStubAdjacencySequences);
        stSet_destruct(stubColumns);
    }
    stList_destruct(freeStubCaps);

    //Now convert to set of final aligned pairs to return.
    stSortedSet *sortedAlignment = stSortedSet_construct3((int(*)(const void *, const void *)) alignedPair_cmpFn,
            (void(*)(void *)) alignedPair_destruct);
    stHashIterator *endAlignmentsIt = stHash_getIterator(endAlignments);
    while ((end = stHash_getNext(endAlignmentsIt)) != NULL) {
        stSortedSetIterator *endAlignmentIt = stSortedSet_getIterator(stHash_search(endAlignments, end));
        stSet *alignedPairsToColumns = stSet_construct3(alignedPair_columnFn, alignedPair_columnEqualsFn, NULL);
        AlignedPair *alignedPair;
        while ((alignedPair = stSortedSet_getNext(endAlignmentIt)) != NULL) {
            AlignedPair *alignedPair2 = stSet_search(alignedPairsToColumns, alignedPair);
            if(alignedPair2 == NULL) {
                stSet_insert(alignedPairsToColumns, alignedPair);
            }
            else {
                assert(alignedPair != alignedPair2);
                stSortedSet_insert(sortedAlignment, stIntTuple_construct5(alignedPair->subsequenceIdentifier, alignedPair->position,
                            alignedPair2->subsequenceIdentifier, alignedPair2->position, alignedPair->strand == alignedPair2->strand));
            }
        }
        stSortedSet_destructIterator(endAlignmentIt);
        stSet_destruct(alignedPairsToColumns);
    }
    stHash_destructIterator(endAlignmentsIt);

    //cleanup
    stHash_destruct(endAlignments);

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

static void computeMissingEndAlignments(Flower *flower, stHash *endAlignments, int64_t spanningTrees,
        int64_t maxSequenceLength, int64_t maximumNumberOfSequencesBeforeSwitchingToFast, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
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
                        makeEndAlignment(end, spanningTrees, maxSequenceLength,
                                maximumNumberOfSequencesBeforeSwitchingToFast, gapGamma,
                                pairwiseAlignmentBandingParameters));
            } else {
                stHash_insert(endAlignments, end, stSortedSet_construct());
            }
        }
    }
    flower_destructEndIterator(endIterator);
    stSortedSet_destruct(endsToAlign);
}


stSortedSet *makeFlowerAlignment(Flower *flower, int64_t spanningTrees, int64_t maxSequenceLength,
        int64_t maximumNumberOfSequencesBeforeSwitchingToFast, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters, bool pruneOutStubAlignments) {
    stHash *endAlignments = stHash_construct2(NULL, (void(*)(void *)) stSortedSet_destruct);
    computeMissingEndAlignments(flower, endAlignments, spanningTrees, maxSequenceLength,
            maximumNumberOfSequencesBeforeSwitchingToFast, gapGamma, pairwiseAlignmentBandingParameters);
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

stSortedSet *makeFlowerAlignment3(Flower *flower, stList *listOfEndAlignmentFiles, int64_t spanningTrees,
        int64_t maxSequenceLength, int64_t maximumNumberOfSequencesBeforeSwitchingToFast, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters, bool pruneOutStubAlignments) {
    stHash *endAlignments = stHash_construct2(NULL, (void(*)(void *)) stSortedSet_destruct);
    if(listOfEndAlignmentFiles != NULL) {
        loadEndAlignments(flower, endAlignments, listOfEndAlignmentFiles);
    }
    computeMissingEndAlignments(flower, endAlignments, spanningTrees, maxSequenceLength,
            maximumNumberOfSequencesBeforeSwitchingToFast, gapGamma, pairwiseAlignmentBandingParameters);
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
