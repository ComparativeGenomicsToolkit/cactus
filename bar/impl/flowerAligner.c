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
 * Gets an ordered list of pairs from the end alignment for the given adjacency sequence.
 */
stList *getInducedAlignment(stSortedSet *endAlignment,
        AdjacencySequence *adjacencySequence) {
    stList *inducedAlignment = stList_construct();
    if (adjacencySequence->strand) {
        AlignedPair *alignedPair = alignedPair_construct(
                adjacencySequence->subsequenceIdentifier, adjacencySequence->start - 1,
                adjacencySequence->strand, 0, 0, 0, 0);
        AlignedPair *alignedPair2 = stSortedSet_searchGreaterThan(endAlignment,
                alignedPair);
        if (alignedPair2 != NULL) {
            stSortedSetIterator *it = stSortedSet_getIteratorFrom(endAlignment,
                    alignedPair2);
            while ((alignedPair2 = stSortedSet_getNext(it)) != NULL) {
                if (alignedPair2->subsequenceIdentifier == adjacencySequence->subsequenceIdentifier) {
                    if (alignedPair2->position >= adjacencySequence->start
                            + adjacencySequence->length) {
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
        AlignedPair *alignedPair = alignedPair_construct(
                adjacencySequence->subsequenceIdentifier, adjacencySequence->start + 1,
                adjacencySequence->strand, 0, 0, 0, 0);
        AlignedPair *alignedPair2 = stSortedSet_searchLessThan(endAlignment,
                alignedPair);
        if (alignedPair2 != NULL) {
            stSortedSetIterator *it = stSortedSet_getIteratorFrom(endAlignment,
                    alignedPair2);
            stSortedSet_getNext(it);
            stSortedSet_getNext(it); //shift it up..
            while ((alignedPair2 = stSortedSet_getPrevious(it)) != NULL) {
                if (alignedPair2->subsequenceIdentifier == adjacencySequence->subsequenceIdentifier) {
                    if (alignedPair2->position <= adjacencySequence->start
                            - adjacencySequence->length) {
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
    for (int32_t i = 0; i < stList_length(inducedAlignment); i++) {
        AlignedPair *alignedPair = stList_get(inducedAlignment, i);
        (void)alignedPair;
        assert(alignedPair->subsequenceIdentifier == adjacencySequence->subsequenceIdentifier);
        assert(alignedPair->strand == adjacencySequence->strand);
        if (adjacencySequence->strand) {
            assert(alignedPair->position >= adjacencySequence->start);
            assert(
                    alignedPair->position < adjacencySequence->start
                            + adjacencySequence->length);
        } else {
            assert(alignedPair->position <= adjacencySequence->start);
            assert(
                    alignedPair->position > adjacencySequence->start
                            - adjacencySequence->length);
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
    for (int32_t i = 0; i < stList_length(inducedAlignment1); i++) {
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
    for (int32_t i = stList_length(inducedAlignment1) - 1; i >= 0; i--) {
        AlignedPair *alignedPair = stList_get(inducedAlignment1, i);
        totalScore += alignedPair->score;
        iA[i] = totalScore;
    }
    return iA;
}

/*
 * Chooses a point along the adjacency sequence at which to filter the two alignments,
 */
static int64_t getCutOff(stList *inducedAlignment1, stList *inducedAlignment2,
        int32_t *cutOff1, int32_t *cutOff2) {
    int64_t *cScore1 = cumulateScoreForward(inducedAlignment1);
    int64_t *cScore2 = cumulateScoreBackward(inducedAlignment2);

    //Check the score arrays for sanity..
    for (int32_t i = 1; i < stList_length(inducedAlignment1); i++) {
        assert(cScore1[i - 1] < cScore1[i]);
    }
    for (int32_t i = 1; i < stList_length(inducedAlignment2); i++) {
        assert(cScore2[i - 1] > cScore2[i]);
    }

    //Find the score threshold by iterating along each alignment.
    *cutOff1 = 0;
    *cutOff2 = 0;
    int64_t maxScore = -1;
    if (stList_length(inducedAlignment2) > 0) {
        maxScore = cScore2[0];
    }
    int32_t j = 0;
    int32_t pPos1 = INT32_MIN, pPos2 = INT32_MIN;
    for (int32_t i = 0; i < stList_length(inducedAlignment1); i++) {
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
                assert(
                        cScore1[stList_length(inducedAlignment1) - 1]
                                >= maxScore);
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

static void pruneAlignmentsP(stList *inducedAlignment,
        stSortedSet *endAlignment, int32_t start, int32_t end,
        stSortedSet *pairsToDelete) {
    for (int32_t i = start; i < end; i++) {
        AlignedPair *alignedPair = stList_get(inducedAlignment, i);
        if (stSortedSet_search(endAlignment, alignedPair) != NULL) { //can be missing if we are pruning the reverse strand alignment at the same time
            assert(
                    stSortedSet_search(endAlignment, alignedPair->reverse)
                            != NULL);
            stSortedSet_remove(endAlignment, alignedPair);
            stSortedSet_remove(endAlignment, alignedPair->reverse);
            if (stSortedSet_search(pairsToDelete, alignedPair) == NULL) { // &&
                assert(
                        stSortedSet_search(pairsToDelete, alignedPair->reverse)
                                == NULL);
                stSortedSet_insert(pairsToDelete, alignedPair);
                stSortedSet_insert(pairsToDelete, alignedPair->reverse);
            } else {
                assert(
                        stSortedSet_search(pairsToDelete, alignedPair->reverse)
                                != NULL);
            }
        }
    }
}

/*
 * Chooses a point along the adjacency sequence at which to filter the two alignments,
 * then filters the aligned pairs by this point.
 */
static void pruneAlignments(Cap *cap, stList *inducedAlignment1,
        stList *inducedAlignment2, stSortedSet *endAlignment1,
        stSortedSet *endAlignment2) {
    int32_t cutOff1 = 0, cutOff2 = 0;
    getCutOff(inducedAlignment1, inducedAlignment2, &cutOff1, &cutOff2);
    stSortedSet *pairsToDelete = stSortedSet_construct2(
            (void(*)(void *)) alignedPair_destruct);
    //Now do the actual filtering of the alignments.
    pruneAlignmentsP(inducedAlignment1, endAlignment1, cutOff1,
            stList_length(inducedAlignment1), pairsToDelete);
    pruneAlignmentsP(inducedAlignment2, endAlignment2, 0, cutOff2,
            pairsToDelete);
    stSortedSet_destruct(pairsToDelete);
}

/*
 * Following functions deal with trimming the ends of stub alignments.
 */
static stSortedSet *getFreeStubAdjacencySequences(Flower *flower) {
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    stSortedSet *freeStubAdjacencySequences = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stInt64Tuple_cmpFn,
            (void(*)(void *)) stInt64Tuple_destruct);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        if (end_isFree(end) && end_isStubEnd(end)) {
            Cap *cap;
            End_InstanceIterator *capIterator = end_getInstanceIterator(end);
            while ((cap = end_getNext(capIterator)) != NULL) {
                if (cap_getAdjacency(cap) != NULL) {
                    cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                    Cap *adjacentCap = cap_getAdjacency(cap);
                    assert(
                            cap_getCoordinate(cap) != cap_getCoordinate(
                                    adjacentCap));
                    if (cap_getCoordinate(cap) > cap_getCoordinate(adjacentCap)) {
                        Cap *cap2 = cap;
                        cap = adjacentCap;
                        adjacentCap = cap2;
                    }
                    assert(
                            cap_getCoordinate(cap) < cap_getCoordinate(
                                    adjacentCap));
                    assert(cap_getSequence(cap) != NULL);
                    assert(!cap_getSide(cap));
                    assert(cap_getSide(adjacentCap));
                    Sequence *sequence = cap_getSequence(cap);
                    stSortedSet_insert(
                            freeStubAdjacencySequences,
                            stInt64Tuple_construct(3,
                                    sequence_getName(sequence),
                                    cap_getCoordinate(cap) + 1,
                                    cap_getCoordinate(adjacentCap)));
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

static int32_t findFirstNonStubAlignment(stList *inducedAlignment,
        bool reverse, stSortedSet *freeStubAdjacencySequences) {
    AlignedPair *pAlignedPair = NULL;
    int32_t j = -1;
    for (int32_t i = reverse ? stList_length(inducedAlignment) - 1 : 0; i
            < stList_length(inducedAlignment) && i >= 0; i += reverse ? -1 : 1) {
        AlignedPair *alignedPair = stList_get(inducedAlignment, i);
        assert(
                pAlignedPair == NULL || pAlignedPair->subsequenceIdentifier
                        == alignedPair->subsequenceIdentifier);
        if (pAlignedPair == NULL || pAlignedPair->position
                != alignedPair->position) {
            pAlignedPair = alignedPair;
            j = i;
        }
        stInt64Tuple *k = stInt64Tuple_construct(3,
                alignedPair->reverse->subsequenceIdentifier,
                (int64_t) alignedPair->reverse->position, INT64_MAX);
        stInt64Tuple *l = stSortedSet_searchLessThanOrEqual(
                freeStubAdjacencySequences, k);
        stInt64Tuple_destruct(k);
        if (l != NULL && alignedPair->reverse->subsequenceIdentifier
                == stInt64Tuple_getPosition(l, 0)
                && alignedPair->reverse->position < stInt64Tuple_getPosition(l,
                        2)) {
            assert(
                    stInt64Tuple_getPosition(l, 1)
                            <= alignedPair->reverse->position);
            continue;
        }
        assert(j != -1);
        return j;
    }
    return (reverse ? -1 : stList_length(inducedAlignment));
}

static stSortedSet *pruneStubAlignments_freeStubAdjacencySequences = NULL;
static void pruneStubAlignments(Cap *cap, stList *inducedAlignment1,
        stList *inducedAlignment2, stSortedSet *endAlignment1,
        stSortedSet *endAlignment2) {
    assert(cap != NULL);
    End *end = cap_getEnd(cap);
    assert(cap_getAdjacency(cap) != NULL);
    End *adjacentEnd = cap_getEnd(cap_getAdjacency(cap));
    assert(end != NULL);
    assert(adjacentEnd != NULL);
    int32_t cutOff1 = stList_length(inducedAlignment1) - 1;
    int32_t cutOff2 = 0;
    if (end_isStubEnd(adjacentEnd) && end_isFree(adjacentEnd)) {
        cutOff1 = findFirstNonStubAlignment(inducedAlignment1, 1,
                pruneStubAlignments_freeStubAdjacencySequences);
        assert(stList_length(inducedAlignment2) == 0);
        cutOff2 = stList_length(inducedAlignment2);
    }
    if (end_isStubEnd(end) && end_isFree(end)) {
        assert(stList_length(inducedAlignment1) == 0);
        cutOff1 = -1;
        cutOff2 = findFirstNonStubAlignment(inducedAlignment2, 0,
                pruneStubAlignments_freeStubAdjacencySequences);
    }
    stSortedSet *pairsToDelete = stSortedSet_construct2(
            (void(*)(void *)) alignedPair_destruct);
    //Now do the actual filtering of the alignments.
    pruneAlignmentsP(inducedAlignment1, endAlignment1, cutOff1 + 1,
            stList_length(inducedAlignment1), pairsToDelete);
    pruneAlignmentsP(inducedAlignment2, endAlignment2, 0, cutOff2,
            pairsToDelete);
    stSortedSet_destruct(pairsToDelete);
}

static stHash *capScoresFn_Hash = NULL;
void getScore(Cap *cap, stList *inducedAlignment1, stList *inducedAlignment2,
        stSortedSet *endAlignment1, stSortedSet *endAlignment2) {

    int32_t i, j;
    int64_t *maxScore = st_malloc(sizeof(int64_t));
    maxScore[0] = getCutOff(inducedAlignment1, inducedAlignment2, &i, &j);
    assert(cap != NULL);
    assert(stHash_search(capScoresFn_Hash, cap) == NULL);
    stHash_insert(capScoresFn_Hash, cap, maxScore);
}

static int sortCapsFn(const void *cap1, const void *cap2) {
    assert(stHash_search(capScoresFn_Hash, (void *) cap1) != NULL);
    assert(stHash_search(capScoresFn_Hash, (void *) cap2) != NULL);
    int64_t i = ((int64_t *) stHash_search(capScoresFn_Hash, (void *) cap1))[0]
            - ((int64_t *) stHash_search(capScoresFn_Hash, (void *) cap2))[0];
    return (i > 0) ? 1 : ((i < 0) ? -1 : 0);
}

static int makeFlowerAlignmentP(Cap *cap, stHash *endAlignments,
        void(*fn)(Cap *, stList *, stList *, stSortedSet *, stSortedSet *)) {
    stSortedSet *endAlignment1 = stHash_search(endAlignments,
            end_getPositiveOrientation(cap_getEnd(cap)));
    assert(endAlignment1 != NULL);

    Cap *adjacentCap = cap_getAdjacency(cap);
    assert(adjacentCap != NULL);
    assert(cap_getSide(adjacentCap));
    assert(cap_getStrand(adjacentCap));
    adjacentCap = cap_getReverse(adjacentCap);
    stSortedSet *endAlignment2 = stHash_search(endAlignments,
            end_getPositiveOrientation(cap_getEnd(adjacentCap)));
    assert(endAlignment2 != NULL);

    AdjacencySequence *adjacencySequence1 = adjacencySequence_construct(cap,
            INT32_MAX);
    AdjacencySequence *adjacencySequence2 = adjacencySequence_construct(
            adjacentCap, INT32_MAX);
    assert(adjacencySequence1->length == adjacencySequence2->length);
    assert(adjacencySequence1->subsequenceIdentifier == adjacencySequence2->subsequenceIdentifier);
    assert(adjacencySequence1->strand == !adjacencySequence2->strand);
    assert(
            adjacencySequence2->start == adjacencySequence1->start
                    + adjacencySequence1->length - 1);

    stList *inducedAlignment1 = getInducedAlignment(endAlignment1,
            adjacencySequence1);
    stList *inducedAlignment2 = getInducedAlignment(endAlignment2,
            adjacencySequence2);
    stList_reverse(inducedAlignment2);

    fn(cap, inducedAlignment1, inducedAlignment2, endAlignment1, endAlignment2);

    //Cleanup.
    adjacencySequence_destruct(adjacencySequence1);
    adjacencySequence_destruct(adjacencySequence2);
    stList_destruct(inducedAlignment1);
    stList_destruct(inducedAlignment2);
    return 1;
}

int32_t getMaxAdjacencyLength(Flower *flower) {
    assert(flower_getGroupNumber(flower) <= 1);
    Flower_CapIterator *capIt = flower_getCapIterator(flower);
    Cap *cap;
    int32_t maxAdjacencyLength = 0;
    while ((cap = flower_getNextCap(capIt)) != NULL) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        int32_t adjacencyLength = abs(
                cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap)) - 1;
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
    assert(flower_getGroupNumber(flower) <= 1);
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    int32_t maxInstanceNumber = 0;
    End *dominantEnd = NULL;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if(end_getInstanceNumber(end) > maxInstanceNumber) {
            maxInstanceNumber = end_getInstanceNumber(end);
            dominantEnd = end;
        }
    }
    flower_destructEndIterator(endIt);
    if(dominantEnd == NULL) {
        return NULL;
    }
    assert(end_getOrientation(dominantEnd));
    if(end_getInstanceNumber(dominantEnd)*2 < flower_getCapNumber(flower)) {
        return NULL;
    }
    Cap *cap;
    Flower_CapIterator *capIt = flower_getCapIterator(flower);
    while((cap = flower_getNextCap(capIt)) != NULL) {
        assert(cap_getAdjacency(cap) != NULL);
        if(end_getPositiveOrientation(cap_getEnd(cap)) != dominantEnd &&
                end_getPositiveOrientation(cap_getEnd(cap_getAdjacency(cap))) != dominantEnd) {
            flower_destructCapIterator(capIt);
            return NULL;
        }
    }
    flower_destructCapIterator(capIt);
    return dominantEnd;
}

stSortedSet *makeFlowerAlignment(Flower *flower, int32_t spanningTrees,
        int32_t maxSequenceLength, int32_t maximumNumberOfSequencesBeforeSwitchingToFast, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters,
        bool pruneOutStubAlignments) {
    //return stSortedSet_construct3((int(*)(const void *, const void *)) alignedPair_cmpFn,
    //        (void(*)(void *)) alignedPair_destruct);

    //Make the end alignments, representing each as an adjacency alignment.
    End *end;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    stHash *endAlignments = stHash_construct2(NULL,
            (void(*)(void *)) stSortedSet_destruct);
    End *dominantEnd = getDominantEnd(flower); //an end to which all adjacencies are incident
    if (dominantEnd != NULL && getMaxAdjacencyLength(flower) <= 2
            * maxSequenceLength) {
        stHash_insert(
                endAlignments,
                dominantEnd,
                makeEndAlignment(dominantEnd, spanningTrees,
                        2 * maxSequenceLength, maximumNumberOfSequencesBeforeSwitchingToFast, gapGamma,
                        pairwiseAlignmentBandingParameters));
        while ((end = flower_getNextEnd(endIterator)) != NULL) {
            if (end != dominantEnd) {
                assert(stHash_search(endAlignments, end) == NULL);
                stHash_insert(endAlignments, end, stSortedSet_construct());
            }
        }
    } else {
        while ((end = flower_getNextEnd(endIterator)) != NULL) {
            assert(stHash_search(endAlignments, end) == NULL);
            stHash_insert(
                    endAlignments,
                    end,
                    makeEndAlignment(end, spanningTrees, maxSequenceLength,
                            maximumNumberOfSequencesBeforeSwitchingToFast, gapGamma,
                            pairwiseAlignmentBandingParameters));
        }
    }
    flower_destructEndIterator(endIterator);

    //Prune end alignments
    endIterator = flower_getEndIterator(flower);
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
                if (makeFlowerAlignmentP(cap, endAlignments, getScore)) {
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

        makeFlowerAlignmentP(cap, endAlignments, pruneAlignments);
        assert(cap_getAdjacency(cap) != NULL);
        if ((end_isFree(cap_getEnd(cap)) && end_isStubEnd(cap_getEnd(cap)))
                || (end_isFree(cap_getEnd(cap_getAdjacency(cap)))
                        && end_isStubEnd(cap_getEnd(cap_getAdjacency(cap))))) {
            stList_append(freeStubCaps, cap);
        }
    }
    stList_destruct(caps);
    stHash_destruct(capScoresFn_Hash);

    if (pruneOutStubAlignments) {
        pruneStubAlignments_freeStubAdjacencySequences
                = getFreeStubAdjacencySequences(flower);
        while (stList_length(freeStubCaps) > 0) {
            makeFlowerAlignmentP(stList_pop(freeStubCaps), endAlignments,
                    pruneStubAlignments);
        }
        stSortedSet_destruct(pruneStubAlignments_freeStubAdjacencySequences);
    }
    stList_destruct(freeStubCaps);

    //Now convert to set of final aligned pairs to return.
    stSortedSet *sortedAlignment = stSortedSet_construct3(
            (int(*)(const void *, const void *)) alignedPair_cmpFn,
            (void(*)(void *)) alignedPair_destruct);
    stList *endAlignmentsList = stHash_getValues(endAlignments);
    while (stList_length(endAlignmentsList) > 0) {
        stSortedSet *endAlignment = stList_pop(endAlignmentsList);
        while (stSortedSet_size(endAlignment) > 0) {
            AlignedPair *alignedPair = stSortedSet_getFirst(endAlignment);
            stSortedSet_remove(endAlignment, alignedPair);
            assert(
                    stSortedSet_search(endAlignment, alignedPair->reverse)
                            != NULL);
            stSortedSet_remove(endAlignment, alignedPair->reverse);
            assert(stSortedSet_search(sortedAlignment, alignedPair) == NULL);
            assert(
                    stSortedSet_search(sortedAlignment, alignedPair->reverse)
                            == NULL);
            stSortedSet_insert(sortedAlignment, alignedPair);
            stSortedSet_insert(sortedAlignment, alignedPair->reverse);
        }
    }
    stList_destruct(endAlignmentsList);
    stHash_destruct(endAlignments);

    return sortedAlignment;
}

