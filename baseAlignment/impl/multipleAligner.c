#include "multipleAligner.h"
#include "sonLib.h"
#include "stPosetAlignment.h"
#include "pairwiseAligner.h"
#include <stdlib.h>

/*
 * Functions to align a bunch of sequences, creating a global alignment.
 */

/*
 * Constructs a random spanning tree linking all the nodes in items into one component.
 */
void constructSpanningTree(int32_t numberOfSequences,
        stSortedSet *pairwiseAlignments) {
    stList *list = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < numberOfSequences; i++) {
        stList_append(list, stIntTuple_construct(1, i));
    }
    while (stList_length(list) > 1) {
        stIntTuple *i = st_randomChoice(list);
        stList_removeItem(list, i);
        int32_t j = stIntTuple_getPosition(i, 0);
        stIntTuple_destruct(i);
        stIntTuple *k = st_randomChoice(list);
        int32_t l = stIntTuple_getPosition(k, 0);
        assert(l != j);
        stIntTuple *m = j < l ? stIntTuple_construct(2, j, l)
                : stIntTuple_construct(2, l, j);
        if (!stSortedSet_search(pairwiseAlignments, m)) {
            stSortedSet_insert(pairwiseAlignments, m);
        } else {
            stIntTuple_destruct(m);
        }
    }
    stList_destruct(list);
}

static int32_t *calculateIndelProbs(stList *alignedPairs,
        int32_t sequenceLength, int32_t sequenceIndex) {
    int32_t *indelProbs = st_malloc(sizeof(int32_t) * (sequenceLength + 1));
    indelProbs[0] = sequenceLength;
    for (int32_t i = 0; i < sequenceLength; i++) {
        indelProbs[i + 1] = PAIR_ALIGNMENT_PROB_1;
    }
    for (int32_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int32_t j = stIntTuple_getPosition(alignedPair, sequenceIndex + 1);
        assert(j >= 0);
        assert(j < sequenceLength);
        int32_t score = stIntTuple_getPosition(alignedPair, 0);
        assert(score >= 0);
        assert(score <= PAIR_ALIGNMENT_PROB_1);
        indelProbs[j + 1] -= score;
        if (indelProbs[j + 1] < 0) {
            indelProbs[j + 1] = 0;
        }
    }
    return indelProbs;
}

static void calculateIndelProbsAndAddToHash(stList *alignedPairs,
        int32_t sequence1Length, int32_t sequence1, int32_t sequence2,
        stHash *indelProbsHash, int32_t sequenceIndex) {
    stIntTuple *i = stIntTuple_construct(2, sequence1, sequence2);
    assert(stHash_search(indelProbsHash, i) == NULL);
    stHash_insert(indelProbsHash, i, calculateIndelProbs(alignedPairs,
            sequence1Length, sequenceIndex));
}

static int32_t getIndelProb(int32_t sequence1, int32_t position1,
        int32_t sequence2, stHash *indelProbsHash) {
    stIntTuple *i = stIntTuple_construct(2, sequence1, sequence2);
    int32_t *indelProbs = stHash_search(indelProbsHash, i);
    assert(indelProbs != NULL);
    int32_t sequenceLength = indelProbs[0];
    assert(position1 >= 0);
    assert(position1 <= sequenceLength);
    stIntTuple_destruct(i);
    int32_t j = indelProbs[position1 + 1];
    assert(j >= 0);
    assert(j <= PAIR_ALIGNMENT_PROB_1);
    return j;
}

/*static int64_t getMaxMatchProbs(stList *alignedPairs, int32_t sequenceLength,
 int32_t sequenceIndex) {
 int32_t *maxMatchProbs = st_malloc(sizeof(int32_t) * sequenceLength);
 for (int32_t i = 0; i < sequenceLength; i++) {
 maxMatchProbs[i] = 0;
 }
 for (int32_t i = 0; i < stList_length(alignedPairs); i++) {
 stIntTuple *alignedPair = stList_get(alignedPairs, i);
 int32_t j = stIntTuple_getPosition(alignedPair, sequenceIndex + 1);
 int32_t score = stIntTuple_getPosition(alignedPair, 0);
 assert(score >= 0);
 assert(score <= PAIR_ALIGNMENT_PROB_1);
 if (score > maxMatchProbs[j]) {
 maxMatchProbs[j] = score;
 }
 }
 int64_t j = 0;
 for (int32_t i = 0; i < sequenceLength; i++) {
 j += maxMatchProbs[i];
 }
 free(maxMatchProbs);
 return j;
 }*/

stList *makeAlignment(stList *sequences, int32_t spanningTrees,
        void *modelParameters) {
    //Get the set of pairwise alignments (by constructing spanning trees)
    stSortedSet *pairwiseAlignments = stSortedSet_construct3((int(*)(
            const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);

    int32_t sequenceNo = stList_length(sequences);
    int32_t maxAlignmentPairs = (sequenceNo * sequenceNo - sequenceNo) / 2;
    if (spanningTrees * sequenceNo < maxAlignmentPairs) {
        for (int32_t i = 0; i < spanningTrees; i++) {
            constructSpanningTree(sequenceNo, pairwiseAlignments);
        }
    } else {
        for (int32_t i = 0; i < stList_length(sequences); i++) {
            for (int32_t j = i + 1; j < stList_length(sequences); j++) {
                stSortedSet_insert(pairwiseAlignments, stIntTuple_construct(2,
                        i, j));
            }
        }
    }

    //Indel probs hash
    stHash *indelProbsHash = stHash_construct3(
            (uint32_t(*)(const void *)) stIntTuple_hashKey, (int(*)(
                    const void *, const void *)) stIntTuple_equalsFn, (void(*)(
                    void *)) stIntTuple_destruct, free);

    //Construct the alignments
    //and sort them by weight
    stSortedSetIterator *pairwiseAlignmentsIterator = stSortedSet_getIterator(
            pairwiseAlignments);
    stIntTuple *pairwiseAlignment;
    stList *alignedPairs = stList_construct();
    while ((pairwiseAlignment = (stIntTuple *) stSortedSet_getNext(
            pairwiseAlignmentsIterator)) != NULL) {
        int32_t sequence1 = stIntTuple_getPosition(pairwiseAlignment, 0);
        int32_t sequence2 = stIntTuple_getPosition(pairwiseAlignment, 1);
        stList *alignedPairs2 = getAlignedPairs_Fast(
                stList_get(sequences, sequence1), stList_get(sequences,
                        sequence2), modelParameters);

        //Make indel probs
        calculateIndelProbsAndAddToHash(alignedPairs2, strlen(stList_get(
                sequences, sequence1)), sequence1, sequence2, indelProbsHash, 0);
        calculateIndelProbsAndAddToHash(alignedPairs2, strlen(stList_get(
                sequences, sequence2)), sequence2, sequence1, indelProbsHash, 1);

        //Now deal with the match probs..
        while (stList_length(alignedPairs2) > 0) {
            stIntTuple *alignedPair = (stIntTuple *) stList_pop(alignedPairs2);
#ifdef BEN_DEBUG
            assert(stIntTuple_length(alignedPair) == 3);
#endif
            stList_append(alignedPairs, stIntTuple_construct(5,
            /* score */stIntTuple_getPosition(alignedPair, 0),
            /*seq 1 */sequence1, stIntTuple_getPosition(alignedPair, 1),
            /*seq 2 */sequence2, stIntTuple_getPosition(alignedPair, 2)));
            stIntTuple_destruct(alignedPair);
        }
        stList_destruct(alignedPairs2);
    }
    stSortedSet_destructIterator(pairwiseAlignmentsIterator);

    //Sort the pairs (by weight)
    stList_sort(alignedPairs,
            (int(*)(const void *, const void *)) stIntTuple_cmpFn);

    //Greedily construct poset and filter pairs..
    stPosetAlignment *posetAlignment = stPosetAlignment_construct(
            stList_length(sequences));
    stList *acceptedAlignedPairs = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    int32_t pScore = INT32_MAX;
    while (stList_length(alignedPairs) > 0) {
        stIntTuple *alignedPair = stList_pop(alignedPairs);
        int32_t score = stIntTuple_getPosition(alignedPair, 0);
#ifdef BEN_DEBUG
        assert(score > 0);
        assert(score <= pScore);
#endif
        pScore = score;
        int32_t sequence1 = stIntTuple_getPosition(alignedPair, 1);
        int32_t position1 = stIntTuple_getPosition(alignedPair, 2);
        int32_t sequence2 = stIntTuple_getPosition(alignedPair, 3);
        int32_t position2 = stIntTuple_getPosition(alignedPair, 4);

        //score >= 0.1 * indelProb1 && score >= 0.1 * indelProb2 &&
        if (stPosetAlignment_isPossible(posetAlignment, sequence1, position1,
                sequence2, position2)) {
            int32_t indelProb1 = getIndelProb(sequence1, position1, sequence2,
                    indelProbsHash);
            int32_t indelProb2 = getIndelProb(sequence2, position2, sequence1,
                    indelProbsHash);
            int32_t gamma = 0.5;
            if (score >= gamma * indelProb1 && score >= gamma * indelProb2) {
                stPosetAlignment_add(posetAlignment, sequence1, position1,
                        sequence2, position2);
                //Add a converted version to the accepted aligned pairs.
                stList_append(acceptedAlignedPairs, alignedPair);
            } else {
                //stPosetAlignment_add(posetAlignment, sequence1, position1,
                //        sequence2, position2);
                //Add a converted version to the accepted aligned pairs.
                //stList_append(acceptedAlignedPairs, alignedPair);

                stIntTuple_destruct(alignedPair);
            }
        } else {
            stIntTuple_destruct(alignedPair);
        }
    }
    stList_destruct(alignedPairs);
    stPosetAlignment_destruct(posetAlignment);
    stSortedSet_destruct(pairwiseAlignments);
    stHash_destruct(indelProbsHash);

    //Return the accepted pairs
    return acceptedAlignedPairs;
}
