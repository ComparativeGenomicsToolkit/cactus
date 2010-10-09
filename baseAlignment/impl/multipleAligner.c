#include <stdlib.h>

#include "multipleAligner.h"
#include "sonLib.h"
#include "stPosetAlignment.h"
#include "pairwiseAligner.h"

/*
 * Functions to align a bunch of sequences, creating a global alignment.
 */

static int32_t matrix_index(int32_t i, int32_t j, int32_t sequenceNo) {
#ifdef BEN_DEBUG
    assert(i >= 0 && i < sequenceNo);
    assert(j >= 0 && j < sequenceNo);
#endif
    return i * sequenceNo + j;
}

static int32_t matrix_get(int32_t *matrix, int32_t nodeNo, int32_t i, int32_t j) {
    return i < j ? matrix[matrix_index(i, j, nodeNo)] : matrix[matrix_index(j,
            i, nodeNo)];
}

static void matrix_set(int32_t *matrix, int32_t nodeNo, int32_t i, int32_t j,
        int32_t distance) {
    if (i < j) {
        matrix[matrix_index(i, j, nodeNo)] = distance;
    } else {
        matrix[matrix_index(j, i, nodeNo)] = distance;
    }
}

static int32_t *matrix_construct(int32_t nodeNo) {
    int32_t *matrix = st_malloc(sizeof(int32_t) * nodeNo * nodeNo);
    //initialise distance matrix
    for (int32_t i = 0; i < nodeNo; i++) {
        matrix_set(matrix, nodeNo, i, i, 0);
        for (int32_t j = i + 1; j < nodeNo; j++) {
            matrix_set(matrix, nodeNo, i, j, INT32_MAX);
        }
    }
    return matrix;
}

static int32_t addDistances(int32_t i, int32_t j) {
    return i == INT32_MAX || j == INT32_MAX ? INT32_MAX : i + j;
}

static void updatePathDistances(int32_t *distanceMatrix, int32_t *alignedMatrix,
        int32_t sequenceNo, int32_t i, int32_t j, int32_t distance) {
#ifdef BEN_DEBUG
    assert(matrix_get(alignedMatrix, i, j) == INT32_MAX);
#endif
    matrix_set(alignedMatrix, sequenceNo, i, j, 0);
    if (matrix_get(distanceMatrix, sequenceNo, i, j) > distance) {
        matrix_set(distanceMatrix, sequenceNo, i, j, distance);
        for (int32_t k = 0; k < sequenceNo; k++) {
            if (i != k && j != k) {
                int32_t l = addDistances(matrix_get(distanceMatrix, sequenceNo, i, k), distance);
                if (l < matrix_get(distanceMatrix, sequenceNo, j, k)) {
                    matrix_set(distanceMatrix, sequenceNo, j, k, l);
                }
                l = addDistances(matrix_get(distanceMatrix, sequenceNo, j, k), distance);
                if (l < matrix_get(distanceMatrix, sequenceNo, i, k)) {
                    matrix_set(distanceMatrix, sequenceNo, i, k, l);
                }
            }
        }
    }
}

static stIntTuple *getMostDistanceUnalignedPair(int32_t *distanceMatrix,
        int32_t *alignedMatrix, int32_t sequenceNo) {
    /*
     * Gets the pair not-yet aligned with maximum pairwise distance. Breaks ties arbitrarily.
     */
    int32_t sequence1 = INT32_MAX;
    int32_t sequence2 = INT32_MAX;
    int32_t distance = -1;
    for (int32_t i = 0; i < sequenceNo - 1; i++) {
        for (int32_t j = i + 1; j < sequenceNo; j++) {
            if (matrix_get(alignedMatrix, sequenceNo, i, j) == INT32_MAX) { //has not yet been aligned
                int32_t k = matrix_get(distanceMatrix, sequenceNo, i, j);
                if (k > distance) {
                    sequence1 = i;
                    sequence2 = j;
                    distance = k;
                }
            }
        }
    }
    return stIntTuple_construct(3, sequence1, sequence2, distance);
}

/*static int32_t *getIndelProbs(stList *alignedPairs, int32_t sequenceLength, int32_t sequenceIndex) {
    int32_t *indelProbs = st_malloc(sizeof(int32_t)*sequenceLength);
    for(int32_t i=0; i<sequenceLength; i++) {
        indelProbs[i] = 1000;
    }
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int32_t j = stIntTuple_getPosition(alignedPair, sequenceIndex+1);
        int32_t score = stIntTuple_getPosition(alignedPair, 0);
        indelProbs[j] -= score;
        if(indelProbs[j] < 0) {
            indelProbs[j] = 0;
        }
    }
    return indelProbs;
}*/

static int64_t getMaxMatchProbs(stList *alignedPairs, int32_t sequenceLength, int32_t sequenceIndex) {
    int32_t *maxMatchProbs = st_malloc(sizeof(int32_t)*sequenceLength);
    for(int32_t i=0; i<sequenceLength; i++) {
        maxMatchProbs[i] = 0;
    }
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int32_t j = stIntTuple_getPosition(alignedPair, sequenceIndex+1);
        int32_t score = stIntTuple_getPosition(alignedPair, 0);
        assert(score >= 0);
        assert(score <= 1000);
        if(score > maxMatchProbs[j]) {
            maxMatchProbs[j] = score;
        }
    }
    int64_t j = 0;
    for(int32_t i=0; i<sequenceLength; i++) {
        j += maxMatchProbs[i];
    }
    free(maxMatchProbs);
    return j;
}


static int32_t getCertaintyOfAlignmentDistance(stList *alignedPairs, int32_t sequenceLength1, int32_t sequenceLength2) {
    double maxMatchProbs = getMaxMatchProbs(alignedPairs, sequenceLength1, 0) + getMaxMatchProbs(alignedPairs, sequenceLength2, 1);
    double alignmentCertainty = sequenceLength1 + sequenceLength2 == 0 ? 1.0 : maxMatchProbs / (sequenceLength1 + sequenceLength2);
    int32_t i = alignmentCertainty * 1000;
    return 1000 - (i < 0 ? 0 : i > 1000 ? 1000 : i);
}

stList *makeAlignment(stList *sequences, int32_t alignmentsPerSequence,
        void *modelParameters) {
    int32_t sequenceNo = stList_length(sequences);
    int32_t *distanceMatrix = matrix_construct(sequenceNo);
    int32_t *alignedMatrix = matrix_construct(sequenceNo);

    int32_t maximumPossibleNumberOfAlignedPairs = (sequenceNo * sequenceNo
            - sequenceNo) / 2;
    int32_t numberOfAlignedPairs = sequenceNo * alignmentsPerSequence;
    if (numberOfAlignedPairs > maximumPossibleNumberOfAlignedPairs) {
        numberOfAlignedPairs = maximumPossibleNumberOfAlignedPairs;
    }

#ifdef BEN_DEBUG
    int32_t maxPathDistance = INT32_MAX;
#endif
    stList *alignedPairs = stList_construct();
    while (numberOfAlignedPairs-- > 0) {
        stIntTuple *pairwiseAlignment = getMostDistanceUnalignedPair(
                distanceMatrix, alignedMatrix, sequenceNo);
        int32_t sequence1 = stIntTuple_getPosition(pairwiseAlignment, 0);
        int32_t sequence2 = stIntTuple_getPosition(pairwiseAlignment, 1);
#ifdef BEN_DEBUG
        assert(sequence1 != INT32_MAX && sequence2 != INT32_MAX);
        assert(maxPathDistance >= stIntTuple_getPosition(pairwiseAlignment, 2));
        maxPathDistance = stIntTuple_getPosition(pairwiseAlignment, 2);
#endif
        char *string1 = stList_get(sequences, sequence1);
        char *string2 = stList_get(sequences, sequence2);
        int32_t sequence1Length = strlen(string1);
        int32_t sequence2Length = strlen(string2);
        stList *alignedPairs2 = getAlignedPairs(string1, string2, modelParameters);
        //int32_t *indelProbs1 = getIndelProbs(alignedPairs2, sequence1Length, 0);
        //int32_t *indelProbs2 = getIndelProbs(alignedPairs2, sequence2Length, 1);
        int32_t distance = getCertaintyOfAlignmentDistance(alignedPairs2, sequence1Length, sequence2Length);
        //free(indelProbs1);
        //free(indelProbs2);

#ifdef BEN_DEBUG
        assert(distance >= 0);
#endif
        updatePathDistances(distanceMatrix, alignedMatrix, sequenceNo, sequence1,
                sequence2, distance);

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
        if (stPosetAlignment_isPossible(posetAlignment, sequence1, position1,
                sequence2, position2)) {
            stPosetAlignment_add(posetAlignment, sequence1, position1,
                    sequence2, position2);
            //Add a converted version to the accepted aligned pairs.
            stList_append(acceptedAlignedPairs, alignedPair);
        } else {
            stIntTuple_destruct(alignedPair);
        }
    }
    stList_destruct(alignedPairs);
    stPosetAlignment_destruct(posetAlignment);
    free(distanceMatrix);
    free(alignedMatrix);

    //Return the accepted pairs
    return acceptedAlignedPairs;
}
