/*
 * chaining.c
 *
 *  Created on: 8 Mar 2012
 *      Author: benedictpaten
 */
#include <stdio.h>
#include <stdlib.h>

#include "sonLib.h"
#include "chaining.h"

APairArray *aPairArray_construct(stList *blastPairs) {
    APairArray *aPairArray = st_malloc(sizeof(APairArray));
    aPairArray->length = stList_length(blastPairs);
    aPairArray->aPairs = st_malloc(sizeof(APair) * aPairArray->length);
    int32_t pX = INT32_MIN, pY = INT32_MIN;
    for (int32_t i = 0; i < aPairArray->length; i++) {
        stIntTuple *blastPair = stList_get(blastPairs, i);
        APair aPair = aPairArray->aPairs[i];
        aPair.x = stIntTuple_getPosition(blastPair, 0);
        aPair.y = stIntTuple_getPosition(blastPair, 1);
        //Check coordinates are sorted
        assert(pX <= aPair.x);
        if (pX == aPair.x) {
            assert(pY <= aPair.y);
        }
        pX = aPair.x;
        pY = aPair.y;
        //
        aPair.fScore = 1;
        aPair.bScore = 1;
    }
    return aPairArray;
}

void aPairArray_destruct(APairArray *aPairArray) {
    free(aPairArray->aPairs);
    free(aPairArray);
}

int aPair_sortByYCoordinate(APair *pair1, APair *pair2) {
    return pair1->y > pair2->y ? 1 : (pair1->y < pair2->y ? -1 : 0);
}

void assignScoreForwards(APair *pair, stSortedSet *predecessors) {
    APair *pair2 = stSortedSet_searchLessThan(predecessors, pair);
    if (pair2 != NULL) { //Nothing earlier, so insert
        pair->fScore += pair2->fScore;
    }
}

void insertForwards(APair *pair, stSortedSet *predecessors) {
    APair *pair2 = stSortedSet_searchLessThanOrEqual(predecessors, pair);
    if (pair2 == NULL || pair2->fScore < pair->fScore) {
        pair2 = stSortedSet_searchGreaterThanOrEqual(predecessors, pair);
        while (pair2->fScore <= pair->fScore) {
            stSortedSet_remove(predecessors, pair2);
        }
        stSortedSet_insert(predecessors, pair);
    }
}

void aPairArray_calculateForwardScores(APairArray *aPairArray) {
    int32_t i = 0;
    stSortedSet *predecessors = stSortedSet_construct3(
            (int(*)(const void *, const void *)) aPair_sortByYCoordinate, NULL);
    while (i < aPairArray->length) {
        APair *pair = &aPairArray->aPairs[i++];
        assignScoreForwards(pair, predecessors);
        int32_t j = i;
        while (j < aPairArray->length) {
            APair *pair2 = &aPairArray->aPairs[j++];
            if (pair2->x != pair->x) {
                break;
            }
            assignScoreForwards(pair2, predecessors);
        }
        //j and i now represent an interval of APairs with the same y coordinate
        insertForwards(pair, predecessors);
        do {
            APair *pair2 = &aPairArray->aPairs[i++];
            insertForwards(pair2, predecessors);
        } while (i < j);
    }
    stSortedSet_destruct(predecessors);
}

void assignScoreBackwards(APair *pair, stSortedSet *successors) {
    APair *pair2 = stSortedSet_searchGreaterThan(successors, pair);
    if (pair2 != NULL) { //Nothing earlier, so insert
        pair->bScore += pair2->bScore;
    }
}

void insertBackwards(APair *pair, stSortedSet *successors) {
    APair *pair2 = stSortedSet_searchGreaterThanOrEqual(successors, pair);
    if (pair2 == NULL || pair2->bScore < pair->bScore) {
        pair2 = stSortedSet_searchLessThanOrEqual(successors, pair);
        while (pair2->bScore <= pair->bScore) {
            stSortedSet_remove(successors, pair2);
        }
        stSortedSet_insert(successors, pair);
    }
}

void aPairArray_calculateBackwardScores(APairArray *aPairArray) {
    int32_t i = aPairArray->length - 1;
    stSortedSet *successors = stSortedSet_construct3(
            (int(*)(const void *, const void *)) aPair_sortByYCoordinate, NULL);
    while (i >= 0) {
        APair *pair = &aPairArray->aPairs[i--];
        assignScoreBackwards(pair, successors);
        int32_t j = i;
        while (j >= 0) {
            APair *pair2 = &aPairArray->aPairs[j--];
            if (pair2->x != pair->x) {
                break;
            }
            assignScoreForwards(pair2, successors);
        }
        //j and i now represent an interval of APairs with the same y coordinate
        insertForwards(pair, successors);
        do {
            APair *pair2 = &aPairArray->aPairs[i--];
            insertForwards(pair2, successors);
        } while (i > j);
    }
    stSortedSet_destruct(successors);
}

double aPairArray_maxScore(APairArray *aPairArray) {
    double maxScore = INT64_MIN;
    for (int32_t i = 0; i < aPairArray->length; i++) {
        APair *aPair = &aPairArray->aPairs[i];
        if (aPair->fScore + aPair->bScore > maxScore) {
            maxScore = aPair->fScore + aPair->bScore;
        }
    }
    return maxScore;
}

stList *filterToRemoveOverlap(stList *overlappingPairs) {
    stList *nonOverlappingPairs = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);

    //Traverse backwards
    stSortedSet *set = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn, NULL);
    int32_t pX = INT32_MAX, pY = INT32_MIN;
    for (int32_t i = stList_length(overlappingPairs) - 1; i >= 0; i--) {
        stIntTuple *pair = stList_get(overlappingPairs, i);
        int32_t x = stIntTuple_getPosition(pair, 0);
        int32_t y = stIntTuple_getPosition(pair, 1);
        if (x < pX && y < pY) {
            stSortedSet_insert(set, pair);
        }
        pX = x;
        pY = y;
    }

    //Traverse forwards to final set of pairs
    pX = INT32_MIN;
    pY = INT32_MIN;
    for (int32_t i = 0; i < stList_length(overlappingPairs); i++) {
        stIntTuple *pair = stList_get(overlappingPairs, i);
        int32_t x = stIntTuple_getPosition(pair, 0);
        int32_t y = stIntTuple_getPosition(pair, 1);
        if (x > pX && y > pY && stSortedSet_search(set, pair) == NULL) {
            stList_append(nonOverlappingPairs, stIntTuple_construct(2, x, y));
        }
        pX = x;
        pY = y;
    }
    stSortedSet_destruct(set);

    return nonOverlappingPairs;
}

stList *aPairArray_selectHighScoringPairs(APairArray *aPairArray, double alpha) {
    //Get max score
    double permissableScore = aPairArray_maxScore(aPairArray) * alpha;
    stList *finalPairs = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < aPairArray->length; i++) {
        APair *aPair = &aPairArray->aPairs[i];
        if (aPair->fScore + aPair->bScore >= permissableScore) {
            stList_append(finalPairs,
                    stIntTuple_construct(2, aPair->x, aPair->y));
        }
    }
    return finalPairs;
}

stList *getAnchorChain(stList *sortedBlastPairs, double alpha) {
    //Initialise the array
    APairArray *aPairArray = aPairArray_construct(sortedBlastPairs);

    //Get forward scores
    aPairArray_calculateForwardScores(aPairArray);

    //Get backward scores
    aPairArray_calculateBackwardScores(aPairArray);

    //Get sub-set of anchors on good chains
    stList *highScoringPairs = aPairArray_selectHighScoringPairs(aPairArray,
            alpha);

    //Now filter so that pairs do not overlap
    stList *chosenPairs = filterToRemoveOverlap(highScoringPairs);
    stList_destruct(highScoringPairs);

    //Cleanup
    aPairArray_destruct(aPairArray);

    return chosenPairs;
}
