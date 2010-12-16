/*
 * pairwiseExpectedAlignmentAccuracy.c
 *
 *  Created on: 14 Dec 2010
 *      Author: benedictpaten
 */

#include "multipleAligner.h"
#include "sonLib.h"
#include "stPosetAlignment.h"
#include "pairwiseAligner.h"
#include <stdlib.h>

static int cmpFn(int32_t i, int32_t j) {
    return i > j ? 1 : (i < j ? -1 : 0);
}

static int getExpectedAlignmentAccuracyScoresP(const void *a, const void *b) {
#ifdef BEN_DEBUG
    assert(stIntTuple_length((stIntTuple *)a) == 3);
    assert(stIntTuple_length((stIntTuple *)b) == 3);
#endif
    int32_t i = cmpFn(stIntTuple_getPosition((stIntTuple *)a, 1), stIntTuple_getPosition((stIntTuple *)b, 1)); //compare by x and y coordinates
    return i != 0 ? i : cmpFn(stIntTuple_getPosition((stIntTuple *)a, 2), stIntTuple_getPosition((stIntTuple *)b, 2));
}

static int getExpectedAlignmentAccuracyScoresP2(const void *a, const void *b) {
#ifdef BEN_DEBUG
    assert(stIntTuple_length((stIntTuple *)a) == 3);
    assert(stIntTuple_length((stIntTuple *)b) == 3);
#endif
    return cmpFn(stIntTuple_getPosition((stIntTuple *)a, 2), stIntTuple_getPosition((stIntTuple *)b, 2)); //compare by y coordinate only
}

static int64_t *getCummulativeIndelProbs(int32_t *indelProbs, int32_t length) {
    int64_t *cummulativeIndelProbs = st_malloc(sizeof(int64_t) * (length+1));
    cummulativeIndelProbs[0] = 0.0;
    for(int32_t i=0; i<length; i++) {
        cummulativeIndelProbs[i+1] = indelProbs[i] + cummulativeIndelProbs[i];
    }
    return cummulativeIndelProbs;
}

static int64_t getIndelScore(int64_t *cummulativeIndelProbs, int32_t from, int32_t to, int32_t length) {
#ifdef BEN_DEBUG
    assert(from <= to);
    assert(from >= 0);
    assert(to <= length);
#endif
    return from == to ? 0.0 : cummulativeIndelProbs[to] - cummulativeIndelProbs[from];
}

static int32_t *getMinYArray(stList *alignedPairs, int32_t seqLength1) {
    int32_t index = seqLength1;
    int32_t *minY = st_malloc(sizeof(int32_t) * (seqLength1+1));
    int32_t j=INT32_MAX;
    for(int32_t i=stList_length(alignedPairs)-1; i>=0; i--) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int32_t x = stIntTuple_getPosition(alignedPair, 1);
        int32_t y = stIntTuple_getPosition(alignedPair, 2);
        while(x < index) {
            minY[index--] = j; //the minimum value of y seen on that line.
        }
        if(y < j) {
            j = y;
        }
    }
    while(index >= 0) {
        minY[index--] = j;
    }
    return minY;
}

static stHash *getExpectedAlignmentAccuracyScores_Forward(stList *alignedPairs, int32_t *indelProbs1, int32_t *indelProbs2, int32_t seqLength1, int32_t seqLength2) {
    //Make the basic datastructures
    stHash *scores = stHash_construct3((uint32_t(*)(const void *))stIntTuple_hashKey, (int(*)(const void *,const void *))stIntTuple_equalsFn, NULL, free);
    alignedPairs = stList_copy(alignedPairs, NULL);
    stList_sort(alignedPairs, getExpectedAlignmentAccuracyScoresP);
    int64_t *cummulativeIndelProbs1 = getCummulativeIndelProbs(indelProbs1, seqLength1);
    int64_t *cummulativeIndelProbs2 = getCummulativeIndelProbs(indelProbs2, seqLength2);
    //Now the stuff for loop
    stSortedSet *pLine = NULL, *line = stSortedSet_construct3(getExpectedAlignmentAccuracyScoresP2, NULL); //Lines sorted by y coordinate only
    int32_t *minY = getMinYArray(alignedPairs, seqLength1);
    int32_t index = -1;
    for(int32_t i=0; i<stList_length(alignedPairs); i++) { //Main loop
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int32_t x = stIntTuple_getPosition(alignedPair, 1);
        int32_t y = stIntTuple_getPosition(alignedPair, 2);

#ifdef BEN_DEBUG
        assert(x < seqLength1);
        assert(y < seqLength2);
        assert(x >= index);
#endif

        if(x > index) {
            index = x; //newline (which may have skipped some previous lines)

            //Prune unnecessary pairs from the line
            stSortedSetIterator *it = stSortedSet_getIterator(line);
            stIntTuple *alignedPair2;
            stList *list = stList_construct();
            while((alignedPair2 = stSortedSet_getNext(it)) != NULL) {
                if(stIntTuple_getPosition(alignedPair2, 2) < minY[x]) {
                    stList_append(list, alignedPair2);
                }
                else {
                    break;
                }
            }
            stSortedSet_destructIterator(it);
            for(int32_t i=0; i<stList_length(list)-1; i++) {
                stSortedSet_remove(line, stList_get(list, i));
            }
            stList_destruct(list);

            if(pLine != NULL) {
                stSortedSet_destruct(pLine);
            }
            pLine = line;
            line = stSortedSet_copyConstruct(line, NULL);
        }

        stIntTuple *alignedPair2 = stSortedSet_searchLessThan(pLine, alignedPair);
#ifdef BEN_DEBUG
        if(alignedPair2 != NULL) {
            assert(stHash_search(scores, alignedPair2) != NULL); //has score
            assert(stIntTuple_getPosition(alignedPair2, 1) < x); //precedes pair in both dimensions
            assert(stIntTuple_getPosition(alignedPair2, 2) < y);
        }
#endif
        //Calculate the score of the alignment
        int64_t score = stIntTuple_getPosition(alignedPair, 0) * 2; //Factor of 2 to make parity with indels.
        int64_t totalScore =  score + (alignedPair2 == NULL ? //First case when there is no preceding aligned pair
                getIndelScore(cummulativeIndelProbs1, 0, x, seqLength1) + getIndelScore(cummulativeIndelProbs2, 0, y, seqLength2) :
                *(int64_t *)stHash_search(scores, alignedPair2) +
                getIndelScore(cummulativeIndelProbs1, stIntTuple_getPosition(alignedPair2, 1)+1, x, seqLength1) +
                getIndelScore(cummulativeIndelProbs2, stIntTuple_getPosition(alignedPair2, 2)+1, y, seqLength2));

        //Add score to hash..
        int64_t *iA = st_malloc(sizeof(int64_t));
        iA[0] = totalScore;
        stHash_insert(scores, alignedPair, iA);
        //Now update line..
        if(indelProbs1[x] + indelProbs2[y] < score) { //It is more probable that it is a match
            alignedPair2 = stSortedSet_searchLessThanOrEqual(line, alignedPair);
            if(alignedPair2 == NULL) { //There is no pair on the line with y coordinate less than or equal, so we can just add this one
                stSortedSet_insert(line, alignedPair);
            }
            else { //Check if the previous point eclipses the new point (has a score which is greater and so makes the new point pointless)
                int32_t x2 = stIntTuple_getPosition(alignedPair2, 1);
                int32_t y2 = stIntTuple_getPosition(alignedPair2, 2);
#ifdef BEN_DEBUG
                assert(x2 <= x); //precedes pair in one or both dimensions
                assert(y2 <= y);
                assert(!(x == x2 && y == y2));
#endif
                int64_t totalScore2 = *(int64_t *)stHash_search(scores, alignedPair2) +
                        getIndelScore(cummulativeIndelProbs1, x2+1, x+1, seqLength1) +
                        getIndelScore(cummulativeIndelProbs2, y2+1, y+1, seqLength2);
                if(totalScore2 < totalScore) { //The previous point does not eclipse the new point, so we can insert this one..
                    stSortedSet_insert(line, alignedPair); //possibly replacing alignedPair2 if y == y2.
                    //Now remove any points eclipsed by this new point
                    stSortedSetIterator *it = stSortedSet_getIteratorFrom(line, alignedPair);
                    alignedPair2 = stSortedSet_getNext(it);
#ifdef BEN_DEBUG
                    assert(alignedPair2 == alignedPair);
#endif
                    stList *list = stList_construct();
                    while((alignedPair2 = stSortedSet_getNext(it)) != NULL) {
                        x2 = stIntTuple_getPosition(alignedPair2, 1);
                        y2 = stIntTuple_getPosition(alignedPair2, 2);
#ifdef BEN_DEBUG
                        assert(x2 <= x);
                        assert(y2 > y);
                        assert(stHash_search(scores, alignedPair2) != NULL);
#endif
                        if(*(int64_t *)stHash_search(scores, alignedPair2) + getIndelScore(cummulativeIndelProbs1, x2+1, x+1, seqLength1) <
                                totalScore + getIndelScore(cummulativeIndelProbs2, y+1, y2+1, seqLength2)) { //is eclipsed
                            stList_append(list, alignedPair2);
                        }
                        else {
                            break;
                        }
                    }
                    stSortedSet_destructIterator(it);
                    //Remove the eclipsed pairs from the line
                    while(stList_length(list) > 0) {
                        stSortedSet_remove(line, stList_pop(list));
                    }
                    stList_destruct(list);
                }
            }
        }
    }
    //Cleanup
    stList_destruct(alignedPairs);
    free(cummulativeIndelProbs1);
    free(cummulativeIndelProbs2);
    free(minY);
    if(pLine != NULL) {
        stSortedSet_destruct(pLine);
    }
    stSortedSet_destruct(line);

    return scores;
}

static int32_t *reverseArray(int32_t *iA, int32_t length) {
    int32_t *iA2 = st_malloc(sizeof(int32_t) * length);
    for(int32_t i=0; i<length/2; i++) {
        iA2[i] = iA[length-1-i];
        iA2[length-1-i] = iA[i];
    }
    if((length % 2) != 0) { //add in the odd number if required.
        iA2[length/2] = iA[length/2];
    }
    return iA2;
}

static stIntTuple *reversePair(stIntTuple *alignedPair, int32_t seqLength1, int32_t seqLength2) {
#ifdef BEN_DEBUG
    assert(stIntTuple_length(alignedPair) == 3);
#endif
    return stIntTuple_construct(3, stIntTuple_getPosition(alignedPair, 0), seqLength1 - 1 - stIntTuple_getPosition(alignedPair, 1), seqLength2 - 1 - stIntTuple_getPosition(alignedPair, 2));
}

stHash *getExpectedAlignmentAccuracyScores(stList *alignedPairs, int32_t *indelProbs1, int32_t *indelProbs2, int32_t seqLength1, int32_t seqLength2) {
#ifdef BEN_DEBUG
    assert(seqLength1 >= 0);
    assert(seqLength2 >= 0);
#endif
    //Get the forward scores
    stHash *forwardScores = getExpectedAlignmentAccuracyScores_Forward(alignedPairs, indelProbs1, indelProbs2, seqLength1, seqLength2);
    int32_t *reverseIndelProbs1 = reverseArray(indelProbs1, seqLength1);
    int32_t *reverseIndelProbs2 = reverseArray(indelProbs2, seqLength2);
    //Make the reverse pairs
    stList *reverseAlignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stList_append(reverseAlignedPairs, reversePair(stList_get(alignedPairs, i), seqLength1, seqLength2));
    }
    //Get the reverse scores.
    stHash *reverseScores = getExpectedAlignmentAccuracyScores_Forward(reverseAlignedPairs, reverseIndelProbs1, reverseIndelProbs2, seqLength1, seqLength2);
    //Combine the forward and reverse scores.
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        stIntTuple *reverseAlignedPair = stList_get(reverseAlignedPairs, i);
        int64_t *forwardScore = (int64_t *)stHash_remove(forwardScores, alignedPair);
#ifdef BEN_DEBUG
        assert(forwardScore != NULL);
        assert(stHash_search(reverseScores, reverseAlignedPair) != NULL);
        assert(seqLength1 > 0);
        assert(seqLength2 > 0);
#endif
        int64_t reverseScore = *(int64_t *)stHash_search(reverseScores, reverseAlignedPair);
        //stHash_remove(forwardScores, )
        double *d = st_malloc(sizeof(int64_t));
        d[0] = ((double)(forwardScore[0] + reverseScore - 2 * stIntTuple_getPosition(alignedPair, 0)) / (seqLength1 + seqLength2)) / PAIR_ALIGNMENT_PROB_1; //We subtract the aligned pair score count so as not to count it twice.
#ifdef BEN_DEBUG
        assert(d[0] >= -0.00001);
        assert(d[0] <= 1.00001);
#endif
        stHash_insert(forwardScores, alignedPair, d);
        free(forwardScore); //cleanup old memory
    }
    //Cleanup
    free(reverseIndelProbs1);
    free(reverseIndelProbs2);
    stHash_destruct(reverseScores);
    stList_destruct(reverseAlignedPairs);
#ifdef BEN_DEBUG
    assert(stHash_size(forwardScores) == stList_length(alignedPairs));
#endif
    return forwardScores;
}
