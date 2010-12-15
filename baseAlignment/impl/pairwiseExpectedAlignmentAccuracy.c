/*
 * pairwiseExpectedAlignmentAccuracy.c
 *
 *  Created on: 14 Dec 2010
 *      Author: benedictpaten
 */

static int cmpFn(int32_t i, int32_t j) {
    return i > j ? 1 : (i < j ? -1 : 0);
}

static int getExpectedAlignmentAccuracyScoresP(const void *a, const void *b) {
    int32_t i = cmpFn(stIntTuple_getPosition((stIntTuple *)a, 1), x2 = stIntTuple_getPosition((stIntTuple *)b, 1)); //compare by x and y coordinates
    return i != 0 ? i : cmpFn(stIntTuple_getPosition((stIntTuple *)a, 2), x2 = stIntTuple_getPosition((stIntTuple *)b, 2));
}

static int getExpectedAlignmentAccuracyScoresP2(const void *a, const void *b) {
    return cmpFn(stIntTuple_getPosition((stIntTuple *)a, 2), x2 = stIntTuple_getPosition((stIntTuple *)b, 2)); //compare by y coordinate only
}

static int64_t getCummulativeIndelProbs(int32_t indelProbs, int32_t length) {

}

static int64_t getIndelScore(int64_t indelProbs, int32_t from, int32_t to, int32_t length) {

}

static stHash *getExpectedAlignmentAccuracyScores_Forward(stList *alignedPairs, int32_t *indelProbs1, int32_t *indelProbs2, int32_t seqLength1, int32_t seqLength2) {
    stHash *scores = stHash_construct2(NULL, free);
    stList_sort(alignedPairs, getExpectedAlignmentAccuracyScoresP);
    stSortedSet *pLine = NULL, *line = stSortedSet_construct3(getExpectedAlignmentAccuracyScoresP2, NULL);
    int32_t i = 0, index = 0, x, y, score;
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        score = stIntTuple_getPosition(alignedPair, 0);
        x = stIntTuple_getPosition(alignedPair, 1);
        y = stIntTuple_getPosition(alignedPair, 2);
        if(x > index) {
            index++; //newline
            if(pLine != NULL) {
                stSortedSet_destruct(pLine);
            }
            pLine = line;
            line = stSortedSet_copyConstruct(line, NULL);
        }
        else {
            stIntTuple *alignedPair2 = stSortedSet_searchLessThan(pLine, alignedPair);
            int64_t score = stIntTuple_getPosition(alignedPair, 0) * 2;
            int64_t totalScore =  score + (alignedPair2 == NULL ?
                    getIndelScore(indelProbs1, 0, x) + getIndelScore(indelProbs2, 0, y) :
                    *(int64_t *)stHash_search(scores, alignedPair2) +
                    getIndelScore(indelProbs1, stIntTuple_getPosition(alignedPair2, 1)+1, x) +
                    getIndelScore(indelProbs2, stIntTuple_getPosition(alignedPair2, 2)+1, y));

            //Add score to hash..
            int64_t iA = st_malloc(sizeof(int64_t) * 1);
            iA[0] = score;
            stHash_insert(scores, alignedPair, iA);
            //Now update line..
            if(indelProbs1[x] + indelProbs[y] < score) { //It is more probable that it is a match
                alignedPair2 = stSortedSet_searchLessThanOrEqual(pLine, alignedPair);
                if(alignedPair2 == NULL) { //There is no pair in the set, so we can just add this one
                    stSortedSet_insert(line, alignedPair);
                }
                else { //Check if the previous point eclipses the new point
                    int32_t x2 = stIntTuple_getPosition(alignedPair2, 1);
                    int32_t y2 = stIntTuple_getPosition(alignedPair2, 2);
                    int64_t totalScore2 = *(int64_t *)stHash_search(scores, alignedPair2) +
                            (x2 < x ? getIndelScore(indelProb1, x2+1, x+1) : 0) +
                            (y2 < y ? getIndelScore(indelProb1, y2+1, y+1) : 0);
                    if(totalScore2 < totalScore) { //The previous point does not eclipse the new point, so we can insert this one..
                        stSortedSet_insert(line, alignedPair); //possibly replacing alignedPair2 if y == y2.
                        //Now remove any points eclipsed by this new point
                        stSortedSetIterator *it = stSortedSet_getIteratorFrom(line, alignedPair);
                        stSortedSet_getNext(it);
                        stList *list = stList_construct();
                        while(alignedPair2 = stSortedSet_getNext(it) != NULL) {
                            if(*(int64_t *)stHash_search(scores, alignedPair2) + getIndelScore(indelProb1, stIntTuple_getPosition(alignedPair2, 1)+1, x+1) <
                                    totalScore + getIndelScore(indelProb2, y+1, stIntTuple_getPosition(alignedPair2, 2)+1)) { //is eclipsed
                                stList_append(list, alignedPair2);
                            }
                            else {
                                break;
                            }
                        }
                        stSortedSet_destructIterator(it);
                        while(stList_length(list) > 0) {
                            stSortedSet_remove(list, stList_pop(list));
                        }
                        stList_destruct(list);
                    }
                }
            }
        }
    }
    return scores;
}
