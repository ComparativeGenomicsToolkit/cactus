#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseExpectedAlignmentAccuracy.h"
#include "pairwiseAligner.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int32_t fn(int32_t x, int32_t y, int32_t seqLength1) {
    return x * (seqLength1 + 1) + y;
}

static int64_t max(int64_t i, int64_t j) {
    return i > j ? i : j;
}

stHash *getExpectedAlignmentAccuracyScores2(stList *alignedPairs, int32_t *indelProbs1,
        int32_t *indelProbs2, int32_t seqLength1, int32_t seqLength2) {
    int64_t *fMatrix = st_calloc((seqLength1 + 1) * (seqLength2 + 1), sizeof(int64_t));
    int64_t *bMatrix = st_calloc((seqLength1 + 1) * (seqLength2 + 1), sizeof(int64_t));
    int64_t *fMatrix2 = st_calloc((seqLength1 + 1) * (seqLength2 + 1), sizeof(int64_t));
    int64_t *bMatrix2 = st_calloc((seqLength1 + 1) * (seqLength2 + 1), sizeof(int64_t));
    //Initialise matrices
    fMatrix[fn(0, 0, seqLength1)] = 0;
    bMatrix[fn(seqLength1, seqLength2, seqLength1)] = 0;
    for(int32_t x=1; x<seqLength1+1; x++) {
        fMatrix[fn(x, 0, seqLength1)] = indelProbs1[x-1] + fMatrix[fn(x-1, 0, seqLength1)];
    }
    for(int32_t y=1; y<seqLength2+1; y++) {
        fMatrix[fn(0, y, seqLength1)] = indelProbs2[y-1] + fMatrix[fn(0, y-1, seqLength1)];
    }
    for(int32_t x=seqLength1-1; x>=0; x--) {
        bMatrix[fn(x, 0, seqLength1)] = indelProbs1[x] + bMatrix[fn(x+1, 0, seqLength1)];
    }
    for(int32_t y=seqLength2-1; y>=0; y--) {
        bMatrix[fn(0, y, seqLength1)] = indelProbs1[y] + bMatrix[fn(0, y+1, seqLength1)];
    }
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        fMatrix[fn(stIntTuple_getPosition(alignedPair, 1)+1,
                   stIntTuple_getPosition(alignedPair, 2)+1, seqLength1)] = stIntTuple_getPosition(alignedPair, 0) * 2;
        bMatrix[fn(stIntTuple_getPosition(alignedPair, 1),
                           stIntTuple_getPosition(alignedPair, 2), seqLength1)] = stIntTuple_getPosition(alignedPair, 0) * 2;
    }
    //Forwards
    for(int32_t x=1; x<seqLength1+1; x++) {
        for(int32_t y=1; y<seqLength2+1; y++) {
            int64_t i = fMatrix[fn(x-1, y-1, seqLength1)] + fMatrix[fn(x, y, seqLength1)];
            int64_t j = fMatrix[fn(x-1, y, seqLength1)] + indelProbs1[x-1];
            int64_t k = fMatrix[fn(x, y-1, seqLength1)] + indelProbs2[y-1];
            fMatrix[fn(x, y, seqLength1)] = max(max(i, j), k);
            fMatrix2[fn(x, y, seqLength1)] = i;
        }
    }
    //Backwards
    for(int32_t x=seqLength1; x>=0; x--) {
        for(int32_t y=seqLength2; y>=0; y--) {
            int64_t i = bMatrix[fn(x+1, y+1, seqLength1)] + bMatrix[fn(x, y, seqLength1)];
            int64_t j = bMatrix[fn(x+1, y, seqLength1)] + indelProbs1[x];
            int64_t k = bMatrix[fn(x, y+1, seqLength1)] + indelProbs2[y];
            bMatrix[fn(x, y, seqLength1)] = max(max(i, j), k);
            bMatrix2[fn(x, y, seqLength1)] = i;
        }
    }
    //Now do the summations..
    stHash *scores = stHash_construct2(NULL, free);
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int64_t j = fMatrix2[fn(stIntTuple_getPosition(alignedPair, 1)+1,
                   stIntTuple_getPosition(alignedPair, 2)+1, seqLength1)] +
                    bMatrix2[fn(stIntTuple_getPosition(alignedPair, 1),
                   stIntTuple_getPosition(alignedPair, 2), seqLength1)] -
                   stIntTuple_getPosition(alignedPair, 0);
        double *d = st_malloc(sizeof(double));
        d[0] = (double)j / (seqLength1 + seqLength2 * PAIR_ALIGNMENT_PROB_1);
        stHash_insert(scores, alignedPair, d);
    }
    return scores;
}

int32_t *calculateIndelProbs(stList *alignedPairs,
        int32_t sequenceLength, int32_t sequenceIndex);

char *getRandomSequence(int32_t length);

char *evolveSequence(const char *startSequence);

void test_getExpectedAlignmentAccuracyScores(CuTest *testCase) {
    for(int32_t test=0; test<100; test++) {
        //Generate the inputs
        int32_t seqLength1 = st_randomInt(0, 100);
        char *seq1 = getRandomSequence(seqLength1);
        char *seq2 = evolveSequence(seq1);
        int32_t seqLength2 = strlen(seq2);
        stList *alignedPairs = getAlignedPairs(seq1, seq2);
        int32_t *indelProbs1 = calculateIndelProbs(alignedPairs, seqLength1, 0);
        int32_t *indelProbs2 = calculateIndelProbs(alignedPairs, seqLength2, 1);
        //Run the algorithms
        stHash *scores1 = getExpectedAlignmentAccuracyScores(alignedPairs, indelProbs1, indelProbs2, seqLength1, seqLength2);
        stHash *scores2 = getExpectedAlignmentAccuracyScores2(alignedPairs, indelProbs1, indelProbs2, seqLength1, seqLength2);

        CuAssertTrue(testCase, stHash_size(scores1) == stList_length(alignedPairs));
        CuAssertTrue(testCase, stHash_size(scores2) == stList_length(alignedPairs));
        for(int32_t i=0; i<stList_length(alignedPairs); i++) {
            stIntTuple *alignedPair = stList_get(alignedPairs, i);
            CuAssertTrue(testCase, stHash_search(scores1, alignedPair) != NULL);
            CuAssertTrue(testCase, stHash_search(scores2, alignedPair) != NULL);
            CuAssertDblEquals(testCase, *(double *)stHash_search(scores2, alignedPair), *(double *)stHash_search(scores1, alignedPair), 0.000001);
        }

        //Cleanup
        free(indelProbs1);
        free(indelProbs2);
        stList_destruct(alignedPairs);
        stHash_destruct(scores1);
    }
}

CuSuite* pairwiseExpectedAlignmentAccuracyTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_getExpectedAlignmentAccuracyScores);
    return suite;
}
