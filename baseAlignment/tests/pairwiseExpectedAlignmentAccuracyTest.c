#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseExpectedAlignmentAccuracy.h"
#include "pairwiseAligner.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int32_t fn(int32_t x, int32_t y, int32_t seqLength1, int32_t seqLength2) {
    assert(x >= 0);
    assert(y >= 0);
    assert(x < seqLength1+1);
    assert(y < seqLength2+1);
    return x * (seqLength2 + 1) + y;
}

static int64_t max(int64_t i, int64_t j) {
    return i > j ? i : j;
}

int64_t *getMatrix(int32_t seqLength1, int32_t seqLength2) {
    int64_t *matrix = st_malloc((seqLength1 + 1) * (seqLength2 + 1) * sizeof(int64_t));
    for(int32_t i=0; i<(seqLength1+1)*(seqLength2+1); i++) {
        matrix[i] = 0;
    }
    return matrix;
}

stHash *getExpectedAlignmentAccuracyScores2(stList *alignedPairs, int32_t *indelProbs1,
        int32_t *indelProbs2, int32_t seqLength1, int32_t seqLength2) {
    int64_t *fMatrix = getMatrix(seqLength1, seqLength2);
    int64_t *bMatrix = getMatrix(seqLength1, seqLength2);
    int64_t *fMatrix2 = getMatrix(seqLength1, seqLength2);
    int64_t *bMatrix2 = getMatrix(seqLength1, seqLength2);
    //Initialise matrices
    fMatrix[fn(0, 0, seqLength1, seqLength2)] = 0;
    bMatrix[fn(seqLength1, seqLength2, seqLength1, seqLength2)] = 0;
    for(int32_t x=1; x<seqLength1+1; x++) {
        assert(fMatrix[fn(x, 0, seqLength1, seqLength2)] == 0);
        assert(indelProbs1[x-1] >= 0);
        fMatrix[fn(x, 0, seqLength1, seqLength2)] = indelProbs1[x-1] + fMatrix[fn(x-1, 0, seqLength1, seqLength2)];
    }
    for(int32_t y=1; y<seqLength2+1; y++) {
        assert(fMatrix[fn(0, y, seqLength1, seqLength2)] == 0);
        assert(indelProbs2[y-1] >= 0);
        fMatrix[fn(0, y, seqLength1, seqLength2)] = indelProbs2[y-1] + fMatrix[fn(0, y-1, seqLength1, seqLength2)];
    }
    for(int32_t x=seqLength1-1; x>=0; x--) {
        assert(bMatrix[fn(x, seqLength2, seqLength1, seqLength2)] == 0);
        bMatrix[fn(x, seqLength2, seqLength1, seqLength2)] = indelProbs1[x] + bMatrix[fn(x+1, seqLength2, seqLength1, seqLength2)];
    }
    for(int32_t y=seqLength2-1; y>=0; y--) {
        assert(bMatrix2[fn(seqLength1, y, seqLength1, seqLength2)] == 0);
        bMatrix[fn(seqLength1, y, seqLength1, seqLength2)] = indelProbs2[y] + bMatrix[fn(seqLength1, y+1, seqLength1, seqLength2)];
    }
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);

        assert(alignedPair != NULL);
        assert(stIntTuple_length(alignedPair) == 3);
        assert(fMatrix2[fn(stIntTuple_getPosition(alignedPair, 1)+1,
                           stIntTuple_getPosition(alignedPair, 2)+1, seqLength1, seqLength2)] == 0);
        assert(fMatrix[fn(stIntTuple_getPosition(alignedPair, 1)+1,
                   stIntTuple_getPosition(alignedPair, 2)+1, seqLength1, seqLength2)] == 0);

        fMatrix[fn(stIntTuple_getPosition(alignedPair, 1)+1,
                   stIntTuple_getPosition(alignedPair, 2)+1, seqLength1, seqLength2)] = (int64_t)stIntTuple_getPosition(alignedPair, 0) * 2;
        bMatrix[fn(stIntTuple_getPosition(alignedPair, 1),
                           stIntTuple_getPosition(alignedPair, 2), seqLength1, seqLength2)] = (int64_t)stIntTuple_getPosition(alignedPair, 0) * 2;
    }
    //Forwards
    for(int32_t x=1; x<seqLength1+1; x++) {
        for(int32_t y=1; y<seqLength2+1; y++) {
            assert(fMatrix[fn(x, y, seqLength1, seqLength2)] >= 0);
            assert(fMatrix[fn(x, y, seqLength1, seqLength2)] <= (int64_t)PAIR_ALIGNMENT_PROB_1 * 2);
            assert(fMatrix[fn(x-1, y-1, seqLength1, seqLength2)] >= 0);
            assert(fMatrix[fn(x-1, y, seqLength1, seqLength2)] >= 0);
            assert(fMatrix[fn(x, y-1, seqLength1, seqLength2)] >= 0);

            int64_t i = fMatrix[fn(x-1, y-1, seqLength1, seqLength2)] + fMatrix[fn(x, y, seqLength1, seqLength2)];
            int64_t j = fMatrix[fn(x-1, y, seqLength1, seqLength2)] + indelProbs1[x-1];
            int64_t k = fMatrix[fn(x, y-1, seqLength1, seqLength2)] + indelProbs2[y-1];
            int64_t l = fMatrix[fn(x-1, y-1, seqLength1, seqLength2)] + indelProbs1[x-1] + indelProbs2[y-1];

            assert(l <= j && l <= k);
            assert(i >= 0);
            assert(j >= 0);
            assert(k >= 0);

            fMatrix[fn(x, y, seqLength1, seqLength2)] = max(max(i, j), max(k, l));
            fMatrix2[fn(x, y, seqLength1, seqLength2)] = i;
            assert(fMatrix[fn(x, y, seqLength1, seqLength2)] >= 0);
            assert(fMatrix2[fn(x, y, seqLength1, seqLength2)] >= 0);
        }
    }
    //Backwards
    for(int32_t x=seqLength1-1; x>=0; x--) {
        for(int32_t y=seqLength2-1; y>=0; y--) {
            assert(bMatrix[fn(x, y, seqLength1, seqLength2)] >= 0);

            int64_t i = bMatrix[fn(x+1, y+1, seqLength1, seqLength2)] + bMatrix[fn(x, y, seqLength1, seqLength2)];
            int64_t j = bMatrix[fn(x+1, y, seqLength1, seqLength2)] + indelProbs1[x];
            int64_t k = bMatrix[fn(x, y+1, seqLength1, seqLength2)] + indelProbs2[y];
            int64_t l = bMatrix[fn(x+1, y+1, seqLength1, seqLength2)] + indelProbs1[x] + indelProbs2[y];

            assert(l <= j && l <= k);

            bMatrix[fn(x, y, seqLength1, seqLength2)] = max(max(i, j), k);
            bMatrix2[fn(x, y, seqLength1, seqLength2)] = i;
            assert(bMatrix[fn(x, y, seqLength1, seqLength2)] >= 0);
            assert(bMatrix2[fn(x, y, seqLength1, seqLength2)] >= 0);
        }
    }
    //Now do the summations..
    stHash *scores = stHash_construct2(NULL, free);
    for(int32_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        assert(alignedPair != NULL);
        assert(stIntTuple_length(alignedPair) == 3);
        int64_t j = fMatrix2[fn(stIntTuple_getPosition(alignedPair, 1)+1,
                   stIntTuple_getPosition(alignedPair, 2)+1, seqLength1, seqLength2)];
        int64_t k =
                    bMatrix2[fn(stIntTuple_getPosition(alignedPair, 1),
                   stIntTuple_getPosition(alignedPair, 2), seqLength1, seqLength2)];
        int64_t l = 2 * stIntTuple_getPosition(alignedPair, 0);
        double *d = st_malloc(sizeof(double));
        d[0] = ((double)(j + k - l) / (seqLength1 + seqLength2)) / PAIR_ALIGNMENT_PROB_1;
#ifdef BEN_DEBUG
        assert(d[0] >= -0.00001);
        assert(d[0] <= 1.00001);
#endif
        stHash_insert(scores, alignedPair, d);
    }
    //Cleanup
    free(fMatrix);
    free(bMatrix);
    free(fMatrix2);
    free(bMatrix2);
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
        st_uglyf("I have the sequences to get maximum expected accuracy values on :%s: :%s: %i %i\n", seq1, seq2, seqLength1, seqLength2);
        stList *alignedPairs = getAlignedPairs(seq1, seq2);
        st_uglyf("I have %i aligned sequence pairs\n", stList_length(alignedPairs));
        int32_t i=stList_length(alignedPairs);
        int32_t *indelProbs1 = calculateIndelProbs(alignedPairs, seqLength1, 0);
        int32_t *indelProbs2 = calculateIndelProbs(alignedPairs, seqLength2, 1);
        //Run the algorithms
        stHash *scores2 = getExpectedAlignmentAccuracyScores2(alignedPairs, indelProbs1, indelProbs2, seqLength1, seqLength2);
        stHash *scores1 = getExpectedAlignmentAccuracyScores(alignedPairs, indelProbs1, indelProbs2, seqLength1, seqLength2);


        CuAssertTrue(testCase, i == stList_length(alignedPairs));
        CuAssertTrue(testCase, stHash_size(scores1) == stList_length(alignedPairs));
        CuAssertTrue(testCase, stHash_size(scores2) == stList_length(alignedPairs));
        for(i=0; i<stList_length(alignedPairs); i++) {
            stIntTuple *alignedPair = stList_get(alignedPairs, i);
            st_uglyf("I have the aligned pair %i %i %i %f %f\n", stIntTuple_getPosition(alignedPair, 0), stIntTuple_getPosition(alignedPair, 1), stIntTuple_getPosition(alignedPair, 2), *(double *)stHash_search(scores2, alignedPair), *(double *)stHash_search(scores1, alignedPair));
            CuAssertTrue(testCase, stHash_search(scores1, alignedPair) != NULL);
            CuAssertTrue(testCase, stHash_search(scores2, alignedPair) != NULL);
            CuAssertDblEquals(testCase, *(double *)stHash_search(scores2, alignedPair), *(double *)stHash_search(scores1, alignedPair), 0.000000001);
        }

        //Cleanup
        free(indelProbs1);
        free(indelProbs2);
        stList_destruct(alignedPairs);
        stHash_destruct(scores1);
        stHash_destruct(scores2);
    }
}

CuSuite* pairwiseExpectedAlignmentAccuracyTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_getExpectedAlignmentAccuracyScores);
    return suite;
}
