#include "CuTest.h"
#include "sonLib.h"
#include "multipleAligner.h"
#include "stPosetAlignment.h"
#include <stdlib.h>
#include <string.h>

/*
 * Test the spanning tree generator.
 */

void constructSpanningTree(int32_t numberOfSequences, stSortedSet *pairwiseAlignments);

void test_getSpanningTree(CuTest *testCase) {
    for(int32_t test=0; test<100; test++) {
        int32_t numberOfSequences = st_randomInt(0, 15);
        int32_t spanningTrees = st_randomInt(1, 3);
        stSortedSet *pairwiseAlignments = stSortedSet_construct3((int(*)(
                    const void *, const void *)) stIntTuple_cmpFn,
                    (void(*)(void *)) stIntTuple_destruct);
        for (int i = 0; i < spanningTrees; i++) {
            constructSpanningTree(numberOfSequences, pairwiseAlignments);
        }
        //Check we have all the sequences in at least one pair and that they are joined by a tree.
        stList *list = stList_construct2(numberOfSequences);
        stList *tuples = stList_construct();
        for(int32_t i=0; i<numberOfSequences; i++) {
            stList *list2 = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
            stIntTuple *pos = stIntTuple_construct(1, i);
            stList_append(list2, pos);
            stList_set(list, i, list2);
            stList_append(tuples, pos);
        }

        //Now do merging..
        stSortedSetIterator *it = stSortedSet_getIterator(pairwiseAlignments);
        stIntTuple *pair;
        while((pair = stSortedSet_getNext(it)) != NULL) {
            CuAssertTrue(testCase, stIntTuple_length(pair) == 2);
            int32_t i = stIntTuple_getPosition(pair, 0);
            int32_t j = stIntTuple_getPosition(pair, 1);
            stList *list2 = stList_get(list, i);
            stList *list3 = stList_get(list, j);
            if(list2 != list3) {
                stList_appendAll(list2, list3);
                while(stList_length(list3) > 0) {
                    stIntTuple *pos = stList_pop(list3);
                    stList_set(list, stIntTuple_getPosition(pos, 0), list2);
                }
                stList_destruct(list3);
            }
        }

        //Now check
        for(int32_t i=0; i<numberOfSequences; i++) {
            stList *list2 = stList_get(list, i);
            CuAssertTrue(testCase, stList_length(list2) == numberOfSequences);
            CuAssertTrue(testCase, stList_contains(list2, stList_get(tuples, i)));
        }

        //Cleanup
        if(numberOfSequences > 0) {
            stList_destruct(stList_get(list, 0));
        }
        stList_destruct(list);
        stList_destruct(tuples);
    }
}

/*
 * Test the multiple alignment code with multiple examples.
 */

//Declarations from pairwiseAlignerTest

char *getRandomSequence(int32_t length);

char *evolveSequence(const char *startSequence);

stList *getRandomSequences(int32_t sequenceNumber, int32_t approxLength) {
    stList *sequences = stList_construct3(0, free);
    char *firstSequence = getRandomSequence(approxLength);
    for(int32_t i=0; i<sequenceNumber; i++) {
        stList_append(sequences, evolveSequence(firstSequence));
    }
    return sequences;
}

void test_multipleAlignerRandom(CuTest *testCase) {
    for(int32_t test=0; test<100; test++) {
        stList *randomSequences = getRandomSequences(st_randomInt(0, 10), st_randomInt(0, 100));
        int32_t spanningTrees = st_randomInt(0, 5);
        stPosetAlignment *posetAlignment = stPosetAlignment_construct(stList_length(randomSequences));

        for(int32_t i=0; i<stList_length(randomSequences); i++) {
            st_uglyf("Sequence to align: %s\n", stList_get(randomSequences, i));
        }

        stList *alignedPairs = makeAlignment(randomSequences, spanningTrees, &spanningTrees, 1);
        //Check the aligned pairs.
        stListIterator *iterator = stList_getIterator(alignedPairs);
        stIntTuple *alignedPair;
        while((alignedPair = stList_getNext(iterator)) != NULL) {
            CuAssertTrue(testCase, stIntTuple_length(alignedPair) == 5);
            int32_t score = stIntTuple_getPosition(alignedPair, 0);
            int32_t seqX = stIntTuple_getPosition(alignedPair, 1);
            int32_t x = stIntTuple_getPosition(alignedPair, 2);
            int32_t seqY = stIntTuple_getPosition(alignedPair, 3);
            int32_t y = stIntTuple_getPosition(alignedPair, 4);
            st_uglyf("Got aligned pair, score: %i x seq: %i x pos: %i x seq: %i y pos: %i\n", score, seqX, x, seqY, y);
            CuAssertTrue(testCase, score > 0);
            CuAssertTrue(testCase, score <= 1000);
            CuAssertTrue(testCase, seqX >= 0);
            CuAssertTrue(testCase, seqX < stList_length(randomSequences));
            CuAssertTrue(testCase, x >= 0);
            CuAssertTrue(testCase, x < strlen(stList_get(randomSequences, seqX)));
            CuAssertTrue(testCase, seqY >= 0);
            CuAssertTrue(testCase, seqY < stList_length(randomSequences));
            CuAssertTrue(testCase, y >= 0);
            CuAssertTrue(testCase, y < strlen(stList_get(randomSequences, seqY)));
            //Check we can form an alignment
            CuAssertTrue(testCase, stPosetAlignment_add(posetAlignment, seqX, x, seqY, y));
        }
        stList_destructIterator(iterator);
        stList_destruct(randomSequences);
        stPosetAlignment_destruct(posetAlignment);
    }
}

CuSuite* multipleAlignerTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_multipleAlignerRandom);
    SUITE_ADD_TEST(suite, test_getSpanningTree);
    return suite;
}
