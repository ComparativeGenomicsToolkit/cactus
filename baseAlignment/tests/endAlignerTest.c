#include "flowersShared.h"
#include "endAligner.h"
#include "adjacencySequences.h"

void test_alignedPair_cmpFn(CuTest *testCase) {
    stList *list = stList_construct3(0, (void (*)(void *))alignedPair_destruct);

    //(int (*)(const void *, const void *))alignedPair_cmpFn,

    Name seq1 = 5;
    Name seq2 = 10;

    AlignedPair *aP1 = alignedPair_construct(seq1, 2, 1,
                        seq1, 7, 1, 90);
    AlignedPair *aP2 = alignedPair_construct(seq1, 2, 1,
            seq2, 4, 0, 10);
    AlignedPair *aP3 = alignedPair_construct(seq1, 2, 1,
                seq2, 4, 1, 75);
    AlignedPair *aP4 = alignedPair_construct(seq1, 3, 1,
                        seq2, 4, 0, 10);
    AlignedPair *aP5 = alignedPair_construct(seq1, 4, 1,
                            seq2, 4, 1, 90);

    AlignedPair *ordering[5] = { aP2, aP1, aP4, aP5, aP3 };
    for(int32_t i=0; i<5; i++) {
        AlignedPair *aP = ordering[i];
        st_uglyf("I got %i %i\n", aP, aP->reverse);
        stList_append(list, aP);
        stList_append(list, aP->reverse);
    }

    stList_sort(list, (int (*)(const void *, const void *))alignedPair_cmpFn);

    AlignedPair *correctOrdering[10] = { aP1, aP2, aP3, aP4, aP5, aP1->reverse, aP2->reverse, aP4->reverse, aP3->reverse, aP5->reverse };
    for(int32_t i=0; i<10; i++) {
        st_uglyf("Checking %i %i\n", correctOrdering[i], stList_get(list, i));
        CuAssertTrue(testCase, correctOrdering[i] == stList_get(list, i));
    }


    stList_destruct(list);
}

int32_t isInAdjacencySequence(AlignedPair *alignedPair, AdjacencySequence *adjacencySequence) {
    if (alignedPair->sequence == adjacencySequence->sequenceName) {
        if (alignedPair->strand == adjacencySequence->strand) {
            if (alignedPair->strand) {
                if (alignedPair->position >= adjacencySequence->start
                        && alignedPair->position < adjacencySequence->start
                                + adjacencySequence->length) {
                    return 1;
                }
            } else {
                if (alignedPair->position <= adjacencySequence->start
                        && alignedPair->position > adjacencySequence->start
                                - adjacencySequence->length) {
                    return 1;
                }
            }

        }
    }
    return 0;
}

/*
 * Checks that the position referred to is in an adjacency coming from the end.
 */
int32_t isInAdjacency(AlignedPair *alignedPair, End *end, int32_t maxLength) {
    Cap *cap;
    End_InstanceIterator *it = end_getInstanceIterator(end);
    while ((cap = end_getNext(it)) != NULL) {
        if (cap_getSide(cap)) {
            cap = cap_getReverse(cap);
        }
        AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap,
                maxLength);
        int32_t i = isInAdjacencySequence(alignedPair, adjacencySequence);
        adjacencySequence_destruct(adjacencySequence);
        if(i) {
            return 1;
        }
    }
    return 0;
}

static void testMakeEndAlignments(CuTest *testCase) {
    setup();
    End *ends[3] = { end1, end2, end3 };
    int32_t maxLength = 4;
    for (int32_t endIndex = 0; endIndex < 3; endIndex++) {
        End *end = ends[endIndex];
        stSortedSet *endAlignment = makeEndAlignment(end, 5, maxLength, 1,
                &endIndex);

        stSortedSetIterator *iterator = stSortedSet_getIterator(endAlignment);
        AlignedPair *alignedPair;
        //Check pairs are part of valid sequences from end

        while ((alignedPair = stSortedSet_getNext(iterator)) != NULL) {
            CuAssertTrue(testCase, alignedPair->score > 0); //Check score is valid.
            CuAssertTrue(testCase, alignedPair->score <= 1000);
            CuAssertTrue(testCase, stSortedSet_search(endAlignment, alignedPair->reverse) != NULL); //Check other end is in.
            //Check coordinates are in sequence..
            CuAssertTrue(testCase, isInAdjacency(alignedPair, end, maxLength));
        }
        stSortedSet_destructIterator(iterator);
        stSortedSet_destruct(endAlignment);
    }
    teardown();
}

CuSuite* endAlignerTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testMakeEndAlignments);
    SUITE_ADD_TEST(suite, test_alignedPair_cmpFn);
    return suite;
}
