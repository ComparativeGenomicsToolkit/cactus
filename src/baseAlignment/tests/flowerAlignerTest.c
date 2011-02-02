#include "flowersShared.h"
#include "flowerAligner.h"
#include "endAligner.h"
#include "adjacencySequences.h"
#include "pairwiseAligner.h"

stList *getInducedAlignment(stSortedSet *endAlignment, AdjacencySequence *adjacencySequence);

static int getRandomPosition(AdjacencySequence *adjacencySequence) {
    if(adjacencySequence->strand) {
        return st_randomInt(adjacencySequence->start, adjacencySequence->start + adjacencySequence->length);
    }
    else {
        return st_randomInt(adjacencySequence->start - adjacencySequence->length + 1, adjacencySequence->start + 1);
    }
}

int32_t isInAdjacencySequence(AlignedPair *alignedPair, AdjacencySequence *adjacencySequence);

stList *getinducedAlignment2(stSortedSet *endAlignment, AdjacencySequence *adjacencySequence) {
    stList *inducedAlignment = stList_construct();
    stSortedSetIterator *it = stSortedSet_getIterator(endAlignment);
    AlignedPair *alignedPair;
    while((alignedPair = stSortedSet_getNext(it)) != NULL) {
        if(isInAdjacencySequence(alignedPair, adjacencySequence)) {
            stList_append(inducedAlignment, alignedPair);
        }
    }
    stSortedSet_destructIterator(it);
    stList_sort(inducedAlignment, (int (*)(const void *, const void *))alignedPair_cmpFn);
    if(!adjacencySequence->strand) {
        stList_reverse(inducedAlignment);
    }
    return inducedAlignment;
}

void test_getInducedAlignment(CuTest *testCase) {
    for(int32_t test=0; test<100; test++) {
        setup();

        stSortedSet *sortedAlignment = stSortedSet_construct3((int (*)(const void *, const void *))alignedPair_cmpFn,
                       (void (*)(void *))alignedPair_destruct);


        stList *adjacencySequences = stList_construct3(0, (void (*)(void *))adjacencySequence_destruct);
        Cap *caps[] = { cap1, cap_getReverse(cap4),
                cap5, cap_getReverse(cap8),
                cap9 };
        for(int32_t i=0; i<5; i++) {
            stList_append(adjacencySequences, adjacencySequence_construct(caps[i], INT32_MAX));
        }

        //Make random aligned pairs
        while(st_random() > 0.001) {
            AdjacencySequence *aS1 = st_randomChoice(adjacencySequences);
            AdjacencySequence *aS2 = st_randomChoice(adjacencySequences);
            if(aS1 != aS2) {
                AlignedPair *alignedPair =
                        alignedPair_construct(aS1->sequenceName, getRandomPosition(aS1), aS1->strand,
                                              aS2->sequenceName, getRandomPosition(aS2), aS2->strand,
                                              st_randomInt(0, PAIR_ALIGNMENT_PROB_1));
                stSortedSet_insert(sortedAlignment, alignedPair);
                stSortedSet_insert(sortedAlignment, alignedPair->reverse);
            }
        }

        for(int32_t i=0; i<stList_length(adjacencySequences); i++) {
            AdjacencySequence *adjacencySequence = stList_get(adjacencySequences, i);
            stList *inducedAlignment = getInducedAlignment(sortedAlignment, adjacencySequence);
            stList *inducedAlignment2 = getinducedAlignment2(sortedAlignment, adjacencySequence);

            /*st_uglyf("The lengths are %i %i\n", stList_length(inducedAlignment), stList_length(inducedAlignment2));
            st_uglyf("Adj %i %i %i %i\n", adjacencySequence->sequenceName, adjacencySequence->start, adjacencySequence->length, adjacencySequence->strand);

            for(int32_t j=0; j<stList_length(inducedAlignment); j++) {
                AlignedPair *aP = stList_get(inducedAlignment, j);
                st_uglyf("Hello %i %i %i\n", aP->sequence, aP->position, aP->strand);
            }

            for(int32_t j=0; j<stList_length(inducedAlignment2); j++) {
                AlignedPair *aP = stList_get(inducedAlignment2, j);
                st_uglyf("Goodbye %i %i %i\n", aP->sequence, aP->position, aP->strand);
            }*/

            CuAssertTrue(testCase, stList_length(inducedAlignment) == stList_length(inducedAlignment2));
            for(int32_t j=0; j<stList_length(inducedAlignment); j++) {
                CuAssertTrue(testCase, stList_get(inducedAlignment, j) == stList_get(inducedAlignment2, j));
            }

            stList_destruct(inducedAlignment);
            stList_destruct(inducedAlignment2);
        }

        //cleanup
        stSortedSet_destruct(sortedAlignment);
        teardown();
    }
}

/*
 * Just runs the flower alignment through, doesn't really check its okay.
 */
void test_flowerAlignerRandom(CuTest *testCase) {
    setup();
    int32_t maxLength = 5;
    stSortedSet *flowerAlignment = makeFlowerAlignment(flower, 5, maxLength,  0.5, 1, pairwiseParameters);
    //Check the aligned pairs are all good..
    stSortedSetIterator *iterator = stSortedSet_getIterator(flowerAlignment);
    AlignedPair *alignedPair;
    while((alignedPair = stSortedSet_getNext(iterator)) != NULL) {
        CuAssertTrue(testCase, alignedPair->score > 0); //Check score is valid
        CuAssertTrue(testCase, alignedPair->score <= PAIR_ALIGNMENT_PROB_1);
        CuAssertTrue(testCase, stSortedSet_search(flowerAlignment, alignedPair->reverse) != NULL); //Check other end is in.
    }
    stSortedSet_destructIterator(iterator);
    stSortedSet_destruct(flowerAlignment);

    teardown();
}

CuSuite* flowerAlignerTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_getInducedAlignment);
    SUITE_ADD_TEST(suite, test_flowerAlignerRandom);
    return suite;
}
