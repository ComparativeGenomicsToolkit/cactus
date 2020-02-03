/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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

int64_t isInAdjacencySequence(AlignedPair *alignedPair, AdjacencySequence *adjacencySequence);

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
    for(int64_t test=0; test<100; test++) {
        setup(testCase);

        stSortedSet *sortedAlignment = stSortedSet_construct3((int (*)(const void *, const void *))alignedPair_cmpFn,
                       (void (*)(void *))alignedPair_destruct);


        stList *adjacencySequences = stList_construct3(0, (void (*)(void *))adjacencySequence_destruct);
        Cap *caps[] = { cap1, cap_getReverse(cap4),
                cap5, cap_getReverse(cap8),
                cap9 };
        for(int64_t i=0; i<5; i++) {
            stList_append(adjacencySequences, adjacencySequence_construct(caps[i], INT64_MAX));
        }

        //Make random aligned pairs
        while(st_random() > 0.001) {
            AdjacencySequence *aS1 = st_randomChoice(adjacencySequences);
            AdjacencySequence *aS2 = st_randomChoice(adjacencySequences);
            if(aS1 != aS2) {
                AlignedPair *alignedPair =
                        alignedPair_construct(aS1->subsequenceIdentifier, getRandomPosition(aS1), aS1->strand,
                                              aS2->subsequenceIdentifier, getRandomPosition(aS2), aS2->strand,
                                              st_randomInt(0, PAIR_ALIGNMENT_PROB_1), st_randomInt(0, PAIR_ALIGNMENT_PROB_1));
                stSortedSet_insert(sortedAlignment, alignedPair);
                stSortedSet_insert(sortedAlignment, alignedPair->reverse);
            }
        }

        for(int64_t i=0; i<stList_length(adjacencySequences); i++) {
            AdjacencySequence *adjacencySequence = stList_get(adjacencySequences, i);
            stList *inducedAlignment = getInducedAlignment(sortedAlignment, adjacencySequence);
            stList *inducedAlignment2 = getinducedAlignment2(sortedAlignment, adjacencySequence);

            /*st_logInfo("The lengths are %" PRIi64 " %" PRIi64 "\n", stList_length(inducedAlignment), stList_length(inducedAlignment2));
            st_logInfo("Adj %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 "\n", adjacencySequence->sequenceName, adjacencySequence->start, adjacencySequence->length, adjacencySequence->strand);

            for(int64_t j=0; j<stList_length(inducedAlignment); j++) {
                AlignedPair *aP = stList_get(inducedAlignment, j);
                st_logInfo("Hello %" PRIi64 " %" PRIi64 " %" PRIi64 "\n", aP->sequence, aP->position, aP->strand);
            }

            for(int64_t j=0; j<stList_length(inducedAlignment2); j++) {
                AlignedPair *aP = stList_get(inducedAlignment2, j);
                st_logInfo("Goodbye %" PRIi64 " %" PRIi64 " %" PRIi64 "\n", aP->sequence, aP->position, aP->strand);
            }*/

            CuAssertTrue(testCase, stList_length(inducedAlignment) == stList_length(inducedAlignment2));
            for(int64_t j=0; j<stList_length(inducedAlignment); j++) {
                CuAssertTrue(testCase, stList_get(inducedAlignment, j) == stList_get(inducedAlignment2, j));
            }

            stList_destruct(inducedAlignment);
            stList_destruct(inducedAlignment2);
        }

        //cleanup
        stSortedSet_destruct(sortedAlignment);
        teardown(testCase);
    }
}

/*
 * Just runs the flower alignment through, doesn't really check its okay.
 */
void test_flowerAlignerRandom(CuTest *testCase) {
    setup(testCase);
    int64_t maxLength = 5;
    StateMachine *sM = stateMachine5_construct(fiveState);
    stSortedSet *flowerAlignment = makeFlowerAlignment(sM, flower, 5, maxLength, 1, 0.5, pairwiseParameters, st_random() > 0.5);
    stateMachine_destruct(sM);
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

    teardown(testCase);
}

CuSuite* flowerAlignerTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_getInducedAlignment);
    SUITE_ADD_TEST(suite, test_flowerAlignerRandom);
    return suite;
}
