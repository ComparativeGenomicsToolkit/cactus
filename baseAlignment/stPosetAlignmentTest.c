/*
 * stPosetAlignmentTest.c
 *
 *  Created on: 1 June 2010
 *      Author: benedictpaten
 */

#include <stdlib.h>

#include "sonLib.h"
#include "stPosetAlignment.h"
#include "CuTest.h"

/*
 * The following test code builds a random alignment and checking as each
 * extra pair of aligned positions is added to the alignment that the alignment
 * remains partially ordered.
 */
static const int32_t MAX_SEQUENCE_NUMBER = 10;
static const int32_t MAX_SEQUENCE_SIZE = 100;
static const int32_t MAX_ALIGNMENTS = 100;
static const int32_t MAX_ALIGNED_PAIRS = 10000;


static int32_t sequenceNumber;
static stPosetAlignment *posetAlignment = NULL;

static void teardown() {
    if(posetAlignment != NULL) {
        stPosetAlignment_destruct(posetAlignment);
    }
}

static void setup() {
    teardown();
    sequenceNumber = st_randomInt(0, MAX_SEQUENCE_NUMBER);
    posetAlignment = stPosetAlignment_construct(sequenceNumber);
}


/*
 * Tests the constructor and the basic setup.
 */
static void test_stPosetAlignment_getSequenceNumber(CuTest *testCase) {
    setup();
    CuAssertTrue(testCase, sequenceNumber == stPosetAlignment_getSequenceNumber(posetAlignment));
    teardown();
}


/*
 * This builds an adjacency list structure for the the sequences. Every sequence-position
 * has a column in the hash with which it can be aligned with.
 */
static stHash *buildAdjacencyList(stList *pairs, int32_t sequenceNumber) {
    stHash *hash = stHash_construct3((uint32_t (*)(const void *))stIntTuple_hashKey,
            (int (*)(const void *, const void *))stIntTuple_equalsFn,
            (void (*)(void *))stIntTuple_destruct, (void (*)(void *))stSortedSet_destruct);
    for(int32_t seq=0; seq<sequenceNumber; seq++) {
        for(int32_t position=0; position<MAX_SEQUENCE_SIZE; position++) {
            stIntTuple *seqPos = stIntTuple_construct(2, seq, position);
            stSortedSet *column = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, NULL);
            stSortedSet_insert(column, seqPos);
            stHash_insert(hash, seqPos, column);
        }
    }
    stListIterator *it = stList_getIterator(pairs);
    stIntTuple *pair;
    while((pair = stList_getNext(it)) != NULL) {
       stIntTuple *seqPos1 = stIntTuple_construct(2, stIntTuple_getPosition(pair, 0), stIntTuple_getPosition(pair, 1));
       stIntTuple *seqPos2 = stIntTuple_construct(2, stIntTuple_getPosition(pair, 2), stIntTuple_getPosition(pair, 3));
       stSortedSet *column1 = stHash_search(hash, seqPos1);
       stSortedSet *column2 = stHash_search(hash, seqPos2);
       stSortedSetIterator *it2 = stSortedSet_getIterator(column2);
       stIntTuple *seqPos3;
       while((seqPos3 = stSortedSet_getNext(it2)) != NULL) {
           stSortedSet_insert(column1, seqPos3);
           stHash_insert(hash, seqPos3, column1);
       }
       stSortedSet_destructIterator(it2);
       stSortedSet_destruct(column2);
       //Cleanup loop.
       stIntTuple_destruct(seqPos1);
       stIntTuple_destruct(seqPos2);
    }
    stList_destructIterator(it);
    return hash;
}

/*
 * Function does the actual depth first search to detect if the thing has an acyclic ordering.
 */
static int32_t dfs(stHash *adjacencyList, stIntTuple *seqPos,
                               stSortedSet *started, stSortedSet *done) {
    if(stSortedSet_search(started, seqPos) != NULL) {
        if(stSortedSet_search(done, seqPos) == NULL) {
            //We have detected a cycle
            return 1;
        }
        //We have already explored this area, but no cycle.
        return 0;
    }
    stSortedSet_insert(started, seqPos);
    stSortedSet *column = stHash_search(adjacencyList, seqPos);
    stSortedSetIterator *it = stSortedSet_getIterator(column);
    stIntTuple *seqPos2;
    int32_t i=0;
    while((seqPos2 = stSortedSet_getNext(it)) != NULL) {
        i = i || dfs(adjacencyList, seqPos2, started, done);
    }
    stSortedSet_destructIterator(it);
    stSortedSet_insert(done, seqPos);
    return i;
}

/*
 * Uses the functions above to build an adjacency list, then by DFS attempts to create
 * a valid topological sort, returning non-zero if the graph contains a cycle.
 */
static int32_t containsACycle(stList *pairs, int32_t sequenceNumber) {
    //Build an adjacency list structure..
    stHash *adjacencyList = buildAdjacencyList(pairs, sequenceNumber);

    //Do a topological sort of the adjacency list
    stSortedSet *started = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, (void (*)(void *))stIntTuple_destruct);
    stSortedSet *done = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, NULL);
    int32_t cyclic = 0;
    for(int32_t seq=0; seq<sequenceNumber; seq++) {
        stIntTuple *seqPos = stIntTuple_construct(2, seq, 0);
        cyclic = cyclic || dfs(adjacencyList, seqPos, started, done);
        stIntTuple_destruct(seqPos);
    }

    //cleanup
    stHash_destruct(adjacencyList);
    stSortedSet_destruct(started);
    stSortedSet_destruct(done);

    return cyclic;
}


static void test_stPosetAlignment_addAndIsPossible(CuTest *testCase) {
    for(int32_t trial=0; trial<100; trial++) {
        setup();

        //Make random number of sequences.
        stList *sequenceLengths = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
        for(int32_t i=0; i<sequenceNumber; i++) {
            stList_append(sequenceLengths, stIntTuple_construct(1, st_randomInt(0, MAX_SEQUENCE_SIZE)));
        }

        //Propose random alignment pairs...
        stList *pairs = stList_construct3(0, (void(*)(void *))stIntTuple_destruct);
        int32_t maxAlignedPairs = st_randomInt(0, MAX_ALIGNMENTS);
        for(int32_t i=0; i<maxAlignedPairs; i++) {
            int32_t seq1 = st_randomInt(0, sequenceNumber);
            int32_t position1 = st_randomInt(0, stIntTuple_getPosition(stList_get(sequenceLengths, seq1), 0));
            int32_t seq2 = st_randomInt(0, sequenceNumber);
            int32_t position2 = st_randomInt(0, stIntTuple_getPosition(stList_get(sequenceLengths, seq2), 0));
            if(seq1 != seq2) {
                stList_append(pairs, stIntTuple_construct(4, seq1, position1, seq2, position2));
                if(stPosetAlignment_isPossible(posetAlignment, seq1, position1, seq2, position2)) {
                    //For each accepted pair check it doesn't create a cycle.
                    CuAssertTrue(testCase, !containsACycle(pairs, sequenceNumber));
                    CuAssertTrue(testCase, stPosetAlignment_add(posetAlignment, seq1, position1, seq2, position2));
                }
                else {
                    //For each rejected pair check it creates a cycle..
                    CuAssertTrue(testCase, containsACycle(pairs, sequenceNumber));
                    CuAssertTrue(testCase, !stPosetAlignment_isPossible(posetAlignment, seq1, position1, seq2, position2));
                    stList_pop(pairs); //remove the pair which created the cycle.
                    CuAssertTrue(testCase, !containsACycle(pairs, sequenceNumber)); //Check we're back to being okay..
                }
            }
        }

        //Cleanup
        stList_destruct(sequenceLengths);
        stList_destruct(pairs);
        teardown();
        st_logInfo("Passed a random ordering test with %i sequences and %i aligned pairs\n", sequenceNumber, maxAlignedPairs);
    }
}

CuSuite* stPosetAlignmentTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stPosetAlignment_addAndIsPossible);
    SUITE_ADD_TEST(suite, test_stPosetAlignment_getSequenceNumber);
    return suite;
}
