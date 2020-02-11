/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stPinchIterator.h"
#include "pairwiseAlignment.h"
#include <math.h>

static void testIterator(CuTest *testCase, stPinchIterator *pinchIterator, stList *randomPairwiseAlignments) {
    for (int64_t trim = 0; trim < 10; trim++) {
        stPinchIterator_setTrim(pinchIterator, trim);
        //Test get next
        for (int64_t i = 0; i < stList_length(randomPairwiseAlignments); i++) {
            struct PairwiseAlignment *pairwiseAlignment = stList_get(randomPairwiseAlignments, i);
            int64_t contigX = -1, contigY = -1;
            sscanf(pairwiseAlignment->contig1, "%" PRIi64 "", &contigX);
            sscanf(pairwiseAlignment->contig2, "%" PRIi64 "", &contigY);
            int64_t x = pairwiseAlignment->start1;
            int64_t y = pairwiseAlignment->start2;
            for (int64_t j = 0; j < pairwiseAlignment->operationList->length; j++) {
                struct AlignmentOperation *op = pairwiseAlignment->operationList->list[j];
                if (op->opType == PAIRWISE_MATCH) {
                    if (op->length > 2 * trim) {
                        stPinch *pinch = stPinchIterator_getNext(pinchIterator);
                        CuAssertTrue(testCase, pinch != NULL);
                        CuAssertIntEquals(testCase, contigX, pinch->name1);
                        CuAssertIntEquals(testCase, contigY, pinch->name2);
                        CuAssertIntEquals(testCase, (pairwiseAlignment->strand1 ? x : x - op->length) + trim, pinch->start1);
                        CuAssertIntEquals(testCase, (pairwiseAlignment->strand2 ? y : y - op->length) + trim, pinch->start2);
                        CuAssertIntEquals(testCase, op->length - 2 * trim, pinch->length);
                        CuAssertTrue(testCase, pinch->length > 0);
                        CuAssertIntEquals(testCase, pairwiseAlignment->strand1 == pairwiseAlignment->strand2, pinch->strand);
                    }
                }
                if (op->opType != PAIRWISE_INDEL_Y) {
                    x += pairwiseAlignment->strand1 ? op->length : -op->length;
                }
                if (op->opType != PAIRWISE_INDEL_X) {
                    y += pairwiseAlignment->strand2 ? op->length : -op->length;
                }
            }
        }
        CuAssertPtrEquals(testCase, NULL, stPinchIterator_getNext(pinchIterator));
        CuAssertPtrEquals(testCase, NULL, stPinchIterator_getNext(pinchIterator));
        CuAssertPtrEquals(testCase, NULL, stPinchIterator_getNext(pinchIterator));
        stPinchIterator_reset(pinchIterator);
    }
}

static stList *getRandomPairwiseAlignments() {
    stList *pairwiseAlignments = stList_construct3(0, (void(*)(void *)) destructPairwiseAlignment);
    int64_t randomAlignmentNumber = st_randomInt(0, 10);
    for (int64_t i = 0; i < randomAlignmentNumber; i++) {
        char *contig1 = stString_print("%" PRIi64 "", i);
        char *contig2 = stString_print("%" PRIi64 "", i * 10);
        int64_t start1 = st_randomInt(100000, 1000000);
        int64_t start2 = st_randomInt(100000, 1000000);
        int64_t strand1 = st_random() > 0.5;
        int64_t strand2 = st_random() > 0.5;
        int64_t end1 = start1;
        int64_t end2 = start2;
        struct List *operationList = constructEmptyList(0, NULL);
        while (st_random() > 0.1) {
            int64_t length = st_randomInt(0, 10);
            int64_t type = st_randomInt(0, 3);
            assert(type < 3);
            listAppend(operationList, constructAlignmentOperation(type, length, 0));
            if (type != PAIRWISE_INDEL_Y) {
                end1 += strand1 ? length : -length;
            }
            if (type != PAIRWISE_INDEL_X) {
                end2 += strand2 ? length : -length;
            }
        }
        stList_append(pairwiseAlignments,
                constructPairwiseAlignment(contig1, start1, end1, strand1, contig2, start2, end2, strand2, 0.0, operationList));
        free(contig1);
        free(contig2);
    }
    return pairwiseAlignments;
}

static void testPinchIteratorFromFile(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        stList *pairwiseAlignments = getRandomPairwiseAlignments();
        st_logInfo("Doing a random pinch iterator from file test %" PRIi64 " with %" PRIi64 " alignments\n", test, stList_length(pairwiseAlignments));
        //Put alignments in a file
        char *tempFile = "tempFileForPinchIteratorTest.cig";
        FILE *fileHandle = fopen(tempFile, "w");
        for (int64_t i = 0; i < stList_length(pairwiseAlignments); i++) {
            cigarWrite(fileHandle, stList_get(pairwiseAlignments, i), 0);
        }
        fclose(fileHandle);
        //Get an iterator
        stPinchIterator *pinchIterator = stPinchIterator_constructFromFile(tempFile);
        //Now test it
        testIterator(testCase, pinchIterator, pairwiseAlignments);
        //Cleanup
        stPinchIterator_destruct(pinchIterator);
        stFile_rmtree(tempFile);
        stList_destruct(pairwiseAlignments);
    }
}

static void testPinchIteratorFromList(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        stList *pairwiseAlignments = getRandomPairwiseAlignments();
        st_logInfo("Doing a random pinch iterator from list test %" PRIi64 " with %" PRIi64 " alignments\n", test, stList_length(pairwiseAlignments));
        //Get an iterator
        stPinchIterator *pinchIterator = stPinchIterator_constructFromList(pairwiseAlignments);
        //Now test it
        testIterator(testCase, pinchIterator, pairwiseAlignments);
        //Cleanup
        stPinchIterator_destruct(pinchIterator);
        stList_destruct(pairwiseAlignments);
    }
}

CuSuite* pinchIteratorTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testPinchIteratorFromFile);
    SUITE_ADD_TEST(suite, testPinchIteratorFromList);
    return suite;
}
