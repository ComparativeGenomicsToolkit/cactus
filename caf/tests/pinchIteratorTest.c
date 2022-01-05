/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "stPinchIterator.h"
#include "pairwiseAlignment.h"
#include "paf.h"
#include <math.h>

static void testIterator(CuTest *testCase, stPinchIterator *pinchIterator, stList *randomPairwiseAlignments) {
    for (int64_t trim = 0; trim < 10; trim++) {
        stPinchIterator_setTrim(pinchIterator, trim);
        //Test get next
        stPinch pinchToFillOut;
        for (int64_t i = 0; i < stList_length(randomPairwiseAlignments); i++) {
            Paf *paf = stList_get(randomPairwiseAlignments, i);
            int64_t contigX = -1, contigY = -1;
            sscanf(paf->query_name, "%" PRIi64 "", &contigX);
            sscanf(paf->target_name, "%" PRIi64 "", &contigY);
            int64_t x = paf->same_strand ? paf->query_start : paf->query_end;
            int64_t y = paf->target_start;
            Cigar *c = paf->cigar;
            while(c != NULL) {
                if (c->op == match) {
                    if (c->length > 2 * trim) {
                        stPinch *pinch = stPinchIterator_getNext(pinchIterator, &pinchToFillOut);
                        CuAssertTrue(testCase, pinch != NULL);
                        CuAssertIntEquals(testCase, contigX, pinch->name1);
                        CuAssertIntEquals(testCase, contigY, pinch->name2);
                        CuAssertIntEquals(testCase, (paf->same_strand ? x : x - c->length) + trim, pinch->start1);
                        CuAssertIntEquals(testCase, y + trim, pinch->start2);
                        CuAssertIntEquals(testCase, c->length - 2 * trim, pinch->length);
                        CuAssertTrue(testCase, pinch->length > 0);
                        CuAssertIntEquals(testCase, paf->same_strand, pinch->strand);
                    }
                }
                if (c->op != query_delete) {
                    x += paf->same_strand ? c->length : -c->length;
                }
                if (c->op != query_insert) {
                    y += c->length;
                }
                c = c->next;
            }
        }
        CuAssertPtrEquals(testCase, NULL, stPinchIterator_getNext(pinchIterator, &pinchToFillOut));
        CuAssertPtrEquals(testCase, NULL, stPinchIterator_getNext(pinchIterator, &pinchToFillOut));
        CuAssertPtrEquals(testCase, NULL, stPinchIterator_getNext(pinchIterator, &pinchToFillOut));
        stPinchIterator_reset(pinchIterator);
    }
}

static stList *getRandomPairwiseAlignments() {
    stList *pafs = stList_construct3(0, (void(*)(void *)) paf_destruct);
    int64_t randomAlignmentNumber = st_randomInt(0, 10);
    for (int64_t i = 0; i < randomAlignmentNumber; i++) {
        Paf *paf = st_calloc(1, sizeof(Paf));
        paf->query_name = stString_print("%" PRIi64 "", i);
        paf->target_name = stString_print("%" PRIi64 "", i * 10);
        paf->query_start = st_randomInt(100000, 1000000);
        paf->target_start = st_randomInt(100000, 1000000);
        paf->same_strand = st_random() > 0.5;
        int64_t i = paf->query_start, j = paf->target_start;
        Cigar **pc = &(paf->cigar);
        do {
            Cigar *c = st_calloc(1, sizeof(Cigar));
            c->length = st_randomInt(1, 10);
            c->op = st_random() > 0.3 ? (st_random() > 0.5 ? match : query_insert): query_delete;
            if (c->op != query_delete) {
                i += c->length;
            }
            if (c->op != query_insert) {
                j += c->length;
            }
            *pc = c;
            pc = &(c->next);
        } while(st_random() > 0.1 || paf->query_start == i || paf->target_start == j);
        paf->query_end = i;
        paf->target_end = j;
        paf->query_length = i;
        paf->target_length = j;
        paf_check(paf);
        stList_append(pafs, paf);
    }
    return pafs;
}

static void testPinchIteratorFromFile(CuTest *testCase) {
    for (int64_t test = 0; test < 100; test++) {
        stList *pairwiseAlignments = getRandomPairwiseAlignments();
        st_logInfo("Doing a random pinch iterator from file test %" PRIi64 " with %" PRIi64 " alignments\n", test, stList_length(pairwiseAlignments));
        //Put alignments in a file
        char *tempFile = "tempFileForPinchIteratorTest.cig";
        FILE *fileHandle = fopen(tempFile, "w");
        assert(fileHandle != NULL);
        write_pafs(fileHandle, pairwiseAlignments);
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

CuSuite* pinchIteratorTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testPinchIteratorFromFile);
    return suite;
}
