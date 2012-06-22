/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"
#include <unistd.h>

static CactusDisk *cactusDisk = NULL;
static Flower *flower;
static Event *event;
static EventTree *eventTree;
static MetaSequence *metaSequence;
static Sequence *sequence;
static const char *sequenceString = "ACTGGCACTG";
static const char *headerString = ">one";

static bool nestedTest = 0;

void cactusSequenceTestTeardown() {
	if(!nestedTest && cactusDisk != NULL) {
		testCommon_deleteTemporaryCactusDisk(cactusDisk);
		cactusDisk = NULL;
		flower = NULL;
		event = NULL;
		eventTree = NULL;
		metaSequence = NULL;
		sequence = NULL;
	}
}

void cactusSequenceTestSetup() {
	if(!nestedTest) {
		cactusSequenceTestTeardown();
		cactusDisk = testCommon_getTemporaryCactusDisk();
		flower = flower_construct(cactusDisk);
		eventTree = eventTree_construct2(flower);
		event = eventTree_getRootEvent(eventTree);
		metaSequence = metaSequence_construct(1, 10, sequenceString,
						   headerString, event_getName(event), cactusDisk);
		sequence = sequence_construct(metaSequence, flower);
	}
}

void testSequence_construct(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertTrue(testCase, sequence != NULL);
	cactusSequenceTestTeardown();
}

void testSequence_getMetaSequence(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertTrue(testCase, metaSequence == sequence_getMetaSequence(sequence));
	cactusSequenceTestTeardown();
}

void testSequence_getStart(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertIntEquals(testCase, 1, sequence_getStart(sequence));
	cactusSequenceTestTeardown();
}

void testSequence_getLength(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertIntEquals(testCase, 10, sequence_getLength(sequence));
	cactusSequenceTestTeardown();
}

void testSequence_getName(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertTrue(testCase, sequence_getName(sequence) != NULL_NAME);
	CuAssertTrue(testCase, flower_getSequence(flower, sequence_getName(sequence)) == sequence);
	cactusSequenceTestTeardown();
}

void testSequence_getEvent(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertTrue(testCase, sequence_getEvent(sequence) == event);
	cactusSequenceTestTeardown();
}

void testSequence_getString(CuTest* testCase) {
	cactusSequenceTestSetup();
	int32_t i, j;
	for(i=1; i<11; i++) {
		for(j=11-i; j>=0; j--) {
			CuAssertStrEquals(testCase, metaSequence_getString(metaSequence, i, j, 1), sequence_getString(sequence, i, j, 1));
			CuAssertStrEquals(testCase, metaSequence_getString(metaSequence, i, j, 0), sequence_getString(sequence, i, j, 0));
		}
	}
	cactusSequenceTestTeardown();
}

static char *getRandomDNASequence() {
    int32_t stringLength = st_randomInt(0, 100000);
    char *string = st_malloc(sizeof(char) * (stringLength + 1));
    string[stringLength] = '\0';
    for (int32_t j = 0; j < stringLength; j++) {
        char cA[10] = { 'a', 'c', 'g', 't', 'A', 'C', 'G', 'T', 'n', 'N' };
        string[j] = cA[st_randomInt(0, 10)];
    }
    return string;
}

void testSequence_addAndGetBigStringsP(CuTest* testCase, bool preCacheSequences) {
    for(int32_t i=0; i<100; i++) {
        cactusSequenceTestSetup();
        //Create a bunch of sequences
        stList *strings = stList_construct3(0, free);
        stList *sequences = stList_construct();
        int32_t coordinateStart = st_randomInt(0, 100);
        do {
            char *string = getRandomDNASequence();
            stList_append(strings, string);
            metaSequence = metaSequence_construct(coordinateStart, strlen(string), string,
                                               "Hello I am header", event_getName(event), cactusDisk);
            stList_append(sequences, sequence_construct(metaSequence, flower));
        } while(st_random() > 0.5);

        if(preCacheSequences) {
            for(int32_t j=0; j<stList_length(sequences); j++) { //Ensures the sequences have caps representing their ends.
                Sequence *sequence = stList_get(sequences, j);
                End *end = end_construct(st_random() > 0.5, flower);
                assert(end_getSide(end));
                Cap *cap = cap_construct2(end_getReverse(end), sequence_getStart(sequence)-1, 1, sequence);
                Cap *adjacentCap = cap_construct2(end, sequence_getStart(sequence) + sequence_getLength(sequence), 1, sequence);
                cap_makeAdjacent(cap, adjacentCap);
            }
            stList *flowerList = stList_construct();
            stList_append(flowerList, flower);
            cactusDisk_preCacheStrings(cactusDisk, flowerList);
            stList_destruct(flowerList);
        }
        //Do different requests for portions of the strings
        while(st_random() > 0.01) {
            int32_t j = st_randomInt(0, stList_length(strings));
            Sequence *sequence = stList_get(sequences, j);
            if(sequence_getLength(sequence) == 0) {
                continue;
            }
            char *string = stList_get(strings, j);
            CuAssertIntEquals(testCase, strlen(string), sequence_getLength(sequence));
            //Choose a random interval to request
            int64_t start = st_randomInt(0, strlen(string));
            int64_t length = st_randomInt(0, strlen(string)-start);
            CuAssertTrue(testCase, start >= 0);
            CuAssertTrue(testCase, start + length <= strlen(string));
            char *subString = stString_getSubString(string, start, length);
            CuAssertIntEquals(testCase, length, strlen(subString));
            bool strand = st_random() > 0.5;
            if(!strand) {
                char *subString2 = cactusMisc_reverseComplementString(subString);
                CuAssertIntEquals(testCase, strlen(subString2), strlen(subString));
                free(subString);
                subString = subString2;
            }
            char *subSequence = NULL;
            if(preCacheSequences) {
                subSequence = cactusDisk_getStringFromCache(cactusDisk, sequence_getMetaSequence(sequence)->stringName,
                        start, length, strand);
            }
            else {
                subSequence = sequence_getString(sequence, coordinateStart + start, length, strand);
            }
            assert(subSequence != NULL);
            CuAssertStrEquals(testCase, subString, subSequence);
            free(subString);
            free(subSequence);
        }
        stList_destruct(sequences);
        cactusSequenceTestTeardown();
        stList_destruct(strings);
    }
}

void testSequence_addAndGetBigStrings(CuTest* testCase) {
    testSequence_addAndGetBigStringsP(testCase, 0);
}

void testSequence_addAndGetBigStrings_preCacheSequences(CuTest* testCase) {
    testSequence_addAndGetBigStringsP(testCase, 1);
}

void testSequence_getHeader(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertStrEquals(testCase, headerString, sequence_getHeader(sequence));
	cactusSequenceTestTeardown();
}

void testSequence_getFlower(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertTrue(testCase, sequence_getFlower(sequence) == flower);
	cactusSequenceTestTeardown();
}

void testSequence_serialisation(CuTest* testCase) {
	cactusSequenceTestSetup();
	int64_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(sequence,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))sequence_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	sequence_destruct(sequence);
	void *vA2 = vA;
	sequence = sequence_loadFromBinaryRepresentation(&vA2, flower);
	nestedTest = 1;
	testSequence_getMetaSequence(testCase);
	testSequence_getStart(testCase);
	testSequence_getLength(testCase);
	testSequence_getName(testCase);
	testSequence_getEvent(testCase);
	testSequence_getString(testCase);
	testSequence_getHeader(testCase);
	testSequence_getFlower(testCase);
	nestedTest = 0;
	cactusSequenceTestTeardown();
}

CuSuite* cactusSequenceTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testSequence_getMetaSequence);
	SUITE_ADD_TEST(suite, testSequence_getStart);
	SUITE_ADD_TEST(suite, testSequence_getLength);
	SUITE_ADD_TEST(suite, testSequence_getName);
	SUITE_ADD_TEST(suite, testSequence_getEvent);
	SUITE_ADD_TEST(suite, testSequence_getString);
	SUITE_ADD_TEST(suite, testSequence_addAndGetBigStrings);
	SUITE_ADD_TEST(suite, testSequence_addAndGetBigStrings_preCacheSequences);
	SUITE_ADD_TEST(suite, testSequence_getHeader);
	SUITE_ADD_TEST(suite, testSequence_getFlower);
	SUITE_ADD_TEST(suite, testSequence_serialisation);
	SUITE_ADD_TEST(suite, testSequence_construct);
	return suite;
}
