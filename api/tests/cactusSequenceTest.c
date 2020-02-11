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

static void cactusSequenceTestTeardown(CuTest* testCase) {
    if(!nestedTest && cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
        cactusDisk = NULL;
        flower = NULL;
        event = NULL;
        eventTree = NULL;
        metaSequence = NULL;
        sequence = NULL;
    }
}

static void cactusSequenceTestSetup(CuTest* testCase) {
	if(!nestedTest) {
		cactusSequenceTestTeardown(testCase);
		cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
		flower = flower_construct(cactusDisk);
		eventTree = eventTree_construct2(cactusDisk);
		event = eventTree_getRootEvent(eventTree);
		metaSequence = metaSequence_construct(1, 10, sequenceString,
						   headerString, event_getName(event), cactusDisk);
		sequence = sequence_construct(metaSequence, flower);
	}
}

void testSequence_construct(CuTest* testCase) {
	cactusSequenceTestSetup(testCase);
	CuAssertTrue(testCase, sequence != NULL);
	cactusSequenceTestTeardown(testCase);
}

void testSequence_getMetaSequence(CuTest* testCase) {
	cactusSequenceTestSetup(testCase);
	CuAssertTrue(testCase, metaSequence == sequence_getMetaSequence(sequence));
	cactusSequenceTestTeardown(testCase);
}

void testSequence_getStart(CuTest* testCase) {
	cactusSequenceTestSetup(testCase);
	CuAssertIntEquals(testCase, 1, sequence_getStart(sequence));
	cactusSequenceTestTeardown(testCase);
}

void testSequence_getLength(CuTest* testCase) {
	cactusSequenceTestSetup(testCase);
	CuAssertIntEquals(testCase, 10, sequence_getLength(sequence));
	cactusSequenceTestTeardown(testCase);
}

void testSequence_getName(CuTest* testCase) {
	cactusSequenceTestSetup(testCase);
	CuAssertTrue(testCase, sequence_getName(sequence) != NULL_NAME);
	CuAssertTrue(testCase, flower_getSequence(flower, sequence_getName(sequence)) == sequence);
	cactusSequenceTestTeardown(testCase);
}

void testSequence_getEvent(CuTest* testCase) {
	cactusSequenceTestSetup(testCase);
	CuAssertTrue(testCase, sequence_getEvent(sequence) == event);
	cactusSequenceTestTeardown(testCase);
}

void testSequence_getString(CuTest* testCase) {
	cactusSequenceTestSetup(testCase);
	int64_t i, j;
	for(i=1; i<11; i++) {
		for(j=11-i; j>=0; j--) {
			CuAssertStrEquals(testCase, metaSequence_getString(metaSequence, i, j, 1), sequence_getString(sequence, i, j, 1));
			CuAssertStrEquals(testCase, metaSequence_getString(metaSequence, i, j, 0), sequence_getString(sequence, i, j, 0));
		}
	}
	cactusSequenceTestTeardown(testCase);
}

static char *getRandomDNASequence(int64_t minSequenceLength, int64_t maxSequenceLength) {
    int64_t stringLength = st_randomInt(minSequenceLength, maxSequenceLength);
    char *string = st_malloc(sizeof(char) * (stringLength + 1));
    string[stringLength] = '\0';
    for (int64_t j = 0; j < stringLength; j++) {
        char cA[10] = { 'a', 'c', 'g', 't', 'A', 'C', 'G', 'T', 'n', 'N' };
        string[j] = cA[st_randomInt(0, 10)];
    }
    return string;
}

void testSequence_addAndGetBigStringsP(CuTest* testCase,
        bool preCacheSequences, bool reopenCactusDisk,
        int64_t minSequenceLength, int64_t maxSequenceLength, int64_t testNo) {
    for(int64_t i=0; i<testNo; i++) {
        cactusSequenceTestSetup(testCase);
        //Create a bunch of sequences
        stList *strings = stList_construct3(0, free);
        stList *sequenceNames = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
        int64_t coordinateStart = st_randomInt(0, 100);
        Name flowerName = flower_getName(flower);
        do {
            char *string = getRandomDNASequence(minSequenceLength, maxSequenceLength);
            stList_append(strings, string);
            stList_append(sequenceNames, stIntTuple_construct1(sequence_getName(sequence_construct(metaSequence_construct(coordinateStart, strlen(string), string,
                                               "Hello I am header", event_getName(event), cactusDisk), flower))));

        } while(st_random() > 0.5);

        if(reopenCactusDisk) {
            cactusDisk_write(cactusDisk);
            cactusDisk_destruct(cactusDisk);
            stKVDatabaseConf *conf = testCommon_getTemporaryKVDatabaseConf(testCase->name);
            cactusDisk = cactusDisk_construct(conf, false, true);
            flower = cactusDisk_getFlower(cactusDisk, flowerName);
            testCommon_deleteTemporaryKVDatabase(testCase->name);
        }

        stList *sequences = stList_construct(); //Get the meta-sequences from the reopened cactusDisk.
        for(int64_t j=0; j<stList_length(sequenceNames); j++) {
            stList_append(sequences, flower_getSequence(flower, stIntTuple_get(stList_get(sequenceNames, j), 0)));
        }

        if(preCacheSequences) {
            for(int64_t j=0; j<stList_length(sequences); j++) { //Ensures the sequences have caps representing their ends.
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
        bool first = 1;
        while(st_random() > 0.1) {
            int64_t j = st_randomInt(0, stList_length(strings));
            Sequence *sequence = stList_get(sequences, j);
            if(sequence_getLength(sequence) == 0) {
                continue;
            }
            char *string = stList_get(strings, j);
            CuAssertIntEquals(testCase, strlen(string), sequence_getLength(sequence));
            //Choose a random interval to request
            int64_t start = first ? 0 : st_randomInt(0, strlen(string));
            int64_t length = first ? strlen(string) : st_randomInt(0, strlen(string)-start);
            st_logInfo("Getting interval start: %" PRIi64 " length %" PRIi64 " of sequence of length %" PRIi64 " \n", start, length, sequence_getLength(sequence));
            first = 0;
            CuAssertTrue(testCase, start >= 0);
            CuAssertTrue(testCase, start + length <= strlen(string));
            char *subString = stString_getSubString(string, start, length);
            CuAssertIntEquals(testCase, length, strlen(subString));
            bool strand = st_random() > 0.5;
            if(!strand) {
                char *subString2 = stString_reverseComplementString(subString);
                CuAssertIntEquals(testCase, strlen(subString2), strlen(subString));
                free(subString);
                subString = subString2;
            }
            char *subSequence = NULL;
            if(preCacheSequences) {
                subSequence = cactusDisk_getStringFromCache(cactusDisk, sequence_getMetaSequence(sequence)->stringName,
                        start, length, strand);
            }
            if(preCacheSequences || subSequence == NULL) {
                subSequence = sequence_getString(sequence, coordinateStart + start, length, strand);
            }
            CuAssertTrue(testCase, subSequence != NULL);
            for(int64_t k=0; k<length; k++) {
                CuAssertIntEquals(testCase, subString[k], subSequence[k]);
            }
            //CuAssertStrEquals(testCase, subString, subSequence);
            free(subString);
            free(subSequence);
        }
        stList_destruct(sequences);
        stList_destruct(sequenceNames);
        cactusSequenceTestTeardown(testCase);
        stList_destruct(strings);
    }
}

void testSequence_addAndGetBigStrings(CuTest* testCase) {
    testSequence_addAndGetBigStringsP(testCase, 0, 0, 0, 100000, 100);
}

void testSequence_addAndGetBigStrings_preCacheSequences(CuTest* testCase) {
    testSequence_addAndGetBigStringsP(testCase, 1, 0, 0, 100000, 100);
}

void testSequence_addAndGetBigStrings_reopenCactusDisk(CuTest* testCase) {
    testSequence_addAndGetBigStringsP(testCase, 0, 1, 0, 100000, 100);
}

void testSequence_addAndGetBigStrings_preCacheSequences_reopenCactusDisk(CuTest* testCase) {
    testSequence_addAndGetBigStringsP(testCase, 1, 1, 0, 100000, 100);
}

void testSequence_addAndGetBigStrings_massive(CuTest* testCase) {
    testSequence_addAndGetBigStringsP(testCase, 1, 1, 5000000, 10000000, 5);
}

void testSequence_getHeader(CuTest* testCase) {
	cactusSequenceTestSetup(testCase);
	CuAssertStrEquals(testCase, headerString, sequence_getHeader(sequence));
	cactusSequenceTestTeardown(testCase);
}

void testSequence_getFlower(CuTest* testCase) {
	cactusSequenceTestSetup(testCase);
	CuAssertTrue(testCase, sequence_getFlower(sequence) == flower);
	cactusSequenceTestTeardown(testCase);
}

void testSequence_serialisation(CuTest* testCase) {
	cactusSequenceTestSetup(testCase);
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
	cactusSequenceTestTeardown(testCase);
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
	SUITE_ADD_TEST(suite, testSequence_addAndGetBigStrings_reopenCactusDisk);
	SUITE_ADD_TEST(suite, testSequence_addAndGetBigStrings_preCacheSequences_reopenCactusDisk);
	SUITE_ADD_TEST(suite, testSequence_addAndGetBigStrings_massive);
	SUITE_ADD_TEST(suite, testSequence_getHeader);
	SUITE_ADD_TEST(suite, testSequence_getFlower);
	SUITE_ADD_TEST(suite, testSequence_serialisation);
	SUITE_ADD_TEST(suite, testSequence_construct);
	return suite;
}
