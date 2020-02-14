/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk = NULL;
static Flower *flower;

static Event *rootEvent;
static Event *internalEvent;
static Event *leafEvent1;
static Event *leafEvent2;

static EventTree *eventTree;

static bool nestedTest = 0;

static void cactusEventTreeTestTeardown(CuTest* testCase) {
	if(!nestedTest && cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
		cactusDisk = NULL;
		eventTree = NULL;
		rootEvent = NULL;
		internalEvent = NULL;
		leafEvent1 = NULL;
		leafEvent2 = NULL;
	}
}

static void cactusEventTreeTestSetup(CuTest* testCase) {
	if(!nestedTest) {
		cactusEventTreeTestTeardown(testCase);
		cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
		flower = flower_construct(cactusDisk);

		eventTree = eventTree_construct2(cactusDisk);
		rootEvent = eventTree_getRootEvent(eventTree);
		internalEvent = event_construct3("INTERNAL", 0.5, rootEvent, eventTree);
		leafEvent1 = event_construct3("LEAF1", 0.2, internalEvent, eventTree);
		leafEvent2 = event_construct3("LEAF2", 1.3, internalEvent, eventTree);
		event_setOutgroupStatus(leafEvent1, 1);
	}
}

void testEventTree_construct(CuTest* testCase) {
	nestedTest = 0;
	cactusEventTreeTestSetup(testCase);
	CuAssertTrue(testCase, eventTree != NULL);
	cactusEventTreeTestTeardown(testCase);
}

static int64_t unaryEventFunction(Event *event) {
	assert(event != NULL);
	assert(event_getChildNumber(event) == 1);
	return 1;
}

void testEventTree_copyConstruct(CuTest* testCase) {
	cactusEventTreeTestSetup(testCase);
	flower_construct(cactusDisk);
	EventTree *eventTree2 = eventTree_copyConstruct(eventTree, unaryEventFunction);
	CuAssertIntEquals(testCase, eventTree_getEventNumber(eventTree), eventTree_getEventNumber(eventTree2));
	CuAssertTrue(testCase, event_getName(eventTree_getEvent(eventTree2, event_getName(rootEvent))) == event_getName(rootEvent));
	CuAssertTrue(testCase, event_getName(eventTree_getEvent(eventTree2, event_getName(internalEvent))) == event_getName(internalEvent));
	CuAssertTrue(testCase, event_getName(eventTree_getEvent(eventTree2, event_getName(leafEvent1))) == event_getName(leafEvent1));
	CuAssertTrue(testCase, event_getName(eventTree_getEvent(eventTree2, event_getName(leafEvent2))) == event_getName(leafEvent2));
	CuAssertTrue(testCase, event_isOutgroup(eventTree_getEvent(eventTree2, event_getName(leafEvent1))));
	CuAssertTrue(testCase, !event_isOutgroup(eventTree_getEvent(eventTree2, event_getName(leafEvent2))));
	cactusEventTreeTestTeardown(testCase);
}

void testEventTree_getRootEvent(CuTest* testCase) {
	cactusEventTreeTestSetup(testCase);
	CuAssertTrue(testCase, eventTree_getRootEvent(eventTree) == rootEvent);
	cactusEventTreeTestTeardown(testCase);
}

void testEventTree_getEvent(CuTest* testCase) {
	cactusEventTreeTestSetup(testCase);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree, event_getName(rootEvent)) == rootEvent);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree, event_getName(internalEvent)) == internalEvent);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree, event_getName(leafEvent1)) == leafEvent1);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree, event_getName(leafEvent2)) == leafEvent2);
	cactusEventTreeTestTeardown(testCase);
}

void testEventTree_getCommonAncestor(CuTest* testCase) {
	cactusEventTreeTestSetup(testCase);
	//self
	CuAssertTrue(testCase, eventTree_getCommonAncestor(leafEvent1, leafEvent1) == leafEvent1);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(leafEvent2, leafEvent2) == leafEvent2);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(internalEvent, internalEvent) == internalEvent);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(rootEvent, rootEvent) == rootEvent);

	CuAssertTrue(testCase, eventTree_getCommonAncestor(leafEvent1, leafEvent2) == internalEvent);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(leafEvent2, leafEvent1) == internalEvent);

	CuAssertTrue(testCase, eventTree_getCommonAncestor(leafEvent1, internalEvent) == internalEvent);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(internalEvent, leafEvent1) == internalEvent);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(leafEvent2, internalEvent) == internalEvent);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(internalEvent, leafEvent2) == internalEvent);

	CuAssertTrue(testCase, eventTree_getCommonAncestor(leafEvent1, rootEvent) == rootEvent);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(rootEvent, leafEvent1) == rootEvent);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(leafEvent2, rootEvent) == rootEvent);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(rootEvent, leafEvent2) == rootEvent);

	CuAssertTrue(testCase, eventTree_getCommonAncestor(internalEvent, rootEvent) == rootEvent);
	CuAssertTrue(testCase, eventTree_getCommonAncestor(rootEvent, internalEvent) == rootEvent);

	cactusEventTreeTestTeardown(testCase);
}

void testEventTree_getEventNumber(CuTest* testCase) {
	cactusEventTreeTestSetup(testCase);
	CuAssertIntEquals(testCase, 4, eventTree_getEventNumber(eventTree));
	cactusEventTreeTestTeardown(testCase);
}

void testEventTree_getFirst(CuTest* testCase) {
	cactusEventTreeTestSetup(testCase);
	CuAssertTrue(testCase, eventTree_getFirst(eventTree) == rootEvent);
	cactusEventTreeTestTeardown(testCase);
}

void testEventTree_iterator(CuTest* testCase) {
	cactusEventTreeTestSetup(testCase);
	EventTree_Iterator *iterator = eventTree_getIterator(eventTree);
	CuAssertTrue(testCase, eventTree_getNext(iterator) == rootEvent);
	CuAssertTrue(testCase, eventTree_getNext(iterator) == internalEvent);
	EventTree_Iterator *iterator2 = eventTree_copyIterator(iterator);
	CuAssertTrue(testCase, eventTree_getNext(iterator) == leafEvent1);
	CuAssertTrue(testCase, eventTree_getNext(iterator) == leafEvent2);
	CuAssertTrue(testCase, eventTree_getNext(iterator) == NULL);
	CuAssertTrue(testCase, eventTree_getNext(iterator2) == leafEvent1);
	CuAssertTrue(testCase, eventTree_getNext(iterator2) == leafEvent2);
	CuAssertTrue(testCase, eventTree_getNext(iterator2) == NULL);

	CuAssertTrue(testCase, eventTree_getPrevious(iterator2) == leafEvent2);
	CuAssertTrue(testCase, eventTree_getPrevious(iterator2) == leafEvent1);
	CuAssertTrue(testCase, eventTree_getPrevious(iterator2) == internalEvent);
	CuAssertTrue(testCase, eventTree_getPrevious(iterator2) == rootEvent);
	CuAssertTrue(testCase, eventTree_getPrevious(iterator2) == NULL);

	eventTree_destructIterator(iterator);
	eventTree_destructIterator(iterator2);

	cactusEventTreeTestTeardown(testCase);
}

void testEventTree_makeNewickString(CuTest* testCase) {
	cactusEventTreeTestSetup(testCase);
	CuAssertStrEquals(testCase, "((LEAF1:0.2,LEAF2:1.3)INTERNAL:0.5)ROOT:9.22337e+18;", eventTree_makeNewickString(eventTree));
	cactusEventTreeTestTeardown(testCase);
}

void testEventTree_serialisation(CuTest* testCase) {
	cactusEventTreeTestSetup(testCase);
	int64_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(eventTree,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))eventTree_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	eventTree_destruct(eventTree);
	void *vA2 = vA;
	eventTree = eventTree_loadFromBinaryRepresentation(&vA2, cactusDisk);
	rootEvent = eventTree_getRootEvent(eventTree);
	internalEvent = event_getChild(rootEvent, 0);
	leafEvent1 = event_getChild(internalEvent, 0);
	leafEvent2 = event_getChild(internalEvent, 1);
	free(vA);
	nestedTest = 1;
	testEventTree_copyConstruct(testCase);
	testEventTree_getRootEvent(testCase);
	testEventTree_getEvent(testCase);
	testEventTree_getCommonAncestor(testCase);
	testEventTree_getEventNumber(testCase);
	testEventTree_getFirst(testCase);
	testEventTree_makeNewickString(testCase);
	testEventTree_iterator(testCase);
	nestedTest = 0;
	cactusEventTreeTestTeardown(testCase);
}

CuSuite* cactusEventTreeTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testEventTree_copyConstruct);
	SUITE_ADD_TEST(suite, testEventTree_getRootEvent);
	SUITE_ADD_TEST(suite, testEventTree_getEvent);
	SUITE_ADD_TEST(suite, testEventTree_getCommonAncestor);
	SUITE_ADD_TEST(suite, testEventTree_getEventNumber);
	SUITE_ADD_TEST(suite, testEventTree_getFirst);
	SUITE_ADD_TEST(suite, testEventTree_iterator);
	SUITE_ADD_TEST(suite, testEventTree_makeNewickString);
	SUITE_ADD_TEST(suite, testEventTree_serialisation);
	SUITE_ADD_TEST(suite, testEventTree_construct);
	return suite;
}
