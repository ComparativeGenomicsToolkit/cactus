#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk = NULL;
static Flower *flower;

static MetaEvent *rootMetaEvent;
static MetaEvent *internalMetaEvent;
static MetaEvent *leafMetaEvent1;
static MetaEvent *leafMetaEvent2;

static Event *rootEvent;
static Event *internalEvent;
static Event *leafEvent1;
static Event *leafEvent2;

static EventTree *eventTree;

static bool nestedTest = 0;

static void cactusEventTreeTestTeardown() {
	if(!nestedTest && cactusDisk != NULL) {
		cactusDisk_destruct(cactusDisk);
		testCommon_deleteTemporaryCactusDisk();
		cactusDisk = NULL;
		eventTree = NULL;
		rootEvent = NULL;
		internalEvent = NULL;
		leafEvent1 = NULL;
		leafEvent2 = NULL;
	}
}

static void cactusEventTreeTestSetup() {
	if(!nestedTest) {
		cactusEventTreeTestTeardown();
		cactusDisk = cactusDisk_construct(testCommon_getTemporaryCactusDisk());
		flower = flower_construct(cactusDisk);

		rootMetaEvent = metaEvent_construct("ROOT", cactusDisk);
		internalMetaEvent = metaEvent_construct("INTERNAL", cactusDisk);
		leafMetaEvent1 = metaEvent_construct("LEAF1", cactusDisk);
		leafMetaEvent2 = metaEvent_construct("LEAF2", cactusDisk);
		eventTree = eventTree_construct(rootMetaEvent, flower);

		rootEvent = eventTree_getEvent(eventTree, metaEvent_getName(rootMetaEvent));
		internalEvent = event_construct(internalMetaEvent, 0.5, rootEvent, eventTree);
		leafEvent1 = event_construct(leafMetaEvent1, 0.2, internalEvent, eventTree);
		leafEvent2 = event_construct(leafMetaEvent2, 1.3, internalEvent, eventTree);
	}
}

void testEventTree_construct(CuTest* testCase) {
	nestedTest = 0;
	cactusEventTreeTestSetup();
	CuAssertTrue(testCase, eventTree != NULL);
	cactusEventTreeTestTeardown();
}

static int32_t unaryEventFunction(Event *event) {
	assert(event != NULL);
	assert(event_getChildNumber(event) == 1);
	return 1;
}

void testEventTree_copyConstruct(CuTest* testCase) {
	cactusEventTreeTestSetup();
	Flower *flower2 = flower_construct(cactusDisk);
	EventTree *eventTree2 = eventTree_copyConstruct(eventTree, flower2, unaryEventFunction);
	CuAssertIntEquals(testCase, eventTree_getEventNumber(eventTree), eventTree_getEventNumber(eventTree2));
	CuAssertTrue(testCase, event_getMetaEvent(eventTree_getEvent(eventTree2, event_getName(rootEvent))) == event_getMetaEvent(rootEvent));
	CuAssertTrue(testCase, event_getMetaEvent(eventTree_getEvent(eventTree2, event_getName(internalEvent))) == event_getMetaEvent(internalEvent));
	CuAssertTrue(testCase, event_getMetaEvent(eventTree_getEvent(eventTree2, event_getName(leafEvent1))) == event_getMetaEvent(leafEvent1));
	CuAssertTrue(testCase, event_getMetaEvent(eventTree_getEvent(eventTree2, event_getName(leafEvent2))) == event_getMetaEvent(leafEvent2));
	cactusEventTreeTestTeardown();
}

void testEventTree_getRootEvent(CuTest* testCase) {
	cactusEventTreeTestSetup();
	CuAssertTrue(testCase, eventTree_getRootEvent(eventTree) == rootEvent);
	cactusEventTreeTestTeardown();
}

void testEventTree_getEvent(CuTest* testCase) {
	cactusEventTreeTestSetup();
	CuAssertTrue(testCase, eventTree_getEvent(eventTree, metaEvent_getName(rootMetaEvent)) == rootEvent);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree, metaEvent_getName(internalMetaEvent)) == internalEvent);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree, metaEvent_getName(leafMetaEvent1)) == leafEvent1);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree, metaEvent_getName(leafMetaEvent2)) == leafEvent2);
	cactusEventTreeTestTeardown();
}

void testEventTree_getCommonAncestor(CuTest* testCase) {
	cactusEventTreeTestSetup();
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

	cactusEventTreeTestTeardown();
}

void testEventTree_getFlower(CuTest* testCase) {
	cactusEventTreeTestSetup();
	CuAssertTrue(testCase, eventTree_getFlower(eventTree) == flower);
	cactusEventTreeTestTeardown();
}

void testEventTree_getEventNumber(CuTest* testCase) {
	cactusEventTreeTestSetup();
	CuAssertIntEquals(testCase, 4, eventTree_getEventNumber(eventTree));
	cactusEventTreeTestTeardown();
}

void testEventTree_getFirst(CuTest* testCase) {
	cactusEventTreeTestSetup();
	CuAssertTrue(testCase, eventTree_getFirst(eventTree) == rootEvent);
	cactusEventTreeTestTeardown();
}

void testEventTree_iterator(CuTest* testCase) {
	cactusEventTreeTestSetup();
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

	cactusEventTreeTestTeardown();
}

void testEventTree_makeNewickString(CuTest* testCase) {
	cactusEventTreeTestSetup();
	CuAssertStrEquals(testCase, "((4:0.2,5:1.3)3:0.5)2:2.14748e+09;", eventTree_makeNewickString(eventTree));
	cactusEventTreeTestTeardown();
}

void testEventTree_addSiblingUnaryEvent(CuTest *testCase) {
	cactusEventTreeTestSetup();
	//Create two sibling flowers with the basic event tree..
	//then try adding events from on into the other.
	Group *group1 = group_construct2(flower);
	Group *group2 = group_construct2(flower);
	Flower *flower2 = flower_construct(cactusDisk);
	Flower *flower3 = flower_construct(cactusDisk);
	flower_setParentGroup(flower2, group1);
	flower_setParentGroup(flower3, group2);
	EventTree *eventTree2 = eventTree_copyConstruct(flower_getEventTree(flower), flower2, NULL);
	MetaEvent *unaryInternalMetaEvent1 = metaEvent_construct("UNARY1", cactusDisk);
	MetaEvent *unaryInternalMetaEvent2 = metaEvent_construct("UNARY2", cactusDisk);
	MetaEvent *unaryInternalMetaEvent3 = metaEvent_construct("UNARY3", cactusDisk);
	Event *parentUnaryEvent1 = event_construct2(unaryInternalMetaEvent1, 0.1, internalEvent, leafEvent1, eventTree);
	Event *parentUnaryEvent2 = event_construct2(unaryInternalMetaEvent2, 0.1, parentUnaryEvent1, leafEvent1, eventTree);
	Event *parentUnaryEvent3 = event_construct2(unaryInternalMetaEvent3, 0.1, internalEvent, leafEvent2, eventTree);
	//now event tree contains the added unary events.
	EventTree *eventTree3 = eventTree_copyConstruct(flower_getEventTree(flower), flower3, NULL);
	//add a couple of denovo events into the new event tree
	MetaEvent *unaryInternalMetaEvent4 = metaEvent_construct("UNARY4", cactusDisk);
	MetaEvent *unaryInternalMetaEvent5 = metaEvent_construct("UNARY5", cactusDisk);
	Event *internalEventChild = eventTree_getEvent(eventTree3, event_getName(internalEvent));
	Event *unaryEvent1 = eventTree_getEvent(eventTree3, event_getName(parentUnaryEvent1));
	Event *unaryEvent2 = eventTree_getEvent(eventTree3, event_getName(parentUnaryEvent2));
	Event *unaryEvent3 = eventTree_getEvent(eventTree3, event_getName(parentUnaryEvent3));
	Event *unaryEvent4 = event_construct2(unaryInternalMetaEvent4, 0.1,
			internalEventChild, unaryEvent3, eventTree3);
	Event *unaryEvent5 = event_construct2(unaryInternalMetaEvent5, 0.1,
				internalEventChild, unaryEvent4, eventTree3);
	//Now event-tree 2 does not contain the unary events but event-tree 3 does..

	CuAssertTrue(testCase, eventTree_getEvent(eventTree2, event_getName(unaryEvent1)) == NULL);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree2, event_getName(unaryEvent2)) == NULL);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree2, event_getName(unaryEvent3)) == NULL);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree2, event_getName(unaryEvent4)) == NULL);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree2, event_getName(unaryEvent5)) == NULL);

	eventTree_addSiblingUnaryEvent(eventTree2, unaryEvent1);
	Event *unaryEvent12 = eventTree_getEvent(eventTree2, event_getName(unaryEvent1));
	CuAssertTrue(testCase, unaryEvent12 != NULL);
	CuAssertTrue(testCase, event_getName(event_getParent(unaryEvent12)) == event_getName(internalEvent));
	CuAssertTrue(testCase, event_getChildNumber(unaryEvent12) == 1);
	CuAssertTrue(testCase, event_getName(event_getChild(unaryEvent12, 0)) == event_getName(leafEvent1));

	eventTree_addSiblingUnaryEvent(eventTree2, unaryEvent2);
	Event *unaryEvent22 = eventTree_getEvent(eventTree2, event_getName(unaryEvent2));
	CuAssertTrue(testCase, unaryEvent22 != NULL);
	CuAssertTrue(testCase, event_getName(event_getParent(unaryEvent22)) == event_getName(unaryEvent1));
	CuAssertTrue(testCase, event_getChildNumber(unaryEvent22) == 1);
	CuAssertTrue(testCase, event_getName(event_getChild(unaryEvent22, 0)) == event_getName(leafEvent1));

	eventTree_addSiblingUnaryEvent(eventTree2, unaryEvent3);
	Event *unaryEvent32 = eventTree_getEvent(eventTree2, event_getName(unaryEvent3));
	CuAssertTrue(testCase, unaryEvent32 != NULL);
	CuAssertTrue(testCase, event_getName(event_getParent(unaryEvent32)) == event_getName(internalEvent));
	CuAssertTrue(testCase, event_getChildNumber(unaryEvent32) == 1);
	CuAssertTrue(testCase, event_getName(event_getChild(unaryEvent32, 0)) == event_getName(leafEvent2));

	eventTree_addSiblingUnaryEvent(eventTree2, unaryEvent4);
	Event *unaryEvent42 = eventTree_getEvent(eventTree2, event_getName(unaryEvent4));
	CuAssertTrue(testCase, unaryEvent42 != NULL);
	CuAssertTrue(testCase, event_getName(event_getParent(unaryEvent42)) == event_getName(internalEvent));
	CuAssertTrue(testCase, event_getChildNumber(unaryEvent42) == 1);
	CuAssertTrue(testCase, event_getName(event_getChild(unaryEvent42, 0)) == event_getName(unaryEvent3));

	eventTree_addSiblingUnaryEvent(eventTree2, unaryEvent5);
	Event *unaryEvent52 = eventTree_getEvent(eventTree2, event_getName(unaryEvent5));
	CuAssertTrue(testCase, unaryEvent52 != NULL);
	CuAssertTrue(testCase, event_getName(event_getParent(unaryEvent52)) == event_getName(internalEvent));
	CuAssertTrue(testCase, event_getChildNumber(unaryEvent52) == 1);
	CuAssertTrue(testCase, event_getName(event_getChild(unaryEvent52, 0)) == event_getName(unaryEvent4));

	//uglyf("Event-tree-1 %s\n", eventTree_makeNewickString(eventTree));
	//uglyf("Event-tree-2 %s\n", eventTree_makeNewickString(eventTree3));
	//uglyf("Event-tree-3 %s\n", eventTree_makeNewickString(eventTree2));

	cactusEventTreeTestTeardown();
}

void testEventTree_serialisation(CuTest* testCase) {
	cactusEventTreeTestSetup();
	int32_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(eventTree,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))eventTree_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	eventTree_destruct(eventTree);
	void *vA2 = vA;
	eventTree = eventTree_loadFromBinaryRepresentation(&vA2, flower);
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
	testEventTree_getFlower(testCase);
	testEventTree_getEventNumber(testCase);
	testEventTree_getFirst(testCase);
	testEventTree_makeNewickString(testCase);
	testEventTree_iterator(testCase);
	testEventTree_addSiblingUnaryEvent(testCase);
	nestedTest = 0;
	cactusEventTreeTestTeardown();
}

CuSuite* cactusEventTreeTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testEventTree_copyConstruct);
	SUITE_ADD_TEST(suite, testEventTree_getRootEvent);
	SUITE_ADD_TEST(suite, testEventTree_getEvent);
	SUITE_ADD_TEST(suite, testEventTree_getCommonAncestor);
	SUITE_ADD_TEST(suite, testEventTree_getFlower);
	SUITE_ADD_TEST(suite, testEventTree_getEventNumber);
	SUITE_ADD_TEST(suite, testEventTree_getFirst);
	SUITE_ADD_TEST(suite, testEventTree_iterator);
	SUITE_ADD_TEST(suite, testEventTree_makeNewickString);
	SUITE_ADD_TEST(suite, testEventTree_addSiblingUnaryEvent);
	SUITE_ADD_TEST(suite, testEventTree_serialisation);
	SUITE_ADD_TEST(suite, testEventTree_construct);
	return suite;
}
