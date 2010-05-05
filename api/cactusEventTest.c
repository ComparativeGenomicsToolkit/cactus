#include "cactusGlobalsPrivate.h"

static NetDisk *netDisk = NULL;
static Net *net;

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

static void cactusEventTestTeardown() {
	if(!nestedTest && netDisk != NULL) {
		netDisk_destruct(netDisk);
		testCommon_deleteTemporaryNetDisk();
		netDisk = NULL;
		eventTree = NULL;
		rootEvent = NULL;
		internalEvent = NULL;
		leafEvent1 = NULL;
		leafEvent2 = NULL;
	}
}

static void cactusEventTestSetup() {
	if(!nestedTest) {
		cactusEventTestTeardown();
		netDisk = netDisk_construct(testCommon_getTemporaryNetDisk());
		net = net_construct(netDisk);
		internalMetaEvent = metaEvent_construct("INTERNAL", netDisk);
		leafMetaEvent1 = metaEvent_construct("LEAF1", netDisk);
		leafMetaEvent2 = metaEvent_construct("LEAF2", netDisk);
		eventTree = net_getEventTree(net);
		rootEvent = eventTree_getRootEvent(eventTree);
		rootMetaEvent = event_getMetaEvent(rootEvent);
		internalEvent = event_construct(internalMetaEvent, 0.5, rootEvent, eventTree);
		leafEvent1 = event_construct(leafMetaEvent1, 0.2, internalEvent, eventTree);
		leafEvent2 = event_construct(leafMetaEvent2, 1.3, internalEvent, eventTree);
	}
}

void testEvent_construct(CuTest* testCase) {
	nestedTest = 0;
	cactusEventTestSetup();
	CuAssertTrue(testCase, eventTree != NULL);
	CuAssertTrue(testCase, rootEvent != NULL);
	CuAssertTrue(testCase, internalEvent != NULL);
	CuAssertTrue(testCase, leafEvent1 != NULL);
	CuAssertTrue(testCase, leafEvent2 != NULL);
	cactusEventTestTeardown();
}

void testEvent_getParent(CuTest* testCase) {
	cactusEventTestSetup();
	CuAssertTrue(testCase, event_getParent(rootEvent) == NULL);
	CuAssertTrue(testCase, event_getParent(internalEvent) == rootEvent);
	CuAssertTrue(testCase, event_getParent(leafEvent1) == internalEvent);
	CuAssertTrue(testCase, event_getParent(leafEvent2) == internalEvent);
	cactusEventTestTeardown();
}

void testEvent_getName(CuTest* testCase) {
	cactusEventTestSetup();
	CuAssertTrue(testCase, event_getName(rootEvent) == metaEvent_getName(rootMetaEvent));
	CuAssertTrue(testCase, event_getName(internalEvent) == metaEvent_getName(internalMetaEvent));
	CuAssertTrue(testCase, event_getName(leafEvent1) == metaEvent_getName(leafMetaEvent1));
	CuAssertTrue(testCase, event_getName(leafEvent2) == metaEvent_getName(leafMetaEvent2));
	cactusEventTestTeardown();
}

void testEvent_getMetaEvent(CuTest* testCase) {
	cactusEventTestSetup();
	CuAssertTrue(testCase, event_getMetaEvent(rootEvent) == rootMetaEvent);
	CuAssertTrue(testCase, event_getMetaEvent(internalEvent) == internalMetaEvent);
	CuAssertTrue(testCase, event_getMetaEvent(leafEvent1) == leafMetaEvent1);
	CuAssertTrue(testCase, event_getMetaEvent(leafEvent2) == leafMetaEvent2);
	cactusEventTestTeardown();
}

void testEvent_getHeader(CuTest* testCase) {
	cactusEventTestSetup();
	CuAssertTrue(testCase, event_getHeader(rootEvent) == metaEvent_getHeader(rootMetaEvent));
	CuAssertTrue(testCase, event_getHeader(internalEvent) == metaEvent_getHeader(internalMetaEvent));
	CuAssertTrue(testCase, event_getHeader(leafEvent1) == metaEvent_getHeader(leafMetaEvent1));
	CuAssertTrue(testCase, event_getHeader(leafEvent2) == metaEvent_getHeader(leafMetaEvent2));
	cactusEventTestTeardown();
}

void testEvent_getBranchLength(CuTest* testCase) {
	cactusEventTestSetup();
	CuAssertDblEquals(testCase, INT32_MAX, event_getBranchLength(rootEvent), 1.000);
	CuAssertDblEquals(testCase, 0.5, event_getBranchLength(internalEvent), 0.001);
	CuAssertDblEquals(testCase, 0.2, event_getBranchLength(leafEvent1), 0.001);
	CuAssertDblEquals(testCase, 1.3, event_getBranchLength(leafEvent2), 0.001);
	cactusEventTestTeardown();
}

void testEvent_getSubTreeBranchLength(CuTest* testCase) {
	cactusEventTestSetup();
	CuAssertDblEquals(testCase, 2.0, event_getSubTreeBranchLength(rootEvent), 0.001);
	CuAssertDblEquals(testCase, 1.5, event_getSubTreeBranchLength(internalEvent), 0.001);
	CuAssertDblEquals(testCase, 0.0, event_getSubTreeBranchLength(leafEvent1), 0.001);
	CuAssertDblEquals(testCase, 0.0, event_getSubTreeBranchLength(leafEvent2), 0.001);
	cactusEventTestTeardown();
}

void testEvent_getSubTreeEventNumber(CuTest* testCase) {
	cactusEventTestSetup();
	CuAssertTrue(testCase, event_getSubTreeEventNumber(rootEvent) == 3);
	CuAssertTrue(testCase, event_getSubTreeEventNumber(internalEvent) == 2);
	CuAssertTrue(testCase, event_getSubTreeEventNumber(leafEvent1) == 0);
	CuAssertTrue(testCase, event_getSubTreeEventNumber(leafEvent2) == 0);
	cactusEventTestTeardown();
}

void testEvent_getChildNumber(CuTest* testCase) {
	cactusEventTestSetup();
	CuAssertTrue(testCase, event_getChildNumber(rootEvent) == 1);
	CuAssertTrue(testCase, event_getChildNumber(internalEvent) == 2);
	CuAssertTrue(testCase, event_getChildNumber(leafEvent1) == 0);
	CuAssertTrue(testCase, event_getChildNumber(leafEvent2) == 0);
	cactusEventTestTeardown();
}

void testEvent_getChild(CuTest* testCase) {
	cactusEventTestSetup();
	CuAssertTrue(testCase, event_getChild(rootEvent, 0) == internalEvent);
	CuAssertTrue(testCase, event_getChild(internalEvent, 0) == leafEvent1);
	CuAssertTrue(testCase, event_getChild(internalEvent, 1) == leafEvent2);
	cactusEventTestTeardown();
}

void testEvent_getEventTree(CuTest* testCase) {
	cactusEventTestSetup();
	CuAssertTrue(testCase, event_getEventTree(rootEvent) == eventTree);
	CuAssertTrue(testCase, event_getEventTree(internalEvent) == eventTree);
	CuAssertTrue(testCase, event_getEventTree(leafEvent1) == eventTree);
	CuAssertTrue(testCase, event_getEventTree(leafEvent2) == eventTree);
	cactusEventTestTeardown();
}

void testEvent_isAncestor(CuTest* testCase) {
	cactusEventTestSetup();
	//negatives
	CuAssertTrue(testCase, !event_isAncestor(leafEvent1, leafEvent2));
	CuAssertTrue(testCase, !event_isAncestor(leafEvent2, leafEvent1));
	CuAssertTrue(testCase, !event_isAncestor(internalEvent, leafEvent1));
	CuAssertTrue(testCase, !event_isAncestor(internalEvent, leafEvent2));
	CuAssertTrue(testCase, !event_isAncestor(rootEvent, leafEvent1));
	CuAssertTrue(testCase, !event_isAncestor(rootEvent, leafEvent2));
	CuAssertTrue(testCase, !event_isAncestor(rootEvent, internalEvent));
	//selfs
	CuAssertTrue(testCase, !event_isAncestor(leafEvent1, leafEvent1));
	CuAssertTrue(testCase, !event_isAncestor(leafEvent2, leafEvent2));
	CuAssertTrue(testCase, !event_isAncestor(internalEvent, internalEvent));
	CuAssertTrue(testCase, !event_isAncestor(rootEvent, rootEvent));
	//positives
	CuAssertTrue(testCase, event_isAncestor(leafEvent1, internalEvent));
	CuAssertTrue(testCase, event_isAncestor(leafEvent2, internalEvent));
	CuAssertTrue(testCase, event_isAncestor(leafEvent1, rootEvent));
	CuAssertTrue(testCase, event_isAncestor(leafEvent2, rootEvent));
	CuAssertTrue(testCase, event_isAncestor(internalEvent, rootEvent));

	cactusEventTestTeardown();
}

void testEvent_isDescendant(CuTest* testCase) {
	cactusEventTestSetup();
	//negatives
	CuAssertTrue(testCase, !event_isDescendant(leafEvent1, leafEvent2));
	CuAssertTrue(testCase, !event_isDescendant(leafEvent2, leafEvent1));
	CuAssertTrue(testCase, !event_isDescendant(leafEvent1, rootEvent));
	CuAssertTrue(testCase, !event_isDescendant(leafEvent2, rootEvent));
	CuAssertTrue(testCase, !event_isDescendant(leafEvent1, internalEvent));
	CuAssertTrue(testCase, !event_isDescendant(leafEvent2, internalEvent));
	CuAssertTrue(testCase, !event_isDescendant(internalEvent, rootEvent));
	//selfs
	CuAssertTrue(testCase, !event_isDescendant(leafEvent1, leafEvent1));
	CuAssertTrue(testCase, !event_isDescendant(leafEvent2, leafEvent2));
	CuAssertTrue(testCase, !event_isDescendant(internalEvent, internalEvent));
	CuAssertTrue(testCase, !event_isDescendant(rootEvent, rootEvent));
	//positives
	CuAssertTrue(testCase, event_isDescendant(internalEvent, leafEvent1));
	CuAssertTrue(testCase, event_isDescendant(internalEvent, leafEvent2));
	CuAssertTrue(testCase, event_isDescendant(rootEvent, internalEvent));
	CuAssertTrue(testCase, event_isDescendant(rootEvent, leafEvent1));
	CuAssertTrue(testCase, event_isDescendant(rootEvent, leafEvent2));

	cactusEventTestTeardown();
}

void testEvent_isSibling(CuTest* testCase) {
	cactusEventTestSetup();
	//negatives
	CuAssertTrue(testCase, !event_isSibling(leafEvent1, rootEvent));
	CuAssertTrue(testCase, !event_isSibling(leafEvent2, rootEvent));
	CuAssertTrue(testCase, !event_isSibling(leafEvent1, internalEvent));
	CuAssertTrue(testCase, !event_isSibling(leafEvent2, internalEvent));
	CuAssertTrue(testCase, !event_isSibling(internalEvent, rootEvent));
	CuAssertTrue(testCase, !event_isSibling(internalEvent, leafEvent1));
	CuAssertTrue(testCase, !event_isSibling(internalEvent, leafEvent2));
	CuAssertTrue(testCase, !event_isSibling(rootEvent, internalEvent));
	CuAssertTrue(testCase, !event_isSibling(rootEvent, leafEvent1));
	CuAssertTrue(testCase, !event_isSibling(rootEvent, leafEvent2));
	//selfs
	CuAssertTrue(testCase, !event_isSibling(leafEvent1, leafEvent1));
	CuAssertTrue(testCase, !event_isSibling(leafEvent2, leafEvent2));
	CuAssertTrue(testCase, !event_isSibling(internalEvent, internalEvent));
	CuAssertTrue(testCase, !event_isSibling(rootEvent, rootEvent));

	//positives
	CuAssertTrue(testCase, event_isSibling(leafEvent1, leafEvent2));
	CuAssertTrue(testCase, event_isSibling(leafEvent2, leafEvent1));

	cactusEventTestTeardown();
}

void testEvent_serialisation(CuTest* testCase) {
	cactusEventTestSetup();
	int32_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(leafEvent1,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))event_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	event_destruct(leafEvent1);
	void *vA2 = vA;
	leafEvent1 = event_loadFromBinaryRepresentation(&vA2, eventTree);
	free(vA);
	nestedTest = 1;
	testEvent_getParent(testCase);
	testEvent_getName(testCase);
	testEvent_getMetaEvent(testCase);
	testEvent_getHeader(testCase);
	testEvent_getBranchLength(testCase);
	testEvent_getSubTreeBranchLength(testCase);
	testEvent_getSubTreeEventNumber(testCase);
	testEvent_getChildNumber(testCase);
	//testEvent_getChild(testCase); -- won't work as doesn't preserve order of leaves, which is okay -- here's a replacement.
	CuAssertTrue(testCase, event_getChild(rootEvent, 0) == internalEvent);
	CuAssertTrue(testCase, event_getChild(internalEvent, 0) == leafEvent2);
	CuAssertTrue(testCase, event_getChild(internalEvent, 1) == leafEvent1);

	testEvent_getEventTree(testCase);
	testEvent_isAncestor(testCase);
	testEvent_isDescendant(testCase);
	testEvent_isSibling(testCase);
	nestedTest = 0;
	cactusEventTestTeardown();
}

CuSuite* cactusEventTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testEvent_getParent);
	SUITE_ADD_TEST(suite, testEvent_getName);
	SUITE_ADD_TEST(suite, testEvent_getMetaEvent);
	SUITE_ADD_TEST(suite, testEvent_getHeader);
	SUITE_ADD_TEST(suite, testEvent_getBranchLength);
	SUITE_ADD_TEST(suite, testEvent_getSubTreeBranchLength);
	SUITE_ADD_TEST(suite, testEvent_getSubTreeEventNumber);
	SUITE_ADD_TEST(suite, testEvent_getChildNumber);
	SUITE_ADD_TEST(suite, testEvent_getChild);
	SUITE_ADD_TEST(suite, testEvent_getEventTree);
	SUITE_ADD_TEST(suite, testEvent_isAncestor);
	SUITE_ADD_TEST(suite, testEvent_isDescendant);
	SUITE_ADD_TEST(suite, testEvent_isSibling);
	SUITE_ADD_TEST(suite, testEvent_serialisation);
	SUITE_ADD_TEST(suite, testEvent_construct); //put this last to ensure cleanup.
	return suite;
}
