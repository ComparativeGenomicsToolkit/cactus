#include "cactusGlobalsPrivate.h"

static NetDisk *netDisk = NULL;
static Net *net;
static MetaEvent *metaEvent;
static EventTree *eventTree;
static MetaSequence *metaSequence;
static Sequence *sequence;
static const char *sequenceString = "ACTGGCACTG";
static const char *headerString = ">one";

static bool nestedTest = 0;

void cactusSequenceTestTeardown() {
	if(!nestedTest && netDisk != NULL) {
		netDisk_destruct(netDisk);
		testCommon_deleteTemporaryNetDisk();
		netDisk = NULL;
		net = NULL;
		metaEvent = NULL;
		eventTree = NULL;
		metaSequence = NULL;
		sequence = NULL;
	}
}

void cactusSequenceTestSetup() {
	if(!nestedTest) {
		cactusSequenceTestTeardown();
		netDisk = netDisk_construct(testCommon_getTemporaryNetDisk());
		net = net_construct(netDisk);
		metaEvent = metaEvent_construct("ROOT", netDisk);
		eventTree = eventTree_construct(metaEvent, net);
		metaSequence = metaSequence_construct(1, 10, sequenceString,
						   headerString, metaEvent_getName(metaEvent), netDisk);
		sequence = sequence_construct(metaSequence, net);
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
	CuAssertTrue(testCase, net_getSequence(net, sequence_getName(sequence)) == sequence);
	cactusSequenceTestTeardown();
}

void testSequence_getEvent(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertTrue(testCase, sequence_getEvent(sequence) == eventTree_getEvent(eventTree, metaEvent_getName(metaEvent)));
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

void testSequence_getHeader(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertStrEquals(testCase, headerString, sequence_getHeader(sequence));
	cactusSequenceTestTeardown();
}

void testSequence_getNet(CuTest* testCase) {
	cactusSequenceTestSetup();
	CuAssertTrue(testCase, sequence_getNet(sequence) == net);
	cactusSequenceTestTeardown();
}

void testSequence_serialisation(CuTest* testCase) {
	cactusSequenceTestSetup();
	int32_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(sequence,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))sequence_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	sequence_destruct(sequence);
	void *vA2 = vA;
	sequence = sequence_loadFromBinaryRepresentation(&vA2, net);
	nestedTest = 1;
	testSequence_getMetaSequence(testCase);
	testSequence_getStart(testCase);
	testSequence_getLength(testCase);
	testSequence_getName(testCase);
	testSequence_getEvent(testCase);
	testSequence_getString(testCase);
	testSequence_getHeader(testCase);
	testSequence_getNet(testCase);
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
	SUITE_ADD_TEST(suite, testSequence_getHeader);
	SUITE_ADD_TEST(suite, testSequence_getNet);
	SUITE_ADD_TEST(suite, testSequence_serialisation);
	SUITE_ADD_TEST(suite, testSequence_construct);
	return suite;
}
