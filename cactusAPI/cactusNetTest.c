#include "cactusGlobalsPrivate.h"

/*
 * Global variables for test.
 */

static NetDisk *netDisk = NULL;
static Net *net;
static MetaEvent *metaEvent;
static EventTree *eventTree;
static MetaSequence *metaSequence;
static MetaSequence *metaSequence2;
static Sequence *sequence;
static Sequence *sequence2;
static End *end;
static End *end2;
static Block *block;
static Block *block2;
static Group *group;
static Group *group2;
static Face *face;
static Face *face2;
static Chain *chain;
static Chain *chain2;

/*
 * Setup/teardown functions.
 */

static void cactusNetTestTeardown() {
	if(netDisk != NULL) {
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

static void cactusNetTestSetup() {
	cactusNetTestTeardown();
	netDisk = netDisk_construct(testCommon_getTemporaryNetDisk());
	net = net_construct(netDisk);
	metaEvent = metaEvent_construct("ROOT", netDisk);
	eventTree = eventTree_construct(metaEvent, net);
}

static void sequenceSetup() {
	metaSequence = metaSequence_construct(0, 10, "ACTGACTGAC",
			">one", metaEvent_getName(metaEvent), netDisk);
	sequence = sequence_construct(metaSequence, net);
	metaSequence2 = metaSequence_construct(0, 10, "ACTGACTGAC",
				">two", metaEvent_getName(metaEvent), netDisk);
	sequence2 = sequence_construct(metaSequence2, net);
}

static void endsSetup() {
	end = end_construct(1, net);
	end2 = end_construct(1, net);
}

static void blocksSetup() {
	block = block_construct(1, net);
	block2 = block_construct(2, net);
}

static void chainsSetup() {
	chain = chain_construct(net);
	chain2 = chain_construct(net);
}

static void groupsSetup() {
	group = group_construct(net, net_construct(netDisk));
	group2 = group_construct(net, net_construct(netDisk));
}

static void facesSetup() {
	face = face_construct(net);
	face2 = face_construct(net);
}

/*
 * This tests all the retrieval functions for each type of object.
 */
static void testObjectRetrieval(CuTest* testCase,
		void (*setupFn)(),
		int32_t (*getObjectNumberFn)(Net *net),
		void *(*getFirstObjectFn)(Net *net),
		Name (*objectGetNameFn)(void *),
		void *(*getObjectFn)(Net *net, Name name),
		void *(*constructIterator)(Net *net),
		void (*destructIterator)(void *iterator),
		void *(*getNext)(void *iterator),
		void *(*getPrevious)(void *iterator),
		void *(*copyIterator)(void *iterator),
		void **object, void **object2) {
	cactusNetTestSetup();
	/*
	 * Test number function
	 */
	CuAssertTrue(testCase, getObjectNumberFn(net) == 0);
	setupFn();
	CuAssertTrue(testCase, getObjectNumberFn(net) == 2);

	/*
	 * Test get first function.
	 */
	CuAssertTrue(testCase, getFirstObjectFn(net) == *object);

	/*
	 * Test get function
	 */
	CuAssertTrue(testCase, getObjectFn(net, objectGetNameFn(*object)) == *object);
	CuAssertTrue(testCase, getObjectFn(net, objectGetNameFn(*object2)) == *object2);

	/*
	 * Test iterator.
	 */
	void *iterator = constructIterator(net);
	CuAssertTrue(testCase, getNext(iterator) == *object);
	CuAssertTrue(testCase, getNext(iterator) == *object2);
	CuAssertTrue(testCase, getNext(iterator) == NULL);
	void *iterator2 = copyIterator(iterator);
	CuAssertTrue(testCase, getPrevious(iterator) == *object2);
	CuAssertTrue(testCase, getPrevious(iterator) == *object);
	CuAssertTrue(testCase, getPrevious(iterator) == NULL);
	destructIterator(iterator);
	CuAssertTrue(testCase, getPrevious(iterator2) == *object2);
	CuAssertTrue(testCase, getPrevious(iterator2) == *object);
	CuAssertTrue(testCase, getPrevious(iterator2) == NULL);
	destructIterator(iterator2);
	cactusNetTestTeardown();
}

/*
 * Now all the actual tests.
 */

void testNet_constructAndDestruct(CuTest* testCase) {
	cactusNetTestSetup();
	CuAssertTrue(testCase, net != NULL);
	cactusNetTestTeardown();
}

void testNet_getName(CuTest* testCase) {
	cactusNetTestSetup();
	CuAssertTrue(testCase, net_getName(net) != NULL_NAME);
	CuAssertTrue(testCase, netDisk_getNet(netDisk, net_getName(net)) == net);
	cactusNetTestTeardown();
}

void testNet_getNetDisk(CuTest* testCase) {
	cactusNetTestSetup();
	CuAssertTrue(testCase, net_getNetDisk(net) == netDisk);
	cactusNetTestTeardown();
}

void testNet_getEventTree(CuTest* testCase) {
	cactusNetTestSetup();
	CuAssertTrue(testCase, net_getEventTree(net) == eventTree);
	cactusNetTestTeardown();
}

void testNet_sequence(CuTest* testCase) {
	testObjectRetrieval(testCase, sequenceSetup,
			(int32_t (*)(Net *net))net_getSequenceNumber,
			(void *(*)(Net *))net_getFirstSequence,
			(Name (*)(void *))sequence_getName,
			(void *(*)(Net *, Name name))net_getSequence,
			(void *(*)(Net *))net_getSequenceIterator,
			(void (*)(void *))net_destructSequenceIterator,
			(void *(*)(void *))net_getNextSequence,
			(void *(*)(void *))net_getPreviousSequence,
			(void *(*)(void *))net_copySequenceIterator,
			(void **)(&sequence), (void **)(&sequence2));
}

void testNet_end(CuTest* testCase) {
	testObjectRetrieval(testCase, endsSetup,
			(int32_t (*)(Net *net))net_getEndNumber,
			(void *(*)(Net *))net_getFirstEnd,
			(Name (*)(void *))end_getName,
			(void *(*)(Net *, Name name))net_getEnd,
			(void *(*)(Net *))net_getEndIterator,
			(void (*)(void *))net_destructEndIterator,
			(void *(*)(void *))net_getNextEnd,
			(void *(*)(void *))net_getPreviousEnd,
			(void *(*)(void *))net_copyEndIterator,
			(void **)(&end), (void **)(&end2));
}

void testNet_block(CuTest* testCase) {
	testObjectRetrieval(testCase, blocksSetup,
			(int32_t (*)(Net *net))net_getBlockNumber,
			(void *(*)(Net *))net_getFirstBlock,
			(Name (*)(void *))block_getName,
			(void *(*)(Net *, Name name))net_getBlock,
			(void *(*)(Net *))net_getBlockIterator,
			(void (*)(void *))net_destructBlockIterator,
			(void *(*)(void *))net_getNextBlock,
			(void *(*)(void *))net_getPreviousBlock,
			(void *(*)(void *))net_copyBlockIterator,
			(void **)(&block), (void **)(&block2));
}

void testNet_chain(CuTest* testCase) {
	testObjectRetrieval(testCase, chainsSetup,
			(int32_t (*)(Net *net))net_getChainNumber,
			(void *(*)(Net *))net_getFirstChain,
			(Name (*)(void *))chain_getName,
			(void *(*)(Net *, Name name))net_getChain,
			(void *(*)(Net *))net_getChainIterator,
			(void (*)(void *))net_destructChainIterator,
			(void *(*)(void *))net_getNextChain,
			(void *(*)(void *))net_getPreviousChain,
			(void *(*)(void *))net_copyChainIterator,
			(void **)(&chain), (void **)(&chain2));
}

void testNet_group(CuTest* testCase) {
	testObjectRetrieval(testCase, groupsSetup,
			(int32_t (*)(Net *net))net_getGroupNumber,
			(void *(*)(Net *))net_getFirstGroup,
			(Name (*)(void *))group_getName,
			(void *(*)(Net *, Name name))net_getGroup,
			(void *(*)(Net *))net_getGroupIterator,
			(void (*)(void *))net_destructGroupIterator,
			(void *(*)(void *))net_getNextGroup,
			(void *(*)(void *))net_getPreviousGroup,
			(void *(*)(void *))net_copyGroupIterator,
			(void **)(&group), (void **)(&group2));
}

void testNet_face(CuTest* testCase) {
	testObjectRetrieval(testCase, facesSetup,
			(int32_t (*)(Net *net))net_getFaceNumber,
			(void *(*)(Net *))net_getFirstFace,
			(Name (*)(void *))face_getName,
			(void *(*)(Net *, Name name))net_getFace,
			(void *(*)(Net *))net_getFaceIterator,
			(void (*)(void *))net_destructFaceIterator,
			(void *(*)(void *))net_getNextFace,
			(void *(*)(void *))net_getPreviousFace,
			(void *(*)(void *))net_copyFaceIterator,
			(void **)(&face), (void **)(&face2));
}

CuSuite* cactusNetTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testNet_getName);
	SUITE_ADD_TEST(suite, testNet_getNetDisk);
	SUITE_ADD_TEST(suite, testNet_getEventTree);
	SUITE_ADD_TEST(suite, testNet_sequence);
	SUITE_ADD_TEST(suite, testNet_end);
	SUITE_ADD_TEST(suite, testNet_block);
	SUITE_ADD_TEST(suite, testNet_group);
	SUITE_ADD_TEST(suite, testNet_chain);
	SUITE_ADD_TEST(suite, testNet_face);
	SUITE_ADD_TEST(suite, testNet_constructAndDestruct);
	return suite;
}
