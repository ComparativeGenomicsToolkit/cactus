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
static Segment *segment;
static Segment *segment2;
static Cap *cap;
static Cap *cap2;
static Reference *reference;

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
		reference = NULL;
	}
}

static void cactusNetTestSetup() {
	cactusNetTestTeardown();
	netDisk = netDisk_construct(testCommon_getTemporaryNetDisk());
	net = net_construct(netDisk);
	metaEvent = metaEvent_construct("ROOT", netDisk);
	eventTree = eventTree_construct(metaEvent, net);
	assert(net_getReference(net) == NULL);
	reference = reference_construct(net);
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

static void capsSetup() {
	endsSetup();
	cap = cap_construct(end, eventTree_getRootEvent(eventTree));
	cap2 = cap_construct(end2, eventTree_getRootEvent(eventTree));
}

static void blocksSetup() {
	block = block_construct(1, net);
	block2 = block_construct(2, net);
}

static void segmentsSetup() {
	blocksSetup();
	segment = segment_construct(block, eventTree_getRootEvent(eventTree));
	segment2 = segment_construct(block2, eventTree_getRootEvent(eventTree));
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

void testNet_getReference(CuTest* testCase) {
	cactusNetTestSetup();
	CuAssertTrue(testCase, net_getReference(net) == reference);
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

void testNet_cap(CuTest* testCase) {
	testObjectRetrieval(testCase, capsSetup,
			(int32_t (*)(Net *net))net_getCapNumber,
			(void *(*)(Net *))net_getFirstCap,
			(Name (*)(void *))cap_getName,
			(void *(*)(Net *, Name name))net_getCap,
			(void *(*)(Net *))net_getCapIterator,
			(void (*)(void *))net_destructCapIterator,
			(void *(*)(void *))net_getNextCap,
			(void *(*)(void *))net_getPreviousCap,
			(void *(*)(void *))net_copyCapIterator,
			(void **)(&cap), (void **)(&cap2));
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

void testNet_segment(CuTest* testCase) {
	testObjectRetrieval(testCase, segmentsSetup,
			(int32_t (*)(Net *net))net_getSegmentNumber,
			(void *(*)(Net *))net_getFirstSegment,
			(Name (*)(void *))segment_getName,
			(void *(*)(Net *, Name name))net_getSegment,
			(void *(*)(Net *))net_getSegmentIterator,
			(void (*)(void *))net_destructSegmentIterator,
			(void *(*)(void *))net_getNextSegment,
			(void *(*)(void *))net_getPreviousSegment,
			(void *(*)(void *))net_copySegmentIterator,
			(void **)(&segment), (void **)(&segment2));
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

void testNet_getEndNumber(CuTest *testCase) {
	/*
	 * Tests the different end number functions.
	 */
	cactusNetTestSetup();
	CuAssertTrue(testCase, net_getEndNumber(net) == 0);
	CuAssertTrue(testCase, net_getBlockEndNumber(net) == 0);
	CuAssertTrue(testCase, net_getStubEndNumber(net) == 0);
	CuAssertTrue(testCase, net_getFreeStubEndNumber(net) == 0);
	CuAssertTrue(testCase, net_getAttachedStubEndNumber(net) == 0);
	int32_t blockNumber = 10;
	int32_t freeStubEndNumber = 5;
	int32_t attachedStubEndNumber = 3;
	int32_t i;
	for(i=0; i<blockNumber; i++) {
		block_construct(1, net);
	}
	for(i=0; i<freeStubEndNumber; i++) {
		end_construct(0, net);
	}
	for(i=0; i<attachedStubEndNumber; i++) {
		end_construct(1, net);
	}

	CuAssertTrue(testCase, net_getEndNumber(net) == blockNumber*2 + freeStubEndNumber + attachedStubEndNumber);
	CuAssertTrue(testCase, net_getBlockEndNumber(net) == blockNumber*2);
	CuAssertTrue(testCase, net_getStubEndNumber(net) == freeStubEndNumber + attachedStubEndNumber);
	CuAssertTrue(testCase, net_getFreeStubEndNumber(net) == freeStubEndNumber);
	CuAssertTrue(testCase, net_getAttachedStubEndNumber(net) == attachedStubEndNumber);
	cactusNetTestTeardown();
}


/*void testNet_mergeNets(CuTest *testCase) {
	cactusNetTestSetup();
	//construct the nets to merge...
	Net *net1 = net_construct(netDisk);
	Net *net2 = net_construct(netDisk);

	//construct nets that are children of the nets to merge.
	Net *net3 = net_construct(netDisk);
	Net *net4 = net_construct(netDisk);
	Net *net5 = net_construct(netDisk);

	//Make event tree
	MetaEvent *internalMetaEvent = metaEvent_construct("INTERNAL", netDisk);
	MetaEvent *leafMetaEvent1 = metaEvent_construct("LEAF1", netDisk);
	MetaEvent *leafMetaEvent2 = metaEvent_construct("LEAF2", netDisk);
	Event *internalEvent = event_construct(internalMetaEvent, 0.5, eventTree_getRootEvent(eventTree), eventTree);
	Event *leafEvent1 = event_construct(leafMetaEvent1, 0.2, internalEvent, eventTree);
	event_construct(leafMetaEvent2, 1.3, internalEvent, eventTree);

	//Copy the event tree into the children (this is tested by the merge event function)..
	EventTree *eventTree1 = eventTree_copyConstruct(eventTree, net1, NULL);
	eventTree_copyConstruct(eventTree, net2, NULL);
	MetaEvent *unaryInternalMetaEvent1 = metaEvent_construct("UNARY1", netDisk);
	MetaEvent *unaryInternalMetaEvent2 = metaEvent_construct("UNARY2", netDisk);
	MetaEvent *unaryInternalMetaEvent3 = metaEvent_construct("UNARY3", netDisk);

	Event *internalEventChild = eventTree_getEvent(eventTree1, event_getName(internalEvent));
	Event *unaryEvent1 = event_construct2(unaryInternalMetaEvent1, 0.1,
				internalEventChild, eventTree_getEvent(eventTree1, event_getName(leafEvent1)), eventTree1);
	Event *unaryEvent2 = event_construct2(unaryInternalMetaEvent2, 0.1,
				internalEventChild, unaryEvent1, eventTree1);
	event_construct2(unaryInternalMetaEvent3, 0.1,
				internalEventChild, unaryEvent2, eventTree1);
	CuAssertTrue(testCase, eventTree_getEventNumber(eventTree1) == 7);

	//Make some sequences
	MetaSequence *metaSequence1 = metaSequence_construct(0, 5, "ACTGG", "one", event_getName(unaryEvent1), netDisk);
	MetaSequence *metaSequence2 = metaSequence_construct(0, 5, "CCCCC", "two", event_getName(unaryEvent2), netDisk);
	MetaSequence *metaSequence3 = metaSequence_construct(0, 5, "TTTTT", "three", event_getName(leafEvent1), netDisk);
	Sequence *sequence1 = sequence_construct(metaSequence1, net1);
	Sequence *sequence2 = sequence_construct(metaSequence2, net1);
	Sequence *sequence3 = sequence_construct(metaSequence3, net2);

	//Make children and parent nets and some groups..
	Group *group3 = group_construct(net1, net3);
	Group *group4 = group_construct(net2, net4);
	Group *group5 = group_construct(net2, net5);
	Group *group6 = group_construct2(net1);
	Group *group7 = group_construct2(net2);

	//Make some stubs to put in the ends..
	End *end1 = end_construct(0, net1);
	End *end2 = end_construct(0, net1);
	End *end3 = end_construct(0, net2);
	end_setGroup(end1, group6);
	end_setGroup(end2, group6);
	end_setGroup(end3, group7);

	//Make a few caps
	Cap *cap1 = cap_construct2(end1, 0, 1, 1, sequence1);
	Cap *cap2 = cap_construct2(end1, 0, 1, 1, sequence2);
	cap_makeParentAndChild(cap2, cap1);
	Cap *cap3 = cap_construct2(end3, 0, 1, 1, sequence3);

	//Make some blocks..
	Block *block1 = block_construct(0, net1);
	Block *block2 = block_construct(0, net2);
	Block *block3 = block_construct(0, net2);
	end_setGroup(block_get5End(block3), group7);

	//Make a segment
	Segment *segment1 = segment_construct(block1, eventTree_getEvent(eventTree1, metaEvent_getName(leafMetaEvent1)));

	//Make some chains...
	Chain *chain1 = chain_construct(net1);
	Chain *chain2 = chain_construct(net2);

	//Make a couple of links in the chain
	Link *link1 = link_construct(end1, end2, group6, chain1);
	Link *link2 = link_construct(end3, block_get5End(block3), group7, chain2);

	Name netName1 = net_getName(net1);
	Net *net6 = net_mergeNets(net1, net2);

	//Check the events
	EventTree *eventTree3 = net_getEventTree(net6);
	CuAssertTrue(testCase, eventTree_getEventNumber(eventTree3) == 7);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent1)) != NULL);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent2)) != NULL);
	CuAssertTrue(testCase, eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent3)) != NULL);

	CuAssertTrue(testCase, event_getName(event_getParent(eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent1)))) == metaEvent_getName(unaryInternalMetaEvent2));
	CuAssertTrue(testCase, event_getName(event_getParent(eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent2)))) == metaEvent_getName(unaryInternalMetaEvent3));
	CuAssertTrue(testCase, event_getName(event_getParent(eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent3)))) == metaEvent_getName(internalMetaEvent));
	CuAssertTrue(testCase, event_getName(event_getParent(eventTree_getEvent(eventTree3, metaEvent_getName(leafMetaEvent1)))) == metaEvent_getName(unaryInternalMetaEvent1));

	//Check the sequences
	CuAssertTrue(testCase, net_getSequence(net6, metaSequence_getName(metaSequence2)) != NULL);
	CuAssertTrue(testCase, net_getSequence(net6, metaSequence_getName(metaSequence3)) != NULL);
	CuAssertTrue(testCase, net_getSequence(net6, metaSequence_getName(metaSequence1)) != NULL);

	CuAssertTrue(testCase, sequence_getEvent(net_getSequence(net6, metaSequence_getName(metaSequence1))) == eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent1)));
	CuAssertTrue(testCase, sequence_getEvent(net_getSequence(net6, metaSequence_getName(metaSequence2))) == eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent2)));
	CuAssertTrue(testCase, sequence_getEvent(net_getSequence(net6, metaSequence_getName(metaSequence3))) == eventTree_getEvent(eventTree3, metaEvent_getName(leafMetaEvent1)));

	//Check the groups..
	CuAssertTrue(testCase, net_getGroupNumber(net6) == 5);
	CuAssertTrue(testCase, net_getGroup(net6, group_getName(group3)) != NULL);
	CuAssertTrue(testCase, net_getGroup(net6, group_getName(group4)) != NULL);
	CuAssertTrue(testCase, net_getGroup(net6, group_getName(group5)) != NULL);
	CuAssertTrue(testCase, net_getGroup(net6, group_getName(group6)) != NULL);
	CuAssertTrue(testCase, net_getGroup(net6, group_getName(group7)) != NULL);
	CuAssertTrue(testCase, net_getGroup(net6, group_getName(group3)) == net_getParentGroup(net3));
	CuAssertTrue(testCase, net_getGroup(net6, group_getName(group4)) == net_getParentGroup(net4));
	CuAssertTrue(testCase, net_getGroup(net6, group_getName(group5)) == net_getParentGroup(net5));

	//Check the ends
	CuAssertTrue(testCase, net_getEndNumber(net6) == 9);
	CuAssertTrue(testCase, net_getEnd(net6, end_getName(end1)) == end1);
	CuAssertTrue(testCase, net_getEnd(net6, end_getName(end2)) == end2);
	CuAssertTrue(testCase, net_getEnd(net6, end_getName(end3)) == end3);
	CuAssertTrue(testCase, net_getEnd(net6, end_getName(block_get5End(block1))) == block_get5End(block1));
	CuAssertTrue(testCase, net_getEnd(net6, end_getName(block_get5End(block2))) == block_get5End(block2));
	CuAssertTrue(testCase, net_getEnd(net6, end_getName(block_get5End(block3))) == block_get5End(block3));
	CuAssertTrue(testCase, net_getEnd(net6, end_getName(block_get3End(block1))) == block_get3End(block1));
	CuAssertTrue(testCase, net_getEnd(net6, end_getName(block_get3End(block2))) == block_get3End(block2));
	CuAssertTrue(testCase, net_getEnd(net6, end_getName(block_get3End(block3))) == block_get3End(block3));

	//Check we have the right index of caps.
	CuAssertTrue(testCase, net_getCapNumber(net6) == 5);
	CuAssertTrue(testCase, net_getCap(net6, cap_getName(cap1)) == cap1);
	CuAssertTrue(testCase, net_getCap(net6, cap_getName(cap2)) == cap2);
	CuAssertTrue(testCase, net_getCap(net6, cap_getName(cap3)) == cap3);

	//Check the the events have been reassigned for the caps
	CuAssertTrue(testCase, cap_getEvent(cap1) == eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent1)));
	CuAssertTrue(testCase, cap_getEvent(cap2) == eventTree_getEvent(eventTree3, metaEvent_getName(unaryInternalMetaEvent2)));
	CuAssertTrue(testCase, cap_getEvent(cap3) == eventTree_getEvent(eventTree3, metaEvent_getName(leafMetaEvent1)));

	//Check the caps have the right sequences.
	CuAssertTrue(testCase, cap_getSequence(cap1) == net_getSequence(net6, metaSequence_getName(metaSequence1)));
	CuAssertTrue(testCase, cap_getSequence(cap2) == net_getSequence(net6, metaSequence_getName(metaSequence2)));
	CuAssertTrue(testCase, cap_getSequence(cap3) == net_getSequence(net6, metaSequence_getName(metaSequence3)));

	//Check the blocks
	CuAssertTrue(testCase, net_getBlockNumber(net6) == 3);
	CuAssertTrue(testCase, net_getBlock(net6, block_getName(block1)) == block1);
	CuAssertTrue(testCase, net_getBlock(net6, block_getName(block2)) == block2);
	CuAssertTrue(testCase, net_getBlock(net6, block_getName(block3)) == block3);

	//Check the segments
	CuAssertTrue(testCase, net_getSegmentNumber(net6) == 1);
	CuAssertTrue(testCase, net_getSegment(net6, segment_getName(segment1)) == segment1);
	CuAssertTrue(testCase, segment_getSequence(segment1) == NULL);
	CuAssertTrue(testCase, segment_getEvent(segment1) == eventTree_getEvent(eventTree3, event_getName(leafEvent1)));

	//Check the chains
	CuAssertTrue(testCase, net_getChainNumber(net6) == 2);
	CuAssertTrue(testCase, net_getChain(net6, chain_getName(chain1)) == chain1);
	CuAssertTrue(testCase, net_getChain(net6, chain_getName(chain2)) == chain2);
	CuAssertTrue(testCase, link_getChain(link1) == chain1);
	CuAssertTrue(testCase, link_getChain(link2) == chain2);
	CuAssertTrue(testCase, link_get5End(link1) == end1);
	CuAssertTrue(testCase, link_get3End(link1) == end2);
	CuAssertTrue(testCase, link_get5End(link2) == end3);
	CuAssertTrue(testCase, link_get3End(link2) == block_get5End(block3));

	//Check net1 is no longer in the database anywhere...
	CuAssertTrue(testCase, netDisk_getNet(netDisk, netName1) == NULL);

	//Check the merged net is in the database.
	CuAssertTrue(testCase, netDisk_getNet(netDisk, net_getName(net6)) == net6);

	cactusNetTestTeardown();
}*/

void testNet_builtBlocks(CuTest *testCase) {
	cactusNetTestSetup();

	CuAssertTrue(testCase, !net_builtBlocks(net));
	net_setBuiltBlocks(net, 0);
	CuAssertTrue(testCase, !net_builtBlocks(net));
	net_setBuiltBlocks(net, 1);
	CuAssertTrue(testCase, net_builtBlocks(net));
	net_setBuiltBlocks(net, 0);
	CuAssertTrue(testCase, !net_builtBlocks(net));

	cactusNetTestTeardown();
}

void testNet_builtTrees(CuTest *testCase) {
	cactusNetTestSetup();

	CuAssertTrue(testCase, !net_builtTrees(net));
	net_setBuiltTrees(net, 0);
	CuAssertTrue(testCase, !net_builtTrees(net));
	net_setBuiltTrees(net, 1);
	CuAssertTrue(testCase, net_builtTrees(net));
	net_setBuiltTrees(net, 0);
	CuAssertTrue(testCase, !net_builtTrees(net));

	cactusNetTestTeardown();
}

void testNet_builtFaces(CuTest *testCase) {
	cactusNetTestSetup();

	CuAssertTrue(testCase, !net_builtFaces(net));
	net_setBuildFaces(net, 0);
	CuAssertTrue(testCase, !net_builtFaces(net));
	net_setBuildFaces(net, 1);
	CuAssertTrue(testCase, net_builtFaces(net));
	net_setBuildFaces(net, 0);
	CuAssertTrue(testCase, !net_builtFaces(net));

	cactusNetTestTeardown();
}

void testNet_isLeaf(CuTest *testCase) {
	cactusNetTestSetup();
	CuAssertTrue(testCase, net_isLeaf(net));
	Group *group = group_construct2(net);
	CuAssertTrue(testCase, net_isLeaf(net));
	group_makeNestedNet(group);
	CuAssertTrue(testCase, !net_isLeaf(net));
	cactusNetTestTeardown();
}

void testNet_isTerminal(CuTest *testCase) {
	cactusNetTestSetup();
	CuAssertTrue(testCase, net_isTerminal(net));
	group_construct2(net);
	CuAssertTrue(testCase, net_isTerminal(net));
	end_construct(0, net);
	CuAssertTrue(testCase, net_isTerminal(net));
	block_construct(1, net);
	CuAssertTrue(testCase, !net_isTerminal(net));
	cactusNetTestTeardown();
}

CuSuite* cactusNetTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testNet_getName);
	SUITE_ADD_TEST(suite, testNet_getNetDisk);
	SUITE_ADD_TEST(suite, testNet_getEventTree);
	SUITE_ADD_TEST(suite, testNet_sequence);
	SUITE_ADD_TEST(suite, testNet_cap);
	SUITE_ADD_TEST(suite, testNet_end);
	SUITE_ADD_TEST(suite, testNet_getEndNumber);
	SUITE_ADD_TEST(suite, testNet_segment);
	SUITE_ADD_TEST(suite, testNet_block);
	SUITE_ADD_TEST(suite, testNet_group);
	SUITE_ADD_TEST(suite, testNet_chain);
	SUITE_ADD_TEST(suite, testNet_face);
	SUITE_ADD_TEST(suite, testNet_getReference);
	//SUITE_ADD_TEST(suite, testNet_mergeNets);
	SUITE_ADD_TEST(suite, testNet_builtBlocks);
	SUITE_ADD_TEST(suite, testNet_builtTrees);
	SUITE_ADD_TEST(suite, testNet_builtFaces);
	SUITE_ADD_TEST(suite, testNet_isLeaf);
	SUITE_ADD_TEST(suite, testNet_isTerminal);
	SUITE_ADD_TEST(suite, testNet_constructAndDestruct);
	return suite;
}
