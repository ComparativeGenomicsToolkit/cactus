#include "cactusEndsTestShared.h"

static bool nestedTest = 0;

void cactusEndTestSetup() {
	if(!nestedTest) {
		cactusEndsTestSharedSetup();
	}
}

void cactusEndTestTeardown() {
	if(!nestedTest) {
		cactusEndsTestSharedTeardown();
	}
}

void testEnd_construct(CuTest* testCase) {
	nestedTest = 0;
	cactusEndTestSetup();
	CuAssertTrue(testCase, end != NULL);
	cactusEndTestTeardown();
}

int32_t testEnd_copyConstructP(Event *event) {
	assert(event != NULL);
	return 1;
}

void testEnd_copyConstruct(CuTest* testCase) {
	cactusEndTestSetup();
	Net *net2 = net_construct(netDisk);
	eventTree_copyConstruct(eventTree, net2, testEnd_copyConstructP);
	sequence_construct(metaSequence, net2);

	End *end2 = end_copyConstruct(end, net2);
	CuAssertTrue(testCase, end_getName(end2) != NULL_NAME);
	CuAssertTrue(testCase, end_getName(end2) ==  end_getName(end));
	CuAssertTrue(testCase, net_getEnd(net2, end_getName(end2)) == end2);
	CuAssertTrue(testCase, cap_getName(end_getInstance(end2, cap_getName(rootCap))) == cap_getName(rootCap));
	CuAssertTrue(testCase, cap_getName(end_getInstance(end2, cap_getName(leaf1Cap))) == cap_getName(leaf1Cap));
	CuAssertTrue(testCase, cap_getName(end_getInstance(end2, cap_getName(leaf2Cap))) == cap_getName(leaf2Cap));
	cactusEndTestTeardown();
}

void testEnd_getName(CuTest* testCase) {
	cactusEndTestSetup();
	CuAssertTrue(testCase, end_getName(end) != NULL_NAME);
	CuAssertTrue(testCase, net_getEnd(net, end_getName(end)) == end);
	CuAssertTrue(testCase, net_getEnd(net, end_getName(end_getReverse(end))) == end);
	cactusEndTestTeardown();
}

void testEnd_getOrientation(CuTest* testCase) {
	cactusEndTestSetup();
	CuAssertTrue(testCase, end_getOrientation(end));
	CuAssertTrue(testCase, !end_getOrientation(end_getReverse(end)));
	cactusEndTestTeardown();
}

void testEnd_getReverse(CuTest* testCase) {
	cactusEndTestSetup();
	CuAssertTrue(testCase, end_getReverse(end) != NULL);
	CuAssertTrue(testCase, end_getReverse(end_getReverse(end)) == end);
	cactusEndTestTeardown();
}

void testEnd_getSide(CuTest *testCase) {
	cactusEndTestSetup();
	Block *block = block_construct(10, net);
	End *_5End = block_get5End(block);
	End *_3End = block_get3End(block);

	CuAssertTrue(testCase, end_getSide(end));
	CuAssertTrue(testCase, !end_getSide(end_getReverse(end)));

	CuAssertTrue(testCase, end_getSide(_5End));
	CuAssertTrue(testCase, !end_getSide(_3End));
	CuAssertTrue(testCase, end_getSide(end_getReverse(_3End)));
	CuAssertTrue(testCase, !end_getSide(end_getReverse(_5End)));

	cactusEndTestTeardown();
}

void testEnd_getNet(CuTest* testCase) {
	cactusEndTestSetup();
	CuAssertTrue(testCase, end_getNet(end) == net);
	cactusEndTestTeardown();
}

void testEnd_getBlock(CuTest* testCase) {
	cactusEndTestSetup();
	Block *block = block_construct(10, net);
	End *leftEnd = block_get5End(block);
	End *rightEnd = block_get3End(block);

	CuAssertTrue(testCase, end_getBlock(end) == NULL);
	CuAssertTrue(testCase, end_getBlock(end_getReverse(end)) == NULL);

	CuAssertTrue(testCase, end_getBlock(leftEnd) == block);
	CuAssertTrue(testCase, end_getBlock(end_getReverse(leftEnd)) == block_getReverse(block));
	CuAssertTrue(testCase, block_getOrientation(block) == end_getOrientation(leftEnd));

	CuAssertTrue(testCase, end_getBlock(rightEnd) == block);
	CuAssertTrue(testCase, end_getBlock(end_getReverse(rightEnd)) == block_getReverse(block));
	CuAssertTrue(testCase, block_getOrientation(block) == end_getOrientation(rightEnd));

	cactusEndTestTeardown();
}

void testEnd_getOtherBlockEnd(CuTest *testCase) {
	cactusEndTestSetup();
	Block *block = block_construct(10, net);
	End *leftEnd = block_get5End(block);
	End *rightEnd = block_get3End(block);

	CuAssertTrue(testCase, end_getOtherBlockEnd(end) == NULL);
	CuAssertTrue(testCase, end_getOtherBlockEnd(end_getReverse(end)) == NULL);

	CuAssertTrue(testCase, end_getOtherBlockEnd(leftEnd) == rightEnd);
	CuAssertTrue(testCase, end_getOtherBlockEnd(rightEnd) == leftEnd);

	CuAssertTrue(testCase, end_getOtherBlockEnd(end_getReverse(leftEnd)) == end_getReverse(rightEnd));
	CuAssertTrue(testCase, end_getOtherBlockEnd(end_getReverse(rightEnd)) == end_getReverse(leftEnd));

	cactusEndTestTeardown();
}

void testEnd_getGroup(CuTest* testCase) {
	cactusEndTestSetup();
	Net *net2 = net_construct(netDisk);
	eventTree_copyConstruct(eventTree, net2, testEnd_copyConstructP);
	sequence_construct(metaSequence, net2);
	End *end2 = end_copyConstruct(end, net2);
	CuAssertTrue(testCase, end_getGroup(end) == NULL);
	Group *group = group_construct(net, net2);
	CuAssertTrue(testCase, end_getGroup(end) == group);
	CuAssertTrue(testCase, end_getGroup(end2) == NULL);
	cactusEndTestTeardown();
}

void testEnd_setGroup(CuTest* testCase) {
	cactusEndTestSetup();
	Net *net2 = net_construct(netDisk);
	Group *group2 = group_construct2(net2);
	End *end2 = end_construct(1, net2);
	End *end3 = end_construct(1, net2);
	CuAssertTrue(testCase, group_getEndNumber(group2) == 0);
	CuAssertTrue(testCase, end_getGroup(end2) == NULL);
	CuAssertTrue(testCase, end_getGroup(end3) == NULL);
	end_setGroup(end2, group2);
	CuAssertTrue(testCase, group_getEndNumber(group2) == 1);
	CuAssertTrue(testCase, end_getGroup(end2) == group2);
	CuAssertTrue(testCase, group_getEnd(group2, end_getName(end2)) == end2);
	CuAssertTrue(testCase, end_getGroup(end3) == NULL);
	end_setGroup(end3, group2);
	CuAssertTrue(testCase, group_getEndNumber(group2) == 2);
	CuAssertTrue(testCase, end_getGroup(end2) == group2);
	CuAssertTrue(testCase, group_getEnd(group2, end_getName(end2)) == end2);
	CuAssertTrue(testCase, end_getGroup(end3) == group2);
	CuAssertTrue(testCase, group_getEnd(group2, end_getName(end3)) == end3);
	end_setGroup(end3, NULL);
	end_setGroup(end2, group2);
	CuAssertTrue(testCase, group_getEndNumber(group2) == 1);
	CuAssertTrue(testCase, end_getGroup(end2) == group2);
	CuAssertTrue(testCase, group_getEnd(group2, end_getName(end2)) == end2);
	CuAssertTrue(testCase, end_getGroup(end3) == NULL);
	cactusEndTestTeardown();
}

void testEnd_getInstanceNumber(CuTest* testCase) {
	cactusEndTestSetup();
	End *end2 = end_construct(0, net);
	CuAssertTrue(testCase, end_getInstanceNumber(end2) == 0);
	CuAssertTrue(testCase, end_getInstanceNumber(end) == 4);
	cactusEndTestTeardown();
}

void testEnd_getInstance(CuTest* testCase) {
	cactusEndTestSetup();
	CuAssertTrue(testCase, end_getInstance(end, cap_getName(rootCap)) == cap_getReverse(rootCap));
	CuAssertTrue(testCase, end_getInstance(end, cap_getName(leaf1Cap)) == cap_getReverse(leaf1Cap));
	CuAssertTrue(testCase, end_getInstance(end, cap_getName(leaf2Cap)) == leaf2Cap);
	cactusEndTestTeardown();
}

void testEnd_getFirst(CuTest* testCase) {
	cactusEndTestSetup();
	CuAssertTrue(testCase, end_getFirst(end) == cap_getReverse(rootCap));
	cactusEndTestTeardown();
}

void testEnd_getSetRootInstance(CuTest* testCase) {
	cactusEndTestSetup();
	CuAssertTrue(testCase, end_getRootInstance(end) == cap_getReverse(rootCap));
	CuAssertTrue(testCase, end_getRootInstance(end_getReverse(end)) == rootCap);

	End *end2 = end_construct(0, net);
	CuAssertTrue(testCase, end_getRootInstance(end2) == NULL);
	CuAssertTrue(testCase, end_getRootInstance(end_getReverse(end2)) == NULL);
	cactusEndTestTeardown();
}

void testEnd_instanceIterator(CuTest* testCase) {
	cactusEndTestSetup();
	End_InstanceIterator *iterator = end_getInstanceIterator(end);
	CuAssertTrue(testCase, iterator != NULL);
	CuAssertTrue(testCase, end_getNext(iterator) == cap_getReverse(rootCap));
	CuAssertTrue(testCase, end_getNext(iterator) == cap_getReverse(leaf1Cap));

	End_InstanceIterator *iterator2 = end_copyInstanceIterator(iterator);

	CuAssertTrue(testCase, end_getNext(iterator) == leaf2Cap);
	CuAssertTrue(testCase, end_getNext(iterator) == leaf3Cap);
	CuAssertTrue(testCase, end_getNext(iterator) == NULL);
	CuAssertTrue(testCase, end_getPrevious(iterator) == leaf3Cap);
	CuAssertTrue(testCase, end_getPrevious(iterator) == leaf2Cap);
	CuAssertTrue(testCase, end_getPrevious(iterator) == cap_getReverse(leaf1Cap));
	CuAssertTrue(testCase, end_getPrevious(iterator) == cap_getReverse(rootCap));
	CuAssertTrue(testCase, end_getPrevious(iterator) == NULL);

	CuAssertTrue(testCase, end_getNext(iterator2) == leaf2Cap);
	CuAssertTrue(testCase, end_getNext(iterator2) == leaf3Cap);
	CuAssertTrue(testCase, end_getNext(iterator2) == NULL);
	CuAssertTrue(testCase, end_getPrevious(iterator2) == leaf3Cap);
	CuAssertTrue(testCase, end_getPrevious(iterator2) == leaf2Cap);
	CuAssertTrue(testCase, end_getPrevious(iterator2) == cap_getReverse(leaf1Cap));
	CuAssertTrue(testCase, end_getPrevious(iterator2) == cap_getReverse(rootCap));
	CuAssertTrue(testCase, end_getPrevious(iterator2) == NULL);

	end_destructInstanceIterator(iterator);
	end_destructInstanceIterator(iterator2);

	iterator = end_getInstanceIterator(end_getReverse(end));
	CuAssertTrue(testCase, end_getNext(iterator) == rootCap);
	CuAssertTrue(testCase, end_getNext(iterator) == leaf1Cap);
	CuAssertTrue(testCase, end_getNext(iterator) == cap_getReverse(leaf2Cap));
	CuAssertTrue(testCase, end_getNext(iterator) == cap_getReverse(leaf3Cap));
	CuAssertTrue(testCase, end_getNext(iterator) == NULL);
	CuAssertTrue(testCase, end_getPrevious(iterator) == cap_getReverse(leaf3Cap));
	CuAssertTrue(testCase, end_getPrevious(iterator) == cap_getReverse(leaf2Cap));
	CuAssertTrue(testCase, end_getPrevious(iterator) == leaf1Cap);
	CuAssertTrue(testCase, end_getPrevious(iterator) == rootCap);
	CuAssertTrue(testCase, end_getPrevious(iterator) == NULL);

	end_destructInstanceIterator(iterator);

	cactusEndTestTeardown();
}

void testEnd_isBlockOrStubEnd(CuTest* testCase) {
	cactusEndTestSetup();
	CuAssertTrue(testCase, end_isStubEnd(end));
	CuAssertTrue(testCase, !end_isBlockEnd(end));
	Block *block = block_construct(2, net);
	End *leftEnd = block_get5End(block);
	End *rightEnd = block_get3End(block);
	CuAssertTrue(testCase, end_isBlockEnd(leftEnd));
	CuAssertTrue(testCase, end_isBlockEnd(rightEnd));
	CuAssertTrue(testCase, !end_isStubEnd(leftEnd));
	CuAssertTrue(testCase, !end_isStubEnd(rightEnd));
	cactusEndTestTeardown();
}

void testEnd_isAttachedOrFree(CuTest* testCase) {
	cactusEndTestSetup();
	End *end2 = end_construct(1, net);
	End *end3 = end_construct(0, net);
	Block *block = block_construct(2, net);
	End *end4 = block_get5End(block);
	End *end5 = block_get3End(block);
	CuAssertTrue(testCase, end_isAttached(end2));
	CuAssertTrue(testCase, !end_isAttached(end3));
	CuAssertTrue(testCase, !end_isFree(end2));
	CuAssertTrue(testCase, end_isFree(end3));

	CuAssertTrue(testCase, !end_isAttached(end4));
	CuAssertTrue(testCase, !end_isAttached(end5));
	CuAssertTrue(testCase, end_isFree(end4));
	CuAssertTrue(testCase, end_isFree(end5));
	cactusEndTestTeardown();
}

void testEnd_serialisation(CuTest* testCase) {
	cactusEndTestSetup();
	Name rootInstanceName = cap_getName(rootCap);
	Name leaf1InstanceName = cap_getName(leaf1Cap);
	Name leaf2InstanceName = cap_getName(leaf2Cap);
	Name leaf3InstanceName = cap_getName(leaf3Cap);
	int32_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(end,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))end_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	end_destruct(end);
	void *vA2 = vA;
	end = end_loadFromBinaryRepresentation(&vA2, net);
	rootCap = cap_getReverse(end_getInstance(end, rootInstanceName));
	leaf1Cap = cap_getReverse(end_getInstance(end, leaf1InstanceName));
	leaf2Cap = end_getInstance(end, leaf2InstanceName);
	leaf3Cap = end_getInstance(end, leaf3InstanceName);
	CuAssertTrue(testCase, leaf3Cap != NULL);
	free(vA);
	nestedTest = 1;
	testEnd_copyConstruct(testCase);
	testEnd_getName(testCase);
	testEnd_getOrientation(testCase);
	testEnd_getReverse(testCase);
	testEnd_getSide(testCase);
	testEnd_getNet(testCase);
	testEnd_getBlock(testCase);
	testEnd_getOtherBlockEnd(testCase);
	testEnd_getGroup(testCase);
	testEnd_setGroup(testCase);
	testEnd_getInstanceNumber(testCase);
	testEnd_getInstance(testCase);
	testEnd_getFirst(testCase);
	testEnd_getSetRootInstance(testCase);
	testEnd_instanceIterator(testCase);
	testEnd_isBlockOrStubEnd(testCase);
	testEnd_isAttachedOrFree(testCase);
	nestedTest = 0;
	cactusEndTestTeardown();
}

CuSuite* cactusEndTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testEnd_copyConstruct);
	SUITE_ADD_TEST(suite, testEnd_getName);
	SUITE_ADD_TEST(suite, testEnd_getOrientation);
	SUITE_ADD_TEST(suite, testEnd_getReverse);
	SUITE_ADD_TEST(suite, testEnd_getSide);
	SUITE_ADD_TEST(suite, testEnd_getNet);
	SUITE_ADD_TEST(suite, testEnd_getBlock);
	SUITE_ADD_TEST(suite, testEnd_getOtherBlockEnd);
	SUITE_ADD_TEST(suite, testEnd_getGroup);
	SUITE_ADD_TEST(suite, testEnd_setGroup);
	SUITE_ADD_TEST(suite, testEnd_getInstanceNumber);
	SUITE_ADD_TEST(suite, testEnd_getInstance);
	SUITE_ADD_TEST(suite, testEnd_getFirst);
	SUITE_ADD_TEST(suite, testEnd_getSetRootInstance);
	SUITE_ADD_TEST(suite, testEnd_instanceIterator);
	SUITE_ADD_TEST(suite, testEnd_isBlockOrStubEnd);
	SUITE_ADD_TEST(suite, testEnd_isAttachedOrFree);
	SUITE_ADD_TEST(suite, testEnd_serialisation);
	SUITE_ADD_TEST(suite, testEnd_construct);
	return suite;
}
