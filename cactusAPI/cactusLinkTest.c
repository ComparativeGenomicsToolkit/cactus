#include "cactusChainsTestShared.h"

static bool nestedTest = 0;

void cactusLinkTestSetup() {
	if(!nestedTest) {
		cactusChainsSharedTestSetup();
	}
}

void cactusLinkTestTeardown() {
	if(!nestedTest) {
		cactusChainsSharedTestTeardown();
	}
}

void testLink_construct(CuTest* testCase) {
	nestedTest = 0;
	cactusLinkTestSetup();
	CuAssertTrue(testCase, link1 != NULL);
	CuAssertTrue(testCase, link2 != NULL);
	cactusLinkTestTeardown();
}

void testLink_getNextLink(CuTest* testCase) {
	cactusLinkTestSetup();
	CuAssertTrue(testCase, link_getNextLink(link1) == link2);
	CuAssertTrue(testCase, link_getNextLink(link2) == NULL);
	cactusLinkTestTeardown();
}

void testLink_getPreviousLink(CuTest* testCase) {
	cactusLinkTestSetup();
	CuAssertTrue(testCase, link_getPreviousLink(link2) == link1);
	CuAssertTrue(testCase, link_getPreviousLink(link1) == NULL);
	cactusLinkTestTeardown();
}

void testLink_getGroup(CuTest* testCase) {
	cactusLinkTestSetup();
	CuAssertTrue(testCase, link_getGroup(link1) == group1);
	CuAssertTrue(testCase, link_getGroup(link2) == group2);
	cactusLinkTestTeardown();
}

void testLink_getLeft(CuTest* testCase) {
	cactusLinkTestSetup();
	CuAssertTrue(testCase, link_get5End(link1) == end1);
	CuAssertTrue(testCase, link_get5End(link2) == block_get3End(block));
	cactusLinkTestTeardown();
}

void testLink_getRight(CuTest* testCase) {
	cactusLinkTestSetup();
	CuAssertTrue(testCase, link_get3End(link1) == block_get5End(block));
	CuAssertTrue(testCase, link_get3End(link2) == end2);
	cactusLinkTestTeardown();
}

void testLink_getChain(CuTest* testCase) {
	cactusLinkTestSetup();
	CuAssertTrue(testCase, link_getChain(link1) == chain);
	CuAssertTrue(testCase, link_getChain(link2) == chain);
	cactusLinkTestTeardown();
}

void testLink_getIndex(CuTest* testCase) {
	cactusLinkTestSetup();
	CuAssertIntEquals(testCase, 0, link_getIndex(link1));
	CuAssertIntEquals(testCase, 1, link_getIndex(link2));
	cactusLinkTestTeardown();
}

void testLink_split(CuTest *testCase) {
	cactusLinkTestSetup();
	CuAssertTrue(testCase, net_getChainNumber(net) == 1);
	CuAssertTrue(testCase, chain_getLength(chain) == 2);
	CuAssertTrue(testCase, group_getLink(group1) == link1);
	link_split(link1);
	CuAssertTrue(testCase, group_getLink(group1) == NULL);
	CuAssertTrue(testCase, net_getChainNumber(net) == 1);
	chain = net_getFirstChain(net);
	CuAssertTrue(testCase, chain_getLength(chain) == 1);
	Link *link3 = chain_getLink(chain, 0);
	CuAssertTrue(testCase, link_get5End(link3) == block_get3End(block));
	CuAssertTrue(testCase, link_get3End(link3) == end2);
	CuAssertTrue(testCase, group_getLink(group2) == link3);
	link_split(chain_getLink(chain, 0));
	CuAssertTrue(testCase, net_getChainNumber(net) == 0);
	CuAssertTrue(testCase, group_getLink(group1) == NULL);
	CuAssertTrue(testCase, group_getLink(group2) == NULL);
	cactusLinkTestTeardown();
}

void testLink_serialisation(CuTest* testCase) {
	cactusLinkTestSetup();
	int32_t i;
	void *vA = binaryRepresentation_makeBinaryRepresentation(link2,
			(void (*)(void *, void (*)(const void *, size_t, size_t)))link_writeBinaryRepresentation, &i);
	CuAssertTrue(testCase, i > 0);
	link_destruct(link2);
	void *vA2 = vA;
	link2 = link_loadFromBinaryRepresentation(&vA2, chain);
	nestedTest = 1;
	testLink_getNextLink(testCase);
	testLink_getPreviousLink(testCase);
	testLink_getGroup(testCase);
	testLink_getLeft(testCase);
	testLink_getRight(testCase);
	testLink_getChain(testCase);
	testLink_getIndex(testCase);
	//testLink_split(testCase); //can't do that test, because it disrupts stuff..
	nestedTest = 0;
	cactusLinkTestTeardown();
}

CuSuite* cactusLinkTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testLink_getNextLink);
	SUITE_ADD_TEST(suite, testLink_getPreviousLink);
	SUITE_ADD_TEST(suite, testLink_getGroup);
	SUITE_ADD_TEST(suite, testLink_getLeft);
	SUITE_ADD_TEST(suite, testLink_getRight);
	SUITE_ADD_TEST(suite, testLink_getChain);
	SUITE_ADD_TEST(suite, testLink_getIndex);
	SUITE_ADD_TEST(suite, testLink_serialisation);
	SUITE_ADD_TEST(suite, testLink_split);
	SUITE_ADD_TEST(suite, testLink_construct);
	return suite;
}
