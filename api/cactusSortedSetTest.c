#include "cactusGlobalsPrivate.h"

static struct avl_table *sortedSet = NULL;
static int32_t size = 9;
static int32_t input[] = { 1, 5, -1, 10, 3, 12, 3, -10, -10 };
static int32_t sortedInput[] = { -10, -1, 1, 3, 5, 10, 12 };
static int32_t sortedSize = 7;

static int32_t sortedSetTest_intCmp(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return *((int32_t *)o1) - *((int32_t *)o2);
}

static void sortedSetTest_destructElement(void *o, void *a) {
	assert(a == NULL);
	destructInt((int32_t *)o);
}

static void cactusSortedSetTestTeardown() {
	if(sortedSet != NULL) {
		sortedSet_destruct(sortedSet, sortedSetTest_destructElement);
		sortedSet = NULL;
	}
}

static void cactusSortedSetTestSetup() {
	sortedSet = sortedSet_construct(sortedSetTest_intCmp);
}

void testSortedSet_construct(CuTest* testCase) {
	cactusSortedSetTestSetup();
	CuAssertTrue(testCase, sortedSet != NULL);
	cactusSortedSetTestTeardown();
}

void testSortedSet(CuTest* testCase) {
	cactusSortedSetTestSetup();
	int32_t i;
	CuAssertIntEquals(testCase, 0, sortedSet_getLength(sortedSet));
	for(i=0; i<size; i++) {
		sortedSet_insert(sortedSet, constructInt(input[i]));
	}
	CuAssertIntEquals(testCase, sortedSize, sortedSet_getLength(sortedSet));
	CuAssertIntEquals(testCase, sortedInput[0], *(int32_t *)sortedSet_getFirst(sortedSet));
	CuAssertIntEquals(testCase, sortedInput[sortedSize-1], *(int32_t *)sortedSet_getLast(sortedSet));
	for(i=0; i<sortedSize; i++) {
		CuAssertIntEquals(testCase, sortedSize-i, sortedSet_getLength(sortedSet));
		CuAssertTrue(testCase, *((int32_t *)sortedSet_find(sortedSet, &sortedInput[i])) == sortedInput[i]);
		sortedSet_delete(sortedSet, &sortedInput[i]);
		CuAssertTrue(testCase, sortedSet_find(sortedSet, &sortedInput[i]) == NULL);
	}
	cactusSortedSetTestTeardown();
}

void testIterator_construct(CuTest* testCase) {
	cactusSortedSetTestSetup();
	struct avl_traverser *iterator = iterator_construct(sortedSet);
	CuAssertTrue(testCase, iterator != NULL);
	iterator_destruct(iterator);
	cactusSortedSetTestTeardown();
}

void testIterator(CuTest* testCase) {
	cactusSortedSetTestSetup();
	int32_t i;
	for(i=0; i<size; i++) {
		sortedSet_insert(sortedSet, constructInt(input[i]));
	}
	struct avl_traverser *iterator = iterator_construct(sortedSet);
	CuAssertTrue(testCase, iterator != NULL);

	for(i=0; i<sortedSize; i++) {
		CuAssertIntEquals(testCase, sortedInput[i], *(int32_t *)iterator_getNext(iterator));
	}
	CuAssertTrue(testCase, iterator_getNext(iterator) == NULL);
	struct avl_traverser *iterator2 = iterator_copy(iterator);
	CuAssertTrue(testCase, iterator2 != NULL);
	for(i=0; i<sortedSize; i++) {
		CuAssertIntEquals(testCase, sortedInput[sortedSize - 1 - i], *(int32_t *)iterator_getPrevious(iterator));
		CuAssertIntEquals(testCase, sortedInput[sortedSize - 1 - i], *(int32_t *)iterator_getPrevious(iterator2));
	}
	CuAssertTrue(testCase, iterator_getPrevious(iterator) == NULL);
	CuAssertTrue(testCase, iterator_getPrevious(iterator2) == NULL);
	iterator_destruct(iterator);
	iterator_destruct(iterator2);
	cactusSortedSetTestTeardown();
}

CuSuite* cactusSortedSetTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testSortedSet_construct);
	SUITE_ADD_TEST(suite, testSortedSet);
	SUITE_ADD_TEST(suite, testIterator);
	return suite;
}
