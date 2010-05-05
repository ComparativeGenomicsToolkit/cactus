#include "cactusChainsTestShared.h"

static Hash *hash;
static Hash *hash2;
static int32_t *one, *two, *three, *four, *five, *six;

static uint32_t hashKey(void *o) {
	return *((int32_t *)o);
}

static int32_t hashEqualsKey(void *o, void *o2) {
	return *((int32_t *)o) == *((int32_t *)o2);
}

static void destructKey(void *o) {
	destructInt(o);
}

static void destructValue(void *o) {
	destructInt(o);
}

static void testSetup() {
	//compare by value of memory address
	hash = hash_construct();
	//compare by value of ints.
	hash2 = hash_construct3(hashKey, hashEqualsKey, destructKey, destructValue);
	one = constructInt(0);
	two = constructInt(1);
	three = constructInt(2);
	four = constructInt(3);
	five = constructInt(4);
	six = constructInt(5);

	hash_insert(hash, one, two);
	hash_insert(hash, three, four);
	hash_insert(hash, five, six);

	hash_insert(hash2, one, two);
	hash_insert(hash2, three, four);
	hash_insert(hash2, five, six);
}

static void testTeardown() {
	hash_destruct(hash);
	hash_destruct(hash2);
}

void testHash_construct(CuTest* testCase) {
	testSetup();
	/* Do nothing */
	testTeardown();
}

void testHash_search(CuTest* testCase) {
	testSetup();

	int32_t *i = constructInt(0);

	//Check search by memory address
	CuAssertTrue(testCase, hash_search(hash, one) == two);
	CuAssertTrue(testCase, hash_search(hash, three) == four);
	CuAssertTrue(testCase, hash_search(hash, five) == six);
	//Check not present
	CuAssertTrue(testCase, hash_search(hash, six) == NULL);
	CuAssertTrue(testCase, hash_search(hash, i) == NULL);

	//Check search by memory address
	CuAssertTrue(testCase, hash_search(hash2, one) == two);
	CuAssertTrue(testCase, hash_search(hash2, three) == four);
	CuAssertTrue(testCase, hash_search(hash2, five) == six);
	//Check not present
	CuAssertTrue(testCase, hash_search(hash2, six) == NULL);
	//Check is searching by memory.
	CuAssertTrue(testCase, hash_search(hash2, i) == two);

	destructInt(i);

	testTeardown();
}

void testHash_remove(CuTest* testCase) {
	testSetup();

	CuAssertTrue(testCase, hash_remove(hash, one) == two);
	CuAssertTrue(testCase, hash_search(hash, one) == NULL);

	CuAssertTrue(testCase, hash_remove(hash2, one) == two);
	CuAssertTrue(testCase, hash_search(hash2, one) == NULL);

	hash_insert(hash2, one, two);
	CuAssertTrue(testCase, hash_search(hash2, one) == two);

	testTeardown();
}

void testHash_insert(CuTest* testCase) {
	/*
	 * Tests inserting already present keys.
	 */
	testSetup();

	CuAssertTrue(testCase, hash_search(hash, one) == two);
	hash_insert(hash, one, two);
	CuAssertTrue(testCase, hash_search(hash, one) == two);
	hash_insert(hash, one, three);
	CuAssertTrue(testCase, hash_search(hash, one) == three);
	hash_insert(hash, one, two);
	CuAssertTrue(testCase, hash_search(hash, one) == two);

	testTeardown();
}

void testHash_size(CuTest *testCase) {
	/*
	 * Tests the size function of the hash.
	 */
	testSetup();

	CuAssertTrue(testCase, hash_size(hash) == 3);
	CuAssertTrue(testCase, hash_size(hash2) == 3);
	Hash *hash3 = hash_construct();
	CuAssertTrue(testCase, hash_size(hash3) == 0);
	hash_destruct(hash3);

	testTeardown();
}

CuSuite* cactusHashTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testHash_search);
	SUITE_ADD_TEST(suite, testHash_remove);
	SUITE_ADD_TEST(suite, testHash_insert);
	SUITE_ADD_TEST(suite, testHash_size);
	SUITE_ADD_TEST(suite, testHash_construct);
	return suite;
}
