#include "cactusGlobalsPrivate.h"

static TCBDB *database = NULL;

static void cactusDatabaseTestTeardown() {
	if(database != NULL) {
		database_destruct(database);
		system("rm rf tempDatabase");
		database = NULL;
	}
}

static void cactusDatabaseTestSetup() {
	cactusDatabaseTestTeardown();
	database = database_construct("tempDatabase");
}

void testDatabase_construct(CuTest* testCase) {
	cactusDatabaseTestSetup();
	CuAssertTrue(testCase, database != NULL);
	cactusDatabaseTestTeardown();
}

void testDatabase(CuTest* testCase) {
	cactusDatabaseTestSetup();
	int32_t numberOfRecords = 100000;
	Name j;
	void *vA;
	void *vA2;

	/*
	 * Test the starting state.
	 */
	CuAssertIntEquals(testCase, 0, database_getNumberOfRecords(database));
	/*
	 * Try adding the values.
	 */
	for(j=0; j<numberOfRecords; j++) {
		CuAssertTrue(testCase, database_writeRecord(database, j, &j, sizeof(Name)) == 0);
		CuAssertIntEquals(testCase, j+1, database_getNumberOfRecords(database));
	}
	/*
	 * Try adding the values again, to check we don't get multiple keys.
	 */
	for(j=0; j<numberOfRecords; j++) {
		CuAssertTrue(testCase, database_writeRecord(database, j, &j, sizeof(Name)) == 0);
	}
	CuAssertIntEquals(testCase, numberOfRecords, database_getNumberOfRecords(database));
	/*
	 * Check we can retrieve all the names.
	 */
	for(j=0; j<numberOfRecords; j++) {
		vA = database_getRecord(database, j);
		vA2 = vA;
		CuAssertTrue(testCase, binaryRepresentation_getName(&vA2) == j);
		free(vA);
	}
	/*
	 * Check we can remove all the records.
	 */
	for(j=0; j<numberOfRecords; j++) {
		CuAssertTrue(testCase, database_removeRecord(database, j) == 0);
		CuAssertTrue(testCase, database_getRecord(database, j) == NULL);
		CuAssertIntEquals(testCase, numberOfRecords - 1 - j, database_getNumberOfRecords(database));
	}

	cactusDatabaseTestTeardown();
}

void testDatabaseIterator(CuTest* testCase) {
	cactusDatabaseTestSetup();
	int32_t numberOfRecords = 100000;
	Name j;

	for(j=0; j<numberOfRecords; j++) {
		CuAssertTrue(testCase, database_writeRecord(database, j, &j, sizeof(Name)) == 0);
		CuAssertIntEquals(testCase, j+1, database_getNumberOfRecords(database));
	}
	BDBCUR *iterator = databaseIterator_construct(database);

	for(j=0; j<numberOfRecords; j++) {
		CuAssertIntEquals(testCase, j, databaseIterator_getNext(iterator));
	}
	CuAssertTrue(testCase, databaseIterator_getNext(iterator) == NULL_NAME);

	databaseIterator_destruct(iterator);

	cactusDatabaseTestTeardown();
}


void testDatabase_transaction(CuTest* testCase) {
	cactusDatabaseTestSetup();
	int32_t numberOfRecords = 100000;
	Name j;

	/*
	 * Start the transaction.
	 */
	CuAssertTrue(testCase, database_startTransaction(database) == 0);
	/*
	 * Try adding the values.
	 */
	for(j=0; j<numberOfRecords; j++) {
		CuAssertTrue(testCase, database_writeRecord(database, j, &j, sizeof(Name)) == 0);
		CuAssertIntEquals(testCase, j+1, database_getNumberOfRecords(database));
	}

	/*
	 * Now abort the transaction.
	 */
	CuAssertTrue(testCase, database_abortTransaction(database) == 0);
	CuAssertIntEquals(testCase, 0, database_getNumberOfRecords(database));

	/*
	 * Start the transaction again.
	 */
	database_startTransaction(database);
	/*
	 * Try adding the values, again
	 */
	for(j=0; j<numberOfRecords; j++) {
		CuAssertTrue(testCase, database_writeRecord(database, j, &j, sizeof(Name)) == 0);
		CuAssertIntEquals(testCase, j+1, database_getNumberOfRecords(database));
	}
	/*
	 * Finish the transaction by commiting it.
	 */
	CuAssertTrue(testCase, database_commitTransaction(database) == 0);
	cactusDatabaseTestTeardown();
}


CuSuite* cactusDatabaseTestSuite(void) {
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, testDatabase);
	SUITE_ADD_TEST(suite, testDatabaseIterator);
	SUITE_ADD_TEST(suite, testDatabase_transaction);
	SUITE_ADD_TEST(suite, testDatabase_construct);
	return suite;
}
