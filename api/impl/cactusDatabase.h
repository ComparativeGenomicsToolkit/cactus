#ifndef CACTUS_DATABASE_H_
#define CACTUS_DATABASE_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Database functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a sorted-set database object.
 */
TCBDB *database_construct(const char *name);

/*
 * Destructs a database.
 */
void database_destruct(TCBDB *database);

/*
 * Returns number of records in database.
 */
int32_t database_getNumberOfRecords(TCBDB *database);

/*
 * Gets a record from the database, given the key. The record is in newly allocated memory, and must be freed.
 */
void *database_getRecord(TCBDB *database, Name key);

/*
 * Writes a key value record to the database. Returns zero if successful.
 */
int32_t database_writeRecord(TCBDB *database, Name key, const void *value, int32_t sizeOfRecord);

/*
 * Removes a record from the database. Returns zero if successful.
 */
int32_t database_removeRecord(TCBDB *database, Name key);

/*
 * Constructs an iterator over the sorted database records.
 */
BDBCUR *databaseIterator_construct(TCBDB *database);

/*
 * Gets the next element from the database iterator.
 */
Name databaseIterator_getNext(BDBCUR *iterator);

/*
 * Destructs a database iterator.
 */
void databaseIterator_destruct(BDBCUR *iterator);

/*
 * Starts a transaction with the database. Returns zero for success.
 */
int32_t database_startTransaction(TCBDB *database);

/*
 * Commits the transaction to the database. Returns zero for success.
 */
int32_t database_commitTransaction(TCBDB *database);

/*
 * Aborts the transaction to the database. Returns zero for success.
 */
int32_t database_abortTransaction(TCBDB *database);

#endif
