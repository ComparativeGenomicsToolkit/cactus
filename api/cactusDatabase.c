#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Database functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int database_constructP(const char *vA1, int size1, const char *vA2, int size2, void *a) {
	assert(size1 == sizeof(Name));
	assert(size2 == sizeof(Name));
	assert(a == NULL);
	return cactusMisc_nameCompare(*(Name *)vA1, *(Name *)vA2);
}

TCBDB *database_construct(const char *name) {
	int32_t ecode;
	TCBDB *database;
	database = tcbdbnew();
	tcbdbsetcmpfunc(database, database_constructP, NULL);
	if(!tcbdbopen(database, name, BDBOWRITER | BDBOCREAT)) {
	   ecode = tcbdbecode(database);
	   exitOnFailure(1, "Opening database: %s with error: %s\n", name, tcbdberrmsg(ecode));
	}
	return database;
}

void database_destruct(TCBDB *database) {
	int32_t ecode;
	if(!tcbdbclose(database)){
		ecode = tcbdbecode(database);
		exitOnFailure(1, "Closing database error: %s\n", tcbdberrmsg(ecode));
	}
	tcbdbdel(database);
}

int32_t database_getNumberOfRecords(TCBDB *database) {
	return tcbdbrnum(database);
}

void *database_getRecord(TCBDB *database, Name key) {
	//Return value must be freed.
	int32_t i; //the size is ignored
	return tcbdbget(database, &key, sizeof(Name), &i);
}

int32_t database_writeRecord(TCBDB *database, Name key, const void *value, int32_t sizeOfRecord) {
	int32_t ecode;
	if(!tcbdbput(database, &key, sizeof(Name), value, sizeOfRecord)){
		ecode = tcbdbecode(database);
		fprintf(stderr, "Writing key/value to database error: %s\n", tcbdberrmsg(ecode));
		return 1;
	}
	return 0;
}

int32_t database_removeRecord(TCBDB *database, Name key) {
	int32_t ecode;
	if(!tcbdbout(database, &key, sizeof(Name))){
		ecode = tcbdbecode(database);
		fprintf(stderr, "Removing key/value to database error: %s\n", tcbdberrmsg(ecode));
		return 1;
	}
	return 0;
}

BDBCUR *databaseIterator_construct(TCBDB *database) {
	BDBCUR *iterator;
	iterator = tcbdbcurnew(database);
	tcbdbcurfirst(iterator);
	return iterator;
}

Name databaseIterator_getNext(BDBCUR *iterator) {
	const void *vA;
	int32_t i = 0;
	Name name = NULL_NAME;
	vA = tcbdbcurkey3(iterator, &i);
	if(vA != NULL) {
		assert(i == sizeof(Name));
		name = *(Name *)vA;
		tcbdbcurnext(iterator);
	}
	return name;
}

void databaseIterator_destruct(BDBCUR *iterator) {
	tcbdbcurdel(iterator);
}

int32_t database_startTransaction(TCBDB *database) {
	int32_t ecode;
	if(!tcbdbtranbegin(database)){
		ecode = tcbdbecode(database);
		fprintf(stderr, "Tried to start a transaction but got error: %s\n", tcbdberrmsg(ecode));
		return 1;
	}
	return 0;
}

int32_t database_commitTransaction(TCBDB *database) {
	int32_t ecode;
	if(!tcbdbtrancommit(database)){
		ecode = tcbdbecode(database);
		fprintf(stderr, "Tried to commit a transaction but got error: %s\n", tcbdberrmsg(ecode));
		return 1;
	}
	return 0;
}

int32_t database_abortTransaction(TCBDB *database) {
	int32_t ecode;
	if(!tcbdbtranabort(database)) {
		ecode = tcbdbecode(database);
		fprintf(stderr, "Tried to abort a transaction but got error: %s\n", tcbdberrmsg(ecode));
		return 1;
	}
	return 0;
}
