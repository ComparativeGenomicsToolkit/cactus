#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static int32_t cactusDisk_constructNetsP(const void *o1, const void *o2) {
	return cactusMisc_nameCompare(net_getName((Net *)o1), net_getName((Net *)o2));
}

static int32_t cactusDisk_constructMetaSequencesP(const void *o1, const void *o2) {
	return cactusMisc_nameCompare(metaSequence_getName((MetaSequence *)o1), metaSequence_getName((MetaSequence *)o2));
}

static int32_t cactusDisk_constructMetaEventsP(const void *o1, const void *o2) {
	return cactusMisc_nameCompare(metaEvent_getName((MetaEvent *)o1), metaEvent_getName((MetaEvent *)o2));
}

CactusDisk *cactusDisk_construct(const char *cactusDiskFile) {
	CactusDisk *cactusDisk;
	cactusDisk = st_malloc(sizeof(CactusDisk));
	int32_t i;

	//construct lists of in memory objects
	cactusDisk->metaEvents = stSortedSet_construct3(cactusDisk_constructMetaEventsP, NULL);
	cactusDisk->metaSequences = stSortedSet_construct3(cactusDisk_constructMetaSequencesP, NULL);
	cactusDisk->nets = stSortedSet_construct3(cactusDisk_constructNetsP, NULL);

	//the files to write the databases in
	cactusDisk->netsDatabaseName = pathJoin(cactusDiskFile, "nets");
	cactusDisk->metaDataDatabaseName = pathJoin(cactusDiskFile, "metaData");
	cactusDisk->iDDatabaseName = pathJoin(cactusDiskFile, "uniqueIDs");

	st_logInfo("Constructing the databases: %s, %s, %s\n",
			cactusDisk->iDDatabaseName, cactusDisk->netsDatabaseName, cactusDisk->metaDataDatabaseName);

	//create the net disk directory if it doesn't already exist.
	i = mkdir(cactusDiskFile, S_IRWXU);
	st_logInfo("Tried to create the base net disk directory with exit value: %i\n", i);

	//open the sequences database
	cactusDisk->netsDatabase = database_construct(cactusDisk->netsDatabaseName);
	cactusDisk->metaDataDatabase = database_construct(cactusDisk->metaDataDatabaseName);
	cactusDisk->iDDatabase = database_construct(cactusDisk->iDDatabaseName);

	//construct the string file
	cactusDisk->stringFile = pathJoin(cactusDiskFile, "strings");
	cactusDisk->stringFileLength = 0;

	//initialise the unique ids.
	//cactusDisk_getUniqueID(cactusDisk);
	cactusDisk->uniqueNumber = 0;
	cactusDisk->maxUniqueNumber = 0;

	return cactusDisk;
}

void cactusDisk_destruct(CactusDisk *cactusDisk){
	Net *net;
	MetaSequence *metaSequence;
	MetaEvent *metaEvent;

	while((net = cactusDisk_getFirstNetInMemory(cactusDisk)) != NULL) {
		net_destruct(net, FALSE);
	}
	stSortedSet_destruct(cactusDisk->nets);

	while((metaSequence = cactusDisk_getFirstMetaSequenceInMemory(cactusDisk)) != NULL) {
		metaSequence_destruct(metaSequence);
	}
	stSortedSet_destruct(cactusDisk->metaSequences);

	while((metaEvent = cactusDisk_getFirstMetaEventInMemory(cactusDisk)) != NULL) {
		metaEvent_destruct(metaEvent);
	}
	stSortedSet_destruct(cactusDisk->metaEvents);

	//close DBs
	database_destruct(cactusDisk->metaDataDatabase);
	database_destruct(cactusDisk->netsDatabase);
	database_destruct(cactusDisk->iDDatabase);

	//free string names
	free(cactusDisk->metaDataDatabaseName);
	free(cactusDisk->netsDatabaseName);
	free(cactusDisk->iDDatabaseName);
	free(cactusDisk->stringFile);

	free(cactusDisk);
}

void cactusDisk_write(CactusDisk *cactusDisk) {
	CactusDisk_NetIterator *netIterator;
	struct avl_traverser *metaDataIterator;
	void *vA;
	int32_t recordSize;
	Net *net;
	MetaSequence *metaSequence;
	MetaEvent *metaEvent;

	exitOnFailure(database_startTransaction(cactusDisk->netsDatabase), "Failed to start a transaction for the database: %s\n", cactusDisk->netsDatabaseName);
	exitOnFailure(database_startTransaction(cactusDisk->metaDataDatabase), "Failed to start a transaction for the database: %s\n", cactusDisk->metaDataDatabaseName);

	netIterator = cactusDisk_getNetsInMemoryIterator(cactusDisk);
	while((net = cactusDisk_getNextNet(netIterator)) != NULL) {
		vA = binaryRepresentation_makeBinaryRepresentation(net,
				(void (*)(void *, void (*)(const void * ptr, size_t size, size_t count)))net_writeBinaryRepresentation, &recordSize);
		exitOnFailure(database_writeRecord(cactusDisk->netsDatabase, net_getName(net), vA, recordSize),
				"Failed to write the net: %s to disk\n",
				cactusMisc_nameToStringStatic(net_getName(net)));
		free(vA);
	}
	cactusDisk_destructNetsInMemoryIterator(netIterator);

	metaDataIterator = stSortedSet_getIterator(cactusDisk->metaSequences);
	while((metaSequence = stSortedSet_getNext(metaDataIterator)) != NULL) {
		vA = binaryRepresentation_makeBinaryRepresentation(metaSequence,
				(void (*)(void *, void (*)(const void * ptr, size_t size, size_t count)))metaSequence_writeBinaryRepresentation, &recordSize);
		exitOnFailure(database_writeRecord(cactusDisk->metaDataDatabase, metaSequence_getName(metaSequence), vA, recordSize),
				"Failed to write meta sequence: %s to disk\n",
				cactusMisc_nameToStringStatic(metaSequence_getName(metaSequence)));
		free(vA);
	}
	stSortedSet_destructIterator(metaDataIterator);

	metaDataIterator = stSortedSet_getIterator(cactusDisk->metaEvents);
	while((metaEvent = stSortedSet_getNext(metaDataIterator)) != NULL) {
		vA = binaryRepresentation_makeBinaryRepresentation(metaEvent,
				(void (*)(void *, void (*)(const void *, size_t, size_t)))metaEvent_writeBinaryRepresentation, &recordSize);
		exitOnFailure(database_writeRecord(cactusDisk->metaDataDatabase, metaEvent_getName(metaEvent), vA, recordSize),
				"Failed to write meta event: %s to disk\n",
				cactusMisc_nameToStringStatic(metaEvent_getName(metaEvent)));
		free(vA);
	}
	stSortedSet_destructIterator(metaDataIterator);

	exitOnFailure(database_commitTransaction(cactusDisk->netsDatabase), "Failed to commit a transaction for the database: %s\n", cactusDisk->netsDatabaseName);
	exitOnFailure(database_commitTransaction(cactusDisk->metaDataDatabase), "Failed to commit a transaction for the database: %s\n", cactusDisk->metaDataDatabaseName);
}

void *cactusDisk_getObject(CactusDisk *cactusDisk, TCBDB *database, void *(*getObjectInMemory)(CactusDisk *, Name ),
		void *(*loadFromBinaryRepresentation)(void **, CactusDisk *),
		Name objectName) {
	void *cA;
	void *cA2;
	void *object;
	//try in memory list first.
	if((object = getObjectInMemory(cactusDisk, objectName)) != NULL) {
		return object;
	}
	//else try the database.
	cA = database_getRecord(database, objectName);
	if(cA == NULL) {
		return NULL;
	}
	else {
		cA2 = cA;
		object = loadFromBinaryRepresentation(&cA2, cactusDisk);
		free(cA);
		return object;
	}
}

Net *cactusDisk_getNet(CactusDisk *cactusDisk, Name netName) {
	return cactusDisk_getObject(cactusDisk, cactusDisk->netsDatabase,
			(void *(*)(CactusDisk *, Name ))cactusDisk_getNetInMemory,
			(void *(*)(void **, CactusDisk *))net_loadFromBinaryRepresentation,
			netName);
}

int32_t cactusDisk_getNetNumberOnDisk(CactusDisk *cactusDisk) {
	return database_getNumberOfRecords(cactusDisk->netsDatabase);
}

CactusDisk_NetNameIterator *cactusDisk_getNetNamesOnDiskIterator(CactusDisk *cactusDisk) {
	return databaseIterator_construct(cactusDisk->netsDatabase);
}

Name cactusDisk_getNextNetName(CactusDisk_NetNameIterator *netIterator) {
	return databaseIterator_getNext(netIterator);
}

void cactusDisk_destructNetNamesOnDiskIterator(CactusDisk_NetNameIterator *netIterator) {
	databaseIterator_destruct(netIterator);
}

int32_t cactusDisk_getNetNumberInMemory(CactusDisk *cactusDisk) {
	return stSortedSet_size(cactusDisk->nets);
}

CactusDisk_NetIterator *cactusDisk_getNetsInMemoryIterator(CactusDisk *cactusDisk) {
	return stSortedSet_getIterator(cactusDisk->nets);
}

Net *cactusDisk_getNextNet(CactusDisk_NetIterator *netIterator) {
	return stSortedSet_getNext(netIterator);
}

Net *cactusDisk_getPreviousNet(CactusDisk_NetIterator *netIterator) {
	return stSortedSet_getPrevious(netIterator);
}

CactusDisk_NetIterator *cactusDisk_copyNetIterator(CactusDisk_NetIterator *netIterator) {
	return stSortedSet_copyIterator(netIterator);
}

void cactusDisk_destructNetsInMemoryIterator(CactusDisk_NetIterator *netIterator) {
	stSortedSet_destructIterator(netIterator);
}

/*
 * Private functions.
 */

void cactusDisk_addNet(CactusDisk *cactusDisk, Net *net) {
	assert(stSortedSet_search(cactusDisk->nets, net) == NULL);
	stSortedSet_insert(cactusDisk->nets, net);
}

void cactusDisk_deleteNetFromDisk(CactusDisk *cactusDisk, Name netName) {
	if(database_getRecord(cactusDisk->netsDatabase, netName) != NULL) {
		exitOnFailure(database_removeRecord(cactusDisk->netsDatabase, netName),
				"Failed to remove the net: %s from the net disk database\n", cactusMisc_nameToStringStatic(netName));
	}
}

void cactusDisk_unloadNet(CactusDisk *cactusDisk, Net *net) {
	assert(cactusDisk_getNetInMemory(cactusDisk, net_getName(net)) != NULL);
	stSortedSet_remove(cactusDisk->nets, net);
}

Net *cactusDisk_getNetInMemory(CactusDisk *cactusDisk, Name netName) {
	static Net net;
	net.name = netName;
	return stSortedSet_search(cactusDisk->nets, &net);
}

Net *cactusDisk_getFirstNetInMemory(CactusDisk *cactusDisk) {
	return stSortedSet_getFirst(cactusDisk->nets);
}

/*
 * Functions on meta sequences.
 */

void cactusDisk_addMetaSequence(CactusDisk *cactusDisk, MetaSequence *metaSequence) {
	assert(stSortedSet_search(cactusDisk->metaSequences, metaSequence) == NULL);
	stSortedSet_insert(cactusDisk->metaSequences, metaSequence);
}

void cactusDisk_deleteMetaSequenceFromDisk(CactusDisk *cactusDisk, Name metaSequenceName) {
	if(database_getRecord(cactusDisk->metaDataDatabase, metaSequenceName) != NULL) {
		exitOnFailure(database_removeRecord(cactusDisk->metaDataDatabase, metaSequenceName),
				"Failed to remove the meta sequence: %s from the net disk database\n", cactusMisc_nameToStringStatic(metaSequenceName));
	}
}
void cactusDisk_unloadMetaSequence(CactusDisk *cactusDisk, MetaSequence *metaSequence) {
	stSortedSet_remove(cactusDisk->metaSequences, metaSequence);
}

MetaSequence *cactusDisk_getFirstMetaSequenceInMemory(CactusDisk *cactusDisk) {
	return stSortedSet_getFirst(cactusDisk->metaSequences);
}

MetaSequence *cactusDisk_getMetaSequenceInMemory(CactusDisk *cactusDisk, Name metaSequenceName) {
	static MetaSequence metaSequence;
	metaSequence.name = metaSequenceName;
	return stSortedSet_search(cactusDisk->metaSequences, &metaSequence);
}

MetaSequence *cactusDisk_getMetaSequence(CactusDisk *cactusDisk, Name metaSequenceName) {
return cactusDisk_getObject(cactusDisk, cactusDisk->metaDataDatabase,
		(void *(*)(CactusDisk *, Name ))cactusDisk_getMetaSequenceInMemory,
		(void *(*)(void **, CactusDisk *))metaSequence_loadFromBinaryRepresentation,
		metaSequenceName);
}

/*
 * Functions on meta events.
 */

void cactusDisk_addMetaEvent(CactusDisk *cactusDisk, MetaEvent *metaEvent) {
	assert(stSortedSet_search(cactusDisk->metaEvents, metaEvent) == NULL);
	stSortedSet_insert(cactusDisk->metaEvents, metaEvent);
}

void cactusDisk_deleteMetaEventFromDisk(CactusDisk *cactusDisk, Name metaEventName) {
	exitOnFailure(database_removeRecord(cactusDisk->metaDataDatabase, metaEventName),
			"Failed to remove the meta event: %s from the net disk database\n", cactusMisc_nameToStringStatic(metaEventName));
}

void cactusDisk_unloadMetaEvent(CactusDisk *cactusDisk, MetaEvent *metaEvent) {
	stSortedSet_remove(cactusDisk->metaEvents, metaEvent);
}

MetaEvent *cactusDisk_getFirstMetaEventInMemory(CactusDisk *cactusDisk) {
	return stSortedSet_getFirst(cactusDisk->metaEvents);
}

MetaEvent *cactusDisk_getMetaEventInMemory(CactusDisk *cactusDisk, Name metaEventName) {
	static MetaEvent metaEvent;
	metaEvent.name = metaEventName;
	return stSortedSet_search(cactusDisk->metaEvents, &metaEvent);
}

MetaEvent *cactusDisk_getMetaEvent(CactusDisk *cactusDisk, Name metaEventName) {
return cactusDisk_getObject(cactusDisk, cactusDisk->metaDataDatabase,
		(void *(*)(CactusDisk *, Name ))cactusDisk_getMetaEventInMemory,
		(void *(*)(void **, CactusDisk *))metaEvent_loadFromBinaryRepresentation,
		metaEventName);
}

/*
 * Functions on strings stored by the net disk.
 */

int64_t cactusDisk_addString(CactusDisk *cactusDisk, const char *string, int32_t length) {
	int64_t fileOffset;
	FILE *fileHandle;

	assert(length == (int32_t)strlen(string));
	fileOffset = cactusDisk->stringFileLength;
	cactusDisk->stringFileLength += length;
	fileHandle = fopen(cactusDisk->stringFile, "a");
	fprintf(fileHandle, "%s", string);
	fclose(fileHandle);
	return fileOffset;
}

char *cactusDisk_getString(CactusDisk *cactusDisk, int64_t offset, int32_t start, int32_t length, int32_t strand) {
	FILE *fileHandle;
	char *cA;
	char *cA2;

	fileHandle = fopen(cactusDisk->stringFile, "r");
	fseek(fileHandle, offset+start, SEEK_SET);
	cA = st_malloc(sizeof(char)*(length+1));
	fread(cA, sizeof(char), length, fileHandle);
	cA[length] = '\0';
	fclose(fileHandle);

	if(!strand) {
		cA2 = cactusMisc_reverseComplementString(cA);
		free(cA);
		return cA2;
	}
	return cA;
}

void cactusDisk_getBlockOfUniqueIDs(CactusDisk *cactusDisk) {
	void *vA;
	Name keyName;
	exitOnFailure(database_startTransaction(cactusDisk->iDDatabase), "Failed to start transaction to get a block of unique names\n");
	keyName = 0;
	vA = database_getRecord(cactusDisk->iDDatabase, keyName);
	if(vA == NULL) {
		cactusDisk->uniqueNumber = 0;
	}
	else {
		cactusDisk->uniqueNumber = *((Name *)vA);
	}
	cactusDisk->maxUniqueNumber = cactusDisk->uniqueNumber + CACTUSDISK_NAME_INCREMENT;
	exitOnFailure(database_writeRecord(cactusDisk->iDDatabase, keyName, &cactusDisk->maxUniqueNumber, sizeof(Name)), "Failed to update the new max unique name (%s) in the database\n", cactusMisc_nameToStringStatic(cactusDisk->maxUniqueNumber));
	exitOnFailure(database_commitTransaction(cactusDisk->iDDatabase), "Failed to commit a transaction to get a block of unique names\n");
	free(vA);
}

int64_t cactusDisk_getUniqueID(CactusDisk *cactusDisk) {
	assert(cactusDisk->uniqueNumber <= cactusDisk->maxUniqueNumber);
	if(cactusDisk->uniqueNumber == cactusDisk->maxUniqueNumber) {
		cactusDisk_getBlockOfUniqueIDs(cactusDisk);
	}
	return cactusDisk->uniqueNumber++;
}
