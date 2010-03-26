#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t netDisk_constructNetsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(net_getName((Net *)o1), net_getName((Net *)o2));
}

int32_t netDisk_constructMetaSequencesP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(metaSequence_getName((MetaSequence *)o1), metaSequence_getName((MetaSequence *)o2));
}

int32_t netDisk_constructMetaEventsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(metaEvent_getName((MetaEvent *)o1), metaEvent_getName((MetaEvent *)o2));
}

NetDisk *netDisk_construct(const char *netDiskFile) {
	NetDisk *netDisk;
	netDisk = malloc(sizeof(NetDisk));
	int32_t i;

	//construct lists of in memory objects
	netDisk->metaEvents = sortedSet_construct(netDisk_constructMetaEventsP);
	netDisk->metaSequences = sortedSet_construct(netDisk_constructMetaSequencesP);
	netDisk->nets = sortedSet_construct(netDisk_constructNetsP);

	//the files to write the databases in
	netDisk->netsDatabaseName = pathJoin(netDiskFile, "nets");
	netDisk->metaDataDatabaseName = pathJoin(netDiskFile, "metaData");
	netDisk->iDDatabaseName = pathJoin(netDiskFile, "uniqueIDs");

	logInfo("Constructing the databases: %s, %s, %s\n",
			netDisk->iDDatabaseName, netDisk->netsDatabaseName, netDisk->metaDataDatabaseName);

	//create the net disk directory if it doesn't already exist.
	i = mkdir(netDiskFile, S_IRWXU);
	logInfo("Tried to create the base net disk directory with exit value: %i\n", i);

	//open the sequences database
	netDisk->netsDatabase = database_construct(netDisk->netsDatabaseName);
	netDisk->metaDataDatabase = database_construct(netDisk->metaDataDatabaseName);
	netDisk->iDDatabase = database_construct(netDisk->iDDatabaseName);

	//construct the string file
	netDisk->stringFile = pathJoin(netDiskFile, "strings");
	netDisk->stringFileLength = 0;

	//initialise the unique ids.
	//netDisk_getUniqueID(netDisk);
	netDisk->uniqueNumber = 0;
	netDisk->maxUniqueNumber = 0;

	return netDisk;
}

void netDisk_destruct(NetDisk *netDisk){
	Net *net;
	MetaSequence *metaSequence;
	MetaEvent *metaEvent;

	while((net = netDisk_getFirstNetInMemory(netDisk)) != NULL) {
		net_destruct(net, FALSE);
	}
	sortedSet_destruct(netDisk->nets, NULL);

	while((metaSequence = netDisk_getFirstMetaSequenceInMemory(netDisk)) != NULL) {
		metaSequence_destruct(metaSequence);
	}
	sortedSet_destruct(netDisk->metaSequences, NULL);

	while((metaEvent = netDisk_getFirstMetaEventInMemory(netDisk)) != NULL) {
		metaEvent_destruct(metaEvent);
	}
	sortedSet_destruct(netDisk->metaEvents, NULL);

	//close DBs
	database_destruct(netDisk->metaDataDatabase);
	database_destruct(netDisk->netsDatabase);
	database_destruct(netDisk->iDDatabase);

	//free string names
	free(netDisk->metaDataDatabaseName);
	free(netDisk->netsDatabaseName);
	free(netDisk->iDDatabaseName);
	free(netDisk->stringFile);

	free(netDisk);
}

void netDisk_write(NetDisk *netDisk) {
	NetDisk_NetIterator *netIterator;
	struct avl_traverser *metaDataIterator;
	void *vA;
	int32_t recordSize;
	Net *net;
	MetaSequence *metaSequence;
	MetaEvent *metaEvent;

	exitOnFailure(database_startTransaction(netDisk->netsDatabase), "Failed to start a transaction for the database: %s\n", netDisk->netsDatabaseName);
	exitOnFailure(database_startTransaction(netDisk->metaDataDatabase), "Failed to start a transaction for the database: %s\n", netDisk->metaDataDatabaseName);

	netIterator = netDisk_getNetsInMemoryIterator(netDisk);
	while((net = netDisk_getNextNet(netIterator)) != NULL) {
		vA = binaryRepresentation_makeBinaryRepresentation(net,
				(void (*)(void *, void (*)(const void * ptr, size_t size, size_t count)))net_writeBinaryRepresentation, &recordSize);
		exitOnFailure(database_writeRecord(netDisk->netsDatabase, net_getName(net), vA, recordSize),
				"Failed to write the net: %s to disk\n",
				netMisc_nameToStringStatic(net_getName(net)));
		free(vA);
	}
	netDisk_destructNetsInMemoryIterator(netIterator);

	metaDataIterator = iterator_construct(netDisk->metaSequences);
	while((metaSequence = iterator_getNext(metaDataIterator)) != NULL) {
		vA = binaryRepresentation_makeBinaryRepresentation(metaSequence,
				(void (*)(void *, void (*)(const void * ptr, size_t size, size_t count)))metaSequence_writeBinaryRepresentation, &recordSize);
		exitOnFailure(database_writeRecord(netDisk->metaDataDatabase, metaSequence_getName(metaSequence), vA, recordSize),
				"Failed to write meta sequence: %s to disk\n",
				netMisc_nameToStringStatic(metaSequence_getName(metaSequence)));
		free(vA);
	}
	iterator_destruct(metaDataIterator);

	metaDataIterator = iterator_construct(netDisk->metaEvents);
	while((metaEvent = iterator_getNext(metaDataIterator)) != NULL) {
		vA = binaryRepresentation_makeBinaryRepresentation(metaEvent,
				(void (*)(void *, void (*)(const void *, size_t, size_t)))metaEvent_writeBinaryRepresentation, &recordSize);
		exitOnFailure(database_writeRecord(netDisk->metaDataDatabase, metaEvent_getName(metaEvent), vA, recordSize),
				"Failed to write meta event: %s to disk\n",
				netMisc_nameToStringStatic(metaEvent_getName(metaEvent)));
		free(vA);
	}
	iterator_destruct(metaDataIterator);

	exitOnFailure(database_commitTransaction(netDisk->netsDatabase), "Failed to commit a transaction for the database: %s\n", netDisk->netsDatabaseName);
	exitOnFailure(database_commitTransaction(netDisk->metaDataDatabase), "Failed to commit a transaction for the database: %s\n", netDisk->metaDataDatabaseName);
}

void *netDisk_getObject(NetDisk *netDisk, TCBDB *database, void *(*getObjectInMemory)(NetDisk *, Name ),
		void *(*loadFromBinaryRepresentation)(void **, NetDisk *),
		Name objectName) {
	void *cA;
	void *cA2;
	void *object;
	//try in memory list first.
	if((object = getObjectInMemory(netDisk, objectName)) != NULL) {
		return object;
	}
	//else try the database.
	cA = database_getRecord(database, objectName);
	if(cA == NULL) {
		return NULL;
	}
	else {
		cA2 = cA;
		object = loadFromBinaryRepresentation(&cA2, netDisk);
		free(cA);
		return object;
	}
}

Net *netDisk_getNet(NetDisk *netDisk, Name netName) {
	return netDisk_getObject(netDisk, netDisk->netsDatabase,
			(void *(*)(NetDisk *, Name ))netDisk_getNetInMemory,
			(void *(*)(void **, NetDisk *))net_loadFromBinaryRepresentation,
			netName);
}

int32_t netDisk_getNetNumberOnDisk(NetDisk *netDisk) {
	return database_getNumberOfRecords(netDisk->netsDatabase);
}

NetDisk_NetNameIterator *netDisk_getNetNamesOnDiskIterator(NetDisk *netDisk) {
	return databaseIterator_construct(netDisk->netsDatabase);
}

Name netDisk_getNextNetName(NetDisk_NetNameIterator *netIterator) {
	return databaseIterator_getNext(netIterator);
}

void netDisk_destructNetNamesOnDiskIterator(NetDisk_NetNameIterator *netIterator) {
	databaseIterator_destruct(netIterator);
}

int32_t netDisk_getNetNumberInMemory(NetDisk *netDisk) {
	return sortedSet_getLength(netDisk->nets);
}

NetDisk_NetIterator *netDisk_getNetsInMemoryIterator(NetDisk *netDisk) {
	return iterator_construct(netDisk->nets);
}

Net *netDisk_getNextNet(NetDisk_NetIterator *netIterator) {
	return iterator_getNext(netIterator);
}

Net *netDisk_getPreviousNet(NetDisk_NetIterator *netIterator) {
	return iterator_getPrevious(netIterator);
}

NetDisk_NetIterator *netDisk_copyNetIterator(NetDisk_NetIterator *netIterator) {
	return iterator_copy(netIterator);
}

void netDisk_destructNetsInMemoryIterator(NetDisk_NetIterator *netIterator) {
	iterator_destruct(netIterator);
}

/*
 * Private functions.
 */

void netDisk_addNet(NetDisk *netDisk, Net *net) {
	assert(sortedSet_find(netDisk->nets, net) == NULL);
	sortedSet_insert(netDisk->nets, net);
}

void netDisk_deleteNetFromDisk(NetDisk *netDisk, Name netName) {
	exitOnFailure(database_removeRecord(netDisk->netsDatabase, netName),
			"Failed to remove the net: %s from the net disk database\n", netMisc_nameToStringStatic(netName));
}

void netDisk_unloadNet(NetDisk *netDisk, Net *net) {
	assert(netDisk_getNetInMemory(netDisk, net_getName(net)) != NULL);
	sortedSet_delete(netDisk->nets, net);
}

Net *netDisk_getNetInMemory(NetDisk *netDisk, Name netName) {
	static Net net;
	net.name = netName;
	return sortedSet_find(netDisk->nets, &net);
}

Net *netDisk_getFirstNetInMemory(NetDisk *netDisk) {
	return sortedSet_getFirst(netDisk->nets);
}

/*
 * Functions on meta sequences.
 */

void netDisk_addMetaSequence(NetDisk *netDisk, MetaSequence *metaSequence) {
	assert(sortedSet_find(netDisk->metaSequences, metaSequence) == NULL);
	sortedSet_insert(netDisk->metaSequences, metaSequence);
}

void netDisk_deleteMetaSequenceFromDisk(NetDisk *netDisk, Name metaSequenceName) {
	exitOnFailure(database_removeRecord(netDisk->metaDataDatabase, metaSequenceName),
			"Failed to remove the meta sequence: %s from the net disk database\n", netMisc_nameToStringStatic(metaSequenceName));
}

void netDisk_unloadMetaSequence(NetDisk *netDisk, MetaSequence *metaSequence) {
	sortedSet_delete(netDisk->metaSequences, metaSequence);
}

MetaSequence *netDisk_getFirstMetaSequenceInMemory(NetDisk *netDisk) {
	return sortedSet_getFirst(netDisk->metaSequences);
}

MetaSequence *netDisk_getMetaSequenceInMemory(NetDisk *netDisk, Name metaSequenceName) {
	static MetaSequence metaSequence;
	metaSequence.name = metaSequenceName;
	return sortedSet_find(netDisk->metaSequences, &metaSequence);
}

MetaSequence *netDisk_getMetaSequence(NetDisk *netDisk, Name metaSequenceName) {
return netDisk_getObject(netDisk, netDisk->metaDataDatabase,
		(void *(*)(NetDisk *, Name ))netDisk_getMetaSequenceInMemory,
		(void *(*)(void **, NetDisk *))metaSequence_loadFromBinaryRepresentation,
		metaSequenceName);
}

/*
 * Functions on meta events.
 */

void netDisk_addMetaEvent(NetDisk *netDisk, MetaEvent *metaEvent) {
	assert(sortedSet_find(netDisk->metaEvents, metaEvent) == NULL);
	sortedSet_insert(netDisk->metaEvents, metaEvent);
}

void netDisk_deleteMetaEventFromDisk(NetDisk *netDisk, Name metaEventName) {
	exitOnFailure(database_removeRecord(netDisk->metaDataDatabase, metaEventName),
			"Failed to remove the meta event: %s from the net disk database\n", netMisc_nameToStringStatic(metaEventName));
}

void netDisk_unloadMetaEvent(NetDisk *netDisk, MetaEvent *metaEvent) {
	sortedSet_delete(netDisk->metaEvents, metaEvent);
}

MetaEvent *netDisk_getFirstMetaEventInMemory(NetDisk *netDisk) {
	return sortedSet_getFirst(netDisk->metaEvents);
}

MetaEvent *netDisk_getMetaEventInMemory(NetDisk *netDisk, Name metaEventName) {
	static MetaEvent metaEvent;
	metaEvent.name = metaEventName;
	return sortedSet_find(netDisk->metaEvents, &metaEvent);
}

MetaEvent *netDisk_getMetaEvent(NetDisk *netDisk, Name metaEventName) {
return netDisk_getObject(netDisk, netDisk->metaDataDatabase,
		(void *(*)(NetDisk *, Name ))netDisk_getMetaEventInMemory,
		(void *(*)(void **, NetDisk *))metaEvent_loadFromBinaryRepresentation,
		metaEventName);
}

/*
 * Functions on strings stored by the net disk.
 */

int64_t netDisk_addString(NetDisk *netDisk, const char *string, int32_t length) {
	int64_t fileOffset;
	FILE *fileHandle;

	assert(length == (int32_t)strlen(string));
	fileOffset = netDisk->stringFileLength;
	netDisk->stringFileLength += length;
	fileHandle = fopen(netDisk->stringFile, "a");
	fprintf(fileHandle, "%s", string);
	fclose(fileHandle);
	return fileOffset;
}

char *netDisk_getString(NetDisk *netDisk, int64_t offset, int32_t start, int32_t length, int32_t strand) {
	FILE *fileHandle;
	char *cA;
	char *cA2;

	fileHandle = fopen(netDisk->stringFile, "r");
	fseek(fileHandle, offset+start, SEEK_SET);
	cA = malloc(sizeof(char)*(length+1));
	fread(cA, sizeof(char), length, fileHandle);
	cA[length] = '\0';
	fclose(fileHandle);

	if(!strand) {
		cA2 = netMisc_reverseComplementString(cA);
		free(cA);
		return cA2;
	}
	return cA;
}

void netDisk_getBlockOfUniqueIDs(NetDisk *netDisk) {
	void *vA;
	Name keyName;
	exitOnFailure(database_startTransaction(netDisk->iDDatabase), "Failed to start transaction to get a block of unique names\n");
	keyName = 0;
	vA = database_getRecord(netDisk->iDDatabase, keyName);
	if(vA == NULL) {
		netDisk->uniqueNumber = 0;
	}
	else {
		netDisk->uniqueNumber = *((Name *)vA);
	}
	netDisk->maxUniqueNumber = netDisk->uniqueNumber + NETDISK_NAME_INCREMENT;
	exitOnFailure(database_writeRecord(netDisk->iDDatabase, keyName, &netDisk->maxUniqueNumber, sizeof(Name)), "Failed to update the new max unique name (%s) in the database\n", netMisc_nameToStringStatic(netDisk->maxUniqueNumber));
	exitOnFailure(database_commitTransaction(netDisk->iDDatabase), "Failed to commit a transaction to get a block of unique names\n");
	free(vA);
}

int64_t netDisk_getUniqueID(NetDisk *netDisk) {
	assert(netDisk->uniqueNumber <= netDisk->maxUniqueNumber);
	if(netDisk->uniqueNumber == netDisk->maxUniqueNumber) {
		netDisk_getBlockOfUniqueIDs(netDisk);
	}
	return netDisk->uniqueNumber++;
}
