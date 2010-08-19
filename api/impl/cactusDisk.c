#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic flower disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static int32_t cactusDisk_constructFlowersP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(flower_getName((Flower *) o1),
            flower_getName((Flower *) o2));
}

static int32_t cactusDisk_constructMetaSequencesP(const void *o1,
        const void *o2) {
    return cactusMisc_nameCompare(metaSequence_getName((MetaSequence *) o1),
            metaSequence_getName((MetaSequence *) o2));
}

static int32_t cactusDisk_constructMetaEventsP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(metaEvent_getName((MetaEvent *) o1),
            metaEvent_getName((MetaEvent *) o2));
}

CactusDisk *cactusDisk_construct(const char *cactusDiskFile) {
    CactusDisk *cactusDisk;
    cactusDisk = st_malloc(sizeof(CactusDisk));
    int32_t i;

    //construct lists of in memory objects
    cactusDisk->metaEvents = stSortedSet_construct3(
            cactusDisk_constructMetaEventsP, NULL);
    cactusDisk->metaSequences = stSortedSet_construct3(
            cactusDisk_constructMetaSequencesP, NULL);
    cactusDisk->flowers = stSortedSet_construct3(cactusDisk_constructFlowersP,
            NULL);

    //the files to write the databases in
    cactusDisk->flowersDatabaseName = pathJoin(cactusDiskFile, "flowers");
    cactusDisk->metaDataDatabaseName = pathJoin(cactusDiskFile, "metaData");
    cactusDisk->iDDatabaseName = pathJoin(cactusDiskFile, "uniqueIDs");

    st_logInfo("Constructing the databases: %s, %s, %s\n",
            cactusDisk->iDDatabaseName, cactusDisk->flowersDatabaseName,
            cactusDisk->metaDataDatabaseName);

    //create the flower disk directory if it doesn't already exist.
    i = mkdir(cactusDiskFile, S_IRWXU);
    st_logInfo(
            "Tried to create the base flower disk directory with exit value: %i\n",
            i);

    //open the sequences database
    cactusDisk->flowersDatabase = database_construct(
            cactusDisk->flowersDatabaseName);
    cactusDisk->metaDataDatabase = database_construct(
            cactusDisk->metaDataDatabaseName);
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

void cactusDisk_destruct(CactusDisk *cactusDisk) {
    Flower *flower;
    MetaSequence *metaSequence;
    MetaEvent *metaEvent;

    while ((flower = cactusDisk_getFirstFlowerInMemory(cactusDisk)) != NULL) {
        flower_destruct(flower, FALSE);
    }
    stSortedSet_destruct(cactusDisk->flowers);

    while ((metaSequence = cactusDisk_getFirstMetaSequenceInMemory(cactusDisk))
            != NULL) {
        metaSequence_destruct(metaSequence);
    }
    stSortedSet_destruct(cactusDisk->metaSequences);

    while ((metaEvent = cactusDisk_getFirstMetaEventInMemory(cactusDisk))
            != NULL) {
        metaEvent_destruct(metaEvent);
    }
    stSortedSet_destruct(cactusDisk->metaEvents);

    //close DBs
    database_destruct(cactusDisk->metaDataDatabase);
    database_destruct(cactusDisk->flowersDatabase);
    database_destruct(cactusDisk->iDDatabase);

    //free string names
    free(cactusDisk->metaDataDatabaseName);
    free(cactusDisk->flowersDatabaseName);
    free(cactusDisk->iDDatabaseName);
    free(cactusDisk->stringFile);

    free(cactusDisk);
}

void cactusDisk_write(CactusDisk *cactusDisk) {
    CactusDisk_FlowerIterator *flowerIterator;
    struct avl_traverser *metaDataIterator;
    void *vA;
    int32_t recordSize;
    Flower *flower;
    MetaSequence *metaSequence;
    MetaEvent *metaEvent;

    exitOnFailure(database_startTransaction(cactusDisk->flowersDatabase),
            "Failed to start a transaction for the database: %s\n",
            cactusDisk->flowersDatabaseName);
    exitOnFailure(database_startTransaction(cactusDisk->metaDataDatabase),
            "Failed to start a transaction for the database: %s\n",
            cactusDisk->metaDataDatabaseName);

    flowerIterator = cactusDisk_getFlowersInMemoryIterator(cactusDisk);
    while ((flower = cactusDisk_getNextFlower(flowerIterator)) != NULL) {
        vA = binaryRepresentation_makeBinaryRepresentation(flower,
                (void(*)(void *, void(*)(const void * ptr, size_t size,
                        size_t count))) flower_writeBinaryRepresentation,
                &recordSize);
        exitOnFailure(database_writeRecord(cactusDisk->flowersDatabase,
                flower_getName(flower), vA, recordSize),
                "Failed to write the flower: %s to disk\n",
                cactusMisc_nameToStringStatic(flower_getName(flower)));
        free(vA);
    }
    cactusDisk_destructFlowersInMemoryIterator(flowerIterator);

    metaDataIterator = stSortedSet_getIterator(cactusDisk->metaSequences);
    while ((metaSequence = stSortedSet_getNext(metaDataIterator)) != NULL) {
        vA = binaryRepresentation_makeBinaryRepresentation(metaSequence,
                (void(*)(void *, void(*)(const void * ptr, size_t size,
                        size_t count))) metaSequence_writeBinaryRepresentation,
                &recordSize);
        exitOnFailure(database_writeRecord(cactusDisk->metaDataDatabase,
                metaSequence_getName(metaSequence), vA, recordSize),
                "Failed to write meta sequence: %s to disk\n",
                cactusMisc_nameToStringStatic(
                        metaSequence_getName(metaSequence)));
        free(vA);
    }
    stSortedSet_destructIterator(metaDataIterator);

    metaDataIterator = stSortedSet_getIterator(cactusDisk->metaEvents);
    while ((metaEvent = stSortedSet_getNext(metaDataIterator)) != NULL) {
        vA
                = binaryRepresentation_makeBinaryRepresentation(
                        metaEvent,
                        (void(*)(void *, void(*)(const void *, size_t, size_t))) metaEvent_writeBinaryRepresentation,
                        &recordSize);
        exitOnFailure(database_writeRecord(cactusDisk->metaDataDatabase,
                metaEvent_getName(metaEvent), vA, recordSize),
                "Failed to write meta event: %s to disk\n",
                cactusMisc_nameToStringStatic(metaEvent_getName(metaEvent)));
        free(vA);
    }
    stSortedSet_destructIterator(metaDataIterator);

    exitOnFailure(database_commitTransaction(cactusDisk->flowersDatabase),
            "Failed to commit a transaction for the database: %s\n",
            cactusDisk->flowersDatabaseName);
    exitOnFailure(database_commitTransaction(cactusDisk->metaDataDatabase),
            "Failed to commit a transaction for the database: %s\n",
            cactusDisk->metaDataDatabaseName);
}

void *cactusDisk_getObject(CactusDisk *cactusDisk, TCBDB *database,
        void *(*getObjectInMemory)(CactusDisk *, Name),
        void *(*loadFromBinaryRepresentation)(void **, CactusDisk *),
        Name objectName) {
    void *cA;
    void *cA2;
    void *object;
    //try in memory list first.
    if ((object = getObjectInMemory(cactusDisk, objectName)) != NULL) {
        return object;
    }
    //else try the database.
    cA = database_getRecord(database, objectName);
    if (cA == NULL) {
        return NULL;
    } else {
        cA2 = cA;
        object = loadFromBinaryRepresentation(&cA2, cactusDisk);
        free(cA);
        return object;
    }
}

Flower *cactusDisk_getFlower(CactusDisk *cactusDisk, Name flowerName) {
    return cactusDisk_getObject(
            cactusDisk,
            cactusDisk->flowersDatabase,
            (void *(*)(CactusDisk *, Name)) cactusDisk_getFlowerInMemory,
            (void *(*)(void **, CactusDisk *)) flower_loadFromBinaryRepresentation,
            flowerName);
}

int32_t cactusDisk_getFlowerNumberOnDisk(CactusDisk *cactusDisk) {
    return database_getNumberOfRecords(cactusDisk->flowersDatabase);
}

CactusDisk_FlowerNameIterator *cactusDisk_getFlowerNamesOnDiskIterator(
        CactusDisk *cactusDisk) {
    return databaseIterator_construct(cactusDisk->flowersDatabase);
}

Name cactusDisk_getNextFlowerName(CactusDisk_FlowerNameIterator *flowerIterator) {
    return databaseIterator_getNext(flowerIterator);
}

void cactusDisk_destructFlowerNamesOnDiskIterator(
        CactusDisk_FlowerNameIterator *flowerIterator) {
    databaseIterator_destruct(flowerIterator);
}

int32_t cactusDisk_getFlowerNumberInMemory(CactusDisk *cactusDisk) {
    return stSortedSet_size(cactusDisk->flowers);
}

CactusDisk_FlowerIterator *cactusDisk_getFlowersInMemoryIterator(
        CactusDisk *cactusDisk) {
    return stSortedSet_getIterator(cactusDisk->flowers);
}

Flower *cactusDisk_getNextFlower(CactusDisk_FlowerIterator *flowerIterator) {
    return stSortedSet_getNext(flowerIterator);
}

Flower *cactusDisk_getPreviousFlower(CactusDisk_FlowerIterator *flowerIterator) {
    return stSortedSet_getPrevious(flowerIterator);
}

CactusDisk_FlowerIterator *cactusDisk_copyFlowerIterator(
        CactusDisk_FlowerIterator *flowerIterator) {
    return stSortedSet_copyIterator(flowerIterator);
}

void cactusDisk_destructFlowersInMemoryIterator(
        CactusDisk_FlowerIterator *flowerIterator) {
    stSortedSet_destructIterator(flowerIterator);
}

/*
 * Private functions.
 */

void cactusDisk_addFlower(CactusDisk *cactusDisk, Flower *flower) {
    assert(stSortedSet_search(cactusDisk->flowers, flower) == NULL);
    stSortedSet_insert(cactusDisk->flowers, flower);
}

void cactusDisk_deleteFlowerFromDisk(CactusDisk *cactusDisk, Name flowerName) {
    if (database_getRecord(cactusDisk->flowersDatabase, flowerName) != NULL) {
        exitOnFailure(
                database_removeRecord(cactusDisk->flowersDatabase, flowerName),
                "Failed to remove the flower: %s from the flower disk database\n",
                cactusMisc_nameToStringStatic(flowerName));
    }
}

void cactusDisk_unloadFlower(CactusDisk *cactusDisk, Flower *flower) {
    assert(cactusDisk_getFlowerInMemory(cactusDisk, flower_getName(flower)) != NULL);
    stSortedSet_remove(cactusDisk->flowers, flower);
}

Flower *cactusDisk_getFlowerInMemory(CactusDisk *cactusDisk, Name flowerName) {
    static Flower flower;
    flower.name = flowerName;
    return stSortedSet_search(cactusDisk->flowers, &flower);
}

Flower *cactusDisk_getFirstFlowerInMemory(CactusDisk *cactusDisk) {
    return stSortedSet_getFirst(cactusDisk->flowers);
}

/*
 * Functions on meta sequences.
 */

void cactusDisk_addMetaSequence(CactusDisk *cactusDisk,
        MetaSequence *metaSequence) {
    assert(stSortedSet_search(cactusDisk->metaSequences, metaSequence) == NULL);
    stSortedSet_insert(cactusDisk->metaSequences, metaSequence);
}

void cactusDisk_deleteMetaSequenceFromDisk(CactusDisk *cactusDisk,
        Name metaSequenceName) {
    if (database_getRecord(cactusDisk->metaDataDatabase, metaSequenceName)
            != NULL) {
        exitOnFailure(
                database_removeRecord(cactusDisk->metaDataDatabase,
                        metaSequenceName),
                "Failed to remove the meta sequence: %s from the flower disk database\n",
                cactusMisc_nameToStringStatic(metaSequenceName));
    }
}
void cactusDisk_unloadMetaSequence(CactusDisk *cactusDisk,
        MetaSequence *metaSequence) {
    stSortedSet_remove(cactusDisk->metaSequences, metaSequence);
}

MetaSequence *cactusDisk_getFirstMetaSequenceInMemory(CactusDisk *cactusDisk) {
    return stSortedSet_getFirst(cactusDisk->metaSequences);
}

MetaSequence *cactusDisk_getMetaSequenceInMemory(CactusDisk *cactusDisk,
        Name metaSequenceName) {
    static MetaSequence metaSequence;
    metaSequence.name = metaSequenceName;
    return stSortedSet_search(cactusDisk->metaSequences, &metaSequence);
}

MetaSequence *cactusDisk_getMetaSequence(CactusDisk *cactusDisk,
        Name metaSequenceName) {
    return cactusDisk_getObject(
            cactusDisk,
            cactusDisk->metaDataDatabase,
            (void *(*)(CactusDisk *, Name)) cactusDisk_getMetaSequenceInMemory,
            (void *(*)(void **, CactusDisk *)) metaSequence_loadFromBinaryRepresentation,
            metaSequenceName);
}

/*
 * Functions on meta events.
 */

void cactusDisk_addMetaEvent(CactusDisk *cactusDisk, MetaEvent *metaEvent) {
    assert(stSortedSet_search(cactusDisk->metaEvents, metaEvent) == NULL);
    stSortedSet_insert(cactusDisk->metaEvents, metaEvent);
}

void cactusDisk_deleteMetaEventFromDisk(CactusDisk *cactusDisk,
        Name metaEventName) {
    exitOnFailure(
            database_removeRecord(cactusDisk->metaDataDatabase, metaEventName),
            "Failed to remove the meta event: %s from the flower disk database\n",
            cactusMisc_nameToStringStatic(metaEventName));
}

void cactusDisk_unloadMetaEvent(CactusDisk *cactusDisk, MetaEvent *metaEvent) {
    stSortedSet_remove(cactusDisk->metaEvents, metaEvent);
}

MetaEvent *cactusDisk_getFirstMetaEventInMemory(CactusDisk *cactusDisk) {
    return stSortedSet_getFirst(cactusDisk->metaEvents);
}

MetaEvent *cactusDisk_getMetaEventInMemory(CactusDisk *cactusDisk,
        Name metaEventName) {
    static MetaEvent metaEvent;
    metaEvent.name = metaEventName;
    return stSortedSet_search(cactusDisk->metaEvents, &metaEvent);
}

MetaEvent *cactusDisk_getMetaEvent(CactusDisk *cactusDisk, Name metaEventName) {
    return cactusDisk_getObject(
            cactusDisk,
            cactusDisk->metaDataDatabase,
            (void *(*)(CactusDisk *, Name)) cactusDisk_getMetaEventInMemory,
            (void *(*)(void **, CactusDisk *)) metaEvent_loadFromBinaryRepresentation,
            metaEventName);
}

/*
 * Functions on strings stored by the flower disk.
 */

int64_t cactusDisk_addString(CactusDisk *cactusDisk, const char *string,
        int32_t length) {
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

char *cactusDisk_getString(CactusDisk *cactusDisk, int64_t offset,
        int32_t start, int32_t length, int32_t strand) {
    FILE *fileHandle;
    char *cA;
    char *cA2;

    fileHandle = fopen(cactusDisk->stringFile, "r");
    fseek(fileHandle, offset + start, SEEK_SET);
    cA = st_malloc(sizeof(char) * (length + 1));
    fread(cA, sizeof(char), length, fileHandle);
    cA[length] = '\0';
    fclose(fileHandle);

    if (!strand) {
        cA2 = cactusMisc_reverseComplementString(cA);
        free(cA);
        return cA2;
    }
    return cA;
}

void cactusDisk_getBlockOfUniqueIDs(CactusDisk *cactusDisk) {
    void *vA;
    Name keyName;
    exitOnFailure(database_startTransaction(cactusDisk->iDDatabase),
            "Failed to start transaction to get a block of unique names\n");
    keyName = 0;
    vA = database_getRecord(cactusDisk->iDDatabase, keyName);
    if (vA == NULL) {
        cactusDisk->uniqueNumber = 0;
    } else {
        cactusDisk->uniqueNumber = *((Name *) vA);
    }
    cactusDisk->maxUniqueNumber = cactusDisk->uniqueNumber
            + CACTUSDISK_NAME_INCREMENT;
    exitOnFailure(database_writeRecord(cactusDisk->iDDatabase, keyName,
            &cactusDisk->maxUniqueNumber, sizeof(Name)),
            "Failed to update the new max unique name (%s) in the database\n",
            cactusMisc_nameToStringStatic(cactusDisk->maxUniqueNumber));
    exitOnFailure(database_commitTransaction(cactusDisk->iDDatabase),
            "Failed to commit a transaction to get a block of unique names\n");
    free(vA);
}

int64_t cactusDisk_getUniqueID(CactusDisk *cactusDisk) {
    assert(cactusDisk->uniqueNumber <= cactusDisk->maxUniqueNumber);
    if (cactusDisk->uniqueNumber == cactusDisk->maxUniqueNumber) {
        cactusDisk_getBlockOfUniqueIDs(cactusDisk);
    }
    return cactusDisk->uniqueNumber++;
}
