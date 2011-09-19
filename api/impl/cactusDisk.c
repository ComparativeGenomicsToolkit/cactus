/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"
#include <unistd.h>

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic flower disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

#define CACTUS_DISK_NAME_INCREMENT 16384
#define CACTUS_DISK_BUCKET_NUMBER 65536
#define CACTUS_DISK_PARAMETER_KEY -100000

const char *CACTUS_DISK_EXCEPTION_ID = "CACTUS_DISK_EXCEPTION_ID";

static int32_t cactusDisk_constructFlowersP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(flower_getName((Flower *) o1), flower_getName((Flower *) o2));
}

static int32_t cactusDisk_constructMetaSequencesP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(metaSequence_getName((MetaSequence *) o1), metaSequence_getName((MetaSequence *) o2));
}

static void cactusDisk_writeBinaryRepresentation(CactusDisk *cactusDisk,
        void(*writeFn)(const void * ptr, size_t size, size_t count)) {
    binaryRepresentation_writeElementType(CODE_CACTUS_DISK, writeFn);
    binaryRepresentation_writeBool(cactusDisk->storeSequencesInAFile, writeFn);
    if (cactusDisk->storeSequencesInAFile) {
        assert(cactusDisk->sequencesFileName != NULL);
        binaryRepresentation_writeString(cactusDisk->sequencesFileName, writeFn);
    }
    binaryRepresentation_writeElementType(CODE_CACTUS_DISK, writeFn);
}

static void cactusDisk_loadFromBinaryRepresentation(void **binaryString, CactusDisk *cactusDisk, stKVDatabaseConf *conf) {
    cactusDisk->sequencesFileHandle = NULL;
    cactusDisk->sequencesFileName = NULL;
    cactusDisk->absSequencesFileName = NULL;
    assert(binaryRepresentation_peekNextElementType(*binaryString) == CODE_CACTUS_DISK);
    binaryRepresentation_popNextElementType(binaryString);
    cactusDisk->storeSequencesInAFile = binaryRepresentation_getBool(binaryString);
    if (cactusDisk->storeSequencesInAFile) {
        cactusDisk->sequencesFileName = binaryRepresentation_getString(binaryString);
        if (stKVDatabaseConf_getDir(conf) == NULL) {
            stThrowNew(CACTUS_DISK_EXCEPTION_ID,
                    "The database conf does not contain a directory in which the sequence file is to be found!\n");
        }
        cactusDisk->absSequencesFileName = stString_print("%s/%s", stKVDatabaseConf_getDir(conf),
                cactusDisk->sequencesFileName);
    }
    assert(binaryRepresentation_peekNextElementType(*binaryString) == CODE_CACTUS_DISK);
    binaryRepresentation_popNextElementType(binaryString);
}

/*
 * The following two functions compress and decompress the data in the cactus disk..
 */

static void *compress(void *data, int32_t *dataSize) {
    //Compression
    int64_t compressedSize;
    void *data2 = stCompression_compress(data, *dataSize, &compressedSize, -1);
    free(data);
    *dataSize = compressedSize;
    return data2;
}

static void *decompress(void *data, int64_t *dataSize) {
    //Decompression
    int64_t uncompressedSize;
    void *data2 = stCompression_decompress(data, *dataSize, &uncompressedSize);
    free(data);
    *dataSize = uncompressedSize;
    return data2;
}

static void *getRecord(CactusDisk *cactusDisk, Name objectName, char *type) {
    void *cA = NULL;
    int64_t recordSize = 0;
    if (stCache_containsRecord(cactusDisk->cache, objectName, 0, INT64_MAX)) { //If we already have the record, we won't update it.
        cA = stCache_getRecord(cactusDisk->cache, objectName, 0, INT64_MAX, &recordSize);
    } else {
        stTry {
            cA = stKVDatabase_getRecord2(cactusDisk->database, objectName,
                    &recordSize);
        }
        stCatch(except)
        {
            stThrowNewCause(
                    except,
                    ST_KV_DATABASE_EXCEPTION_ID,
                    "An unknown database error occurred when getting a %s",
                    type);
        }
        stTryEnd;
        if (cA == NULL) {
            return NULL;
        }
        stCache_setRecord(cactusDisk->cache, objectName, 0, recordSize, cA); //Add the compressed record to the cache.
    }
    //Decompression
    void *cA2 = decompress(cA, &recordSize);
    return cA2;
}

static bool containsRecord(CactusDisk *cactusDisk, Name objectName) {
    return stCache_containsRecord(cactusDisk->cache, objectName, 0, INT64_MAX) || stKVDatabase_containsRecord(
            cactusDisk->database, objectName);
}

static CactusDisk *cactusDisk_constructPrivate(stKVDatabaseConf *conf, bool create, const char *sequencesFileName) {
    CactusDisk *cactusDisk;
    cactusDisk = st_malloc(sizeof(CactusDisk));

    //construct lists of in memory objects
    cactusDisk->metaSequences = stSortedSet_construct3(cactusDisk_constructMetaSequencesP, NULL);
    cactusDisk->flowers = stSortedSet_construct3(cactusDisk_constructFlowersP, NULL);
    cactusDisk->flowerNamesMarkedForDeletion
            = stSortedSet_construct3((int(*)(const void *, const void *)) strcmp, free);

    //Now open the database
    cactusDisk->database = stKVDatabase_construct(conf, create);
    cactusDisk->cache = stCache_construct();
    cactusDisk->stringCache = stCache_construct();

    //initialise the unique ids.
    int32_t seed = (time(NULL) << 16) | (getpid() & 65535); //Likely to be unique
    st_logDebug("The cactus disk is seeding the random number generator with the value %i\n", seed);
    st_randomSeed(seed);
    cactusDisk->uniqueNumber = 0;
    cactusDisk->maxUniqueNumber = 0;

    //Now load any stuff..
    if (containsRecord(cactusDisk, CACTUS_DISK_PARAMETER_KEY)) {
        if (create) {
            stThrowNew(CACTUS_DISK_EXCEPTION_ID, "Tried to create a cactus disk, but the cactus disk already exists");
        }
        if (sequencesFileName != NULL) {
            stThrowNew(CACTUS_DISK_EXCEPTION_ID,
                    "A sequences file name is specified, but the cactus disk is not being created");
        }
        void *record = getRecord(cactusDisk, CACTUS_DISK_PARAMETER_KEY, "cactus_disk parameters");
        void *record2 = record;
        cactusDisk_loadFromBinaryRepresentation(&record, cactusDisk, conf);
        free(record2);
    } else {
        assert(create);
        if (sequencesFileName == NULL) {
            cactusDisk->storeSequencesInAFile = 0;
            cactusDisk->sequencesFileName = NULL;
            cactusDisk->sequencesFileHandle = NULL;
            cactusDisk->absSequencesFileName = NULL;
        } else {
            if (stKVDatabaseConf_getDir(conf) == NULL) {
                stThrowNew(CACTUS_DISK_EXCEPTION_ID,
                        "The database conf does not contain a directory in which the sequence file is to be found!\n");
            }
            cactusDisk->storeSequencesInAFile = 1;
            cactusDisk->sequencesFileName = stString_copy(sequencesFileName);
            cactusDisk->absSequencesFileName = stString_print("%s/%s", stKVDatabaseConf_getDir(conf),
                    cactusDisk->sequencesFileName);
            cactusDisk->sequencesFileHandle = fopen(cactusDisk->absSequencesFileName, "w");
            assert(cactusDisk->sequencesFileHandle != NULL);
            fclose(cactusDisk->sequencesFileHandle); //Flush it first time.
            cactusDisk->sequencesFileHandle = NULL;
        }
    }

    return cactusDisk;
}

CactusDisk *cactusDisk_construct(stKVDatabaseConf *conf, bool create) {
    return cactusDisk_constructPrivate(conf, create, NULL);
}

CactusDisk *cactusDisk_construct2(stKVDatabaseConf *conf, const char *sequencesFileName) {
    return cactusDisk_constructPrivate(conf, 1, sequencesFileName);
}

void cactusDisk_destruct(CactusDisk *cactusDisk) {
    Flower *flower;
    MetaSequence *metaSequence;

    while ((flower = stSortedSet_getFirst(cactusDisk->flowers)) != NULL) {
        flower_destruct(flower, FALSE);
    }
    stSortedSet_destruct(cactusDisk->flowers);

    stSortedSet_destruct(cactusDisk->flowerNamesMarkedForDeletion);

    while ((metaSequence = stSortedSet_getFirst(cactusDisk->metaSequences)) != NULL) {
        metaSequence_destruct(metaSequence);
    }
    stSortedSet_destruct(cactusDisk->metaSequences);

    //close DB
    stKVDatabase_destruct(cactusDisk->database);

    //Close the sequences files.
    if (cactusDisk->storeSequencesInAFile) {
        free(cactusDisk->sequencesFileName);
        free(cactusDisk->absSequencesFileName);
        if (cactusDisk->sequencesFileHandle != NULL) {
            fclose(cactusDisk->sequencesFileHandle);
        }
    }
#ifdef BEN_DEBUG
    else {
        assert(cactusDisk->sequencesFileName == NULL);
        assert(cactusDisk->sequencesFileHandle == NULL);
        assert(cactusDisk->absSequencesFileName == NULL);
    }
#endif

    stCache_destruct(cactusDisk->cache); //Get rid of the cache
    stCache_destruct(cactusDisk->stringCache);
    free(cactusDisk);
}

void cactusDisk_write(CactusDisk *cactusDisk) {
    Flower *flower;
    int32_t recordSize;

    stList *updateRequests = stList_construct3(0, (void(*)(void *)) stKVDatabaseBulkRequest_destruct);
    stList *removeRequests = stList_construct3(0, (void(*)(void *)) stInt64Tuple_destruct);

    st_logDebug("Starting to write the cactus to disk\n");

    stSortedSetIterator *it = stSortedSet_getIterator(cactusDisk->flowers);
    //Sort flowers to update.
    while ((flower = stSortedSet_getNext(it)) != NULL) {
        void
                *vA =
                        binaryRepresentation_makeBinaryRepresentation(
                                flower,
                                (void(*)(void *, void(*)(const void * ptr, size_t size, size_t count))) flower_writeBinaryRepresentation,
                                &recordSize);
        //Compression
        vA = compress(vA, &recordSize);
        if (containsRecord(cactusDisk, flower_getName(flower))) {
            int64_t recordSize2;
            void *vA2 = stCache_getRecord(cactusDisk->cache, flower_getName(flower), 0, INT64_MAX, &recordSize2);
            if (!stCache_recordsIdentical(vA, recordSize, vA2, recordSize2)) { //Only rewrite if we actually did something
                stList_append(updateRequests,
                        stKVDatabaseBulkRequest_constructUpdateRequest(flower_getName(flower), vA, recordSize));
            }
            free(vA2);
        } else {
            stList_append(updateRequests,
                    stKVDatabaseBulkRequest_constructInsertRequest(flower_getName(flower), vA, recordSize));
        }
        free(vA);
    }
    stSortedSet_destructIterator(it);

    st_logDebug("Got the flowers to update\n");

    //Remove nets that are marked for deletion..
    it = stSortedSet_getIterator(cactusDisk->flowerNamesMarkedForDeletion);
    char *nameString;
    while ((nameString = stSortedSet_getNext(it)) != NULL) {
        Name name = cactusMisc_stringToName(nameString);
        if (containsRecord(cactusDisk, name)) {
            stList_append(updateRequests, stKVDatabaseBulkRequest_constructUpdateRequest(name, &name, 0)); //We set it to null in the first atomic operation.
            stList_append(removeRequests, stInt64Tuple_construct(1, name));
        }
    }
    stSortedSet_destructIterator(it);

    st_logDebug("Avoided updating nets marked for deletion\n");

    //Sort sequences to insert
    it = stSortedSet_getIterator(cactusDisk->metaSequences);
    MetaSequence *metaSequence;
    while ((metaSequence = stSortedSet_getNext(it)) != NULL) {
        if (!containsRecord(cactusDisk, metaSequence_getName(metaSequence))) { //We do not edit meta sequences, so we do not update it..
            void
                    *vA =
                            binaryRepresentation_makeBinaryRepresentation(
                                    metaSequence,
                                    (void(*)(void *, void(*)(const void * ptr, size_t size, size_t count))) metaSequence_writeBinaryRepresentation,
                                    &recordSize);
            //Compression
            vA = compress(vA, &recordSize);
            stList_append(updateRequests,
                    stKVDatabaseBulkRequest_constructInsertRequest(metaSequence_getName(metaSequence), vA, recordSize));
            free(vA);
        }
    }
    stSortedSet_destructIterator(it);

    st_logDebug("Got the sequences we are going to add to the database.");

    if (!containsRecord(cactusDisk, CACTUS_DISK_PARAMETER_KEY)) { //We only write the parameters once.
        //Finally the database info.
        void
                *cactusDiskParameters =
                        binaryRepresentation_makeBinaryRepresentation(
                                cactusDisk,
                                (void(*)(void *, void(*)(const void * ptr, size_t size, size_t count))) cactusDisk_writeBinaryRepresentation,
                                &recordSize);
        //Compression
        cactusDiskParameters = compress(cactusDiskParameters, &recordSize);
        stList_append(
                updateRequests,
                stKVDatabaseBulkRequest_constructInsertRequest(CACTUS_DISK_PARAMETER_KEY, cactusDiskParameters,
                        recordSize));
        free(cactusDiskParameters);
    }

    st_logDebug("Checked if need to write the initial parameters\n");

    if (stList_length(updateRequests) > 0) {
        stTry {
            int maxRecords = 100000;
            while(stList_length(updateRequests) > maxRecords) {
                    stList *list = stList_construct3(0, (void(*)(void *)) stKVDatabaseBulkRequest_destruct);
                    while(stList_length(list) < maxRecords) {
                        stList_append(list, stList_pop(updateRequests));
                    }
                    stKVDatabase_bulkSetRecords(cactusDisk->database,
                                        list);
                    stList_destruct(list);
                }

            stKVDatabase_bulkSetRecords(cactusDisk->database,
                    updateRequests);
        }
        stCatch(except)
        {
            stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                    "Failed when trying to set records in updating the cactus disk");
        }stTryEnd;
    }

    st_logDebug("Updated the database with inserts\n");

    if (stList_length(removeRequests) > 0) {
        stTry {
            stKVDatabase_bulkRemoveRecords(cactusDisk->database,
                    removeRequests);
        }
        stCatch(except)
        {
            stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                    "Failed when trying to remove records in updating the cactus disk");
        }stTryEnd;
    }

    st_logDebug("Now removed flowers we don't need\n");

    stList_destruct(updateRequests);
    stList_destruct(removeRequests);

    st_logDebug("Finished writing to the database\n");
}

static Cap *getNextStub(Cap *cap) {
    Cap *cap2 = cap_getAdjacency(cap);
    Cap *cap3 = cap_getOtherSegmentCap(cap2);
    return cap3 == NULL ? cap2 : getNextStub(cap3);
}

Flower *cactusDisk_getFlower(CactusDisk *cactusDisk, Name flowerName) {
    static Flower flower;
    flower.name = flowerName;
    Flower *flower2;
    if ((flower2 = stSortedSet_search(cactusDisk->flowers, &flower)) != NULL) {
        return flower2;
    }
    void *cA = getRecord(cactusDisk, flowerName, "flower");

    if (cA == NULL) {
        return NULL;
    }
    void *cA2 = cA;
    flower2 = flower_loadFromBinaryRepresentation(&cA2, cactusDisk);
    free(cA);
    return flower2;
}

MetaSequence *cactusDisk_getMetaSequence(CactusDisk *cactusDisk, Name metaSequenceName) {
    static MetaSequence metaSequence;
    metaSequence.name = metaSequenceName;
    MetaSequence *metaSequence2;
    if ((metaSequence2 = stSortedSet_search(cactusDisk->metaSequences, &metaSequence)) != NULL) {
        return metaSequence2;
    }
    void *cA = getRecord(cactusDisk, metaSequenceName, "metaSequence");
    if (cA == NULL) {
        return NULL;
    }
    void *cA2 = cA;
    metaSequence2 = metaSequence_loadFromBinaryRepresentation(&cA2, cactusDisk);
    free(cA);
    return metaSequence2;
}

/*
 * Private functions.
 */

bool cactusDisk_flowerIsLoaded(CactusDisk *cactusDisk, Name flowerName) {
    static Flower flower;
    flower.name = flowerName;
    return stSortedSet_search(cactusDisk->flowers, &flower) != NULL;
}

void cactusDisk_addFlower(CactusDisk *cactusDisk, Flower *flower) {
    assert(stSortedSet_search(cactusDisk->flowers, flower) == NULL);
    stSortedSet_insert(cactusDisk->flowers, flower);
}

void cactusDisk_removeFlower(CactusDisk *cactusDisk, Flower *flower) {
    assert(cactusDisk_flowerIsLoaded(cactusDisk, flower_getName(flower)));
    stSortedSet_remove(cactusDisk->flowers, flower);
}

void cactusDisk_deleteFlowerFromDisk(CactusDisk *cactusDisk, Flower *flower) {
    char *nameString = cactusMisc_nameToString(flower_getName(flower));
    if (stSortedSet_search(cactusDisk->flowerNamesMarkedForDeletion, nameString) == NULL) {
        stSortedSet_insert(cactusDisk->flowerNamesMarkedForDeletion, nameString);
    } else {
        free(nameString);
    }
}

/*
 * Functions on meta sequences.
 */

void cactusDisk_addMetaSequence(CactusDisk *cactusDisk, MetaSequence *metaSequence) {
    assert(stSortedSet_search(cactusDisk->metaSequences, metaSequence) == NULL);
    stSortedSet_insert(cactusDisk->metaSequences, metaSequence);
}

void cactusDisk_removeMetaSequence(CactusDisk *cactusDisk, MetaSequence *metaSequence) {
    assert(stSortedSet_search(cactusDisk->metaSequences, metaSequence) != NULL);
    stSortedSet_remove(cactusDisk->metaSequences, metaSequence);
}

/*
 * Functions on strings stored by the flower disk.
 */

Name cactusDisk_addString(CactusDisk *cactusDisk, const char *string) {
    Name name;
    //bool done = 0;
    if (cactusDisk->storeSequencesInAFile) {
        if(cactusDisk->sequencesFileHandle != NULL) {
            fclose(cactusDisk->sequencesFileHandle);
        }
        cactusDisk->sequencesFileHandle = fopen(cactusDisk->absSequencesFileName, "a");
        assert(cactusDisk->sequencesFileHandle != NULL);
        name = ftell(cactusDisk->sequencesFileHandle) + 1;
        fprintf(cactusDisk->sequencesFileHandle, ">%s", string);
        fclose(cactusDisk->sequencesFileHandle);
        cactusDisk->sequencesFileHandle = NULL;
        return name;
    } else {
        name = cactusDisk_getUniqueID(cactusDisk);
        stTry {
            stKVDatabase_insertRecord(cactusDisk->database, name, string,
                    (strlen(string) + 1) * sizeof(char));
        }
        stCatch(except)
        {
            stThrowNewCause(
                    except,
                    ST_KV_DATABASE_EXCEPTION_ID,
                    "An unknown database error occurred when we tried to add a string to the cactus disk");
        }stTryEnd;
    }
    return name;
}

char *cactusDisk_getString(CactusDisk *cactusDisk, Name name, int32_t start, int32_t length, int32_t strand,
        int32_t totalSequenceLength) {
    assert(length >= 0);
    if (length == 0) {
        return stString_copy("");
    }
    //bool done = 0;
    char *string = NULL;

    //First try getting it from the cache
    if (stCache_containsRecord(cactusDisk->stringCache, name, start, sizeof(char) * length)) {
        int64_t recordSize;
        string = stCache_getRecord(cactusDisk->stringCache, name, start, sizeof(char) * length, &recordSize);
        assert(string != NULL);
        assert(recordSize == length);
        string = realloc(string, length + 1);
        assert(string != NULL);
    } else {
        if (cactusDisk->storeSequencesInAFile) {
            if (cactusDisk->sequencesFileHandle == NULL) {
                cactusDisk->sequencesFileHandle = fopen(cactusDisk->absSequencesFileName, "r");
                assert(cactusDisk->sequencesFileHandle != NULL);
            }
            fseek(cactusDisk->sequencesFileHandle, name + start, SEEK_SET);
            string = st_malloc(sizeof(char) * (length + 1));
            fread(string, sizeof(char), length, cactusDisk->sequencesFileHandle);
#ifdef BEN_DEBUG
            for (int32_t i = 0; i < length; i++) {
                assert(string[i] != '>');
            }
#endif
        } else {
            stTry {
                string = stKVDatabase_getPartialRecord(
                        cactusDisk->database, name, start * sizeof(char),
                        (length + 1) * sizeof(char),
                        (totalSequenceLength + 1) * sizeof(char));
            }
            stCatch(except)
            {
                stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                        "An unknown database error occurred when getting a sequence string");
            }stTryEnd;
        }
        stCache_setRecord(cactusDisk->stringCache, name, start, sizeof(char) * length, string);
    }
    string[length] = '\0';
    if (!strand) {
        char *string2 = cactusMisc_reverseComplementString(string);
        free(string);
        return string2;
    }
    return string;
}

/*
 * Function to get unique ID.
 */

void cactusDisk_getBlockOfUniqueIDs(CactusDisk *cactusDisk) {
    bool done = 0;
    int32_t collisionCount = 0;
    while (!done) {
        stTry {
            Name keyName = st_randomInt(-CACTUS_DISK_BUCKET_NUMBER, 0);
            assert(keyName >= -CACTUS_DISK_BUCKET_NUMBER);
            assert(keyName < 0);
            int64_t bucketSize = INT64_MAX / CACTUS_DISK_BUCKET_NUMBER;
            int64_t minimumValue = bucketSize * (abs(keyName) - 1) + 1; //plus one for the reserved '0' value.
            int64_t maximumValue = minimumValue + (bucketSize - 1);
            assert(minimumValue >= 1);
            assert(maximumValue <= INT64_MAX);
            assert(minimumValue < maximumValue);
            if (stKVDatabase_containsRecord(cactusDisk->database, keyName)) {
                cactusDisk->maxUniqueNumber = stKVDatabase_incrementInt64(
                        cactusDisk->database, keyName,
                        CACTUS_DISK_NAME_INCREMENT);
                cactusDisk->uniqueNumber = cactusDisk->maxUniqueNumber
                - CACTUS_DISK_NAME_INCREMENT;
                if(cactusDisk->uniqueNumber <= 0 || cactusDisk->uniqueNumber < minimumValue || cactusDisk->uniqueNumber > maximumValue) {
                    st_errAbort("Got a non positive unique number %lli %lli %lli %lli", cactusDisk->uniqueNumber, cactusDisk->maxUniqueNumber, minimumValue, maximumValue);
                }
                assert(cactusDisk->uniqueNumber >= minimumValue);
                assert(cactusDisk->uniqueNumber <= maximumValue);
                assert(cactusDisk->uniqueNumber > 0);
            } else {
                stTry {
                    stKVDatabase_insertInt64(cactusDisk->database,
                            keyName, minimumValue);
                }
                stCatch(except)
                {
                    collisionCount++;
                    if (collisionCount >= 10) {
                        stThrowNewCause(
                                except,
                                ST_KV_DATABASE_EXCEPTION_ID,
                                "Repeated unknown database errors occurred when we tried to get a unique ID, collision count %i",
                                collisionCount);
                    } else {
                        st_logDebug(
                                "Got an exception when trying to insert a record: %s",
                                except);
                    }
                }stTryEnd;
                continue;
            }
            if (cactusDisk->maxUniqueNumber >= maximumValue) {
                st_errAbort(
                        "We have exhausted a bucket, which seems really unlikely");
            }
            done = 1;
        }
        stCatch(except) {
            collisionCount++;
            if (collisionCount >= 10) {
                stThrowNewCause(
                        except,
                        ST_KV_DATABASE_EXCEPTION_ID,
                        "Repeated unknown database errors occurred when we tried to get a unique ID, collision count %i",
                        collisionCount);
            } else {
                st_logDebug(
                        "Got an exception when trying to insert a record: %s",
                        except);
            }
        }stTryEnd;
    }
}

int64_t cactusDisk_getUniqueID(CactusDisk *cactusDisk) {
    assert(cactusDisk->uniqueNumber <= cactusDisk->maxUniqueNumber);
    if (cactusDisk->uniqueNumber == cactusDisk->maxUniqueNumber) {
        cactusDisk_getBlockOfUniqueIDs(cactusDisk);
    }
    return cactusDisk->uniqueNumber++;
}

