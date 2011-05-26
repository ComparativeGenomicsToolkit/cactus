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
    return cactusMisc_nameCompare(flower_getName((Flower *) o1),
            flower_getName((Flower *) o2));
}

static int32_t cactusDisk_constructMetaSequencesP(const void *o1,
        const void *o2) {
    return cactusMisc_nameCompare(metaSequence_getName((MetaSequence *) o1),
            metaSequence_getName((MetaSequence *) o2));
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

static void cactusDisk_loadFromBinaryRepresentation(void **binaryString,
        CactusDisk *cactusDisk) {
    cactusDisk->sequencesFileHandle = NULL;
    cactusDisk->sequencesFileName = NULL;
    cactusDisk->absSequencesFileName = NULL;
    assert(binaryRepresentation_peekNextElementType(*binaryString)
            == CODE_CACTUS_DISK);
    binaryRepresentation_popNextElementType(binaryString);
    cactusDisk->storeSequencesInAFile = binaryRepresentation_getBool(
            binaryString);
    if (cactusDisk->storeSequencesInAFile) {
        cactusDisk->sequencesFileName = binaryRepresentation_getString(
                binaryString);
    }
    assert(binaryRepresentation_peekNextElementType(*binaryString)
            == CODE_CACTUS_DISK);

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
    //bool done = 0;
    void *cA = NULL;
    int64_t recordSize = 0;
    //while (!done) {
    stTry {
        cA = stKVDatabase_getRecord2(cactusDisk->database, objectName, &recordSize);
    }
    stCatch(except)
    {
        stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                "An unknown database error occurred when getting a %s", type);
    }
    stTryEnd;
    //}
    if (cA == NULL) {
        return NULL;
    }
    //Decompression
    void *cA2 = decompress(cA, &recordSize);
    return cA2;
}

static CactusDisk *cactusDisk_constructPrivate(stKVDatabaseConf *conf,
        bool create, bool preCacheSequences, const char *sequencesFileName) {
    CactusDisk *cactusDisk;
    cactusDisk = st_malloc(sizeof(CactusDisk));

    //construct lists of in memory objects
    cactusDisk->metaSequences = stSortedSet_construct3(
            cactusDisk_constructMetaSequencesP, NULL);
    cactusDisk->flowers = stSortedSet_construct3(cactusDisk_constructFlowersP,
            NULL);
    cactusDisk->flowerNamesMarkedForDeletion = stSortedSet_construct3((int(*)(
            const void *, const void *)) strcmp, free);

    //Now open the database
    //preCacheSequences = 0;
    cactusDisk->preCacheSequences = preCacheSequences;
    cactusDisk->database = stKVDatabase_construct(conf, create);
    stKVDatabase_makeMemCache(cactusDisk->database, 0, preCacheSequences ? 0
            : 5000);

    //initialise the unique ids.
    //cactusDisk_getUniqueID(cactusDisk);
    int32_t seed = (time(NULL) << 16) | (getpid() & 65535); //Likely to be unique
    st_logDebug(
            "The cactus disk is seeding the random number generator with the value %i\n",
            seed);
    st_randomSeed(seed);
    cactusDisk->uniqueNumber = 0;
    cactusDisk->maxUniqueNumber = 0;

    //Now load any stuff..
    if (stKVDatabase_containsRecord(cactusDisk->database,
            CACTUS_DISK_PARAMETER_KEY)) {
        if (create) {
            stThrowNew(CACTUS_DISK_EXCEPTION_ID,
                    "Tried to create a cactus disk, but the cactus disk already exists");
        }
        if (sequencesFileName != NULL) {
            stThrowNew(CACTUS_DISK_EXCEPTION_ID,
                    "A sequences file name is specified, but the cactus disk is not being created");
        }
        void *record = getRecord(cactusDisk, CACTUS_DISK_PARAMETER_KEY,
                "cactus_disk parameters");
        void *record2 = record;
        cactusDisk_loadFromBinaryRepresentation(&record, cactusDisk);
        free(record2);
        if (cactusDisk->storeSequencesInAFile) {
            assert(cactusDisk->sequencesFileName != NULL);
            if (stKVDatabaseConf_getDir(conf) == NULL) {
                stThrowNew(
                        CACTUS_DISK_EXCEPTION_ID,
                        "The database conf does not contain a directory in which the sequence file is to be found!\n");
            }

            cactusDisk->absSequencesFileName = stString_print("%s/%s",
                    stKVDatabaseConf_getDir(conf),
                    cactusDisk->sequencesFileName);
            cactusDisk->sequencesFileHandle = fopen(
                    cactusDisk->absSequencesFileName, "r");
            assert(cactusDisk->sequencesFileHandle != NULL);
        }
    } else {
        assert(create);
        if (sequencesFileName == NULL) {
            cactusDisk->storeSequencesInAFile = 0;
            cactusDisk->sequencesFileName = NULL;
            cactusDisk->sequencesFileHandle = NULL;
        } else {
            if (stKVDatabaseConf_getDir(conf) == NULL) {
                stThrowNew(
                        CACTUS_DISK_EXCEPTION_ID,
                        "The database conf does not contain a directory in which the sequence file is to be found!\n");
            }
            cactusDisk->storeSequencesInAFile = 1;
            cactusDisk->sequencesFileName = stString_copy(sequencesFileName);
            cactusDisk->absSequencesFileName = stString_print("%s/%s",
                    stKVDatabaseConf_getDir(conf),
                    cactusDisk->sequencesFileName);
            cactusDisk->sequencesFileHandle = fopen(cactusDisk->absSequencesFileName, "w");
            assert(cactusDisk->sequencesFileHandle != NULL);
            fclose(cactusDisk->sequencesFileHandle); //Flush it first time.
            cactusDisk->sequencesFileHandle = fopen(
                    cactusDisk->absSequencesFileName, "r");
            assert(cactusDisk->sequencesFileHandle != NULL);
        }
    }

    return cactusDisk;
}

CactusDisk *cactusDisk_construct2(stKVDatabaseConf *conf, bool create,
        bool preCacheSequences) {
    return cactusDisk_constructPrivate(conf, create, preCacheSequences, NULL);
}

CactusDisk *cactusDisk_construct(stKVDatabaseConf *conf, bool create) {
    return cactusDisk_constructPrivate(conf, create, 0, NULL);
}

CactusDisk *cactusDisk_construct3(stKVDatabaseConf *conf,
        bool preCacheSequences, const char *sequencesFileName) {
    return cactusDisk_constructPrivate(conf, 1, preCacheSequences,
            sequencesFileName);
}

void cactusDisk_destruct(CactusDisk *cactusDisk) {
    Flower *flower;
    MetaSequence *metaSequence;

    while ((flower = stSortedSet_getFirst(cactusDisk->flowers)) != NULL) {
        flower_destruct(flower, FALSE);
    }
    stSortedSet_destruct(cactusDisk->flowers);

    stSortedSet_destruct(cactusDisk->flowerNamesMarkedForDeletion);

    while ((metaSequence = stSortedSet_getFirst(cactusDisk->metaSequences))
            != NULL) {
        metaSequence_destruct(metaSequence);
    }
    stSortedSet_destruct(cactusDisk->metaSequences);

    //close DB
    stKVDatabase_destruct(cactusDisk->database);

    //Close the sequences files.
    if (cactusDisk->storeSequencesInAFile) {
        free(cactusDisk->sequencesFileName);
        free(cactusDisk->absSequencesFileName);
        fclose(cactusDisk->sequencesFileHandle);
    }
#ifdef BEN_DEBUG
    else {
        assert(cactusDisk->sequencesFileName == NULL);
        assert(cactusDisk->sequencesFileHandle == NULL);
        assert(cactusDisk->absSequencesFileName == NULL);
    }
#endif

    free(cactusDisk);
}

static void cactusDisk_writeP(stList *list) {
    for (int32_t i = 0; i < stList_length(list); i += 3) {
        free(stList_get(list, i + 1));
        stIntTuple_destruct(stList_get(list, i + 2));
    }
    stList_destruct(list);
}

void cactusDisk_write(CactusDisk *cactusDisk) {
    stList *flowersToUpdate = stList_construct();
    stList *flowersToInsert = stList_construct();
    stList *flowersToRemove = stList_construct();
    stList *sequencesToInsert = stList_construct();

    stSortedSetIterator *it = stSortedSet_getIterator(cactusDisk->flowers);
    Flower *flower;
    int32_t recordSize;

    while ((flower = stSortedSet_getNext(it)) != NULL) {
        void *vA = binaryRepresentation_makeBinaryRepresentation(flower,
                (void(*)(void *, void(*)(const void * ptr, size_t size,
                        size_t count))) flower_writeBinaryRepresentation,
                &recordSize);
        //Compression
        vA = compress(vA, &recordSize);
        if (stKVDatabase_containsRecord(cactusDisk->database, flower_getName(
                flower))) {
            stList_append(flowersToUpdate, flower);
            stList_append(flowersToUpdate, vA);
            stList_append(flowersToUpdate, stIntTuple_construct(1, recordSize));

        } else {
            stList_append(flowersToInsert, flower);
            stList_append(flowersToInsert, vA);
            stList_append(flowersToInsert, stIntTuple_construct(1, recordSize));
        }
    }
    stSortedSet_destructIterator(it);

    it = stSortedSet_getIterator(cactusDisk->metaSequences);
    MetaSequence *metaSequence;
    while ((metaSequence = stSortedSet_getNext(it)) != NULL) {
        if (!stKVDatabase_containsRecord(cactusDisk->database,
                metaSequence_getName(metaSequence))) { //We do not edit meta sequences, so we do not update it..
            void
                    *vA =
                            binaryRepresentation_makeBinaryRepresentation(
                                    metaSequence,
                                    (void(*)(void *, void(*)(const void * ptr,
                                            size_t size, size_t count))) metaSequence_writeBinaryRepresentation,
                                    &recordSize);
            //Compression
            vA = compress(vA, &recordSize);
            stList_append(sequencesToInsert, metaSequence);
            stList_append(sequencesToInsert, vA);
            stList_append(sequencesToInsert,
                    stIntTuple_construct(1, recordSize));
        }
    }
    stSortedSet_destructIterator(it);

    //Remove nets that are marked for deletion..
    it = stSortedSet_getIterator(cactusDisk->flowerNamesMarkedForDeletion);
    char *nameString;
    while ((nameString = stSortedSet_getNext(it)) != NULL) {
        Name name = cactusMisc_stringToName(nameString);
        if (stKVDatabase_containsRecord(cactusDisk->database, name)) {
            stList_append(flowersToRemove, nameString);
        }
    }
    stSortedSet_destructIterator(it);

    //Finally the database info.
    int32_t cactusDiskParametersRecordSize;
    void *cactusDiskParameters = binaryRepresentation_makeBinaryRepresentation(
            cactusDisk, (void(*)(void *, void(*)(const void * ptr, size_t size,
                    size_t count))) cactusDisk_writeBinaryRepresentation,
            &cactusDiskParametersRecordSize);
    //Compression
    cactusDiskParameters = compress(cactusDiskParameters,
            &cactusDiskParametersRecordSize);
    bool containsCactusDiskParameters = stKVDatabase_containsRecord(
            cactusDisk->database, CACTUS_DISK_PARAMETER_KEY);

    stKVDatabase_startTransaction(cactusDisk->database);
    stTry {
        if(containsCactusDiskParameters) {
            stKVDatabase_updateRecord(cactusDisk->database, CACTUS_DISK_PARAMETER_KEY, cactusDiskParameters, cactusDiskParametersRecordSize);
        }
        else {
            stKVDatabase_insertRecord(cactusDisk->database, CACTUS_DISK_PARAMETER_KEY, cactusDiskParameters, cactusDiskParametersRecordSize);
        }

        for (int32_t i = 0; i < stList_length(flowersToInsert); i += 3) {
            Flower *flower = stList_get(flowersToInsert, i);
            void *vA = stList_get(flowersToInsert, i + 1);
            stIntTuple *recordSize = stList_get(flowersToInsert, i + 2);
            stKVDatabase_insertRecord(cactusDisk->database, flower_getName(flower), vA, stIntTuple_getPosition(
                            recordSize, 0));
        }

        for (int32_t i = 0; i < stList_length(flowersToUpdate); i += 3) {
            Flower *flower = stList_get(flowersToUpdate, i);
            void *vA = stList_get(flowersToUpdate, i + 1);
            stIntTuple *recordSize = stList_get(flowersToUpdate, i + 2);
            stKVDatabase_updateRecord(cactusDisk->database, flower_getName(flower), vA, stIntTuple_getPosition(
                            recordSize, 0));
        }

        for (int32_t i = 0; i < stList_length(flowersToRemove); i++) {
            Name name = cactusMisc_stringToName(stList_get(flowersToRemove, i));
            stKVDatabase_removeRecord(cactusDisk->database, name);
        }

        for (int32_t i = 0; i < stList_length(sequencesToInsert); i += 3) {
            metaSequence = stList_get(sequencesToInsert, i);
            void *vA = stList_get(sequencesToInsert, i + 1);
            stIntTuple *recordSize = stList_get(sequencesToInsert, i + 2);
            stKVDatabase_insertRecord(cactusDisk->database, metaSequence_getName(metaSequence), vA,
                    stIntTuple_getPosition(recordSize, 0));
        }

        stKVDatabase_commitTransaction(cactusDisk->database);
    }
    stCatch(except)
    {
        stKVDatabase_abortTransaction(cactusDisk->database);
        stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                "An unknown database error occurred when updating flowers and metasequences on the cactus disk");
    }
    stTryEnd;

    cactusDisk_writeP(flowersToUpdate);
    cactusDisk_writeP(flowersToInsert);
    stList_destruct(flowersToRemove);
    cactusDisk_writeP(sequencesToInsert);
    free(cactusDiskParameters);
}

static Cap *getNextStub(Cap *cap) {
    Cap *cap2 = cap_getAdjacency(cap);
    Cap *cap3 = cap_getOtherSegmentCap(cap2);
    return cap3 == NULL ? cap2 : getNextStub(cap3);
}

static void preCacheSequences(CactusDisk *cactusDisk, Flower *flower) {
    /*
     * Gets all the sequences in the flower up front (works with the c
     */
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end)) {
            End_InstanceIterator *instanceIterator = end_getInstanceIterator(
                    end);
            Cap *cap, *cap2;
            while ((cap = end_getNext(instanceIterator)) != NULL) {
                Sequence *sequence = cap_getSequence(cap);
                cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                if (sequence != NULL && !cap_getSide(cap)) {
                    cap2 = getNextStub(cap);
#ifdef BEN_DEBUG
                    assert(end_isStubEnd(cap_getEnd(cap2)));
                    assert(cap2 != NULL);
                    assert(cap_getStrand(cap2));
                    assert(sequence == cap_getSequence(cap2));
                    assert(cap_getSide(cap2));
#endif
                    int32_t start = cap_getCoordinate(cap) + 1;
                    int32_t stop = cap_getCoordinate(cap2);
                    int32_t length = stop - start;
#ifdef BEN_DEBUG
                    assert(length >= 0);
#endif
                    if (length > 0) {
                        char *cA = sequence_getString(sequence, start, length,
                                1); //Strand doesn't matter as we store and access the strings in the same orientation
                        free(cA);
                    }
                }
            }
            end_destructInstanceIterator(instanceIterator);
        }
    }
    flower_destructEndIterator(endIt);
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
    if (cactusDisk->preCacheSequences) {
        preCacheSequences(cactusDisk, flower2);
    }
    return flower2;
}

MetaSequence *cactusDisk_getMetaSequence(CactusDisk *cactusDisk,
        Name metaSequenceName) {
    static MetaSequence metaSequence;
    metaSequence.name = metaSequenceName;
    MetaSequence *metaSequence2;
    if ((metaSequence2 = stSortedSet_search(cactusDisk->metaSequences,
            &metaSequence)) != NULL) {
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
    if (stSortedSet_search(cactusDisk->flowerNamesMarkedForDeletion, nameString)
            == NULL) {
        stSortedSet_insert(cactusDisk->flowerNamesMarkedForDeletion, nameString);
    } else {
        free(nameString);
    }
}

/*
 * Functions on meta sequences.
 */

void cactusDisk_addMetaSequence(CactusDisk *cactusDisk,
        MetaSequence *metaSequence) {
    assert(stSortedSet_search(cactusDisk->metaSequences, metaSequence) == NULL);
    stSortedSet_insert(cactusDisk->metaSequences, metaSequence);
}

void cactusDisk_removeMetaSequence(CactusDisk *cactusDisk,
        MetaSequence *metaSequence) {
    assert(stSortedSet_search(cactusDisk->metaSequences, metaSequence) != NULL);
    stSortedSet_remove(cactusDisk->metaSequences, metaSequence);
}

/*
 * Functions on strings stored by the flower disk.
 */

Name cactusDisk_addString(CactusDisk *cactusDisk, const char *string) {
    Name name;
    bool done = 0;
    if (cactusDisk->storeSequencesInAFile) {
        fclose(cactusDisk->sequencesFileHandle);
        cactusDisk->sequencesFileHandle = fopen(cactusDisk->absSequencesFileName, "a");
        assert(cactusDisk->sequencesFileHandle != NULL);
        name = ftell(cactusDisk->sequencesFileHandle) + 1;
        fprintf(cactusDisk->sequencesFileHandle, ">%s", string);
        fclose(cactusDisk->sequencesFileHandle);
        cactusDisk->sequencesFileHandle = fopen(cactusDisk->absSequencesFileName, "r");
        assert(cactusDisk->sequencesFileHandle != NULL);
        return name;
    } else {
        name = cactusDisk_getUniqueID(cactusDisk);
        while (!done) {
            stKVDatabase_startTransaction(cactusDisk->database);
            stTry {
                stKVDatabase_insertRecord(cactusDisk->database, name, string, (strlen(string) + 1) * sizeof(char));
                stKVDatabase_commitTransaction(cactusDisk->database);
                done = 1;
            }
            stCatch(except)
            {
                if (stExcept_getId(except) == ST_KV_DATABASE_RETRY_TRANSACTION_EXCEPTION_ID) {
                    st_logDebug(
                            "We have caught a retry transaction exception when adding a string to the database, we will try again\n");
                    stExcept_free(except);
                    stKVDatabase_abortTransaction(cactusDisk->database);
                } else {
                    stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                            "An unknown database error occurred when we tried to add a string to the cactus disk");
                }
            }stTryEnd;
        }
    }
    return name;
}

char *cactusDisk_getString(CactusDisk *cactusDisk, Name name, int32_t start,
        int32_t length, int32_t strand, int32_t totalSequenceLength) {
    //bool done = 0;
    char *string = NULL;
    if (cactusDisk->storeSequencesInAFile) {
        fseek(cactusDisk->sequencesFileHandle, name + start, SEEK_SET);
        string = st_malloc(sizeof(char) * (length + 1));
        fread(string, sizeof(char), length, cactusDisk->sequencesFileHandle);
#ifdef BEN_DEBUG
        for(int32_t i=0; i<length; i++) {
            assert(string[i] != '>');
        }
#endif
    } else {
        stTry {
            string = stKVDatabase_getPartialRecord(cactusDisk->database, name, start * sizeof(char), (length + 1)
                    * sizeof(char), (totalSequenceLength + 1) * sizeof(char));
        }
        stCatch(except)
        {
            stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                    "An unknown database error occurred when getting a sequence string");
        }stTryEnd;
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
        stKVDatabase_startTransaction(cactusDisk->database);
        stTry {
            Name keyName = st_randomInt(-CACTUS_DISK_BUCKET_NUMBER, 0);
            assert(keyName >= -CACTUS_DISK_BUCKET_NUMBER);
            assert(keyName < 0);
            void *vA = stKVDatabase_getRecord(cactusDisk->database, keyName);
            int64_t bucketSize = INT64_MAX / CACTUS_DISK_BUCKET_NUMBER;
            int64_t minimumValue = bucketSize * (abs(keyName) - 1) + 1; //plus one for the reserved '0' value.
            int64_t maximumValue = minimumValue + (bucketSize - 1);
            bool recordExists = 0;
            if (vA == NULL) {
                cactusDisk->uniqueNumber = minimumValue;
            } else {
                recordExists = 1;
                cactusDisk->uniqueNumber = *((Name *) vA);
                free(vA);
            }
            cactusDisk->maxUniqueNumber = cactusDisk->uniqueNumber + CACTUS_DISK_NAME_INCREMENT;
            if (cactusDisk->maxUniqueNumber >= maximumValue) {
                st_errAbort("We have exhausted a bucket, which seems really unlikely");
            }
            if (recordExists) {
                stKVDatabase_updateRecord(cactusDisk->database, keyName, &cactusDisk->maxUniqueNumber, sizeof(Name));
            } else {
                stKVDatabase_insertRecord(cactusDisk->database, keyName, &cactusDisk->maxUniqueNumber, sizeof(Name));
            }
            stKVDatabase_commitTransaction(cactusDisk->database);
            done = 1;
        }
        stCatch(except)
        {
            stKVDatabase_abortTransaction(cactusDisk->database);
            if (stExcept_getId(except) == ST_KV_DATABASE_RETRY_TRANSACTION_EXCEPTION_ID || collisionCount++
                    < 10) {
                st_logDebug(
                        "We have caught a retry transaction exception when allocating a new id, we will try again\n");
                stExcept_free(except);
            } else {
                stThrowNewCause(
                        except,
                        ST_KV_DATABASE_EXCEPTION_ID,
                        "An unknown database error occurred when we tried to get a unique ID, collision count %i",
                        collisionCount);
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

static void flowerIterator(Flower *flower, void(*fn)(Flower *, void *),
        void *extra) {
    fn(flower, extra);
    Group *group;
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_getNestedFlower(group) != NULL) {
            flowerIterator(group_getNestedFlower(group), fn, extra);
        }
    }
    flower_destructGroupIterator(groupIt);
}

static void copyAcrossFlowersP(Flower *flower, CactusDisk *newCactusDisk) {
    //Find the meta sequences associated with the flower..
    Flower_SequenceIterator *sequenceIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    int32_t recordSize;
    while ((sequence = flower_getNextSequence(sequenceIt)) != NULL) {
        MetaSequence *metaSequence = sequence_getMetaSequence(sequence);
        if (cactusDisk_getMetaSequence(newCactusDisk, metaSequence_getName(
                metaSequence)) == NULL) {
            //We need to add it..
            //Now copy across the actual flower.
            void
                    *vA =
                            binaryRepresentation_makeBinaryRepresentation(
                                    flower,
                                    (void(*)(void *, void(*)(const void * ptr,
                                            size_t size, size_t count))) metaSequence_writeBinaryRepresentation,
                                    &recordSize);
            //Compression
            vA = compress(vA, &recordSize);
            stKVDatabase_insertRecord(newCactusDisk->database,
                    metaSequence_getName(metaSequence), vA, recordSize);
        }
    }
    flower_destructSequenceIterator(sequenceIt);
    //Now copy across the actual flower.
    void *vA = binaryRepresentation_makeBinaryRepresentation(flower,
            (void(*)(void *, void(*)(const void * ptr, size_t size,
                    size_t count))) flower_writeBinaryRepresentation,
            &recordSize);
    //Compression
    vA = compress(vA, &recordSize);
    stKVDatabase_insertRecord(newCactusDisk->database, flower_getName(flower),
            vA, recordSize);
}

static void copyAcrossFlowers(Flower *flower, CactusDisk *newCactusDisk) {
    stKVDatabase_startTransaction(newCactusDisk->database);
    stTry {
        flowerIterator(flower, (void(*)(Flower *, void *))copyAcrossFlowersP,
                newCactusDisk);
        stKVDatabase_commitTransaction(newCactusDisk->database);
    }
    stCatch(except)
    {
        stKVDatabase_abortTransaction(newCactusDisk->database);
        stThrowNewCause(
                except,
                ST_KV_DATABASE_EXCEPTION_ID,
                "An unknown database error occurred when we tried to copy out a flower subtree\n");
    }
    stTryEnd;
}

void cactusDisk_splitFlowers(Flower *flower, CactusDisk *newCactusDisk) {
    Group *parentGroup = flower_getParentGroup(flower);
    flower_setParentGroup(flower, NULL);
    copyAcrossFlowers(flower, newCactusDisk);
    flower_setParentGroup(flower, parentGroup);
}

void cactusDisk_mergeFlowers(Flower *flower, CactusDisk *oldCactusDisk) {
    Flower *oldFlower = cactusDisk_getFlower(oldCactusDisk, flower_getName(flower));
    assert(oldFlower != NULL);
    assert(flower_getParentGroup(flower) == NULL);
    flower_setParentGroup(flower, flower_getParentGroup(oldFlower));
    copyAcrossFlowers(flower, oldCactusDisk);
}

void cactusDisk_relabelIDs(CactusDisk *cactusDisk) {


}

