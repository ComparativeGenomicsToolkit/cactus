/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

// For fsync declaration (technically a POSIX extension).
#define _POSIX_C_SOURCE 200809L

#include "cactusGlobalsPrivate.h"
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#define CACTUS_DISK_NAME_INCREMENT 16384
#define CACTUS_DISK_BUCKET_NUMBER 65536
#define CACTUS_DISK_PARAMETER_KEY -100000
#define CACTUS_DISK_SEQUENCE_CHUNK_SIZE 500

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

static char *getStringFromDisk(FILE *fileHandle, int64_t name, int64_t start, int64_t length) {
    int64_t k = fseek(fileHandle, name + start, SEEK_SET);
    if (k != 0) {
        st_errAbort("Could not fseek to start of desired sequence: %" PRIi64 "\n", name + start);
    }
    char *string = st_malloc(sizeof(char) * (length + 1));
    int64_t bytesRead = fread(string, sizeof(char), length, fileHandle);
    if (bytesRead != length) {
        st_errAbort("Read only %" PRIi64 " bytes of string of length %" PRIi64 " when caching substrings from DB\n",
                bytesRead, length);
    }
    string[length] = '\0';
#ifndef NDEBUG
    for (int64_t j = 0; j < length; j++) {
        assert(string[j] != '>');
        assert(string[j] != ' ');
    }
#endif
    return string;
}

Name cactusDisk_addString(CactusDisk *cactusDisk, const char *string) {
    /*
     * Adds a string to the database.
     */
    if (cactusDisk->storeSequencesInAFile) {
        if (cactusDisk->sequencesWriteFileHandle == NULL) {
            //We do not allow the read file handle to be open at the same time.
            if (cactusDisk->sequencesReadFileHandle != NULL) {
                fclose(cactusDisk->sequencesReadFileHandle);
                cactusDisk->sequencesReadFileHandle = NULL;
            }
            cactusDisk->sequencesWriteFileHandle = fopen(cactusDisk->absSequencesFileName, "a");
            assert(cactusDisk->sequencesWriteFileHandle != NULL);
        }
        else {
            //The read file handle should not be open at the same time.
            assert(cactusDisk->sequencesReadFileHandle == NULL);
        }
        Name name = ftell(cactusDisk->sequencesWriteFileHandle) + 1;

        //Extra temporary cheesy code to avoid potential overflow in fprintf
        int64_t chunkSize = 1000000000; //1 gig approx chunks, to avoid a possible overflow issue with fprintf
        int64_t length = strlen(string);
        if (length > chunkSize) {
            fprintf(cactusDisk->sequencesWriteFileHandle, ">");
            for (int64_t i = 0; i < length;) {
                int64_t j = i + chunkSize <= length ? chunkSize : length - i;
                char *string2 = memcpy(st_malloc(sizeof(char) * (j + 1)), string + i, sizeof(char) * j);
                string2[j] = '\0';
                int64_t k = fprintf(cactusDisk->sequencesWriteFileHandle, "%s", string2);
                (void) k;
                assert(k == j);
                free(string2);
                i += j;
            }
        } else {
            //Replacing this line
            int64_t k = fprintf(cactusDisk->sequencesWriteFileHandle, ">%s", string);
            (void) k;
            assert(k == length + 1);
        }

#ifndef NDEBUG
        // Extra fsync may not be necessary.
        fsync(fileno(cactusDisk->sequencesWriteFileHandle));
        fclose(cactusDisk->sequencesWriteFileHandle);
        cactusDisk->sequencesWriteFileHandle = NULL;
        cactusDisk->sequencesReadFileHandle = fopen(cactusDisk->absSequencesFileName, "r");
        char *string2 = getStringFromDisk(cactusDisk->sequencesReadFileHandle, name, 0, length);
        for (int64_t i = 0; i < length; i++) {
            assert(string[i] == string2[i]);
        }
        free(string2);
#endif
        return name;
    } else {
        int64_t stringSize = strlen(string);
        int64_t intervalSize = ceil((double) stringSize / CACTUS_DISK_SEQUENCE_CHUNK_SIZE);
        Name name = cactusDisk_getUniqueIDInterval(cactusDisk, intervalSize);
        stList *insertRequests = stList_construct3(0, (void (*)(void *)) stKVDatabaseBulkRequest_destruct);
        for (int64_t i = 0; i * CACTUS_DISK_SEQUENCE_CHUNK_SIZE < stringSize; i++) {
            int64_t j =
                    (i + 1) * CACTUS_DISK_SEQUENCE_CHUNK_SIZE < stringSize ?
                            CACTUS_DISK_SEQUENCE_CHUNK_SIZE : stringSize - i * CACTUS_DISK_SEQUENCE_CHUNK_SIZE;
            char *subString = stString_getSubString(string, i * CACTUS_DISK_SEQUENCE_CHUNK_SIZE, j);
            stList_append(insertRequests, stKVDatabaseBulkRequest_constructInsertRequest(name + i, subString, j + 1));
            free(subString);
        }
        stTry
            {
                stKVDatabase_bulkSetRecords(cactusDisk->database, insertRequests);
            }
            stCatch(except)
                {
                    stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                            "An unknown database error occurred when we tried to add a string to the cactus disk");
                }stTryEnd
        ;
        stList_destruct(insertRequests);
        return name;
    }
}

/*
 * Functions used to precache the sequences in the database for a given set of flowers.
 */

typedef struct _substring {
    /*
     * Struct used to represent a substring of a string to be retrieved from the database.
     */
    Name name;
    int64_t start;
    int64_t length;
} Substring;

/*
 * Basic methods on substrings.
 */

static Substring *substring_construct(Name name, int64_t start, int64_t length) {
    Substring *substring = st_malloc(sizeof(Substring));
    substring->name = name;
    substring->start = start;
    substring->length = length;
    return substring;
}

static Substring *substring_clone(Substring *substring) {
    return substring_construct(substring->name, substring->start, substring->length);
}

static void substring_destruct(Substring *substring) {
    free(substring);
}

static int substring_cmp(Substring *substring1, Substring *substring2) {
    int i = cactusMisc_nameCompare(substring1->name, substring2->name);
    if (i != 0) {
        return i;
    }
    i = substring1->start < substring2->start ? -1 : (substring1->start > substring2->start ? 1 : 0);
    if (i != 0) {
        return i;
    }
    return substring1->length < substring2->length ? -1 : (substring1->length > substring2->length ? 1 : 0);
}

static stList *mergeSubstrings(stList *substrings, int64_t proximityToMerge) {
    /*
     * Merge set of substrings into fewer substrings, if they overlap by less than proximityToMerge
     */
    stList *mergedSubstrings = stList_construct3(0, (void (*)(void *)) substring_destruct);
    if (stList_length(substrings) == 0) {
        return mergedSubstrings;
    }
    stList_sort(substrings, (int (*)(const void *, const void *)) substring_cmp);
    Substring *pSubsequence = substring_clone(stList_get(substrings, 0));
    stList_append(mergedSubstrings, pSubsequence);
    for (int64_t i = 1; i < stList_length(substrings); i++) {
        Substring *substring = stList_get(substrings, i);
        if (pSubsequence->name == substring->name
                && pSubsequence->start + pSubsequence->length + proximityToMerge >= substring->start) { //Merge
            if (pSubsequence->start + pSubsequence->length < substring->start + substring->length) {
                pSubsequence->length = substring->start + substring->length - pSubsequence->start;
            }
        } else {
            pSubsequence = substring_clone(substring);
            stList_append(mergedSubstrings, pSubsequence);
        }
    }
    return mergedSubstrings;
}

static void cacheSubstringsFromDB(CactusDisk *cactusDisk, stList *substrings) {
    /*
     * Caches the given set of substrings in the cactusDisk cache.
     */
    if (cactusDisk->storeSequencesInAFile) {
        if (cactusDisk->sequencesReadFileHandle == NULL) {
            if(cactusDisk->sequencesWriteFileHandle != NULL) {
                fsync(fileno(cactusDisk->sequencesWriteFileHandle));
                fclose(cactusDisk->sequencesWriteFileHandle);
                cactusDisk->sequencesWriteFileHandle = NULL;
            }
            cactusDisk->sequencesReadFileHandle = fopen(cactusDisk->absSequencesFileName, "r");
            assert(cactusDisk->sequencesReadFileHandle != NULL);
        }
        else {
            assert(cactusDisk->sequencesWriteFileHandle == NULL);
        }
        for (int64_t i = 0; i < stList_length(substrings); i++) {
            Substring *substring = stList_get(substrings, i);
            char *string = getStringFromDisk(cactusDisk->sequencesReadFileHandle, substring->name, substring->start,
                    substring->length);
            stCache_setRecord(cactusDisk->stringCache, substring->name, substring->start, substring->length, string);
#ifndef NDEBUG
            int64_t bytesRead;
            char *string2 = stCache_getRecord(cactusDisk->stringCache, substring->name, substring->start,
                    substring->length, &bytesRead);
            assert(bytesRead == substring->length);
            for (int64_t j = 0; j < substring->length; j++) {
                assert(string2[j] == string[j]);
            }
            free(string2);
#endif
            free(string);
        }
    } else {
        stList *getRequests = stList_construct3(0, free);
        for (int64_t i = 0; i < stList_length(substrings); i++) {
            Substring *substring = stList_get(substrings, i);
            int64_t intervalSize = (substring->length + substring->start - 1) / CACTUS_DISK_SEQUENCE_CHUNK_SIZE
                    - substring->start / CACTUS_DISK_SEQUENCE_CHUNK_SIZE + 1;
            Name shiftedName = substring->name + substring->start / CACTUS_DISK_SEQUENCE_CHUNK_SIZE;
            for (int64_t j = 0; j < intervalSize; j++) {
                int64_t *k = st_malloc(sizeof(int64_t));
                k[0] = shiftedName + j;
                stList_append(getRequests, k);
            }
        }
        if (stList_length(getRequests) == 0) {
            stList_destruct(getRequests);
            return;
        }
        stList *records = NULL;
        stTry
            {
                records = stKVDatabase_bulkGetRecords(cactusDisk->database, getRequests);
            }
            stCatch(except)
                {
                    stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                            "An unknown database error occurred when getting a sequence string");
                }stTryEnd
        ;
        assert(records != NULL);
        assert(stList_length(records) == stList_length(getRequests));
        stList_destruct(getRequests);
        stListIterator *recordsIt = stList_getIterator(records);
        for (int64_t i = 0; i < stList_length(substrings); i++) {
            Substring *substring = stList_get(substrings, i);
            int64_t intervalSize = (substring->length + substring->start - 1) / CACTUS_DISK_SEQUENCE_CHUNK_SIZE
                    - substring->start / CACTUS_DISK_SEQUENCE_CHUNK_SIZE + 1;
            stList *strings = stList_construct();
            while (intervalSize-- > 0) {
                int64_t recordSize;
                stKVDatabaseBulkResult *result = stList_getNext(recordsIt);
                assert(result != NULL);
                char *string = stKVDatabaseBulkResult_getRecord(result, &recordSize);
                assert(string != NULL);
                assert(strlen(string) == recordSize - 1);
                stList_append(strings, string);
                assert(recordSize <= CACTUS_DISK_SEQUENCE_CHUNK_SIZE + 1);
            }
            assert(stList_length(strings) > 0);
            char *joinedString = stString_join2("", strings);
            stCache_setRecord(cactusDisk->stringCache, substring->name,
                    (substring->start / CACTUS_DISK_SEQUENCE_CHUNK_SIZE) * CACTUS_DISK_SEQUENCE_CHUNK_SIZE,
                    strlen(joinedString), joinedString);
            free(joinedString);
            stList_destruct(strings);
        }
        assert(stList_getNext(recordsIt) == NULL);
        stList_destructIterator(recordsIt);
        stList_destruct(records);
    }
}

void cactusDisk_preCacheStrings2(CactusDisk *cactusDisk, stList *substrings) {
    /*
     * Precaches the given substrings, so that they are all in memory.
     */
    //Now do some simple merging to reduce granularity
    stList *mergedSubstrings = mergeSubstrings(substrings, CACTUS_DISK_SEQUENCE_CHUNK_SIZE);
    //Now cache the sequences
    cacheSubstringsFromDB(cactusDisk, mergedSubstrings);
    stList_destruct(mergedSubstrings);
}

static stList *getSubstringsForFlowers(stList *flowers) {
    /*
     * Get the set of substrings for sequence intervals in the given set of flowers.
     */
    stList *substrings = stList_construct3(0, (void (*)(void *)) substring_destruct);
    for (int64_t i = 0; i < stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        Flower_EndIterator *endIt = flower_getEndIterator(flower);
        End *end;
        while ((end = flower_getNextEnd(endIt)) != NULL) {
            if (end_isStubEnd(end)) {
                End_InstanceIterator *instanceIt = end_getInstanceIterator(end);
                Cap *cap;
                while ((cap = end_getNext(instanceIt)) != NULL) {
                    Sequence *sequence;
                    if ((sequence = cap_getSequence(cap)) != NULL) {
                        cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                        if (!cap_getSide(cap)) { //We have a sequence interval of interest
                            Cap *adjacentCap = cap_getAdjacency(cap);
                            assert(adjacentCap != NULL);
                            int64_t length = cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) - 1;
                            assert(length >= 0);
                            if (length > 0) {
                                stList_append(substrings,
                                        substring_construct(sequence_getMetaSequence(sequence)->stringName,
                                                cap_getCoordinate(cap) + 1 - sequence_getStart(sequence), length));
                            }
                        }
                    }
                }
                end_destructInstanceIterator(instanceIt);
            }
        }
        flower_destructEndIterator(endIt);
    }
    return substrings;
}

void cactusDisk_preCacheStrings(CactusDisk *cactusDisk, stList *flowers) {
    /*
     * Precaches the sequences in the set of flowers, so that they are all in memory.
     */
    stList *substrings = getSubstringsForFlowers(flowers);
    cactusDisk_preCacheStrings2(cactusDisk, substrings);
    stList_destruct(substrings);
}

static stList *getSubstringsForFlowerSegments(stList *flowers) {
    /*
     * Get the set of substrings representing the strings in the segments of the given flowers.
     */
    stList *substrings = stList_construct3(0, (void (*)(void *)) substring_destruct);
    for (int64_t i = 0; i < stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        Flower_EndIterator *blockIt = flower_getBlockIterator(flower);
        Block *block;
        while ((block = flower_getNextBlock(blockIt)) != NULL) {
            Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
            Segment *segment;
            while ((segment = block_getNext(instanceIt)) != NULL) {
                Sequence *sequence;
                if ((sequence = segment_getSequence(segment)) != NULL) {
                    segment = segment_getStrand(segment) ? segment : segment_getReverse(segment);
                    assert(segment_getLength(segment) > 0);
                    stList_append(substrings,
                            substring_construct(sequence_getMetaSequence(sequence)->stringName,
                                    segment_getStart(segment) - sequence_getStart(sequence),
                                    segment_getLength(segment)));
                }
            }
            block_destructInstanceIterator(instanceIt);
        }
        flower_destructBlockIterator(blockIt);
    }
    return substrings;
}

void cactusDisk_preCacheSegmentStrings(CactusDisk *cactusDisk, stList *flowers) {
    /*
     * Precaches the sequences in the set blocks of the given flowers, so that they are all in memory.
     */
    stList *substrings = getSubstringsForFlowerSegments(flowers);
    cactusDisk_preCacheStrings2(cactusDisk, substrings);
    stList_destruct(substrings);
}

char *cactusDisk_getStringFromCache(CactusDisk *cactusDisk, Name name, int64_t start, int64_t length, int64_t strand) {
    /*
     * Gets a sequence from the cache.
     */
    char *string = NULL;
    if (stCache_containsRecord(cactusDisk->stringCache, name, start, sizeof(char) * length)) {
        int64_t recordSize;
        string = stCache_getRecord(cactusDisk->stringCache, name, start, sizeof(char) * length, &recordSize);
        assert(string != NULL);
        assert(recordSize == length);
        string = st_realloc(string, sizeof(char) * (length + 1));
        string[length] = '\0';
        if (!strand) {
            char *string2 = stString_reverseComplementString(string);
            free(string);
            string = string2;
        }
    }
    return string;
}

char *cactusDisk_getString(CactusDisk *cactusDisk, Name name, int64_t start, int64_t length, int64_t strand,
        int64_t totalSequenceLength) {
    /*
     * Gets a string from the database.
     *
     */
    assert(length >= 0);
    if (length == 0) {
        return stString_copy("");
    }
    //First try getting it from the cache
    char *string = cactusDisk_getStringFromCache(cactusDisk, name, start, length, strand);
    if (string == NULL) { //If not in the cache, add it to the cache and then get it from the cache.
        stList *list = stList_construct3(0, (void (*)(void *)) substring_destruct);
        stList_append(list, substring_construct(name, start, length));
        cacheSubstringsFromDB(cactusDisk, list);
        stList_destruct(list);
        string = cactusDisk_getStringFromCache(cactusDisk, name, start, length, strand);
    }
    assert(string != NULL);
    return string;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic flower disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

const char *CACTUS_DISK_EXCEPTION_ID = "CACTUS_DISK_EXCEPTION_ID";

static int cactusDisk_constructFlowersP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(flower_getName((Flower *) o1), flower_getName((Flower *) o2));
}

static int cactusDisk_constructMetaSequencesP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(metaSequence_getName((MetaSequence *) o1), metaSequence_getName((MetaSequence *) o2));
}

static void cactusDisk_writeBinaryRepresentation(CactusDisk *cactusDisk,
        void (*writeFn)(const void * ptr, size_t size, size_t count)) {
    binaryRepresentation_writeElementType(CODE_CACTUS_DISK, writeFn);
    binaryRepresentation_writeBool(cactusDisk->storeSequencesInAFile, writeFn);
    if (cactusDisk->storeSequencesInAFile) {
        assert(cactusDisk->sequencesFileName != NULL);
        binaryRepresentation_writeString(cactusDisk->sequencesFileName, writeFn);
    }
    binaryRepresentation_writeElementType(CODE_CACTUS_DISK, writeFn);
}

static void cactusDisk_loadFromBinaryRepresentation(void **binaryString, CactusDisk *cactusDisk, stKVDatabaseConf *conf) {
    cactusDisk->sequencesReadFileHandle = NULL;
    cactusDisk->sequencesWriteFileHandle = NULL; //I think these lines are not needed.
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

static void *compress(void *data, int64_t *dataSize) {
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
    *dataSize = uncompressedSize;
    return data2;
}

static stList *getRecords(CactusDisk *cactusDisk, stList *objectNames, char *type) {
    if (stList_length(objectNames) == 0) {
        return stList_construct3(0, NULL);
    }
    stList *records = NULL;
    stTry
        {
            records = stKVDatabase_bulkGetRecords(cactusDisk->database, objectNames);
        }
        stCatch(except)
            {
                stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                        "An unknown database error occurred when getting a bulk set of %s", type);
            }stTryEnd
    ;
    assert(records != NULL);
    assert(stList_length(objectNames) == stList_length(records));
    stList_setDestructor(records, free);
    for (int64_t i = 0; i < stList_length(objectNames); i++) {
        Name objectName = *((int64_t *) stList_get(objectNames, i));
        int64_t recordSize;
        void *record;
        stKVDatabaseBulkResult *result = stList_get(records, i);
        assert(result != NULL);
        if (!stCache_containsRecord(cactusDisk->cache, objectName, 0, INT64_MAX)) {
            record = stKVDatabaseBulkResult_getRecord(result, &recordSize);
            assert(recordSize >= 0);
            assert(record != NULL);
            record = decompress(record, &recordSize);
            stCache_setRecord(cactusDisk->cache, objectName, 0, recordSize, record);
        } else {
            record = stCache_getRecord(cactusDisk->cache, objectName, 0, INT64_MAX, &recordSize);
            assert(recordSize >= 0);
            assert(record != NULL);
        }
        stKVDatabaseBulkResult_destruct(result);
        stList_set(records, i, record);
    }
    return records;
}

static void *getRecord(CactusDisk *cactusDisk, Name objectName, char *type) {
    void *cA = NULL;
    int64_t recordSize = 0;
    if (stCache_containsRecord(cactusDisk->cache, objectName, 0, INT64_MAX)) { //If we already have the record, we won't update it.
        cA = stCache_getRecord(cactusDisk->cache, objectName, 0, INT64_MAX, &recordSize);
    } else {
        stTry
            {
                cA = stKVDatabase_getRecord2(cactusDisk->database, objectName, &recordSize);
            }
            stCatch(except)
                {
                    stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                            "An unknown database error occurred when getting a %s", type);
                }stTryEnd
        ;
        if (cA == NULL) {
            return NULL;
        }
        stCache_setRecord(cactusDisk->cache, objectName, 0, recordSize, cA); //Add the compressed record to the cache.
    }
    //Decompression
    void *cA2 = decompress(cA, &recordSize);
    free(cA);
    return cA2;
}

static bool containsRecord(CactusDisk *cactusDisk, Name objectName) {
    return stCache_containsRecord(cactusDisk->cache, objectName, 0, INT64_MAX)
            || stKVDatabase_containsRecord(cactusDisk->database, objectName);
}

static CactusDisk *cactusDisk_constructPrivate(stKVDatabaseConf *conf, bool create, const char *sequencesFileName) {
    //sequencesFileName = NULL; //Disable the ability to store the sequences on disk.
    CactusDisk *cactusDisk = st_calloc(1, sizeof(CactusDisk));

    //construct lists of in memory objects
    cactusDisk->metaSequences = stSortedSet_construct3(cactusDisk_constructMetaSequencesP, NULL);
    cactusDisk->flowers = stSortedSet_construct3(cactusDisk_constructFlowersP, NULL);
    cactusDisk->flowerNamesMarkedForDeletion = stSortedSet_construct3((int (*)(const void *, const void *)) strcmp,
            free);
    cactusDisk->updateRequests = stList_construct3(0, (void (*)(void *)) stKVDatabaseBulkRequest_destruct);

    //Now open the database
    cactusDisk->database = stKVDatabase_construct(conf, create);
    cactusDisk->cache = stCache_construct();
    cactusDisk->stringCache = stCache_construct();

    //initialise the unique ids.
    int64_t seed = (clock() << 24) | (time(NULL) << 16) | (getpid() & 65535); //Likely to be unique
    st_logDebug("The cactus disk is seeding the random number generator with the value %" PRIi64 "\n", seed);
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
            cactusDisk->sequencesReadFileHandle = NULL;
            cactusDisk->sequencesWriteFileHandle = NULL;
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
            //Make sure the file exists
            cactusDisk->sequencesReadFileHandle = fopen(cactusDisk->absSequencesFileName, "w");
            assert(cactusDisk->sequencesReadFileHandle != NULL);
            fclose(cactusDisk->sequencesReadFileHandle); //Flush it first time.
            cactusDisk->sequencesReadFileHandle = NULL;
            cactusDisk->sequencesWriteFileHandle = NULL;
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
        if (cactusDisk->sequencesReadFileHandle != NULL) {
            fclose(cactusDisk->sequencesReadFileHandle);
        }
        if (cactusDisk->sequencesWriteFileHandle != NULL) {
            fsync(fileno(cactusDisk->sequencesWriteFileHandle));
            fclose(cactusDisk->sequencesWriteFileHandle);
        }
    } else {
        assert(cactusDisk->sequencesFileName == NULL);
        assert(cactusDisk->sequencesReadFileHandle == NULL);
        assert(cactusDisk->absSequencesFileName == NULL);
    }

    stCache_destruct(cactusDisk->cache); //Get rid of the cache
    stCache_destruct(cactusDisk->stringCache);

    stList_destruct(cactusDisk->updateRequests);

    free(cactusDisk);
}

void cactusDisk_addUpdateRequest(CactusDisk *cactusDisk, Flower *flower) {
    int64_t recordSize;
    void *vA = binaryRepresentation_makeBinaryRepresentation(flower,
            (void (*)(void *, void (*)(const void * ptr, size_t size, size_t count))) flower_writeBinaryRepresentation,
            &recordSize);
    //Compression
    vA = compress(vA, &recordSize);
    if (containsRecord(cactusDisk, flower_getName(flower))) {
        int64_t recordSize2;
        void *vA2 = stCache_getRecord(cactusDisk->cache, flower_getName(flower), 0, INT64_MAX, &recordSize2);
        if (!stCache_recordsIdentical(vA, recordSize, vA2, recordSize2)) { //Only rewrite if we actually did something
            stList_append(cactusDisk->updateRequests,
                    stKVDatabaseBulkRequest_constructUpdateRequest(flower_getName(flower), vA, recordSize));
        }
        free(vA2);
    } else {
        stList_append(cactusDisk->updateRequests,
                stKVDatabaseBulkRequest_constructInsertRequest(flower_getName(flower), vA, recordSize));
    }
    free(vA);
}

void cactusDisk_write(CactusDisk *cactusDisk) {
    Flower *flower;
    int64_t recordSize;

    stList *removeRequests = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    st_logDebug("Starting to write the cactus to disk\n");

    stSortedSetIterator *it = stSortedSet_getIterator(cactusDisk->flowers);
    //Sort flowers to update.
    while ((flower = stSortedSet_getNext(it)) != NULL) {
        cactusDisk_addUpdateRequest(cactusDisk, flower);
    }
    stSortedSet_destructIterator(it);

    st_logDebug("Got the flowers to update\n");

    //Remove nets that are marked for deletion..
    it = stSortedSet_getIterator(cactusDisk->flowerNamesMarkedForDeletion);
    char *nameString;
    while ((nameString = stSortedSet_getNext(it)) != NULL) {
        Name name = cactusMisc_stringToName(nameString);
        if (containsRecord(cactusDisk, name)) {
            stList_append(cactusDisk->updateRequests, stKVDatabaseBulkRequest_constructUpdateRequest(name, &name, 0)); //We set it to null in the first atomic operation.
            stList_append(removeRequests, stIntTuple_construct1(name));
        }
    }
    stSortedSet_destructIterator(it);

    st_logDebug("Avoided updating nets marked for deletion\n");

    // Insert and/or update meta-sequences.
    it = stSortedSet_getIterator(cactusDisk->metaSequences);
    MetaSequence *metaSequence;
    while ((metaSequence = stSortedSet_getNext(it)) != NULL) {
        void *vA =
                binaryRepresentation_makeBinaryRepresentation(metaSequence,
                        (void (*)(void *, void (*)(const void * ptr, size_t size, size_t count))) metaSequence_writeBinaryRepresentation,
                        &recordSize);
        //Compression
        vA = compress(vA, &recordSize);
        if (!containsRecord(cactusDisk, metaSequence_getName(metaSequence))) {
            stList_append(cactusDisk->updateRequests,
                    stKVDatabaseBulkRequest_constructInsertRequest(metaSequence_getName(metaSequence), vA, recordSize));
        } else {
            stList_append(cactusDisk->updateRequests,
                    stKVDatabaseBulkRequest_constructUpdateRequest(metaSequence_getName(metaSequence), vA, recordSize));
        }
        free(vA);
    }
    stSortedSet_destructIterator(it);

    st_logDebug("Got the sequences we are going to add to the database.\n");

    if (!containsRecord(cactusDisk, CACTUS_DISK_PARAMETER_KEY)) { //We only write the parameters once.
        //Finally the database info.
        void *cactusDiskParameters =
                binaryRepresentation_makeBinaryRepresentation(cactusDisk,
                        (void (*)(void *, void (*)(const void * ptr, size_t size, size_t count))) cactusDisk_writeBinaryRepresentation,
                        &recordSize);
        //Compression
        cactusDiskParameters = compress(cactusDiskParameters, &recordSize);
        stList_append(cactusDisk->updateRequests,
                stKVDatabaseBulkRequest_constructInsertRequest(CACTUS_DISK_PARAMETER_KEY, cactusDiskParameters,
                        recordSize));
        free(cactusDiskParameters);
    }

    st_logDebug("Checked if need to write the initial parameters\n");

    if (stList_length(cactusDisk->updateRequests) > 0) {
        st_logDebug("Going to write %" PRIi64 " updates\n", stList_length(cactusDisk->updateRequests));
        stTry
            {
                st_logDebug("Writing %" PRIi64 " updates\n", stList_length(cactusDisk->updateRequests));
                assert(stList_length(cactusDisk->updateRequests) > 0);
                stKVDatabase_bulkSetRecords(cactusDisk->database, cactusDisk->updateRequests);
            }
            stCatch(except)
                {
                    stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                            "Failed when trying to set records in updating the cactus disk");
                }stTryEnd
        ;
    }

    st_logDebug("Updated the database with inserts\n");

    if (stList_length(removeRequests) > 0) {
        stTry
            {
                stKVDatabase_bulkRemoveRecords(cactusDisk->database, removeRequests);
            }
            stCatch(except)
                {
                    stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                            "Failed when trying to remove records in updating the cactus disk");
                }stTryEnd
        ;
    }

    st_logDebug("Now removed flowers we don't need\n");

    stList_destruct(cactusDisk->updateRequests);
    cactusDisk->updateRequests = stList_construct3(0, (void (*)(void *)) stKVDatabaseBulkRequest_destruct);
    stList_destruct(removeRequests);

    st_logDebug("Finished writing to the database\n");
}

stList *cactusDisk_getFlowers(CactusDisk *cactusDisk, stList *flowerNames) {
    stList *records = getRecords(cactusDisk, flowerNames, "flowers");
    assert(stList_length(flowerNames) == stList_length(records));
    stList *flowers = stList_construct();
    for (int64_t i = 0; i < stList_length(flowerNames); i++) {
        Name flowerName = *((int64_t *) stList_get(flowerNames, i));
        static Flower flower;
        flower.name = flowerName;
        Flower *flower2;
        if ((flower2 = stSortedSet_search(cactusDisk->flowers, &flower)) == NULL) {
            void *record = stList_get(records, i);
            assert(record != NULL);
            void *cA = record;
            flower2 = flower_loadFromBinaryRepresentation(&cA, cactusDisk);
            assert(flower2 != NULL);
        }
        stList_append(flowers, flower2);
    }
    stList_destruct(records);
    return flowers;
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
 * Function to get unique ID.
 */

void cactusDisk_getBlockOfUniqueIDs(CactusDisk *cactusDisk, int64_t intervalSize) {
    intervalSize = intervalSize < CACTUS_DISK_NAME_INCREMENT ? CACTUS_DISK_NAME_INCREMENT : intervalSize;
    bool done = 0;
    int64_t collisionCount = 0;
    while (!done) {
        stTry
            {
                Name keyName = st_randomInt(-CACTUS_DISK_BUCKET_NUMBER, 0);
                assert(keyName >= -CACTUS_DISK_BUCKET_NUMBER);
                assert(keyName < 0);
                int64_t bucketSize = INT64_MAX / CACTUS_DISK_BUCKET_NUMBER;
                int64_t minimumValue = bucketSize * (llabs(keyName) - 1) + 1; //plus one for the reserved '0' value.
                int64_t maximumValue = minimumValue + (bucketSize - 1);
                assert(minimumValue >= 1);
                assert(maximumValue <= INT64_MAX);
                assert(minimumValue < maximumValue);
                if (stKVDatabase_containsRecord(cactusDisk->database, keyName)) {
                    cactusDisk->maxUniqueNumber = stKVDatabase_incrementInt64(cactusDisk->database, keyName,
                            intervalSize);
                    cactusDisk->uniqueNumber = cactusDisk->maxUniqueNumber - intervalSize;
                    if (cactusDisk->uniqueNumber <= 0 || cactusDisk->uniqueNumber < minimumValue
                            || cactusDisk->uniqueNumber > maximumValue) {
                        st_errAbort("Got a non positive unique number %lli %lli %lli %lli", cactusDisk->uniqueNumber,
                                cactusDisk->maxUniqueNumber, minimumValue, maximumValue);
                    }
                    assert(cactusDisk->uniqueNumber >= minimumValue);
                    assert(cactusDisk->uniqueNumber <= maximumValue);
                    assert(cactusDisk->uniqueNumber > 0);
                } else {
                    stTry
                        {
                            stKVDatabase_insertInt64(cactusDisk->database, keyName, minimumValue);
                        }
                        stCatch(except)
                            {
                                collisionCount++;
                                if (collisionCount >= 10) {
                                    stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                                            "Repeated unknown database errors occurred when we tried to get a unique ID, collision count %" PRIi64 "",
                                            collisionCount);
                                } else {
                                    st_logDebug("Got an exception when trying to insert a uid record: %s",
                                            stExcept_getMsg(except));
                                }
                            }stTryEnd
                    ;
                    continue;
                }
                if (cactusDisk->maxUniqueNumber >= maximumValue) {
                    st_errAbort("We have exhausted a bucket, which seems really unlikely");
                }
                done = 1;
            }
            stCatch(except)
                {
                    collisionCount++;
                    if (collisionCount >= 10) {
                        stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                                "Repeated unknown database errors occurred when we tried to get a unique ID, collision count %" PRIi64 "",
                                collisionCount);
                    } else {
                        st_logDebug("Got an exception when trying to insert a record: %s", stExcept_getMsg(except));
                    }
                }stTryEnd
        ;
    }
}

int64_t cactusDisk_getUniqueIDInterval(CactusDisk *cactusDisk, int64_t intervalSize) {
    assert(cactusDisk->uniqueNumber <= cactusDisk->maxUniqueNumber);
    if (cactusDisk->uniqueNumber + intervalSize > cactusDisk->maxUniqueNumber) {
        cactusDisk_getBlockOfUniqueIDs(cactusDisk, intervalSize);
    }
    Name uniqueNumber = cactusDisk->uniqueNumber;
    cactusDisk->uniqueNumber += intervalSize;
    return uniqueNumber;
}

int64_t cactusDisk_getUniqueID(CactusDisk *cactusDisk) {
    return cactusDisk_getUniqueIDInterval(cactusDisk, 1);
}

