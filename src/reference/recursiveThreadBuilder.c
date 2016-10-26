/*
 * recursiveFileBuilder.c
 *
 *  Created on: 16 Mar 2012
 *      Author: benedictpaten
 */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "cactus.h"
#include "sonLib.h"

static void *compress(char *string, int64_t *dataSize) {
    void *data = stCompression_compress(string, strlen(string) + 1, dataSize, 1); //going with least, fastest compression-1);
    free(string);
    return data;
}

static char *decompress(void *data, int64_t dataSize) {
    int64_t uncompressedSize;
    char *string = stCompression_decompress(data, dataSize, &uncompressedSize);
    assert(strlen(string)+1 == uncompressedSize);
    free(data);
    return string;
}

static void cacheNonNestedRecords(stCache *cache, stList *caps, char *(*segmentWriteFn)(Segment *),
        char *(*terminalAdjacencyWriteFn)(Cap *)) {
    /*
     * Caches the set of terminal adjacency and segment records present in the threads.
     */
    for (int64_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        int64_t recordSize;
        while (1) {
            Cap *adjacentCap = cap_getAdjacency(cap);
            assert(adjacentCap != NULL);
            Group *group = end_getGroup(cap_getEnd(cap));
            assert(group != NULL);
            if (group_isLeaf(group)) { //Record must not be in the database already
                void *data = compress(terminalAdjacencyWriteFn(cap), &recordSize);
                assert(!stCache_containsRecord(cache, cap_getName(cap), 0, INT64_MAX));
                stCache_setRecord(cache, cap_getName(cap), 0, recordSize, data);
                free(data);
            }
            if ((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
                break;
            }
            Segment *segment = cap_getSegment(adjacentCap);
            assert(!stCache_containsRecord(cache, segment_getName(segment), 0, INT64_MAX));
            void *data = compress(segmentWriteFn(segment), &recordSize);
            stCache_setRecord(cache, segment_getName(segment), 0, recordSize, data);
            free(data);
        }
    }
}

static stList *getNestedRecordNames(stList *caps) {
    /*
     * Gets the names of non-terminal adjacencies as a list of cap names.
     */
    stList *getRequests = stList_construct3(0, free);
    for (int64_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        while (1) {
            Cap *adjacentCap = cap_getAdjacency(cap);
            assert(adjacentCap != NULL);
            Group *group = end_getGroup(cap_getEnd(cap));
            assert(group != NULL);
            if (!group_isLeaf(group)) { //Record must be in the database already
                int64_t *j = st_malloc(sizeof(int64_t));
                j[0] = cap_getName(cap);
                stList_append(getRequests, j);
            }
            if ((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
                break;
            }
        }
    }
    return getRequests;
}

static void cacheNestedRecords(stKVDatabase *database, stCache *cache, stList *caps) {
    /*
     * Caches all the non-terminal adjacencies by retrieving them from the database.
     */
    stList *getRequests = getNestedRecordNames(caps);
    if (stList_length(caps) > 10000) {
        st_logCritical("Going to request %" PRIi64 " records from the database: %" PRIi64 "\n", stList_length(caps));
    }
    //Do the retrieval of the records
    stList *records = NULL;
    stTry {
            records = stKVDatabase_bulkGetRecords(database, getRequests);
        }stCatch(except)
            {
                stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                        "An unknown database error occurred when we tried to bulk get records from the database");
            }stTryEnd;
    assert(records != NULL);
    assert(stList_length(records) == stList_length(getRequests));
    //Now cache the resulting records
    while (stList_length(records) > 0) {
        stKVDatabaseBulkResult *result = stList_pop(records);
        int64_t *recordName = stList_pop(getRequests);
        int64_t recordSize;
        void *record = stKVDatabaseBulkResult_getRecord(result, &recordSize);
        assert(record != NULL);
        assert(!stCache_containsRecord(cache, *recordName, 0, INT64_MAX));
        stCache_setRecord(cache, *recordName, 0, recordSize, record);
        stKVDatabaseBulkResult_destruct(result); //Cleanup the memory as we go.
        free(recordName);
    }
    assert(stList_length(getRequests) == 0);
    stList_destruct(getRequests);
    stList_destruct(records);
}

static stCache *cacheRecords(stKVDatabase *database, stList *caps, char *(*segmentWriteFn)(Segment *),
        char *(*terminalAdjacencyWriteFn)(Cap *)) {
    /*
     * Cache all the elements needed to construct the set of threads.
     */
    stCache *cache = stCache_construct();
    cacheNestedRecords(database, cache, caps);
    cacheNonNestedRecords(cache, caps, segmentWriteFn, terminalAdjacencyWriteFn);
    return cache;
}

static void deleteNestedRecords(stKVDatabase *database, stList *caps) {
    /*
     * Removes the non-terminal adjacencies from the database.
     */
    stList *deleteRequests = getNestedRecordNames(caps);
    for (int64_t i = 0; i < stList_length(deleteRequests); i++) {
        int64_t *record = stList_get(deleteRequests, i);
        stList_set(deleteRequests, i, stIntTuple_construct1( record[0])); //Hack
        free(record);
    }
    stList_setDestructor(deleteRequests, (void(*)(void *)) stIntTuple_destruct);
    //Do the deletion of the records
    stTry {
            stKVDatabase_bulkRemoveRecords(database, deleteRequests);
        }stCatch(except)
            {
                stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                        "An unknown database error occurred when we tried to bulk remove records from the database");
            }stTryEnd;
    stList_destruct(deleteRequests);
}

static char *getThread(stCache *cache, Cap *startCap) {
    /*
     * Iterate through, first calculating the length of the final record, then concatenating the results.
     */
    Cap *cap = startCap;
    stList *strings = stList_construct3(0, free);
    while (1) { //Calculate the size of the entries in the DB that represent the thread.
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        int64_t recordSize;
        assert(stCache_containsRecord(cache, cap_getName(cap), 0, INT64_MAX));
        void *data = stCache_getRecord(cache, cap_getName(cap), 0, INT64_MAX, &recordSize);
        stList_append(strings, decompress(data, recordSize));
        if ((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
            break;
        }
        assert(stCache_containsRecord(cache, segment_getName(cap_getSegment(adjacentCap)), 0, INT64_MAX));
        data = stCache_getRecord(cache, segment_getName(cap_getSegment(adjacentCap)), 0, INT64_MAX, &recordSize);
        stList_append(strings, decompress(data, recordSize));
    }
    char *string = stString_join2("", strings);
    stList_destruct(strings);
    return string;
}

void buildRecursiveThreads(stKVDatabase *database, stList *caps, char *(*segmentWriteFn)(Segment *),
        char *(*terminalAdjacencyWriteFn)(Cap *)) {
    //Cache records
    stCache *cache = cacheRecords(database, caps, segmentWriteFn, terminalAdjacencyWriteFn);

    //Build new threads
    stList *records = stList_construct3(0, (void(*)(void *)) stKVDatabaseBulkRequest_destruct);
    for (int64_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        char *string = getThread(cache, cap);
        assert(string != NULL);
        int64_t recordSize;
        void *data = compress(string, &recordSize);
        stList_append(records, stKVDatabaseBulkRequest_constructInsertRequest(cap_getName(cap), data, recordSize));
        free(data);
    }

    //Delete old records and insert new records
    deleteNestedRecords(database, caps);
    stTry {
            stKVDatabase_bulkSetRecords(database, records);
        }stCatch(except)
            {
                stThrowNewCause(except, ST_KV_DATABASE_EXCEPTION_ID,
                        "An unknown database error occurred when we tried to bulk insert records from the database");
            }stTryEnd;

    //Cleanup
    stCache_destruct(cache);
    stList_destruct(records);
}

stList *buildRecursiveThreadsInList(stKVDatabase *database, stList *caps, char *(*segmentWriteFn)(Segment *),
        char *(*terminalAdjacencyWriteFn)(Cap *)) {
    stList *threadStrings = stList_construct3(0, free);

    //Cache records
    stCache *cache = cacheRecords(database, caps, segmentWriteFn, terminalAdjacencyWriteFn);

    //Build new threads
    for (int64_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        stList_append(threadStrings, getThread(cache, cap));
    }

    stCache_destruct(cache);

    return threadStrings;
}

