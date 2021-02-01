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
#include "recursiveThreadBuilder.h"

// OpenMP
#if defined(_OPENMP)
#include <omp.h>
#endif

//static omp_lock_t writelock;

RecordHolder *recordHolder_construct() {
    //omp_init_lock(&writelock);
    return stHash_construct2(NULL, free);
}

void recordHolder_destruct(RecordHolder *rh) {
    stHash_destruct(rh);
    //omp_destroy_lock(&writelock);
}

/*static void *compress(char *string, int64_t *dataSize) {
    *dataSize = (1+strlen(string)) * sizeof(char);
    return string;

    //void *data = stCompression_compress(string, strlen(string) + 1, dataSize, 1); //going with least, fastest compression-1);
    //free(string);
    //return data;
}

static char *decompress(void *data, int64_t dataSize) {
    return data;
    int64_t uncompressedSize;
    char *string = stCompression_decompress(data, dataSize, &uncompressedSize);
    assert(strlen(string)+1 == uncompressedSize);
    free(data);
    return string;
}

static char *recordHolder_get(RecordHolder *rh, Name name) {
    return stHash_search(rh, name);
}*/

static void recordHolder_add(RecordHolder *rh, Name name, char *string) {
    //omp_set_lock(&writelock);
    assert(stHash_search(rh, (void *)name) == NULL);
    stHash_insert(rh, (void *)name, string);
    //omp_unset_lock(&writelock);
}

static char *recordHolder_remove(RecordHolder *rh, Name name) {
    //omp_set_lock(&writelock);
    char *string = stHash_remove(rh, (void *)name);
    //omp_unset_lock(&writelock);
    return string;
}

static void cacheNonNestedRecords(RecordHolder *rh, stList *caps, char *(*segmentWriteFn)(Segment *),
        char *(*terminalAdjacencyWriteFn)(Cap *)) {
    /*
     * Caches the set of terminal adjacency and segment records present in the threads.
     */
    for (int64_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        //int64_t recordSize;
        while (1) {
            Cap *adjacentCap = cap_getAdjacency(cap);
            assert(adjacentCap != NULL);
            Group *group = end_getGroup(cap_getEnd(cap));
            assert(group != NULL);
            if (group_isLeaf(group)) { //Record must not be in the database already
                recordHolder_add(rh, cap_getName(cap), terminalAdjacencyWriteFn(cap));
                /*assert(!recordHolder_get(rh, cap_getName(cap)));
                void *data = compress(terminalAdjacencyWriteFn(cap), &recordSize);
                assert(!stCache_containsRecord(cache, cap_getName(cap), 0, INT64_MAX));
                stCache_setRecord(cache, cap_getName(cap), 0, recordSize, data);
                free(data);*/
            }
            if ((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
                break;
            }
            Segment *segment = cap_getSegment(adjacentCap);
            recordHolder_add(rh, segment_getName(segment), segmentWriteFn(segment));
            /*assert(!stCache_containsRecord(cache, segment_getName(segment), 0, INT64_MAX));
            void *data = compress(segmentWriteFn(segment), &recordSize);
            stCache_setRecord(cache, segment_getName(segment), 0, recordSize, data);
            free(data);*/
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

static void cacheNestedRecords(stKVDatabase *database, RecordHolder *rh, stList *caps) {
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
        recordHolder_add(rh, *recordName, stString_copy(record));
        /*assert(!stCache_containsRecord(cache, *recordName, 0, INT64_MAX));
        stCache_setRecord(cache, *recordName, 0, recordSize, record);*/
        stKVDatabaseBulkResult_destruct(result); //Cleanup the memory as we go.
        free(recordName);
    }
    assert(stList_length(getRequests) == 0);
    stList_destruct(getRequests);
    stList_destruct(records);
}

static RecordHolder *cacheRecords(stKVDatabase *database, stList *caps, char *(*segmentWriteFn)(Segment *),
        char *(*terminalAdjacencyWriteFn)(Cap *)) {
    /*
     * Cache all the elements needed to construct the set of threads.
     */
    RecordHolder *rh = recordHolder_construct(); //stCache_construct();
    cacheNestedRecords(database, rh, caps);
    cacheNonNestedRecords(rh, caps, segmentWriteFn, terminalAdjacencyWriteFn);
    return rh;
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

static char *getThread(RecordHolder *rh, Cap *startCap, bool deleteUsedRecords) {
    /*
     * Iterate through, first calculating the length of the final record, then concatenating the results.
     */
    Cap *cap = startCap;
    stList *strings = stList_construct3(0, free);
    while (1) { //Calculate the size of the entries in the DB that represent the thread.
        char *s = recordHolder_remove(rh, cap_getName(cap));
        assert(s != NULL);
        stList_append(strings, s);

        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);

        /*int64_t recordSize;
        assert(stCache_containsRecord(cache, cap_getName(cap), 0, INT64_MAX));
        void *data = stCache_getRecord(cache, cap_getName(cap), 0, INT64_MAX, &recordSize);
        stList_append(strings, decompress(data, recordSize));
        if(deleteUsedRecords) {
            recordHolder_delete(rh, cap_getName(cap));
            //stCache_delete(cache, cap_getName(cap));
        }*/
        if ((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
            break;
        }
        s = recordHolder_remove(rh, segment_getName(cap_getSegment(adjacentCap)));
        assert(s != NULL);
        stList_append(strings, s);

        /*assert(stCache_containsRecord(cache, segment_getName(cap_getSegment(adjacentCap)), 0, INT64_MAX));
        data = stCache_getRecord(cache, segment_getName(cap_getSegment(adjacentCap)), 0, INT64_MAX, &recordSize);
        stList_append(strings, decompress(data, recordSize));
        if(deleteUsedRecords) {
            //stCache_deleteRecord(cache, segment_getName(cap_getSegment(adjacentCap)));
        }*/
    }
    char *string = stString_join2("", strings);
    stList_destruct(strings);
    return string;
}

void buildRecursiveThreads(stKVDatabase *database, stList *caps, char *(*segmentWriteFn)(Segment *),
        char *(*terminalAdjacencyWriteFn)(Cap *)) {
    //Cache records
    RecordHolder *rh = cacheRecords(database, caps, segmentWriteFn, terminalAdjacencyWriteFn);

    //Build new threads
    stList *records = stList_construct3(stList_length(caps), (void(*)(void *)) stKVDatabaseBulkRequest_destruct);
//#pragma omp parallel for
    for (int64_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        char *string = getThread(rh, cap, 0);
        assert(string != NULL);
        //int64_t recordSize;
        //void *data = compress(string, &recordSize);
        stList_set(records, i, stKVDatabaseBulkRequest_constructInsertRequest(cap_getName(cap),
                                                                              string, sizeof(char)*(strlen(string)+1))); //data, recordSize));
        //free(data);
        free(string);
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
    recordHolder_destruct(rh);
    //stCache_destruct(cache);
    stList_destruct(records);
}

stList *buildRecursiveThreadsInListP(RecordHolder *rh, stList *caps, char *(*segmentWriteFn)(Segment *),
                                    char *(*terminalAdjacencyWriteFn)(Cap *), bool deleteUsedRecords) {
    //Build new threads
    stList *threadStrings = stList_construct3(stList_length(caps), free);
//#pragma omp parallel for
    for (int64_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        stList_set(threadStrings, i, getThread(rh, cap, deleteUsedRecords));
    }
    return threadStrings;
}

stList *buildRecursiveThreadsInList(stKVDatabase *database, stList *caps, char *(*segmentWriteFn)(Segment *),
        char *(*terminalAdjacencyWriteFn)(Cap *)) {
    //Cache records
    RecordHolder *rh = cacheRecords(database, caps, segmentWriteFn, terminalAdjacencyWriteFn);
    stList *threadStrings = buildRecursiveThreadsInListP(rh, caps, segmentWriteFn, terminalAdjacencyWriteFn, 0);
    recordHolder_destruct(rh);
    return threadStrings;
}

void buildRecursiveThreadsNoDb(RecordHolder *rh, stList *caps, char *(*segmentWriteFn)(Segment *),
                           char *(*terminalAdjacencyWriteFn)(Cap *)) {
    //Cache records
    cacheNonNestedRecords(rh, caps, segmentWriteFn, terminalAdjacencyWriteFn);

    //Build new threads and add to cache
//#pragma omp parallel for
    for (int64_t i = 0; i < stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        char *string = getThread(rh, cap, 1);
        assert(string != NULL);
        /*int64_t recordSize;
        void *data = compress(string, &recordSize); // this frees the string*/
        recordHolder_add(rh, cap_getName(cap), string);
    }
}

stList *buildRecursiveThreadsInListNoDb(RecordHolder *rh, stList *caps, char *(*segmentWriteFn)(Segment *),
                                    char *(*terminalAdjacencyWriteFn)(Cap *)) {
    cacheNonNestedRecords(rh, caps, segmentWriteFn, terminalAdjacencyWriteFn);
    return buildRecursiveThreadsInListP(rh, caps, segmentWriteFn, terminalAdjacencyWriteFn, 1);
}


