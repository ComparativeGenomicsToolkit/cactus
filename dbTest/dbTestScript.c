#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "sonLib.h"

void usage() {
    fprintf(stderr, "dpTestScript, version 0.1\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --databaseConf : The database connection script\n");
    fprintf(stderr, "-c --firstKey : First key.\n");
    fprintf(stderr, "-d --keyNumber : Total number of keys.\n");
    fprintf(stderr, "-e --addRecords : Add records instead of getting them.\n");
    fprintf(stderr, "-f --setRecords : After adding/getting records, set the records.\n");
    fprintf(stderr, "-g --minRecordSize : Min size of record.\n");
    fprintf(stderr, "-h --maxRecordSize : Min size of record.\n");
    fprintf(stderr, "-i --create : Make the database.\n");
}

void *getRandomRecord(int64_t minRecordSize, int64_t maxRecordSize, int64_t *recordSize) {
    *recordSize = maxRecordSize;
    char *cA = st_malloc(sizeof(char) * (*recordSize));
    for (int64_t i = 0; i < *recordSize; i++) {
        cA[i] = st_randomInt(0, 128);
    }
    return cA;
}

int main(int argc, char *argv[]) {
    /*
     * Script for adding a reference genome to a flower.
     */

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * databaseString = NULL;
    int64_t firstKey = INT64_MIN;
    int64_t keyNumber = INT64_MIN;
    bool addRecords = 0, setRecords = 0, create = 0;
    int64_t minRecordSize = 0, maxRecordSize = 100000000;
    int32_t i;

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "databaseConf", required_argument, 0, 'b' }, {
                "firstKey", required_argument, 0, 'c' }, { "keyNumber", required_argument, 0, 'd' }, { "addRecords", no_argument, 0, 'e' },
                { "setRecords", no_argument, 0, 'f' }, { "minRecordSize", required_argument, 0, 'g' }, { "maxRecordSize",
                        required_argument, 0, 'h' }, { "maxRecordSize", no_argument, 0, 'i' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:efg:h:i", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                databaseString = stString_copy(optarg);
                break;
            case 'c':
                i = sscanf(optarg, "%" PRIi64 "", &firstKey);
                if (i != 1) {
                    st_errAbort("Did not parse a valid firstKey number: %s", optarg);
                }
                break;
            case 'd':
                i = sscanf(optarg, "%" PRIi64 "", &keyNumber);
                if (i != 1) {
                    st_errAbort("Did not parse a valid keyNumber: %s", optarg);
                }
                break;
            case 'e':
                addRecords = 1;
                break;
            case 'f':
                setRecords = 1;
                break;
            case 'g':
                i = sscanf(optarg, "%" PRIi64 "", &minRecordSize);
                if (i != 1) {
                    st_errAbort("Did not parse a valid minRecordSize number: %s", optarg);
                }
                break;
            case 'h':
                i = sscanf(optarg, "%" PRIi64 "", &maxRecordSize);
                if (i != 1) {
                    st_errAbort("Did not parse a valid maxRecordSize number: %s", optarg);
                }
                break;
            case 'i':
                create = 1;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(databaseString != NULL);
    assert(firstKey != INT64_MIN);
    assert(keyNumber != INT64_MIN);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(databaseString);
    stKVDatabase *database = stKVDatabase_construct(kvDatabaseConf, create);
    st_logInfo("Set up the database\n");

    //////////////////////////////////////////////
    //Seed the random number generator
    //////////////////////////////////////////////

    int32_t seed = (time(NULL) << 16) | (getpid() & 65535); //Likely to be unique
    st_logDebug("Seeding the random number generator with the value %i\n", seed);
    st_randomSeed(seed);

    ///////////////////////////////////////////////////////////////////////////
    // Do meat of manipulating the database
    ///////////////////////////////////////////////////////////////////////////

    st_logInfo("Modifying records to the database in the range %" PRIi64 " to %" PRIi64 "\n", firstKey, firstKey+keyNumber);

    if (addRecords) {
        st_logInfo("Adding records to the database\n");
        stList *recordNames = stList_construct3(0, free);
        for (int64_t i = 0; i < keyNumber; i++) {
            int64_t recordSize;
            void *vA = getRandomRecord(minRecordSize, maxRecordSize, &recordSize);
            stList_append(recordNames, stKVDatabaseBulkRequest_constructInsertRequest(firstKey + i, vA, recordSize));
        } stTry {
            stKVDatabase_bulkSetRecords(database, recordNames);
        }
        stCatch(except)
        {
            stThrowNewCause(
                    except,
                    ST_KV_DATABASE_EXCEPTION_ID,
                    "An unknown database error occurred when getting a bulk add");
        }
        stTryEnd;
        stList_destruct(recordNames);
    } else {
        st_logInfo("Getting records from the database\n");
        stList *recordNames = stList_construct3(0, free);
        for (int64_t i = 0; i < keyNumber; i++) {
            int64_t *iA = st_malloc(sizeof(int64_t));
            iA[0] = firstKey + i;
            stList_append(recordNames, iA);
        } stTry {
            stList_destruct(stKVDatabase_bulkGetRecords(database, recordNames));
        }
        stCatch(except)
        {
            stThrowNewCause(
                    except,
                    ST_KV_DATABASE_EXCEPTION_ID,
                    "An unknown database error occurred when getting a bulk set");
        }
        stTryEnd;
        stList_destruct(recordNames);
    }

    if (setRecords) {
        st_logInfo("Setting records in the database\n");
        stList *recordNames = stList_construct3(0, free);
        for (int64_t i = 0; i < keyNumber; i++) {
            int64_t recordSize;
            void *vA = getRandomRecord(minRecordSize, maxRecordSize, &recordSize);
            stList_append(recordNames, stKVDatabaseBulkRequest_constructUpdateRequest(firstKey + i, vA, recordSize));
        } stTry {
            stKVDatabase_bulkSetRecords(database, recordNames);
        }
        stCatch(except)
        {
            stThrowNewCause(
                    except,
                    ST_KV_DATABASE_EXCEPTION_ID,
                    "An unknown database error occurred when getting a bulk add");
        }
        stTryEnd;
        stList_destruct(recordNames);
    }

    st_logDebug("Finished!\n");

    ///////////////////////////////////////////////////////////////////////////
    // Close the database down
    ///////////////////////////////////////////////////////////////////////////

    stKVDatabase_destruct(database);
}
