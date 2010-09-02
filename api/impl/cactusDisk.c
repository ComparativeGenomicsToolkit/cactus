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

CactusDisk *cactusDisk_construct(stKVDatabaseConf *conf, bool create) {
    CactusDisk *cactusDisk;
    cactusDisk = st_malloc(sizeof(CactusDisk));

    //construct lists of in memory objects
    cactusDisk->metaSequences = stSortedSet_construct3(
            cactusDisk_constructMetaSequencesP, NULL);
    cactusDisk->flowers = stSortedSet_construct3(cactusDisk_constructFlowersP, NULL);
    cactusDisk->flowerNamesMarkedForDeletion = stSortedSet_construct2(free);

    //Now open the database
    cactusDisk->database = stKVDatabase_construct(conf, create);

    //initialise the unique ids.
    //cactusDisk_getUniqueID(cactusDisk);
    st_randomSeed(time(NULL));
    cactusDisk->uniqueNumber = 0;
    cactusDisk->maxUniqueNumber = 0;

    return cactusDisk;
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

    free(cactusDisk);
}

void cactusDisk_write(CactusDisk *cactusDisk) {
    void *vA;
    int32_t recordSize;
    Flower *flower;
    MetaSequence *metaSequence;

    stKVDatabase_startTransaction(cactusDisk->database);

    stSortedSetIterator *it = stSortedSet_getIterator(cactusDisk->flowers);
    while ((flower = stSortedSet_getNext(it)) != NULL) {
        vA = binaryRepresentation_makeBinaryRepresentation(flower,
                (void(*)(void *, void(*)(const void * ptr, size_t size,
                        size_t count))) flower_writeBinaryRepresentation,
                &recordSize);
        if(stKVDatabase_getRecord(cactusDisk->database, flower_getName(flower)) != NULL) {
            stKVDatabase_updateRecord(cactusDisk->database, flower_getName(flower),
                    vA, recordSize);
        }
        else {
            stKVDatabase_insertRecord(cactusDisk->database, flower_getName(flower),
                                vA, recordSize);
        }
        free(vA);
    }
    stSortedSet_destructIterator(it);

    it = stSortedSet_getIterator(cactusDisk->metaSequences);
    while ((metaSequence = stSortedSet_getNext(it)) != NULL) {
        vA = binaryRepresentation_makeBinaryRepresentation(metaSequence,
                (void(*)(void *, void(*)(const void * ptr, size_t size,
                        size_t count))) metaSequence_writeBinaryRepresentation,
                &recordSize);
        if(stKVDatabase_getRecord(cactusDisk->database, metaSequence_getName(metaSequence)) != NULL) {
            stKVDatabase_updateRecord(cactusDisk->database, metaSequence_getName(metaSequence),
                    vA, recordSize);
        }
        else {
            stKVDatabase_insertRecord(cactusDisk->database, metaSequence_getName(metaSequence),
                                vA, recordSize);
        }
        free(vA);
    }
    stSortedSet_destructIterator(it);

    //Remove nets that are marked for deletion..
    it = stSortedSet_getIterator(cactusDisk->flowerNamesMarkedForDeletion);
    char *nameString;
    while ((nameString = stSortedSet_getNext(it)) != NULL) {
        Name name = cactusMisc_stringToName(nameString);
        if((vA = stKVDatabase_getRecord(cactusDisk->database, name)) != NULL) {
            stKVDatabase_removeRecord(cactusDisk->database, name);
            free(vA);
        }
    }
    stSortedSet_destructIterator(it);

    stKVDatabase_commitTransaction(cactusDisk->database);
}

Flower *cactusDisk_getFlower(CactusDisk *cactusDisk, Name flowerName) {
    static Flower flower;
    flower.name = flowerName;
    Flower *flower2;
    if ((flower2 = stSortedSet_search(cactusDisk->flowers, &flower)) != NULL) {
        return flower2;
    }
    void *cA = stKVDatabase_getRecord(cactusDisk->database, flowerName);
    if (cA == NULL) {
        return NULL;
    }
    void *cA2 = cA;
    flower2 = flower_loadFromBinaryRepresentation(&cA2, cactusDisk);
    free(cA);
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
    void *cA = stKVDatabase_getRecord(cactusDisk->database, metaSequenceName);
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

void cactusDisk_addFlower(CactusDisk *cactusDisk, Flower *flower) {
    assert(stSortedSet_search(cactusDisk->flowers, flower) == NULL);
    stSortedSet_insert(cactusDisk->flowers, flower);
}

void cactusDisk_removeFlower(CactusDisk *cactusDisk, Flower *flower) {
    assert(stSortedSet_search(cactusDisk->flowers, flower) != NULL);
    stSortedSet_remove(cactusDisk->flowers, flower);
}

void cactusDisk_deleteFlowerFromDisk(CactusDisk *cactusDisk, Flower *flower) {
    char *nameString = cactusMisc_nameToString(flower_getName(flower));
    if(stSortedSet_search(cactusDisk->flowerNamesMarkedForDeletion, nameString) == NULL) {
        stSortedSet_insert(cactusDisk->flowerNamesMarkedForDeletion, nameString);
    }
    else {
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
    Name name = cactusDisk_getUniqueID(cactusDisk);
    stKVDatabase_startTransaction(cactusDisk->database);
    stKVDatabase_insertRecord(cactusDisk->database, name, string, (strlen(string)+1)*sizeof(char));
    stKVDatabase_commitTransaction(cactusDisk->database);
    return name;
}

char *cactusDisk_getString(CactusDisk *cactusDisk, Name name,
        int32_t start, int32_t length, int32_t strand) {
    char *string = stKVDatabase_getPartialRecord(cactusDisk->database, name, start*sizeof(char), (length+1)*sizeof(char));
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

#define CACTUS_DISK_NAME_INCREMENT 16384
#define CACTUS_DISK_BUCKET_NUMBER 65536

void cactusDisk_getBlockOfUniqueIDs(CactusDisk *cactusDisk) {
    stKVDatabase_startTransaction(cactusDisk->database);
    Name keyName = st_randomInt(-CACTUS_DISK_BUCKET_NUMBER, 0);
    assert(keyName >= -CACTUS_DISK_BUCKET_NUMBER);
    assert(keyName < 0);
    void *vA = stKVDatabase_getRecord(cactusDisk->database, keyName);
    int64_t bucketSize = INT64_MAX / CACTUS_DISK_BUCKET_NUMBER;
    int64_t minimumValue = bucketSize * (abs(keyName)-1) + 1; //plus one for the reserved '0' value.
    int64_t maximumValue = minimumValue + (bucketSize - 1);
    bool recordExists = 0;
    if (vA == NULL) {
        cactusDisk->uniqueNumber = minimumValue;
    } else {
        recordExists = 1;
        cactusDisk->uniqueNumber = *((Name *) vA);
    }
    cactusDisk->maxUniqueNumber = cactusDisk->uniqueNumber
            + CACTUS_DISK_NAME_INCREMENT;
    if(cactusDisk->maxUniqueNumber >= maximumValue) {
        st_errAbort("We have exhausted a bucket, which seems really unlikely\n");
    }
    if(recordExists) {
        stKVDatabase_updateRecord(cactusDisk->database, keyName,
                &cactusDisk->maxUniqueNumber, sizeof(Name));
    }
    else {
        stKVDatabase_insertRecord(cactusDisk->database, keyName,
                        &cactusDisk->maxUniqueNumber, sizeof(Name));
    }
    stKVDatabase_commitTransaction(cactusDisk->database);
    free(vA);
}

int64_t cactusDisk_getUniqueID(CactusDisk *cactusDisk) {
    assert(cactusDisk->uniqueNumber <= cactusDisk->maxUniqueNumber);
    if (cactusDisk->uniqueNumber == cactusDisk->maxUniqueNumber) {
        cactusDisk_getBlockOfUniqueIDs(cactusDisk);
    }
    return cactusDisk->uniqueNumber++;
}
