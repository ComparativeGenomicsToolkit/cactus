/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_DISK_PRIVATE_H_
#define CACTUS_DISK_PRIVATE_H_

#include "cactusGlobals.h"

struct _cactusDisk {
    stKVDatabase *database;
    stSortedSet *metaSequences;
    stSortedSet *flowers;
    stSortedSet *flowerNamesMarkedForDeletion;
    stList *updateRequests;
    stCache *cache;
    stCache *stringCache;
    Name uniqueNumber;
    Name maxUniqueNumber;
    bool storeSequencesInAFile;
    char *sequencesFileName;
    char *absSequencesFileName;
    FILE *sequencesReadFileHandle;
    FILE *sequencesWriteFileHandle;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Cactus disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Returns non-zero if the given flower is loaded in memory.
 */
bool cactusDisk_flowerIsLoaded(CactusDisk *cactusDisk, Name flowerName);

/*
 * Adds a newly constructed flower to the memory of the cactusDisk.
 */
void cactusDisk_addFlower(CactusDisk *cactusDisk, Flower *flower);

/*
 * Registers the flower is being freed from memory.
 */
void cactusDisk_removeFlower(CactusDisk *cactusDisk, Flower *flower);

/*
 * Registers the flower should be removed from the disk.
 */
void cactusDisk_deleteFlowerFromDisk(CactusDisk *cactusDisk, Flower *flower);

/*
 * Functions on meta sequences.
 */

/*
 * Adds a newly constructed meta sequence to the memory of the cactusDisk.
 */
void cactusDisk_addMetaSequence(CactusDisk *cactusDisk,
        MetaSequence *metaSequence);

/*
 * Registers the meta sequence is being freed from memory.
 */
void cactusDisk_removeMetaSequence(CactusDisk *cactusDisk,
        MetaSequence *metaSequence);

bool cactusDisk_storedInFile(CactusDisk *cactusDisk);



/*
 * Functions on strings stored by the cactus disk.
 */

/*
 * Adds the sequence string to the database.
 */
Name cactusDisk_addString(CactusDisk *cactusDisk, const char *string);

/*
 * Retrieves a string from the bucket of sequence.
 */
char *cactusDisk_getString(CactusDisk *cactusDisk, Name name,
        int64_t start, int64_t length, int64_t strand, int64_t totalSequenceLength);

/*
 * Gets the string from a cache.
 */
char *cactusDisk_getStringFromCache(CactusDisk *cactusDisk, Name name, int64_t start, int64_t length, int64_t strand);

#endif
