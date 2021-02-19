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
    stSortedSet *sequences;
    stSortedSet *flowers;
    stSortedSet *flowerNamesMarkedForDeletion;
    stList *updateRequests;
    stCache *cache;
    stCache *stringCache;
    EventTree *eventTree;
    Name uniqueNumber;
    Name maxUniqueNumber;

    // Structures used if cactus disk is to only be held in memory
    bool inMemory;
    stHash *allStrings; // If the strings are being all stored in memory, a map of names to strings
    Name currentName; // Used as a counter for issuing names
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
void cactusDisk_addSequence(CactusDisk *cactusDisk,
        Sequence *sequence);

/*
 * Registers the meta sequence is being freed from memory.
 */
void cactusDisk_removeSequence(CactusDisk *cactusDisk,
        Sequence *sequence);

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

/*
 * Set the event tree for this disk. (Hopefully this only happens once.)
 */
void cactusDisk_setEventTree(CactusDisk *cactusDisk, EventTree *eventTree);

#endif
