/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_DISK_H_
#define CACTUS_DISK_H_

#include "cactusGlobals.h"

// General database exception id
extern const char *CACTUS_DISK_EXCEPTION_ID;

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic cactus disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a cactus disk to load flowers. If 'create' is non-zero the cactus disk is
 * to be created. An exception will be thrown if the 'create' is non-zero and the cactus
 * disk already exists.
 */
CactusDisk *cactusDisk_construct(stKVDatabaseConf *conf, bool create);

/*
 * Create a cactus disk, with the option to store the sequences in a file.
 */
CactusDisk *cactusDisk_construct2(stKVDatabaseConf *conf,
        const char *sequencesFileName);

/*
 * Destructs the cactus disk and all open flowers and sequences, and
 * then disconnects from the cactus DB.
 */
void cactusDisk_destruct(CactusDisk *cactusDisk);

/*
 * Retrieves the next unique ID.
 */
int64_t cactusDisk_getUniqueID(CactusDisk *cactusDisk);

/*
 * Retrieves a contiguous interval of unique ids starting from the return value to return value + intervalSize (exclusive).
 */
int64_t cactusDisk_getUniqueIDInterval(CactusDisk *cactusDisk, int64_t intervalSize);

/*
 * Writes the updated state of the parts of the cactus disk in memory to disk.
 *
 */
void cactusDisk_write(CactusDisk *cactusDisk);

/*
 * This is used to serialise a flower before a call to a cactusDisk_write, it is exposed for use in the cactus_caf code.
 */
void cactusDisk_addUpdateRequest(CactusDisk *cactusDisk, Flower *flower);

/*
 * Gets a flower the cactusDisk contains. If the flower is not in memory it will be loaded. If not in memory or on disk, returns NULL.
 */
Flower *cactusDisk_getFlower(CactusDisk *cactusDisk, Name flowerName);

/*
 * Gets a list of flowers for the given list of flower names.
 */
stList *cactusDisk_getFlowers(CactusDisk *cactusDisk, stList *flowerNames);

/*
 * Functions on meta sequences.
 */

/*
 * Gets the meta sequence for an object.
 */
MetaSequence *cactusDisk_getMetaSequence(CactusDisk *cactusDisk, Name metaSequenceName);

/*
 * Precaches all the sequences in a given set of flowers into the cache.
 */
void cactusDisk_preCacheStrings(CactusDisk *cactusDisk, stList *flowers);

/*
 * Precaches all the segments sequences in the given set of flowers.
 */
void cactusDisk_preCacheSegmentStrings(CactusDisk *cactusDisk, stList *flowers);

#endif
