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
 * Constructs a cactus disk to hold the flower hierarchy.
 */
CactusDisk *cactusDisk_construct();

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
 * Gets a flower the cactusDisk contains. If the flower is not in memory it will be loaded. If not in memory or on disk, returns NULL.
 */
Flower *cactusDisk_getFlower(CactusDisk *cactusDisk, Name flowerName);

/*
 * Gets a sequence
 */
Sequence *cactusDisk_getSequence(CactusDisk *cactusDisk, Name sequenceName);

/*
 * Get the event tree.
 */
EventTree *cactusDisk_getEventTree(CactusDisk *cactusDisk);

#endif
