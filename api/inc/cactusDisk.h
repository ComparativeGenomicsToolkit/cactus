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
 * Deterministic naming inside parallel loops.
 *
 * Names normally come off one shared counter, so an object's name depends on
 * how the threads happened to interleave and differs from run to run.  To avoid
 * that, a parallel loop hands each iteration its own interval of names, chosen
 * from the iteration index rather than from the thread:
 *
 *   Name base = cactusDisk_reserveNames(cactusDisk, total);   // before the loop
 *   #pragma omp parallel for
 *   for (i ...) {
 *       cactusDisk_pushNameInterval(cactusDisk, base + offset[i], size[i]);
 *       ... construct objects ...
 *       cactusDisk_popNameInterval(cactusDisk);
 *   }
 *
 * The intervals must not overlap, and the caller has to size them: sizing is a
 * per-loop matter, since only the caller knows what an iteration can build.
 * Anything an iteration allocates beyond its interval falls back to the shared
 * counter -- correct, but not reproducible, so it is logged as a warning.
 *
 * The interval is held in thread-local state, so it applies to the thread that
 * pushed it and only until it pops.
 */
Name cactusDisk_reserveNames(CactusDisk *cactusDisk, int64_t nameNumber);

void cactusDisk_pushNameInterval(CactusDisk *cactusDisk, Name start, int64_t nameNumber);

void cactusDisk_popNameInterval(CactusDisk *cactusDisk);

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
