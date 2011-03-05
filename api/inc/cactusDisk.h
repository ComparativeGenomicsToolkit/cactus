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
 * As with the cactusDisk_construct, but with the option to pre-cache all the sequence contained
 * within a flower when retrieved, to minimise the amount of I/O.
 */
CactusDisk *cactusDisk_construct2(stKVDatabaseConf *conf, bool create, bool preCacheSequences);

/*
 * Create a cactus disk, with the option to store the sequences in a file.
 */
CactusDisk *cactusDisk_construct3(stKVDatabaseConf *conf,
        bool preCacheSequences, const char *sequencesFileName);

/*
 * Destructs the cactus disk, and all open flowers and sequences.
 */
void cactusDisk_destruct(CactusDisk *cactusDisk);

/*
 * Retrieves the next unique ID.
 */
int64_t cactusDisk_getUniqueID(CactusDisk *cactusDisk);

/*
 * Writes the updated state of the parts of the cactus disk in memory to disk.
 *
 */
void cactusDisk_write(CactusDisk *cactusDisk);

/*
 * Gets a flower the cactusDisk contains. If the flower is not in memory it will be loaded. If not in memory or on disk, returns NULL.
 */
Flower *cactusDisk_getFlower(CactusDisk *cactusDisk, Name flowerName);

/*
 * Functions on meta sequences.
 */

/*
 * Gets the meta sequence for an object.
 */
MetaSequence *cactusDisk_getMetaSequence(CactusDisk *cactusDisk, Name metaSequenceName);


void splitFlowers(Flower *flower, CactusDisk *newCactusDisk);

/*
 * Copy the flower and its descendants and associated stuff into the given cactus disk.
 */
void splitFlowers(Flower *flower, CactusDisk *newCactusDisk);

/*
 * Merge the flower and its descendants and associated
 */
void mergeFlowers(Flower *flower, CactusDisk *oldCactusDisk);

#endif
