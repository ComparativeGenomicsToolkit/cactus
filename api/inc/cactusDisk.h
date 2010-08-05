#ifndef CACTUS_DISK_H_
#define CACTUS_DISK_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic cactus disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a cactus disk to load flowers.
 */
CactusDisk *cactusDisk_construct(const char *cactusDiskFile);

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
 * Returns the total number of flowers on disk.
 */
int32_t cactusDisk_getFlowerNumberOnDisk(CactusDisk *cactusDisk);

/*
 * Gets an iterator to iterate through all the flower names of flowers on disk.
 */
CactusDisk_FlowerNameIterator *cactusDisk_getFlowerNamesOnDiskIterator(CactusDisk *cactusDisk);

/*
 * Gets the next flower name from the iterator, returns 'NULL_NAME' when done.
 */
Name cactusDisk_getNextFlowerName(CactusDisk_FlowerNameIterator *flowerIterator);

/*
 * Destructs the iterator.
 */
void cactusDisk_destructFlowerNamesOnDiskIterator(CactusDisk_FlowerNameIterator *flowerIterator);

/*
 * Returns the number of flowers currently in memory.
 */
int32_t cactusDisk_getFlowerNumberInMemory(CactusDisk *cactusDisk);


/*
 * Gets an iterator to iterate through the flowers currently in memory.
 */
CactusDisk_FlowerIterator *cactusDisk_getFlowersInMemoryIterator(CactusDisk *cactusDisk);

/*
 * Gets the next flower from the iterator.
 */
Flower *cactusDisk_getNextFlower(CactusDisk_FlowerIterator *flowerIterator);

/*
 * Gets the previous flower from the iterator.
 */
Flower *cactusDisk_getPreviousFlower(CactusDisk_FlowerIterator *flowerIterator);

/*
 * Duplicates the iterator.
 */
CactusDisk_FlowerIterator *cactusDisk_copyFlowerIterator(CactusDisk_FlowerIterator *flowerIterator);

/*
 * Destructs the iterator.
 */
void cactusDisk_destructFlowersInMemoryIterator(CactusDisk_FlowerIterator *flowerIterator);

/*
 * Functions on meta sequences.
 */

/*
 * Gets the meta sequence for an object.
 */
MetaSequence *cactusDisk_getMetaSequence(CactusDisk *cactusDisk, Name metaSequenceName);

/*
 * Functions on meta events.
 */

/*
 * Gets a meta event from either disk or memory.
 */
MetaEvent *cactusDisk_getMetaEvent(CactusDisk *cactusDisk, Name metaEventName);

#endif
