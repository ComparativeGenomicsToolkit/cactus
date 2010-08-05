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
Flower *cactusDisk_getNet(CactusDisk *cactusDisk, Name flowerName);

/*
 * Returns the total number of flowers on disk.
 */
int32_t cactusDisk_getNetNumberOnDisk(CactusDisk *cactusDisk);

/*
 * Gets an iterator to iterate through all the flower names of flowers on disk.
 */
CactusDisk_NetNameIterator *cactusDisk_getNetNamesOnDiskIterator(CactusDisk *cactusDisk);

/*
 * Gets the next flower name from the iterator, returns 'NULL_NAME' when done.
 */
Name cactusDisk_getNextNetName(CactusDisk_NetNameIterator *flowerIterator);

/*
 * Destructs the iterator.
 */
void cactusDisk_destructNetNamesOnDiskIterator(CactusDisk_NetNameIterator *flowerIterator);

/*
 * Returns the number of flowers currently in memory.
 */
int32_t cactusDisk_getNetNumberInMemory(CactusDisk *cactusDisk);


/*
 * Gets an iterator to iterate through the flowers currently in memory.
 */
CactusDisk_NetIterator *cactusDisk_getNetsInMemoryIterator(CactusDisk *cactusDisk);

/*
 * Gets the next flower from the iterator.
 */
Flower *cactusDisk_getNextNet(CactusDisk_NetIterator *flowerIterator);

/*
 * Gets the previous flower from the iterator.
 */
Flower *cactusDisk_getPreviousNet(CactusDisk_NetIterator *flowerIterator);

/*
 * Duplicates the iterator.
 */
CactusDisk_NetIterator *cactusDisk_copyNetIterator(CactusDisk_NetIterator *flowerIterator);

/*
 * Destructs the iterator.
 */
void cactusDisk_destructNetsInMemoryIterator(CactusDisk_NetIterator *flowerIterator);

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
