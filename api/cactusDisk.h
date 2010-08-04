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
 * Constructs a cactus disk to load nets.
 */
CactusDisk *cactusDisk_construct(const char *cactusDiskFile);

/*
 * Destructs the cactus disk, and all open nets and sequences.
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
 * Gets a net the cactusDisk contains. If the net is not in memory it will be loaded. If not in memory or on disk, returns NULL.
 */
Net *cactusDisk_getNet(CactusDisk *cactusDisk, Name netName);

/*
 * Returns the total number of nets on disk.
 */
int32_t cactusDisk_getNetNumberOnDisk(CactusDisk *cactusDisk);

/*
 * Gets an iterator to iterate through all the net names of nets on disk.
 */
CactusDisk_NetNameIterator *cactusDisk_getNetNamesOnDiskIterator(CactusDisk *cactusDisk);

/*
 * Gets the next net name from the iterator, returns 'NULL_NAME' when done.
 */
Name cactusDisk_getNextNetName(CactusDisk_NetNameIterator *netIterator);

/*
 * Destructs the iterator.
 */
void cactusDisk_destructNetNamesOnDiskIterator(CactusDisk_NetNameIterator *netIterator);

/*
 * Returns the number of nets currently in memory.
 */
int32_t cactusDisk_getNetNumberInMemory(CactusDisk *cactusDisk);


/*
 * Gets an iterator to iterate through the nets currently in memory.
 */
CactusDisk_NetIterator *cactusDisk_getNetsInMemoryIterator(CactusDisk *cactusDisk);

/*
 * Gets the next net from the iterator.
 */
Net *cactusDisk_getNextNet(CactusDisk_NetIterator *netIterator);

/*
 * Gets the previous net from the iterator.
 */
Net *cactusDisk_getPreviousNet(CactusDisk_NetIterator *netIterator);

/*
 * Duplicates the iterator.
 */
CactusDisk_NetIterator *cactusDisk_copyNetIterator(CactusDisk_NetIterator *netIterator);

/*
 * Destructs the iterator.
 */
void cactusDisk_destructNetsInMemoryIterator(CactusDisk_NetIterator *netIterator);

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
