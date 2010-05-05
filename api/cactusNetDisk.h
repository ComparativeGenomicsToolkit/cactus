#ifndef CACTUS_NET_DISK_H_
#define CACTUS_NET_DISK_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a net disk to load nets.
 */
NetDisk *netDisk_construct(const char *netDiskFile);

/*
 * Destructs the net disk, and all open nets and sequences.
 */
void netDisk_destruct(NetDisk *netDisk);

/*
 * Retrieves the next unique ID.
 */
int64_t netDisk_getUniqueID(NetDisk *netDisk);

/*
 * Writes the updated state of the parts of the net disk in memory to disk.
 *
 */
void netDisk_write(NetDisk *netDisk);

/*
 * Gets a net the netDisk contains. If the net is not in memory it will be loaded. If not in memory or on disk, returns NULL.
 */
Net *netDisk_getNet(NetDisk *netDisk, Name netName);

/*
 * Returns the total number of nets on disk.
 */
int32_t netDisk_getNetNumberOnDisk(NetDisk *netDisk);

/*
 * Gets an iterator to iterate through all the net names of nets on disk.
 */
NetDisk_NetNameIterator *netDisk_getNetNamesOnDiskIterator(NetDisk *netDisk);

/*
 * Gets the next net name from the iterator, returns 'NULL_NAME' when done.
 */
Name netDisk_getNextNetName(NetDisk_NetNameIterator *netIterator);

/*
 * Destructs the iterator.
 */
void netDisk_destructNetNamesOnDiskIterator(NetDisk_NetNameIterator *netIterator);

/*
 * Returns the number of nets currently in memory.
 */
int32_t netDisk_getNetNumberInMemory(NetDisk *netDisk);


/*
 * Gets an iterator to iterate through the nets currently in memory.
 */
NetDisk_NetIterator *netDisk_getNetsInMemoryIterator(NetDisk *netDisk);

/*
 * Gets the next net from the iterator.
 */
Net *netDisk_getNextNet(NetDisk_NetIterator *netIterator);

/*
 * Gets the previous net from the iterator.
 */
Net *netDisk_getPreviousNet(NetDisk_NetIterator *netIterator);

/*
 * Duplicates the iterator.
 */
NetDisk_NetIterator *netDisk_copyNetIterator(NetDisk_NetIterator *netIterator);

/*
 * Destructs the iterator.
 */
void netDisk_destructNetsInMemoryIterator(NetDisk_NetIterator *netIterator);

/*
 * Functions on meta sequences.
 */

/*
 * Gets the meta sequence for an object.
 */
MetaSequence *netDisk_getMetaSequence(NetDisk *netDisk, Name metaSequenceName);

/*
 * Functions on meta events.
 */

/*
 * Gets a meta event from either disk or memory.
 */
MetaEvent *netDisk_getMetaEvent(NetDisk *netDisk, Name metaEventName);

#endif
