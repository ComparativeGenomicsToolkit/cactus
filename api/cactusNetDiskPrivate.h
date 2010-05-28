#ifndef CACTUS_NET_DISK_PRIVATE_H_
#define CACTUS_NET_DISK_PRIVATE_H_

#include "cactusGlobals.h"

#define NETDISK_NAME_INCREMENT 1000000

struct _netDisk {
	char *stringFile;
	int64_t stringFileLength;
	char *iDDatabaseName;
	char *netsDatabaseName;
	char *metaDataDatabaseName;
	TCBDB *netsDatabase;
	TCBDB *metaDataDatabase;
	TCBDB *iDDatabase;
	st_SortedSet *metaSequences;
	st_SortedSet *metaEvents;
	st_SortedSet *nets;
	Name uniqueNumber;
	Name maxUniqueNumber;
};


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Net disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Adds a newly constructed net to the memory of the netDisk.
 */
void netDisk_addNet(NetDisk *netDisk, Net *net);

/*
 * Removes the net from the disk (if it is on disk);
 */
void netDisk_deleteNetFromDisk(NetDisk *netDisk, Name netName);

/*
 * Registers the net is being freed from memory.
 */
void netDisk_unloadNet(NetDisk *netDisk, Net *net);

/*
 * Gets the first net in the list of nets in memory, or returns NULL if the list is empty.
 */
Net *netDisk_getFirstNetInMemory(NetDisk *netDisk);

/*
 * Returns a net in memory, or NULL, if not in memory.
 */
Net *netDisk_getNetInMemory(NetDisk *netDisk, Name netName);

/*
 * Functions on meta sequences.
 */

/*
 * Adds a newly constructed meta sequence to the memory of the netDisk.
 */
void netDisk_addMetaSequence(NetDisk *netDisk, MetaSequence *metaSequence);

/*
 * Removes the meta sequence from the disk.
 */
void netDisk_deleteMetaSequenceFromDisk(NetDisk *netDisk, Name metaSequenceName);

/*
 * Registers the meta sequence is being freed from memory.
 */
void netDisk_unloadMetaSequence(NetDisk *netDisk, MetaSequence *metaSequence);

/*
 * Gets a meta sequence from the set held in memory.
 */
MetaSequence *netDisk_getMetaSequenceInMemory(NetDisk *netDisk, Name metaSequenceName);

/*
 * Gets the first meta sequence in memory, or NULL, of none remaining.
 */
MetaSequence *netDisk_getFirstMetaSequenceInMemory(NetDisk *netDisk);

/*
 * Functions on meta events.
 */

/*
 * Adds a newly constructed meta sequence to the memory of the netDisk.
 */
void netDisk_addMetaEvent(NetDisk *netDisk, MetaEvent *metaEvent);

/*
 * Deletes the meta event from the disk.
 */
void netDisk_deleteMetaEventFromDisk(NetDisk *netDisk, Name metaEventName);

/*
 * Registers the meta event is being freed from memory.
 */
void netDisk_unloadMetaEvent(NetDisk *netDisk, MetaEvent *metaEvent);

/*
 * Gets a meta event from the pool of meta events in memory, or returns null.
 */
MetaEvent *netDisk_getMetaEventInMemory(NetDisk *netDisk, Name metaEventName);

/*
 * Get first meta event in memory, or NULL if none remaining.
 */
MetaEvent *netDisk_getFirstMetaEventInMemory(NetDisk *netDisk);

/*
 * Functions on strings stored by the net disk.
 */

/*
 * Adds the sequence string to the bucket of sequence.
 *
 * This function is NOT thread safe.
 */
int64_t netDisk_addString(NetDisk *netDisk, const char *string, int32_t length);

/*
 * Retrieves a string from the bucket of sequence.
 */
char *netDisk_getString(NetDisk *netDisk, int64_t offset, int32_t start, int32_t length, int32_t strand);


#endif
