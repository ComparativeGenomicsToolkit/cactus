#ifndef CACTUS_DISK_PRIVATE_H_
#define CACTUS_DISK_PRIVATE_H_

#include "cactusGlobals.h"

#define CACTUSDISK_NAME_INCREMENT 1000000

struct _cactusDisk {
	char *stringFile;
	int64_t stringFileLength;
	char *iDDatabaseName;
	char *netsDatabaseName;
	char *metaDataDatabaseName;
	TCBDB *netsDatabase;
	TCBDB *metaDataDatabase;
	TCBDB *iDDatabase;
	stSortedSet *metaSequences;
	stSortedSet *metaEvents;
	stSortedSet *nets;
	Name uniqueNumber;
	Name maxUniqueNumber;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Cactus disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Adds a newly constructed net to the memory of the cactusDisk.
 */
void cactusDisk_addNet(CactusDisk *cactusDisk, Net *net);

/*
 * Removes the net from the disk (if it is on disk);
 */
void cactusDisk_deleteNetFromDisk(CactusDisk *cactusDisk, Name netName);

/*
 * Registers the net is being freed from memory.
 */
void cactusDisk_unloadNet(CactusDisk *cactusDisk, Net *net);

/*
 * Gets the first net in the list of nets in memory, or returns NULL if the list is empty.
 */
Net *cactusDisk_getFirstNetInMemory(CactusDisk *cactusDisk);

/*
 * Returns a net in memory, or NULL, if not in memory.
 */
Net *cactusDisk_getNetInMemory(CactusDisk *cactusDisk, Name netName);

/*
 * Functions on meta sequences.
 */

/*
 * Adds a newly constructed meta sequence to the memory of the cactusDisk.
 */
void cactusDisk_addMetaSequence(CactusDisk *cactusDisk, MetaSequence *metaSequence);

/*
 * Removes the meta sequence from the disk.
 */
void cactusDisk_deleteMetaSequenceFromDisk(CactusDisk *cactusDisk, Name metaSequenceName);

/*
 * Registers the meta sequence is being freed from memory.
 */
void cactusDisk_unloadMetaSequence(CactusDisk *cactusDisk, MetaSequence *metaSequence);

/*
 * Gets a meta sequence from the set held in memory.
 */
MetaSequence *cactusDisk_getMetaSequenceInMemory(CactusDisk *cactusDisk, Name metaSequenceName);

/*
 * Gets the first meta sequence in memory, or NULL, of none remaining.
 */
MetaSequence *cactusDisk_getFirstMetaSequenceInMemory(CactusDisk *cactusDisk);

/*
 * Functions on meta events.
 */

/*
 * Adds a newly constructed meta sequence to the memory of the cactusDisk.
 */
void cactusDisk_addMetaEvent(CactusDisk *cactusDisk, MetaEvent *metaEvent);

/*
 * Deletes the meta event from the disk.
 */
void cactusDisk_deleteMetaEventFromDisk(CactusDisk *cactusDisk, Name metaEventName);

/*
 * Registers the meta event is being freed from memory.
 */
void cactusDisk_unloadMetaEvent(CactusDisk *cactusDisk, MetaEvent *metaEvent);

/*
 * Gets a meta event from the pool of meta events in memory, or returns null.
 */
MetaEvent *cactusDisk_getMetaEventInMemory(CactusDisk *cactusDisk, Name metaEventName);

/*
 * Get first meta event in memory, or NULL if none remaining.
 */
MetaEvent *cactusDisk_getFirstMetaEventInMemory(CactusDisk *cactusDisk);

/*
 * Functions on strings stored by the cactus disk.
 */

/*
 * Adds the sequence string to the bucket of sequence.
 *
 * This function is NOT thread safe.
 */
int64_t cactusDisk_addString(CactusDisk *cactusDisk, const char *string, int32_t length);

/*
 * Retrieves a string from the bucket of sequence.
 */
char *cactusDisk_getString(CactusDisk *cactusDisk, int64_t offset, int32_t start, int32_t length, int32_t strand);


#endif
