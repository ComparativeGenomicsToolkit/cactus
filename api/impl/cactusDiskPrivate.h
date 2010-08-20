#ifndef CACTUS_DISK_PRIVATE_H_
#define CACTUS_DISK_PRIVATE_H_

#include "cactusGlobals.h"

#define CACTUSDISK_NAME_INCREMENT 10000
#define CACTUSDISK_BUCKET_NUMBER 10000

struct _cactusDisk {
    char *stringFile;
    int64_t stringFileLength;
    stKVDatabase *database;
    stSortedSet *metaSequences;
    stSortedSet *flowers;
    stSortedSet *flowerNamesMarkedForDeletion;
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
 * Adds a newly constructed flower to the memory of the cactusDisk.
 */
void cactusDisk_addFlower(CactusDisk *cactusDisk, Flower *flower);

/*
 * Registers the flower is being freed from memory.
 */
void cactusDisk_removeFlower(CactusDisk *cactusDisk, Flower *flower);

/*
 * Registers the flower should be removed from the disk.
 */
void cactusDisk_deleteFlowerFromDisk(CactusDisk *cactusDisk, Flower *flower);

/*
 * Functions on meta sequences.
 */

/*
 * Adds a newly constructed meta sequence to the memory of the cactusDisk.
 */
void cactusDisk_addMetaSequence(CactusDisk *cactusDisk,
        MetaSequence *metaSequence);

/*
 * Registers the meta sequence is being freed from memory.
 */
void cactusDisk_removeMetaSequence(CactusDisk *cactusDisk,
        MetaSequence *metaSequence);

/*
 * Functions on strings stored by the cactus disk.
 */

/*
 * Adds the sequence string to the bucket of sequence.
 *
 * This function is NOT thread safe.
 */
int64_t cactusDisk_addString(CactusDisk *cactusDisk, const char *string,
        int32_t length);

/*
 * Retrieves a string from the bucket of sequence.
 */
char *cactusDisk_getString(CactusDisk *cactusDisk, int64_t offset,
        int32_t start, int32_t length, int32_t strand);

#endif
