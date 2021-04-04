/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_DISK_PRIVATE_H_
#define CACTUS_DISK_PRIVATE_H_

#include "cactusGlobals.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

struct _cactusDisk {
    stSortedSet *sequences;
    stSortedSet *flowers;
    EventTree *eventTree;
#if defined(_OPENMP)
    omp_lock_t writelock; // This lock used to gate access to concurrently accessed variables
#endif
    stHash *allStrings; // If the strings are being all stored in memory, a map of names to strings
    Name currentName; // Used as a counter for issuing names
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
 * Functions on meta sequences.
 */

/*
 * Adds a newly constructed meta sequence to the memory of the cactusDisk.
 */
void cactusDisk_addSequence(CactusDisk *cactusDisk,
        Sequence *sequence);

/*
 * Registers the meta sequence is being freed from memory.
 */
void cactusDisk_removeSequence(CactusDisk *cactusDisk,
        Sequence *sequence);

/*
 * Functions on strings stored by the cactus disk.
 */

/*
 * Adds the sequence string to the database.
 */
Name cactusDisk_addString(CactusDisk *cactusDisk, const char *string);

/*
 * Retrieves a string from the bucket of sequence.
 */
char *cactusDisk_getString(CactusDisk *cactusDisk, Name name,
        int64_t start, int64_t length, int64_t strand, int64_t totalSequenceLength);

/*
 * Set the event tree for this disk. (Hopefully this only happens once.)
 */
void cactusDisk_setEventTree(CactusDisk *cactusDisk, EventTree *eventTree);

#endif
