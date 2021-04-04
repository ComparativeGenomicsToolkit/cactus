/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

// OpenMP
#if defined(_OPENMP)
#include <omp.h>
#endif

/*
 * Functions on meta sequences.
 */

void cactusDisk_addSequence(CactusDisk *cactusDisk, Sequence *sequence) {
#if defined(_OPENMP)
    omp_set_lock(&(cactusDisk->writelock));
#endif
    assert(stSortedSet_search(cactusDisk->sequences, sequence) == NULL);
    stSortedSet_insert(cactusDisk->sequences, sequence);
#if defined(_OPENMP)
    omp_unset_lock(&(cactusDisk->writelock));
#endif
}

void cactusDisk_removeSequence(CactusDisk *cactusDisk, Sequence *sequence) {
#if defined(_OPENMP)
    omp_set_lock(&(cactusDisk->writelock));
#endif
    assert(stSortedSet_search(cactusDisk->sequences, sequence) != NULL);
    stSortedSet_remove(cactusDisk->sequences, sequence);
#if defined(_OPENMP)
    omp_unset_lock(&(cactusDisk->writelock));
#endif
}

/*
 * Functions used to precache the sequences in the database for a given set of flowers.
 */

typedef struct _substring {
    /*
     * Struct used to represent a substring of a string to be retrieved from the database.
     */
    Name name;
    int64_t start;
    int64_t length;
} Substring;

/*
 * Basic methods on substrings.
 */

static Substring *substring_construct(Name name, int64_t start, int64_t length) {
    Substring *substring = st_malloc(sizeof(Substring));
    substring->name = name;
    substring->start = start;
    substring->length = length;
    return substring;
}

static Substring *substring_clone(Substring *substring) {
    return substring_construct(substring->name, substring->start, substring->length);
}

static void substring_destruct(Substring *substring) {
    free(substring);
}

static int substring_cmp(Substring *substring1, Substring *substring2) {
    int i = cactusMisc_nameCompare(substring1->name, substring2->name);
    if (i != 0) {
        return i;
    }
    i = substring1->start < substring2->start ? -1 : (substring1->start > substring2->start ? 1 : 0);
    if (i != 0) {
        return i;
    }
    return substring1->length < substring2->length ? -1 : (substring1->length > substring2->length ? 1 : 0);
}

Name cactusDisk_addString(CactusDisk *cactusDisk, const char *string) {
    /*
     * Adds a string to the database.
     */
    Name name = cactusDisk_getUniqueID(cactusDisk);
#if defined(_OPENMP)
    omp_set_lock(&(cactusDisk->writelock));
#endif
    stHash_insert(cactusDisk->allStrings, (void *)name, stString_copy(string)); // Cheeky 64bit to pointer conversion
#if defined(_OPENMP)
    omp_unset_lock(&(cactusDisk->writelock));
#endif
    return name;
}

char *cactusDisk_getString(CactusDisk *cactusDisk, Name name, int64_t start, int64_t length, int64_t strand,
        int64_t totalSequenceLength) {
    /*
     * Gets a string from the database.
     *
     */
    assert(length >= 0);
    if (length == 0) {
        return stString_copy("");
    }

#if defined(_OPENMP)
    omp_set_lock(&(cactusDisk->writelock));
#endif
    char *string = stHash_search(cactusDisk->allStrings, (void *)name); // Cheeky 64bit int to pointer conversion
#if defined(_OPENMP)
    omp_unset_lock(&(cactusDisk->writelock));
#endif

    assert(string != NULL);
    string = stString_getSubString(string, start, length);
    if(!strand) {
        char *reverseComplement = stString_reverseComplementString(string);
        free(string);
        return reverseComplement;
    }
    return string;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic flower disk functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static int cactusDisk_constructFlowersP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(flower_getName((Flower *) o1), flower_getName((Flower *) o2));
}

static int cactusDisk_constructSequencesP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(sequence_getName((Sequence *) o1), sequence_getName((Sequence *) o2));
}

/*
 * The following two functions compress and decompress the data in the cactus disk..
 */

static CactusDisk *cactusDisk_construct() {
    CactusDisk *cactusDisk = st_calloc(1, sizeof(CactusDisk));
    cactusDisk->sequences = stSortedSet_construct3(cactusDisk_constructSequencesP, NULL);
    cactusDisk->flowers = stSortedSet_construct3(cactusDisk_constructFlowersP, NULL);
    cactusDisk->eventTree = NULL;
    cactusDisk->allStrings = stHash_construct2(NULL, free);
    cactusDisk->currentName = 1; // Start the naming of objects from 1
#if defined(_OPENMP)
        omp_init_lock(&(cactusDisk->writelock));
#endif
    return cactusDisk;
}

void cactusDisk_destruct(CactusDisk *cactusDisk) {
    Flower *flower;
    while ((flower = stSortedSet_getFirst(cactusDisk->flowers)) != NULL) {
        flower_destruct(flower, FALSE);
    }
    stSortedSet_destruct(cactusDisk->flowers);

    Sequence *sequence;
    while ((sequence = stSortedSet_getFirst(cactusDisk->sequences)) != NULL) {
        sequence_destruct(sequence);
    }
    stSortedSet_destruct(cactusDisk->sequences);
    stHash_destruct(cactusDisk->allStrings); // cleanup the library of strings we hold in memory

#if defined(_OPENMP)
    omp_destroy_lock(&(cactusDisk->writelock));
#endif

    free(cactusDisk);
}

Flower *cactusDisk_getFlower(CactusDisk *cactusDisk, Name flowerName) {
    Flower flower;
    flower.name = flowerName;
#if defined(_OPENMP)
        omp_set_lock(&(cactusDisk->writelock));
#endif
        Flower *flower2 = stSortedSet_search(cactusDisk->flowers, &flower);
#if defined(_OPENMP)
        omp_unset_lock(&(cactusDisk->writelock));
#endif
    return flower2;
}

Sequence *cactusDisk_getSequence(CactusDisk *cactusDisk, Name sequenceName) {
    Sequence sequence;
    sequence.name = sequenceName;
#if defined(_OPENMP)
    omp_set_lock(&(cactusDisk->writelock));
#endif
    Sequence *sequence2 = stSortedSet_search(cactusDisk->sequences, &sequence);
#if defined(_OPENMP)
    omp_unset_lock(&(cactusDisk->writelock));
#endif
    return sequence2;
}

/*
 * Private functions.
 */

void cactusDisk_addFlower(CactusDisk *cactusDisk, Flower *flower) {
#if defined(_OPENMP)
        omp_set_lock(&(cactusDisk->writelock));
#endif
        assert(stSortedSet_search(cactusDisk->flowers, flower) == NULL);
        stSortedSet_insert(cactusDisk->flowers, flower);
#if defined(_OPENMP)
        omp_unset_lock(&(cactusDisk->writelock));
#endif
}

void cactusDisk_removeFlower(CactusDisk *cactusDisk, Flower *flower) {
#if defined(_OPENMP)
        omp_set_lock(&(cactusDisk->writelock));
#endif
        assert(stSortedSet_search(cactusDisk->flowers, flower) != NULL);
        stSortedSet_remove(cactusDisk->flowers, flower);
#if defined(_OPENMP)
        omp_unset_lock(&(cactusDisk->writelock));
#endif
}

void cactusDisk_setEventTree(CactusDisk *cactusDisk, EventTree *eventTree) {
    cactusDisk->eventTree = eventTree;
}

/*
 * Function to get unique ID.
 */

int64_t cactusDisk_getUniqueIDInterval(CactusDisk *cactusDisk, int64_t intervalSize) {
#if defined(_OPENMP)
    omp_set_lock(&(cactusDisk->writelock));
#endif
    Name n = cactusDisk->currentName;
    cactusDisk->currentName += intervalSize;
#if defined(_OPENMP)
    omp_unset_lock(&(cactusDisk->writelock));
#endif
    return n;
}

int64_t cactusDisk_getUniqueID(CactusDisk *cactusDisk) {
    return cactusDisk_getUniqueIDInterval(cactusDisk, 1);
}

EventTree *cactusDisk_getEventTree(CactusDisk *cactusDisk) {
    return cactusDisk->eventTree;
}
