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
 * Functions for strings
 *
 * Sequence is held two bases to the byte rather than one. Every sequence reaching this point is
 * drawn from {A,C,G,T,N} in either case: cactus_sanitizeFastaHeaders folds the IUPAC ambiguity
 * codes to N and rejects anything else, and the ancestral base caller in blockMLString.c emits the
 * same alphabet. Ten symbols fit in a nibble with room to spare, so a genome costs half of what it
 * used to, for the whole life of the process.
 */

// Code 0 means "not a base we can store", so an unexpected character is caught rather than
// silently read back as an A.
static const uint8_t cactusDisk_baseToCode[256] = {
    ['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4, ['N'] = 5,
    ['a'] = 6, ['c'] = 7, ['g'] = 8, ['t'] = 9, ['n'] = 10
};

static const char cactusDisk_codeToBase[16] = {
    0, 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', 0, 0, 0, 0, 0
};

typedef struct _packedString {
    int64_t length;
    uint8_t bases[]; // Base i is in the low nibble of bases[i/2] if i is even, the high nibble if odd
} PackedString;

Name cactusDisk_addString(CactusDisk *cactusDisk, const char *string) {
    /*
     * Adds a string to the database.
     */
    Name name = cactusDisk_getUniqueID(cactusDisk);
    int64_t length = strlen(string);
    PackedString *packedString = st_malloc(sizeof(PackedString) + (length + 1) / 2);
    packedString->length = length;
    for (int64_t i = 0; i < length; i++) {
        uint8_t code = cactusDisk_baseToCode[(unsigned char) string[i]];
        if (code == 0) {
            st_errAbort("Sequence contains '%c' (0x%02x) at position %" PRIi64 ", but only A, C, G, T "
                        "and N, in either case, can be stored", string[i], (unsigned char) string[i], i);
        }
        if (i % 2 == 0) {
            packedString->bases[i / 2] = code; // Also zeroes the high nibble, so a trailing odd base
                                               // leaves the byte in a defined state
        } else {
            packedString->bases[i / 2] |= code << 4;
        }
    }
#if defined(_OPENMP)
    omp_set_lock(&(cactusDisk->writelock));
#endif
    stHash_insert(cactusDisk->allStrings, (void *)name, packedString); // Cheeky 64bit to pointer conversion
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
    PackedString *packedString = stHash_search(cactusDisk->allStrings, (void *)name); // Cheeky 64bit int to pointer conversion
#if defined(_OPENMP)
    omp_unset_lock(&(cactusDisk->writelock));
#endif

    assert(packedString != NULL);
    assert(start >= 0);
    assert(start + length <= packedString->length);
    char *string = st_malloc(length + 1);
    for (int64_t i = 0; i < length; i++) {
        int64_t j = start + i;
        uint8_t byte = packedString->bases[j / 2];
        string[i] = cactusDisk_codeToBase[j % 2 == 0 ? (byte & 0x0f) : (byte >> 4)];
    }
    string[length] = '\0';
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

CactusDisk *cactusDisk_construct() {
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
        flower_destruct(flower, FALSE, FALSE);
    }
    stSortedSet_destruct(cactusDisk->flowers);

    Sequence *sequence;
    while ((sequence = stSortedSet_getFirst(cactusDisk->sequences)) != NULL) {
        sequence_destruct(sequence);
    }
    stSortedSet_destruct(cactusDisk->sequences);
    stHash_destruct(cactusDisk->allStrings); // cleanup the library of strings we hold in memory

    if(cactusDisk->eventTree != NULL) {
        eventTree_destruct(cactusDisk->eventTree);
    }

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
