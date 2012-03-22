/*
 * Copyright (C) 2009-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cactus.h"
#include "sonLib.h"
#include "pairwiseAlignment.h"

/*
 * Iterates through a cigar file and converts the coordinates into cactus okay coordinates.
 */

static stHash *makeSequenceHeaderToSequenceHash(CactusDisk *cactusDisk) {
    stHash *sequenceHeaderToSequenceHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    Flower *flower = cactusDisk_getFlower(cactusDisk, 0);
    assert(flower != NULL);
    Sequence *sequence;
    Flower_SequenceIterator *seqIt = flower_getSequenceIterator(flower);
    while ((sequence = flower_getNextSequence(seqIt)) != NULL) {
        stHash_insert(sequenceHeaderToSequenceHash, cactusMisc_nameToString(sequence_getName(sequence)), sequence);
    }
    flower_destructSequenceIterator(seqIt);
    return sequenceHeaderToSequenceHash;
}

static void convertCoordinates(struct PairwiseAlignment *pairwiseAlignment, FILE *outputCigarFileHandle,
        stHash *sequenceHeaderToSequenceHash) {
    Sequence *sequence1 = stHash_search(sequenceHeaderToSequenceHash, pairwiseAlignment->contig1);
    Sequence *sequence2 = stHash_search(sequenceHeaderToSequenceHash, pairwiseAlignment->contig2);
    if (sequence1 == NULL) {
        st_errAbort("Could not match contig name in alignment to cactus sequence: %s", pairwiseAlignment->contig1);
    }
    if (sequence2 == NULL) {
        st_errAbort("Could not match contig name in alignment to cactus sequence: %s", pairwiseAlignment->contig2);
    }
    //Fix the names
    free(pairwiseAlignment->contig1);
    pairwiseAlignment->contig1 = cactusMisc_nameToString(sequence_getName(sequence1));
    free(pairwiseAlignment->contig2);
    pairwiseAlignment->contig2 = cactusMisc_nameToString(sequence_getName(sequence2));
    //Now fix the coordinates by adding one
    pairwiseAlignment->start1 += 2;
    pairwiseAlignment->start2 += 2;
    pairwiseAlignment->end1 += 2;
    pairwiseAlignment->end2 += 2;
    if (pairwiseAlignment->start1 < sequence_getStart(sequence1) || pairwiseAlignment->end1 > sequence_getStart(sequence1) + sequence_getLength(sequence1)) {
        st_errAbort("Coordinates of pairwise alignment appear incorrect: %i %i %i", pairwiseAlignment->start1,
                pairwiseAlignment->end1, sequence_getLength(sequence1));
    }
    if (pairwiseAlignment->start2 < sequence_getStart(sequence1) || pairwiseAlignment->end2 > sequence_getStart(sequence1) + sequence_getLength(sequence2)) {
        st_errAbort("Coordinates of pairwise alignment appear incorrect: %i %i %i", pairwiseAlignment->start2,
                pairwiseAlignment->end2, sequence_getLength(sequence2));
    }
}

int main(int argc, char *argv[]) {
    /*
     * Gets the cigars, iterates through them converting them
     */
    assert(argc == 5);
    st_setLogLevelFromString(argv[1]);
    st_logDebug("Set up logging\n");

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[2]);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    stHash *sequenceHeaderToSequenceHash = makeSequenceHeaderToSequenceHash(cactusDisk);
    st_logDebug("Set up the flower disk\n");

    FILE *inputCigarFileHandle = fopen(argv[3], "r");
    FILE *outputCigarFileHandle = fopen(argv[4], "w");
    struct PairwiseAlignment *pairwiseAlignment;
    while ((pairwiseAlignment = cigarRead(inputCigarFileHandle)) != NULL) {
        convertCoordinates(pairwiseAlignment, outputCigarFileHandle, sequenceHeaderToSequenceHash);
        destructPairwiseAlignment(pairwiseAlignment);
    }
    fclose(inputCigarFileHandle);
    fclose(outputCigarFileHandle);

    return 0;
    //Cleanup
    cactusDisk_destruct(cactusDisk);
    stHash_destruct(sequenceHeaderToSequenceHash);

    st_logDebug("Am finished\n");
    return 0;
}
