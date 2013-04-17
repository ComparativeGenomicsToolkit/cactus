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

static stHash *makeSequenceHeaderToCapHash(CactusDisk *cactusDisk) {
    stHash *sequenceHeaderToCapsHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    Flower *flower = cactusDisk_getFlower(cactusDisk, 0);
    assert(flower != NULL);
    Cap *cap;
    Flower_CapIterator *capIt = flower_getCapIterator(flower);
    while ((cap = flower_getNextCap(capIt)) != NULL) {
        cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
        assert(cap_getAdjacency(cap) != NULL);
        if (!cap_getSide(cap)) {
            Sequence *sequence = cap_getSequence(cap);
            assert(sequence != NULL);
            stList *sequenceHeaderTokens = stString_split(sequence_getHeader(sequence));
            if (stList_length(sequenceHeaderTokens) == 0) {
                st_errAbort("Sequence has absent header: %s", sequence_getHeader(sequence));
            }
            char *sequenceNameString = stString_copy(stList_get(sequenceHeaderTokens, 0));
            stList_destruct(sequenceHeaderTokens);
            if (stHash_search(sequenceHeaderToCapsHash, sequenceNameString) != NULL) {
                st_errAbort("Could not make a unique map of fasta headers to sequence names: '%s'", sequenceNameString);
            }
            stHash_insert(sequenceHeaderToCapsHash, sequenceNameString, cap);
        }
    }
    flower_destructCapIterator(capIt);
    return sequenceHeaderToCapsHash;
}

static void convertCoordinates(struct PairwiseAlignment *pairwiseAlignment, FILE *outputCigarFileHandle,
        stHash *sequenceHeaderToCapHash) {
    Cap *cap1 = stHash_search(sequenceHeaderToCapHash, pairwiseAlignment->contig1);
    Cap *cap2 = stHash_search(sequenceHeaderToCapHash, pairwiseAlignment->contig2);
    if (cap1 == NULL) {
        st_errAbort("Could not match contig name in alignment to cactus cap: '%s'", pairwiseAlignment->contig1);
    }
    if (cap2 == NULL) {
        st_errAbort("Could not match contig name in alignment to cactus cap: '%s'", pairwiseAlignment->contig2);
    }
    //Fix the names
    free(pairwiseAlignment->contig1);
    pairwiseAlignment->contig1 = cactusMisc_nameToString(cap_getName(cap1));
    free(pairwiseAlignment->contig2);
    pairwiseAlignment->contig2 = cactusMisc_nameToString(cap_getName(cap2));
    //Now fix the coordinates by adding one
    pairwiseAlignment->start1 += 2;
    pairwiseAlignment->start2 += 2;
    pairwiseAlignment->end1 += 2;
    pairwiseAlignment->end2 += 2;
    if (pairwiseAlignment->start1 <= cap_getCoordinate(cap1) || pairwiseAlignment->end1 > cap_getCoordinate(cap_getAdjacency(cap1))) {
        st_errAbort("Coordinates of pairwise alignment appear incorrect: %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 "", pairwiseAlignment->start1, pairwiseAlignment->end1,
                cap_getCoordinate(cap1), cap_getCoordinate(cap_getAdjacency(cap1)));
    }
    if (pairwiseAlignment->start2 <= cap_getCoordinate(cap2) || pairwiseAlignment->end2 > cap_getCoordinate(cap_getAdjacency(cap2))) {
            st_errAbort("Coordinates of pairwise alignment appear incorrect: %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 "", pairwiseAlignment->start2, pairwiseAlignment->end2,
                    cap_getCoordinate(cap2), cap_getCoordinate(cap_getAdjacency(cap2)));
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
    stHash *sequenceHeaderToCapHash = makeSequenceHeaderToCapHash(cactusDisk);
    st_logDebug("Set up the flower disk and built hash\n");

    FILE *inputCigarFileHandle = fopen(argv[3], "r");
    FILE *outputCigarFileHandle = fopen(argv[4], "w");
    st_logDebug("Opened files for writing\n");

    struct PairwiseAlignment *pairwiseAlignment;
    while ((pairwiseAlignment = cigarRead(inputCigarFileHandle)) != NULL) {
        convertCoordinates(pairwiseAlignment, outputCigarFileHandle, sequenceHeaderToCapHash);
        cigarWrite(outputCigarFileHandle, pairwiseAlignment, 0);
        destructPairwiseAlignment(pairwiseAlignment);
    }
    st_logDebug("Finished converting alignments\n");

    fclose(inputCigarFileHandle);
    fclose(outputCigarFileHandle);

    return 0;
    //Cleanup
    cactusDisk_destruct(cactusDisk);
    stHash_destruct(sequenceHeaderToCapHash);

    st_logDebug("Am finished\n");
    return 0;
}
