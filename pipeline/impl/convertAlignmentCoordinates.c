/*
 * Copyright (C) 2009-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus.h"
#include "sonLib.h"
#include "pairwiseAlignment.h"
#include "bioioC.h"

void stripUniqueIdsFromSequences(Flower *flower) {
    Flower_SequenceIterator *flowerIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(flowerIt)) != NULL) {
        const char *header;
        char *firstToken, *newHeader;
        // Strip the ID token from the header (should be the first
        // |-separated token) and complain if there isn't one.
        header = sequence_getHeader(sequence);
        stList *tokens = fastaDecodeHeader(header);
        assert(stList_length(tokens) > 1);
        firstToken = stList_removeFirst(tokens);
        assert(!strncmp(firstToken, "id=", 3));
        free(firstToken);
        newHeader = fastaEncodeHeader(tokens);
        sequence_setHeader(sequence, newHeader);
        stList_destruct(tokens);
    }
    flower_destructSequenceIterator(flowerIt);
}

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

void convertAlignmentCoordinates(char *inputAlignmentFile, char *outputAlignmentFile, CactusDisk *cactusDisk) {
    stHash *sequenceHeaderToCapHash = makeSequenceHeaderToCapHash(cactusDisk);
    st_logDebug("Set up the flower disk and built hash\n");

    FILE *inputCigarFileHandle = fopen(inputAlignmentFile, "r");
    FILE *outputCigarFileHandle = fopen(outputAlignmentFile, "w");
    st_logDebug("Opened files for writing\n");

    struct PairwiseAlignment *pairwiseAlignment;
    while ((pairwiseAlignment = cigarRead(inputCigarFileHandle)) != NULL) {
        convertCoordinates(pairwiseAlignment, outputCigarFileHandle, sequenceHeaderToCapHash);
        cigarWrite(outputCigarFileHandle, pairwiseAlignment, 0);
        destructPairwiseAlignment(pairwiseAlignment);
    }
    st_logDebug("Finished converting alignments\n");

    //Cleanup
    fclose(inputCigarFileHandle);
    fclose(outputCigarFileHandle);
    stHash_destruct(sequenceHeaderToCapHash);
}
