/*
 * Copyright (C) 2009-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus.h"
#include "sonLib.h"
#include "paf.h"
#include "bioioC.h"

void stripUniqueIdsFromLeafSequences(Flower *flower) {
    Flower_SequenceIterator *flowerIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(flowerIt)) != NULL) {
        // Strip the ID token from the header (should be the first
        // |-separated token) and complain if there isn't one.
        const char *header = sequence_getHeader(sequence);
        stList *tokens = fastaDecodeHeader(header);
        if(stList_length(tokens) > 1 && !strncmp(stList_get(tokens, 0), "id=", 3)) {
            free(stList_removeFirst(tokens));
            char *newHeader = fastaEncodeHeader(tokens);
            sequence_setHeader(sequence, newHeader);
        }
        stList_destruct(tokens);
    }
    flower_destructSequenceIterator(flowerIt);
}

/*
 * Iterates through a cigar file and converts the coordinates into cactus okay coordinates.
 */

static stHash *makeSequenceHeaderToCapHash(Flower *flower) {
    stHash *sequenceHeaderToCapsHash = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
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

static void convertCoordinates(Paf *paf, FILE *outputCigarFileHandle,
                               stHash *sequenceHeaderToCapHash) {
    Cap *cap1 = stHash_search(sequenceHeaderToCapHash, paf->query_name);
    Cap *cap2 = stHash_search(sequenceHeaderToCapHash, paf->target_name);
    if (cap1 == NULL) {
        st_errAbort("Could not match contig name in alignment to cactus cap: '%s'", paf->query_name);
    }
    if (cap2 == NULL) {
        st_errAbort("Could not match contig name in alignment to cactus cap: '%s'", paf->target_name);
    }
    //Fix the names
    free(paf->query_name);
    paf->query_name = cactusMisc_nameToString(cap_getName(cap1));
    free(paf->target_name);
    paf->target_name = cactusMisc_nameToString(cap_getName(cap2));
    //Now fix the coordinates by adding one
    paf->query_start += 2;
    paf->target_start += 2;
    paf->query_end += 2;
    paf->target_end += 2;
    paf->query_length += 2;
    paf->target_length += 2;
    paf_check(paf);
    if (paf->query_start <= cap_getCoordinate(cap1) || paf->query_end > cap_getCoordinate(cap_getAdjacency(cap1))) {
        st_errAbort("Coordinates of pairwise alignment appear incorrect: %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 "",
                    paf->query_start, paf->query_end,
                    cap_getCoordinate(cap1), cap_getCoordinate(cap_getAdjacency(cap1)));
    }
    if (paf->target_start <= cap_getCoordinate(cap2) || paf->target_end > cap_getCoordinate(cap_getAdjacency(cap2))) {
        st_errAbort("Coordinates of pairwise alignment appear incorrect: %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 "",
                    paf->target_start, paf->target_end,
                    cap_getCoordinate(cap2), cap_getCoordinate(cap_getAdjacency(cap2)));
    }
}

void convertAlignmentCoordinates(char *inputAlignmentFile, char *outputAlignmentFile, Flower *flower) {
    stHash *sequenceHeaderToCapHash = makeSequenceHeaderToCapHash(flower);
    st_logDebug("Set up the flower disk and built hash\n");

    FILE *inputAlignmentFileHandle = fopen(inputAlignmentFile, "r");
    FILE *outputAlignmentFileHandle = fopen(outputAlignmentFile, "w");
    st_logDebug("Opened files for writing\n");

    Paf *paf;
    while ((paf = paf_read(inputAlignmentFileHandle)) != NULL) {
        convertCoordinates(paf, outputAlignmentFileHandle, sequenceHeaderToCapHash);
        paf_check(paf);
        paf_write(paf, outputAlignmentFileHandle);
        paf_destruct(paf);
    }
    st_logDebug("Finished converting alignments\n");

    //Cleanup
    fclose(inputAlignmentFileHandle);
    fclose(outputAlignmentFileHandle);
    stHash_destruct(sequenceHeaderToCapHash);
}
