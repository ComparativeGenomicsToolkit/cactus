/*
 * Copyright (C) 2009-2013 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <getopt.h>

#include "cactus.h"
#include "sonLib.h"
#include "endAligner.h"
#include "flowerAligner.h"
#include "commonC.h"
#include "stCaf.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"

void usage() {
    fprintf(stderr, "cactus_baseReAligner [options] seq1[fasta] seq2[fasta], version 0.2\n");
    fprintf(stderr,
            "Realigns a set of pairwise alignments, as cigars, read from the command line and written back to the command line\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-l --gapGamma : (float [0, 1]) The gap gamma (as in the AMAP function)\n");
    fprintf(stderr,
            "-o --splitMatrixBiggerThanThis : (int >= 0)  No dp matrix bigger than this number squared will be computed.\n");
    fprintf(stderr, "-r --digaonalExpansion : (int >= 0 and even) Number of x-y diagonals to expand around anchors\n");
    fprintf(stderr, "-t --constraintDiagonalTrim : (int >= 0) Amount to trim from ends of each anchor\n");
    fprintf(stderr,
            "-w --alignAmbiguityCharacters : Align ambiguity characters (anything not ACTGactg) as a wildcard\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

struct PairwiseAlignment *convertAlignedPairsToPairwiseAlignment(char *seqName1, char *seqName2, int64_t score,
        stList *alignedPairs) {
    //Make pairwise alignment
    int64_t pX = -1, pY = -1, mL = 0;
    struct List *opList = constructEmptyList(0, (void(*)(void *)) destructAlignmentOperation);
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(alignedPair, 1);
        int64_t y = stIntTuple_get(alignedPair, 2);
        assert(x - pX > 0);
        assert(y - pY > 0);
        if (x - pX > 1) { //There is an indel.
            if (mL > 0) {
                listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, mL, 0));
                mL = 0;
            }
            listAppend(opList, constructAlignmentOperation(PAIRWISE_INDEL_X, x - pX -1, 0));
        } else if (y - pY > 1) {
            if (mL > 0) {
                listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, mL, 0));
                mL = 0;
            }
            listAppend(opList, constructAlignmentOperation(PAIRWISE_INDEL_Y, y - pY - 1, 0));
        }
        mL++;
        pX = x;
        pY = y;
    }
    if (mL > 0) {
        listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, mL, 0));
    }
    struct PairwiseAlignment *pA = constructPairwiseAlignment(seqName1, 0, pX+1, 1, seqName2, 0, pY+1, 1, score, opList);
    return pA;
}

void rebasePairwiseAlignmentCoordinates(int64_t *start, int64_t *end, int64_t *strand, int64_t coordinateShift,
        bool flipStrand) {
    *start += coordinateShift;
    *end += coordinateShift;
    if (flipStrand) {
        *strand = !(*strand);
    }
}

char *getSubSequence(char *seq, int64_t start, int64_t end, bool strand) {
    return strand ? stString_getSubString(seq, start, end - start) : cactusMisc_reverseComplementString(
            stString_getSubString(seq, end, start - end));
}

stHash *sequences = NULL;
void addToSequencesHash(const char *header, const char *sequence, int64_t length) {
    if (stHash_search(sequences, (char *) header) != NULL) {
        assert(stString_eq(sequence, stHash_search(sequences, header)));
    } else {
        stHash_insert(sequences, stString_copy(header), stString_copy(sequence));
    }
}

int main(int argc, char *argv[]) {
    char * logLevelString = NULL;
    float gapGamma = 0.5;
    int64_t i, j;
    PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters = pairwiseAlignmentBandingParameters_construct();

    /*
     * Parse the options.
     */
    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "help", no_argument, 0,
                'h' }, { "gapGamma", required_argument, 0, 'l' }, { "splitMatrixBiggerThanThis", required_argument, 0,
                'o' }, { "diagonalExpansion", required_argument, 0, 'r' }, { "constraintDiagonalTrim",
                required_argument, 0, 't' }, { "alignAmbiguityCharacters", no_argument, 0, 'w' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:hl:o:r:t:wx:y:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                st_setLogLevelFromString(logLevelString);
                break;
            case 'h':
                usage();
                return 0;
            case 'l':
                i = sscanf(optarg, "%f", &gapGamma);
                assert(i == 1);
                assert(gapGamma >= 0.0);
                break;
            case 'o':
                i = sscanf(optarg, "%" PRIi64 "", &j);
                assert(i == 1);
                assert(j >= 0);
                pairwiseAlignmentBandingParameters->splitMatrixBiggerThanThis = (int64_t) j * j;
                break;
            case 'r':
                i = sscanf(optarg, "%" PRIi64 "", &pairwiseAlignmentBandingParameters->diagonalExpansion);
                assert(i == 1);
                assert(pairwiseAlignmentBandingParameters->diagonalExpansion >= 0);
                assert(pairwiseAlignmentBandingParameters->diagonalExpansion % 2 == 0);
            case 't':
                i = sscanf(optarg, "%" PRIi64 "", &pairwiseAlignmentBandingParameters->constraintDiagonalTrim);
                assert(i == 1);
                assert(pairwiseAlignmentBandingParameters->constraintDiagonalTrim >= 0);
                break;
            case 'w':
                pairwiseAlignmentBandingParameters->alignAmbiguityCharacters = 1;
                break;
            default:
                usage();
                return 1;
        }
    }

    assert(optind + 2 == argc);
    char *seqFile1 = argv[optind++], *seqFile2 = argv[optind++];

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Starting realigning pairwise alignments\n");

    //Read in input sequences
    sequences = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    FILE *seqFileHandle = fopen(seqFile1, "r");
    fastaReadToFunction(seqFileHandle, addToSequencesHash);
    fclose(seqFileHandle);
    seqFileHandle = fopen(seqFile2, "r");
    fastaReadToFunction(seqFileHandle, addToSequencesHash);
    fclose(seqFileHandle);
    //Now do the business of processing the sequences.
    struct PairwiseAlignment *pA;
    FILE *fileHandleIn = stdin;
    FILE *fileHandleOut = stdout;
    while ((pA = cigarRead(fileHandleIn)) != NULL) {
        //Get sequences
        char *seqX = stHash_search(sequences, pA->contig1);
        char *seqY = stHash_search(sequences, pA->contig2);
        assert(seqX != NULL && seqY != NULL);
        //Convert to an alignment on the forward strand starting at 0
        bool flipStrand1 = !pA->strand1, flipStrand2 = !pA->strand2;
        int64_t coordinateShift1 = (pA->strand1 ? pA->start1 : pA->end1);
        int64_t coordinateShift2 = (pA->strand2 ? pA->start2 : pA->end2);
        char *subSeqX = getSubSequence(seqX, pA->start1, pA->end1, pA->strand1);
        char *subSeqY = getSubSequence(seqY, pA->start2, pA->end2, pA->strand2);
        rebasePairwiseAlignmentCoordinates(&(pA->start1), &(pA->end1), &(pA->strand1), -coordinateShift1, flipStrand1);
        rebasePairwiseAlignmentCoordinates(&(pA->start2), &(pA->end2), &(pA->strand2), -coordinateShift2, flipStrand2);
        //Convert input alignment into anchor pairs
        stList *anchorPairs = convertPairwiseForwardStrandAlignmentToAnchorPairs(pA,
                pairwiseAlignmentBandingParameters->constraintDiagonalTrim);
        //Get posterior prob pairs
        stList *alignedPairs = getAlignedPairsUsingAnchors(subSeqX, subSeqY, anchorPairs,
                pairwiseAlignmentBandingParameters, 1, 1);
        //Convert to partial ordered set of pairs
        stList_destruct(filterPairwiseAlignmentToMakePairsOrdered(alignedPairs, gapGamma));
        //Convert back to cigar
        struct PairwiseAlignment *rPA = convertAlignedPairsToPairwiseAlignment(pA->contig1, pA->contig2, pA->score,
                alignedPairs);
        //Rebase realigned-pA.
        rebasePairwiseAlignmentCoordinates(&(rPA->start1), &(rPA->end1), &(rPA->strand1), coordinateShift1, flipStrand1);
        rebasePairwiseAlignmentCoordinates(&(rPA->start2), &(rPA->end2), &(rPA->strand2), coordinateShift2, flipStrand2);
        //Write out alignment
        cigarWrite(fileHandleOut, rPA, 0);
        //Clean up
        stList_destruct(alignedPairs);
        destructPairwiseAlignment(pA);
        destructPairwiseAlignment(rPA);
        free(subSeqX);
        free(subSeqY);
    }
    stHash_destruct(sequences);

    st_logInfo("Finished realigning pairwise alignments, exiting.\n");

    //while(1);

    return 0;
}

