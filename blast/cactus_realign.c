/*
 * Copyright (C) 2009-2013 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>

#include "cactus.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "commonC.h"

void usage() {
    fprintf(stderr, "cactus_baseReAligner [options] seq1[fasta] seq2[fasta], version 0.2\n");
    fprintf(stderr, "Realigns a set of pairwise alignments, as cigars, read from the command line and written back to the command line\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-l --gapGamma : (float [0, 1]) The gap gamma (as in the AMAP function)\n");
    fprintf(stderr, "-o --splitMatrixBiggerThanThis : (int >= 0)  No dp matrix bigger than this number squared will be computed.\n");
    fprintf(stderr, "-r --diagonalExpansion : (int >= 0 and even) Number of x-y diagonals to expand around anchors\n");
    fprintf(stderr, "-t --constraintDiagonalTrim : (int >= 0) Amount to trim from ends of each anchor\n");
    fprintf(stderr, "-w --alignAmbiguityCharacters : Align ambiguity characters (anything not ACTGactg) as a wildcard\n");
    fprintf(stderr, "-x --dummy : Do everything but realign, used for testing.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

struct PairwiseAlignment *convertAlignedPairsToPairwiseAlignment(char *seqName1, char *seqName2, int64_t score, int64_t length1,
        int64_t length2, stList *alignedPairs) {
    //Make pairwise alignment
    int64_t pX = -1, pY = -1, mL = 0;
    //Create an end matched pair, which is used to ensure the alignment has the correct end indels.
    struct List *opList = constructEmptyList(0, (void(*)(void *)) destructAlignmentOperation);
    stList_append(alignedPairs, stIntTuple_construct2(length1, length2));
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(alignedPair, 0);
        int64_t y = stIntTuple_get(alignedPair, 1);
        assert(x - pX > 0);
        assert(y - pY > 0);
        if(x - pX > 0 && y - pY > 0) { //This is a hack for filtering
            if (x - pX > 1) { //There is an indel.
                if (mL > 0) {
                    listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, mL, 0));
                    mL = 0;
                }
                listAppend(opList, constructAlignmentOperation(PAIRWISE_INDEL_X, x - pX - 1, 0));
            }
            if (y - pY > 1) {
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
    }
    //Deal with a trailing match, but exclude the final match
    if (mL > 1) {
        listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, mL - 1, 0));
    }
    stIntTuple_destruct(stList_pop(alignedPairs));
    //Construct the alignment
    struct PairwiseAlignment *pA = constructPairwiseAlignment(seqName1, 0, length1, 1, seqName2, 0, length2, 1, score, opList);
    return pA;
}

void rebasePairwiseAlignmentCoordinates(int64_t *start, int64_t *end, int64_t *strand, int64_t coordinateShift, bool flipStrand) {
    *start += coordinateShift;
    *end += coordinateShift;
    if (flipStrand) {
        *strand = *strand ? 0 : 1;
        int64_t i = *end;
        *end = *start;
        *start = i;
    }
}

char *getSubSequence(char *seq, int64_t start, int64_t end, bool strand) {
    if (strand) {
        return stString_getSubString(seq, start, end - start);
    }
    seq = stString_getSubString(seq, end, start - end);
    char *rSeq = cactusMisc_reverseComplementString(seq);
    free(seq);
    return rSeq;
}

stHash *sequences = NULL;
void addToSequencesHash(const char *header, const char *sequence, int64_t length) {
    stList *tokens = stString_split(header);
    char *firstToken = stList_get(tokens, 0);
    if (stHash_search(sequences, (char *) firstToken) != NULL) {
        st_logInfo("Got a repeat header: %s with sequence length: %" PRIi64 " vs. the existing hashed sequence of length: %" PRIi64 ", complete header: %s\n", (char *) firstToken, length, strlen(stHash_search(sequences, (char *) firstToken)), header);
        if(length > strlen(stHash_search(sequences, (char *) firstToken))) { //The new sequence is a more complete version of the original sequence (can happen with overlapping fragments).
            st_logInfo("Replacing sequence\n");
#ifndef NDEBUG
            //Check old sequence is substring of new string
            char *cA = stString_getSubString(sequence, 0, strlen(stHash_search(sequences, (char *) firstToken)));
            st_uglyf("Differences %s\n", cA);
            st_uglyf("DIfferences %s\n", stHash_search(sequences, (char *) firstToken));
            assert(stString_eq(cA, stHash_search(sequences, (char *) firstToken)));
            free(cA);
#endif
            //Remove the old sequence and cleanup
            free(stHash_removeAndFreeKey(sequences, firstToken)); //The old first token is not cleaned up currently - a memory leak
            //Now insert the new sequence
            stHash_insert(sequences, stString_copy(firstToken), stString_copy(sequence));
        }
    } else {
        st_logInfo("Adding sequence for header: %s, with length %" PRIi64 ", complete header: %s\n", (char *) firstToken, strlen(sequence), header);
        stHash_insert(sequences, stString_copy(firstToken), stString_copy(sequence));
    }
    stList_destruct(tokens);
}

void *convertToAnchorPair(void *aPair, void *extraArg) {
    stIntTuple *i = stIntTuple_construct2(stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2));
    stIntTuple_destruct(aPair);
    return i;
}

bool matchFn(void *aPair, void *seqs) {
    char x = ((char **)seqs)[0][stIntTuple_get(aPair, 0)];
    char y = ((char **)seqs)[1][stIntTuple_get(aPair, 1)];
    return x == y && toupper(x) != 'N';
}

bool gapGammaFilter(void *aPair, void *gapGamma) {
    bool b = ((double)stIntTuple_get(aPair, 0)) / PAIR_ALIGNMENT_PROB_1 >= *((float *)gapGamma);
    if(!b) { //Cleanup.
        stIntTuple_destruct(aPair);
    }
    return b;
}

int main(int argc, char *argv[]) {
    char * logLevelString = NULL;
    float gapGamma = 0.9;
    int64_t i, j;
    PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters = pairwiseAlignmentBandingParameters_construct();
    pairwiseAlignmentBandingParameters->constraintDiagonalTrim = 0;
    pairwiseAlignmentBandingParameters->splitMatrixBiggerThanThis = 10;
    pairwiseAlignmentBandingParameters->diagonalExpansion = 4;
    bool dummy = 0;

    /*
     * Parse the options.
     */
    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "help", no_argument, 0, 'h' }, { "gapGamma",
                required_argument, 0, 'l' }, { "splitMatrixBiggerThanThis", required_argument, 0, 'o' }, { "diagonalExpansion",
                required_argument, 0, 'r' }, { "constraintDiagonalTrim", required_argument, 0, 't' }, { "alignAmbiguityCharacters",
                no_argument, 0, 'w' }, { "dummy", no_argument, 0, 'x' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:hl:o:r:t:wx", long_options, &option_index);

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
            case 'x':
                dummy = 1;
                break;
            default:
                usage();
                return 1;
        }
    }

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Starting realigning pairwise alignments\n");

    //Read in input sequences
    sequences = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    assert(optind < argc);
    while (optind < argc) {
        FILE *seqFileHandle = fopen(argv[optind++], "r");
        fastaReadToFunction(seqFileHandle, addToSequencesHash);
        fclose(seqFileHandle);
    }

    //Now do the business of processing the sequences.
    struct PairwiseAlignment *pA;
    FILE *fileHandleIn = stdin;
    FILE *fileHandleOut = stdout;
    while ((pA = cigarRead(fileHandleIn)) != NULL) {
        st_logInfo("Processing alignment for sequences: %s and %s\n", pA->contig1, pA->contig2);
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
        checkPairwiseAlignment(pA);
        //Convert input alignment into anchor pairs
        stList *anchorPairs = convertPairwiseForwardStrandAlignmentToAnchorPairs(pA,
                pairwiseAlignmentBandingParameters->constraintDiagonalTrim);
        //Filter anchorPairs to remove anchor pairs that include mismatches
        char *seqs[2] = { subSeqX, subSeqY };
        stList *filteredAnchoredPairs = stList_filter2(anchorPairs, matchFn, seqs);
        //Get posterior prob pairs
        stList *alignedPairs = NULL;
        if (dummy) {
            alignedPairs = convertPairwiseForwardStrandAlignmentToAnchorPairs(pA, 0);
        } else {
            alignedPairs = getAlignedPairsUsingAnchors(subSeqX, subSeqY, filteredAnchoredPairs, pairwiseAlignmentBandingParameters, 1, 1);
            //Convert to partial ordered set of pairs
            if(gapGamma < 0.55) { //Shouldn't be needed if we only take pairs with > 50% posterior prob
                stList_destruct(filterPairwiseAlignmentToMakePairsOrdered(alignedPairs, gapGamma));
            }
            else { //Just filter by gap gamma (the alterative is quite expensive)
                stList *l = stList_filter2(alignedPairs, gapGammaFilter, &gapGamma);
                stList_setDestructor(alignedPairs, NULL);
                stList_setDestructor(l, (void (*)(void *))stIntTuple_destruct);
                stList_destruct(alignedPairs);
                alignedPairs = l;
            }
            //Convert to ordered list of sequence coordinate pairs
            stList_mapReplace(alignedPairs, convertToAnchorPair, NULL);
            stList_sort(alignedPairs, (int(*)(const void *, const void *)) stIntTuple_cmpFn); //Ensure we have an monotonically increasing ordering
        }
        //Convert back to cigar
        struct PairwiseAlignment *rPA = convertAlignedPairsToPairwiseAlignment(pA->contig1, pA->contig2, pA->score, pA->end1, pA->end2,
                alignedPairs);
        //Rebase realigned-pA.
        rebasePairwiseAlignmentCoordinates(&(rPA->start1), &(rPA->end1), &(rPA->strand1), coordinateShift1, flipStrand1);
        rebasePairwiseAlignmentCoordinates(&(rPA->start2), &(rPA->end2), &(rPA->strand2), coordinateShift2, flipStrand2);
        checkPairwiseAlignment(rPA);
        //Write out alignment
        cigarWrite(fileHandleOut, rPA, 0);
        //Clean up
        stList_destruct(filteredAnchoredPairs);
        stList_destruct(anchorPairs);
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

