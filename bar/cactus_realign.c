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
    fprintf(stderr, "cactus_relign [options] seq1[fasta] seq2[fasta], version 0.2\n");
    fprintf(stderr,
            "Realigns a set of pairwise alignments, as cigars, read from the command line and written back to the command line\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-l --gapGamma : (float [0, 1]) The gap gamma (as in the AMAP function)\n");
    fprintf(stderr,
            "-o --splitMatrixBiggerThanThis : (int >= 0)  No dp matrix bigger than this number squared will be computed.\n");
    fprintf(stderr, "-r --diagonalExpansion : (int >= 0 and even) Number of x-y diagonals to expand around anchors\n");
    fprintf(stderr, "-t --constraintDiagonalTrim : (int >= 0) Amount to trim from ends of each anchor\n");
    fprintf(stderr,
            "-w --alignAmbiguityCharacters : Align ambiguity characters (anything not ACTGactg) as a wildcard\n");
    fprintf(stderr,
            "-x --rescoreOriginalAlignment : Rescore the original alignment. The output cigar is the same alignment.\n");
    fprintf(stderr, "-i --rescoreByIdentity : Set score equal to alignment identity, treating indels as mismatches.\n");
    fprintf(stderr,
            "-j --rescoreByPosteriorProb : Set score equal to avg. posterior match probability, treating indels as residues with 0 match probability.\n");
    fprintf(stderr, "-k --rescoreByIdentityIgnoringGaps : Set score equal to alignment identity, ignoring indels.\n");
    fprintf(stderr,
            "-m --rescoreByPosteriorProbIgnoringGaps : Set score equal to avg. posterior match probability, ignoring gaps.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
    fprintf(stderr,
            "-s --splitIndelsLongerThanThis : Split alignments with consecutive runs of indels that are longer than this.\n");
    fprintf(stderr,
            "-u --outputPosteriorProbs [FILE] : Outputs the posterior match probs of positions in the alignment to the given tab separated file, each line being X-coordinate, Y-coordinate, posterior-match prob.\n");
    fprintf(stderr,
                "-v --outputExpectations [FILE] : Instead of realigning, switches to calculating expectations, dumping out expectations as matrix in the given file.\n");
    fprintf(stderr,
                "-y --loadHmm [FILE] : Loads HMM from given file.\n");
}

struct PairwiseAlignment *convertAlignedPairsToPairwiseAlignment(char *seqName1, char *seqName2, double score,
        int64_t length1, int64_t length2, stList *alignedPairs) {
    //Make pairwise alignment
    int64_t pX = -1, pY = -1, mL = 0;
    //Create an end matched pair, which is used to ensure the alignment has the correct end indels.
    struct List *opList = constructEmptyList(0, (void (*)(void *)) destructAlignmentOperation);
    stList_append(alignedPairs, stIntTuple_construct2(length1, length2));
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(alignedPair, 0);
        int64_t y = stIntTuple_get(alignedPair, 1);
        assert(x - pX > 0);
        assert(y - pY > 0);
        if (x - pX > 0 && y - pY > 0) { //This is a hack for filtering
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
    struct PairwiseAlignment *pA = constructPairwiseAlignment(seqName1, 0, length1, 1, seqName2, 0, length2, 1, score,
            opList);
    return pA;
}

// Check if a pairwise alignment has an indel longer than
// maxIndelLength (to avoid creating/destructing a list if
// unnecessary)
bool hasLongIndel(struct PairwiseAlignment *pA, int64_t maxIndelLength) {
    int64_t i;
    int64_t curIndelRunLength = 0;
    for (i = 0; i < pA->operationList->length; i++) {
        struct AlignmentOperation *op = pA->operationList->list[i];
        if (op->opType == PAIRWISE_MATCH) {
            curIndelRunLength = 0;
        } else {
            curIndelRunLength += op->length;
            if (curIndelRunLength > maxIndelLength) {
                return TRUE;
            }
        }
    }
    return FALSE;
}

// Split a pairwise alignment into two or more pairwise alignments if
// it has a long indel.
stList *splitPairwiseAlignment(const struct PairwiseAlignment *pA, const int64_t maxIndelLength) {
    stList *ret = stList_construct3(0, free);
    int64_t i = 0;
    int64_t j = 0;
    int64_t curPos1 = pA->start1;
    int64_t curPos2 = pA->start2;
    int64_t curIndelRunLength = 0;
    int64_t curStart1 = pA->start1;
    int64_t curStart2 = pA->start2;
    int64_t curEnd1 = 0;
    int64_t curEnd2 = 0;
    struct List *curOperationList = constructEmptyList(0, (void (*)(void *)) destructAlignmentOperation);
    // Temporary list of a run of indel operations so that we don't
    // end alignments with indels.
    struct List *indelOpList = constructEmptyList(0, (void (*)(void *)) destructAlignmentOperation);
    for (i = 0; i < pA->operationList->length; i++) {
        struct AlignmentOperation *op = pA->operationList->list[i];
        switch (op->opType) {
        case PAIRWISE_MATCH:
            if (curIndelRunLength > maxIndelLength && curOperationList->length != 0) {
                // The last indel run was too long, so discard those
                // last indels and append the alignment so far to the return
                // value.
                if (curOperationList->length != 0) {
                    stList_append(ret,
                            constructPairwiseAlignment(stString_copy(pA->contig1), curStart1, curEnd1, pA->strand1,
                                    stString_copy(pA->contig2), curStart2, curEnd2, pA->strand2, pA->score,
                                    curOperationList));
                }
                curOperationList = constructEmptyList(0, (void (*)(void *)) destructAlignmentOperation);
                indelOpList = constructEmptyList(0, (void (*)(void *)) destructAlignmentOperation);
                curStart1 = curPos1;
                curStart2 = curPos2;
                curEnd1 = curStart1;
                curEnd2 = curStart2;
            } else if (curOperationList->length == 0) {
                // Indel run at the start of the alignment
                indelOpList = constructEmptyList(0, (void (*)(void *)) destructAlignmentOperation);
                curStart1 = curPos1;
                curStart2 = curPos2;
                curEnd1 = curStart1;
                curEnd2 = curStart2;
            }
            curIndelRunLength = 0;
            // Since we're keeping the indel run, (or we've already
            // split and cleared the overly-long indel list), add the
            // indel run to the op list and clear the indel op buffer.
            for (j = 0; j < indelOpList->length; j++) {
                struct AlignmentOperation *indelOp = indelOpList->list[j];
                listAppend(curOperationList, indelOp);
            }
            indelOpList = constructEmptyList(0, (void (*)(void *)) destructAlignmentOperation);

            curPos1 += pA->strand1 ? op->length : -op->length;
            curPos2 += pA->strand2 ? op->length : -op->length;
            curEnd1 = curPos1;
            curEnd2 = curPos2;
            listAppend(curOperationList, constructAlignmentOperation(op->opType, op->length, op->score));
            break;
        case PAIRWISE_INDEL_X:
            curIndelRunLength += op->length;
            curPos1 += pA->strand1 ? op->length : -op->length;
            listAppend(indelOpList, constructAlignmentOperation(op->opType, op->length, op->score));
            break;
        case PAIRWISE_INDEL_Y:
            curIndelRunLength += op->length;
            curPos2 += pA->strand2 ? op->length : -op->length;
            listAppend(indelOpList, constructAlignmentOperation(op->opType, op->length, op->score));
            break;
        }
    }
    assert(curPos1 == pA->end1);
    assert(curPos2 == pA->end2);
    if (curOperationList->length != 0) {
        // Append the remaining pairwise alignment
        stList_append(ret,
                constructPairwiseAlignment(stString_copy(pA->contig1), curStart1, curEnd1, pA->strand1,
                        stString_copy(pA->contig2), curStart2, curEnd2, pA->strand2, pA->score, curOperationList));
    }
    // Check all alignments
    for (i = 0; i < stList_length(ret); i++) {
        checkPairwiseAlignment(stList_get(ret, i));
    }
    return ret;
}

void rebasePairwiseAlignmentCoordinates(int64_t *start, int64_t *end, int64_t *strand, int64_t coordinateShift,
bool flipStrand) {
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
        st_logInfo(
                "Got a repeat header: %s with sequence length: %" PRIi64 " vs. the existing hashed sequence of length: %" PRIi64 ", complete header: %s\n",
                (char *) firstToken, length, strlen(stHash_search(sequences, (char *) firstToken)), header);
        if (length > strlen(stHash_search(sequences, (char *) firstToken))) { //The new sequence is a more complete version of the original sequence (can happen with overlapping fragments).
            st_logInfo("Replacing sequence\n");
#ifndef NDEBUG
            //Check old sequence is substring of new string
            char *cA = stString_getSubString(sequence, 0, strlen(stHash_search(sequences, (char *) firstToken)));
            assert(stString_eq(cA, stHash_search(sequences, (char * ) firstToken)));
            free(cA);
#endif
            //Remove the old sequence and cleanup
            free(stHash_removeAndFreeKey(sequences, firstToken)); //The old first token is not cleaned up currently - a memory leak
            //Now insert the new sequence
            stHash_insert(sequences, stString_copy(firstToken), stString_copy(sequence));
        }
    } else {
        st_logInfo("Adding sequence for header: %s, with length %" PRIi64 ", complete header: %s\n",
                (char *) firstToken, strlen(sequence), header);
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
    char x = toupper(((char **) seqs)[0][stIntTuple_get(aPair, 0)]);
    char y = toupper(((char **) seqs)[1][stIntTuple_get(aPair, 1)]);
    return x == y && x != 'N';
}

bool gapGammaFilter(void *aPair, void *gapGamma) {
    bool b = ((double) stIntTuple_get(aPair, 0)) / PAIR_ALIGNMENT_PROB_1 >= *((float *) gapGamma);
    if (!b) { //Cleanup.
        stIntTuple_destruct(aPair);
    }
    return b;
}

/*
 * Functions to rescore an alignment by identity / or some proxy to it.
 */

static int64_t getNumberOfMatchingAlignedPairs(char *subSeqX, char *subSeqY, stList *alignedPairs) {
    /*
     * Gives the average identity of matches in the alignment, treating indels as mismatches.
     */
    int64_t matches = 0;
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        int64_t x = stIntTuple_get(aPair, 1), y = stIntTuple_get(aPair, 2);
        matches += toupper(subSeqX[x]) == toupper(subSeqY[y]) && toupper(subSeqX[x]) != 'N';
    }
    return matches;
}

double scoreByIdentity(char *subSeqX, char *subSeqY, int64_t lX, int64_t lY, stList *alignedPairs) {
    /*
     * Gives the average identity of matches in the alignment, treating indels as mismatches.
     */
    int64_t matches = getNumberOfMatchingAlignedPairs(subSeqX, subSeqY, alignedPairs);
    return 100.0 * ((lX + lY) == 0 ? 0 : (2.0 * matches) / (lX + lY));
}

double scoreByIdentityIgnoringGaps(char *subSeqX, char *subSeqY, stList *alignedPairs) {
    /*
     * Gives the average identity of matches in the alignment, ignoring indels.
     */
    int64_t matches = getNumberOfMatchingAlignedPairs(subSeqX, subSeqY, alignedPairs);
    return 100.0 * matches / (double) stList_length(alignedPairs);
}

static double totalScore(stList *alignedPairs) {
    double score = 0.0;
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        score += stIntTuple_get(aPair, 0);
    }
    return score;
}

double scoreByPosteriorProbability(int64_t lX, int64_t lY, stList *alignedPairs) {
    /*
     * Gives the average posterior match probability per base of the two sequences, treating bases in indels as having 0 match probability.
     */
    return 100.0 * ((lX + lY) == 0 ? 0 : (2.0 * totalScore(alignedPairs)) / ((lX + lY) * PAIR_ALIGNMENT_PROB_1));
}

double scoreByPosteriorProbabilityIgnoringGaps(stList *alignedPairs) {
    /*
     * Gives the average posterior match probability per base of the two sequences, ignoring indels.
     */
    return 100.0 * totalScore(alignedPairs) / ((double) stList_length(alignedPairs) * PAIR_ALIGNMENT_PROB_1);
}

void writePosteriorProbs(char *posteriorProbsFile, stList *alignedPairs) {
    /*
     * Writes the posterior match probabibilities to a tab separated file, each line being X coordinate, Y coordinate, Match probability
     */
    FILE *fH = fopen(posteriorProbsFile, "w");
    for(int64_t i=0;i<stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        fprintf(fH, "%" PRIi64 "\t%" PRIi64 "\t%f\n", stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2), ((double)stIntTuple_get(aPair, 0))/PAIR_ALIGNMENT_PROB_1);
    }
    fclose(fH);
}

stList *scoreAnchorPairs(stList *anchorPairs, stList *alignedPairs) {
    /*
     * Selects the aligned pairs contained in anchor pairs.
     */
    stSortedSet *anchorPairsSet = stList_getSortedSet(anchorPairs, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    assert(stList_length(anchorPairs) == stSortedSet_size(anchorPairsSet));
    stList *scoredAnchorPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);

    for(int64_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        stIntTuple *j = stIntTuple_construct2(stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2));
        if(stSortedSet_search(anchorPairsSet, j) != NULL) {
            stList_append(scoredAnchorPairs, stIntTuple_construct3(stIntTuple_get(aPair, 0), stIntTuple_get(aPair, 1), stIntTuple_get(aPair, 2)));
            stSortedSet_remove(anchorPairsSet, j);
        }
        stIntTuple_destruct(j);
    }

    //The following should not really be needed, and may be masking a bug/numerical precision issues
    stSortedSetIterator *it = stSortedSet_getIterator(anchorPairsSet);
    stIntTuple *pair;
    while((pair = stSortedSet_getNext(it))) {
        stList_append(scoredAnchorPairs, stIntTuple_construct3(0, stIntTuple_get(pair, 0), stIntTuple_get(pair, 1)));
    }
    stSortedSet_destructIterator(it);

    stSortedSet_destruct(anchorPairsSet);
    assert(stList_length(anchorPairs) == stList_length(scoredAnchorPairs));

    return scoredAnchorPairs;
}

int main(int argc, char *argv[]) {
    char * logLevelString = NULL;
    float gapGamma = 0.9;
    int64_t i, j;
    PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters = pairwiseAlignmentBandingParameters_construct();
    pairwiseAlignmentBandingParameters->constraintDiagonalTrim = 0;
    pairwiseAlignmentBandingParameters->splitMatrixBiggerThanThis = 10;
    pairwiseAlignmentBandingParameters->diagonalExpansion = 4;
    bool rescoreOriginalAlignment = 0;
    bool rescoreByIdentity = 0;
    bool rescoreByPosteriorProbability = 0;
    bool rescoreByIdentityIgnoringGaps = 0;
    bool rescoreByPosteriorProbabilityIgnoringGaps = 0;
    // -1 here signals don't split indels at all
    int64_t splitIndelsLongerThanThis = -1;
    char *posteriorProbsFile = NULL;
    char *expectationsFile = NULL;
    char *hmmFile = NULL;
    Hmm *hmmExpectations = NULL;
    /*
     * Parse the options.
     */

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' },
                { "help", no_argument, 0, 'h' }, { "gapGamma", required_argument, 0, 'l' }, {
                        "splitMatrixBiggerThanThis", required_argument, 0, 'o' }, { "diagonalExpansion",
                required_argument, 0, 'r' }, { "constraintDiagonalTrim", required_argument, 0, 't' }, {
                        "alignAmbiguityCharacters", no_argument, 0, 'w' }, { "rescoreOriginalAlignment", no_argument, 0,
                        'x' }, { "rescoreByIdentity", no_argument, 0, 'i' }, { "rescoreByPosteriorProb", no_argument, 0,
                        'j' }, { "rescoreByPosteriorProbIgnoringGaps", no_argument, 0, 'm' }, {
                        "rescoreByIdentityIgnoringGaps", no_argument, 0, 'k' },
                        { "splitIndelsLongerThanThis",required_argument, 0, 's' },
                { "outputPosteriorProbs", required_argument, 0, 'u' },
                { "outputExpectations", required_argument, 0, 'v' },
                { "loadHmm", required_argument, 0, 'y' },
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:hl:o:r:t:s:wxijkmuv:y:", long_options, &option_index);

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
            break;
        case 't':
            i = sscanf(optarg, "%" PRIi64 "", &pairwiseAlignmentBandingParameters->constraintDiagonalTrim);
            assert(i == 1);
            assert(pairwiseAlignmentBandingParameters->constraintDiagonalTrim >= 0);
            break;
        case 'w':
            pairwiseAlignmentBandingParameters->alignAmbiguityCharacters = 1;
            break;
        case 'x':
            rescoreOriginalAlignment = 1;
            break;
        case 'i':
            rescoreByIdentity = 1;
            break;
        case 'j':
            rescoreByPosteriorProbability = 1;
            break;
        case 'k':
            rescoreByIdentityIgnoringGaps = 1;
            break;
        case 'm':
            rescoreByPosteriorProbabilityIgnoringGaps = 1;
            break;
        case 's':
            i = sscanf(optarg, "%" PRIi64, &splitIndelsLongerThanThis);
            assert(i == 1);
            assert(splitIndelsLongerThanThis >= 0);
            break;
        case 'u':
            posteriorProbsFile = stString_copy(optarg);
            break;
        case 'v':
            expectationsFile = stString_copy(optarg);
            hmmExpectations = hmm_constructEmpty(0.000000000001); //The tiny pseudo count prevents overflow
            break;
        case 'y':
            hmmFile = stString_copy(optarg);
            break;
        default:
            usage();
            return 1;
        }
    }

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Starting realigning pairwise alignments\n");

    //Load the model, if specified
    if(hmmFile != NULL) {
        st_logInfo("Loading the hmm from file %s\n", hmmFile);
        Hmm *hmm = hmm_loadFromFile(hmmFile);
        loadTheGlobalHmm(hmm);
        //hmm_normalise(hmm);
        hmm_destruct(hmm);
    }

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
        if(expectationsFile != NULL) {
            st_logInfo("Computing expectations\n");
            getExpectationsUsingAnchors(hmmExpectations, subSeqX, subSeqY, filteredAnchoredPairs,
                                pairwiseAlignmentBandingParameters, 1, 1);
        }
        else {
            //Get posterior prob pairs
            stList *alignedPairs = getAlignedPairsUsingAnchors(subSeqX, subSeqY, filteredAnchoredPairs,
                    pairwiseAlignmentBandingParameters, 1, 1);
            //Convert to partial ordered set of pairs
            if (rescoreOriginalAlignment) {
                stList *rescoredPairs = scoreAnchorPairs(anchorPairs, alignedPairs);
                stList_destruct(alignedPairs);
                alignedPairs = rescoredPairs;
            } else if (1) { //Shouldn't be needed if we only take pairs with > 50% posterior prob
                stList_destruct(filterPairwiseAlignmentToMakePairsOrdered(alignedPairs, gapGamma));
            } else { //Just filter by gap gamma (the alterative is quite expensive)
                stList *l = stList_filter2(alignedPairs, gapGammaFilter, &gapGamma);
                stList_setDestructor(alignedPairs, NULL);
                stList_setDestructor(l, (void (*)(void *)) stIntTuple_destruct);
                stList_destruct(alignedPairs);
                alignedPairs = l;
            }
            //Rescore
            if (rescoreByPosteriorProbability) {
                pA->score = scoreByPosteriorProbability(strlen(subSeqX), strlen(subSeqY), alignedPairs);
            } else if (rescoreByPosteriorProbabilityIgnoringGaps) {
                pA->score = scoreByPosteriorProbabilityIgnoringGaps(alignedPairs);
            } else if (rescoreByIdentity) {
                pA->score = scoreByIdentity(subSeqX, subSeqY, strlen(subSeqX), strlen(subSeqY), alignedPairs);
            } else if (rescoreByIdentityIgnoringGaps) {
                pA->score = scoreByIdentityIgnoringGaps(subSeqX, subSeqY, alignedPairs);
            }
            //Output the posterior match probs, if needed
            if(posteriorProbsFile != NULL) {
                writePosteriorProbs(posteriorProbsFile, alignedPairs);
            }
            //Convert to ordered list of sequence coordinate pairs
            stList_mapReplace(alignedPairs, convertToAnchorPair, NULL);
            stList_sort(alignedPairs, (int (*)(const void *, const void *)) stIntTuple_cmpFn); //Ensure we have an monotonically increasing ordering
            //Convert back to cigar
            struct PairwiseAlignment *rPA = convertAlignedPairsToPairwiseAlignment(pA->contig1, pA->contig2, pA->score,
                    pA->end1, pA->end2, alignedPairs);
            //Rebase realigned-pA.
            rebasePairwiseAlignmentCoordinates(&(rPA->start1), &(rPA->end1), &(rPA->strand1), coordinateShift1,
                    flipStrand1);
            rebasePairwiseAlignmentCoordinates(&(rPA->start2), &(rPA->end2), &(rPA->strand2), coordinateShift2,
                    flipStrand2);
            checkPairwiseAlignment(rPA);
            //Write out alignment
            if (splitIndelsLongerThanThis != -1) {
                // Write multiple split alignments
                stList *pAs = splitPairwiseAlignment(rPA, splitIndelsLongerThanThis);
                for (i = 0; i < stList_length(pAs); i++) {
                    cigarWrite(fileHandleOut, stList_get(pAs, i), 0);
                }
                stList_destruct(pAs);
            } else {
                // Write just one unsplit alignment
                cigarWrite(fileHandleOut, rPA, 0);
            }

            //Clean up
            stList_destruct(alignedPairs);
            destructPairwiseAlignment(rPA);
        }
        destructPairwiseAlignment(pA);
        stList_destruct(filteredAnchoredPairs);
        stList_destruct(anchorPairs);
        free(subSeqX);
        free(subSeqY);
    }
    stHash_destruct(sequences);

    if(expectationsFile != NULL) {
        st_logInfo("Writing out expectations to file %s\n", expectationsFile);
        FILE *fH = fopen(expectationsFile, "w");
        hmm_write(hmmExpectations, fH);
        hmm_destruct(hmmExpectations);
        fclose(fH);
    }

    st_logInfo("Finished realigning pairwise alignments, exiting.\n");

    //while(1);

    return 0;
}
