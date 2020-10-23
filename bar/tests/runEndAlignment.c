/*
 * This is a hook directly into the multiple alignment part of the end aligner.  It
 * can be run directly on fasta inputs, in order to develop and debug new ideas
 * without going through cactus_bar which can only read from the cactus_disk server.
 * Otherwise, it attempts to take in the same sorts of parameters as cactus_bar
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <getopt.h>
#include <sys/mman.h>
#include <stdio.h>
#include <stdio.h>
#include <time.h>

#include "cactus.h"
#include "sonLib.h"
#include "endAligner.h"
#include "flowerAligner.h"
#include "rescue.h"
#include "commonC.h"
#include "stCaf.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stateMachine.h"
#include "multipleAligner.h"
#include "poaAligner.h"

void usage() {
    fprintf(stderr, "cactus_runEndAlignment [input-fasta], version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(
            stderr,
            "-i --spanningTrees (int >= 0) : The number of spanning trees construct in forming the set of pairwise alignments to include. If the number of pairwise alignments is less than the product of the total number of sequences and the number of spanning trees then all pairwise alignments will be included.\n");
    fprintf(
            stderr,
            "-j --maximumLength (int  >= 0 ) : The maximum length of a sequence to align, only the prefix and suffix maximum length bases are aligned\n");
    fprintf(stderr, "-k --useBanding : Use banding to speed up the alignments\n");
    fprintf(stderr, "-l --gapGamma : (float >= 0) The gap gamma (as in the AMAP function)\n");
    fprintf(stderr, "-L --matchGamma : (float [0, 1]) The match gamma (the avg. weight or greater to be allowed in the alignment)\n");
    fprintf(stderr, "-o --splitMatrixBiggerThanThis : (int >= 0)  No dp matrix bigger than this number squared will be computed.\n");
    fprintf(stderr,
            "-p --anchorMatrixBiggerThanThis : (int >= 0)  Any matrix bigger than this number squared will be broken apart with banding.\n");
    fprintf(
            stderr,
            "-q --repeatMaskMatrixBiggerThanThis : (int >= 0) Any matrix bigger than this after initial banding will be broken apart without repeat masking of the sequences\n");
    fprintf(stderr, "-r --digaonalExpansion : (int >= 0 and even) Number of x-y diagonals to expand around anchors\n");
    fprintf(stderr, "-t --constraintDiagonalTrim : (int >= 0) Amount to trim from ends of each anchor\n");

    fprintf(stderr, "-u --minimumDegree : (int >= 0) Minimum number of sequences in a block to be included in the output graph\n");

    fprintf(stderr, "-w --alignAmbiguityCharacters : Align ambiguity characters (anything not ACTGactg) as a wildcard\n");

    fprintf(stderr, "-y --pruneOutStubAlignments : Prune out alignments of sequences that terminates in free stubs stubs\n");

    fprintf(stderr, "-A --minimumIngroupDegree : Number of ingroup sequences required in a block.\n");

    fprintf(stderr, "-B --minimumOutgroupDegree : Number of outgroup sequences required in a block.\n");

    fprintf(stderr, "-D --precomputedAlignments : Precomputed end alignments.\n");

    fprintf(stderr, "-E --endAlignmentsToPrecomputeOutputFile [fileName] : If this output file is provided then bar will read stdin first to parse the flower, then to parse the names of the end alignments to precompute. The results will be placed in this file.\n");

    fprintf(stderr,
            "-F --useProgressiveMerging : Use progressive merging instead of poset merging for constructing multiple sequence alignments.\n");

    fprintf(stderr, "-G --calculateWhichEndsToComputeSeparately : Decide which end alignments to compute separately.\n");

    fprintf(stderr, "-I --largeEndSize : The size of sequences in an end at which point to compute it separately.\n");

    fprintf(stderr, "-J --ingroupCoverageFile : Binary coverage file containing ingroup regions that are covered by outgroups. These regions will be 'rescued' into single-degree blocks if they haven't been aligned to anything after the bar phase finished.\n");

    fprintf(stderr, "-K --minimumSizeToRescue : Unaligned but covered segments must be at least this size to be rescued.\n");

    fprintf(stderr, "-M --minimumCoverageToRescue : Unaligned segments must have at least this proportion of their bases covered by an outgroup to be rescued.\n");

    fprintf(stderr, "-P --partialOrderAlignmentLength (int >= 0): Use partial order aligner instead of Pecan for multiple alignment subproblems, on blocks up to given length (0=disable POA).\n");

    fprintf(stderr, "-n --noPecan: Dont do pecan, just poa\n");
    
    fprintf(stderr, "-h --help : Print this help screen\n");
}

static int64_t minimumIngroupDegree = 0, minimumOutgroupDegree = 0, minimumDegree = 0, minimumNumberOfSpecies = 0;

/* read each fasta sequence into a SeqFrag record with dummy End names */
static void add_fasta_fragment(void* destination, const char *name, const char *seq, int64_t length) {
    stList* fragList = (stList*)destination;
    SeqFrag* frag = (SeqFrag*)malloc(sizeof(SeqFrag));
    frag->seq = stString_copy(seq);
    frag->length = length;
    frag->leftEndId = 0; //(int64_t)name;
    frag->rightEndId = 1; //(int64_t)name + 1;
    stList_append(fragList, frag);
    fprintf(stderr, "added fragment with name %s and length %ld\n", name, length);
}

static void print_results(MultipleAlignment *mA) {

    stListIterator *it = stList_getIterator(mA->alignedPairs);
    stIntTuple *mAP;
    while ((mAP = stList_getNext(it)) != NULL) {
        fprintf(stderr, "tuple %ld\t%ld\t%ld\t%ld\n", stIntTuple_get(mAP, 1), stIntTuple_get(mAP, 2), stIntTuple_get(mAP, 3), stIntTuple_get(mAP, 4)); 
    }
}

int main(int argc, char *argv[]) {

    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * in_fasta = NULL;
    int64_t i, j;
    int64_t spanningTrees = 10;
    int64_t maximumLength = 1500;
    bool useProgressiveMerging = 0;
    float matchGamma = 0.5;
    bool useBanding = 0;
    int64_t k;
    stList *listOfEndAlignmentFiles = NULL;
    char *endAlignmentsToPrecomputeOutputFile = NULL;
    bool calculateWhichEndsToComputeSeparately = 0;
    int64_t largeEndSize = 1000000;
    int64_t chainLengthForBigFlower = 1000000;
    int64_t longChain = 2;
    char *ingroupCoverageFilePath = NULL;
    int64_t minimumSizeToRescue = 1;
    double minimumCoverageToRescue = 0.0;
    int64_t poaWindow = 10000;
    bool doPecan = true;

    PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters = pairwiseAlignmentBandingParameters_construct();

    /*
     * Setup the input parameters for cactus core.
     */
    bool pruneOutStubAlignments = 0;

    /*
     * Parse the options.
     */
    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "cactusDisk", required_argument, 0, 'b' },
                                                { "help", no_argument, 0, 'h' }, { "spanningTrees", required_argument, 0, 'i' },
                                                { "maximumLength", required_argument, 0, 'j' }, { "useBanding", no_argument, 0, 'k' },
                                                { "gapGamma", required_argument, 0, 'l' }, { "matchGamma", required_argument, 0, 'L' },
                                                { "splitMatrixBiggerThanThis", required_argument, 0, 'o' },
                                                { "anchorMatrixBiggerThanThis", required_argument, 0, 'p' },
                                                { "repeatMaskMatrixBiggerThanThis", required_argument, 0, 'q' },
                                                { "diagonalExpansion", required_argument, 0, 'r' },
                                                { "constraintDiagonalTrim", required_argument, 0, 't' },
                                                { "minimumDegree", required_argument, 0, 'u' },
                                                { "alignAmbiguityCharacters", no_argument, 0, 'w' },
                                                { "pruneOutStubAlignments", no_argument, 0, 'y' },
                                                { "minimumIngroupDegree", required_argument, 0, 'A' },
                                                { "minimumOutgroupDegree", required_argument, 0, 'B' },
                                                { "precomputedAlignments", required_argument, 0, 'D' },
                                                { "endAlignmentsToPrecomputeOutputFile", required_argument, 0, 'E' },
                                                { "useProgressiveMerging",  no_argument, 0, 'F' },
                                                { "calculateWhichEndsToComputeSeparately", no_argument, 0, 'G' },
                                                { "largeEndSize", required_argument, 0, 'I' },
                                                {"ingroupCoverageFile", required_argument, 0, 'J'},
                                                {"minimumSizeToRescue", required_argument, 0, 'K'},
                                                {"minimumCoverageToRescue", required_argument, 0, 'M'},
                                                { "minimumNumberOfSpecies", required_argument, 0, 'N' },
                                                {"inputFasta", required_argument, 0, 'f'},
                                                {"partialOrderAlignmentLength", required_argument, 0, 'P'},
                                                {"noPecan", no_argument, 0, 'n'},
                                                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:hi:j:kl:o:p:q:r:t:u:wy:A:B:D:E:FGI:J:K:L:M:N:f:P:n", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                st_setLogLevelFromString(logLevelString);
                break;
            case 'b':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'f':
                in_fasta = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case 'i':
                i = sscanf(optarg, "%" PRIi64 "", &spanningTrees);
                (void) i;
                assert(i == 1);
                assert(spanningTrees >= 0);
                break;
            case 'j':
                i = sscanf(optarg, "%" PRIi64 "", &maximumLength);
                assert(i == 1);
                assert(maximumLength >= 0);
                break;
            case 'k':
                useBanding = !useBanding;
                break;
            case 'l':
                i = sscanf(optarg, "%f", &pairwiseAlignmentBandingParameters->gapGamma);
                assert(i == 1);
                assert(pairwiseAlignmentBandingParameters->gapGamma >= 0.0);
                break;
            case 'L':
                i = sscanf(optarg, "%f", &matchGamma);
                assert(i == 1);
                assert(matchGamma >= 0.0);
                break;
            case 'o':
                i = sscanf(optarg, "%" PRIi64 "", &k);
                assert(i == 1);
                assert(k >= 0);
                pairwiseAlignmentBandingParameters->splitMatrixBiggerThanThis = (int64_t) k * k;
                break;
            case 'p':
                i = sscanf(optarg, "%" PRIi64 "", &k);
                assert(i == 1);
                assert(k >= 0);
                pairwiseAlignmentBandingParameters->anchorMatrixBiggerThanThis = (int64_t) k * k;
                break;
            case 'q':
                i = sscanf(optarg, "%" PRIi64 "", &k);
                assert(i == 1);
                assert(k >= 0);
                pairwiseAlignmentBandingParameters->repeatMaskMatrixBiggerThanThis = (int64_t) k * k;
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
            case 'u':
                i = sscanf(optarg, "%" PRIi64 "", &minimumDegree);
                assert(i == 1);
                break;
            case 'w':
                pairwiseAlignmentBandingParameters->alignAmbiguityCharacters = 1;
                break;
            case 'y':
                pruneOutStubAlignments = 1;
                break;
            case 'A':
                i = sscanf(optarg, "%" PRIi64 "", &minimumIngroupDegree);
                assert(i == 1);
                break;
            case 'B':
                i = sscanf(optarg, "%" PRIi64 "", &minimumOutgroupDegree);
                assert(i == 1);
                break;
            case 'D':
                listOfEndAlignmentFiles = stString_split(optarg);
                break;
            case 'E':
                endAlignmentsToPrecomputeOutputFile = stString_copy(optarg);
                break;
            case 'F':
                useProgressiveMerging = 1;
                break;
            case 'G':
                calculateWhichEndsToComputeSeparately = 1;
                break;
            case 'I':
                i = sscanf(optarg, "%" PRIi64 "", &largeEndSize);
                assert(i == 1);
                break;
            case 'J':
                ingroupCoverageFilePath = stString_copy(optarg);
                break;
            case 'K':
                i = sscanf(optarg, "%" PRIi64, &minimumSizeToRescue);
                assert(i == 1);
                break;
            case 'M':
                i = sscanf(optarg, "%lf", &minimumCoverageToRescue);
                assert(i == 1);
                break;
            case 'N':
                i = sscanf(optarg, "%" PRIi64, &minimumNumberOfSpecies);
                if (i != 1) {
                    st_errAbort("Error parsing minimumNumberOfSpecies parameter");
                }
                break;
            case 'P':
                i = sscanf(optarg, "%" PRIi64 "", &poaWindow);
                if (i != 1) {
                    st_errAbort("Error parsing poaWindow parameter");
                }
                break;
            case 'n':
                doPecan = false;
                break;
            default:
                usage();
                return 1;
        }
    }

    if (in_fasta == NULL) {
        fprintf(stderr, "-f option is mandatory\n");
        exit(1);
    }
    st_setLogLevelFromString(logLevelString);

    StateMachine *sM = stateMachine5_construct(fiveState);

    fprintf(stderr, "opening %s\n", in_fasta);
    FILE* fastaFile = fopen(in_fasta, "r");
    stList* seqFrags = stList_construct3(0, free);

    fastaReadToFunction(fastaFile, seqFrags, add_fasta_fragment);

    // cap the lengths
    stListIterator *it = stList_getIterator(seqFrags);
    SeqFrag* frag;
    while ((frag = stList_getNext(it)) != NULL) {
        if (frag->length > maximumLength) {
            frag->length = maximumLength;
            frag->seq[maximumLength] = '\0';
        }
    }
    
    fclose(fastaFile);

    fprintf(stderr, "doing pecan alignment\n");
    clock_t t = clock();
    MultipleAlignment *mA = NULL;
    if (doPecan) {
        mA = makeAlignment(sM, seqFrags, spanningTrees, 100000000, useProgressiveMerging, pairwiseAlignmentBandingParameters->gapGamma, pairwiseAlignmentBandingParameters);
    }
    t = clock() - t;
    double pecanTime = ((double)t)/CLOCKS_PER_SEC;

    if (doPecan) {
        print_results(mA);
    }

    fprintf(stderr, "doing poa alignment\n");
    t = clock();
    MultipleAlignment *pA = makePartialOrderAlignment(sM, seqFrags, pairwiseAlignmentBandingParameters->gapGamma, pairwiseAlignmentBandingParameters, poaWindow);
    t = clock() - t;
    double abpoaTime = ((double)t)/CLOCKS_PER_SEC;

    print_results(pA);

    fprintf(stderr, "Pecan Time = %lf seconds   abPOA Time = %lf seconds\n", pecanTime, abpoaTime);

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup
    ///////////////////////////////////////////////////////////////////////////

    stList_destruct(seqFrags);

    stateMachine_destruct(sM);
    //destructCactusCoreInputParameters(cCIP);
    free(cactusDiskDatabaseString);
    if (listOfEndAlignmentFiles != NULL) {
        stList_destruct(listOfEndAlignmentFiles);
    }
    if (logLevelString != NULL) {
        free(logLevelString);
    }
    st_logInfo("Finished with the flower disk for this flower.\n");

    //while(1);
    return 0;
}
