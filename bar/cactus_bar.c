/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#define _POSIX_C_SOURCE 200809L

#include <assert.h>
#include <getopt.h>
#include <sys/mman.h>
#include <stdio.h>

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

void usage() {
    fprintf(stderr, "cactus_baseAligner [flower-names, ordered by order they should be processed], version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --cactusDisk : The location of the flower disk directory\n");
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

    fprintf(stderr, "-P --partialOrderAlignmentWindow (int >= 0): Use partial order aligner instead of Pecan for multiple alignment subproblems, on blocks up to given length (0=disable POA).\n");

    fprintf(stderr, "-h --help : Print this help screen\n");
}

stPinch *getNextAlignedPairAlignment(stSortedSetIterator *it) {
    AlignedPair *alignedPair = stSortedSet_getNext(it);
    if (alignedPair == NULL) {
        return NULL;
    }
    static stPinch pinch;
    stPinch_fillOut(&pinch, alignedPair->subsequenceIdentifier, alignedPair->reverse->subsequenceIdentifier, alignedPair->position,
            alignedPair->reverse->position, 1, alignedPair->strand == alignedPair->reverse->strand);
    return &pinch;
}

static int64_t minimumIngroupDegree = 0, minimumOutgroupDegree = 0, minimumDegree = 0, minimumNumberOfSpecies = 0;
static Flower *flower;

bool blockFilterFn(stPinchBlock *pinchBlock) {
    return !stCaf_containsRequiredSpecies(pinchBlock, flower, minimumIngroupDegree, minimumOutgroupDegree, minimumDegree, minimumNumberOfSpecies);
}

int main(int argc, char *argv[]) {

    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
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
    // toggle from pecan to abpoa for multiple alignment, by setting to non-zero
    // Note that poa uses about N^2, so maximum value is generally in 10s of kb
    int64_t poaWindow = 0;

    PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters = pairwiseAlignmentBandingParameters_construct();

    /*
     * Setup the input parameters for cactus core.
     */
    bool pruneOutStubAlignments = 0;

    /*
     * Parse the options.
     */
    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "cactusDisk", required_argument, 0, 'b' }, { "help", no_argument, 0, 'h' }, { "spanningTrees", required_argument, 0, 'i' },
                { "maximumLength", required_argument, 0, 'j' }, { "useBanding", no_argument, 0, 'k' },
                { "gapGamma", required_argument, 0, 'l' }, { "matchGamma", required_argument, 0, 'L' },
                { "splitMatrixBiggerThanThis", required_argument, 0, 'o' }, { "anchorMatrixBiggerThanThis",
                        required_argument, 0, 'p' }, { "repeatMaskMatrixBiggerThanThis", required_argument, 0, 'q' }, {
                        "diagonalExpansion", required_argument, 0, 'r' }, { "constraintDiagonalTrim", required_argument, 0, 't' }, {
                        "minimumDegree", required_argument, 0, 'u' }, { "alignAmbiguityCharacters", no_argument, 0, 'w' }, {
                        "pruneOutStubAlignments", no_argument, 0, 'y' }, {
                        "minimumIngroupDegree", required_argument, 0, 'A' }, { "minimumOutgroupDegree", required_argument, 0, 'B' },
                { "precomputedAlignments", required_argument, 0, 'D' }, {
                        "endAlignmentsToPrecomputeOutputFile", required_argument, 0, 'E' }, { "useProgressiveMerging",
                        no_argument, 0, 'F' }, { "calculateWhichEndsToComputeSeparately", no_argument, 0, 'G' }, { "largeEndSize",
                        required_argument, 0, 'I' },
                        {"ingroupCoverageFile", required_argument, 0, 'J'},
                        {"minimumSizeToRescue", required_argument, 0, 'K'},
                        {"minimumCoverageToRescue", required_argument, 0, 'M'},
                        { "minimumNumberOfSpecies", required_argument, 0, 'N' },
                        {"partialOrderAlignmentWindow", required_argument, 0, 'P'},                        
                        { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:hi:j:kl:o:p:q:r:t:u:wy:A:B:D:E:FGI:J:K:L:M:N:P:", long_options, &option_index);

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
                    st_errAbort("Error parsing poaLength parameter");
                }
                break;
            default:
                usage();
                return 1;
        }
    }

    st_setLogLevelFromString(logLevelString);

    /*
     * Load the flowerdisk
     */
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, false, true); //We precache the sequences
    st_logInfo("Set up the flower disk\n");

    /*
     * Load the hmm
     */
    StateMachine *sM = stateMachine5_construct(fiveState);

    /*
     * For each flower.
     */
    if (calculateWhichEndsToComputeSeparately) {
        stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
        if (stList_length(flowers) != 1) {
            st_errAbort("We are breaking up a flower's end alignments for precomputation but we have %" PRIi64 " flowers.\n", stList_length(flowers));
        }
        stSortedSet *endsToAlignSeparately = getEndsToAlignSeparately(stList_get(flowers, 0), maximumLength, largeEndSize);
        assert(stSortedSet_size(endsToAlignSeparately) != 1);
        stSortedSetIterator *it = stSortedSet_getIterator(endsToAlignSeparately);
        End *end;
        while ((end = stSortedSet_getNext(it)) != NULL) {
            fprintf(stdout, "%s\t%" PRIi64 "\t%" PRIi64 "\n", cactusMisc_nameToStringStatic(end_getName(end)), end_getInstanceNumber(end), getTotalAdjacencyLength(end));
        }
        return 0; //avoid cleanup costs
        stSortedSet_destructIterator(it);
        stSortedSet_destruct(endsToAlignSeparately);
    } else if (endAlignmentsToPrecomputeOutputFile != NULL) {
        /*
         * In this case we will align a set of end and save the alignments in a file.
         */
        stList *names = flowerWriter_parseNames(stdin);
        Flower *flower = cactusDisk_getFlower(cactusDisk, *((Name *)stList_get(names, 0)));
        FILE *fileHandle = fopen(endAlignmentsToPrecomputeOutputFile, "w");
        if (fileHandle == NULL) {
            st_errnoAbort("Opening end alignment file %s failed", endAlignmentsToPrecomputeOutputFile);
        }
        for(int64_t i=1; i<stList_length(names); i++) {
            End *end = flower_getEnd(flower, *((Name *)stList_get(names, i)));
            if (end == NULL) {
                st_errAbort("The end %" PRIi64 " was not found in the flower\n", *((Name *)stList_get(names, i)));
            }
            stSortedSet *endAlignment = makeEndAlignment(sM, end, spanningTrees, maximumLength, useProgressiveMerging,
                                                         matchGamma, pairwiseAlignmentBandingParameters, poaWindow);
            writeEndAlignmentToDisk(end, endAlignment, fileHandle);
            stSortedSet_destruct(endAlignment);
        }
        fclose(fileHandle);
        return 0; //avoid cleanup costs
        stList_destruct(names);
        st_logInfo("Finished precomputing end alignments\n");
    } else {
        /*
         * Compute complete flower alignments, possibly loading some precomputed alignments.
         */
        bedRegion *bedRegions = NULL;
        size_t numBeds = 0;
        if (ingroupCoverageFilePath != NULL) {
            // Pre-load the mmap for the coverage file.
            FILE *coverageFile = fopen(ingroupCoverageFilePath, "rb");
            if (coverageFile == NULL) {
                st_errnoAbort("Opening coverage file %s failed",
                              ingroupCoverageFilePath);
            }
            fseek(coverageFile, 0, SEEK_END);
            int64_t coverageFileLen = ftell(coverageFile);
            assert(coverageFileLen >= 0);
            assert(coverageFileLen % sizeof(bedRegion) == 0);
            if (coverageFileLen == 0) {
                // mmap doesn't like length-0 mappings, for obvious
                // reasons. Pretend that the coverage file doesn't
                // exist in this case, since it contains no data.
                ingroupCoverageFilePath = NULL;
            } else {
                // Establish a memory mapping for the file.
                bedRegions = mmap(NULL, coverageFileLen, PROT_READ, MAP_SHARED,
                                  fileno(coverageFile), 0);
                if (bedRegions == MAP_FAILED) {
                    st_errnoAbort("Failure mapping coverage file");
                }

                numBeds = coverageFileLen / sizeof(bedRegion);
            }
            fclose(coverageFile);
        }

        stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
        if (listOfEndAlignmentFiles != NULL && stList_length(flowers) != 1) {
            st_errAbort("We have precomputed alignments but %" PRIi64 " flowers to align.\n", stList_length(flowers));
        }
        cactusDisk_preCacheStrings(cactusDisk, flowers);
        for (j = 0; j < stList_length(flowers); j++) {
            flower = stList_get(flowers, j);
            st_logInfo("Processing a flower\n");

            stSortedSet *alignedPairs = makeFlowerAlignment3(sM, flower, listOfEndAlignmentFiles, spanningTrees, maximumLength,
                    useProgressiveMerging, matchGamma, pairwiseAlignmentBandingParameters, pruneOutStubAlignments, poaWindow);
            st_logInfo("Created the alignment: %" PRIi64 " pairs\n", stSortedSet_size(alignedPairs));
            stPinchIterator *pinchIterator = stPinchIterator_constructFromAlignedPairs(alignedPairs, getNextAlignedPairAlignment);

            /*
             * Run the cactus caf functions to build cactus.
             */
            stPinchThreadSet *threadSet = stCaf_setup(flower);
            stCaf_anneal(threadSet, pinchIterator, NULL);
            if (minimumDegree < 2) {
                stCaf_makeDegreeOneBlocks(threadSet);
            }
            if (minimumIngroupDegree > 0 || minimumOutgroupDegree > 0 || minimumDegree > 1) {
                stCaf_melt(flower, threadSet, blockFilterFn, 0, 0, 0, INT64_MAX);
            }

            if (ingroupCoverageFilePath != NULL) {
                // Rescue any sequence that is covered by outgroups
                // but currently unaligned into single-degree blocks.
                stPinchThreadSetIt pinchIt = stPinchThreadSet_getIt(threadSet);
                stPinchThread *thread;
                while ((thread = stPinchThreadSetIt_getNext(&pinchIt)) != NULL) {
                    Cap *cap = flower_getCap(flower,
                                             stPinchThread_getName(thread));
                    assert(cap != NULL);
                    Sequence *sequence = cap_getSequence(cap);
                    assert(sequence != NULL);
                    rescueCoveredRegions(thread, bedRegions, numBeds,
                                         sequence_getName(sequence),
                                         minimumSizeToRescue,
                                         minimumCoverageToRescue);
                }
                stCaf_joinTrivialBoundaries(threadSet);
            }

            stCaf_finish(flower, threadSet, chainLengthForBigFlower, longChain, INT64_MAX, INT64_MAX); //Flower now destroyed.
            stPinchThreadSet_destruct(threadSet);
            st_logInfo("Ran the cactus core script.\n");

            /*
             * Cleanup
             */
            //Clean up the sorted set after cleaning up the iterator
            stPinchIterator_destruct(pinchIterator);
            stSortedSet_destruct(alignedPairs);

            st_logInfo("Finished filling in the alignments for the flower\n");
        }
        stList_destruct(flowers);
        //st_errAbort("Done\n");
        /*
         * Write and close the cactusdisk.
         */
        cactusDisk_write(cactusDisk);
        return 0; //Exit without clean up is quicker, enable cleanup when doing memory leak detection.
        if (bedRegions != NULL) {
            // Clean up our mapping.
            munmap(bedRegions, numBeds * sizeof(bedRegion));
        }
    }


    ///////////////////////////////////////////////////////////////////////////
    // Cleanup
    ///////////////////////////////////////////////////////////////////////////

    stateMachine_destruct(sM);
    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
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

