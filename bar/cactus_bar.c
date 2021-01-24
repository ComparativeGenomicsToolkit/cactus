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
#include <time.h>

#include "cactus.h"
#include "sonLib.h"
#include "endAligner.h"
#include "poaBarAligner.h"
#include "flowerAligner.h"
#include "rescue.h"
#include "commonC.h"
#include "stCaf.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stateMachine.h"
#include "../submodules/cPecan/inc/multipleAligner.h"

void usage() {
    fprintf(stderr, "cactus_baseAligner [flower-names, ordered by order they should be processed], version 0.2\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-p --params : [Required] The cactus config file\n");
    fprintf(stderr, "-d --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-D --precomputedAlignments : Precomputed end alignments.\n");
    fprintf(stderr, "-E --endAlignmentsToPrecomputeOutputFile [fileName] : If this output file is provided then bar will read stdin first to parse the flower, then to parse the names of the end alignments to precompute. The results will be placed in this file.\n");
    fprintf(stderr, "-G --calculateWhichEndsToComputeSeparately : Decide which end alignments to compute separately.\n");
    //fprintf(stderr, "-J --ingroupCoverageFile : Binary coverage file containing ingroup regions that are covered by outgroups. These regions will be 'rescued' into single-degree blocks if they haven't been aligned to anything after the bar phase finished.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    stList *listOfEndAlignmentFiles = NULL;
    char *endAlignmentsToPrecomputeOutputFile = NULL;
    bool calculateWhichEndsToComputeSeparately = 0;
    //char *ingroupCoverageFilePath = NULL;
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * paramsFile = NULL;

    /*
     * Parse the options.
     */

    if(argc <= 1) {
        usage();
        return 0;
    }

    while (1) {
        static struct option long_options[] = {
                { "logLevel", required_argument, 0, 'l' },
                { "params", required_argument, 0, 'p' },
                { "cactusDisk", required_argument, 0, 'd' },
                { "help", no_argument, 0, 'h' },
                { "precomputedAlignments", required_argument, 0, 'D' },
                {"endAlignmentsToPrecomputeOutputFile", required_argument, 0, 'E' },
                { "calculateWhichEndsToComputeSeparately", no_argument, 0, 'G' },
                //{"ingroupCoverageFile", required_argument, 0, 'J'},
                        { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:hD:E:G:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = stString_copy(optarg);
                st_setLogLevelFromString(logLevelString);
                break;
            case 'p':
                paramsFile = optarg;
                break;
            case 'd':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;

            //case 'j':
            //    i = sscanf(optarg, "%" PRIi64 "", &maximumLength);
            //    assert(i == 1);
            //    assert(maximumLength >= 0);
            //    break;

            case 'D':
                listOfEndAlignmentFiles = stString_split(optarg);
                break;
            case 'E':
                endAlignmentsToPrecomputeOutputFile = stString_copy(optarg);
                break;
            case 'G':
                calculateWhichEndsToComputeSeparately = 1;
                break;
            //case 'J':
            //    ingroupCoverageFilePath = stString_copy(optarg);
            //    break;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Parse stuff
    //////////////////////////////////////////////

    // Load the params file
    CactusParams *params = cactusParams_load(paramsFile);
    st_logInfo("Loaded the parameters files, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    // Load the cactus disk
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, false, true); //We precache the sequences
    st_logInfo("Set up the flower disk\n");

    int64_t maximumLength = 1500;
    int64_t poaWindow = cactusParams_get_int(params, 2, "bar", "partialOrderAlignmentWindow");

    //////////////////////////////////////////////
    //Now the actual algorithms..
    //////////////////////////////////////////////

    /*
     * For each flower.
     */
    if (calculateWhichEndsToComputeSeparately) {
        if(poaWindow != 0) {
            return 0; // Do not compute ends separately if using the poa aligner, as the poa aligner is so fast
            // this is unnecessary
            // todo: avoid calling with this flag if using poaMode
        }
        stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
        if (stList_length(flowers) != 1) {
            st_errAbort("We are breaking up a flower's end alignments for precomputation but we have %" PRIi64 " flowers.\n", stList_length(flowers));
        }
        int64_t largeEndSize = cactusParams_get_int(params, 2, "caf", "largeEndSize");
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
            assert(poaWindow == 0);
            /*
             * Get parameters for alignment
             */
            StateMachine *sM = stateMachine5_construct(fiveState);
            PairwiseAlignmentParameters *p = pairwiseAlignmentParameters_constructFromCactusParams(params);
            int64_t spanningTrees = cactusParams_get_int(params, 2, "bar", "spanningTrees");
            bool useProgressiveMerging = cactusParams_get_int(params, 2, "bar", "useProgressiveMerging");
            float matchGamma = cactusParams_get_float(params, 2, "bar", "matchGamma");

            stSortedSet *endAlignment = makeEndAlignment(sM, end, spanningTrees, maximumLength, useProgressiveMerging,
                                                         matchGamma, p);
            stateMachine_destruct(sM);
            writeEndAlignmentToDisk(end, endAlignment, fileHandle);
            stSortedSet_destruct(endAlignment);
            pairwiseAlignmentBandingParameters_destruct(p);
        }
        fclose(fileHandle);
        return 0; //avoid cleanup costs
        stList_destruct(names);
        st_logInfo("Finished precomputing end alignments\n");
    } else {
        /*
         * Run the actual bar algorithm
         */
        stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
        bar(flowers, params, cactusDisk, listOfEndAlignmentFiles, 1);
        cactusDisk_write(cactusDisk); // Write and close the cactusdisk.
    }

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    free(cactusDiskDatabaseString);
    if (listOfEndAlignmentFiles != NULL) {
        stList_destruct(listOfEndAlignmentFiles);
    }
    if (logLevelString != NULL) {
        free(logLevelString);
    }
    st_logInfo("Finished with the flower disk for bar.\n");

    //while(1);

    return 0;
}

