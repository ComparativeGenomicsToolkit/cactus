/*
 * Released under the MIT license, see LICENSE.txt
 */

#include <time.h>
#include <getopt.h>
#include "sonLib.h"
#include "cactus.h"
#include "cactus_setup.h"
#include "stCaf.h"
#include "poaBarAligner.h"
#include "cactusReference.h"
#include "addReferenceCoordinates.h"
#include "traverseFlowers.h"
#include "blockMLString.h"
#include "hal.h"
#include "convertAlignmentCoordinates.h"

// OpenMP
#if defined(_OPENMP)
#include <omp.h>
#endif

/*
 * TODOs:
 *
 * cleanup the python
 * cleanup input alignment format
 *
 */

void usage() {
    fprintf(stderr, "cactus_consolidated, version 0.2\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-p --params : [Required] The cactus config file\n");
    fprintf(stderr, "-f --outputFile : [Required] The file to write the combined cactus to hal output\n");
    fprintf(stderr, "-F --outputHalFastaFile : The file to write the sequences in to build the hal file.\n");
    fprintf(stderr, "-G --outputReferenceFile : The file to write the sequences of the reference in (used in the progressive recursion).\n");
    fprintf(stderr, "-s --sequences [Required] [eventName fastaFile/Directory]xN: The sequences\n");
    fprintf(stderr, "-a --alignments : [Required] The alignments file\n");
    fprintf(stderr, "-S --secondaryAlignments : The secondary alignments file\n");
    fprintf(stderr, "-c --constraintAlignments : The constraint alignments file\n");
    fprintf(stderr, "-g --speciesTree : [Required] The species tree, which will form the skeleton of the event tree\n");
    fprintf(stderr, "-o --outgroupEvents : Leaf events in the species tree identified as outgroups\n");
    fprintf(stderr, "-r --referenceEvent : [Required] The name of the reference event\n");
    fprintf(stderr, "-t --runChecks : Run cactus checks after each stage, used for debugging\n");
    fprintf(stderr, "-T --threads : (int > 0) Use up to this many threads [default: all available]\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

static char *convertAlignments(char *alignmentsFile, Flower *flower) {
    char *tempFile = getTempFile();
    convertAlignmentCoordinates(alignmentsFile, tempFile, flower);
    return tempFile;
}

static RecordHolder *getMergedRecordHolders(stHash *recordHolders, Flower *flower) {
    stList *children = stList_construct();
    getChildFlowers(flower, children);
    RecordHolder *rh = recordHolder_construct();
    for(int64_t i=0; i<stList_length(children); i++) {
        RecordHolder *rh2 = stHash_search(recordHolders, stList_get(children, i));
        assert(rh2 != NULL);
        recordHolder_transferAll(rh, rh2);
    }
    stList_destruct(children);
    return rh;
}

static void callBottomUp(Flower *flower, RecordHolder *rh, void *extraArg) {
    bottomUpNoDb(flower, rh, (Name)extraArg, 0, generateJukesCantorMatrix);
}

static void callHalFn(Flower *flower, RecordHolder *rh, void *extraArg) {
    makeHalFormatNoDb(flower, rh, (Name)extraArg, NULL);
}

static RecordHolder *doBottomUpTraversal(stList *flowerLayers,
                                         void (*bottomUpFn)(Flower *, RecordHolder *, void *), void *extraArgs) {
    // Bottom-up reference coordinates phase
    stHash *recordHolders = stHash_construct();
    for(int64_t i=stList_length(flowerLayers)-1; i>0 ; i--) {
        stList *flowers = stList_get(flowerLayers, i);

        // List to keep the RecordHolder for each flower
        stList *recordHoldersForFlowers = stList_construct3(stList_length(flowers), NULL);

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
        for (int64_t j = 0; j < stList_length(flowers); j++) {
            stList_set(recordHoldersForFlowers, j, getMergedRecordHolders(recordHolders, stList_get(flowers, j)));
            bottomUpFn(stList_get(flowers, j), stList_get(recordHoldersForFlowers, j), extraArgs);
        }

        // Make new map of flowers in the layer to RecordHolders
        stHash_destruct(recordHolders);
        recordHolders = stHash_construct();
        for (int64_t j = 0; j < stList_length(flowers); j++) {
            stHash_insert(recordHolders, stList_get(flowers, j), stList_get(recordHoldersForFlowers, j));
        }
        stList_destruct(recordHoldersForFlowers);
    }
    RecordHolder *rh = getMergedRecordHolders(recordHolders, stList_get(stList_get(flowerLayers, 0), 0));
    stHash_destruct(recordHolders);
    return rh;
}

int flower_sizeCmpFn(const void *a, const void *b) {
    // Sort by number of caps the flowers contains
    int64_t i = flower_getCapNumber((Flower *)a), j = flower_getCapNumber((Flower *)b);
    return i < j ? 1 : (i > j ? -1 : 0); // Sort in descending order
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *paramsFile = NULL;
    char *outputFile = NULL;
    char *outputHalFastaFile = NULL;
    char *outputReferenceFile = NULL;
    char *sequenceFilesAndEvents = NULL;
    char *alignmentsFile = NULL;
    char *secondaryAlignmentsFile = NULL;
    char *constraintAlignmentsFile = NULL;
    char *speciesTree = NULL;
    char *outgroupEvents = NULL;
    char *referenceEventString = NULL;
    bool runChecks = 0;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    //sleep(10);
    //assert(0);

    if(argc <= 1) {
        usage();
        return 0;
    }

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                { "params", required_argument, 0, 'p' },
                { "outputFile", required_argument, 0, 'f' },
                { "outputHalFastaFile", required_argument, 0, 'F' },
                { "outputReferenceFile", required_argument, 0, 'G' },
                { "sequences", required_argument, 0, 's' },
                { "alignments", required_argument, 0, 'a' },
                { "secondaryAlignments", required_argument, 0, 'S' },
                { "speciesTree", required_argument, 0, 'g' },
                { "constraintAlignments", required_argument, 0, 'c' },
                { "outgroupEvents", required_argument, 0, 'o' },
                { "help", no_argument, 0, 'h' },
                { "referenceEvent", required_argument, 0, 'r' },
                { "runChecks", no_argument, 0, 't' },
                { "threads", required_argument, 0, 'T' }, 
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int64_t key = getopt_long(argc, argv, "l:p:s:a:S:c:g:o:hr:F:G:tT:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'p':
                paramsFile = optarg;
                break;
            case 'f':
                outputFile = optarg;
                break;
            case 'F':
                outputHalFastaFile = optarg;
                break;
            case 'G':
                outputReferenceFile = optarg;
                break;
            case 's':
                sequenceFilesAndEvents = optarg;
                break;
            case 'a':
                alignmentsFile = optarg;
                break;
            case 'S':
                secondaryAlignmentsFile = optarg;
                break;
            case 'c':
                constraintAlignmentsFile = optarg;
                break;
            case 'g':
                speciesTree = optarg;
                break;
            case 'o':
                outgroupEvents = optarg;
                break;
            case 'r':
                referenceEventString = optarg;
                break;
            case 't':
                runChecks = 1;
                break;
            case 'T':
            {
                int num_threads = 0;
                int si = sscanf(optarg, "%d", &num_threads);
                assert(si == 1 && num_threads > 0);
                omp_set_num_threads(num_threads);
                break;
            }
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    if (paramsFile == NULL) {
        st_errAbort("must supply --params (-p)");
    }
    if (outputFile == NULL) {
        st_errAbort("must supply --outputFile (-f))");
    }
    if (sequenceFilesAndEvents == NULL) {
        st_errAbort("must supply --sequences (-s)");
    }
    if (alignmentsFile == NULL) {
        st_errAbort("must supply --alignments (-a)");
    }
    if (speciesTree == NULL) {
        st_errAbort("must supply --speciesTree (-f)");
    }
    if (referenceEventString == NULL) {
        st_errAbort("must supply --referenceEvent (-r)");
    }

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_logInfo("Params file: %s\n", paramsFile);
    st_logInfo("Output file string : %s\n", outputFile);
    st_logInfo("Output hal fasta file string : %s\n", outputHalFastaFile);
    st_logInfo("Output reference fasta file string : %s\n", outputReferenceFile);
    st_logInfo("Sequence files and events: %s\n", sequenceFilesAndEvents);
    st_logInfo("Alignments file: %s\n", alignmentsFile);
    st_logInfo("Secondary alignments file: %s\n", secondaryAlignmentsFile);
    st_logInfo("Constraint alignments file: %s\n", constraintAlignmentsFile);
    st_logInfo("Species tree: %s\n", speciesTree);
    st_logInfo("Outgroup events: %s\n", outgroupEvents);
    st_logInfo("Reference event: %s\n", referenceEventString);

    //////////////////////////////////////////////
    //Parse stuff
    //////////////////////////////////////////////

    // Load the params file
    CactusParams *params = cactusParams_load(paramsFile);
    st_logInfo("Loaded the parameters files, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    // Load the cactus disk
    CactusDisk *cactusDisk = cactusDisk_construct();

    st_logInfo("Set up the cactus disk, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //////////////////////////////////////////////
    //Call cactus setup
    //////////////////////////////////////////////

    Flower *flower = cactus_setup_first_flower(cactusDisk, params, speciesTree, outgroupEvents, sequenceFilesAndEvents);
    st_logInfo("Established the first Flower in the hierarchy, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    if(runChecks) {
        flower_checkRecursive(flower);
        st_logInfo("Checked the first flower in the hierarchy, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);
    }

    // Get the Name of the reference event - do this early so we don't fail late in the process
    Event *referenceEvent = eventTree_getEventByHeader(flower_getEventTree(flower), referenceEventString);
    if (referenceEvent == NULL) {
        st_errAbort("Reference event %s not found in tree. Check your "
                    "--referenceEventString option", referenceEventString);
    }
    Name referenceEventName = event_getName(referenceEvent);

    //////////////////////////////////////////////
    //Convert alignment coordinates
    //////////////////////////////////////////////

    alignmentsFile = convertAlignments(alignmentsFile, flower);
    if(secondaryAlignmentsFile != NULL) {
        secondaryAlignmentsFile = convertAlignments(secondaryAlignmentsFile, flower);
    }
    if(constraintAlignmentsFile != NULL) {
        constraintAlignmentsFile = convertAlignments(constraintAlignmentsFile, flower);
    }
    st_logInfo("Converted alignment coordinates, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //////////////////////////////////////////////
    //Strip the unique IDs
    //////////////////////////////////////////////

    stripUniqueIdsFromSequences(flower);
    st_logInfo("Stripped the unique IDs, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //////////////////////////////////////////////
    //Call cactus caf
    //////////////////////////////////////////////

    assert(!flower_builtBlocks(flower));
    caf(flower, params, alignmentsFile, secondaryAlignmentsFile, constraintAlignmentsFile);
    assert(flower_builtBlocks(flower));
    st_logInfo("Ran cactus caf, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    if(runChecks) {
        flower_checkRecursive(flower);
        st_logInfo("Checked the flowers in the hierarchy created by CAF, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);
    }

    //////////////////////////////////////////////
    //Call cactus bar
    //////////////////////////////////////////////

    if (cactusParams_get_int(params, 2, "bar", "runBar")) {
        stList *leafFlowers = stList_construct();
        extendFlowers(flower, leafFlowers, 1); // Get nested flowers to complete
        stList_sort(leafFlowers, flower_sizeCmpFn); // Sort by descending order of size, so that we start processing the
// largest flower as quickly as possible
        st_logInfo("Ran extended flowers ready for bar, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);


        bar(leafFlowers, params, cactusDisk, NULL);
        int64_t usePoa = cactusParams_get_int(params, 2, "bar", "partialOrderAlignment");
        st_logInfo("Ran cactus bar (use poa:%i), %" PRIi64 " seconds have elapsed\n", (int)usePoa, time(NULL) - startTime);

        stList_destruct(leafFlowers);

        if(runChecks) {
            flower_checkRecursive(flower);
            st_logInfo("Checked the flowers in the hierarchy created by BAR, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);
        }
    }

    //////////////////////////////////////////////
    //Call cactus reference
    //////////////////////////////////////////////

    // Get the flowers in the tree so that level 0 contains just the root flower,
    // level 1 contains the flowers that are children of the root flower, etc.
    stList *flowerLayers = getFlowerHierarchyInLayers(flower);
    for(int64_t i=0; i<stList_length(flowerLayers); i++) {
        stList_sort(stList_get(flowerLayers, i), flower_sizeCmpFn); // Sort by descending order of size, so that we start processing the
// largest flower as quickly as possible
    }
    st_logInfo("There are %" PRIi64 " layers in the flowers hierarchy\n", stList_length(flowerLayers));

    // Top-down this constructs the reference sequence
    for(int64_t i=0; i<stList_length(flowerLayers); i++) {
        stList *flowerLayer = stList_get(flowerLayers, i);
        st_logInfo("In the %" PRIi64 " layer there are %" PRIi64 " flowers in the flowers hierarchy\n", i,
                   stList_length(flowerLayer));
        cactus_make_reference(flowerLayer, referenceEventString, cactusDisk, params);
    }
    st_logInfo("Ran cactus make reference, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    // Bottom-up reference coordinates phase
    RecordHolder *rh = doBottomUpTraversal(flowerLayers, callBottomUp, (void *)referenceEventName);
    bottomUpNoDb(flower, rh, referenceEventName, 1, generateJukesCantorMatrix);
    assert(recordHolder_size(rh) == 0);
    recordHolder_destruct(rh);
    st_logInfo("Ran cactus make reference bottom up coordinates, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    // Top-down reference coordinates phase
    for(int64_t i=0; i<stList_length(flowerLayers); i++) {
        stList *flowers = stList_get(flowerLayers, i);
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
        for(int64_t j=0; j<stList_length(flowers); j++) {
            topDown(stList_get(flowers, j), referenceEventName);
        }
    }
    st_logInfo("Ran cactus make reference top down coordinates, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    if(runChecks) {
        flower_checkRecursive(flower);
        st_logInfo("Ran cactus check, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);
    }

    //////////////////////////////////////////////
    //Make c2h files, then build hal
    //////////////////////////////////////////////

    rh = doBottomUpTraversal(flowerLayers, callHalFn, (void *)referenceEventName);
    FILE *fileHandle = fopen(outputFile, "w");
    makeHalFormatNoDb(flower, rh, referenceEventName, fileHandle);
    fclose(fileHandle);
    assert(recordHolder_size(rh) == 0);
    recordHolder_destruct(rh);
    st_logInfo("Ran cactus to hal stage, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //////////////////////////////////////////////
    //Get reference sequences
    //////////////////////////////////////////////

    if(outputHalFastaFile != NULL) {
        fileHandle = fopen(outputHalFastaFile, "w");
        printFastaSequences(flower, fileHandle, referenceEventName);
        fclose(fileHandle);
        st_logInfo("Dumped sequences for hal file, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);
    }

    if(outputReferenceFile != NULL) {
        fileHandle = fopen(outputReferenceFile, "w");
        getReferenceSequences(fileHandle, flower, referenceEventString);
        fclose(fileHandle);
        st_logInfo("Dumped reference sequences, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);
    }

    //////////////////////////////////////////////
    //Cleanup
    //////////////////////////////////////////////

    st_system("rm %s", alignmentsFile);
    if(secondaryAlignmentsFile != NULL) {
        st_system("rm %s", secondaryAlignmentsFile);
    }
    if(constraintAlignmentsFile != NULL) {
        st_system("rm %s", constraintAlignmentsFile);
    }
    st_logInfo("Cactus consolidated is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0; // Exit without cleaning

    // Cleanup the memory
    stList_destruct(flowerLayers);
    cactusParams_destruct(params);
    cactusDisk_destruct(cactusDisk);

    st_logInfo("Cactus consolidated cleanup is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}
