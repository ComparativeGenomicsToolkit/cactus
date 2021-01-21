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

/*
 * TODOs:
 *
 * - add glue from cactus_workflow.py
 *
 * setup a test for cactus_consolidate which feeds in inputs and checks for output...
 * integrate with cactus_workflow
 *
 */

void usage() {
    fprintf(stderr, "cactus_consolidated, version 0.2\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-p --params : [Required] The cactus config file\n");
    fprintf(stderr, "-d --cactusDisk : [Required] The database config string\n");
    fprintf(stderr, "-e --outputDisk : [Required] The database config string to store the cactus to hal output\n");
    fprintf(stderr, "-f --outputFile : [Required] The file to write the combined cactus to hal output\n");
    fprintf(stderr, "-s --sequences [Required] [eventName fastaFile/Directory]xN: [Required] The sequences\n");
    fprintf(stderr, "-a --alignments : [Required] The alignments file\n");
    fprintf(stderr, "-S --secondaryAlignments : The secondary alignments file\n");
    fprintf(stderr, "-c --constraintAlignments : The constraint alignments file\n");
    fprintf(stderr, "-g --speciesTree : [Required] The species tree, which will form the skeleton of the event tree\n");
    fprintf(stderr, "-o --outgroupEvents : Leaf events in the species tree identified as outgroups\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

static char *convertAlignments(char *alignmentsFile, CactusDisk *cactusDisk) {
    char *tempFile = getTempFile();
    convertAlignmentCoordinates(alignmentsFile, tempFile, cactusDisk);
    return tempFile;
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *paramsFile = NULL;
    char *cactusDiskDatabaseString = NULL;
    char *outputDiskDatabaseString = NULL;
    char *outputFile = NULL;
    char *sequenceFilesAndEvents = NULL;
    char *alignmentsFile = NULL;
    char *secondaryAlignmentsFile = NULL;
    char *constraintAlignmentsFile = NULL;
    char *speciesTree = NULL;
    char *outgroupEvents = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    if(argc <= 1) {
        usage();
        return 0;
    }

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                { "params", required_argument, 0, 'p' }, { "cactusDisk", required_argument, 0, 'd' },
                { "outputDisk", required_argument, 0, 'e' },{ "outputFile", required_argument, 0, 'f' },
                { "sequences", required_argument, 0, 's' }, { "alignments", required_argument, 0, 'a' },
                { "secondaryAlignments", required_argument, 0, 'S' }, { "speciesTree", required_argument, 0, 'g' },
                { "constraintAlignments", required_argument, 0, 'c' },{ "outgroupEvents", required_argument, 0, 'o' },
                { "help", no_argument, 0, 'h' },{ 0, 0, 0, 0 } };

        int option_index = 0;

        int64_t key = getopt_long(argc, argv, "l:p:d:s:a:S:c:g:o:h", long_options, &option_index);

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
            case 'd':
                cactusDiskDatabaseString = optarg;
                break;
            case 'e':
                outputDiskDatabaseString = optarg;
                break;
            case 'f':
                outputFile = optarg;
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
    if (cactusDiskDatabaseString == NULL) {
        st_errAbort("must supply --cactusDisk (-d))");
    }
    if (outputDiskDatabaseString == NULL) {
        st_errAbort("must supply --outputDisk (-e))");
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

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_logInfo("Params file: %s\n", paramsFile);
    st_logInfo("Cactus disk database string : %s\n", cactusDiskDatabaseString);
    st_logInfo("Output disk database string : %s\n", outputDiskDatabaseString);
    st_logInfo("Output file string : %s\n", outputFile);
    st_logInfo("Sequence files and events: %s\n", sequenceFilesAndEvents);
    st_logInfo("Alignments file: %s\n", alignmentsFile);
    st_logInfo("Secondary alignments file: %s\n", secondaryAlignmentsFile);
    st_logInfo("Constraint alignments file: %s\n", constraintAlignmentsFile);
    st_logInfo("Species tree: %s\n", speciesTree);
    st_logInfo("Outgroup events: %s\n", outgroupEvents);

    //////////////////////////////////////////////
    //Parse stuff
    //////////////////////////////////////////////

    // Load the params file
    CactusParams *params = cactusParams_load(paramsFile);
    st_logInfo("Loaded the parameters files, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    // Load the cactus disk
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, true, true);
    st_logInfo("Set up the cactus disk, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    // Load the output disk
    //stKVDatabaseConf *kvOutputDatabaseConf = stKVDatabaseConf_constructFromString(outputDiskDatabaseString);
    stKVDatabase *outputDatabase = NULL; //stKVDatabase_construct(kvOutputDatabaseConf, 0);
    //stKVDatabaseConf_destruct(kvOutputDatabaseConf);
    st_logInfo("Set up the output disk, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //////////////////////////////////////////////
    //Call cactus setup
    //////////////////////////////////////////////

    Flower *flower = cactus_setup_first_flower(cactusDisk, params, speciesTree, outgroupEvents, sequenceFilesAndEvents);
    //stripUniqueIdsFromMetaSequences(flower); // Not clear if this is needed
    st_logInfo("Established the first Flower in the hierarchy, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    flower_checkRecursive(flower);

    //////////////////////////////////////////////
    //Convert alignment coordinates
    //////////////////////////////////////////////

    alignmentsFile = convertAlignments(alignmentsFile, cactusDisk);
    if(secondaryAlignmentsFile != NULL) {
        secondaryAlignmentsFile = convertAlignments(secondaryAlignmentsFile, cactusDisk);
    }
    if(constraintAlignmentsFile != NULL) {
        constraintAlignmentsFile = convertAlignments(constraintAlignmentsFile, cactusDisk);
    }
    st_logInfo("Converted alignment coordinates, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //////////////////////////////////////////////
    //Call cactus caf
    //////////////////////////////////////////////

    assert(!flower_builtBlocks(flower));
    caf(flower, params, alignmentsFile, secondaryAlignmentsFile, constraintAlignmentsFile);
    assert(flower_builtBlocks(flower));
    st_logInfo("Ran cactus caf, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    flower_checkRecursive(flower);

    //////////////////////////////////////////////
    //Call cactus bar
    //////////////////////////////////////////////

    stList *leafFlowers = stList_construct();
    extendFlowers(flower, leafFlowers, 1); // Get nested flowers to complete

    flower_checkRecursive(flower);

    st_logInfo("Ran extended flowers ready for bar, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);
    bar(leafFlowers, params, cactusDisk, NULL);
    st_logInfo("Ran cactus bar, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    flower_checkRecursive(flower);

    //////////////////////////////////////////////
    //Call cactus reference
    //////////////////////////////////////////////

    // Get the flowers in the tree so that level 0 contains just the root flower,
    // level 1 contains the flowers that are children of the root flower, etc.
    stList *flowerLayers = getFlowerHierarchyInLayers(flower);

    // Top-down this constructs the reference sequence
    for(int64_t i=0; i<stList_length(flowerLayers); i++) {
        cactus_make_reference(stList_get(flowerLayers, i), cactusDisk, params);
    }

    // Get the Name of the reference event
    char *referenceEventString = cactusParams_get_string(params, 2, "reference", "reference");
    Event *referenceEvent = eventTree_getEventByHeader(flower_getEventTree(flower), referenceEventString);
    if (referenceEvent == NULL) {
        st_errAbort("Reference event %s not found in tree. Check your "
                    "--referenceEventString option", referenceEventString);
    }
    Name referenceEventName = event_getName(referenceEvent);

    // Bottom-up reference coordinates phase
    for(int64_t i=stList_length(flowerLayers)-1; i>=0 ; i--) {
        bottomUp(stList_get(flowerLayers, i), outputDatabase, referenceEventName, i==0, generateJukesCantorMatrix);
    }

    // Top-down reference coordinates phase
    for(int64_t i=0; i<stList_length(flowerLayers); i++) {
        stList *flowers = stList_get(flowerLayers, i);
        for(int64_t j=0; j<stList_length(flowers); j++) {
            topDown(stList_get(flowers, j), referenceEventName);
        }
    }

    //////////////////////////////////////////////
    //Make c2h files, then build hal
    //////////////////////////////////////////////

    // Bottom-up reference coordinates phase
    for(int64_t i=stList_length(flowerLayers)-1; i>0 ; i--) {
        stList *flowers = stList_get(flowerLayers, i);
        for (int64_t j = 0; j < stList_length(flowers); j++) {
            makeHalFormat(stList_get(flowers, j), outputDatabase, referenceEventName, NULL);
        }
    }

    // Now write the complete cactus to hal file.
    FILE *fileHandle = fopen(outputFile, "w");
    makeHalFormat(flower, outputDatabase, referenceEventName, fileHandle);
    fclose(fileHandle);

    //////////////////////////////////////////////
    //Cleanup
    //////////////////////////////////////////////

    stList_destruct(flowerLayers);
    cactusParams_destruct(params);
    st_system("rm %s", alignmentsFile);
    if(secondaryAlignmentsFile != NULL) {
        st_system("rm %s", secondaryAlignmentsFile);
    }
    if(constraintAlignmentsFile != NULL) {
        st_system("rm %s", constraintAlignmentsFile);
    }

    // more to do here...

    return 0;
}
