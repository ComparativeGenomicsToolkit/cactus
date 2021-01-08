/*
 * Released under the MIT license, see LICENSE.txt
 */

#include <time.h>
#include <getopt.h>
#include "sonLib.h"
#include "cactus.h"
#include "cactus_setup.h"
#include "stCaf.h"

/*
 * TODOs:
 *
 * integrate with cactus_workflow
 * test cactus_setup
 * test cactus_caf
 * setup a test for cactus_consolidate which feeds in inputs and checks for output...
 *
 * refactor bar
 * refactor reference
 * refactor output
 */

void usage() {
    fprintf(stderr, "cactus_consolidated, version 0.2\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-p --params : [Required] The cactus config file\n");
    fprintf(stderr, "-d --cactusDisk : [Required] The database config string\n");
    fprintf(stderr, "-s --sequences [Required] [eventName fastaFile/Directory]xN: [Required] The sequences\n");
    fprintf(stderr, "-a --alignments : [Required] The alignments file\n");
    fprintf(stderr, "-S --secondaryAlignments : The secondary alignments file\n");
    fprintf(stderr, "-c --constraintAlignments : The constraint alignments file\n");
    fprintf(stderr, "-g --speciesTree : [Required] The species tree, which will form the skeleton of the event tree\n");
    fprintf(stderr, "-o --outgroupEvents : Leaf events in the species tree identified as outgroups\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *paramsFile = NULL;
    char *cactusDiskDatabaseString = NULL;
    char *sequenceFilesAndEvents = NULL;
    char *alignmentsFile = NULL;
    char *secondaryAlignmentsFile = NULL;
    char * constraintAlignmentsFile = NULL;
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
    stKVDatabaseConf *kvDatabaseConf = kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, true, true);
    st_logInfo("Set up the cactus disk, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //////////////////////////////////////////////
    //Call cactus setup
    //////////////////////////////////////////////

    Flower *flower = cactus_setup_first_flower(cactusDisk, params, speciesTree, outgroupEvents, sequenceFilesAndEvents);
    st_logInfo("Established the first Flower in the hierarchy, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //////////////////////////////////////////////
    //Call cactus caf
    //////////////////////////////////////////////

    caf(flower, params, alignmentsFile, secondaryAlignmentsFile, constraintAlignmentsFile);
    st_logInfo("Ran cactus caf, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //////////////////////////////////////////////
    //Call cactus bar
    //////////////////////////////////////////////

    //////////////////////////////////////////////
    //Call cactus reference
    //////////////////////////////////////////////

    //////////////////////////////////////////////
    //Cleanup
    //////////////////////////////////////////////

    cactusParams_destruct(params);

    return 0;
}
