/*
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include "sonLib.h"
#include "cactus.h"
#include "cactus_params_parser.h"

/*
 * TODOs:
 *
 * refactor setup
 * refactor caf
 * refactor bar
 * refactor reference
 * refactor output
 *
 * // Setup a test which feeds in inputs and checks output...
 */


void usage() {
    fprintf(stderr, "cactus_consolidated, version 0.2\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-p --params : [Required] The cactus config file\n");
    fprintf(stderr, "-d --cactusDisk : [Required] The database config string\n");
    fprintf(stderr, "-s --sequences [Required] [eventName fastaFile/Directory]xN: [Required] The sequences\n");
    fprintf(stderr, "-a --alignments : [Required] The alignments file\n");
    fprintf(stderr, "-g --speciesTree : [Required] The species tree, which will form the skeleton of the event tree\n");
    fprintf(stderr, "-o --outgroupEvents : Leaf events in the species tree identified as outgroups\n");
    fprintf(stderr, "-i --makeEventHeadersAlphaNumeric : Remove non alpha-numeric characters from event header names\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int main(int argc, char *argv[]) {
    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *paramsFile = NULL;
    char *cactusDiskDatabaseString = NULL;
    char *sequenceFiles = NULL;
    char *alignmentsFile = NULL;
    char *speciesTree = NULL;
    char *outgroupEvents = NULL;
    bool makeEventHeadersAlphaNumeric = 0;

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
                { "speciesTree", required_argument, 0, 'g' }, { "outgroupEvents", required_argument, 0, 'o' },
                { "help", no_argument, 0, 'h' }, { "makeEventHeadersAlphaNumeric", no_argument, 0, 'j' },
                { 0, 0, 0, 0 } };

        int option_index = 0;

        int64_t key = getopt_long(argc, argv, "l:p:d:s:a:g:o:hj", long_options, &option_index);

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
                sequenceFiles = optarg;
                break;
            case 'a':
                alignmentsFile = optarg;
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
            case 'j':
                makeEventHeadersAlphaNumeric = 1;
                break;
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
    if (sequenceFiles == NULL) {
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
    st_logInfo("Sequences files: %s\n", sequenceFiles);
    st_logInfo("Alignments file: %s\n", alignmentsFile);
    st_logInfo("Species tree: %s\n", speciesTree);
    st_logInfo("Outgroup events: %s\n", outgroupEvents);

    //////////////////////////////////////////////
    //Parse stuff
    //////////////////////////////////////////////

    // Load the params file
    CactusParams *params = cactusParams_load(paramsFile);

    // Load the cactus disk


    //////////////////////////////////////////////
    //Call cactus setup
    //////////////////////////////////////////////

    //////////////////////////////////////////////
    //Call cactus caf
    //////////////////////////////////////////////

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
