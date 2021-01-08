/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include "bioioC.h"
#include "cactus.h"
#include "cactus_setup.h"

void usage() {
    fprintf(stderr, "cactus_setup, version 0.2\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-p --params : [Required] The cactus config file\n");
    fprintf(stderr, "-d --cactusDisk : [Required] The database config string\n");
    fprintf(stderr, "-s --sequences [Required] [eventName fastaFile/Directory]xN: [Required] The sequences\n");
    fprintf(stderr, "-g --speciesTree : [Required] The species tree, which will form the skeleton of the event tree\n");
    fprintf(stderr, "-o --outgroupEvents : Leaf events in the species tree identified as outgroups\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int main(int argc, char *argv[]) {
    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *paramsFile = NULL;
    char *cactusDiskDatabaseString = NULL;
    char *speciesTree = NULL;
    char *sequenceFilesAndEvents = NULL;
    char *outgroupEvents = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    if(argc <= 1) {
        usage();
        return 0;
    }

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "params", required_argument, 0, 'p' }, { "cactusDisk", required_argument, 0, 'd' },
                                                { "sequences", required_argument, 0, 's' },
                                                { "speciesTree", required_argument, 0, 'g' }, { "outgroupEvents", required_argument, 0, 'o' },
                                                { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int64_t key = getopt_long(argc, argv, "l:p:d:s:a:g:o:h", long_options, &option_index);

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
    st_logInfo("Species tree: %s\n", speciesTree);
    st_logInfo("Outgroup events: %s\n", outgroupEvents);

    //////////////////////////////////////////////
    //Load the params and database
    //////////////////////////////////////////////

    // Load the params file
    CactusParams *params = cactusParams_load(paramsFile);
    st_logInfo("Loaded the parameters files\n");

    stKVDatabaseConf *kvDatabaseConf = kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    if (stKVDatabaseConf_getType(kvDatabaseConf) == stKVDatabaseTypeTokyoCabinet || stKVDatabaseConf_getType(kvDatabaseConf)
            == stKVDatabaseTypeKyotoTycoon) {
        assert(stKVDatabaseConf_getDir(kvDatabaseConf) != NULL);
    }
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, true, true);
    st_logInfo("Set up the flower disk\n");

    //////////////////////////////////////////////
    //Populate the first flower
    //////////////////////////////////////////////

    Flower *flower = cactus_setup_first_flower(cactusDisk, params, speciesTree, outgroupEvents, sequenceFilesAndEvents);
    (void)flower; // prevent the unused assert

    ///////////////////////////////////////////////////////////////////////////
    // Write the flower to disk.
    ///////////////////////////////////////////////////////////////////////////

    //flower_check(flower);
    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    return 0;
}
