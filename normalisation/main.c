#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "normal.h"

/*
 * This script is run bottom up on a cactus tree and ensures the tree is 'normalised' as we expect.
 */

void usage() {
    fprintf(stderr, "cactus_normalisation [net names], version 0.1\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    /*
     * Script for adding a reference genome to a net.
     */

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * netDiskName = NULL;
    int32_t j;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "netDisk", required_argument, 0,
                'c' }, { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:h", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'c':
                netDiskName = stString_copy(optarg);
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

    assert(logLevelString == NULL || strcmp(logLevelString, "CRITICAL") == 0 || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
    assert(netDiskName != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if (logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if (logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Netdisk name : %s\n", netDiskName);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    NetDisk *netDisk = netDisk_construct(netDiskName);
    st_logInfo("Set up the net disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Loop on the nets, doing the reference genome (this process must be run bottom up)
    ///////////////////////////////////////////////////////////////////////////

    for (j = optind; j < argc; j++) {
        /*
         * Read the net.
         */
        const char *netName = argv[j];
        st_logInfo("Processing the net named: %s\n", netName);
        Net *net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
        assert(net != NULL);
        st_logInfo("Parsed the net to normalise\n");

        /*
         * Now run the normalisation functions
         */
        chain_promoteChainsThatExtendHigherLevelChains(net);
        if (!net_deleteIfEmpty(net)) { //If we delete the net we need not run the remaining functions..
            makeTerminalNormal(net);
            net_removeIfRedundant(net);
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Write the net(s) back to disk.
    ///////////////////////////////////////////////////////////////////////////

    netDisk_write(netDisk);
    st_logInfo("Updated the net on disk\n");

    ///////////////////////////////////////////////////////////////////////////
    //Clean up.
    ///////////////////////////////////////////////////////////////////////////

    //Destruct stuff
    netDisk_destruct(netDisk);

    st_logInfo("Cleaned stuff up and am finished\n");
    return 0;
}
