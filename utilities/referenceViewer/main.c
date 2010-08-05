#include "referenceViewer.h"

/*
 * The script builds a circos style plot of the reference structure of a net.
 * The format of the output graph is dot format.
 */

static void usage() {
    fprintf(stderr, "cactus_referenceViewer, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
    fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
    fprintf(stderr, "-e --outputFile : The file to write the dot graph file in.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    CactusDisk *netDisk;
    Flower *net;
    FILE *fileHandle;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * netDiskName = NULL;
    char * netName = NULL;
    char * outputFile = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while(1) {
        static struct option long_options[] = {
            { "logLevel", required_argument, 0, 'a' },
            { "netDisk", required_argument, 0, 'c' },
            { "netName", required_argument, 0, 'd' },
            { "outputFile", required_argument, 0, 'e' },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:e:h", long_options, &option_index);

        if(key == -1) {
            break;
        }

        switch(key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'c':
                netDiskName = stString_copy(optarg);
                break;
            case 'd':
                netName = stString_copy(optarg);
                break;
            case 'e':
                outputFile = stString_copy(optarg);
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

    assert(netDiskName != NULL);
    assert(netName != NULL);
    assert(outputFile != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Net disk name : %s\n", netDiskName);
    st_logInfo("Net name : %s\n", netName);
    st_logInfo("Output graph file : %s\n", outputFile);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    netDisk = cactusDisk_construct(netDiskName);
    st_logInfo("Set up the net disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    net = cactusDisk_getNet(netDisk, cactusMisc_stringToName(netName));
    Group *group = flower_getFirstGroup(net);
    if(group != NULL && !group_isLeaf(group)) {
    	net = group_getNestedNet(group);
    }
    st_logInfo("Parsed the top level net of the cactus tree to build\n");

    ///////////////////////////////////////////////////////////////////////////
    // Build the graph.
    ///////////////////////////////////////////////////////////////////////////

    fileHandle = fopen(outputFile, "w");
    graphViz_setupGraphFile(fileHandle);
    assert(flower_getReference(net) != NULL);
    makeReferenceGraph(flower_getReference(net), fileHandle);
    graphViz_finishGraphFile(fileHandle);
    fclose(fileHandle);
    st_logInfo("Written the reference graph to file\n");

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(netDisk);

    return 0;
}
