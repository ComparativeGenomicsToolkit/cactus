#include <assert.h>
#include <getopt.h>

#include "cactus.h"
#include "sonLib.h"
#include "endAligner.h"
#include "netAligner.h"
#include "cactus_core.h"
#include "commonC.h"
#include "pairwiseAlignment.h"

#define SPANNING_TREES 5
#define MAXIMUM_LENGTH 1000

void usage() {
    fprintf(
            stderr,
            "cactus_colinearAligner [net-names, ordered by order they should be processed], version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --netDisk : The location of the net disk directory\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

static stSortedSet *getAlignment_alignedPairs;
static stSortedSetIterator *getAlignment_iterator = NULL;

static struct PairwiseAlignment *getAlignments() {
    AlignedPair *alignedPair = stSortedSet_getNext(getAlignment_iterator);
    if(alignedPair == NULL) {
        return NULL;
    }
    struct List *opList = constructEmptyList(0, NULL);
    listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, 1, 0.0));
    char *cA = netMisc_nameToString(alignedPair->sequence);
    char *cA2 = netMisc_nameToString(alignedPair->reverse->sequence);

    int32_t i = alignedPair->position;
    int32_t j = i + 1;
    if(!alignedPair->strand) {
        j = i;
        i = alignedPair->position+1;
    }

    int32_t k = alignedPair->reverse->position;
    int32_t l = k + 1;
    if(!alignedPair->reverse->strand) {
        l = k;
        k = alignedPair->reverse->position+1;
    }

    struct PairwiseAlignment *pairwiseAlignment =
            constructPairwiseAlignment(
                    cA, i, j, alignedPair->strand,
                    cA2, k, l, alignedPair->reverse->strand,
                    1.0, opList);
    free(cA);
    free(cA2);
    return pairwiseAlignment;
}

static void startAlignmentStack() {
    if (getAlignment_iterator != NULL) {
        stSortedSet_destructIterator(getAlignment_iterator);
        getAlignment_iterator = NULL;
    }
    getAlignment_iterator = stSortedSet_getIterator(getAlignment_alignedPairs);
}

int main(int argc, char *argv[]) {

    char * logLevelString = NULL;
    char * netDiskName = NULL;
    int32_t j;

    /*
     * Parse the options.
     */
    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "netDisk", required_argument, 0,
                'b' }, { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:h", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
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

    if (logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if (logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    /*
     * Setup the input parameters for cactus core.
     */
    CactusCoreInputParameters *cCIP = constructCactusCoreInputParameters();
    //--maxEdgeDegree 10000000 --minimumTreeCoverage 0 --minimumTreeCoverageForBlocks 0
    //--minimumBlockLength 0 --minimumChainLength 0 --trim 0 --alignRepeats 1 --extensionSteps 0

    /*
     * Load the netdisk
     */
    NetDisk *netDisk = netDisk_construct(netDiskName);
    st_logInfo("Set up the net disk\n");

    /*
     * For each net.
     */
    for (j = optind; j < argc; j++) {
        /*
         * Read the net.
         */
        const char *netName = argv[j];
        st_logInfo("Processing the net named: %s\n", netName);
        Net *net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
        assert(net != NULL);
        st_logInfo("Parsed the net to be aligned\n");

        getAlignment_alignedPairs = makeNetAlignment(net, SPANNING_TREES,
                MAXIMUM_LENGTH, &j);
        st_logInfo("Created the alignment: %i pairs\n", stSortedSet_size(getAlignment_alignedPairs));
        //getAlignment_alignedPairs = stSortedSet_construct();
        //assert(0);

        /*
         * Run the cactus core script.
         */
        cactusCorePipeline(net, cCIP, getAlignments,
                startAlignmentStack, 1);
        st_logInfo("Ran the cactus core script.\n");

        /*
         * Cleanup
         */
        stSortedSet_destruct(getAlignment_alignedPairs);
        assert(getAlignment_iterator != NULL);
        stSortedSet_destructIterator(getAlignment_iterator);
        getAlignment_iterator = NULL;
        st_logInfo("Finished filling in the alignments for the net\n");
    }

    /*
     * Write and close the netdisk.
     */
    netDisk_write(netDisk);
    netDisk_destruct(netDisk);
    destructCactusCoreInputParameters(cCIP);
    st_logInfo("Finished with the net disk for this net.\n");

    return 0;
}

