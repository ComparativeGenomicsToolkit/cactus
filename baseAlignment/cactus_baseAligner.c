#include <assert.h>
#include <getopt.h>

#include "cactus.h"
#include "sonLib.h"
#include "endAligner.h"
#include "flowerAligner.h"
#include "cactus_core.h"
#include "commonC.h"
#include "pairwiseAlignment.h"

#define SPANNING_TREES 5
#define MAXIMUM_LENGTH 1000

void usage() {
    fprintf(
            stderr,
            "cactus_baseAligner [flower-names, ordered by order they should be processed], version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

static stSortedSet *getAlignment_alignedPairs;
static stSortedSetIterator *getAlignment_iterator = NULL;

static struct PairwiseAlignment *getAlignments() {
    AlignedPair *alignedPair = stSortedSet_getNext(getAlignment_iterator);
    if(alignedPair == NULL) {
        return NULL;
    }
    struct List *opList = constructEmptyList(0, (void (*)(void *))destructAlignmentOperation);
    listAppend(opList, constructAlignmentOperation(PAIRWISE_MATCH, 1, 0.0));
    char *cA = cactusMisc_nameToString(alignedPair->sequence);
    char *cA2 = cactusMisc_nameToString(alignedPair->reverse->sequence);

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
    char * cactusDiskDatabaseString = NULL;
    int32_t j;

    /*
     * Parse the options.
     */
    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument, 0,
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
                cactusDiskDatabaseString = stString_copy(optarg);
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
    cCIP->deannealingRounds = 0;
    //--maxEdgeDegree 10000000 --minimumTreeCoverage 0 --minimumTreeCoverageForBlocks 0
    //--minimumBlockLength 0 --minimumChainLength 0 --trim 0 --alignRepeats 1 --extensionSteps 0

    /*
     * Load the flowerdisk
     */
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct2(kvDatabaseConf, 0, 1); //We precache the sequences
    st_logInfo("Set up the flower disk\n");

    /*
     * For each flower.
     */
    for (j = optind; j < argc; j++) {
        /*
         * Read the flower.
         */
        const char *flowerName = argv[j];
        st_logInfo("Processing the flower named: %s\n", flowerName);
        Flower *flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(flowerName));
        assert(flower != NULL);
        st_logInfo("Parsed the flower to be aligned\n");

        getAlignment_alignedPairs = makeFlowerAlignment(flower, SPANNING_TREES,
                MAXIMUM_LENGTH, &j);
        st_logInfo("Created the alignment: %i pairs\n", stSortedSet_size(getAlignment_alignedPairs));
        //getAlignment_alignedPairs = stSortedSet_construct();
        //assert(0);

        /*
         * Run the cactus core script.
         */
        cactusCorePipeline(flower, cCIP, getAlignments,
                startAlignmentStack, 1);
        st_logInfo("Ran the cactus core script.\n");

        /*
         * Cleanup
         */
        stSortedSet_destruct(getAlignment_alignedPairs);
        assert(getAlignment_iterator != NULL);
        stSortedSet_destructIterator(getAlignment_iterator);
        getAlignment_iterator = NULL;
        st_logInfo("Finished filling in the alignments for the flower\n");
        flower_unloadParent(flower); //The parent should not have changed.
    }

    /*
     * Write and close the cactusdisk.
     */
    cactusDisk_write(cactusDisk);
    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    destructCactusCoreInputParameters(cCIP);
    st_logInfo("Finished with the flower disk for this flower.\n");

    return 0;
}

