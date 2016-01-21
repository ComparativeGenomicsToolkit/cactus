/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <getopt.h>
#include <time.h>

#include "cactus.h"
#include "sonLib.h"
#include "commonC.h"
#include "stCaf.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "stLastzAlignments.h"
#include "stGiantComponent.h"
#include "stCafPhylogeny.h"

static void usage() {
    fprintf(stderr, "cactus_caf, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --alignments : The input alignments file\n");
    fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-d --lastzArguments : Lastz arguments\n");
    fprintf(stderr, "-h --help : Print this help screen\n");

    fprintf(stderr, "-i --annealingRounds (array of ints, each greater than or equal to 1) : The rounds of annealing\n");
    fprintf(stderr,
            "-o --deannealingRounds (array of ints, each greater than or equal to 1 and each greater than the last) : The rounds of deannealing\n");

    fprintf(
            stderr,
            "-k --trim (array of integers, each greater or equal to zero) : An array giving the trim for each annealing round. If the array is shorter than the annealing rounds then a trim value of 0 is assumed for annealing rounds greater than the length of the trim array\n");

    fprintf(stderr,
            "-m --minimumTreeCoverage : (float [0.0, 1.0]) Minimum tree coverage proportion of a block to be included in the graph\n");

    fprintf(
            stderr,
            "-n --blockTrim : (int >= 0) The number of bases to trim from the ends of each block in a chain before accepting, this filtering is done after choosing the length of chains\n");

    fprintf(stderr, "-p --minimumDegree : (int >= 0) Minimum number of sequences in a block to be included in the output graph\n");

    fprintf(stderr, "-q --minimumIngroupDegree : Number of ingroup sequences required in a block.\n");

    fprintf(stderr, "-r --minimumOutgroupDegree : Number of outgroup sequences required in a block.\n");

    fprintf(stderr, "-s --singleCopyIngroup : Require that in-group sequences have only single coverage\n");

    fprintf(stderr, "-t --singleCopyOutgroup : Require that out-group sequences have only single coverage\n");

    fprintf(stderr, "-v --minimumSequenceLengthForBlast : The minimum length of a sequence to include when blasting\n");

    fprintf(
            stderr,
            "-w --maxAdjacencyComponentSizeRatio : The components equal or less than log(n) * of this size will be allowed in the cactus. Used to fight giant components.\n");

    fprintf(stderr, "-x --constraints : A file of alignments that will be enforced upon the cactus\n");

    fprintf(stderr,
            "-y --minLengthForChromosome : The minimum length required for a sequence to be considered as a candidate to be chromosome.\n");

    fprintf(
            stderr,
            "-z --proportionOfUnalignedBasesForNewChromosome : Proportion of aligned bases to be not contained in an existing chromosome to cause generation of a new chromosome.\n");
    fprintf(
                stderr,
                "-A --maximumMedianSequenceLengthBetweenLinkedEnds : Maximum nedian length of sequences between linked ends to allow before breaking chains.\n");
    fprintf(stderr, "-B --realign : Realign the lastz hits.\n");
    fprintf(stderr, "-C --realignArguments : Arguments for realignment.\n");
    fprintf(stderr, "-D --phylogenyNumTrees : Number of trees to sample when removing ancient homologies. (default 1)\n");
    fprintf(stderr, "-E --phylogenyRootingMethod : Method of rooting trees: either 'outgroupBranch', 'longestBranch', or 'bestRecon' (default outgroupBranch).\n");
    fprintf(stderr, "-F --phylogenyScoringMethod : Method of deciding which sampled tree is best: either 'reconCost' or .\n");
    fprintf(stderr, "-G --phylogenyBreakpointScalingFactor : scale breakpoint distance by this factor while building phylogenies. Default 0.0.\n");
    fprintf(stderr, "-H --phylogenySkipSingleCopyBlocks : Skip building trees for single-copy blocks. Default is not to skip.\n");
    fprintf(stderr, "-I --phylogenyMaxBaseDistance : maximum distance in bases to walk outside of a block gathering feature columns\n");
    fprintf(stderr, "-J --phylogenyMaxBlockDistance : maximum distance in blocks to walk outside of a block gathering feature columns\n");
    fprintf(stderr, "-K --phylogenyDebugFile : path to file to dump block trees and partitions to\n");
    fprintf(stderr, "-L --phylogenyKeepSingleDegreeBlocks : when splitting blocks, allow blocks to be created of only one ingroup.\n");
    fprintf(stderr, "-M --phylogenyTreeBuildingMethod : neighbor joining or neighbor-joining guided by the species tree\n");
    fprintf(stderr, "-N --phylogenyCostPerLossPerBase : join cost per dup per base for guided neighbor-joining (will be multiplied by maxBaseDistance)\n");
    fprintf(stderr, "-O --phylogenyCostPerLossPerBase : join cost per loss per base for guided neighbor-joining (will be multiplied by maxBaseDistance)\n");
    fprintf(stderr, "-P --referenceEventHeader : name of reference event (necessary for phylogeny estimation)\n");
    fprintf(stderr, "-Q --phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce : assume that this support value or greater means a very confident split, and that they will not be changed by the greedy split algorithm. Do all these very confident splits at once, to save a lot of computation time.\n");
    fprintf(stderr, "-R --numTreeBuildingThreads : Number of threads in the tree-building thread pool. Must be greater than 1. Default 2.\n");
    fprintf(stderr, "-S --phylogeny : Run the tree-building code and split ancient homologies away.\n");
    fprintf(stderr, "-T --minimumBlockHomologySupport: Minimum fraction of possible homologies required not to be considered a transitively collapsed megablock.\n");
    fprintf(stderr, "-U --phylogenyNucleotideScalingFactor: Weighting for the nucleotide information in the distance matrix used to build each tree.\n");
    fprintf(stderr, "-V --minimumBlockDegreeToCheckSupport: Minimum degree required to be checked for being a megablock.\n");
}

static int64_t *getInts(const char *string, int64_t *arrayLength) {
    int64_t *iA = st_malloc(sizeof(int64_t) * strlen(string));
    char *cA = stString_copy(string);
    char *cA2 = cA;
    char *cA3;
    *arrayLength = 0;
    while ((cA3 = stString_getNextWord(&cA)) != NULL) {
        int64_t i = sscanf(cA3, "%" PRIi64 "", &iA[(*arrayLength)++]);
        (void) i;
        assert(i == 1);
        free(cA3);
    }
    free(cA2);
    return iA;
}

static stList *getIntsAsList(const char *string) {
    stList *ret = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    stList *strings = stString_splitByString(string, " ");
    for (int64_t i = 0; i < stList_length(strings); i++) {
        int64_t num;
        int64_t success = sscanf(stList_get(strings, i), "%" PRIi64 "", &num);
        if (success != 1) {
            st_errAbort("error when reading integer list from string: %s", string);
        }
        stList_append(ret, stIntTuple_construct1(num));
    }
    stList_destruct(strings);
    return ret;
}

static int64_t minimumIngroupDegree = 0, minimumOutgroupDegree = 0, minimumDegree = 0;
static float minimumTreeCoverage = 0.0;
static Flower *flower = NULL;

static bool blockFilterFn(stPinchBlock *pinchBlock) {
    if ((minimumIngroupDegree > 0 || minimumOutgroupDegree > 0 || minimumDegree > 0) && !stCaf_containsRequiredSpecies(pinchBlock,
            flower, minimumIngroupDegree, minimumOutgroupDegree, minimumDegree)) {
        return 1;
    }
    if (minimumTreeCoverage > 0.0 && stCaf_treeCoverage(pinchBlock, flower) < minimumTreeCoverage) { //Tree coverage
        return 1;
    }
    return 0;
}

/*
 * Functions used for prefiltering the alignments.
 */

/*
 * Filtering by presence of outgroup. This code is efficient and scales linearly with depth.
 */

bool (*filterFn)(stPinchSegment *, stPinchSegment *) = NULL;
stSet *outgroupThreads = NULL;

bool containsOutgroupSegment(stPinchBlock *block) {
    stPinchBlockIt it = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&it)) != NULL) {
        //if(event_isOutgroup(getEvent(segment, flower))) {
        if (stSet_search(outgroupThreads, stPinchSegment_getThread(segment)) != NULL) {
            assert(event_isOutgroup(getEvent(segment, flower)));
            return 1;
        } else {
            assert(!event_isOutgroup(getEvent(segment, flower)));
        }
    }
    return 0;
}

bool isOutgroupSegment(stPinchSegment *segment) {
    if (stSet_search(outgroupThreads, stPinchSegment_getThread(segment)) != NULL) {
        assert(event_isOutgroup(getEvent(segment, flower)));
        return 1;
    }
    assert(!event_isOutgroup(getEvent(segment, flower)));
    return 0;
}

bool filterByOutgroup(stPinchSegment *segment1, stPinchSegment *segment2) {
    stPinchBlock *block1, *block2;
    if ((block1 = stPinchSegment_getBlock(segment1)) != NULL) {
        if ((block2 = stPinchSegment_getBlock(segment2)) != NULL) {
            if (block1 == block2) {
                return stPinchBlock_getLength(block1) == 1 ? 0 : containsOutgroupSegment(block1);
            }
            if (stPinchBlock_getDegree(block1) < stPinchBlock_getDegree(block2)) {
                return containsOutgroupSegment(block1) && containsOutgroupSegment(block2);
            }
            return containsOutgroupSegment(block2) && containsOutgroupSegment(block1);
        }
        return isOutgroupSegment(segment2) && containsOutgroupSegment(block1);
    }
    if ((block2 = stPinchSegment_getBlock(segment2)) != NULL) {
        return isOutgroupSegment(segment1) && containsOutgroupSegment(block2);
    }
    return isOutgroupSegment(segment1) && isOutgroupSegment(segment2);
}

/*
 * Filtering by presence of repeat species in block. This code is inefficient and does not scale.
 */

static bool checkIntersection(stSortedSet *names1, stSortedSet *names2) {
    stSortedSet *n12 = stSortedSet_getIntersection(names1, names2);
    bool b = stSortedSet_size(n12) > 0;
    stSortedSet_destruct(names1);
    stSortedSet_destruct(names2);
    stSortedSet_destruct(n12);
    return b;
}

static stSortedSet *getNames(stPinchSegment *segment) {
    stSortedSet *names = stSortedSet_construct();
    if (stPinchSegment_getBlock(segment) != NULL) {
        stPinchBlock *block = stPinchSegment_getBlock(segment);
        stPinchBlockIt it = stPinchBlock_getSegmentIterator(block);
        while ((segment = stPinchBlockIt_getNext(&it)) != NULL) {
            stSortedSet_insert(names, getEvent(segment, flower));
        }
    } else {
        stSortedSet_insert(names, getEvent(segment, flower));
    }
    return names;
}

bool filterByRepeatSpecies(stPinchSegment *segment1, stPinchSegment *segment2) {
    return checkIntersection(getNames(segment1), getNames(segment2));
}

static uint64_t choose2(uint64_t n) {
#define CHOOSE_TWO_CACHE_LEN 256
    // This is filled with 0's at init time by the compiler.
    static int64_t chooseTwoCache[256];
    if (n <= 1) {
        return 0;
    } else if (n >= CHOOSE_TWO_CACHE_LEN) {
        return n * (n - 1) / 2;
    } else {
        if (chooseTwoCache[n] != 0) {
            return chooseTwoCache[n];
        } else {
            chooseTwoCache[n] = n * (n - 1) / 2;
            return chooseTwoCache[n];
        }
    }
}

// Get the number of possible pairwise alignments that could support
// this block. Ordinarily this is (degree choose 2), but since we
// don't do outgroup self-alignment, it's a bit smaller.
static uint64_t numPossibleSupportingHomologies(stPinchBlock *block, Flower *flower) {
    uint64_t outgroupDegree = 0, ingroupDegree = 0;
    stPinchBlockIt segIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&segIt)) != NULL) {
        Name capName = stPinchSegment_getName(segment);
        Cap *cap = flower_getCap(flower, capName);
        Event *event = cap_getEvent(cap);
        if (event_isOutgroup(event)) {
            outgroupDegree++;
        } else {
            ingroupDegree++;
        }
    }
    assert(outgroupDegree + ingroupDegree == stPinchBlock_getDegree(block));
    // We do the ingroup-ingroup alignments as an all-against-all
    // alignment, so we can see each ingroup-ingroup homology up to
    // twice.
    return choose2(ingroupDegree) * 2 + ingroupDegree * outgroupDegree;
}

static void dumpBlockInfo(stPinchThreadSet *threadSet, const char *fileName) {
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    FILE *file = fopen(fileName, "w");
    if (file == NULL) {
        st_errnoAbort("couldn't open debug file");
    }
    stPinchBlock *block;
    while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
        uint64_t supportingHomologies = stPinchBlock_getNumSupportingHomologies(block);
        uint64_t possibleSupportingHomologies = numPossibleSupportingHomologies(block, flower);
        double support = ((double) supportingHomologies) / possibleSupportingHomologies;
        fprintf(file, "%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%lf\n", stPinchBlock_getDegree(block), stPinchBlock_getLength(block), supportingHomologies, possibleSupportingHomologies, support);
    }
    fclose(file);
}

// for printThreadSetStatistics
static int uint64_cmp(const uint64_t *x, const uint64_t *y) {
    if (*x < *y) {
        return -1;
    } else if (*x == *y) {
        return 0;
    } else {
        return 1;
    }
}

// for printThreadSetStatistics
static int double_cmp(const double *x, const double *y) {
    if (*x < *y) {
        return -1;
    } else if (*x == *y) {
        return 0;
    } else {
        return 1;
    }
}

// Print a set of statistics (avg, median, max, min) for degree and
// support percentage in the pinch graph.
static void printThreadSetStatistics(stPinchThreadSet *threadSet, Flower *flower, FILE *f)
{
    // Naively finds the max, median, and min by sorting: the lists
    // will have "only" millions of elements, so they should fit
    // comfortably into tens of MB of memory.

    uint64_t numBlocks = stPinchThreadSet_getTotalBlockNumber(threadSet);

    uint64_t *blockDegrees = malloc(numBlocks * sizeof(uint64_t));
    double totalDegree = 0.0;
    double *blockSupports = malloc(numBlocks * sizeof(double));
    double totalSupport = 0.0;

    uint64_t totalAlignedBases = 0;

    stPinchThreadSetBlockIt it = stPinchThreadSet_getBlockIt(threadSet);
    uint64_t i = 0;
    stPinchBlock *block;
    while ((block = stPinchThreadSetBlockIt_getNext(&it)) != NULL) {
        blockDegrees[i] = stPinchBlock_getDegree(block);
        totalDegree += stPinchBlock_getDegree(block);
        uint64_t supportingHomologies = stPinchBlock_getNumSupportingHomologies(block);
        uint64_t possibleSupportingHomologies = numPossibleSupportingHomologies(block, flower);
        double support = 0.0;
        if (possibleSupportingHomologies != 0) {
            support = ((double) supportingHomologies) / possibleSupportingHomologies;
        } else {
            assert(supportingHomologies == 0);
        }
        blockSupports[i] = support;
        totalSupport += support;

        totalAlignedBases += stPinchBlock_getLength(block) * stPinchBlock_getDegree(block);

        i++;
    }

    printf("There were %" PRIu64 " blocks in the sequence graph, representing %" PRIi64
           " total aligned bases\n", numBlocks, totalAlignedBases);

    qsort(blockDegrees, numBlocks, sizeof(uint64_t),
          (int (*)(const void *, const void *)) uint64_cmp);
    qsort(blockSupports, numBlocks, sizeof(double),
          (int (*)(const void *, const void *)) double_cmp);
    printf("Block degree stats: min %" PRIu64 ", avg %lf, median %" PRIu64 ", max %" PRIu64 "\n",
           blockDegrees[0], totalDegree/numBlocks, blockDegrees[(numBlocks - 1) / 2],
           blockDegrees[numBlocks - 1]);
    printf("Block support stats: min %lf, avg %lf, median %lf, max %lf\n",
           blockSupports[0], totalSupport/numBlocks, blockSupports[(numBlocks - 1) / 2],
           blockSupports[numBlocks - 1]);
    free(blockDegrees);
    free(blockSupports);
}

int main(int argc, char *argv[]) {
    /*
     * Script for adding alignments to cactus tree.
     */
    int64_t startTime;
    stKVDatabaseConf *kvDatabaseConf;
    CactusDisk *cactusDisk;
    int key, k;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * alignmentsFile = NULL;
    char * constraintsFile = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * lastzArguments = "";
    int64_t minimumSequenceLengthForBlast = 1;

    //Parameters for annealing/melting rounds
    int64_t *annealingRounds = NULL;
    int64_t annealingRoundsLength = 0;
    stList *meltingRounds = NULL;

    //Parameters for melting
    float maximumAdjacencyComponentSizeRatio = 10;
    int64_t blockTrim = 0;
    int64_t alignmentTrimLength = 0;
    int64_t *alignmentTrims = NULL;
    bool singleCopyIngroup = 0;
    bool singleCopyOutgroup = 0;
    int64_t chainLengthForBigFlower = 1000000;
    int64_t longChain = 2;
    int64_t minLengthForChromosome = 1000000;
    float proportionOfUnalignedBasesForNewChromosome = 0.8;
    bool breakChainsAtReverseTandems = 1;
    int64_t maximumMedianSequenceLengthBetweenLinkedEnds = INT64_MAX;
    bool realign = 0;
    char *realignArguments = "";

    //Parameters for removing ancient homologies
    bool doPhylogeny = false;
    int64_t phylogenyNumTrees = 1;
    enum stCaf_RootingMethod phylogenyRootingMethod = BEST_RECON;
    enum stCaf_ScoringMethod phylogenyScoringMethod = COMBINED_LIKELIHOOD;
    double breakpointScalingFactor = 1.0;
    bool phylogenySkipSingleCopyBlocks = 0;
    int64_t phylogenyMaxBaseDistance = 1000;
    int64_t phylogenyMaxBlockDistance = 100;
    bool phylogenyKeepSingleDegreeBlocks = 0;
    enum stCaf_TreeBuildingMethod phylogenyTreeBuildingMethod = GUIDED_NEIGHBOR_JOINING;
    double phylogenyCostPerDupPerBase = 0.2;
    double phylogenyCostPerLossPerBase = 0.2;
    const char *debugFileName = NULL;
    const char *referenceEventHeader = NULL;
    double phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce = 1.0;
    int64_t numTreeBuildingThreads = 2;
    int64_t minimumBlockDegreeToCheckSupport = 10;
    double minimumBlockHomologySupport = 0.7;
    double nucleotideScalingFactor = 1.0;
    stCaf_meltingMethod onlineMeltingMethod = REMOVE_NON_UNDOABLE_CHAINS;
    // Number of alignments to add before each melting step.
    int64_t numAlignmentsPerBatch = 1;
    // Stop after adding this many alignments. -1: never stop.
    int64_t maxNumAlignments = -1;
    // An alignment will be rejected if its fraction of "redundant"
    // pairs (pairs which are already present in the alignment) is
    // higher than this.
    double maxRedundantFraction = 1.0;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "alignments", required_argument, 0, 'b' }, {
                "cactusDisk", required_argument, 0, 'c' }, { "lastzArguments", required_argument, 0, 'd' },
                { "help", no_argument, 0, 'h' }, { "annealingRounds", required_argument, 0, 'i' }, { "trim", required_argument, 0, 'k' }, {
                        "trimChange", required_argument, 0, 'l', }, { "minimumTreeCoverage", required_argument, 0, 'm' }, { "blockTrim",
                        required_argument, 0, 'n' }, { "deannealingRounds", required_argument, 0, 'o' }, { "minimumDegree",
                        required_argument, 0, 'p' }, { "minimumIngroupDegree", required_argument, 0, 'q' }, {
                        "minimumOutgroupDegree", required_argument, 0, 'r' }, {
                        "singleCopyIngroup", no_argument, 0, 's' }, { "singleCopyOutgroup", no_argument, 0, 't' }, {
                        "minimumSequenceLengthForBlast", required_argument, 0, 'v' }, { "maxAdjacencyComponentSizeRatio",
                        required_argument, 0, 'w' }, { "constraints", required_argument, 0, 'x' }, { "minLengthForChromosome",
                        required_argument, 0, 'y' }, { "proportionOfUnalignedBasesForNewChromosome", required_argument, 0, 'z' },
                        { "maximumMedianSequenceLengthBetweenLinkedEnds", required_argument, 0, 'A' },
                        { "realign", no_argument, 0, 'B' }, { "realignArguments", required_argument, 0, 'C' },
                        { "phylogenyNumTrees", required_argument, 0, 'D' },
                        { "phylogenyRootingMethod", required_argument, 0, 'E' },
                        { "phylogenyScoringMethod", required_argument, 0, 'F' },
                        { "phylogenyBreakpointScalingFactor", required_argument, 0, 'G' },
                        { "phylogenySkipSingleCopyBlocks", no_argument, 0, 'H' },
                        { "phylogenyMaxBaseDistance", required_argument, 0, 'I' },
                        { "phylogenyMaxBlockDistance", required_argument, 0, 'J' },
                        { "phylogenyDebugFile", required_argument, 0, 'K' },
                        { "phylogenyKeepSingleDegreeBlocks", no_argument, 0, 'L' },
                        { "phylogenyTreeBuildingMethod", required_argument, 0, 'M' },
                        { "phylogenyCostPerDupPerBase", required_argument, 0, 'N' },
                        { "phylogenyCostPerLossPerBase", required_argument, 0, 'O' },
                        { "referenceEventHeader", required_argument, 0, 'P' },
                        { "phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce", required_argument, 0, 'Q' },
                        { "numTreeBuildingThreads", required_argument, 0, 'R' },
                        { "phylogeny", no_argument, 0, 'S' },
                        { "minimumBlockHomologySupport", required_argument, 0, 'T' },
                        { "phylogenyNucleotideScalingFactor", required_argument, 0, 'U' },
                        { "minimumBlockDegreeToCheckSupport", required_argument, 0, 'V' },
                        { "onlineMeltingMethod", required_argument, 0, 'W' },
                        { "maximumRedundantFraction", required_argument, 0, 'X' },
                        { 0, 0, 0, 0 } };

        int option_index = 0;

        key = getopt_long(argc, argv, "a:b:c:hi:k:m:n:o:p:q:r:stv:w:x:y:z:A:BC:D:E:F:G:HI:J:K:LM:N:O:P:Q:R:ST:U:V:W:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                st_setLogLevelFromString(logLevelString);
                break;
            case 'b':
                alignmentsFile = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'd':
                lastzArguments = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            case 'i':
                annealingRounds = getInts(optarg, &annealingRoundsLength);
                break;
            case 'o':
                meltingRounds = getIntsAsList(optarg);
                break;
            case 'k':
                alignmentTrims = getInts(optarg, &alignmentTrimLength);
                break;
            case 'm':
                k = sscanf(optarg, "%f", &minimumTreeCoverage);
                assert(k == 1);
                break;
            case 'n':
                k = sscanf(optarg, "%" PRIi64 "", &blockTrim);
                assert(k == 1);
                break;
            case 'p':
                k = sscanf(optarg, "%" PRIi64 "", &minimumDegree);
                assert(k == 1);
                break;
            case 'q':
                k = sscanf(optarg, "%" PRIi64 "", &minimumIngroupDegree);
                assert(k == 1);
                break;
            case 'r':
                k = sscanf(optarg, "%" PRIi64 "", &minimumOutgroupDegree);
                assert(k == 1);
                break;
            case 's':
                singleCopyIngroup = 1;
                break;
            case 't':
                singleCopyOutgroup = 1;
                break;
            case 'v':
                k = sscanf(optarg, "%" PRIi64 "", &minimumSequenceLengthForBlast);
                assert(k == 1);
                break;
            case 'w':
                k = sscanf(optarg, "%f", &maximumAdjacencyComponentSizeRatio);
                assert(k == 1);
                break;
            case 'x':
                constraintsFile = stString_copy(optarg);
                break;
            case 'y':
                k = sscanf(optarg, "%" PRIi64 "", &minLengthForChromosome);
                assert(k == 1);
                break;
            case 'z':
                k = sscanf(optarg, "%f", &proportionOfUnalignedBasesForNewChromosome);
                assert(k == 1);
                break;
            case 'A':
                k = sscanf(optarg, "%" PRIi64 "", &maximumMedianSequenceLengthBetweenLinkedEnds);
                assert(k == 1);
                break;
            case 'B':
                realign = 1;
                break;
            case 'C':
                realignArguments = stString_copy(optarg);
                break;
            case 'D':
                k = sscanf(optarg, "%" PRIi64, &phylogenyNumTrees);
                assert(k == 1);
                break;
            case 'E':
                if (!strcmp(optarg, "outgroupBranch")) {
                    phylogenyRootingMethod = OUTGROUP_BRANCH;
                } else if (!strcmp(optarg, "longestBranch")) {
                    phylogenyRootingMethod = LONGEST_BRANCH;
                } else if (!strcmp(optarg, "bestRecon")) {
                    phylogenyRootingMethod = BEST_RECON;
                } else {
                    st_errAbort("Invalid tree rooting method: %s", optarg);
                }
                break;
            case 'F':
                if (!strcmp(optarg, "reconCost")) {
                    phylogenyScoringMethod = RECON_COST;
                } else if (!strcmp(optarg, "nucLikelihood")) {
                    phylogenyScoringMethod = NUCLEOTIDE_LIKELIHOOD;
                } else if (!strcmp(optarg, "reconLikelihood")) {
                    phylogenyScoringMethod = RECON_LIKELIHOOD;
                } else if (!strcmp(optarg, "combinedLikelihood")) {
                    phylogenyScoringMethod = COMBINED_LIKELIHOOD;
                } else {
                    st_errAbort("Invalid tree scoring method: %s", optarg);
                }
                break;
            case 'G':
                k = sscanf(optarg, "%lf", &breakpointScalingFactor);
                assert(k == 1);
                break;
            case 'H':
                phylogenySkipSingleCopyBlocks = true;
                break;
            case 'I':
                k = sscanf(optarg, "%" PRIi64, &phylogenyMaxBaseDistance);
                assert(k == 1);
                break;
            case 'J':
                k = sscanf(optarg, "%" PRIi64, &phylogenyMaxBlockDistance);
                assert(k == 1);
                break;
            case 'K':
                debugFileName = stString_copy(optarg);
                break;
            case 'L':
                phylogenyKeepSingleDegreeBlocks = true;
                break;
            case 'M':
                if (strcmp(optarg, "neighborJoining") == 0) {
                    phylogenyTreeBuildingMethod = NEIGHBOR_JOINING;
                } else if (strcmp(optarg, "guidedNeighborJoining") == 0) {
                    phylogenyTreeBuildingMethod = GUIDED_NEIGHBOR_JOINING;
                } else {
                    st_errAbort("Unknown tree building method: %s", optarg);
                }
                break;
            case 'N':
                k = sscanf(optarg, "%lf", &phylogenyCostPerDupPerBase);
                assert(k == 1);
                break;
            case 'O':
                k = sscanf(optarg, "%lf", &phylogenyCostPerLossPerBase);
                assert(k == 1);
                break;
            case 'P':
                referenceEventHeader = stString_copy(optarg);
                break;
            case 'Q':
                k = sscanf(optarg, "%lf", &phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce);
                assert(k == 1);
                break;
            case 'R':
                k = sscanf(optarg, "%" PRIi64, &numTreeBuildingThreads);
                assert(k == 1);
                break;
            case 'S':
                doPhylogeny = true;
                break;
            case 'T':
                k = sscanf(optarg, "%lf", &minimumBlockHomologySupport);
                assert(k == 1);
                assert(minimumBlockHomologySupport <= 1.0);
                assert(minimumBlockHomologySupport >= 0.0);
                break;
            case 'U':
                k = sscanf(optarg, "%lf", &nucleotideScalingFactor);
                assert(k == 1);
                break;
            case 'V':
                k = sscanf(optarg, "%" PRIi64, &minimumBlockDegreeToCheckSupport);
                assert(k == 1);
                break;
            case 'W':
                // online melting method: in the form "method|numAlignmentsPerBatch|maxNumAlignments"
                ; // so we can declare variables at the "start" of a case statement
                stList *tokens = stString_splitByString(optarg, "|");
                char *method = stList_get(tokens, 0);
                if (stList_length(tokens) > 1) {
                    char *batchingStr = stList_get(tokens, 1);
                    if (sscanf(batchingStr, "%" PRIi64, &numAlignmentsPerBatch) != 1) {
                        st_errAbort("error reading the number of alignments per batch");
                    }
                    if (numAlignmentsPerBatch <= 0) {
                        st_errAbort("number of alignments per batch must be positive");
                    }
                }
                if (stList_length(tokens) > 2) {
                    char *maxNumAlignmentsStr = stList_get(tokens, 2);
                    if (sscanf(maxNumAlignmentsStr, "%" PRIi64, &maxNumAlignments) != 1) {
                        st_errAbort("error reading the maximum number of alignments to add");
                    }
                }
                if (strcmp(method, "preserveNonUndoableChains") == 0) {
                    onlineMeltingMethod = PRESERVE_NON_UNDOABLE_CHAINS;
                } else if (strcmp(method, "removeNonUndoableChains") == 0) {
                    onlineMeltingMethod = REMOVE_NON_UNDOABLE_CHAINS;
                } else if (strcmp(method, "onlyUndo") == 0) {
                    onlineMeltingMethod = ONLY_UNDO;
                } else if (strcmp(method, "onlyRemove") == 0) {
                    onlineMeltingMethod = ONLY_REMOVE;
                } else if (strcmp(method, "none") == 0) {
                    onlineMeltingMethod = NONE;
                } else {
                    st_errAbort("Unrecognized onlineMeltingMethod.");
                }
                stList_destruct(tokens);
                break;
            case 'X':
                if (sscanf(optarg, "%lf", &maxRedundantFraction) != 1) {
                    st_errAbort("error reading the maximumRedundantFraction");
                }
                break;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(cactusDiskDatabaseString != NULL);
    assert(minimumTreeCoverage >= 0.0);
    assert(minimumTreeCoverage <= 1.0);
    assert(blockTrim >= 0);
    assert(annealingRoundsLength >= 0);
    for (int64_t i = 0; i < annealingRoundsLength; i++) {
        assert(annealingRounds[i] >= 0);
    }
    for (int64_t i = 1; i < stList_length(meltingRounds); i++) {
        int64_t prevRoundLen = stIntTuple_get(stList_get(meltingRounds, i - 1), 0);
        int64_t curRoundLen = stIntTuple_get(stList_get(meltingRounds, i), 0);
        (void) prevRoundLen;
        (void) curRoundLen;
        assert(prevRoundLen < curRoundLen);
        assert(prevRoundLen >= 1);
    }
    assert(alignmentTrimLength >= 0);
    for (int64_t i = 0; i < alignmentTrimLength; i++) {
        assert(alignmentTrims[i] >= 0);
    }
    assert(minimumOutgroupDegree >= 0);
    assert(minimumIngroupDegree >= 0);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower disk name : %s\n", cactusDiskDatabaseString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Sort the constraints
    ///////////////////////////////////////////////////////////////////////////

    stPinchIterator *pinchIteratorForConstraints = NULL;
    if (constraintsFile != NULL) {
        pinchIteratorForConstraints = stPinchIterator_constructFromFile(constraintsFile);
        st_logInfo("Created an iterator for the alignment constaints from file: %s\n", constraintsFile);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Do the alignment
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);

    stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    if (alignmentsFile == NULL) {
        cactusDisk_preCacheStrings(cactusDisk, flowers);
    }
    char *tempFile1 = NULL;
    for (int64_t i = 0; i < stList_length(flowers); i++) {
        flower = stList_get(flowers, i);
        if (!flower_builtBlocks(flower)) { // Do nothing if the flower already has defined blocks
            st_logDebug("Processing flower: %lli\n", flower_getName(flower));

            //Set up the graph and cactus, and add the initial alignments
            stOnlineCactus *cactus;
            stPinchThreadSet *threadSet = stCaf_setupForOnlineCactus(flower, &cactus);

            //Build the set of outgroup threads
            outgroupThreads = stCaf_getOutgroupThreads(flower, threadSet);

            bool sortAlignments = 1;
            if (singleCopyIngroup) {
                filterFn = filterByRepeatSpecies;
            }
            else if (singleCopyOutgroup) {
                if (stSet_size(outgroupThreads) == 0) {
                    filterFn = NULL;
                } else {
                    filterFn = filterByOutgroup;
                }
            }

            //Setup the alignments
            stList *alignmentsList = NULL;
            if (alignmentsFile != NULL) {
                assert(i == 0);
                assert(stList_length(flowers) == 1);
                if (sortAlignments) {
                    tempFile1 = getTempFile();
                    stCaf_sortCigarsFileByScoreInDescendingOrder(alignmentsFile, tempFile1);
                    alignmentsFile = tempFile1;
                }
            } else {
                if (tempFile1 == NULL) {
                    tempFile1 = getTempFile();
                }
                alignmentsList = stCaf_selfAlignFlower(flower, minimumSequenceLengthForBlast, lastzArguments, realign, realignArguments, tempFile1);
                if (sortAlignments) {
                    stCaf_sortCigarsByScoreInDescendingOrder(alignmentsList);
                }
                st_logDebug("Ran lastz and have %" PRIi64 " alignments\n", stList_length(alignmentsList));
            }

            if (annealingRoundsLength > 1) {
                st_errAbort("'online' cactus not currently compatible with more than one "
                            "annealing round.");
            }
            int64_t annealingRound = 0;
            int64_t minimumChainLength = annealingRounds[annealingRound];
            int64_t alignmentTrim = annealingRound < alignmentTrimLength ? alignmentTrims[annealingRound] : 0;

            // Add constraints without checking the chain lengths.
            if (pinchIteratorForConstraints != NULL) {
                stCaf_anneal(threadSet, pinchIteratorForConstraints, filterFn);
            }

            // Do the annealing, checking that we are not creating
            // small chains at each step.
            if (annealingRound == 0) {
                stCaf_annealPreventingSmallChains(flower, threadSet, cactus, alignmentsFile, alignmentsList, alignmentTrim,
                                                  filterFn, meltingRounds, onlineMeltingMethod, numAlignmentsPerBatch, maxNumAlignments, maxRedundantFraction, stString_print("%s-cactusDump", debugFileName));
            }

            // Dump the block degree and length distribution to a file
            if (debugFileName != NULL) {
                dumpBlockInfo(threadSet, stString_print("%s-blockStats-preMelting", debugFileName));
            }

            printf("Sequence graph statistics after annealing:\n");
            printThreadSetStatistics(threadSet, flower, stdout);
            // Check for poorly-supported blocks--those that have
            // been transitively aligned together but with very
            // few homologies supporting the transitive
            // alignment. These "megablocks" can snarl up the
            // graph so that a lot of extra gets thrown away in
            // the first melting step.
            stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
            stPinchBlock *block;
            while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
                if (stPinchBlock_getDegree(block) > minimumBlockDegreeToCheckSupport) {
                    uint64_t supportingHomologies = stPinchBlock_getNumSupportingHomologies(block);
                    uint64_t possibleSupportingHomologies = numPossibleSupportingHomologies(block, flower);
                    double support = ((double) supportingHomologies) / possibleSupportingHomologies;
                    if (support < minimumBlockHomologySupport) {
                        fprintf(stdout, "Destroyed a megablock with degree %" PRIi64
                                " and %" PRIi64 " supporting homologies out of a maximum "
                                "of %" PRIi64 " (%lf%%).\n", stPinchBlock_getDegree(block),
                                supportingHomologies, possibleSupportingHomologies, support);
                        stPinchBlock_destruct(block);
                    }
                }
            }

            // Add in the caps for the conventional cactus code and disable online cactus.
            stCaf_disableAndCleanupOnlineCactus(threadSet, cactus);
            stCaf_addCapsToPinchGraph(flower, threadSet);

            // Do the final melting step.
            stCaf_melt(flower, threadSet, NULL, 0, minimumChainLength, breakChainsAtReverseTandems, maximumMedianSequenceLengthBetweenLinkedEnds);

            // FIXME: remove or cleanup.
            FILE *cactusDumpFile = fopen(stString_print("%s-cactusDump", debugFileName), "a");
            fprintf(cactusDumpFile, "BATCH\t\t\t\t\t\t\t\t\t\t%" PRIi64 "\n", (int64_t) clock() / (CLOCKS_PER_SEC / 1000));
            dumpMaxBlockDegreeAndTAB(threadSet, cactusDumpFile);
            fprintf(cactusDumpFile, "BATCH\t\t\t\t\t\t\t\t\t\t%" PRIi64 "\n", (int64_t) clock() / (CLOCKS_PER_SEC / 1000));
            fclose(cactusDumpFile);

            printf("Sequence graph statistics after final melting round:\n");
            printThreadSetStatistics(threadSet, flower, stdout);
            //This does the filtering of blocks that do not have the required species/tree-coverage/degree.
            stCaf_joinTrivialBoundaries(threadSet);
            stCaf_melt(flower, threadSet, blockFilterFn, blockTrim, 0, 0, INT64_MAX);

            if (debugFileName != NULL) {
                dumpBlockInfo(threadSet, stString_print("%s-blockStats-postMelting", debugFileName));
            }

            printf("Sequence graph statistics after filtering:\n");
            printThreadSetStatistics(threadSet, flower, stdout);

            // Build a tree for each block, then use each tree to
            // partition the homologies between the ingroups sequences
            // into those that occur before the speciation with the
            // outgroup and those which occur late.
            if (stSet_size(outgroupThreads) > 0 && doPhylogeny) {
                st_logDebug("Starting to build trees and partition ingroup homologies\n");
                stHash *threadStrings = stCaf_getThreadStrings(flower, threadSet);
                st_logDebug("Got sets of thread strings and set of threads that are outgroups\n");
                FILE *debugFile = NULL;
                if (debugFileName != NULL) {
                    debugFile = fopen(stString_print("%s-phylogeny", debugFileName), "w");
                    if (debugFile == NULL) {
                        st_errnoAbort("could not open debug file");
                    }
                }
                stCaf_PhylogenyParameters params;
                params.treeBuildingMethod = phylogenyTreeBuildingMethod;
                params.rootingMethod = phylogenyRootingMethod;
                params.scoringMethod = phylogenyScoringMethod;
                params.breakpointScalingFactor = breakpointScalingFactor;
                params.nucleotideScalingFactor = nucleotideScalingFactor;
                params.skipSingleCopyBlocks = phylogenySkipSingleCopyBlocks;
                params.keepSingleDegreeBlocks = phylogenyKeepSingleDegreeBlocks;
                params.costPerDupPerBase = phylogenyCostPerDupPerBase;
                params.costPerLossPerBase = phylogenyCostPerLossPerBase;
                params.maxBaseDistance = phylogenyMaxBaseDistance;
                params.maxBlockDistance = phylogenyMaxBlockDistance;
                params.numTrees = phylogenyNumTrees;
                params.ignoreUnalignedBases = 1;
                params.onlyIncludeCompleteFeatureBlocks = 0;
                params.doSplitsWithSupportHigherThanThisAllAtOnce = phylogenyDoSplitsWithSupportHigherThanThisAllAtOnce;
                params.numTreeBuildingThreads = numTreeBuildingThreads;

                assert(params.numTreeBuildingThreads >= 1);

                stCaf_buildTreesToRemoveAncientHomologies(
                    threadSet, threadStrings, outgroupThreads, flower, &params,
                    debugFile, referenceEventHeader);
                if (debugFile != NULL) {
                    fclose(debugFile);
                }
                stHash_destruct(threadStrings);
                st_logDebug("Finished building trees\n");

                // Enforce the block constraints on minimum degree,
                // etc. after splitting.
                stCaf_melt(flower, threadSet, blockFilterFn, 0, 0, 0, INT64_MAX);
            }

            //Sort out case when we allow blocks of degree 1
            if (minimumDegree < 2) {
                st_logDebug("Creating degree 1 blocks\n");
                stCaf_makeDegreeOneBlocks(threadSet);
                stCaf_melt(flower, threadSet, blockFilterFn, blockTrim, 0, 0, INT64_MAX);
            } else if (maximumAdjacencyComponentSizeRatio < INT64_MAX) { //Deal with giant components
                st_logDebug("Breaking up components greedily\n");
                stCaf_breakupComponentsGreedily(threadSet, maximumAdjacencyComponentSizeRatio);
            }

            //Finish up
            stCaf_finish(flower, threadSet, chainLengthForBigFlower, longChain, minLengthForChromosome,
                    proportionOfUnalignedBasesForNewChromosome); //Flower is then destroyed at this point.
            st_logInfo("Ran the cactus core script\n");

            //Cleanup
            stPinchThreadSet_destruct(threadSet);
            stSet_destruct(outgroupThreads);

            if (alignmentsList != NULL) {
                stList_destruct(alignmentsList);
            }
            st_logInfo("Cleaned up from main loop\n");
        } else {
            st_logInfo("We've already built blocks / alignments for this flower\n");
        }
    }
    stList_destruct(flowers);
    if (tempFile1 != NULL) {
        st_system("rm %s", tempFile1);
    }

    if (constraintsFile != NULL) {
        stPinchIterator_destruct(pinchIteratorForConstraints);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Write the flower to disk.
    ///////////////////////////////////////////////////////////////////////////
    st_logDebug("Writing the flowers to disk\n");
    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk and %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);

    return 0; //Exit without clean up is quicker, enable cleanup when doing memory leak detection.

    //Destruct stuff
    startTime = time(NULL);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    free(cactusDiskDatabaseString);
    if (alignmentsFile != NULL) {
        free(alignmentsFile);
    }
    if (logLevelString != NULL) {
        free(logLevelString);
    }
    if (lastzArguments != NULL) {
        free(lastzArguments);
    }
    st_logInfo("Cleaned stuff up and am finished in: %" PRIi64 " seconds\n", time(NULL) - startTime);

    //while(1);
    return 0;
}
