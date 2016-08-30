#include "bioioC.h"
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

const char bases[4] = {'A', 'C', 'T', 'G'};

void usage() {
    fprintf(stderr, "--countThreshold - Minimum number of counts to include seed in global seed count table");
}

int stringDistance(char *a, char *b) {
    assert(strlen(a) == strlen(b));
    int distance = 0;
    for (int i = 0; i < strlen(a); i++) {
        if (a[i] != b[i]) {
            distance++;
        }
    }
    return distance;
}

void parseCounts(int optind, int argc, char **argv, stHash *seedCounts, stHash *seedToPacked) {
    char *seed = NULL;
    unsigned long packed;
    unsigned long count;
    for (int i = optind; i < argc; i++) {
        FILE *countsTable = fopen(argv[i], "r");

        char line[500];

        while (fgets(line, sizeof(line), countsTable) != NULL) {
            if (strlen(line) == 0) continue;
            seed = st_calloc(100, sizeof(char));

            sscanf(line, "%lx/%s %lu\n", &packed, seed, &count);
            if (strlen(seed) == 0) continue;
            //st_logInfo("Parsed seed %s, packed representation %lx, count %lu\n", seed, packed, count);
            //get rid of the colon
            seed[strlen(seed) - 1] = '\0';
            if (!stHash_search(seedCounts, seed)) {
                stHash_insert(seedCounts, seed, (void*)count);
            }
            else {
                unsigned long seedCount = (unsigned long)stHash_search(seedCounts, (void*)seed);
                stHash_insert(seedCounts, seed, (void*)(count + seedCount));
            }
            stHash_insert(seedToPacked, seed, packed);
        }
    }
}

stHash *blurCounts(stHash *seedCounts) {
    stHash *blurredCounts = stHash_construct();
    char *seed;
    stHashIterator *iter = stHash_getIterator(seedCounts);
    while((seed = stHash_getNext(iter))) {
        int seedLength = strlen(seed);
        int count = 0;
        for (int i = 0; i < seedLength; i++) {
            if (seed[i] == 'x') continue;
            char trueBase = seed[i];
            for (int base = 0; base < 4; base++) {
                if (seed[i] == bases[base]) continue;
                seed[i] = bases[base];
                int64_t deviantSeedCount = stHash_search(seedCounts, (void*)seed);
                count += deviantSeedCount;
            }
            seed[i] = trueBase;

        }
        stHash_insert(blurredCounts, seed, count);
    }
    stHash_destructIterator(iter);

    return blurredCounts;
}



int main(int argc, char **argv) {
    st_setLogLevelFromString("INFO");

    int64_t scoreThreshold = 0;
    char *seedScoresFile;
    bool clusterSeeds = false;
    struct option longopts[] = {{"scoreThreshold", required_argument, 0, 'c'},
    {"seedScoresFile", required_argument, 0, 'd'},
    {"clusterSeeds", no_argument, 0, 'e'},
    {0, 0, 0, 0} };
    int flag, k;
    while((flag = getopt_long(argc, argv, "c:d:e:", longopts, NULL)) != -1) {
        switch(flag) {
            case 'c':
                k = sscanf(optarg, "%" PRIi64 "", &scoreThreshold);
                break;
            case 'd':
                seedScoresFile = stString_copy(optarg);
                break;
            case 'e':
                clusterSeeds = true;
                break;
            case '?':
                break;
            default:
                usage();
                return 1;
                break;
        }
    }
    stHash *seedCounts = stHash_construct();
    stHash *seedToPacked = stHash_construct();
    parseCounts(optind, argc, argv, seedCounts, seedToPacked);

    stHash *seedScores;
    if (clusterSeeds) {
        seedScores = blurCounts(seedCounts);
    }
    else {
        seedScores = seedCounts;
    }
    FILE *seedScoresFileHandle = fopen(seedScoresFile, "w");
    stHashIterator *iter = stHash_getIterator(seedScores);
    char *seed;
    while((seed = stHash_getNext(iter)) != NULL) {
        int64_t score = stHash_search(seedScores, seed);
        if (score >= scoreThreshold) {
            fprintf(seedScoresFileHandle, "%lx/%s: %" PRIi64 "\n", stHash_search(seedToPacked, seed), seed, score);
        }
    }

}
