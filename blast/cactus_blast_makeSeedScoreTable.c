#include "bioioC.h"
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

void usage() {
    fprintf(stderr, "--countThreshold - Minimum number of counts to include seed in global seed count table");
}

int main(int argc, char **argv) {
    st_setLogLevelFromString("INFO");

    stHash *globalCounts = stHash_construct();
    int64_t countThreshold = 0;
    char *seedScoresFile;
    struct option longopts[] = { {"--countThreshold", required_argument, 0, 'c' },
                                 {"--seedScoresFile", required_argument, 0, 'd'},
                                 {0, 0, 0, 0} };
    int flag, k;
    while((flag = getopt_long(argc, argv, "c:d:e:", longopts, NULL)) != -1) {
        switch(flag) {
        case 'c':
            k = sscanf(optarg, "%" PRIi64 "", &countThreshold);
            break;
        case 'd':

            seedScoresFile = stString_copy(optarg);
            break;
        case '?':
        default:
            usage();
            return 1;
        }
    }

    char *seed;
    for (int i = optind; i < argc; i++) {
        FILE *countsTable = fopen(argv[i], "r");

        char line[500];

        while (fgets(line, sizeof(line), countsTable)) {
            seed = malloc(sizeof(char)*50);
            int count;
            sscanf(line, "%s %i\n", seed, &count);
            if (!stHash_search(globalCounts, seed)) {
                stHash_insert(globalCounts, seed, count);
            }
            else {
                int globalCount = stHash_search(globalCounts, seed);
                stHash_insert(globalCounts, seed, count + globalCount);
            }
        }
    }

    FILE *seedScoresFileHandle = fopen(seedScoresFile, "w");
    stHashIterator *iter = stHash_getIterator(globalCounts);
    while((seed = stHash_getNext(iter)) != NULL) {
        int count = stHash_search(globalCounts, seed);
        if (count >= countThreshold) {
            fprintf(globalCountsFileHandle, "%s %i\n", seed, stHash_search(globalCounts, seed));
        }
    }
                
}
