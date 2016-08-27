#include "bioioC.h"
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

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

char *extractSeed(char *seedLine) {
    unsigned long packed;
    char *seed = malloc(sizeof(char)*50);
    sscanf(seedLine, "%lx/%s", &packed, seed);
    return seed;
}
stHash *getCounts(int optind, int argc, char **argv) {
    stHash *globalCounts = stHash_construct();
    char *seed;
    for (int i = optind; i < argc; i++) {
        FILE *countsTable = fopen(argv[i], "r");

        char line[500];

        while (fgets(line, sizeof(line), countsTable) != NULL) {
            if (strlen(line) == 0) continue;
            seed = malloc(sizeof(char)*100);
            unsigned long count;
            sscanf(line, "%s %lu\n", seed, &count);
            if (strlen(seed) == 0) continue;
            if (!stHash_search(globalCounts, seed)) {
                stHash_insert(globalCounts, seed, (void*)count);
            }
            else {
                unsigned long globalCount = (unsigned long)stHash_search(globalCounts, (void*)seed);
                stHash_insert(globalCounts, seed, (void*)(count + globalCount));
            }
        }
    }
    return globalCounts;
}
stHash *createClusters(stList *seeds, int nClusters) {
    stHash *clusters = stHash_construct2(NULL, (void (*)(void *)) stList_destruct);
    int nSeeds = stList_length(seeds);
    int i = 0;
    while (i < nClusters) {
        int seedIndex = rand() % nSeeds;
        char *seed = stList_get(seeds, seedIndex);
        if (stHash_search(clusters, seed)) {
            continue;
        }
        else {
            stHash_insert(clusters, seed, stList_construct());
            i++;
        }
    }
    //stHash_destructIterator(it);
    return clusters;
}

char *findBestCluster(stHash *clusters, char *seed) {
    stHashIterator *it = stHash_getIterator(clusters);
    char *clusterSeed = NULL;
    int minDistance = 100;
    char *bestCluster = NULL;
    while((clusterSeed = stHash_getNext(it)) != NULL) {
        if (strlen(clusterSeed) == 0) continue;
        char *seedUnpacked = extractSeed(seed);
        char *clusterSeedUnpacked = extractSeed(clusterSeed);
        int distance = stringDistance(seedUnpacked, clusterSeedUnpacked);
        free(seedUnpacked);
        free(clusterSeedUnpacked);
        if (distance < minDistance) {
            minDistance = distance;
            bestCluster = clusterSeed;
        }
    }
    //st_logInfo("Found cluster %s for %s with distance %i\n", seed, bestCluster, minDistance);
    stHash_destructIterator(it);
    return bestCluster;
}
stHash *extractFromClusters(stHash *clusters, stHash *clusterMultiplicity) {
    stHash *seedToClusterMultiplicity = stHash_construct();
    stHashIterator *it = stHash_getIterator(clusters);
    char *cluster;
    while((cluster = stHash_getNext(it)) != NULL) {
        int multiplicity = stHash_search(clusterMultiplicity, cluster);
        stList *seedsInCluster = stHash_search(clusters, cluster);
        stListIterator *seedIterator = stList_getIterator(seedsInCluster);
        char *seed;
        while((seed = stList_getNext(seedIterator))) {
            stHash_insert(seedToClusterMultiplicity, seed, multiplicity);
        }
        stList_destructIterator(seedIterator);
    }
    stHash_destructIterator(it);
    return seedToClusterMultiplicity;
}


stHash *makeSeedClusters(stHash *seedCounts, int64_t nClusters) {
    stHash *clusterMultiplicity = stHash_construct();
    stHash *clusters = createClusters(stHash_getKeys(seedCounts), nClusters);
    stHashIterator *it = stHash_getIterator(seedCounts);
    char *seed;
    while((seed = stHash_getNext(it)) != NULL) {
        if (strlen(seed) == 0) continue;
        char *bestCluster = findBestCluster(clusters, seed);
        int64_t oldMultiplicity = stHash_search(clusterMultiplicity, bestCluster);
        stHash_insert(clusterMultiplicity, bestCluster, oldMultiplicity + stHash_search(seedCounts, seed));
        stList *seedsInCluster = stHash_search(clusters, bestCluster);
        stList_append(seedsInCluster, seed);
    }
    stHash *seedToClusterCount = extractFromClusters(clusters, clusterMultiplicity);
    return seedToClusterCount;

}
int main(int argc, char **argv) {
    st_setLogLevelFromString("INFO");

    int64_t scoreThreshold = 0;
    int64_t nClusters = 0;
    char *seedScoresFile;
    bool clusterSeeds = false;
    struct option longopts[] = {{"scoreThreshold", required_argument, 0, 'c'},
    {"seedScoresFile", required_argument, 0, 'd'},
    {"nClusters", required_argument, 0, 'e'},
    {"clusterSeeds", no_argument, 0, 'f'},
    {0, 0, 0, 0} };
    int flag, k;
    while((flag = getopt_long(argc, argv, "c:d:e:f:", longopts, NULL)) != -1) {
        switch(flag) {
            case 'c':
                k = sscanf(optarg, "%" PRIi64 "", &scoreThreshold);
                break;
            case 'd':
                seedScoresFile = stString_copy(optarg);
                break;
            case 'e':
                k = sscanf(optarg, "%" PRIi64 "", &nClusters);
                break;
            case 'f':
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
    stHash *seedCounts = getCounts(optind, argc, argv);

    stHash *seedScores;
    if (clusterSeeds) {
        seedScores = makeSeedClusters(seedCounts, nClusters);
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
            fprintf(seedScoresFileHandle, "%s %" PRIi64 "\n", seed, score);
        }
    }

}
