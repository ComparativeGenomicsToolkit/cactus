/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <dirent.h>
#include <math.h>
#include <ctype.h>

#include "bioioC.h"
#include "cactus.h"

void usage() {
    fprintf(stderr, "cactus_analyseAssembly [fastaFile]xN\n");
}

//We want to report number of sequences,
static stList *sequenceLengths;
static stList *repeatBaseCounts;
int64_t nCount;
static const char *fileNameForStats;
void setupStatsCollation(const char *fileName) {
    fileNameForStats = fileName;
    sequenceLengths = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    repeatBaseCounts = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    nCount = 0;
}

void processSequenceForStats(const char *fastaHeader, const char *string, int64_t length) {
    //Collate stats
    int64_t sequenceLength = strlen(string);
    stList_append(sequenceLengths, stIntTuple_construct1(sequenceLength));
    int64_t j=0;
    for(int64_t i=0; i<sequenceLength; i++) {
        bool isN = string[i] == 'N' || string[i] == 'n';
        j += (tolower(string[i]) == string[i] || isN) ? 1 : 0;
        nCount += isN ? 1 : 0;
    }
    stList_append(repeatBaseCounts, stIntTuple_construct1(j));
}

void cleanupAndReportStatsCollection() {
    //Collate stats
    int64_t totalSequences = stList_length(sequenceLengths);
    int64_t totalLength = 0;
    int64_t repeatBaseCount = 0;
    stList_sort(sequenceLengths, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    for(int64_t i=0; i<totalSequences; i++) {
        totalLength += stIntTuple_get(stList_get(sequenceLengths, i), 0);
        repeatBaseCount += stIntTuple_get(stList_get(repeatBaseCounts, i), 0);
    }
    int64_t medianSequenceLength = totalSequences > 0 ? stIntTuple_get(stList_get(sequenceLengths, totalSequences/2), 0) : 0;
    int64_t maxSequenceLength = totalSequences > 0 ? stIntTuple_get(stList_peek(sequenceLengths), 0) : 0;
    int64_t minSequenceLength = totalSequences > 0 ? stIntTuple_get(stList_get(sequenceLengths, 0), 0) : 0;
    int64_t n50 = 0;
    int64_t j=0;
    for(int64_t i=totalSequences-1; i>=0; i--) {
        n50 = stIntTuple_get(stList_get(sequenceLengths, i), 0);
        j += n50;
        if(j >= totalLength/2) {
            break;
        }
    }
    fprintf(stdout, "Input-sample: %s Total-sequences: %" PRIi64 " Total-length: %" PRIi64 " Proportion-repeat-masked: %f ProportionNs: %f Total-Ns: %" PRIi64 " N50: %" PRIi64 " Median-sequence-length: %" PRIi64 " Max-sequence-length: %" PRIi64 " Min-sequence-length: %" PRIi64 "\n",
            fileNameForStats, totalSequences, totalLength, ((double)repeatBaseCount)/totalLength, ((double)nCount)/totalLength, nCount, n50, medianSequenceLength, maxSequenceLength, minSequenceLength);
    //Cleanup
    stList_destruct(sequenceLengths);
    stList_destruct(repeatBaseCounts);
}


int main(int argc, char *argv[]) {
    if(argc == 1) {
        usage();
        return 0;
    }

    for (int64_t j = 1; j < argc; j++) {
        FILE *fileHandle;
        if (strcmp(argv[j], "-") == 0) {
            fileHandle = stdin;
        } else {
            fileHandle = fopen(argv[j], "r");
            if (fileHandle == NULL) {
                st_errnoAbort("Could not open input file %s", argv[j]);
            }
        }
        setupStatsCollation(argv[j]);
        fastaReadToFunction(fileHandle, processSequenceForStats);
        cleanupAndReportStatsCollection();
        fclose(fileHandle);
    }

    return 0;
}
