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
    fprintf(stderr, "cactus_addUniquePrefixes [fastaFile] [EVENT]\n");
}

static stSet* header_set = NULL;

void addUniqueFastaPrefix(void* destination, const char *fastaHeader, const char *string, int64_t length) {

    char* eventName = (char*)destination;

    // these cause weird crashes in halAppendCactusSubtree -- until that's fixed there's no point in trying to support them
    if (length == 0) {
        fprintf(stderr, "Warning: ignoring empty fast a sequence \"%s\" from event \"%s\"\n", fastaHeader, eventName);
        return;
    }

    if (stSet_search(header_set, (void*)fastaHeader) != NULL) {
        fprintf(stderr, "Error: The fasta header \"%s\" appears more than once for event \"%s\". Please ensure fast headers are unique for each input\n", fastaHeader, eventName);
        exit(1);
    }
    stSet_insert(header_set, stString_copy(fastaHeader));

    // we mimic a very basic check from the preprocessor while we're at it
    size_t header_len = strlen(fastaHeader);
    int64_t pipe_pos = -1;
    for (size_t i = 0; i < header_len; ++i) {
        if (fastaHeader[i] == '|') {
            pipe_pos = (int64_t)i;
        } else if (fastaHeader[i] == ' ' || fastaHeader[i] == '\t') {
            fprintf(stderr, "Error: The fasta header \"%s\" from event \"%s\" contains spaces or tabs. These characters will cause issues in space-separated formats like MAF, and may not function properly when viewed in a browser. Please remove these characters from the input headers and try again.\n", fastaHeader, eventName);
            exit(1);
        }
    }
  
    char* newHeader = NULL;
    if (strncmp(fastaHeader, "id=", 3) != 0 || pipe_pos <= 0) { 
        // no prefix found, we add one
        newHeader = (char*)malloc(strlen(fastaHeader) + strlen(eventName) + 8);
        sprintf(newHeader, "id=%s|%s", eventName, fastaHeader);
    }
    const char* header = newHeader ? newHeader : fastaHeader;
    fastaWrite((char*)string, (char*)header, stdout);
}

int main(int argc, char *argv[]) {
    if(argc != 3) {
        usage();
        return 1;
    }

    header_set = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    
    FILE *fileHandle;
    if (strcmp(argv[1], "-") == 0) {
        fileHandle = stdin;
    } else {
        fileHandle = fopen(argv[1], "r");
        if (fileHandle == NULL) {
            st_errAbort("Could not open input file %s", argv[1]);
        }
    }

    fastaReadToFunction(fileHandle, (void*)argv[2], addUniqueFastaPrefix);
    fclose(fileHandle);
    stSet_destruct(header_set);

    return 0;
}
