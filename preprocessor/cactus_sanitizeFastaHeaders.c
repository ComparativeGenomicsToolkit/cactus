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

    // we cut at whitespace (like preprocessor does by default)
    char* trunc_header = stString_copy(fastaHeader);
    size_t header_len = strlen(trunc_header);
    int64_t pipe_pos = -1;
    for (size_t i = 0; i < header_len; ++i) {
        if (trunc_header[i] == '|') {
            pipe_pos = (int64_t)i;
        } else if (isspace(trunc_header[i])) {
            trunc_header[i] = '\0';
            break;
        }
    }

    if (stSet_search(header_set, (void*)trunc_header) != NULL) {
        fprintf(stderr, "Error: The fasta header \"%s\" appears more than once for event \"%s\". Please ensure fast headers are unique for each input\n", trunc_header, eventName);
        exit(1);
    }

    stSet_insert(header_set, stString_copy(trunc_header));

    if (strncmp(fastaHeader, "id=", 3) != 0 || pipe_pos <= 0) { 
        // no prefix found, we add one
        char* prefixed_header = (char*)st_malloc(strlen(trunc_header) + strlen(eventName) + 8);
        sprintf(prefixed_header, "id=%s|%s", eventName, trunc_header);
        free(trunc_header);
        trunc_header = prefixed_header;
    }
    
    fastaWrite((char*)string, trunc_header, stdout);
    free(trunc_header);
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
