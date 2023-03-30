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
    fprintf(stderr, "cactus_santizeFastaHeaders [fastaFile] [EVENT] (-p)\n");
    fprintf(stderr, "-p: add pangenome-specific processing to strip everything before (up to) last occurrence of #\n");
}

static stSet* header_set = NULL;
static bool strip_pounds = false;

void addUniqueFastaPrefix(void* destination, const char *fastaHeader, const char *string, int64_t length) {

    char* eventName = (char*)destination;

    // these cause weird crashes in halAppendCactusSubtree -- until that's fixed there's no point in trying to support them
    if (length == 0) {
        fprintf(stderr, "Warning: ignoring empty fast a sequence \"%s\" from event \"%s\"\n", fastaHeader, eventName);
        return;
    }

    // we cut at whitespace (like preprocessor does by default)
    // optionally cut out up to last #
    int64_t start = 0;
    int64_t last = strlen(fastaHeader);

    for (size_t i = start; i < last; ++i) {
        if (isspace(fastaHeader[i])) {
            last = i;
            break;
        }
        if (strip_pounds && fastaHeader[i] == '#') {
            start = i + 1;
        }
    }

    char* clipped_header = stString_getSubString(fastaHeader, start, last-start);

    if (stSet_search(header_set, (void*)clipped_header) != NULL) {
        fprintf(stderr, "Error: The sanitzied fasta header \"%s\" appears more than once for event \"%s\". Please ensure fast headers are unique for each input\n", clipped_header, eventName);
        exit(1);
    }

    stSet_insert(header_set, stString_copy(clipped_header));
    
    if (strncmp(clipped_header, "id=", 3) != 0 || strchr(clipped_header, '|') == NULL) {
        // no prefix found, we add one
        char* prefixed_header = (char*)st_malloc(strlen(clipped_header) + strlen(eventName) + 8);
        sprintf(prefixed_header, "id=%s|%s", eventName, clipped_header);
        free(clipped_header);
        clipped_header = prefixed_header;
    }
    
    fastaWrite((char*)string, clipped_header, stdout);
    free(clipped_header);
}

int main(int argc, char *argv[]) {
    if(argc != 3 && argc != 4) {
        usage();
        return 1;
    }
    if (argc == 4) {
        if (strcmp(argv[3], "-p") == 0) {
            strip_pounds = true;
        } else {
            usage();
            return 1;
        }
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
