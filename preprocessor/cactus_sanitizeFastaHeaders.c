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
static bool convert_range = false;

void addUniqueFastaPrefix(void* destination, const char *fastaHeader, const char *string, int64_t length) {

    char* eventName = (char*)destination;

    // these cause weird crashes in halAppendCactusSubtree -- until that's fixed there's no point in trying to support them
    if (length == 0) {
        fprintf(stderr, "Warning: ignoring empty fast a sequence \"%s\" from event \"%s\"\n", fastaHeader, eventName);
        return;
    }

    // pipeline does not support nameless sequences
    if (strlen(fastaHeader) == 0) {
        fprintf(stderr, "Error: empty fasta header (> with nothing after) found in event \"%s\"\n", eventName);
        exit(1);
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

    // FROM shared.common.get_faidx_subpath_rename_cmd() which should be kept consistent:
    //
    // transform chr1:10-15 (1-based inclusive) into chr1_sub_9_15 (0-based end open)
    // this is a format that contains no special characters in order to make assembly hubs
    // happy.  But it does require conversion going into vg which wants chr[9-15] and
    // hal2vg is updated to do this autmatically
    if (convert_range) {
        char *colonpos = strrchr(clipped_header, ':');
        char *dashpos = colonpos != NULL ? strrchr(clipped_header, '-') : NULL;
        if (colonpos != NULL && dashpos != NULL && dashpos > colonpos + 1) {
            bool are_numbers = true;
            for (char *i = colonpos + 1; i < dashpos && are_numbers; ++i) {
                are_numbers = isdigit(*i);
            }
            char *last_pos = clipped_header + strlen(clipped_header);
            for (char *i = dashpos + 1; i < last_pos && are_numbers; ++i) {
                are_numbers = isdigit(*i);
            }
            if (are_numbers) {
                int64_t range_start = -1;
                int64_t range_end = -1;
                int ret = sscanf(colonpos + 1, "%" PRIi64 "-%" PRIi64 "", &range_start, &range_end);
                if (ret != 2 || range_start < 0 || range_end < range_start) {
                    fprintf(stderr, "Error: Could not parse :start-end range in \"%s\" for event \"%s\".\n", clipped_header, eventName);
                    exit(1);
                }
                // samtools is 1-based but will treat 0 like 1 so we do the same while
                // converting here to 0-based
                if (range_start > 0) {
                    range_start -= 1;
                }
                char *range_header = (char*)st_calloc(last_pos - clipped_header + 5, sizeof(char));
                strcpy(range_header, clipped_header);
                sprintf(range_header + (colonpos - clipped_header), "_sub_%" PRIi64 "_%" PRIi64 "", range_start, range_end);
                free(clipped_header);
                clipped_header = range_header;
            }
        }
        if (colonpos != NULL) {
            // : characters apparently cause crashes in and of themselves
            int64_t n = strlen(clipped_header);
            for (int64_t i = 0; i < n; ++i) {
                if (clipped_header[i] == ':') {
                    clipped_header[i] = '_';
                }
            }
        }
    }

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
            convert_range = true;
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
    if (fileHandle != stdin) {
        fclose(fileHandle);
    }
    fflush(stdout);
    stSet_destruct(header_set);

    return 0;
}
