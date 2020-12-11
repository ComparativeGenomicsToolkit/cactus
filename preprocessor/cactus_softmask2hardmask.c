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
    fprintf(stderr, "cactus_softmask2hardmask [fastaFile]\n");
    fprintf(stderr, "-m --minLength N: Only mask intervals > Nbp\n");
}

static void hardmask(void* min_length_p, const char* name, const char* seq, int64_t length) {
    int64_t min_length = *(int64_t*)min_length_p;
    char* uc_seq = (char*)seq;
    int64_t start = -1;
    for (int64_t i = 0; i < length; ++i) {
        bool lower = islower(uc_seq[i]);
        if (lower) {
            if (start == -1) {
                start = i;
            }
        }
        if (start != -1) {
            int64_t end = -1;
            if (!lower) {
                end = i - 1;
            } else if (i == length - 1) {
                end = i;
            }
            if (end > start + min_length - 1) {
                for (int64_t j = start; j <= end; ++j) {
                    assert(islower(uc_seq[j]));
                    uc_seq[j] = 'N';
                }                
            }
            if (end != -1) {
                start = -1;
            }
        }
    }
    fastaWrite(uc_seq, (char*)name, stdout);
}

int main(int argc, char *argv[]) {

    int64_t min_length = 0;
    
    while (1) {
        static struct option long_options[] = { { "minLength", required_argument, 0, 'm' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "m:", long_options, &option_index);
        int i = 0;

        if (key == -1) {
            break;
        }

        switch (key) {
        case 'm':
            i = sscanf(optarg, "%" PRIi64 "", &min_length);
            assert(i == 1);
            break;
        default:
            usage();
            return 1;
        }
    }

    if(argc == 1) {
        usage();
        return 0;
    }

    for (int64_t j = optind; j < argc; j++) {
        FILE *fileHandle;
        if (strcmp(argv[j], "-") == 0) {
            fileHandle = stdin;
        } else {
            fileHandle = fopen(argv[j], "r");
            if (fileHandle == NULL) {
                st_errnoAbort("Could not open input file %s", argv[j]);
            }
        }
        fastaReadToFunction(fileHandle, &min_length, hardmask);
    }

    return 0;
}
