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
}

static void hardmask(void* out_file, const char* name, const char* seq, int64_t length) {
    char* uc_seq = (char*)seq;
    for (int64_t i = 0; i < length; ++i) {
        if (islower(uc_seq[i])) {
            uc_seq[i] = 'N';
        }
    }
    fastaWrite(uc_seq, (char*)name, (FILE*)out_file);
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
        fastaReadToFunction(fileHandle, stdout, hardmask);
    }

    return 0;
}
