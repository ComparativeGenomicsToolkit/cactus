/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 *
 * REpeat Detector (Red) doesn't handle some input edge cases very well. These include
 * - tiny contigs
 * - contigs that are very low information -- ie nearly all the same base
 * 
 * This program can be used to filter them out before running Red then add them back in (-x) after
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
    fprintf(stderr, "cactus_redPrefilter [fastaFile]\n");
    fprintf(stderr, "-m --minLength N: Filter contigs < Nbp. DEFAULT=1000\n");
    fprintf(stderr, "-b --maxBaseFrac F:  Filter contigs with proportion of the same base >= F. DEFAULT=1.0\n");    
    fprintf(stderr, "-x --extract:     Extract (instead of remove) filtered sequences\n");
}

static bool extract = false;
static double max_base_frac = 1.0;
static int64_t run_len_threshold = 10; 

static void redfilter(void* min_length_p, const char* name, const char* seq, int64_t length) {
    int64_t min_length = *(int64_t*)min_length_p;
    char* uc_seq = (char*)seq;

    bool too_short = length < min_length;
    // this is a quick hack to filter out really low-entropy sequence
    // todo: would be better to get Red to handle without crashing,
    //       but I don't think they're bound to happen much in real assemblies
    int64_t base_hist[256] = {0};
    for (int64_t i = 0; i < length; ++i) {
        int64_t c = (int64_t)toupper(seq[i]);
        ++base_hist[c];
    }
    bool is_monomer = length == 0;
    char monomer;
    for (int64_t i = 0; i < 256 && !is_monomer; ++i) {
        if ((double)base_hist[i] / (double)length > max_base_frac) {
            is_monomer = true;
            monomer = (char)i;
        }
    }

    if (is_monomer && extract) {
        // softmask the monomer
        int64_t run_len = 0;
        for (int64_t i = 0; i < length; ++i) {
            if (toupper(uc_seq[i]) == monomer) {
                ++run_len;
            } else {
                run_len = 0;
            }
            if (run_len > run_len_threshold) {
                uc_seq[i] = tolower(uc_seq[i]);
            }
        }
    }

    if ((!extract && !too_short && !is_monomer) || (extract && (too_short || is_monomer))) {
        fastaWrite(uc_seq, (char*)name, stdout);
    }
}

int main(int argc, char *argv[]) {

    int64_t min_length = 1000;
    
    while (1) {
        static struct option long_options[] = { { "minLength", required_argument, 0, 'm' },
                                                { "maxBaseFrac", required_argument, 0, 'b' },            
                                                { "extract", no_argument, 0, 'x' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "m:b:xs", long_options, &option_index);
        int i = 0;

        if (key == -1) {
            break;
        }

        switch (key) {
        case 'm':
            i = sscanf(optarg, "%" PRIi64 "", &min_length);
            assert(i == 1);
            break;
        case 'b':
            i = sscanf(optarg, "%lf", &max_base_frac);
            assert(i == 1);
            break;            
        case 'x':
            extract = true;
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
        fastaReadToFunction(fileHandle, &min_length, redfilter);
    }

    return 0;
}
