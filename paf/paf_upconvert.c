/*
 * paf_upconvert: Converts paf coordinates to refer to subsequences of the original sequences.
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <time.h>
#include "bioioC.h"
#include "commonC.h"
#include "sonLib.h"
#include "paf.h"

void usage() {
    fprintf(stderr, "paf_upconvert [fasta_file]xN [options], version 0.1\n");
    fprintf(stderr, "Converts the coordinates of paf alignments to refer to extracted subsequences.\n");
    fprintf(stderr, "-i --inFile : The input paf file. If omitted then reads pafs from stdin\n");
    fprintf(stderr, "-o --outFile : The output paf file. If omitted then pafs will be written to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

void fastaRead_readCoordinates(void *sequences, const char *fastaHeader, const char *sequence, int64_t length) {
    Interval *i = decode_fasta_header((char *)fastaHeader);
    i->end = i->start + strlen(sequence); // Get the length of the actual fragment
    st_logDebug("Adding sequence interval name: %s start:%" PRIi64 " end:%" PRIi64 " length: %" PRIi64 "\n",
            i->name, i->start, i->end, i->length);
    stList_append((stList *)sequences, i);
}

int cmp_overlapping_intervals(const void *i, const void *j) {
    Interval *x = (Interval *)i, *y = (Interval *)j;
    int k = strcmp(x->name, y->name);
    //st_uglyf("Hello %s %s %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64"\n", x->name, y->name, x->start, y->start, x->end, y->end);
    if(k != 0) { // If seqs have different names
        return k;
    }
    if(x->start < y->start) { // If query interval starts before the interval
        return -1;
    }
    if(x->start <= y->end) {
        assert(x->end <= y->end); // Must be contained within
        return 0;
    }
    return 1;
}

void fix_interval(stList *intervals, char **name, int64_t *start, int64_t *end, int64_t *length) {
    Interval qi; qi.name = *name; qi.start = *start; qi.end = *end;

    Interval *i = stList_binarySearch(intervals, &qi, cmp_overlapping_intervals);
    if(i != NULL) { // If this fails the coordinate range is not contained within a sequence interval
        st_logDebug("Found interval seq name: %s start:%" PRIi64 " end:%" PRIi64 "\n",
                i->name, i->start, i->length+i->start);

        // Fix query sequence name
        free(*name);
        *name = stString_print("%s|%" PRIi64 "|%" PRIi64 "", i->name, i->length, i->start); // Fix name
        *start -= i->start; *end -= i->start; *length = i->length; // Fix coordinates
    }
    else {
        st_logDebug("Did not find sequence for interval: seq: %s align start: %" PRIi64 " align end: %" PRIi64 "\n",
                    *name, *start, *end);
    }
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *paf_file = NULL;
    char *output_file = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:o:hi:", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                paf_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Paf file : %s\n", paf_file);

    //////////////////////////////////////////////
    // Parse the sequences
    //////////////////////////////////////////////

    stList *intervals = stList_construct();
    while(optind < argc) {
        char *seq_file = argv[optind++];
        st_logInfo("Parsing sequence file : %s\n", seq_file);
        FILE *seq_file_handle = fopen(seq_file, "r");
        fastaReadToFunction(seq_file_handle, intervals, fastaRead_readCoordinates);
        fclose(seq_file_handle);
    }
    stList_sort(intervals, cmp_intervals);
    st_logInfo("Read %i sequences from sequence files\n", (int)stList_length(intervals));

    //////////////////////////////////////////////
    // Now parse the bed file and extract the sequences, ensuring they are non-overlapping
    //////////////////////////////////////////////

    FILE *input = paf_file == NULL ? stdin : fopen(paf_file, "r");
    FILE *output = output_file == NULL ? stdout : fopen(output_file, "w");
    Paf *paf;
    while((paf = paf_read(input)) != NULL) {
        // fix query and target coordinates
        fix_interval(intervals, &(paf->query_name), &(paf->query_start), &(paf->query_end), &(paf->query_length));
        fix_interval(intervals, &(paf->target_name), &(paf->target_start), &(paf->target_end), &(paf->target_length));

        paf_check(paf); // Check all is okay

        paf_write(paf, output); // Write out the adjust paf
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    stList_destruct(intervals);
    if(paf_file != NULL) {
        fclose(input);
    }
    if(output_file != NULL) {
        fclose(output);
    }

    st_logInfo("Paf upconvert is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);
    return 0;
}
