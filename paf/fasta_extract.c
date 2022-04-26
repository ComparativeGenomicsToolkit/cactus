/*
 * fasta_extract: Extracts sequences from a fasta file according to a bed file.
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

static int64_t flank = 10;
static int64_t min_size = 100;

void usage() {
    fprintf(stderr, "fasta_extract [fasta_file]xN [options], version 0.1\n");
    fprintf(stderr, "Extracts subsequences from a fasta file according to intervals in a bed file.\n"
                    "To encode the subsequence information each faster header is appended |x,"
                    "where x is the start coordinate (0-based) of the extracted sequence in the original sequence\n");
    fprintf(stderr, "-i --bedFile : The input bed file for the regions to extract. If omitted then reads bed intervals from stdin\n");
    fprintf(stderr, "-o --outFile : The output fasta file. If omitted then sequences will be written to stdout\n");
    fprintf(stderr, "-f --flank : How much flanking sequence to include at each end of each extracted sequence, by default: %" PRIi64 "\n", flank);
    fprintf(stderr, "-m --minSize : The minimum size of a sequence (before adding the flanks) to extract, by default: %" PRIi64 "\n", min_size);
    fprintf(stderr, "-n --skipMissing : Skip bed intervals that reference missing sequences instead of causing an error\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

void report_interval(FILE *output, char *seq_name, int64_t start, int64_t end, stHash *sequences, stHash *sequenceLengths) {
    char *sequence = stHash_search(sequences, seq_name);
    assert(sequence != NULL);
    int64_t seq_length = (int64_t)stHash_search(sequenceLengths, seq_name);
    assert(0 <= start); assert(start <= end); assert(end <= seq_length);
    char *s = stString_getSubString(sequence, start, end-start);
    fprintf(output, ">%s|%" PRIi64 "|%" PRIi64 "\n%s\n", seq_name, seq_length, start, s);
    free(s);
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *bed_file = NULL;
    char *output_file = NULL;
    bool skipMissing = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "bedFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "flank", required_argument, 0, 'f' },
                                                { "minSize", required_argument, 0, 'm' },
                                                { "skipMissing", no_argument, 0, 'n' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:o:f:hnm:i:", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                bed_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'f':
                flank = atol(optarg);
                break;
            case 'm':
                min_size = atol(optarg);
                break;
            case 'n':
                skipMissing = 1;
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
    st_logInfo("Bed file : %s\n", bed_file);
    st_logInfo("Flank size : %" PRIi64 "\n", flank);
    st_logInfo("Minimum sequence size (minSize) : %" PRIi64 "\n", min_size);

    //////////////////////////////////////////////
    // Parse the sequences
    //////////////////////////////////////////////

    stHash *sequences = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    while(optind < argc) {
        char *seq_file = argv[optind++];
        st_logInfo("Parsing sequence file : %s\n", seq_file);
        FILE *seq_file_handle = fopen(seq_file, "r");
        fastaReadToFunction(seq_file_handle, sequences, fastaRead_readToMapFunction);
        fclose(seq_file_handle);
    }
    st_logInfo("Read %i sequences from sequence files\n", (int)stHash_size(sequences));

    stHash *sequenceLengths = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    stHashIterator *it = stHash_getIterator(sequences);
    char *seq_name;
    while((seq_name = stHash_getNext(it)) != NULL) {
        // Add the sequence length (storing an integer directly within the hashtable as opposed to a pointer)
        stHash_insert(sequenceLengths, seq_name, (void *)strlen((char *)stHash_search(sequences, seq_name)));
        // Check we get back what we expect
        assert((int64_t)stHash_search(sequenceLengths, seq_name) == strlen((char *)stHash_search(sequences, seq_name)));
    }
    stHash_destructIterator(it);

    //////////////////////////////////////////////
    // Now parse the bed file and extract the sequences, ensuring they are non-overlapping
    //////////////////////////////////////////////

    FILE *input = bed_file == NULL ? stdin : fopen(bed_file, "r");
    FILE *output = output_file == NULL ? stdout : fopen(output_file, "w");
    char *line, *p_seq_name = NULL;
    int64_t p_start, p_end;
    while((line = stFile_getLineFromFile(input)) != NULL) {
        stList *tokens = stString_split(line); free(line); // Tokenize and cleanup the line
        char *seq_name = stList_get(tokens, 0);
        int64_t start = atol(stList_get(tokens, 1));
        int64_t end = atol(stList_get(tokens, 2));

        if(stHash_search(sequences, seq_name) == NULL) { // Check if we have the sequence
            if(skipMissing) { // If we are expecting to skip missing sequences
                stList_destruct(tokens);
                continue;
            }
            st_errAbort("Missing sequence: %s\n", seq_name); // Report an error if the sequence has not been read
        }


        if(end - start >= min_size) { // Meets the minimum threshold length
            st_logDebug("Processing sequence fragment: %s start:%" PRIi64 " end:%" PRIi64 " \n", seq_name, start, end);

            // Get coordinates of interval, adding on flanks
            int64_t seq_length = (int64_t)stHash_search(sequenceLengths, seq_name);
            int64_t i = start - flank > 0 ? start - flank : 0, j = end + flank <= seq_length ? end + flank : seq_length; // Get expanded bounds of sequence with flanks
            assert(0 <= i); assert(i <= start); assert(start <= end); assert(end <= j); assert(j <= seq_length);

            if(p_seq_name != NULL) { // Has a prior interval
                if(strcmp(p_seq_name, seq_name) == 0 && p_end >= i) { // Previous interval overlaps this one
                    stList_destruct(tokens);
                    p_end = p_end > j ? p_end : j; // Get right-most point of new and old interval to expand previous interval
                    continue; // Go to next interval without creating a new one, no need to print it yet
                }
                else { // Does not overlap previous interval
                    // Print out the previous interval
                    report_interval(output, p_seq_name, p_start, p_end, sequences, sequenceLengths);
                    free(p_seq_name);
                }
            }
            // Create a new interval - either because there was no prior interval or because we have printed out the previous interval
            p_seq_name = stString_copy(seq_name);
            p_start = i; p_end = j;
        }
        stList_destruct(tokens);
    }
    if(p_seq_name != NULL) { // Has a final interval
        // Report the last interval
        report_interval(output, p_seq_name, p_start, p_end, sequences, sequenceLengths);
        free(p_seq_name);
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    stHash_destruct(sequences);
    stHash_destruct(sequenceLengths);
    if(bed_file != NULL) {
        fclose(input);
    }
    if(output_file != NULL) {
        fclose(output);
    }

    st_logInfo("Fasta extract is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);
    return 0;
}
