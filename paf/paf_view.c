/*
 * paf_shatter: Break up paf alignments into individual matches
 *
 *  Released under the MIT license, see LICENSE.txt
 *
 * Overview:
 * (1) Load query and target sequences
 * (2) For each input PAF record pretty print the alignment
*/

#include "paf.h"
#include "inc/paf.h"
#include <getopt.h>
#include <time.h>
#include "bioioC.h"

void usage() {
    fprintf(stderr, "paf_view query_sequence_fasta target_sequence_fasta [options], version 0.1\n");
    fprintf(stderr, "Pretty print PAF alignments\n");
    fprintf(stderr, "-i --inputFile : Input paf file to invert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-a --includeAlignment : Include base level alignment in output\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

stHash *read_sequence_file(char *sequence_file) {
    FILE *seq_file_handle = fopen(sequence_file, "r");
    stHash *sequences = fastaReadToMap(seq_file_handle);
    fclose(seq_file_handle);
    st_logInfo("Read %i sequences from sequence file: %s\n", (int)stHash_size(sequences), sequence_file);
    return sequences;
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    bool include_alignment=0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:ha", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                inputFile = optarg;
                break;
            case 'o':
                outputFile = optarg;
                break;
            case 'a':
                include_alignment = 1;
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    if (optind >= argc-1) {
        fprintf(stderr, "Expected two arguments after options\n");
        exit(1);
    }

    char *query_sequence_file = argv[optind];
    char *target_sequence_file = argv[optind+1];

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Input file string : %s\n", inputFile);
    st_logInfo("Output file string : %s\n", outputFile);
    st_logInfo("Query sequence file string : %s\n", query_sequence_file);
    st_logInfo("Target sequence file string : %s\n", target_sequence_file);

    //////////////////////////////////////////////
    // Parse the query and target
    //////////////////////////////////////////////

    stHash *query_sequences = read_sequence_file(query_sequence_file);
    stHash *target_sequences = read_sequence_file(target_sequence_file);

    //////////////////////////////////////////////
    // Shatter the paf records
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

    Paf *paf;
    while((paf = paf_read(input)) != NULL) {
        // Get the query sequence
        char *query_seq = stHash_search(query_sequences, paf->query_name);
        if(query_seq == NULL) {
            fprintf(stderr, "No query sequence named: %s found in the query file: %s\n", paf->query_name, query_sequence_file);
            exit(1);
        }

        // Get the target sequence
        char *target_seq = stHash_search(target_sequences, paf->target_name);
        if(target_seq == NULL) {
            fprintf(stderr, "No target sequence named: %s found in the query file: %s\n", paf->target_name, target_sequence_file);
            exit(1);
        }

        // Now print the alignment
        paf_pretty_print(paf, query_seq, target_seq, output, include_alignment);

        // Cleanup
        paf_destruct(paf);
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if(inputFile != NULL) {
        fclose(input);
    }
    if(outputFile != NULL) {
        fclose(output);
    }
    stHash_destruct(query_sequences);
    stHash_destruct(target_sequences);

    st_logInfo("Paf view is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}

