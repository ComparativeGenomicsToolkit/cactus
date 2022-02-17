/*
 * paf_view: Pretty print pafs
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
    fprintf(stderr, "paf_view [fasta_files]xN [options], version 0.1\n");
    fprintf(stderr, "Pretty print PAF alignments\n");
    fprintf(stderr, "-i --inputFile : Input paf file to invert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-a --includeAlignment : Include base level alignment in output\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
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

    if (optind >= argc) {
        fprintf(stderr, "Expected at least one sequence file\n");
        exit(1);
    }

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Input file string : %s\n", inputFile);
    st_logInfo("Output file string : %s\n", outputFile);

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

    //////////////////////////////////////////////
    // Shatter the paf records
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

    Paf *paf;
    while((paf = paf_read(input)) != NULL) {
        // Get the query sequence
        char *query_seq = stHash_search(sequences, paf->query_name);
        if(query_seq == NULL) {
            fprintf(stderr, "No query sequence named: %s found\n", paf->query_name);
            exit(1);
        }

        // Get the target sequence
        char *target_seq = stHash_search(sequences, paf->target_name);
        if(target_seq == NULL) {
            fprintf(stderr, "No target sequence named: %s found\n", paf->target_name);
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
    stHash_destruct(sequences);

    st_logInfo("Paf view is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}

