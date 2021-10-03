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

void usage() {
    fprintf(stderr, "paf_view query_sequence_fasta target_sequence_fasta [options], version 0.1\n");
    fprintf(stderr, "Pretty print PAF alignments\n");
    fprintf(stderr, "-i --inputFile : Input paf file to invert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
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
    char *query_sequence_file = NULL;
    char *target_sequence_file = NULL;

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
        int64_t key = getopt_long(argc, argv, "l:i:o:h", long_options, &option_index);
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
    st_logInfo("Input file string : %s\n", inputFile);
    st_logInfo("Output file string : %s\n", outputFile);
    st_logInfo("Query sequence file string : %s\n", query_sequence_file);
    st_logInfo("Target sequence file string : %s\n", target_sequence_file);

    //////////////////////////////////////////////
    // Parse the query and target
    //////////////////////////////////////////////

    fastaRead(query_sequence_file);

    //////////////////////////////////////////////
    // Shatter the paf records
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

    Paf *paf;
    while((paf = paf_read(input)) != NULL) {
        stList *matches = paf_shatter(paf);
        for(int64_t i=0; i<stList_length(matches); i++) {
            paf_write(stList_get(matches, i), output);
        }
        stList_destruct(matches);
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

    st_logInfo("Paf shatter is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}

