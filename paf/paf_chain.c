/*
 * paf_chain: Chain paf alignments
 *
 *  Released under the MIT license, see LICENSE.txt
 *
 * Overview:
 * (1) Load local alignment file (PAF)
 * (2) Sort alignments by chromosome and coordinate
 * (3) Chain alignments forward and reverse, assigning each alignment to a chain
 * (4) Output chained alignments file (PAF)
*/

#include "paf.h"
#include "inc/paf.h"
#include <getopt.h>
#include <time.h>

void usage() {
    fprintf(stderr, "paf_chain [options], version 0.1\n");
    fprintf(stderr, "Chains the records in the PAF file into longer, combined chain PAF records\n");
    fprintf(stderr, "-i --inputFile : Input paf file to invert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

static int intcmp(int64_t i, int64_t j) {
    return i > j ? 1 : (i < j ? -1 : 0);
}

static int paf_cmp_by_location(const void *a, const void *b) {
    Paf *p1 = (Paf *)a, *p2 = (Paf *)b;
    int i = strcmp(p1->query_name, p2->query_name);
    if(i == 0) {
        i = strcmp(p1->target_name, p2->target_name);
        if(i == 0) {
            i = intcmp(p1->query_start, p2->query_start);
            if(i == 0) {
                i = intcmp(p1->target_start, p2->target_start);
            }
        }
    }
    return i;
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;

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

    //////////////////////////////////////////////
    // Tile the paf records
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

    stList *pafs = read_pafs(input); // Load local alignments files (PAF)
    stList_sort(pafs, paf_cmp_by_location); // Sort alignments by chromosome and location

    // Chain alignments forward and reverse, assigning each alignment to a chain

    // For each alignment:

    // Find highest scoring chains that alignment could be chained with:

    // Update chain

    // For each alignment, in reverse order to initial chaining:

    // If not already assigned to an output chain:

    // Trace back chain, assigning alignments to that chain and add chain to output list

    // Output chained alignments file
    write_pafs(output, pafs);

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    stList_destruct(pafs);
    if(inputFile != NULL) {
        fclose(input);
    }
    if(outputFile != NULL) {
        fclose(output);
    }

    st_logInfo("Paf chain is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}

