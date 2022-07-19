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

static int64_t max_gap_length = 10000;
static float percentage_to_trim = 0.02;
static int64_t chain_gap_open = 5000;
static int64_t chain_gap_extend = 1;

void usage() {
    fprintf(stderr, "paf_chain [options], version 0.1\n");
    fprintf(stderr, "Chains the records in the PAF file into chains, rescoring them as chains.\nChains are indicated with the cn tag.\n");
    fprintf(stderr, "-i --inputFile : Input paf file to invert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-g --maxGapLength [INT] : The maximum allowable length of a gap in either sequence to chain (default:%" PRIi64 "bp)\n", max_gap_length);
    fprintf(stderr, "-d --chainGapOpen [INT] : The cost of opening a chain gap (default:%" PRIi64 "bp)\n", chain_gap_open);
    fprintf(stderr, "-e --chainGapExtend [INT] : The cost of extending a chain gap (default:%" PRIi64 "bp)\n", chain_gap_extend);
    fprintf(stderr, "-t --trimFraction : Fraction (from 0 to 1) of aligned bases to discount from the ends of the alignments when chaining"
                    "to trim from each end of the alignment when chaining, allowing slightly overlapping alignments to be chained (default:%f)\n", percentage_to_trim);
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int64_t gap_cost(int64_t query_gap_length, int64_t target_gap_length, void *params) {
    assert(query_gap_length >= 0);
    assert(target_gap_length >= 0);
    // These parameters match the ones in lastz except the gap extend is 10 rather than 30
    return query_gap_length + target_gap_length == 0 ? 0 : chain_gap_open + chain_gap_extend * (query_gap_length + target_gap_length);
    //return 1000*(query_gap_length + target_gap_length);
    //int64_t min_indel = llabs(query_gap_length - target_gap_length);
    //int64_t diagonal_gap = query_gap_length < target_gap_length ? query_gap_length : target_gap_length;
    //return min_indel * 30 + 10 * diagonal_gap;
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
                                                { "maxGapLength", required_argument, 0, 'g' },
                                                { "trimFraction", required_argument, 0, 't' },
                                                { "chainGapOpen", required_argument, 0, 'd' },
                                                { "chainGapExtend", required_argument, 0, 'e' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:hg:t:d:e:", long_options, &option_index);
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
            case 'g':
                max_gap_length = atoi(optarg);
                break;
            case 't':
                percentage_to_trim = atof(optarg);
                break;
            case 'd':
                chain_gap_open = atoi(optarg);
                break;
            case 'e':
                chain_gap_extend = atoi(optarg);
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
    st_logInfo("Trim chained alignment ends by : %f %\n", percentage_to_trim);
    st_logInfo("Maximum gap length : %" PRIi64 "\n", max_gap_length);
    st_logInfo("Chain gap open : %" PRIi64 "\n", chain_gap_open);
    st_logInfo("Chain gap extend : %" PRIi64 "\n", chain_gap_extend);

    //////////////////////////////////////////////
    // Tile the paf records
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

    stList *pafs = read_pafs(input); // Load local alignments files (PAF)
    stList *chained_pafs = paf_chain(pafs, gap_cost, NULL, max_gap_length, percentage_to_trim); // Convert to set of chains

    // Output chained alignments file
    write_pafs(output, chained_pafs);

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    // Cleans up the pafs list - but not the pafs themselves, which are destroyed by chaining
    stList_setDestructor(pafs, NULL);
    stList_destruct(pafs);
    stList_destruct(chained_pafs);

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
