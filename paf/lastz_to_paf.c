/*
 * lastz_to_paf: Lastz align to PAF converter:
 * (1) Read in alignments in lastz cigar format (from stdin)
 * (2) Convert to PAF
 * (3) Write out alignments
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "paf.h"
#include "pairwiseAlignment.h"
#include <getopt.h>
#include <time.h>

void usage() {
    fprintf(stderr, "lastz_to_paf [options], version 0.1\n");
    fprintf(stderr, "Converts lastz cigar format aligments to PAF\n");
    fprintf(stderr, "-i --inputFile : Input lastz cigar file to invert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

Paf *cigar_to_paf(struct PairwiseAlignment *pA) {
    checkPairwiseAlignment(pA);

    Paf *paf = st_calloc(1, sizeof(Paf));

    // Set query coordinates
    paf->query_name = stString_copy(pA->contig1);
    paf->query_start = pA->strand1 ? pA->start1 : pA->end1+1;
    paf->query_end = pA->strand1 ? pA->end1 : pA->start1+1;

    // Set target coordinates
    paf->target_name = stString_copy(pA->contig2);
    paf->target_start = pA->strand2 ? pA->start2 : pA->end2+1;
    paf->target_end = pA->strand2 ? pA->end2 : pA->start2+1;

    // Set forward or reverse complement alignment
    paf->same_strand = pA->strand1 == pA->strand2;

    // Set score
    paf->score = pA->score;

    // Set cigar
    Cigar **pCigar = &(paf->cigar);
    for(int64_t i=0; i<pA->operationList->length; i++) {
        struct AlignmentOperation *op =
                pA->operationList->list[pA->strand1 ? i : pA->operationList->length-i-1];
        assert(op->length >= 1);
        Cigar *cigar = st_calloc(1, sizeof(Cigar));
        *pCigar = cigar;
        pCigar = &(cigar->next);
        cigar->length = op->length;
        if(op->opType == PAIRWISE_MATCH) {
            cigar->op = match;
        }
        else if(op->opType == PAIRWISE_INDEL_X) {
            cigar->op = query_insert;
        }
        else {
            assert(op->opType == PAIRWISE_INDEL_Y);
            cigar->op = query_delete;
        }
    }

    paf_check(paf);

    return paf;
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
    // Invert the paf
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

    struct PairwiseAlignment *pairwiseAlignment;
    while((pairwiseAlignment = cigarRead(input)) != NULL) {
        Paf *paf = cigar_to_paf(pairwiseAlignment);
        paf_check(paf);
        paf_write(paf, output);
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

    st_logInfo("Paf invert is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}
