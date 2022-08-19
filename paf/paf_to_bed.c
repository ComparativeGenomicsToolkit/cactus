/*
 * paf_to_bed: Creates a bed file representing the coverages of paf alignments.
 *
 * Released under the MIT license, see LICENSE.txt
 *
 * Rough outline:
 * (1) Load local alignments files (PAF)
 * (3) Create integer array representing counts of alignments to bases in the genome, setting values initially to 0.
 * (4) Iterate through alignments, increase by one the aligned bases count of each base covered by the alignment.
 * (5) Output bed file, representing the coverage of each base in the query sequences
 */

#include "paf.h"
#include "inc/paf.h"
#include <getopt.h>
#include <time.h>

void usage() {
    fprintf(stderr, "paf_to_bed [options], version 0.1\n");
    fprintf(stderr, "Creates a bed file representing the coverage of alignments on the query sequences of the paf alignments\n");
    fprintf(stderr, "-i --inputFile : Input paf file. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output bed file. If not specified outputs to stdout\n");
    fprintf(stderr, "-b --binary : Output 0 for unaligned and 1 for aligned, rather than reporting the actual number of alignments\n");
    fprintf(stderr, "-e --excludeUnaligned : Exclude any interval with 0 alignment coverage from the output\n");
    fprintf(stderr, "-f --excludeAligned : Exclude any interval with > 0 alignment coverage from the output\n");
    fprintf(stderr, "-m --minSize : Exclude any interval shorter than this length from the output\n");
    fprintf(stderr, "-n --includeInverted : Flip the alignments to include the target sequences also.\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

void write_bed(FILE *output, stHash *seq_names_to_alignment_count_arrays,
               bool binary, bool exclude_unaligned, bool exclude_aligned, int64_t min_size) {
    stHashIterator *it = stHash_getIterator(seq_names_to_alignment_count_arrays);
    char *seq_name;
    while((seq_name = stHash_getNext(it)) != NULL) {
        SequenceCountArray *seq_count_array = stHash_search(seq_names_to_alignment_count_arrays, seq_name);
        for(int64_t i=0; i<seq_count_array->length;) {
            for(int64_t j=i+1; j<=seq_count_array->length; j++) {
                if(j == seq_count_array->length ||
                        (binary ? (seq_count_array->counts[i] > 0) != (seq_count_array->counts[j] > 0)
                                : seq_count_array->counts[i] != seq_count_array->counts[j])) {
                    if(j - i >= min_size && (seq_count_array->counts[i] == 0 ? !exclude_unaligned : !exclude_aligned)) {
                        fprintf(output, "%s %" PRIi64 " %" PRIi64 " %i\n",
                                seq_name, i, j, binary ? seq_count_array->counts[i] > 0 : seq_count_array->counts[i]);
                    }
                    i = j;
                    break;
                }
            }
        }
    }
    stHash_destructIterator(it);
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    bool binary = 0;
    bool exclude_unaligned = 0;
    bool exclude_aligned = 0;
    int64_t min_size = 1;
    bool include_inverted_alignments = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "binary", no_argument, 0, 'b' },
                                                { "excludeUnaligned", no_argument, 0, 'e' },
                                                { "excludeAligned", no_argument, 0, 'f' },
                                                { "minSize", required_argument, 0, 'm' },
                                                { "includeInverted", no_argument, 0, 'n' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:hbefm:n", long_options, &option_index);
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
            case 'b':
                binary = 1;
                break;
            case 'e':
                exclude_unaligned = 1;
                break;
            case 'f':
                exclude_aligned = 1;
                break;
            case 'm':
                min_size = atol(optarg);
                break;
            case 'n':
                include_inverted_alignments = 1;
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
    // Calculate the paf coverages
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

    // Create integer array representing counts of alignments to bases in the genome, setting values initially to 0.
    stHash *seq_names_to_alignment_count_arrays = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, free);

    // For each alignment: increase by one the aligned bases count of each base covered by the alignment.
    Paf *paf;
    while((paf = paf_read(input)) != NULL) {
        SequenceCountArray *seq_count_array = get_alignment_count_array(seq_names_to_alignment_count_arrays, paf);
        increase_alignment_level_counts(seq_count_array, paf);

        if(include_inverted_alignments) {
            paf_invert(paf); // Flip the alignment
            seq_count_array = get_alignment_count_array(seq_names_to_alignment_count_arrays, paf);
            increase_alignment_level_counts(seq_count_array, paf);
        }
    }

    // Output local alignments file, sorted by score from best-to-worst
    write_bed(output, seq_names_to_alignment_count_arrays, binary, exclude_unaligned, exclude_aligned, min_size);

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    stHash_destruct(seq_names_to_alignment_count_arrays);
    if(inputFile != NULL) {
        fclose(input);
    }
    if(outputFile != NULL) {
        fclose(output);
    }

    st_logInfo("Paf to bed is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}
