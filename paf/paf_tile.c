/*
 * paf_tile: Greedily assign a "level" to each alignment, so that alignments with the lowest level are tesslating, high-scoring alignments.
 *
 * Released under the MIT license, see LICENSE.txt
 *
 * Rough outline:
 * (1) Load local alignments files (PAF)
 * (2) Sort alignments by score, from best-to-worst
 * (3) Create integer array representing counts of alignments to bases in the genome, setting values initially to 0.
 * (4) Iterate through alignments, from best-to-worst, find maximum?? count, q, of aligned bases to a base covered by the alignment,
 * set the "level" of the alignment to q+1, increase by one the aligned bases count of each base covered by the alignment.
 * (5) Output local alignments file, adding the alignment level tag to each alignment, sorted by score from best-to-worst
 */

#include "paf.h"
#include "inc/paf.h"
#include <getopt.h>
#include <time.h>

void usage() {
    fprintf(stderr, "paf_tile [options], version 0.1\n");
    fprintf(stderr, "Tiles the records in the PAF file\n");
    fprintf(stderr, "-i --inputFile : Input paf file. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int paf_cmp_by_descending_score(const void *a, const void *b) {
    Paf *p1 = (Paf *)a, *p2 = (Paf *)b;
    return p1->score > p2->score ? -1 : (p1->score < p2->score ? 1 : 0);
}

int64_t get_median_alignment_level(uint16_t *counts, Paf *paf) {
    Cigar *c = paf->cigar;
    int64_t i = paf->query_start, max_level=0, matches=0;
    int64_t *level_counts = st_calloc(UINT16_MAX, sizeof(int64_t)); // An array of counts of the number of bases with the given alignment level
    // such that level_counts[i] is the number of bases in the query with level_counts[i] number of alignments to it (at this point in the tiling)
    while(c != NULL) {
        if(c->op != query_delete) {
            if(c->op == match) {
                for(int64_t j=0; j<c->length; j++) {
                    assert(i + j < paf->query_end && i + j >= 0 && i + j < paf->query_length);
                    assert(counts[i + j] < UINT16_MAX); // paranoid check
                    level_counts[counts[i + j]]++;
                    matches++;
                    if(counts[i + j] > max_level) {
                        max_level = counts[i + j];
                    }
                }
            }
            i += c->length;
        }
        c = c->next;
    }
    assert(i == paf->query_end);

    if(matches == 0) { // avoid divide by zero
        free(level_counts);
        return INT16_MAX;
    }

    // Print the alignment levels
    if(st_getLogLevel() >= debug) {
        st_logDebug("Got alignment levels: ");
        for (i = 0; i <= max_level; i++) {
            st_logDebug("%"
            PRIi64
            ":%f ", i, ((float)level_counts[i])/matches);
        }
        char *paf_string = paf_print(paf);
        st_logDebug(" for paf: %s\n", paf_string);
        free(paf_string);
    }

    // Calc the median from the level_counts array
    int64_t j=0;
    for(i=0; i<=max_level; i++) {
        j += level_counts[i];
        if(j >= matches/2.0) {
            free(level_counts);
            assert(i > 0);
            return i;
        }
    }
    assert(0); // This should be unreachable.
    free(level_counts);
    return INT16_MAX;
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
    stList_sort(pafs, paf_cmp_by_descending_score); // Sort alignments by score, from best-to-worst

    // Create integer array representing counts of alignments to bases in the genome, setting values initially to 0.
    stHash *seq_names_to_alignment_count_arrays = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, free);

    // For each alignment: set the "level" of the alignment to q+1, increase by one the aligned bases count of each base covered by the alignment.
    for(int64_t i=0; i<stList_length(pafs); i++) {
        Paf *paf = stList_get(pafs, i);
        SequenceCountArray *seq_count_array = get_alignment_count_array(seq_names_to_alignment_count_arrays, paf);
        increase_alignment_level_counts(seq_count_array, paf);
        paf->tile_level = get_median_alignment_level(seq_count_array->counts, paf); // Store the tile_level
        assert(paf->tile_level > 0); // Tile levels should start at 1
    }

    // Output local alignments file, sorted by score from best-to-worst
    write_pafs(output, pafs);

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    stHash_destruct(seq_names_to_alignment_count_arrays);
    stList_destruct(pafs);
    if(inputFile != NULL) {
        fclose(input);
    }
    if(outputFile != NULL) {
        fclose(output);
    }

    st_logInfo("Paf tile is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}
