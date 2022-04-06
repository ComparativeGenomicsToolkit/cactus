/*
 * paf_consistency: Reads in a set of primary and secondary alignments. Determines an approximate "consistency" score for each
 * secondary alignment based upon how much it is transitively implied by the primary alignments. Secondary
 * alignments with high score are then converted to primary alignments. This is useful for filling in alignments missed
 * by the primary alignments.
 *
 * Released under the MIT license, see LICENSE.txt
 *
 * Rough outline:
 * (1) Load primary alignments (PAF) both forward and inverted.
 * (2) Load the secondary alignments.
 * (3) Iterate through the primary alignments calculating the transitive alignments. For each transitive alignment
 * find any secondary that overlaps this alignment and add to its consistency score.
 * (4) For each secondary alignment determine if it should be output as a primary alignment based upon its consistency
 * score.
 */

#include "paf.h"
#include "inc/paf.h"
#include <getopt.h>
#include <time.h>

void usage() {
    fprintf(stderr, "paf_consistency [PRIMARY_INPUT_PAF_FILE] [options], version 0.1\n");
    fprintf(stderr, "Rescores secondary alignments by consistency with primary alignments\n");
    fprintf(stderr, "-i --secondaryInputFile : Input paf file containing secondary alignments. If not specified reads from stdin\n");
    fprintf(stderr, "-o --secondaryOutputFile : Output paf file to write the secondary alignments in. If not specified outputs to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

static int paf_query_cmp(const void *k, const void *k2) {
    /*
     * Compare pafs by query sequence coordinates
     */
    Paf *p = (Paf *)k, *p2 = (Paf *)k2;
    int i = strcmp(p->query_name, p2->query_name);
    if(i == 0) { // On same query
        i = p->query_start < p2->query_start ? -1 : (p->query_start > p2->query_start ? 1 : 0);
        if(i == 0) { // Have same start coordinate
            i = p->query_end < p2->query_end ? -1 : (p->query_end > p2->query_end ? 1 : 0); // Compare by end coordinate
        }
    }
    return i;
}

typedef struct _sequenceInterval {
    /*
     * Struct to represent a coordinate interval on a sequences.
     */
    char *seq_name;
    int64_t start, end;
    bool strand;
} SequenceInterval;

SequenceInterval get_query_interval(Paf *paf1) {
    /*
     * Get the query interval of the paf as a sequence interval.
     */
    SequenceInterval s;
    s.seq_name = paf1->query_name; // We don't copy this string
    s.start = paf1->query_start;
    s.end = paf1->query_end;
    s.strand = 1; // Assume all query sequences are on the positive strand, only target sequences can flip strand
    // according to alignments
    return s;
}

SequenceInterval overlap(SequenceInterval s, SequenceInterval s2) {
    /*
     * Returns a sequence representing the overlap (if any) between two sequence intervals.
     * To overlap must be on the same sequence and strand.
     * If no overlap seq_name will be initialized to NULL and other fields are undefined.
     */
    SequenceInterval s3;
    s3.seq_name = NULL;
    // Are on the same non-null sequence and strand
    if(s.seq_name != NULL && s2.seq_name != NULL && strcmp(s.seq_name, s2.seq_name) == 0 && s.strand == s2.strand) {
        s3.seq_name = s.seq_name;
        // the following logic will define a valid, non-zero length interval only if start < end
        s3.start = s.start <= s2.start ? s2.start : s.start; // Take the largest start coordinate
        s3.end = s.end <= s2.end ? s.end : s2.end; // Take the smallest end coordinate
        s3.strand = s.strand; // Set the strand to match the inputs
    }
    return s3;
}

bool interval_has_non_zero_length(SequenceInterval s) {
    return s.seq_name != NULL && s.start < s.end;
}

SequenceInterval get_target_interval(Paf *paf, SequenceInterval s) {
    s = overlap(get_query_interval(paf), s); // Make sure that the interval is within the paf query coordinates
    SequenceInterval s2;
    if(interval_has_non_zero_length(s) && paf->target_end > paf->target_start) { // The interval overlaps the paf query coordinates
        // and the target interval has non-zero length
        s2.seq_name = paf->target_name; // Set the sequence name
        s2.strand = paf->same_strand; // Set the strand of the target interval,
        // this will be reversed on the target if the same_strand field is false.

        double f = (float)(paf->target_end - paf->target_start) / (float)(paf->query_end - paf->query_start); // Simple scaling factor
        // to correct for length differences in the query and target intervals
        assert(f > 0);

        // Calculate the overlapping interval target coordinates in a very approximate way
        if(paf->same_strand) {
            s2.start = paf->target_start + (s.start - paf->query_start) * f;
            s2.end = paf->target_start + (s.end - paf->query_start) * f;
            assert(paf->target_start <= s2.start <= s2.end <= paf->target_end);
        }
        else {
            // End and start reverse as we reverse the alignment
            s2.end = paf->target_end - (s.start - paf->query_start) * f;
            s2.start = paf->target_end - (s.end - paf->query_start) * f;
            assert(paf->target_start <= s2.start <= s2.end <= paf->target_end);
        }
    }
    else {
        s2.seq_name = NULL; // Indicate that the lifted interval has zero length
    }
    return s2;
}

static int cmp_paf_to_interval(const void *a, const void *b) {
    SequenceInterval *s = (SequenceInterval *)a;
    Paf *paf = (Paf *)b;
    int i = strcmp(s->seq_name, paf->query_name);
    if(i == 0) {
        i = s->start < paf->query_start ? -1 : (s->start > paf->query_start ? 1 : 0);
        if(i == 0) {
            i = s->end < paf->query_end ? -1 : (s->end > paf->query_end ? 1 : 0);
        }
    }
    return i;
}

int64_t find_index_of_first_overlapping_alignment(stList *paf_alignments, SequenceInterval s) {
    /*
     * Returns index of the first paf alignment overlapping this interval
     */
    return stList_binarySearchFirstIndex(paf_alignments, &s, cmp_paf_to_interval);
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *secondary_input_file = NULL;
    char *secondary_output_file = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "secondaryInputFile", required_argument, 0, 'i' },
                                                { "secondaryOutputFile", required_argument, 0, 'o' },
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
                secondary_input_file = optarg;
                break;
            case 'o':
                secondary_output_file = optarg;
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
        fprintf(stderr, "Expected a primary pafs file\n");
        exit(1);
    }

    char *primary_input_file = argv[optind++]; // The paf file of the primary alignments

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Primary alignments input file string : %s\n", primary_input_file);
    st_logInfo("Secondary alignments input file string : %s\n", secondary_input_file);
    st_logInfo("Secondary alignments output file string : %s\n", secondary_output_file);

    //////////////////////////////////////////////
    // Read in the primary PAFs
    //////////////////////////////////////////////

    FILE *input = fopen(primary_input_file, "r");
    stList *primary_pafs = read_pafs(input); // Load the primary alignments
    fclose(input); input = fopen(primary_input_file, "r");
    stList *p = read_pafs(input); // Load the primary alignments again
    fclose(input);
    while(stList_length(p) > 0) { // Invert all the alignments in p and append them to primary_pafs
        Paf *paf = stList_pop(p);
        paf_invert(paf);
        stList_append(primary_pafs, paf);
    }
    stList_destruct(p); // Cleanup
    stList_sort(primary_pafs, paf_query_cmp); // Sort by query coordinate
    // Now we have all the primary pafs forward and inverted in primary_pafs
    st_logInfo("Read the primary pafs, there are %i total primary pafs \n", (int)stList_length(primary_pafs)/2);

    //////////////////////////////////////////////
    // Load the secondary alignments
    //////////////////////////////////////////////

    input = secondary_input_file == NULL ? stdin : fopen(secondary_input_file, "r");
    stList *secondary_pafs = read_pafs(input); // Load local alignments files (PAF)
    if(secondary_input_file != NULL) {
        fclose(input);
    }
    st_logInfo("Read the secondary pafs, there are %i total secondary pafs \n", (int)stList_length(secondary_pafs));

    //////////////////////////////////////////////
    // Calculate the consistency scores
    //////////////////////////////////////////////

    // An array of alignment scores for the secondary alignments, each score approximates how many bases of
    // the secondary alignment are implied transitively by the primary alignments
    int64_t *secondary_alignment_scores = st_calloc(stList_length(secondary_pafs), sizeof(int64_t));

    // For each possible pair of overlapping alignments, sorted by query sequence coordinates
    for(int64_t i=0; i<stList_length(primary_pafs); i++) {
        Paf *paf1 = stList_get(primary_pafs, i);
        for(int64_t j=i+1; j<stList_length(primary_pafs); j++) {
            Paf *paf2 = stList_get(primary_pafs, j);

            // Get any overlap between the sequences
            SequenceInterval s = overlap(get_query_interval(paf1), get_query_interval(paf2));
            if(interval_has_non_zero_length(s)) { // If the pafs overlap on the query sequence
                // Get the approximate intervals on the target sequences of the overlapping query interval
                // - these projected intervals form an approximate transitive alignment of the target sequences through
                // the shared, overlapping query sequence
                SequenceInterval interval1 = get_target_interval(paf1, s);
                SequenceInterval interval2 = get_target_interval(paf2, s);

                // For each secondary alignment (sorted by query sequence coordinates) that overlaps interval1 on the query sequence
                int64_t k = 0; //find_index_of_first_overlapping_alignment(secondary_pafs, interval1);
                while(k >= 0 && k < stList_length(secondary_pafs)) {
                    Paf *paf3 = stList_get(secondary_pafs, k);
                    SequenceInterval interval3 = overlap(get_query_interval(paf3), interval1); // Approximate interval on query
                    if(interval_has_non_zero_length(interval3)) {

                        // Project interval3 through the secondary paf to calculate the approximate interval on the target sequence
                        SequenceInterval interval4 = get_target_interval(paf3, interval3);

                        // If interval2 and interval4 overlap then we have an agreement in the approximate alignment
                        SequenceInterval interval5 = overlap(interval2, interval4);
                        if(interval_has_non_zero_length(interval5)) {
                            assert(interval5.end > interval5.start);
                            secondary_alignment_scores[k] += interval5.end - interval5.start; // Add the length of the overlap to the score
                        }
                    }
                    else { // No more overlap between a sequence and interval1 is possible on the query sequence
                        //break;
                    }
                    k++;
                }
            }
            else {
                break; // Can be no more overlaps with paf1
            }
        }
    }

    // Adjust the tile-level of secondary alignments if they have significant overlap with the transitively implied primary alignments
    for(int64_t i=0; i<stList_length(secondary_pafs); i++) {
        Paf *paf = stList_get(secondary_pafs, i);
        if(secondary_alignment_scores[i] >= 1) {
            paf->tile_level = 1;
            paf->type = 'P';
            //st_uglyf("%i %s %i\n", (int)secondary_alignment_scores[i], paf->query_name, (int)paf->query_end - paf->query_start);
        }
    }

    // Now write out the secondary PAFs
    FILE *output = secondary_output_file == NULL ? stdout : fopen(secondary_output_file, "w");
    for(int64_t i=0; i<stList_length(secondary_pafs); i++) {
        Paf *paf = stList_get(secondary_pafs, i);
        if(paf->tile_level == 1) {
            paf_write(paf, output);
        }
    }
    //write_pafs(output, secondary_pafs);
    if(secondary_output_file != NULL) {
        fclose(output);
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    free(secondary_alignment_scores);
    stList_destruct(primary_pafs);
    stList_destruct(secondary_pafs);
    st_logInfo("Paf consistency is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}
