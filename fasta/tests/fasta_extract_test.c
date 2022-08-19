#include "paf.h"
#include "CuTest.h"
#include "sonLib.h"
#include "bioioC.h"

/*
 * Test fasta extract
 */

static char *test_fa_file = "./fasta/tests/temp.fa";
static char *test_bed_file = "./fasta/tests/temp.bed";
static char *test_out_file = "./fasta/tests/out.fa";

static void test_fasta_extract(CuTest *testCase) {
    for(int64_t test=0; test<1000; test++) {
        // Write some random sequences to the file
        FILE *fh = fopen(test_fa_file, "w");
        int64_t seq_no = st_randomInt64(1, 10);
        stHash *dna_strings = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
        stHash *dna_strings_with_xs = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, free);
        for(int64_t i=0; i<seq_no; i++) {
            char *r_string = stRandom_getRandomDNAString(st_randomInt64(0, 1000), 0, 0, 0);
            char *header = stString_print("%i", i);
            stHash_insert(dna_strings, header, r_string);
            stHash_insert(dna_strings_with_xs, header, stString_copy(r_string));
            fastaWrite(r_string, header, fh);
        }
        fclose(fh);

        // Choose a random flank length
        int64_t flank = st_randomInt64(0, 10);
        int64_t min_size = st_randomInt64(0, 10);

        // Choose some random intervals and mark them
        fh = fopen(test_bed_file, "w");
        int64_t interval_no = st_randomInt64(0, 100);
        int64_t total_bases = 0;
        for(int64_t i=0; i<interval_no; i++) {
            char *header = stString_print("%i", st_randomInt64(0, seq_no));
            char *r_string = stHash_search(dna_strings_with_xs, header);
            if(strlen(r_string) > 0) {
                int64_t start = st_randomInt64(0, strlen(r_string));
                int64_t end = st_randomInt64(start, strlen(r_string) + 1);
                assert(0 <= start);
                assert(start <= end);
                assert(start < strlen(r_string));
                assert(end <= strlen(r_string));

                // Print the interval to the bed file
                fprintf(fh, "%s %" PRIi64 " %" PRIi64 "\n", header, start, end);

                if(end-start >= min_size) { // If the interval is larger than the minimum size
                    // Mark the bases in the interval and count any new bases
                    for (int64_t j = start - flank; j < end + flank; j++) {
                        if (j >= 0 && j < strlen(r_string)) {
                            if (r_string[j] != 'X') {
                                total_bases++;
                                r_string[j] = 'X';
                            }
                        }
                    }
                }
            }
        }
        fclose(fh);

        st_logDebug("Got fasta_extract_test have %" PRIi64 " sequences, %" PRIi64 " intervals and %" PRIi64 " marked bases\n",
                    seq_no, interval_no, total_bases);

        // Run fasta extract to get subsequences
        st_system("fasta_extract %s -i %s -o %s -f %" PRIi64 " --minSize %" PRIi64 "", test_fa_file, test_bed_file, test_out_file, flank, min_size);

        // Parse extracted sequences and check that the sequences
        // correspond to the correct intervals
        fh = fopen(test_out_file, "r");
        stHash *sub_seqs = fastaReadToMap(fh);
        stHashIterator *it = stHash_getIterator(sub_seqs);
        char *sub_seq_header;
        int64_t total_bases_in_output = 0;
        while((sub_seq_header = stHash_getNext(it)) != NULL) {
            stList *tokens = stString_splitByString(sub_seq_header, "|");
            char *seq_name = stList_get(tokens, 0);
            int64_t sub_seq_start = atol(stList_get(tokens, 2));
            char *seq = stHash_search(dna_strings, seq_name);
            char *sub_seq = stHash_search(sub_seqs, sub_seq_header);

            for(int64_t i=0; i<strlen(sub_seq); i++) {
                CuAssertIntEquals(testCase, sub_seq[i], seq[i + sub_seq_start]); // Check is sequence is as expected
                seq[i + sub_seq_start] = 'X'; // Mark the base so if it occurs in another interval we will know
            }
            total_bases_in_output += strlen(sub_seq);
        }
        stHash_destructIterator(it);
        fclose(fh);

        // Check we saw all the bases
        CuAssertIntEquals(testCase, total_bases, total_bases_in_output);

        // Clean up
        stHash_destruct(sub_seqs);
        stHash_destruct(dna_strings);
        stHash_destruct(dna_strings_with_xs);
        st_system("rm %s %s %s", test_fa_file, test_bed_file, test_out_file);
    }
}

CuSuite* addFastaExtractTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_fasta_extract);
    return suite;
}
