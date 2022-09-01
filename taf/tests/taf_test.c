#include "CuTest.h"
#include "taf.h"

static void test_taf(CuTest *testCase) {
    // Example maf file
    char *example_file = "./tests/chr2_KI270776v1_alt.maf"; //"./tests/chr2_KI270893v1_alt.maf"; //evolverMammals.maf";
    char *temp_copy = "./tests/chr2_KI270776v1_alt.taf.rle"; //evolverMammals.taf";
    bool run_length_encode_bases = 1;

    // Write out the taf file
    FILE *file = fopen(example_file, "r");
    FILE *out_file = fopen(temp_copy, "w");
    Alignment *alignment, *p_alignment = NULL;
    while((alignment = maf_read_block(file)) != NULL) {
        if(p_alignment != NULL) {
            alignment_link_adjacent(p_alignment, alignment);
        }
        taf_write_block(p_alignment, alignment, run_length_encode_bases, out_file);
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment);
        }
        p_alignment = alignment;
    }
    if(p_alignment != NULL) {
        alignment_destruct(p_alignment);
    }
    fclose(file);
    fclose(out_file);

    // Now parse the taf
    file = fopen(example_file, "r");
    FILE *file_copy = fopen(temp_copy, "r");
    LI *li = LI_construct(file_copy);
    Alignment *alignment2, *p_alignment2;
    while((alignment = maf_read_block(file)) != NULL) {
        // Alignment *taf_read_block(Alignment *p_block, bool run_length_encode_bases, LI *li)
        alignment2 = taf_read_block(p_alignment2, run_length_encode_bases, li);
        CuAssertTrue(testCase, alignment2 != NULL);
        // Check that the blocks are the same

        Alignment_Row *row = alignment->row, *row2 = alignment2->row;
        while(row != NULL) {
            CuAssertTrue(testCase, row2 != NULL);
            CuAssertStrEquals(testCase, row->sequence_name, row2->sequence_name);
            CuAssertIntEquals(testCase, row->start, row2->start);
            CuAssertIntEquals(testCase, row->length, row2->length);
            CuAssertIntEquals(testCase, row->sequence_length, row2->sequence_length);
            CuAssertIntEquals(testCase, row->strand, row2->strand);
            CuAssertStrEquals(testCase, row->bases, row2->bases);
            row = row->n_row; row2 = row2->n_row;
        }
        CuAssertTrue(testCase, row2 == NULL);

        alignment_destruct(alignment);
        if(p_alignment2 != NULL) {
            alignment_destruct(p_alignment2);
        }
        p_alignment2 = alignment2;
    }
    CuAssertTrue(testCase, taf_read_block(p_alignment2, run_length_encode_bases, li) == NULL);

    LI_destruct(li);
    fclose(file);
    fclose(out_file);
}

CuSuite* taf_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_taf);
    return suite;
}
