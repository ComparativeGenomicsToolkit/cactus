#include "CuTest.h"
#include "taf.h"

static void test_taf(CuTest *testCase) {
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf";
    char *temp_copy = "./tests/evolverMammals.taf";

    //
    FILE *file = fopen(example_file, "r");
    FILE *out_file = fopen(temp_copy, "w");
    Alignment *alignment, *p_alignment = NULL;
    while((alignment = maf_read_block(file)) != NULL) {
        if(p_alignment != NULL) {
            alignment_link_adjacent(p_alignment, alignment);
        }
        //maf_write_block(alignment, stdout);
        taf_write_block(alignment, out_file);
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

    // Now check that all the maf blocks are the same
    /*file = fopen(example_file, "r");
    FILE *file_copy = fopen(temp_copy, "r");
    Alignment *alignment2;
    while((alignment = maf_read_block(file)) != NULL) {
        alignment2 = taf_read_column(file_copy);
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
        alignment_destruct(alignment2);
    }
    CuAssertTrue(testCase, taf_read_block(file_copy) == NULL);

    fclose(file);
    fclose(out_file);*/
}

CuSuite* taf_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_taf);
    return suite;
}
