#include "CuTest.h"
#include "taf.h"

static void test_maf(CuTest *testCase) {
    // Example maf file
    char *example_file = "./tests/evolverMammals.maf";
    char *temp_copy = "./tests/evolverMammals_copy.maf";

    //
    FILE *file = fopen(example_file, "r");
    FILE *out_file = fopen(temp_copy, "w");
    Alignment *alignment;
    while((alignment = maf_read_block(file)) != NULL) {
        //maf_write_block(alignment, stdout);
        maf_write_block(alignment, out_file);
        alignment_destruct(alignment);
    }
    fclose(file);
    fclose(out_file);

    // Now check that all the maf blocks are the same
    file = fopen(example_file, "r");
    FILE *file_copy = fopen(temp_copy, "r");
    Alignment *alignment2;
    while((alignment = maf_read_block(file)) != NULL) {
        alignment2 = maf_read_block(file_copy);
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
    CuAssertTrue(testCase, maf_read_block(file_copy) == NULL);

    fclose(file);
    fclose(out_file);
}

CuSuite* maf_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_maf);
    return suite;
}
