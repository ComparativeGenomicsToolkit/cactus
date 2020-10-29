/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "flowersShared.h"
#include "randomSequences.h"
#include "poaBarAligner.h"
#include <stdio.h>
#include <ctype.h>

/**
 * Validate MSA. Lengths is an array that is populated with the lengths of the
 * sequences found on the MSA.
 */
void validate_msa(CuTest *testCase, Msa *msa, int64_t *lengths) {
    // Check that each sequence is correctly represented in the msa
    for(int64_t i=0; i<msa->seq_no; i++) {
        int64_t offset = 0;
        for(int64_t j=0; j<msa->column_no; j++) {
            char b = msa_to_base(msa->msa_seq[i][j]);
            if(b != '-') {
                CuAssertTrue(testCase, b == toupper(msa->seqs[i][offset++]));
            }
        }
        lengths[i] = offset;
    }
}

/**
 * Repeatedly generate random sets of closely related strings and test that returned msa is valid
 */
void test_make_partial_order_alignment(CuTest *testCase) {
    for(int64_t test=0; test<100; test++) {
        fprintf(stderr, "Running test_make_partial_order_alignment, test %i\n", (int)test);

        // parent string from which other strings are created

        char *parent_string = getRandomACGTSequence(st_randomInt(1, 100));
        fprintf(stderr, "Parent string (length: %i): %s\n", (int)strlen(parent_string), parent_string);

                // get random string no
        int64_t seq_no = st_randomInt(1, 20);

        // get random strings
        char **seqs = st_malloc(sizeof(char *) * seq_no);
        int *seq_lens = st_malloc(sizeof(int) * seq_no);
        for(int64_t i=0; i<seq_no; i++) {
            seqs[i] = evolveSequence(parent_string);
            seq_lens[i] = strlen(seqs[i]);
            fprintf(stderr, "String %i to align: %s (length: %i)\n", (int)i, seqs[i], (int)seq_lens[i]);
        }

        // generate the alignment
        Msa *msa = msa_make_partial_order_alignment(seqs, seq_lens, seq_no);

        // print the msa
        msa_print(msa, stderr);

        // validate the msa
        int64_t lengths[seq_no];
        validate_msa(testCase, msa, lengths);
        for(int64_t i=0; i<seq_no; i++) {
            CuAssertTrue(testCase, lengths[i] == seq_lens[i]);
        }

        // clean up
        msa_destruct(msa);
        free(parent_string);
    }
}

/**
 * Repeatedly generate random sets of two ends connected by set of strings, check that the resulting msa is valid
 */
void test_make_consistent_partial_order_alignments_two_ends(CuTest *testCase) {
    for(int64_t test=0; test<100; test++) {
        fprintf(stderr, "Running test_make_consistent_partial_order_alignments_two_ends, test %i\n", (int)test);

        // parent string from which other strings are created
        char *parent_string = getRandomACGTSequence(st_randomInt(1, 50));

        // get random string no
        int64_t seq_no = st_randomInt(1, 20);

        // build the two ends
        int64_t end_no = 2;
        int64_t end_lengths[end_no];
        char **end_strings[end_no];
        int *end_string_lengths[end_no];
        int64_t *right_end_indexes[end_no];
        int64_t *right_end_row_indexes[end_no];

        for(int64_t i=0; i<end_no; i++) {
            end_lengths[i] = seq_no;
            end_strings[i] = st_malloc(sizeof(char *) * seq_no);
            end_string_lengths[i] = st_malloc(sizeof(int) * seq_no);
            right_end_indexes[i] = st_malloc(sizeof(int64_t) * seq_no);
            right_end_row_indexes[i] = st_malloc(sizeof(int64_t) * seq_no);
        }

        int64_t j=st_randomInt(0, 1000);
        for(int64_t i=0; i<seq_no; i++) {
            char *c = evolveSequence(parent_string);
            int64_t k = (i + j)%seq_no; // The row index of the corresponding sequence for the second end

            end_strings[0][i] = c;
            end_string_lengths[0][i] = strlen(c);
            right_end_indexes[0][i] = 1;
            right_end_row_indexes[0][i] = k;

            end_strings[1][k] = stString_reverseComplementString(c);
            end_string_lengths[1][k] = strlen(c);
            right_end_indexes[1][k] = 0;
            right_end_row_indexes[1][k] = i;
        }

        // generate the alignments
        Msa **msas = make_consistent_partial_order_alignments(end_no, end_lengths, end_strings, end_string_lengths,
                                                              right_end_indexes, right_end_row_indexes);

        // print the msas
        for(int64_t i=0; i<end_no; i++) {
            fprintf(stderr, "MSA: %i\n", i);
            msa_print(msas[i], stderr);
        }

        // Validate the combination of both MSAs covers the complete sequences
        int64_t lengths1[seq_no], lengths2[seq_no];
        validate_msa(testCase, msas[0], lengths1);
        validate_msa(testCase, msas[1], lengths2);
        for(int64_t i=0; i<seq_no; i++) {
            //int64_t k = (i + j)%seq_no; // The row index of the corresponding sequence for the second end
            CuAssertTrue(testCase, lengths1[i] + lengths2[(i + j)%seq_no] == end_string_lengths[0][i]);
        }

        // clean up
        for(int64_t i=0; i<end_no; i++) {
            msa_destruct(msas[i]);
            free(right_end_indexes[i]);
            free(right_end_row_indexes[i]);
        }
        free(msas);
        free(parent_string);
    }
}

void test_make_flower_alignment_poa(CuTest *testCase) {
    setup(testCase);

    stList *alignment_blocks = make_flower_alignment_poa(flower, false);

    for(int64_t i=0; i<stList_length(alignment_blocks); i++) {
        AlignmentBlock *b = stList_get(alignment_blocks, i);
        //todo: complete
    }

    teardown(testCase);
}

CuSuite* poaBarAlignerTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_make_partial_order_alignment);
    SUITE_ADD_TEST(suite, test_make_consistent_partial_order_alignments_two_ends);
    SUITE_ADD_TEST(suite, test_make_flower_alignment_poa);
    return suite;
}
