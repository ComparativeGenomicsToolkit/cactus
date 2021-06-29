/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "flowersShared.h"
#include "randomSequences.h"
#include "poaBarAligner.h"
#include "stCaf.h"
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
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->wb = 10;
    abpt->wf = 0.01;
    abpoa_post_set_para(abpt);
    for(int64_t test=0; test<100; test++) {
        for (int64_t poa_window_size = 5; poa_window_size < 120; poa_window_size += 15) {
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
            Msa *msa = msa_make_partial_order_alignment(seqs, seq_lens, seq_no, poa_window_size, abpt);

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
    abpoa_free_para(abpt);
}

/**
 * Repeatedly generate random sets of two ends connected by set of strings, check that the resulting msa is valid
 */
void test_make_consistent_partial_order_alignments_two_ends(CuTest *testCase) {
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->wb = 10;
    abpt->wf = 0.01;
    abpoa_post_set_para(abpt);
    
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
        int64_t *overlaps[end_no];

        for(int64_t i=0; i<end_no; i++) {
            end_lengths[i] = seq_no;
            end_strings[i] = st_malloc(sizeof(char *) * seq_no);
            end_string_lengths[i] = st_malloc(sizeof(int) * seq_no);
            right_end_indexes[i] = st_malloc(sizeof(int64_t) * seq_no);
            right_end_row_indexes[i] = st_malloc(sizeof(int64_t) * seq_no);
            overlaps[i] = st_malloc(sizeof(int64_t) * seq_no);
        }

        int64_t j=st_randomInt(0, 1000);
        for(int64_t i=0; i<seq_no; i++) {
            char *c = evolveSequence(parent_string);
            int64_t k = (i + j)%seq_no; // The row index of the corresponding sequence for the second end

            end_strings[0][i] = c;
            end_string_lengths[0][i] = strlen(c);
            right_end_indexes[0][i] = 1;
            right_end_row_indexes[0][i] = k;
            overlaps[0][i] = strlen(c);

            end_strings[1][k] = stString_reverseComplementString(c);
            end_string_lengths[1][k] = strlen(c);
            right_end_indexes[1][k] = 0;
            right_end_row_indexes[1][k] = i;
            overlaps[1][k] = strlen(c);
        }

        // generate the alignments
        Msa **msas = make_consistent_partial_order_alignments(end_no, end_lengths, end_strings, end_string_lengths,
                                                              right_end_indexes, right_end_row_indexes, overlaps, 1000000, abpt);

        // print the msas
        for(int64_t i=0; i<end_no; i++) {
            fprintf(stderr, "MSA: %i\n", (int)i);
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
            free(overlaps[i]);
        }
        free(msas);
        free(parent_string);
    }
    abpoa_free_para(abpt);
}

void test_make_flower_alignment_poa(CuTest *testCase) {
    setup(testCase);

    abpoa_para_t *abpt = abpoa_init_para();
    abpt->wb = 10;
    abpt->wf = 0.01;
    abpoa_post_set_para(abpt);

    fprintf(stderr, "There are %i ends in the flower\n", (int)flower_getEndNumber(flower));
    End *end;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    int64_t i=0; // Index of the end
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        // Now get each string incident with the end
        Cap *cap;
        End_InstanceIterator *capIterator = end_getInstanceIterator(end);
        int64_t j=0;
        while ((cap = end_getNext(capIterator)) != NULL) {
            if (cap_getSide(cap)) {
                cap = cap_getReverse(cap);
            }
            int length;
            char *s = get_adjacency_string(cap, &length, 1);
            Cap *adjacentCap = cap_getAdjacency(cap);
            fprintf(stderr, "For end: %i, cap: %i (% " PRIi64 " to %" PRIi64 ") we have string: %s\n", (int)i, (int)j, cap_getName(cap), cap_getName(adjacentCap), s);
            j++;
        }
        end_destructInstanceIterator(capIterator);
        i++;
    }
    flower_destructEndIterator(endIterator);

    stList *alignment_blocks = make_flower_alignment_poa(flower, 2, 1000000, 5, abpt);

    for(int64_t i=0; i<stList_length(alignment_blocks); i++) {
        AlignmentBlock *b = stList_get(alignment_blocks, i);
        alignmentBlock_print(b, stderr);
    }

    abpoa_free_para(abpt);
    teardown(testCase);
}

void test_alignment_block_iterator(CuTest *testCase) {
    setup(testCase);

    abpoa_para_t *abpt = abpoa_init_para();
    abpt->wb = 10;
    abpt->wf = 0.01;
    abpoa_post_set_para(abpt);

    stList *alignment_blocks = make_flower_alignment_poa(flower, 10000, 1000000, 5, abpt);

    abpoa_free_para(abpt);

    for(int64_t i=0; i<stList_length(alignment_blocks); i++) {
        AlignmentBlock *b = stList_get(alignment_blocks, i);
        alignmentBlock_print(b, stderr);
    }

    stPinchIterator *it = stPinchIterator_constructFromAlignedBlocks(alignment_blocks);

    stPinchIterator_reset(it);
    //stPinchThreadSet *threadSet = stCaf_setup(flower);
    //stCaf_anneal(threadSet, it, NULL);

    stPinch *pinch, pinchToFillOut;
    while((pinch = stPinchIterator_getNext(it, &pinchToFillOut)) != NULL) {
        fprintf(stderr, "Pinch: name1: %" PRIi64 " s1:%i, name2: %" PRIi64 " s2:%i, length:%i, strand:%i\n",
                pinch->name1, (int)pinch->start1, pinch->name2, (int)pinch->start2,
                (int)pinch->length, (int)pinch->strand);
    }

    stPinchIterator_destruct(it);

    teardown(testCase);
}

CuSuite* poaBarAlignerTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_make_partial_order_alignment);
    SUITE_ADD_TEST(suite, test_make_consistent_partial_order_alignments_two_ends);
    SUITE_ADD_TEST(suite, test_make_flower_alignment_poa);
    SUITE_ADD_TEST(suite, test_alignment_block_iterator);
    return suite;
}
