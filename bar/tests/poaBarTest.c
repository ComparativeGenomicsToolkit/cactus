#define _GNU_SOURCE /* mkstemp, under -std=c99 */
/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "flowersShared.h"
#include "randomSequences.h"
#include "poaBarAligner.h"
#ifdef HAVE_MINIPOA
#include "minipoa_c.h"
#endif
#include "stCaf.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>

//#define stderr_logging
/*
 * Every test below runs against both MSA engines.  The parameters mirror what the config would
 * produce, modulo the deliberately tiny band the tests use to keep the random cases quick.
 */
static const BaseAligner TEST_ENGINES[] = {
    BASE_ALIGNER_ABPOA,
#ifdef HAVE_MINIPOA
    BASE_ALIGNER_MINIPOA,
#endif
};
#define TEST_ENGINE_NO ((int64_t)(sizeof(TEST_ENGINES) / sizeof(TEST_ENGINES[0])))

static PoaParameters *test_poa_params(BaseAligner engine) {
    PoaParameters *poaParams = st_calloc(1, sizeof(PoaParameters));
    poaParams->engine = engine;
    if (engine == BASE_ALIGNER_ABPOA) {
        abpoa_para_t *abpt = abpoa_init_para();
        abpt->wb = 10;
        abpt->wf = 0.01;
        abpoa_post_set_para(abpt);
        poaParams->abpt = abpt;
    }
#ifdef HAVE_MINIPOA
    else {
        minipoa_para_t *mpt = minipoa_init_para();
        minipoa_set_band(mpt, 10, 0.01);
        minipoa_set_adaptive_band(mpt, 0);
        minipoa_set_seeding(mpt, 0, 19, 10, 0);
        minipoa_set_progressive(mpt, 0);
        poaParams->mpt = mpt;
    }
#endif
    return poaParams;
}

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
  for (int64_t engine_i = 0; engine_i < TEST_ENGINE_NO; engine_i++) {
    PoaParameters *poaParams = test_poa_params(TEST_ENGINES[engine_i]);
    for(int64_t test=0; test<100; test++) {
        for (int64_t poa_window_size = 5; poa_window_size < 120; poa_window_size += 15) {
#ifdef stderr_logging
            fprintf(stderr, "Running test_make_partial_order_alignment, test %i\n", (int)test);
#endif
            // parent string from which other strings are created

            char *parent_string = getRandomACGTSequence(st_randomInt(1, 100));
#ifdef stderr_logging            
            fprintf(stderr, "Parent string (length: %i): %s\n", (int)strlen(parent_string), parent_string);
#endif

            // get random string no
            int64_t seq_no = st_randomInt(1, 20);

            // get random strings
            char **seqs = st_malloc(sizeof(char *) * seq_no);
            int *seq_lens = st_malloc(sizeof(int) * seq_no);
            for(int64_t i=0; i<seq_no; i++) {
                seqs[i] = evolveSequence(parent_string);
                seq_lens[i] = strlen(seqs[i]);
#ifdef stderr_logging                
                fprintf(stderr, "String %i to align: %s (length: %i)\n", (int)i, seqs[i], (int)seq_lens[i]);
#endif
            }

            // generate the alignment
            Msa *msa = msa_make_partial_order_alignment(seqs, seq_lens, seq_no, poa_window_size, 1000, 0.02, poaParams);

            // print the msa
#ifdef stderr_logging
            msa_print(msa, stderr);
#endif

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
    poaParameters_destruct(poaParams);
  }
}

/**
 * Repeatedly generate random sets of two ends connected by set of strings, check that the resulting msa is valid
 */
void test_make_consistent_partial_order_alignments_two_ends(CuTest *testCase) {
  for (int64_t engine_i = 0; engine_i < TEST_ENGINE_NO; engine_i++) {
    PoaParameters *poaParams = test_poa_params(TEST_ENGINES[engine_i]);

    for(int64_t test=0; test<100; test++) {
#ifdef stderr_logging
        fprintf(stderr, "Running test_make_consistent_partial_order_alignments_two_ends, test %i\n", (int)test);
#endif

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
                                                              right_end_indexes, right_end_row_indexes, overlaps,
                                                              1000000, 100, 0.02, poaParams);

        // print the msas
#ifdef stderr_logging
        for(int64_t i=0; i<end_no; i++) {
            fprintf(stderr, "MSA: %i\n", (int)i);
            msa_print(msas[i], stderr);
        }
#endif

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
    poaParameters_destruct(poaParams);
  }
}

void test_make_flower_alignment_poa(CuTest *testCase) {
  for (int64_t engine_i = 0; engine_i < TEST_ENGINE_NO; engine_i++) {
    setup(testCase);
    PoaParameters *poaParams = test_poa_params(TEST_ENGINES[engine_i]);
#ifdef stderr_logging
    fprintf(stderr, "There are %i ends in the flower\n", (int)flower_getEndNumber(flower));
#endif
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
            int64_t length;
            char *s = get_adjacency_string(cap, &length, 1);
            Cap *adjacentCap = cap_getAdjacency(cap);
#ifdef stderr_logging            
            fprintf(stderr, "For end: %i, cap: %i (% " PRIi64 " to %" PRIi64 ") we have string: %s\n", (int)i, (int)j, cap_getName(cap), cap_getName(adjacentCap), s);
#endif
            j++;
        }
        end_destructInstanceIterator(capIterator);
        i++;
    }
    flower_destructEndIterator(endIterator);

    stList *alignment_blocks = make_flower_alignment_poa(flower, 2, 1000000, 5, 1000, 0.02, poaParams);

    for(int64_t i=0; i<stList_length(alignment_blocks); i++) {
        AlignmentBlock *b = stList_get(alignment_blocks, i);
#ifdef stderr_logging
        alignmentBlock_print(b, stderr);
#endif
    }

    poaParameters_destruct(poaParams);
    teardown(testCase);
  }
}


/*
 * Reference implementation of the mask-filter scan, run over the whole adjacency. This is what
 * get_adjacency_string_and_overlap used to do before it started fetching bounded windows.
 */
static int64_t ref_unmasked_length(const char *seq, int64_t seq_length, int64_t length, bool reversed, int64_t mask_filter) {
    if (mask_filter >= 0) {
        int64_t run_start = -1;
        for (int64_t i = 0; i < length; ++i) {
            char base = reversed ? seq[seq_length - 1 - i] : seq[i];
            if (islower(base) || base == 'N') {
                if (run_start == -1) {
                    run_start = i;
                }
                if (i + 1 - run_start > mask_filter) {
                    return run_start;
                }
            } else {
                run_start = -1;
            }
        }
    }
    return length;
}

/*
 * get_adjacency_string_and_overlap fetches only the prefix it keeps and the suffix the mask filter
 * scans, rather than materialising the whole adjacency. Check that agrees with slicing the whole
 * thing, at every truncation point and on both strands -- the reverse strand takes its window from
 * the far end, which is where this is easy to get wrong.
 */
void test_get_adjacency_string_and_overlap_bounded(CuTest *testCase) {
    setup(testCase);
    End *end;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        Cap *cap;
        End_InstanceIterator *capIterator = end_getInstanceIterator(end);
        while ((cap = end_getNext(capIterator)) != NULL) {
            if (cap_getSide(cap)) {
                cap = cap_getReverse(cap);
            }
            int64_t seq_length;
            char *whole = get_adjacency_string(cap, &seq_length, 1);
            CuAssertTrue(testCase, (int64_t) strlen(whole) == seq_length);

            for (int64_t max_len = 0; max_len <= seq_length + 1; max_len++) {
                for (int64_t mask_filter = -1; mask_filter <= 3; mask_filter++) {
                    int length;
                    int64_t overlap;
                    char *got = get_adjacency_string_and_overlap(cap, &length, &overlap, max_len, mask_filter);

                    int64_t forward = seq_length > max_len ? max_len : seq_length;
                    int64_t backward = forward;
                    if (mask_filter >= 0) {
                        forward = ref_unmasked_length(whole, seq_length, forward, false, mask_filter);
                        backward = ref_unmasked_length(whole, seq_length, forward, true, mask_filter);
                    }
                    int64_t expected_overlap = forward + backward > seq_length ? forward + backward - seq_length : 0;

                    CuAssertTrue(testCase, length == forward);
                    CuAssertTrue(testCase, overlap == expected_overlap);
                    CuAssertTrue(testCase, (int64_t) strlen(got) == forward);
                    CuAssertTrue(testCase, strncmp(got, whole, forward) == 0);
                    free(got);
                }
            }
            free(whole);
        }
        end_destructInstanceIterator(capIterator);
    }
    flower_destructEndIterator(endIterator);
    teardown(testCase);
}

void test_alignment_block_iterator(CuTest *testCase) {
  for (int64_t engine_i = 0; engine_i < TEST_ENGINE_NO; engine_i++) {
    setup(testCase);
    PoaParameters *poaParams = test_poa_params(TEST_ENGINES[engine_i]);

    stList *alignment_blocks = make_flower_alignment_poa(flower, 10000, 1000000, 5, 50, 0.05, poaParams);

    poaParameters_destruct(poaParams);
#ifdef stderr_logging
    for(int64_t i=0; i<stList_length(alignment_blocks); i++) {
        AlignmentBlock *b = stList_get(alignment_blocks, i);
        alignmentBlock_print(b, stderr);
    }
#endif

    stPinchIterator *it = stPinchIterator_constructFromAlignedBlocks(alignment_blocks);

    stPinchIterator_reset(it);
    //stPinchThreadSet *threadSet = stCaf_setup(flower);
    //stCaf_anneal(threadSet, it, NULL);

    stPinch *pinch, pinchToFillOut;
#ifdef stderr_logging
    while((pinch = stPinchIterator_getNext(it, &pinchToFillOut)) != NULL) {
        fprintf(stderr, "Pinch: name1: %" PRIi64 " s1:%i, name2: %" PRIi64 " s2:%i, length:%i, strand:%i\n",
                pinch->name1, (int)pinch->start1, pinch->name2, (int)pinch->start2,
                (int)pinch->length, (int)pinch->strand);
    }
#endif

    stPinchIterator_destruct(it);

    teardown(testCase);
  }
}

/*
 * <bar baseAligner> alone selects the engine; the <bar partialOrderAlignment> boolean it replaced
 * is not read at all, so a copy of an old config that still carries it must not change the result.
 */
static BaseAligner engine_for_config(CuTest *testCase, const char *bar_attrs) {
    const char *tmp = getenv("TMPDIR");
    char path[1024];
    snprintf(path, sizeof path, "%s/cactus_baseAlignerTestXXXXXX", (tmp && *tmp) ? tmp : "/tmp");
    int fd = mkstemp(path);
    CuAssertTrue(testCase, fd >= 0);
    FILE *f = fdopen(fd, "w");
    CuAssertPtrNotNull(testCase, f);
    fprintf(f, "<cactusWorkflowConfig><bar %s><poa/><minipoa/></bar></cactusWorkflowConfig>\n", bar_attrs);
    fclose(f);
    CactusParams *params = cactusParams_load(path);
    BaseAligner engine = baseAligner_constructFromCactusParams(params);
    cactusParams_destruct(params);
    remove(path);
    return engine;
}

void test_baseAligner_selection(CuTest *testCase) {
    // baseAligner decides, and all three values resolve
    CuAssertIntEquals(testCase, BASE_ALIGNER_PECAN, engine_for_config(testCase, "baseAligner=\"pecan\""));
    CuAssertIntEquals(testCase, BASE_ALIGNER_ABPOA, engine_for_config(testCase, "baseAligner=\"abpoa\""));
    CuAssertIntEquals(testCase, BASE_ALIGNER_MINIPOA, engine_for_config(testCase, "baseAligner=\"minipoa\""));

    // absent: abpoa
    CuAssertIntEquals(testCase, BASE_ALIGNER_ABPOA, engine_for_config(testCase, ""));

    // the old boolean is ignored, alone or alongside baseAligner
    CuAssertIntEquals(testCase, BASE_ALIGNER_ABPOA, engine_for_config(testCase, "partialOrderAlignment=\"0\""));
    CuAssertIntEquals(testCase, BASE_ALIGNER_MINIPOA,
                      engine_for_config(testCase, "partialOrderAlignment=\"0\" baseAligner=\"minipoa\""));
}

/*
 * Parse the config cactus actually ships and read every <bar> attribute the C side depends on.
 *
 * Nothing else does this: the build does not look at the XML, and the other tests here write
 * their own.  A stray "--" inside an XML comment got all the way past a clean build and a green
 * suite before an alignment run caught it, which is too late and too indirect.
 */
void test_shipped_config_is_loadable(CuTest *testCase) {
    const char *path = "src/cactus/cactus_progressive_config.xml";
    FILE *f = fopen(path, "r");
    if (f == NULL) {
        // run from somewhere other than the repo root; nothing to check
        return;
    }
    fclose(f);
    CactusParams *params = cactusParams_load((char *)path);
    CuAssertPtrNotNull(testCase, params);

    CuAssertIntEquals(testCase, BASE_ALIGNER_ABPOA, baseAligner_constructFromCactusParams(params));

    // every attribute poaBarAligner.c reads; the getters st_errAbort on a missing one
    PoaParameters *abpoa = poaParameters_constructFromCactusParams(params, BASE_ALIGNER_ABPOA);
    CuAssertPtrNotNull(testCase, abpoa);
    poaParameters_destruct(abpoa);
#ifdef HAVE_MINIPOA
    PoaParameters *minipoa = poaParameters_constructFromCactusParams(params, BASE_ALIGNER_MINIPOA);
    CuAssertPtrNotNull(testCase, minipoa);
    CuAssertTrue(testCase, minipoa->gapOpen > 0);
    CuAssertTrue(testCase, minipoa->gapExt > 0);
    CuAssertTrue(testCase, minipoa->mat[0] > 0);   // A/A must be a match
    CuAssertTrue(testCase, minipoa->mat[1] < 0);   // A/C must be a mismatch
    CuAssertTrue(testCase, minipoa->mat[24] > minipoa->mat[4]); // minipoa requires N/N > N/other
    poaParameters_destruct(minipoa);
#endif
    cactusParams_destruct(params);
}

CuSuite* poaBarAlignerTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_shipped_config_is_loadable);
    SUITE_ADD_TEST(suite, test_baseAligner_selection);
    SUITE_ADD_TEST(suite, test_make_partial_order_alignment);
    SUITE_ADD_TEST(suite, test_make_consistent_partial_order_alignments_two_ends);
    SUITE_ADD_TEST(suite, test_make_flower_alignment_poa);
    SUITE_ADD_TEST(suite, test_get_adjacency_string_and_overlap_bounded);
    SUITE_ADD_TEST(suite, test_alignment_block_iterator);
    return suite;
}
