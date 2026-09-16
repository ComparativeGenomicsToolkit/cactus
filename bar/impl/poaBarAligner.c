/**
 * This is designed as a drop-in replacement for the bar aligner, using the abpoa multiple sequence aligner.
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <limits.h>
#include "abpoa.h"
#ifdef HAVE_MINIPOA
#include "minipoa_c.h"
#endif
#include "poaBarAligner.h"
#include "flowerAligner.h"

#include <stdio.h>
#include <ctype.h>
#include <unistd.h>

/* ============================================================================ */
/* Shared: the CACTUS_BAR_DUMP_DIR switch                                       */
/* ============================================================================ */

/*
 * Set CACTUS_BAR_DUMP_DIR in the environment and every window handed to the base aligner is
 * written there as a FASTA, alongside the substitution matrix and a command line that replays it.
 * That is how a suspect alignment gets out of a 64-thread run and onto a developer's laptop.
 *
 * This used to be a commented-out #define that also deleted its own output on success, so using
 * it meant editing and rebuilding, and the files were gone by the time you looked.  Read once in
 * bar(), before the flower loop, so the OpenMP region never calls getenv.
 */
static const char *bar_dump_dir = NULL;
static int64_t bar_dump_counter = 0;

void bar_dump_dir_init(void) {
    bar_dump_dir = getenv("CACTUS_BAR_DUMP_DIR");
    if (bar_dump_dir != NULL && bar_dump_dir[0] == '\0') {
        bar_dump_dir = NULL; // set-but-empty means off, not a dump into the current directory
    }
    if (bar_dump_dir != NULL) {
        st_logInfo("bar: dumping base-aligner windows to %s\n", bar_dump_dir);
    }
}

// FOR DEBUGGING ONLY: run abpoa from the command line instead of via the API.
// Requires CACTUS_BAR_DUMP_DIR to be set at run time.
//#define CACTUS_ABPOA_FROM_COMMAND_LINE

// OpenMP
#if defined(_OPENMP)
#include <omp.h>
#endif


/* ============================================================================ */
/* Shared: alphabet, and the Msa the backends fill in                           */
/* ============================================================================ */


// char <--> uint8_t conversion copied over from abPOA example
// AaCcGgTtNn ==> 0,1,2,3,4
static unsigned char nst_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
static const char nst_nt256_table[256] = {
       'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};

char msa_to_base(uint8_t n) {
    return (char)nst_nt256_table[n];
}

uint8_t msa_to_byte(char c) {
    return nst_nt4_table[(int)c];
}

static uint8_t rc_table[6] = { 3, 2, 1, 0, 4, 5 };

static inline uint8_t msa_to_rc(uint8_t n) {
    return rc_table[n];
}





void msa_destruct(Msa *msa) {
    for(int64_t i=0; i<msa->seq_no; i++) {
        if (msa->seqs != NULL) {
            free(msa->seqs[i]);
        }
        free(msa->msa_seq[i]);
    }
    free(msa->seqs);
    free(msa->msa_seq);
    free(msa->seq_lens);
    free(msa);
}

void msa_print(Msa *msa, FILE *f) {
    fprintf(f, "MSA. Seq no: %i column no: %i \n", (int)msa->seq_no, (int)msa->column_no);
    for(int64_t i=0; i<msa->seq_no; i++) {
        fprintf(f, "Row:%i [len=%i]\t", (int)i, (int)msa->seq_lens[i]);
        for(int64_t j=0; j<msa->column_no; j++) {
            fprintf(f, "%c", msa_to_base(msa->msa_seq[i][j]));
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}

/* ============================================================================ */
/* Shared: column scores, trimming and window stitching                         */
/* ============================================================================ */


/**
 * flip msa to its reverse complement (for trimming purposees)
 */
static void flip_msa_seq(Msa* msa) {
    if (msa != NULL) {
        int64_t middle = msa->column_no / 2;
        bool odd = msa->column_no % 2 == 1;
        for (int64_t i = 0; i < msa->seq_no; ++i) {
            for (int64_t j = 0; j < middle; ++j) {
                uint8_t buf = msa->msa_seq[i][j];
                msa->msa_seq[i][j] = msa_to_rc(msa->msa_seq[i][msa->column_no - 1 - j]);
                msa->msa_seq[i][msa->column_no - 1 - j] = msa_to_rc(buf);
            }
            if (odd) {
                msa->msa_seq[i][middle] = msa_to_rc(msa->msa_seq[i][middle]);
            }
        }
    }
}

/**
 * Returns an array of floats, one for each corresponding column in the MSA. Each float
 * is the score of the column in the alignment.
 */
static float *make_column_scores(Msa *msa) {
    float *column_scores = st_calloc(msa->column_no, sizeof(float));
    for(int64_t i=0; i<msa->column_no; i++) {
        // Score is simply max(number of aligned bases in the column - 1, 0)
        for(int64_t j=0; j<msa->seq_no; j++) {
            if(msa_to_base(msa->msa_seq[j][i]) != '-') {
                column_scores[i]++;
            }
        }
        if(column_scores[i] >= 1.0) {
            column_scores[i]--;
        }
        assert(column_scores[i] >= 0.0);
    }
    return column_scores;
}

/**
 * Fills in cu_column_scores with the cumulative sum of column scores, from left-to-right, of columns
 * containing a non-gap character in the given "row".
 */
static void sum_column_scores(int64_t row, Msa *msa, float *column_scores, float *cu_column_scores) {
    float cu_score = 0.0; // The cumulative sum of column scores containing bases for the given row
    int64_t j=0; // The index in the DNA string for the given row
    for(int64_t i=0; i<msa->column_no; i++) {
        if(msa_to_base(msa->msa_seq[row][i]) != '-') {
            cu_score += column_scores[i];
            cu_column_scores[j++] = cu_score;
        }
    }
    assert(msa->seq_lens[row] == j); // We should cover all the bases in the DNA sequence
}

/**
 * Removes the suffix of the given row from the MSA and updates the column scores. suffix_start is the beginning
 * suffix to remove.
 */
static void trim_msa_suffix(Msa *msa, float *column_scores, int64_t row, int64_t suffix_start) {
    int64_t seq_index = 0;
    for(int64_t i=0; i<msa->column_no; i++) {
        if(msa_to_base(msa->msa_seq[row][i]) != '-') {
            if(seq_index++ >= suffix_start) {
                msa->msa_seq[row][i] = msa_to_byte('-');
                column_scores[i] = column_scores[i] > 1 ? column_scores[i]-1 : 0;
                assert(column_scores[i] >= 0.0);
            }
        }
    }
}

/**
 * Used to make two MSAs consistent with each other for a shared sequence
 */
static void trim(int64_t row1, Msa *msa1, float *column_scores1,
                 int64_t row2, Msa *msa2, float *column_scores2, int64_t overlap) {
    if(overlap == 0) { // There is no overlap, so no need to trim either MSA
        return;
    }
    assert(overlap > 0); // Otherwise the overlap must be positive

    int64_t seq_len1 = msa1->seq_lens[row1]; // The prefix length of the forward complement sequence in the first MSA
    int64_t seq_len2 = msa2->seq_lens[row2]; // The prefix length of the reverse complement sequence in the second MSA
    // They can be different if either MSA does not include the whole sequence
    assert(overlap <= seq_len1); // The overlap must be less than the length of the prefixes
    assert(overlap <= seq_len2);

    // Get the cumulative cut scores for the columns containing the shared sequence
    float *cu_column_scores1 = st_malloc(msa1->column_no * sizeof(float));
    float *cu_column_scores2 = st_malloc(msa2->column_no * sizeof(float));
    sum_column_scores(row1, msa1, column_scores1, cu_column_scores1);
    sum_column_scores(row2, msa2, column_scores2, cu_column_scores2);

    // The score if we cut all of the overlap in msa1 and keep all of the overlap in msa2
    assert(seq_len2 <= msa2->column_no);
    float max_cut_score = cu_column_scores2[seq_len2-1];
    if(overlap < seq_len1) { // The overlap is less than the length of the first sequence
        assert(seq_len1-overlap-1 >= 0);
        max_cut_score += cu_column_scores1[seq_len1-overlap-1]; // We will keep everything before the overlap
    }
    int64_t max_overlap_cut_point = 0; // the length of the prefix of the overlap of msa1 to keep

    // Now walk through each possible cut point within the overlap
    for(int64_t i=0; i<overlap-1; i++) {
        assert(seq_len2-i-2 >= 0); // Sanity check
        float cut_score = cu_column_scores1[seq_len1-overlap+i] + cu_column_scores2[seq_len2-i-2]; // The score if we keep prefix up to
        // and including column i of MSA1's overlap, and the prefix of msa2 up to and including column seq_len-i-2
        if(cut_score > max_cut_score) {
            max_overlap_cut_point = i + 1;
            max_cut_score = cut_score;
        }
    }

    // The score if we cut all of msa2's overlap and keep all of msa1's
    float f = cu_column_scores1[seq_len1-1];
    if(overlap < seq_len2) {
        assert(seq_len2-overlap-1 >= 0);
        f += cu_column_scores2[seq_len2-overlap-1];
    }

    if(f > max_cut_score) {
        max_cut_score = f;
        max_overlap_cut_point = overlap;
    }

    // Now trim back the two MSAs
    assert(max_overlap_cut_point <= overlap);
    trim_msa_suffix(msa1, column_scores1, row1, seq_len1 - overlap + max_overlap_cut_point);
    trim_msa_suffix(msa2, column_scores2, row2, seq_len2 - max_overlap_cut_point);

    free(cu_column_scores1);
    free(cu_column_scores2);
}

/**
 * recompute the seq_lens of a trimmed msa and clip off empty suffix columns
 * (todo: can this be built into trimming code?)
 */
static void msa_fix_trimmed(Msa* msa) {
    for (int64_t i = 0; i < msa->seq_no; ++i) {
        // recompute the seq_len
        msa->seq_lens[i] = 0;
        for (int64_t j = 0; j < msa->column_no; ++j) {
            if (msa_to_base(msa->msa_seq[i][j]) != '-') {
                ++msa->seq_lens[i];
            }
        }
    }
    // trim empty columns
    int64_t empty_columns = 0;
    for (bool still_empty = true; empty_columns < msa->column_no; ++empty_columns) {
        for (int64_t i = 0; i < msa->seq_no && still_empty; ++i) {
            still_empty = msa_to_base(msa->msa_seq[i][msa->column_no - 1 - empty_columns]) == '-';
        }
        if (!still_empty) {
            break;
        }
    }
    msa->column_no -= empty_columns;
}


/* ============================================================================ */
/* Window dumping: shared between the backends                                  */
/* ============================================================================ */

/*
 * Write the window as FASTA plus its 5x5 matrix.  The matrix file format is the one both abpoa -t
 * and minipoa -m read, so either aligner can be pointed straight at it.
 */
static void dump_window_fasta_and_matrix(Msa *msa, uint8_t **bseqs, const int *mat,
                                         const char *input_path, const char *matrix_path) {
    FILE *mat_file = fopen(matrix_path, "w");
    if (mat_file != NULL) {
        fprintf(mat_file, "\tA\tC\tG\tT\tN\n");
        for (size_t i = 0; i < 5; ++i) {
            fprintf(mat_file, "%c", "ACGTN"[i]);
            for (size_t j = 0; j < 5; ++j) {
                fprintf(mat_file, "\t%d", mat[i * 5 + j]);
            }
            fprintf(mat_file, "\n");
        }
        fclose(mat_file);
    }
    FILE *fa_file = fopen(input_path, "w");
    if (fa_file == NULL) {
        return;
    }
    for (int64_t i = 0; i < msa->seq_no; ++i) {
        fprintf(fa_file, ">%" PRIi64 "\n", i);
        for (int64_t j = 0; j < msa->seq_lens[i]; ++j) {
            fputc(msa_to_base(bseqs[i][j]), fa_file);
        }
        fputc('\n', fa_file);
    }
    fclose(fa_file);
}

/* ============================================================================ */
/* abPOA backend                                                                */
/* ============================================================================ */

abpoa_para_t *abpoaParamaters_constructFromCactusParams(CactusParams *params) {
    abpoa_para_t *abpt = abpoa_init_para();

    // output options
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 0; // generate consensus sequence, set 0 to disable

    // alignment mode. 0:global alignment, 1:local, 2:extension
    // only global works
    abpt->align_mode = ABPOA_GLOBAL_MODE;

    // banding parameters
    abpt->wb = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentBandConstant");
    abpt->wf = cactusParams_get_float(params, 3, "bar", "poa", "partialOrderAlignmentBandFraction");

    // gap scoring model
    abpt->gap_open1 = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentGapOpenPenalty1");
    abpt->gap_ext1 = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentGapExtensionPenalty1");
    abpt->gap_open2 = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentGapOpenPenalty2");
    abpt->gap_ext2 = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentGapExtensionPenalty2");
    
    // seeding paramters
    abpt->disable_seeding = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentDisableSeeding");
    assert(abpt->disable_seeding == 0 || abpt->disable_seeding == 1);
    abpt->k = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentMinimizerK");
    abpt->w = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentMinimizerW");
    abpt->min_w = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentMinimizerMinW");

    // progressive toggle
    abpt->progressive_poa = cactusParams_get_int(params, 3, "bar", "poa", "partialOrderAlignmentProgressiveMode");

    // generate the substitution matrix
    abpt->use_score_matrix = 0;
    abpoa_post_set_para(abpt);

    // optionally override the substitution matrix
    char *submat_string = cactusParams_get_string(params, 3, "bar", "poa", "partialOrderAlignmentSubMatrix");
    if (submat_string && strlen(submat_string) > 0) {
        // Note, this will be used to explicitly override abpoa's subsitution matrix just before aligning
        abpt->use_score_matrix = 1;
        assert(abpt->m == 5);
        int count = 0;
        for (char* val = strtok(submat_string, " "); val != NULL; val = strtok(NULL, " ")) {
            abpt->mat[count++] = atoi(val);
        }
        assert(count == 25);
        int i; abpt->min_mis = 0, abpt->max_mat = 0;
        for (i = 0; i < abpt->m * abpt->m; ++i) {
            if (abpt->mat[i] > abpt->max_mat)
                abpt->max_mat = abpt->mat[i];
            if (-abpt->mat[i] > abpt->min_mis) 
                abpt->min_mis = -abpt->mat[i];
        }
    }
    free(submat_string);    

    return abpt;
}

// It turns out abpoa can write to these, so we make a quick copy before using
static abpoa_para_t *copy_abpoa_params(abpoa_para_t *abpt) {
    abpoa_para_t *abpt_cpy = abpoa_init_para();
    abpt_cpy->out_msa = 1;
    abpt_cpy->out_cons = 0;
    abpt_cpy->align_mode = abpt->align_mode;
    abpt_cpy->wb = abpt->wb;
    abpt_cpy->wf = abpt->wf;
    abpt_cpy->match = abpt->match;
    abpt_cpy->mismatch = abpt->mismatch;
    abpt_cpy->gap_mode = abpt->gap_mode;
    abpt_cpy->gap_open1 = abpt->gap_open1;
    abpt_cpy->gap_ext1 = abpt->gap_ext1;
    abpt_cpy->gap_open2 = abpt->gap_open2;
    abpt_cpy->gap_ext2 = abpt->gap_ext2;
    abpt_cpy->disable_seeding = abpt->disable_seeding;
    abpt_cpy->k = abpt->k;
    abpt_cpy->w = abpt->w;
    abpt_cpy->min_w = abpt->min_w;
    abpt_cpy->progressive_poa = abpt->progressive_poa;
    abpt_cpy->use_score_matrix = 0;
    abpoa_post_set_para(abpt_cpy);
    abpt_cpy->use_score_matrix = abpt->use_score_matrix;
    if (abpt->use_score_matrix == 1) {
        memcpy(abpt_cpy->mat, abpt->mat, abpt->m * abpt->m * sizeof(int));
    }
    abpt_cpy->max_mat = abpt->max_mat;
    abpt_cpy->min_mis = abpt->min_mis;
    return abpt_cpy;
}

// dump the abpoa input to files, and return a command line for running abpoa on them
char* dump_abpoa_input(Msa* msa, abpoa_para_t* abpt, uint8_t **bseqs, char* abpoa_input_path, char* abpoa_matrix_path,
                       char* abpoa_command_path, char* abpoa_output_path) {
    // dump the abpoa input sequences to a FASTA file
    FILE* dump_file = fopen(abpoa_input_path, "w");
    for (int64_t i = 0; i < msa->seq_no; ++i) {
        int64_t seq_len = msa->seq_lens[i];            
        char* buffer = (char*)malloc((seq_len + 1) * sizeof(char));
        for (int64_t j = 0; j < seq_len; ++j) {
            buffer[j] = msa_to_base(bseqs[i][j]);
        }
        buffer[msa->seq_lens[i]] = '\0';
        fprintf(dump_file, ">%ld\n%s\n", i, buffer);
        free(buffer);
    }
    fclose(dump_file);

    // dump the abpoa input matrix to file
    FILE* mat_file = fopen(abpoa_matrix_path, "w");
    fprintf(mat_file, "\tA\tC\tG\tT\tN\n");
    for (size_t i = 0; i < 5; ++i) {
        fprintf(mat_file, "%c", "ACGTN"[i]);
        for (size_t j = 0; j < 5; ++j) {
            fprintf(mat_file, "\t%d", abpt->mat[i * 5 + j]);
        }
        fprintf(mat_file, "\n");
    }
    fclose(mat_file);

    // make a command line
    char* abpoa_command = st_malloc(4096 * sizeof(char));
    sprintf(abpoa_command, "abpoa %s -O %d,%d -E %d,%d -b %d -f %lf -t %s -r 1 -m 0",
            abpoa_input_path,
            abpt->gap_open1,
            abpt->gap_open2,
            abpt->gap_ext1,
            abpt->gap_ext2,
            abpt->wb,
            abpt->wf,
            abpoa_matrix_path);
    if (!abpt->disable_seeding) {
        strcat(abpoa_command, " -S");
        char kw_opts[128];
        sprintf(kw_opts, " -k %d -w %d -n %d", abpt->k, abpt->w, abpt->min_w);
        strcat(abpoa_command, kw_opts);
    }
    if (abpt->progressive_poa) {
        strcat(abpoa_command, " -p");
    }
    strcat(abpoa_command, " > ");
    strcat(abpoa_command, abpoa_output_path);

    // dump the command line
    FILE* cmd_file = fopen(abpoa_command_path, "w");
    fprintf(cmd_file, "%s\n", abpoa_command);
    fclose(cmd_file);

    return abpoa_command;
}

#ifdef CACTUS_ABPOA_FROM_COMMAND_LINE
void abpoa_msa_from_command_line(char* abpoa_command_line, char* abpoa_output_path, uint8_t*** msa_seq, int* col_no) {
    // run abpoa
    st_system(abpoa_command_line);

    // read the result (ascii alignment) back into memory
    size_t n_rows = 0;
    size_t n_cols = 0;
    FILE* msa_file = fopen(abpoa_output_path, "r");
    int64_t buf_size = 500000;
    char* buf = st_malloc(buf_size * sizeof(char));
    
    while (benLine(&buf, &buf_size, msa_file) != -1) {
        if (strlen(buf) && buf[0] != '>') {
            ++n_rows;
        }
    }

    *msa_seq = st_malloc(n_rows * sizeof(uint8_t*));
    fclose(msa_file);
    msa_file = fopen(abpoa_output_path, "r");
    n_rows = 0;
    while (benLine(&buf, &buf_size, msa_file) != -1) {
        if (strlen(buf) && buf[0] != '>') {
            if (n_cols == 0) {
                n_cols = strlen(buf);
            } else {
                assert(n_cols == strlen(buf));
            }
            (*msa_seq)[n_rows] = st_malloc(n_cols * sizeof(uint8_t));
            for (size_t i = 0; i < n_cols; ++i) {
                (*msa_seq)[n_rows][i] = msa_to_byte(buf[i]);
            }
            ++n_rows;
        }
    }
    fclose(msa_file);
    *col_no = (int)n_cols;
    
    free(buf);
}
#endif

/*
 * The abPOA backend.  Fills msa->msa_seq and msa->column_no for one window.
 */
static void run_abpoa_window(Msa *msa, uint8_t **bseqs, PoaParameters *poa_parameters,
                             int64_t max_prog_rows, double max_prog_length_diff) {
    // init abpoa
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = copy_abpoa_params(poa_parameters->abpt);
    if (msa->seq_no > max_prog_rows ||
        // note: these are sorted by length excep in unit tests
        (1. - (double)msa->seq_lens[msa->seq_no-1] / (double)msa->seq_lens[0] > max_prog_length_diff)) {
        abpt->progressive_poa = 0;
    }
    
    // dump the input to file, if asked to at run time
    char abpoa_input_path[1024], abpoa_matrix_path[1024], abpoa_command_path[1024], abpoa_output_path[1024];
    char *abpoa_command_line = NULL;
    if (bar_dump_dir != NULL) {
        // The old name keyed off the Msa pointer, so two windows that reused the same freed
        // allocation silently overwrote each other.  pid + counter is unique for the run.
        int64_t dump_id;
#if defined(_OPENMP)
#pragma omp atomic capture
#endif
        dump_id = ++bar_dump_counter;
        sprintf(abpoa_input_path, "%s/bar_window_%d_%" PRIi64 ".fa", bar_dump_dir, (int)getpid(), dump_id);
        sprintf(abpoa_matrix_path, "%s.mat", abpoa_input_path);
        sprintf(abpoa_command_path, "%s.cmd", abpoa_input_path);
        sprintf(abpoa_output_path, "%s.out", abpoa_input_path);
        abpoa_command_line = dump_abpoa_input(msa, abpt, bseqs,
                                              abpoa_input_path, abpoa_matrix_path, abpoa_command_path, abpoa_output_path);
    }

#ifdef CACTUS_ABPOA_FROM_COMMAND_LINE
    // run abpoa from the command line
    if (abpoa_command_line == NULL) {
        st_errAbort("CACTUS_ABPOA_FROM_COMMAND_LINE needs CACTUS_BAR_DUMP_DIR set in the environment");
    }
    abpoa_msa_from_command_line(abpoa_command_line, abpoa_output_path, &(msa->msa_seq), &(msa->column_no));

    int test_cols = 0;
    uint8_t** test_msa = NULL;
    abpoa_msa(ab, abpt, msa->seq_no, NULL, msa->seq_lens, bseqs, NULL, NULL);
    // abpoa's interface has changed a bit -- instead of passing in pointers to the results, they
    // end up in the ab->abc struct -- we extract them here
    test_msa = ab->abc->msa_base;
    ab->abc->msa_base = NULL;
    test_cols = ab->abc->msa_len;

    // sanity check to make sure we get the same output
    assert(msa->column_no == test_cols);        
    for (int i = 0; i < msa->seq_no; ++i) {
      for (int j = 0; j < test_cols; ++j) {
          //todo: not sure why this doesn't work anymore !!!!
          //assert(test_msa[i][j] == msa->msa_seq[i][j]);
      }
      free(test_msa[i]);
    }
    free(test_msa);
#else
    // perform abpoa-msa
    abpoa_msa(ab, abpt, msa->seq_no, NULL, msa->seq_lens, bseqs, NULL, NULL);
    // abpoa's interface has changed a bit -- instead of passing in pointers to the results, they
    // end up in the ab->abc struct -- we extract them here
    msa->msa_seq = ab->abc->msa_base;
    ab->abc->msa_base = NULL;
    msa->column_no = ab->abc->msa_len;
#endif

    // The dumps are kept.  Deleting them on success made the switch useless for its main
    // job -- looking at a window that aligned badly rather than one that crashed.
    if (abpoa_command_line != NULL) {
        free(abpoa_command_line);
    }

    // free abpoa
    abpoa_free(ab);
    abpoa_free_para(abpt);
}

/* ============================================================================ */
/* minipoa backend                                                              */
/* ============================================================================ */

/*
 * minipoa's substitution matrix comes from <bar><poa>; its gap penalties do not.
 *
 * Sharing the matrix keeps the two engines scoring substitutions identically, and means
 * last_scoring.py's learned matrix reaches both (it writes into <poa>, which local_alignment.py
 * also reads to score FastGA PAFs).
 *
 * The gaps have to be separate.  abPOA's gap model is convex -- min(open1 + L*ext1, open2 +
 * L*ext2) -- so with the shipped 400/30 and 1200/1 its effective extension beyond L~28 is 1.
 * minipoa has a single affine piece, so handing it abPOA's first-piece extension of 30 prices
 * long gaps about 30x above what abPOA charges, and it responds by packing bases into shared
 * columns instead of opening a gap.  That cost about six points of mafComparator accuracy on
 * evolver mammals.
 *
 * last-train reaches minipoa too: it fits a single affine model, which is minipoa's model
 * exactly, so last_scoring.py writes the learned open/extend straight into <minipoa>.  The
 * GapOpen2/GapExtend2 pair it synthesises for abPOA is a stability workaround for that aligner
 * and is deliberately not passed on.
 */
#ifdef HAVE_MINIPOA
static minipoa_para_t *minipoaParameters_constructFromCactusParams(CactusParams *params, PoaParameters *out) {
    minipoa_para_t *mpt = minipoa_init_para();
    if (mpt == NULL) {
        st_errAbort("Failed to allocate minipoa parameters: %s", minipoa_last_error());
    }

    /*
     * minipoa's own matrix, falling back to abPOA's when it is left empty -- that fallback is the
     * only way to get last-train's learned scores into minipoa, since last_scoring.py writes into
     * <poa>.
     */
    char *submat_string = cactusParams_get_string(params, 3, "bar", "minipoa", "minipoaSubMatrix");
    if (submat_string == NULL || strlen(submat_string) == 0) {
        free(submat_string);
        submat_string = cactusParams_get_string(params, 3, "bar", "poa", "partialOrderAlignmentSubMatrix");
    }
    if (submat_string != NULL && strlen(submat_string) > 0) {
        int mat[25];
        int count = 0;
        for (char *val = strtok(submat_string, " "); val != NULL && count < 25; val = strtok(NULL, " ")) {
            mat[count++] = atoi(val);
        }
        if (count != 25) {
            st_errAbort("<bar><poa partialOrderAlignmentSubMatrix> needs 25 values, got %d", count);
        }
        minipoa_set_score_matrix(mpt, mat);
        memcpy(out->mat, mat, sizeof(mat));
    }
    free(submat_string);

    /*
     * Only the first gap piece.  abPOA takes min(open1 + L*ext1, open2 + L*ext2) and minipoa has
     * no second piece at all, so gaps longer than where the two cross (~28bp with the shipped
     * 400/30 and 1200/1) are penalised more heavily here than abPOA would.
     */
    out->gapOpen = cactusParams_get_int(params, 3, "bar", "minipoa", "minipoaGapOpenPenalty");
    out->gapExt = cactusParams_get_int(params, 3, "bar", "minipoa", "minipoaGapExtensionPenalty");
    minipoa_set_gap(mpt, out->gapOpen, out->gapExt);

    out->bandConstant = cactusParams_get_int(params, 3, "bar", "minipoa", "minipoaBandConstant");
    out->bandFraction = cactusParams_get_float(params, 3, "bar", "minipoa", "minipoaBandFraction");
    minipoa_set_band(mpt, out->bandConstant, out->bandFraction);
    out->adaptiveBand = cactusParams_get_int(params, 3, "bar", "minipoa", "minipoaAdaptiveBand");
    minipoa_set_adaptive_band(mpt, out->adaptiveBand);

    out->seeding = !cactusParams_get_int(params, 3, "bar", "minipoa", "minipoaDisableSeeding");
    out->minimizerK = cactusParams_get_int(params, 3, "bar", "minipoa", "minipoaMinimizerK");
    out->minimizerW = cactusParams_get_int(params, 3, "bar", "minipoa", "minipoaMinimizerW");
    out->anchorWindow = cactusParams_get_int(params, 3, "bar", "minipoa", "minipoaAnchorWindow");
    minipoa_set_seeding(mpt, out->seeding, out->minimizerK, out->minimizerW, out->anchorWindow);

    /*
     * Progressive ordering, on by default, because abPOA runs with it on
     * (<poa partialOrderAlignmentProgressiveMode>) and the order sequences are added to a POA
     * graph changes the alignment -- markedly so on diverged input.  Leaving it off here was worth
     * about six points of mafComparator accuracy on the evolver mammals set.
     */
    out->progressive = cactusParams_get_int(params, 3, "bar", "minipoa", "minipoaProgressiveMode");
    out->progressiveMaxRows = cactusParams_get_int(params, 3, "bar", "minipoa", "minipoaProgressiveMaxRows");
    minipoa_set_progressive(mpt, out->progressive);
    return mpt;
}
#else
/*
 * minipoa=off in include.mk.  Selecting it is a configuration error rather than a crash, and the
 * message has to name the build switch -- the config is portable between machines, the build is
 * not.
 */
static void *minipoaParameters_constructFromCactusParams(CactusParams *params, PoaParameters *out) {
    (void)params; (void)out;
    st_errAbort("<bar baseAligner=\"minipoa\"> was selected, but this cactus was built with "
                "minipoa=off (see include.mk). Rebuild with minipoa=on, or choose abpoa or pecan.");
    return NULL;
}
#endif

/*
 * A command line that replays this window through minipoa.  Gap penalties are negated: minipoa
 * maximises, so penalties are negative there, while cactus and abpoa carry them positive.
 */
static char *dump_minipoa_input(Msa *msa, PoaParameters *pp, uint8_t **bseqs, char *input_path,
                                char *matrix_path, char *command_path, char *output_path) {
    dump_window_fasta_and_matrix(msa, bseqs, pp->mat, input_path, matrix_path);

    int f = pp->bandFraction > 0.0 ? (int)(1.0 / pp->bandFraction + 0.5) : 0;
    char *command = st_malloc(4096 * sizeof(char));
    sprintf(command, "minipoa %s -m %s -O -%d -E -%d -b %d -f %d -r 1 -t 1",
            input_path, matrix_path, pp->gapOpen, pp->gapExt, pp->bandConstant, f);
    if (pp->seeding) {
        char kw_opts[128];
        sprintf(kw_opts, " -S -k %d -w %d", pp->minimizerK, pp->minimizerW);
        strcat(command, kw_opts);
        if (pp->anchorWindow > 0) {
            sprintf(kw_opts, " -W %d", pp->anchorWindow);
            strcat(command, kw_opts);
        }
    }
    // -p and -B are on in the shipped config, and progressive ordering in particular changes the
    // alignment, so a replay that omitted them would not reproduce the window it is meant to
    // explain.
    if (pp->progressive) {
        strcat(command, " -p");
    }
    if (pp->adaptiveBand) {
        strcat(command, " -B");
    }
    strcat(command, " > ");
    strcat(command, output_path);

    FILE *cmd_file = fopen(command_path, "w");
    if (cmd_file != NULL) {
        fprintf(cmd_file, "%s\n", command);
        fclose(cmd_file);
    }
    return command;
}

/*
 * The minipoa half of run_poa_window().
 *
 * minipoa_msa() returns a status rather than exiting: its internal band/backtrack traps used to
 * call exit(0) -- a success code -- which from here is indistinguishable from a clean run with
 * truncated output.  A failure here is fatal for the job either way, but it says which window and
 * leaves it on disk when dumping is on, which is the difference between a bug report and a shrug.
 */
static void run_minipoa_window(Msa *msa, uint8_t **bseqs, PoaParameters *poa_parameters) {
#ifndef HAVE_MINIPOA
    (void)msa; (void)bseqs; (void)poa_parameters;
    st_errAbort("minipoa was selected but this cactus was built with minipoa=off (see include.mk)");
#else
    char input_path[1024], matrix_path[1024], command_path[1024], output_path[1024];
    char *command_line = NULL;
    if (bar_dump_dir != NULL) {
        int64_t dump_id;
#if defined(_OPENMP)
#pragma omp atomic capture
#endif
        dump_id = ++bar_dump_counter;
        sprintf(input_path, "%s/bar_window_%d_%" PRIi64 ".fa", bar_dump_dir, (int)getpid(), dump_id);
        sprintf(matrix_path, "%s.mat", input_path);
        sprintf(command_path, "%s.cmd", input_path);
        sprintf(output_path, "%s.out", input_path);
        command_line = dump_minipoa_input(msa, poa_parameters, bseqs, input_path, matrix_path,
                                          command_path, output_path);
    }

    const minipoa_para_t *mpt = (const minipoa_para_t *)poa_parameters->mpt;
    if (poa_parameters->mptNoProgressive != NULL && msa->seq_no > poa_parameters->progressiveMaxRows) {
        mpt = (const minipoa_para_t *)poa_parameters->mptNoProgressive;
    }

    int column_no = 0;
    uint8_t **msa_seq = NULL;
    int ret = minipoa_msa(mpt, (int)msa->seq_no, msa->seq_lens, bseqs, &msa_seq, &column_no);
    if (ret != 0) {
        st_errAbort("minipoa failed on a %" PRIi64 " x %d window: %s. %s",
                    msa->seq_no, msa->seq_lens[0], minipoa_last_error(),
                    command_line != NULL
                        ? command_path
                        : "Set CACTUS_BAR_DUMP_DIR to capture the window that did this.");
    }
    msa->msa_seq = msa_seq;
    msa->column_no = column_no;
    /*
     * A global alignment cannot have fewer columns than its longest input row, so this catches a
     * truncated or empty result -- which is the shape minipoa failures take when they do not
     * return an error.  O(rows), next to nothing against the DP, and worth it: an all-gap MSA
     * produces zero alignment blocks and BAR would otherwise drop the flower in silence.
     */
    int longest = 0;
    for (int64_t i = 0; i < msa->seq_no; i++) {
        if (msa->seq_lens[i] > longest) {
            longest = msa->seq_lens[i];
        }
    }
    if (column_no < longest) {
        st_errAbort("minipoa returned %d columns for a %" PRIi64 " x %d window, which cannot be a "
                    "global alignment of it. %s", column_no, msa->seq_no, longest,
                    command_line != NULL
                        ? command_path
                        : "Set CACTUS_BAR_DUMP_DIR to capture the window that did this.");
    }
    if (command_line != NULL) {
        free(command_line);
    }
#endif
}

/* ============================================================================ */
/* Base aligner selection and dispatch                                          */
/* ============================================================================ */

/*
 * <bar baseAligner="..."> picks the engine.  It is read with cactusParams_has so that a config
 * predating the attribute -- including any a user has saved -- still selects what it used to via
 * the older <bar partialOrderAlignment="0|1"> boolean.
 */
BaseAligner baseAligner_constructFromCactusParams(CactusParams *params) {
    /*
     * Both attributes are optional, in both directions.  An old config has only
     * partialOrderAlignment; a config written from now on may reasonably have only baseAligner --
     * including one produced by following the warning below, which tells the user to delete the
     * legacy attribute.  Reading either unguarded would st_errAbort on the other's config.
     */
    bool hasLegacy = cactusParams_has(params, 2, "bar", "partialOrderAlignment");
    int64_t usePoa = hasLegacy ? cactusParams_get_int(params, 2, "bar", "partialOrderAlignment") : 1;
    if (!cactusParams_has(params, 2, "bar", "baseAligner")) {
        return usePoa ? BASE_ALIGNER_ABPOA : BASE_ALIGNER_PECAN;
    }
    char *name = cactusParams_get_string(params, 2, "bar", "baseAligner");
    BaseAligner engine;
    if (strcmp(name, "pecan") == 0) {
        engine = BASE_ALIGNER_PECAN;
    } else if (strcmp(name, "abpoa") == 0) {
        engine = BASE_ALIGNER_ABPOA;
    } else if (strcmp(name, "minipoa") == 0) {
        engine = BASE_ALIGNER_MINIPOA;
    } else {
        st_errAbort("Unknown <bar baseAligner=\"%s\">; expected pecan, abpoa or minipoa", name);
        engine = BASE_ALIGNER_ABPOA; /* not reached */
    }
    free(name);

    /*
     * Both attributes present and disagreeing is worth saying out loud.  partialOrderAlignment="0"
     * is how the config has always documented "use pecan", so someone who sets it and gets abpoa
     * anyway should not have to discover that from the alignment.
     */
    bool poaImplied = engine != BASE_ALIGNER_PECAN;
    if (hasLegacy && (usePoa != 0) != poaImplied) {
        st_logCritical("Warning: <bar baseAligner=\"%s\"> overrides <bar partialOrderAlignment=\"%" PRIi64
                       "\">, which asks for the opposite. baseAligner wins; remove the other to silence this.\n",
                       baseAligner_toString(engine), usePoa);
    }
    return engine;
}

const char *baseAligner_toString(BaseAligner engine) {
    switch (engine) {
        case BASE_ALIGNER_PECAN: return "pecan";
        case BASE_ALIGNER_ABPOA: return "abpoa";
        case BASE_ALIGNER_MINIPOA: return "minipoa";
    }
    return "unknown";
}

PoaParameters *poaParameters_constructFromCactusParams(CactusParams *params, BaseAligner engine) {
    if (engine == BASE_ALIGNER_PECAN) {
        return NULL;
    }
    PoaParameters *poaParameters = st_calloc(1, sizeof(PoaParameters));
    poaParameters->engine = engine;
    if (engine == BASE_ALIGNER_ABPOA) {
        poaParameters->abpt = abpoaParamaters_constructFromCactusParams(params);
    } else {
        poaParameters->mpt = minipoaParameters_constructFromCactusParams(params, poaParameters);
#ifdef HAVE_MINIPOA
        if (poaParameters->progressive) {
            // Same settings with the guide tree off, for windows too wide to afford a dense NxN
            // distance matrix.  abPOA caps this the same way, via partialOrderAlignmentProgressiveMaxRows.
            PoaParameters scratch = *poaParameters;
            minipoa_para_t *plain = minipoaParameters_constructFromCactusParams(params, &scratch);
            minipoa_set_progressive(plain, 0);
            poaParameters->mptNoProgressive = plain;
        }
#endif
    }
    return poaParameters;
}

void poaParameters_destruct(PoaParameters *poaParameters) {
    if (poaParameters == NULL) {
        return;
    }
    if (poaParameters->abpt != NULL) {
        abpoa_free_para(poaParameters->abpt);
    }
#ifdef HAVE_MINIPOA
    if (poaParameters->mpt != NULL) {
        minipoa_free_para((minipoa_para_t *)poaParameters->mpt);
    }
    if (poaParameters->mptNoProgressive != NULL) {
        minipoa_free_para((minipoa_para_t *)poaParameters->mptNoProgressive);
    }
#endif
    free(poaParameters);
}

/*
 * Align one window with whichever engine was selected.
 *
 * The contract both backends meet: seq_no rows of column_no bytes, one malloc per row plus one
 * for the row array (msa_destruct frees them that way), rows in input order, values 0-4 for ACGTN
 * and 5 for a gap.  Everything around this -- the sliding window, the empty-sequence hack,
 * trimming, stitching, block extraction -- is engine-neutral and shared.
 */
static void run_poa_window(Msa *msa, uint8_t **bseqs, PoaParameters *poa_parameters,
                           int64_t max_prog_rows, double max_prog_length_diff) {
    switch (poa_parameters->engine) {
        case BASE_ALIGNER_ABPOA:
            run_abpoa_window(msa, bseqs, poa_parameters, max_prog_rows, max_prog_length_diff);
            break;
        case BASE_ALIGNER_MINIPOA:
            run_minipoa_window(msa, bseqs, poa_parameters);
            break;
        default:
            st_errAbort("run_poa_window called with base aligner %s, which produces no MSA",
                        baseAligner_toString(poa_parameters->engine));
    }
}






/* ============================================================================ */
/* Shared: windowed MSA construction, engine-neutral                            */
/* ============================================================================ */


Msa *msa_make_partial_order_alignment(char **seqs, int *seq_lens, int64_t seq_no, int64_t window_size,
                                      int64_t max_prog_rows, double max_prog_length_diff, PoaParameters *poa_parameters) {

    assert(seq_no > 0);

    // only one input sequence: no point sending into abpoa; just return it instead
    // (note: current version of abpoa will crash in progressive mode on one sequence)
    // todo: can we filter this out at higher level?
    if (seq_no == 1) {
        Msa *msa = st_malloc(sizeof(Msa));
        msa->seq_no = seq_no;
        msa->seqs = seqs;
        msa->seq_lens = seq_lens;
        msa->column_no = seq_lens[0];
        msa->msa_seq = st_malloc(sizeof(uint8_t*));
        msa->msa_seq[0] = st_malloc(msa->column_no * sizeof(uint8_t));
        for (int64_t i = 0; i < msa->column_no; ++i) {
            msa->msa_seq[0][i] = msa_to_byte(msa->seqs[0][i]);
        }
        return msa;
    }
    
    // we overlap the sliding window, and use the trimming logic to find the best cut point between consecutive windows
    // todo: cli-facing parameter
    float window_overlap_frac = 0.5;
    int64_t window_overlap_size = window_overlap_frac * window_size;
    if (window_overlap_size > 0) {
        --window_overlap_size; // don't want empty window when fully trimmed on each end
    }
    // keep track of what's left to align for the sliding window
    int64_t bases_remaining = 0;
    // keep track of current offsets
    int64_t* seq_offsets = (int64_t*)st_calloc(seq_no, sizeof(int64_t));
    // keep track of empty chunks
    bool* empty_seqs = (bool*)st_calloc(seq_no, sizeof(bool));
    // keep track of overlaps
    int64_t* row_overlaps = (int64_t*)st_calloc(seq_no, sizeof(int64_t));

    // allocate the poa input buffer
    uint8_t **bseqs = (uint8_t**)st_malloc(sizeof(uint8_t*) * seq_no);
    for (int64_t i = 0; i < seq_no; ++i) {
        int64_t row_size = seq_lens[i] < window_size ? seq_lens[i] : window_size;
        bseqs[i] = (uint8_t*)st_malloc(sizeof(uint8_t) * row_size);
        bases_remaining += seq_lens[i];
    }
     
    // collect our windowed outputs here, to be stiched at the end. 
    stList* msa_windows = stList_construct3(0, (void(*)(void *)) msa_destruct);
    
    // remember the previous window
    Msa* prev_msa = NULL;
    
    int64_t prev_bases_remaining = bases_remaining;
    for (int64_t iteration = 0; bases_remaining > 0; ++iteration) {

        // compute the number of bases this msa will overlap with the previous msa per row,
        // assuming that the alignments overlap by window_overlap_size
        if (prev_msa != NULL) {
            for (int64_t i = 0; i < seq_no; ++i) {
                assert(prev_msa->column_no > window_overlap_size);
                row_overlaps[i] = 0;
                for (int64_t j = prev_msa->column_no - window_overlap_size; j < prev_msa->column_no; ++j) {
                    if (msa_to_base(prev_msa->msa_seq[i][j]) != '-') {
                        ++row_overlaps[i];
                    }
                }
                // take the overlaps into account in other other counters
                assert(seq_offsets[i] >= row_overlaps[i]);
                seq_offsets[i] -= row_overlaps[i];
                bases_remaining += row_overlaps[i];
            }
        }

        // Make Msa object
        Msa *msa = st_malloc(sizeof(Msa));
        msa->seq_no = seq_no;
        msa->seqs = NULL;
        msa->seq_lens = st_malloc(sizeof(int) * msa->seq_no);
        
        // load up to window_size of each sequence into the input matrix for poa
        for (int64_t i = 0; i < msa->seq_no; ++i) {
            msa->seq_lens[i] = 0;
            for (int64_t j = seq_offsets[i]; j < seq_lens[i] && msa->seq_lens[i] < window_size; ++j, ++msa->seq_lens[i]) {
                // todo: support iupac characters?
                bseqs[i][msa->seq_lens[i]] = msa_to_byte(seqs[i][j]);
            }
        }

        // poa can't handle empty sequences.  this is a hack to get around that
        int emptyCount = 0;
        for (int64_t i = 0; i < msa->seq_no; ++i) {
            if (msa->seq_lens[i] == 0) {
                empty_seqs[i] = true;
                msa->seq_lens[i] = 1;
                bseqs[i][0] = msa_to_byte('N');
                ++emptyCount;
            } else {
                empty_seqs[i] = false;
            }
        }

        run_poa_window(msa, bseqs, poa_parameters, max_prog_rows, max_prog_length_diff);

        // mask out empty sequences that were phonied in as Ns above
        for (int64_t i = 0; i < msa->seq_no && emptyCount > 0; ++i) {
            if (empty_seqs[i] == true) {
                for (int j = 0; j < msa->column_no; ++j) {
                    if (msa_to_base(msa->msa_seq[i][j]) != '-') {
                        assert(msa_to_base(msa->msa_seq[i][j]) == 'N');
                        msa->msa_seq[i][j] = msa_to_byte('-');
                        --msa->seq_lens[i];
                        assert(msa->seq_lens[i] == 0);
                        --emptyCount;
                        break;
                    }
                }
            }
        }
        assert(emptyCount == 0);

        //if (prev_msa) {
        //    fprintf(stderr, "PREV MSA\n");
        //    msa_print(prev_msa, stderr);
        //}
        //fprintf(stderr, "CUR MSA\n");
        //msa_print(msa, stderr);
        // remember how much we aligned this round
        for (int64_t i = 0; i < msa->seq_no; ++i) {
            //////////////////////////////////////////////////////////////////////////////////////
            // todo: why is this hack necessary?  using it in order for trim to work properly   //
            // after abpoa switched to weirdo 256-bit values  (nst_nt256_table)                //
            for (int64_t j = 0; j < msa->column_no; ++j) {
                msa->msa_seq[i][j] = msa_to_byte(msa_to_base(msa->msa_seq[i][j]));
            }
            bases_remaining -= msa->seq_lens[i];
            seq_offsets[i] += msa->seq_lens[i];
        }

        // todo: there is obviously room for optimization here, as we compute full scores twice for each msa
        //       in addition to flipping the prev_msa back and forth
        //       (not sure if this is at all noticeable on top of abpoa running time though)
        if (prev_msa) {
            // trim() presently assumes we're looking at reverse-complement sequence:
            flip_msa_seq(msa);
            float* prev_column_scores = make_column_scores(prev_msa);
            float* column_scores = make_column_scores(msa);

            // trim with the previous alignment
            for (int64_t i = 0; i < msa->seq_no; ++i) {
                int64_t overlap = msa->seq_lens[i] < row_overlaps[i] ? msa->seq_lens[i] : row_overlaps[i];
                if (overlap > 0) {
                    trim(i, msa, column_scores, i, prev_msa, prev_column_scores, overlap);
                }
            }
            // todo: can this be done as part of trim?
            msa_fix_trimmed(msa);
            msa_fix_trimmed(prev_msa);            
            // flip our msa back to its original strand
            flip_msa_seq(msa);

            free(prev_column_scores);
            free(column_scores);
        }

        // add the msa to our list
        stList_append(msa_windows, msa);
        
        // sanity check        
        assert(prev_bases_remaining > bases_remaining && bases_remaining >= 0);

        prev_msa = msa;
        
        //used only for sanity check
        prev_bases_remaining = bases_remaining; 
    }

    int64_t num_windows = stList_length(msa_windows);
    Msa *output_msa;
    if (num_windows == 1) {
        // if we have only one window, return it
        output_msa = stList_removeFirst(msa_windows);
        output_msa->seqs = seqs;
        free(output_msa->seq_lens); // cleanup old memory
        output_msa->seq_lens = seq_lens;
    } else {
        // otherwise, we stitch all the window msas into a new output msa
        output_msa = st_malloc(sizeof(Msa));
        assert(seq_no > 0);
        output_msa->seq_no = seq_no;
        output_msa->seqs = seqs;
        output_msa->seq_lens = seq_lens;
        output_msa->column_no = 0;
        for (int64_t i = 0; i < num_windows; ++i) {
            Msa* msa_i = (Msa*)stList_get(msa_windows, i);
            output_msa->column_no += msa_i->column_no;
        }
        output_msa->msa_seq = st_malloc(sizeof(uint8_t *) * output_msa->seq_no);
        for (int64_t i = 0; i < output_msa->seq_no; ++i) {
            output_msa->msa_seq[i] = st_malloc(sizeof(uint8_t) * output_msa->column_no);
            int64_t offset = 0;
            for (int64_t j = 0; j < num_windows; ++j) {
                Msa* msa_j = stList_get(msa_windows, j);
                uint8_t* window_row = msa_j->msa_seq[i];
                for (int64_t k = 0; k < msa_j->column_no; ++k) {
                    output_msa->msa_seq[i][offset++] = window_row[k];
                }
            }
            assert(offset == output_msa->column_no);
        }
    } 

    // Clean up
    for (int64_t i = 0; i < seq_no; ++i) {
        free(bseqs[i]);
    }
    free(bseqs);
    free(seq_offsets);
    free(empty_seqs);
    free(row_overlaps);
    stList_destruct(msa_windows);

    return output_msa;
}

/*
 * The partial order alignments below used to run under a nested OpenMP region for flowers
 * with at least 50 ends.  That is gone.  A nested num_threads(k) region creates k new
 * threads per outer thread rather than subdividing the outer team, so it multiplied the
 * number of concurrent abPOA instances -- and with them this phase's peak memory -- while
 * buying no wall-clock time.  Measured on VGP data at equal runtime: 2.2x the memory on
 * MammalsAnc0 (208 ends at most) and 4.9x on AnuraAnc6 (847 ends), the penalty growing
 * with flower size, so it cost most on exactly the largest problems.  The parallelism
 * belongs in the flower loop in bar.c, which already saturates the thread count.
 */

Msa **make_consistent_partial_order_alignments(int64_t end_no, int64_t *end_lengths, char ***end_strings,
        int **end_string_lengths, int64_t **right_end_indexes, int64_t **right_end_row_indexes, int64_t **overlaps,
        int64_t window_size, int64_t max_prog_rows, double max_prog_length_diff, PoaParameters *poa_parameters) {
    // Calculate the initial, potentially inconsistent msas and column scores for each msa
    float *column_scores[end_no];
    Msa **msas = st_malloc(sizeof(Msa *) * end_no);

    for(int64_t i=0; i<end_no; i++) {
        msas[i] = msa_make_partial_order_alignment(end_strings[i], end_string_lengths[i], end_lengths[i], window_size,
                                                   max_prog_rows, max_prog_length_diff, poa_parameters);
        column_scores[i] = make_column_scores(msas[i]);
    }

    // Make the msas consistent with one another
    for(int64_t i=0; i<end_no; i++) { // For each end
        Msa *msa = msas[i];
        for(int64_t j=0; j<msa->seq_no; j++) { //  For each string incident to the ith end
            int64_t right_end_index = right_end_indexes[i][j]; // Find the other end it is incident with
            int64_t right_end_row_index = right_end_row_indexes[i][j]; // And the index of its reverse complement

            // If it hasn't already been trimmed
            if(right_end_index > i || (right_end_index == i /* self loop */ && right_end_row_index > j)) {
                trim(j, msa, column_scores[i],
                        right_end_row_index, msas[right_end_index], column_scores[right_end_index], overlaps[i][j]);
            }
        }
    }

    // Cleanup
    for(int64_t i=0; i<end_no; i++) {
        free(column_scores[i]);
    }

    return msas;
}

/**
 * The follow code is for dealing with the cactus API
 */


void alignmentBlock_destruct(AlignmentBlock *alignmentBlock) {
    AlignmentBlock *a;
    while(alignmentBlock != NULL) {
        a = alignmentBlock;
        alignmentBlock = alignmentBlock->next;
        free(a);
    }
}

/* ============================================================================ */
/* Shared: adjacency strings out of the cactus graph                            */
/* ============================================================================ */


char *get_adjacency_string(Cap *cap, int64_t *length, bool return_string) {
    assert(!cap_getSide(cap));
    Sequence *sequence = cap_getSequence(cap);
    assert(sequence != NULL);
    Cap *cap2 = cap_getAdjacency(cap);
    assert(cap2 != NULL);
    assert(cap_getSide(cap2));
    if (cap_getStrand(cap)) {
        assert(cap_getCoordinate(cap2) > cap_getCoordinate(cap));
        *length = cap_getCoordinate(cap2) - cap_getCoordinate(cap) - 1;
        assert(*length >= 0);
        return return_string ? sequence_getString(sequence, cap_getCoordinate(cap) + 1, *length, 1) : NULL;
    } else {
        assert(cap_getCoordinate(cap) > cap_getCoordinate(cap2));
        *length = cap_getCoordinate(cap) - cap_getCoordinate(cap2) - 1;
        assert(*length >= 0);
        return return_string ? sequence_getString(sequence, cap_getCoordinate(cap2) + 1, *length, 0) : NULL;
    }
}

/*
 * Materialises [offset, offset+length) of the adjacency without building the whole thing. An
 * adjacency can span most of a chromosome, which is gigabytes of ASCII, and BAR holds one per
 * thread.
 */
static char *get_adjacency_substring(Cap *cap, int64_t offset, int64_t length) {
    assert(!cap_getSide(cap));
    assert(offset >= 0 && length >= 0);
    Sequence *sequence = cap_getSequence(cap);
    Cap *cap2 = cap_getAdjacency(cap);
    if (cap_getStrand(cap)) {
        // base i of the adjacency is at coordinate cap+1+i
        return sequence_getString(sequence, cap_getCoordinate(cap) + 1 + offset, length, 1);
    }
    // on the reverse strand the adjacency is the reverse complement of [cap2+1, cap-1], so base i
    // is at coordinate cap-1-i and the window we want starts that far in from the other end
    int64_t total = cap_getCoordinate(cap) - cap_getCoordinate(cap2) - 1;
    assert(offset + length <= total);
    return sequence_getString(sequence, cap_getCoordinate(cap2) + 1 + (total - offset - length), length, 0);
}

/**
 * Used to find where a run of masked (hard or soft) of at least mask_filter bases starts
 * @param seq : The string
 * @param seq_length : The length of the string
 * @param length : The maximum length we want to search in
 * @param reversed : If true, scan from the end of the string
 * @param mask_filter : Cut a string as soon as we hit more than this many hard or softmasked bases (cut is before first masked base)
 * @return length of the filtered string
 */
static int64_t get_unmasked_length(char* seq, int64_t seq_length, int64_t length, bool reversed, int64_t mask_filter) {
    if (mask_filter >= 0) {
        int64_t run_start = -1;
        for (int64_t i = 0; i < length; ++i) {
            char base = reversed ? seq[seq_length - 1 - i] : seq[i];
            if (islower(base) || base == 'N') {
                if (run_start == -1) {
                    // start masked run
                    run_start = i;
                }
                if (i + 1 - run_start > mask_filter) {
                    // our run exceeds the mask_filter, cap before the first masked base
                    return run_start;
                }
            } else {
                run_start = -1;
            }
        }
    }
    return length;
}

/**
 * Used to get a prefix of a given adjacency sequence.
 * @param seq_length
 * @param length
 * @param overlap
 * @param max_seq_length
 * @return
 */
char *get_adjacency_string_and_overlap(Cap *cap, int *length, int64_t *overlap, int64_t max_seq_length, int64_t mask_filter) {
    // Take the adjacency's length without materialising it: only the prefix we keep, and the
    // suffix the mask filter scans, are ever fetched.
    int64_t seq_length;
    get_adjacency_string(cap, &seq_length, 0);
    assert(seq_length >= 0);

    // The prefix, up to max_seq_length
    int64_t length_forward = seq_length > max_seq_length ? max_seq_length : seq_length;
    bool have_whole_adjacency = length_forward == seq_length;
    char *adjacency_string = get_adjacency_substring(cap, 0, length_forward);
    int64_t length_backward = length_forward;

    if (mask_filter >= 0) {
        // apply the mask filter on the forward strand
        length_forward = get_unmasked_length(adjacency_string, length_forward, length_forward, false, mask_filter);
        if (have_whole_adjacency) {
            length_backward = get_unmasked_length(adjacency_string, seq_length, length_forward, true, mask_filter);
        } else {
            // the backward scan reads the last length_forward bases, which are not in the prefix
            char *suffix = get_adjacency_substring(cap, seq_length - length_forward, length_forward);
            length_backward = get_unmasked_length(suffix, length_forward, length_forward, true, mask_filter);
            free(suffix);
        }
    }

    // Cleanup the string
    adjacency_string[length_forward] = '\0'; // Terminate the string at the given length
    char *c = stString_copy(adjacency_string);
    free(adjacency_string);
    adjacency_string = c;

    // Calculate the overlap with the reverse complement
    if (length_forward + length_backward > seq_length) { // There is overlap
        *overlap = length_forward + length_backward - seq_length;
        assert(*overlap >= 0);
    } else { // There is no overlap
        *overlap = 0;
    }

    // Bounded by max_seq_length (<bar bandingLimit>), and abpoa takes its lengths as int
    assert(length_forward <= INT_MAX);
    *length = (int)length_forward;

    return adjacency_string;
}

/* ============================================================================ */
/* Shared: MSA to alignment blocks to pinches                                   */
/* ============================================================================ */


/**
 * Gets the length and sequences present in the next maximal gapless alignment block.
 * @param msa The msa to scan
 * @param start The start of the gapless block
 * @param rows_in_block A boolean array of which sequences are present in the block
 * @param sequences_in_block The number of in the block
 * @return
 */
int64_t get_next_maximal_block_dimensions(Msa *msa, int64_t start, bool *rows_in_block, int64_t *sequences_in_block) {
    assert(start < msa->column_no);

    // Calculate which sequences are in the block
    *sequences_in_block = 0;
    for(int64_t i=0; i<msa->seq_no; i++) {
        rows_in_block[i] = msa_to_base(msa->msa_seq[i][start]) != '-';
        if(rows_in_block[i]) {
            *sequences_in_block += 1;
        }
    }

    // Calculate the maximal block length by looking at successive columns of the MSA and
    // checking they have the same set of sequences present as in the first block
    int64_t end = start;
    while(++end < msa->column_no) {
        for(int64_t i=0; i<msa->seq_no; i++) {
            bool p = msa_to_base(msa->msa_seq[i][end]) != '-'; // Is not a gap
            if(p != rows_in_block[i]) {
                return end;
            }
        }
    }
    return end;
}

/**
 * Make an alignment block for the given interval and sequences
 * @param seq_no The number of sequences in the MSA
 * @param start The start, inclusive, of the block
 * @param length The of the block
 * @param rows_in_block An array specifying which sequences are in the block
 * @param seq_indexes The start coordinates of the sequences in the block
 * @param row_indexes_to_caps The Caps corresponding to the sequences in the block
 * @return The new alignment block
 */
AlignmentBlock *make_alignment_block(int64_t seq_no, int64_t start, int64_t length, bool *rows_in_block,
                                     int64_t *seq_indexes, Cap **row_indexes_to_caps) {
    AlignmentBlock *pB = NULL, *block = NULL;
    for(int64_t i=0; i<seq_no; i++) { // For each row
        if(rows_in_block[i]) { // If the row is in the block
            // Make an alignment block
            AlignmentBlock *b = st_calloc(1, sizeof(AlignmentBlock));
            Cap *cap = row_indexes_to_caps[i];
            assert(!cap_getSide(cap));
            assert(cap_getSequence(cap) != NULL);
            assert(length > 0);

            b->strand = cap_getStrand(cap);
            b->length = length;
            // Calculate the sequence coordinate using Cactus coordinates
            if(b->strand) {
                b->subsequenceIdentifier = cap_getName(cap);
                b->position = cap_getCoordinate(cap) + 1 + seq_indexes[i];
                assert(b->position >= 0);
                assert(b->position + length <= cap_getCoordinate(cap_getAdjacency(cap)));
            }
            else { // In the alignment block all the coordinates are reported with respect to the positive strand sequence
                Cap *adjacentCap = cap_getAdjacency(cap);
                assert(adjacentCap != NULL);
                b->subsequenceIdentifier = cap_getName(adjacentCap);
                b->position = cap_getCoordinate(cap) - seq_indexes[i] - length;
                assert(b->position >= 0);
                assert(b->position + length <= cap_getCoordinate(cap));
                assert(b->position > cap_getCoordinate(adjacentCap));
            }

            // If this is not the first sequence in the block link to the previous sequence in the block
            if (pB != NULL) {
                pB->next = b;
                pB = b;
                assert(b->next == NULL);
            } else { // Otherwise this is the first sequence in the block
                block = b;
                pB = b;
            }
        }
    }
    assert(block != NULL);
    return block;
}

void alignmentBlock_print(AlignmentBlock *ab, FILE *f) {
    fprintf(f, "Alignment block:\n");
    while(ab != NULL) {
        fprintf(f, "\tName: %" PRIi64 "\tPosition: %" PRIi64"\tStrand: %i\tLength: %" PRIi64 "\n",
                ab->subsequenceIdentifier, ab->position, (int)ab->strand, ab->length);
        ab = ab->next;
    }
    fprintf(f, "\n");
}

/**
 * Converts an Msa into a list of AlignmentBlocks.
 * @param msa The msa to convert
 * @param row_indexes_to_caps The Caps for each sequence in the MSA
 * @param alignment_blocks The list to add the alignment blocks to
 */
void create_alignment_blocks(Msa *msa, Cap **row_indexes_to_caps, stList *alignment_blocks) {
    int64_t i=0; // The left most index of the current block
    bool rows_in_block[msa->seq_no]; // An array of bools used to indicate which sequences are present in a block
    int64_t seq_indexes[msa->seq_no]; // The start offsets of the current block
    for(int64_t k=0; k<msa->seq_no; k++) { // Initialize to zero
        seq_indexes[k] = 0;
    }
    int64_t sequences_in_block; // The number of sequences in the block

    //fprintf(stderr, "Start. Col no: %i\n", (int)msa->column_no);
    //msa_print(msa, stderr);

    // Walk through successive gapless blocks
    while(i < msa->column_no) {
        int64_t j = get_next_maximal_block_dimensions(msa, i, rows_in_block, &sequences_in_block);
        assert(j > i);
        assert(j <= msa->column_no);

        // Make the next alignment block
        if(sequences_in_block > 1) { // Only make a block if it contains two or more sequences
            stList_append(alignment_blocks, make_alignment_block(msa->seq_no, i, j - i, rows_in_block,
                                                                 seq_indexes, row_indexes_to_caps));
        }

        // Update the offsets in the sequences in the block, regardless of if we actually
        // created the block
        for(int64_t k=0; k<msa->seq_no; k++) {
            if(rows_in_block[k]) {
                seq_indexes[k] += j - i;
            }
        }

        i = j;
    }
    assert(i == msa->column_no);
}


int caps_comp_by_adjacency_length(const void *a, const void *b) {
    int64_t length1, length2;
    get_adjacency_string((Cap *)a, &length1, 0);
    get_adjacency_string((Cap *)b, &length2, 0);
    if (length1 != length2) {
        return length1 > length2 ? -1 : 1; // sort in descending order of length
    }
    // Equal lengths are common, and the row order decides the alignment, so
    // separate them by name instead of leaving it to whatever qsort does
    return cactusMisc_nameCompare(cap_getName((Cap *)a), cap_getName((Cap *)b));
}

/*
 * Returns end sequences sorted from longest to shortest
 */
void get_end_sequences(End *end, char **end_strings, int *end_string_lengths, int64_t *overlaps,
                       Cap **indices_to_caps, int64_t max_seq_length, int64_t mask_filter) {
    // Make inputs
    Cap *cap;
    End_InstanceIterator *capIterator = end_getInstanceIterator(end);
    int64_t j=0; // Index of the cap in the end's arrays

    // sorting the caps by length from longest to shortest (to make a consistent ordering, also POA seems to
    // create better alignments this way)
    stList *caps = stList_construct();
    while ((cap = end_getNext(capIterator)) != NULL) {
        if (cap_getSide(cap)) {
            cap = cap_getReverse(cap);
        }
        stList_append(caps, cap);
    }
    stList_sort(caps, caps_comp_by_adjacency_length); // sort by descending order of length
    end_destructInstanceIterator(capIterator);

    // Now create the actual end sequences
    for(int64_t i=0; i<stList_length(caps); i++) {
        Cap *cap = stList_get(caps, i);
        assert(!cap_getSide(cap));
        // Get the prefix of the adjacency string and its length and overlap with its reverse complement
        end_strings[j] = get_adjacency_string_and_overlap(cap, &(end_string_lengths[j]),
                                                          &(overlaps[j]), max_seq_length, mask_filter);

        // Populate the caps to end/row indices, and vice versa, data structures
        indices_to_caps[j++] = cap;
    }
    stList_destruct(caps); // cleanup
}

int64_t getMaxSequenceLength(End *end) {
    Cap *cap;
    End_InstanceIterator *capIterator = end_getInstanceIterator(end);
    int64_t max_length=0;
    while ((cap = end_getNext(capIterator)) != NULL) {
        if (cap_getSide(cap)) {
            cap = cap_getReverse(cap);
        }
        int64_t length;
        get_adjacency_string(cap, &length, 0);
        if(length > max_length) {
            max_length = length;
        }
    }
    end_destructInstanceIterator(capIterator);
    return max_length;
}

stList *make_flower_alignment_poa(Flower *flower, int64_t max_seq_length, int64_t window_size, int64_t mask_filter,
                                  int64_t max_prog_rows, double max_prog_length_diff, PoaParameters *poa_parameters) {
    End *dominantEnd = getDominantEnd(flower);
    int64_t seq_no = dominantEnd != NULL ? end_getInstanceNumber(dominantEnd) : -1;
    if(dominantEnd != NULL && getMaxSequenceLength(dominantEnd) < max_seq_length) {
        /*
         * If there is a single end that is connected to all adjacencies that are less than max_seq_length in length,
         * and the adjacencies include no self-aligned (self-loop) sequences
         * just use that alignment
         */
        // Make inputs
        char **end_strings = st_malloc(sizeof(char *) * seq_no);
        int *end_string_lengths = st_malloc(sizeof(int) * seq_no);
        int64_t overlaps[seq_no];
        Cap *indices_to_caps[seq_no];

        get_end_sequences(dominantEnd, end_strings, end_string_lengths, overlaps, indices_to_caps, max_seq_length, mask_filter);
        Msa *msa = msa_make_partial_order_alignment(end_strings, end_string_lengths, seq_no, window_size,
                                                    max_prog_rows, max_prog_length_diff, poa_parameters);

        //Now convert to set of alignment blocks
        stList *alignment_blocks = stList_construct3(0, (void (*)(void *))alignmentBlock_destruct);
        create_alignment_blocks(msa, indices_to_caps, alignment_blocks);

        // Cleanup
        msa_destruct(msa);

        return alignment_blocks;
    }

    // Arrays of ends and connecting the strings necessary to build the POA alignment
    int64_t end_no = flower_getEndNumber(flower); // The number of ends
    int64_t end_lengths[end_no]; // The number of strings incident with each end
    char **end_strings[end_no]; // The actual strings connecting the ends
    int *end_string_lengths[end_no]; // Length of the strings connecting the ends
    int64_t *right_end_indexes[end_no];  // For each string the index of the right end that it is connecting
    int64_t *right_end_row_indexes[end_no]; // For each string the index of the row of its reverse complement
    int64_t *overlaps[end_no]; // For each string the amount it suffix overlaps with its reverse complement

    // Data structures to translate between caps and sequences in above end arrays
    Cap **indices_to_caps[end_no]; // For each string the corresponding Cap
    stHash *caps_to_indices = stHash_construct2(NULL, free); // A hash of caps to their end and row indices

    // Fill out the end information for building the POA alignments arrays
    End *end;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    int64_t i=0; // Index of the end
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        // Initialize the various arrays for the end
        end_lengths[i] = end_getInstanceNumber(end); // The number of strings incident with the end
        end_strings[i] = st_malloc(sizeof(char *)*end_lengths[i]);
        end_string_lengths[i] = st_malloc(sizeof(int)*end_lengths[i]);
        right_end_indexes[i] = st_malloc(sizeof(int64_t)*end_lengths[i]);
        right_end_row_indexes[i] = st_malloc(sizeof(int64_t)*end_lengths[i]);
        indices_to_caps[i] = st_malloc(sizeof(Cap *)*end_lengths[i]);
        overlaps[i] = st_malloc(sizeof(int64_t)*end_lengths[i]);
        get_end_sequences(end, end_strings[i], end_string_lengths[i], overlaps[i], indices_to_caps[i],
                          max_seq_length, mask_filter);
        for(int64_t j=0; j<end_lengths[i]; j++) {
            stHash_insert(caps_to_indices, indices_to_caps[i][j], stIntTuple_construct2(i, j));
        }
        i++;
    }
    flower_destructEndIterator(endIterator);

    // Fill out the end / row indices for each cap
    endIterator = flower_getEndIterator(flower);
    i=0;
    while ((end = flower_getNextEnd(endIterator)) != NULL) { // Note, the use
        // of this iterator assumes that ends are iterated over in a stable order
        for(int64_t j=0; j<end_lengths[i]; j++) {
            Cap *cap = indices_to_caps[i][j];
            assert(!cap_getSide(cap));

            Cap *cap2 = cap_getAdjacency(cap);
            assert(cap2 != NULL);
            cap2 = cap_getReverse(cap2);
            assert(!cap_getSide(cap));
            assert(!cap_getSide(cap2));
            stIntTuple *k = stHash_search(caps_to_indices, cap2);
            assert(k != NULL);

            right_end_indexes[i][j] = stIntTuple_get(k, 0);
            right_end_row_indexes[i][j] = stIntTuple_get(k, 1);
        }
        i++;
    }
    flower_destructEndIterator(endIterator);

    // Now make the consistent MSAs
    Msa **msas = make_consistent_partial_order_alignments(end_no, end_lengths, end_strings, end_string_lengths,
                                                          right_end_indexes, right_end_row_indexes, overlaps, window_size,
                                                          max_prog_rows, max_prog_length_diff, poa_parameters);

    // Temp debug output
    //for(int64_t i=0; i<end_no; i++) {
    //    msa_print(msas[i], stderr);
    //}

    //Now convert to set of alignment blocks
    stList *alignment_blocks = stList_construct3(0, (void (*)(void *))alignmentBlock_destruct);
    for(int64_t i=0; i<end_no; i++) {
        create_alignment_blocks(msas[i], indices_to_caps[i], alignment_blocks);
    }

    // Cleanup
    for(int64_t i=0; i<end_no; i++) {
        msa_destruct(msas[i]);
        free(right_end_indexes[i]);
        free(right_end_row_indexes[i]);
        free(indices_to_caps[i]);
        free(overlaps[i]);
    }
    free(msas);
    stHash_destruct(caps_to_indices);

    // Temp debug output
    //for(int64_t i=0; i<stList_length(alignment_blocks); i++) {
    //    alignmentBlock_print(stList_get(alignment_blocks, i), stderr);
    //}

    return alignment_blocks;
}

/*
 * The following is used for converting the alignment blocks into pinches consumed by the CAF code.
 */

/**
 * Iterator over the list of alignment blocks used to get stPinches in succession.
 */
typedef struct _alignmentBlockIterator {
    stList *alignment_blocks; // The list of alignment blocks
    int64_t i; // Index of the iterator into the alignment_blocks
    AlignmentBlock *current_block; // The current block being considered
} AlignmentBlockIterator;

AlignmentBlockIterator *alignmentBlockIterator_construct(stList *alignment_blocks) {
    AlignmentBlockIterator *alignmentBlockIterator = st_calloc(1, sizeof(AlignmentBlockIterator));
    alignmentBlockIterator->alignment_blocks = alignment_blocks;
    return alignmentBlockIterator;
}

void alignmentBlockIterator_destruct(AlignmentBlockIterator *it) {
    stList_length(it->alignment_blocks);
    free(it);
}

AlignmentBlockIterator *alignmentBlockIterator_start(AlignmentBlockIterator *it) {
    it->i = 0;
    it->current_block = NULL;
    return it;
}

stPinch *alignmentBlockIterator_get_next(AlignmentBlockIterator *it, stPinch *pinchToFillOut) {
    // If there is no current alignment block or the alignment block contains no further pinches
    if(it->current_block == NULL || it->current_block->next == NULL) {
        if(it->i >= stList_length(it->alignment_blocks)) { // We are done
            return NULL;
        }
        it->current_block = stList_get(it->alignment_blocks, it->i++);
    }
    assert(it->current_block->next != NULL); // All alignment blocks should contain at least two sequences

    AlignmentBlock *b = it->current_block;
    assert(b->position >= 0);
    assert(b->next->position >= 0);
    assert(b->length > 0);
    stPinch_fillOut(pinchToFillOut, b->subsequenceIdentifier, b->next->subsequenceIdentifier,
                    b->position, b->next->position, b->length, b->strand == b->next->strand);

    it->current_block = b->next; // Shift to the next sequence to ready the next pinch

    return pinchToFillOut;
}

stPinchIterator *stPinchIterator_constructFromAlignedBlocks(stList *alignment_blocks) {
    stPinchIterator *pinchIterator = st_calloc(1, sizeof(stPinchIterator));
    pinchIterator->alignmentArg = alignmentBlockIterator_construct(alignment_blocks);
    pinchIterator->getNextAlignment = (stPinch *(*)(void *, stPinch *)) alignmentBlockIterator_get_next;
    pinchIterator->destructAlignmentArg = (void(*)(void *)) alignmentBlockIterator_destruct;
    pinchIterator->startAlignmentStack = (void *(*)(void *)) alignmentBlockIterator_start;

    return pinchIterator;
}
