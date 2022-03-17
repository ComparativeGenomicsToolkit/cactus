/**
 * This is designed as a drop-in replacement for the bar aligner, using the abpoa multiple sequence aligner.
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "abpoa.h"
#include "poaBarAligner.h"
#include "flowerAligner.h"

#include <stdio.h>
#include <ctype.h>

// FOR DEBUGGING ONLY: Specify directory where abpoa inputs get dumped
//#define CACTUS_ABPOA_MSA_DUMP_DIR "/home/hickey/dev/cactus/dump"
// FOR DEBUGGING ONLY: Run abpoa from command line instead of via API (only works with CACTUS_ABPOA_MSA_DUMP_DIR defined)
//#define CACTUS_ABPOA_FROM_COMMAND_LINE

// OpenMP
//#if defined(_OPENMP)
//#include <omp.h>
//#endif

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

#ifdef CACTUS_ABPOA_MSA_DUMP_DIR
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
#endif

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

Msa *msa_make_partial_order_alignment(char **seqs, int *seq_lens, int64_t seq_no, int64_t window_size,
                                      abpoa_para_t *poa_parameters) {

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

        // init abpoa
        abpoa_t *ab = abpoa_init();
        abpoa_para_t *abpt = copy_abpoa_params(poa_parameters);
        
#ifdef CACTUS_ABPOA_MSA_DUMP_DIR
        // dump the input to file
        char abpoa_input_path[1024], abpoa_matrix_path[1024], abpoa_command_path[1024], abpoa_output_path[1024];
        sprintf(abpoa_input_path, "%s/ap_in_%ld.fa", CACTUS_ABPOA_MSA_DUMP_DIR, (int64_t)msa);
        sprintf(abpoa_matrix_path, "%s.mat", abpoa_input_path);
        sprintf(abpoa_command_path, "%s.cmd", abpoa_input_path);
        sprintf(abpoa_output_path, "%s.out", abpoa_input_path);
        char* abpoa_command_line = dump_abpoa_input(msa, abpt, bseqs,
                                                    abpoa_input_path, abpoa_matrix_path, abpoa_command_path, abpoa_output_path);
#endif

#ifdef CACTUS_ABPOA_FROM_COMMAND_LINE
        // run abpoa from the command line
        abpoa_msa_from_command_line(abpoa_command_line, abpoa_output_path, &(msa->msa_seq), &(msa->column_no));

        int test_cols = 0;
        uint8_t** test_msa = NULL;
        abpoa_msa(ab, abpt, msa->seq_no, NULL, msa->seq_lens, bseqs, NULL);
        // abpoa's interface has changed a bit -- instead of passing in pointers to the results, they
        // end up in the ab->abc struct -- we extract them here
        msa->msa_seq = ab->abc->msa_base;
        ab->abc->msa_base = NULL;
        msa->column_no = ab->abc->msa_len;

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
        abpoa_msa(ab, abpt, msa->seq_no, NULL, msa->seq_lens, bseqs, NULL);
        // abpoa's interface has changed a bit -- instead of passing in pointers to the results, they
        // end up in the ab->abc struct -- we extract them here
        msa->msa_seq = ab->abc->msa_base;
        ab->abc->msa_base = NULL;
        msa->column_no = ab->abc->msa_len;
#endif

#ifdef CACTUS_ABPOA_MSA_DUMP_DIR
        // we got this far without crashing, so delete the dumped file (they can really pile up otherwise)
        remove(abpoa_input_path);
        remove(abpoa_matrix_path);
        remove(abpoa_command_path);
        remove(abpoa_output_path);
        free(abpoa_command_line);
#endif

        // free abpoa
        abpoa_free(ab);
        abpoa_free_para(abpt);

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

Msa **make_consistent_partial_order_alignments(int64_t end_no, int64_t *end_lengths, char ***end_strings,
        int **end_string_lengths, int64_t **right_end_indexes, int64_t **right_end_row_indexes, int64_t **overlaps,
        int64_t window_size, abpoa_para_t *poa_parameters) {
    // Calculate the initial, potentially inconsistent msas and column scores for each msa
    float *column_scores[end_no];
    Msa **msas = st_malloc(sizeof(Msa *) * end_no);
//#if defined(_OPENMP)
//#pragma omp parallel for schedule(dynamic)
//#endif
    for(int64_t i=0; i<end_no; i++) {
        msas[i] = msa_make_partial_order_alignment(end_strings[i], end_string_lengths[i], end_lengths[i], window_size,
                                                   poa_parameters);
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

char *get_adjacency_string(Cap *cap, int *length, bool return_string) {
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

/**
 * Used to find where a run of masked (hard or soft) of at least mask_filter bases starts
 * @param seq : The string
 * @param seq_length : The length of the string
 * @param length : The maximum length we want to search in
 * @param reversed : If true, scan from the end of the string
 * @param mask_filter : Cut a string as soon as we hit more than this many hard or softmasked bases (cut is before first masked base)
 * @return length of the filtered string
 */
static int get_unmasked_length(char* seq, int64_t seq_length, int64_t length, bool reversed, int64_t mask_filter) {
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
                    return (int)run_start;
                }
            } else {
                run_start = -1;
            }
        }
    }
    return (int)length;
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
    // Get the complete adjacency string
    int seq_length;
    char *adjacency_string = get_adjacency_string(cap, &seq_length, 1);
    assert(seq_length >= 0);

    // Calculate the length of the prefix up to max_seq_length
    *length = seq_length > max_seq_length ? max_seq_length : seq_length;
    assert(*length >= 0);
    int length_backward = *length;

    if (mask_filter >= 0) {
        // apply the mask filter on the forward strand
        *length = get_unmasked_length(adjacency_string, seq_length, *length, false, mask_filter);
        length_backward = get_unmasked_length(adjacency_string, seq_length, *length, true, mask_filter);
    }

    // Cleanup the string
    adjacency_string[*length] = '\0'; // Terminate the string at the given length
    char *c = stString_copy(adjacency_string);
    free(adjacency_string);
    adjacency_string = c;

    // Calculate the overlap with the reverse complement
    if (*length + length_backward > seq_length) { // There is overlap
        *overlap = *length + length_backward - seq_length;
        assert(*overlap >= 0);
    } else { // There is no overlap
        *overlap = 0;
    }

    return adjacency_string;
}

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

void get_end_sequences(End *end, char **end_strings, int *end_string_lengths, int64_t *overlaps,
                       Cap **indices_to_caps, int64_t max_seq_length, int64_t mask_filter) {
    // Make inputs
    Cap *cap;
    End_InstanceIterator *capIterator = end_getInstanceIterator(end);
    int64_t j=0; // Index of the cap in the end's arrays
    while ((cap = end_getNext(capIterator)) != NULL) {
        // Ensure we have the cap in the correct orientation
        if (cap_getSide(cap)) {
            cap = cap_getReverse(cap);
        }
        // Get the prefix of the adjacency string and its length and overlap with its reverse complement
        end_strings[j] = get_adjacency_string_and_overlap(cap, &(end_string_lengths[j]),
                                                          &(overlaps[j]), max_seq_length, mask_filter);

        // Populate the caps to end/row indices, and vice versa, data structures
        indices_to_caps[j] = cap;

        j++;
    }
    end_destructInstanceIterator(capIterator);
}

int64_t getMaxSequenceLength(End *end) {
    Cap *cap;
    End_InstanceIterator *capIterator = end_getInstanceIterator(end);
    int64_t max_length=0;
    while ((cap = end_getNext(capIterator)) != NULL) {
        if (cap_getSide(cap)) {
            cap = cap_getReverse(cap);
        }
        int length;
        get_adjacency_string(cap, &length, 0);
        if(length > max_length) {
            max_length = length;
        }
    }
    end_destructInstanceIterator(capIterator);
    return max_length;
}

stList *make_flower_alignment_poa(Flower *flower, int64_t max_seq_length, int64_t window_size, int64_t mask_filter,
                                  abpoa_para_t * poa_parameters) {
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
        Msa *msa = msa_make_partial_order_alignment(end_strings, end_string_lengths, seq_no, window_size, poa_parameters);

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
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        Cap *cap;
        End_InstanceIterator *capIterator = end_getInstanceIterator(end);
        int64_t j=0;
        while ((cap = end_getNext(capIterator)) != NULL) {
            if (cap_getSide(cap)) {
                cap = cap_getReverse(cap);
            }
            Cap *cap2 = cap_getAdjacency(cap);
            assert(cap2 != NULL);
            cap2 = cap_getReverse(cap2);
            assert(!cap_getSide(cap));
            assert(!cap_getSide(cap2));
            stIntTuple *k = stHash_search(caps_to_indices, cap2);
            assert(k != NULL);

            right_end_indexes[i][j] = stIntTuple_get(k, 0);
            right_end_row_indexes[i][j] = stIntTuple_get(k, 1);

            j++;
        }
        end_destructInstanceIterator(capIterator);
        i++;
    }
    flower_destructEndIterator(endIterator);

    // Now make the consistent MSAs
    Msa **msas = make_consistent_partial_order_alignments(end_no, end_lengths, end_strings, end_string_lengths,
                                                          right_end_indexes, right_end_row_indexes, overlaps, window_size,
                                                          poa_parameters);

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
