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

char msa_to_base(uint8_t n) {
    return "ACGTN-"[n];
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
                column_scores[i] = column_scores[i]-1 > 0 ? column_scores[i]-1 : 0;
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
    assert(overlap > 0); // Otherwise the overlap must be positivie

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

    // Not walk through each possible cut point within the overlap
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

Msa *msa_make_partial_order_alignment(char **seqs, int *seq_lens, int64_t seq_no, int64_t window_size) {

    assert(seq_no > 0);
    
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

    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    // todo: support including modifying abpoa params
    // alignment parameters
    // abpt->align_mode = 0; // 0:global alignment, 1:extension
    // abpt->match = 2;      // match score
    // abpt->mismatch = 4;   // mismatch penalty
    // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
    // abpt->gap_open1 = 4;  // gap open penalty #1
    // abpt->gap_ext1 = 2;   // gap extension penalty #1
    // abpt->gap_open2 = 24; // gap open penalty #2
    // abpt->gap_ext2 = 1;   // gap extension penalty #2
                             // gap_penalty = min{gap_open1 + gap_len * gap_ext1, gap_open2 + gap_len * gap_ext2}
    // abpt->bw = 10;        // extra band used in adaptive banded DP
    // abpt->bf = 0.01; 
     
    // output options
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 0; // generate consensus sequence, set 0 to disable

    abpoa_post_set_para(abpt);

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

        // perform abpoa-msa
        abpoa_msa(ab, abpt, msa->seq_no, NULL, msa->seq_lens, bseqs, NULL, NULL, NULL, NULL, NULL,
                  &(msa->msa_seq), &(msa->column_no));

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
        }

        // add the msa to our list
        stList_append(msa_windows, msa);
        
        // sanity check        
        assert(prev_bases_remaining > bases_remaining && bases_remaining >= 0);

        if (bases_remaining > 0) {
            // reset graph before re-use
            abpoa_reset_graph(ab, abpt, msa->seq_lens[0]); 
        }

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

    abpoa_free(ab, abpt);
    abpoa_free_para(abpt);

    // in debug mode, cactus uses the dreaded -Wall -Werror combo.  This line is a hack to allow compilation with these flags
    if (false) SIMDMalloc(0, 0);

    return output_msa;
}

Msa **make_consistent_partial_order_alignments(int64_t end_no, int64_t *end_lengths, char ***end_strings,
        int **end_string_lengths, int64_t **right_end_indexes, int64_t **right_end_row_indexes, int64_t **overlaps,
        int64_t window_size) {
    // Calculate the initial, potentially inconsistent msas and column scores for each msa
    float *column_scores[end_no];
    Msa **msas = st_malloc(sizeof(Msa *) * end_no);
    for(int64_t i=0; i<end_no; i++) {
        msas[i] = msa_make_partial_order_alignment(end_strings[i], end_string_lengths[i], end_lengths[i], window_size);
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
                        right_end_row_index, msas[right_end_index], column_scores[right_end_index], overlaps[i][j]) ;
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

char *get_adjacency_string(Cap *cap, int *length) {
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
        return sequence_getString(sequence, cap_getCoordinate(cap) + 1, *length, 1);
    } else {
        assert(cap_getCoordinate(cap) > cap_getCoordinate(cap2));
        *length = cap_getCoordinate(cap) - cap_getCoordinate(cap2) - 1;
        assert(*length >= 0);
        return sequence_getString(sequence, cap_getCoordinate(cap2) + 1, *length, 0);
    }
}

/**
 * Used to get a prefix of a given adjacency sequence.
 * @param seq_length
 * @param length
 * @param overlap
 * @param max_seq_length
 * @return
 */
char *get_adjacency_string_and_overlap(Cap *cap, int *length, int64_t *overlap, int64_t max_seq_length) {
    // Get the complete adjacency string
    int seq_length;
    char *adjacency_string = get_adjacency_string(cap, &seq_length);
    assert(seq_length >= 0);

    // Calculate the length of the prefix up to max_seq_length
    *length = seq_length > max_seq_length ? max_seq_length : seq_length;
    assert(*length >= 0);

    // Cleanup the string
    adjacency_string[*length] = '\0'; // Terminate the string at the given length
    char *c = stString_copy(adjacency_string);
    free(adjacency_string);
    adjacency_string = c;

    // Calculate the overlap with the reverse complement
    if(2*(*length) > seq_length) { // There is overlap
        *overlap = 2*(*length) - seq_length;
        assert(*overlap >= 0);
    }
    else { // There is no overlap
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

stList *make_flower_alignment_poa(Flower *flower, int64_t max_seq_length, int64_t window_size) {
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

        // Now get each string incident with the end
        Cap *cap;
        End_InstanceIterator *capIterator = end_getInstanceIterator(end);
        int64_t j=0; // Index of the cap in the end's arrays
        while ((cap = end_getNext(capIterator)) != NULL) {
            assert(j < end_lengths[i]);
            // Ensure we have the cap in the correct orientation
            if (cap_getSide(cap)) {
                cap = cap_getReverse(cap);
            }
            // Get the prefix of the adjacency string and its length and overlap with its reverse complement
            end_strings[i][j] = get_adjacency_string_and_overlap(cap, &(end_string_lengths[i][j]),
                                                                 &(overlaps[i][j]), max_seq_length);

            // Populate the caps to end/row indices, and vice versa, data structures
            indices_to_caps[i][j] = cap;
            stHash_insert(caps_to_indices, cap, stIntTuple_construct2(i, j));

            j++;
        }
        end_destructInstanceIterator(capIterator);
        assert(end_lengths[i] == j);
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
        assert(end_lengths[i] == j);
        i++;
    }
    flower_destructEndIterator(endIterator);

    // Now make the consistent MSAs
    Msa **msas = make_consistent_partial_order_alignments(end_no, end_lengths, end_strings, end_string_lengths,
                                                          right_end_indexes, right_end_row_indexes, overlaps, window_size);

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

stPinch *alignmentBlockIterator_get_next(AlignmentBlockIterator *it) {
    // If there is no current alignment block or the alignment block contains no further pinches
    if(it->current_block == NULL || it->current_block->next == NULL) {
        if(it->i >= stList_length(it->alignment_blocks)) { // We are done
            return NULL;
        }
        it->current_block = stList_get(it->alignment_blocks, it->i++);
    }
    assert(it->current_block->next != NULL); // All alignment blocks should contain at least two sequences

    AlignmentBlock *b = it->current_block;
    static stPinch pinch;
    assert(b->position >= 0);
    assert(b->next->position >= 0);
    assert(b->length > 0);
    stPinch_fillOut(&pinch, b->subsequenceIdentifier, b->next->subsequenceIdentifier,
                    b->position, b->next->position, b->length, b->strand == b->next->strand);

    it->current_block = b->next; // Shift to the next sequence to ready the next pinch

    return &pinch;
}

stPinchIterator *stPinchIterator_constructFromAlignedBlocks(stList *alignment_blocks) {
    stPinchIterator *pinchIterator = st_calloc(1, sizeof(stPinchIterator));
    pinchIterator->alignmentArg = alignmentBlockIterator_construct(alignment_blocks);
    pinchIterator->getNextAlignment = (stPinch *(*)(void *)) alignmentBlockIterator_get_next;
    pinchIterator->destructAlignmentArg = (void(*)(void *)) alignmentBlockIterator_destruct;
    pinchIterator->startAlignmentStack = (void *(*)(void *)) alignmentBlockIterator_start;

    return pinchIterator;
}
