/*
 * This is designed as a drop-in replacement for the multiple aligner from pecan that gets used in the end aligner. 
 * The idea is that this will scale better for larger numbers of input samples, which seems to blow up the memory
 * in pecan.  
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "endAligner.h"
#include "multipleAligner.h"
#include "abpoa.h"
#include "poaAligner.h"
#include "adjacencySequences.h"
#include "pairwiseAligner.h"
#include "sonLib.h"
#include "stPosetAlignment.h"
#include "pairwiseAligner.h"
#include <stdlib.h>
#include <math.h>
#include "stGraph.h"
#include <inttypes.h>

// char <--> uint8_t conversion copied over from abPOA example
// AaCcGgTtNn ==> 0,1,2,3,4
static uint8_t nst_nt4_table[256] = {
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

static inline char toBase(uint8_t n) {
    return "ACGTN-"[n];
}

static inline uint8_t toByte(char c) {
    return nst_nt4_table[(int)c];
}

/** Copy the fragment sequences into the matrix format that abpoa expects, using 
 * offsets and windowSize to subset as needed.  updates bseqs and chunk_lens
 */
static void toPoaInput(stList* seqFrags, int64_t n_seqs, int64_t* seq_lens, int64_t* seq_offsets,
                       int64_t windowSize, uint8_t** bseqs, int* chunk_lens) {
    for (int64_t i = 0; i < n_seqs; ++i) {
        chunk_lens[i] = 0;
        SeqFrag* frag = (SeqFrag*)stList_get(seqFrags, i);
        for (int64_t j = seq_offsets[i]; j < seq_lens[i] && chunk_lens[i] < windowSize; ++j, ++chunk_lens[i]) {
            // todo: support iupac characters?
            bseqs[i][chunk_lens[i]] = toByte(frag->seq[j]);
        }
    }
}

MultipleAlignment *makePartialOrderAlignment(StateMachine *sM, stList *seqFrags,
                                             float matchGamma,
                                             PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters,
                                             int64_t windowSize) {
    
    // transorm the input for abpoa, using its example.c as a guide
    int64_t n_seqs = stList_length(seqFrags);
    
    // lengths of each input fragment
    int64_t *seq_lens = (int64_t*)st_malloc(sizeof(int64_t) * n_seqs);
    int64_t bases_remaining = 0;
    
    // keep track of current lengths
    int *chunk_lens = (int*)st_malloc(sizeof(int) * n_seqs);
    // keep track of current offsets
    int64_t* seq_offsets = (int64_t*)st_calloc(n_seqs, sizeof(int64_t));
    // keep track of empty chunks
    bool* empty_seqs = (bool*)st_calloc(n_seqs, sizeof(bool));

    // initialize the poa input matrix
    uint8_t **bseqs = (uint8_t**)st_malloc(sizeof(uint8_t*) * n_seqs);
    for (int64_t i = 0; i < n_seqs; ++i) {
        SeqFrag* frag = stList_get(seqFrags, i);
        seq_lens[i] = frag->length;
        int64_t chunk_len = windowSize < seq_lens[i] ? windowSize : seq_lens[i];
        bseqs[i] = (uint8_t*)st_malloc(sizeof(uint8_t) * (chunk_len));
        bases_remaining += seq_lens[i];
    }
    
    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    // todo: try to convert pecan parameters into the abpoa scores

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

    // todo: score
    MultipleAlignment *mA = st_calloc(1, sizeof(MultipleAlignment));
    mA->alignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    // todo: refactor the multiple alignment object, as these should be optional
    mA->chosenPairwiseAlignments = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    mA->columns = stSet_construct(); //makeColumns(seqFrags);

    int64_t prev_bases_remaining = bases_remaining;
    while (bases_remaining > 0) {

        // load up to windowSize of each sequence into the input matrix for poa
        toPoaInput(seqFrags, n_seqs, seq_lens, seq_offsets, windowSize, bseqs, chunk_lens);        
        
        // variables to store resulting MSA;
        uint8_t **msa_seq; int msa_l=0;

        // poa can't handle empty sequences.  this is a hack to get around that
        int emptyCount = 0;
        for (int64_t i = 0; i < n_seqs; ++i) {
            if (chunk_lens[i] == 0) {
                empty_seqs[i] = true;
                chunk_lens[i] = 1;
                bseqs[i][0] = toByte('N');
                ++emptyCount;
            } else {
                empty_seqs[i] = false;
            }
        }

        // perform abpoa-msa
        abpoa_msa(ab, abpt, n_seqs, NULL, chunk_lens, bseqs, NULL, NULL, NULL, NULL, NULL, &msa_seq, &msa_l);

        // mask out empty sequences that were phonied in as Ns above
        for (int64_t i = 0; i < n_seqs && emptyCount > 0; ++i) {
            if (empty_seqs[i] == true) {
                for (int j = 0; j < msa_l; ++j) {
                    if (toBase(msa_seq[i][j]) != '-') {
                        assert(toBase(msa_seq[i][j]) == 'N');
                        msa_seq[i][j] = toByte('-');
                        --chunk_lens[i];
                        assert(chunk_lens[i] == 0);
                        --emptyCount;
                        break;
                    }
                }
            }
        }
        assert(emptyCount == 0);

        // todo: lots of room to be smarter about finding the next anchor here by, for example, unaligning any
        // bad-looking columns at the end of the matrix (and updating the various offset structures) or
        // just using a sliding window with overlap, and ignoring alignments near the end
        
        // todo: do we want to use global alignment each time? 

        // add the msa matrix as tuples to our aligned pairs list
        // note that seq_offsets is modified in place here, updating the position along the original sequences
        poaMatrixToAlignedPairs(msa_seq, n_seqs, msa_l, 1.0, seqFrags, seq_offsets, mA->alignedPairs);

        // clean up
        for (int64_t i = 0; i < n_seqs; ++i) {
            free(msa_seq[i]);
        }
        free(msa_seq);

        // remember how much we aligned this round
        for (int64_t i = 0; i < n_seqs; ++i) {
            bases_remaining -= chunk_lens[i];
        }
        // and use for sanity check        
        assert(prev_bases_remaining > bases_remaining && bases_remaining >= 0);

        if (bases_remaining > 0) {
            // reset graph before re-use
            abpoa_reset_graph(ab, abpt, seq_lens[0]); 
        }

        //used only for sanity check
        prev_bases_remaining = bases_remaining; 
    }

    for (int64_t i = 0; i < n_seqs; ++i) {
        free(bseqs[i]);
    }
    free(bseqs);
    free(seq_lens);
    free(chunk_lens);
    free(seq_offsets);
    free(empty_seqs);
    abpoa_free(ab, abpt);
    abpoa_free_para(abpt); 

    // in debug mode, cactus uses the dreaded -Wall -Werror combo.  This line is a hack to allow compilition with these flags
    if (false) SIMDMalloc(0, 0);
    
    return mA;
}

void poaMatrixToAlignedPairs(uint8_t** msaSeq, int64_t numSeqs, int64_t msaWidth, int64_t score, stList* seqFrags, int64_t* offsets, stList* outPairs) {
/*
    fprintf(stdout, ">Multiple_sequence_alignment\n");
    for (int64_t i = 0; i < numSeqs; ++i) {
        for (int64_t j = 0; j < msaWidth; ++j) {
            fprintf(stdout, "%c", toBase(msaSeq[i][j]));
        }
        fprintf(stdout, "\n");
    }
*/

    for (int64_t column = 0; column < msaWidth; ++column) { // For each column
        int64_t anchor = -1, anchorCoordinate;
        for (int64_t seqIdx = 0; seqIdx < numSeqs; ++seqIdx) {
            if(toBase(msaSeq[seqIdx][column]) != '-') { // If is not a gap
                int64_t seqCoordinate = offsets[seqIdx]++;
                if (toBase(msaSeq[seqIdx][column]) != 'N') { // If is not an ambiguous base
                    if (anchor == -1) { // Set as the anchoring base
                        anchor = seqIdx;
                        anchorCoordinate = seqCoordinate;
                    } else { // Otherwise make an aligned pair between the anchor and the base
                        stList_append(outPairs, stIntTuple_construct5(
                                score, anchor, anchorCoordinate,
                                seqIdx, seqCoordinate));
                    }
                }

            }
        }
    }

}

