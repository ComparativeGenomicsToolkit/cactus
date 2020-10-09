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

MultipleAlignment *makePartialOrderAlignment(StateMachine *sM, stList *seqFrags,
                                             int64_t spanningTrees, int64_t maxPairsToConsider,
                                             bool useProgressiveMerging,
                                             float matchGamma,
                                             PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {

    // transorm the input for abpoa, using its example.c as a guide
    int64_t n_seqs = stList_length(seqFrags);
    
    // collect sequence length, trasform ACGT to 0123
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);

    stListIterator *it = stList_getIterator(seqFrags);
    SeqFrag* frag;
    int i = 0;
    int j = 0;    
    while ((frag = stList_getNext(it)) != NULL) {
        fprintf(stderr, "frag = %s\n", frag->seq);
        seq_lens[i] = frag->length;
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * frag->length);
        for (j = 0; j < frag->length; ++j) {
            // todo: support iupac characters!
            bseqs[i][j] = nst_nt4_table[(int)frag->seq[j]];
        }
        ++i;
    }

    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

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

    // variables to store result
    uint8_t **cons_seq; int **cons_cov, *cons_l, cons_n=0;
    uint8_t **msa_seq; int msa_l=0;

    // perform abpoa-msa
    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);

    // print
    fprintf(stdout, "=== output to variables ===\n");
    for (i = 0; i < cons_n; ++i) {
        fprintf(stdout, ">Consensus_sequence\n");
        for (j = 0; j < cons_l[i]; ++j)
            fprintf(stdout, "%c", "ACGTN"[cons_seq[i][j]]);
        fprintf(stdout, "\n");
    }

    // todo: score
    MultipleAlignment* mA = poaMatrixToMultipleAlignment(msa_seq, n_seqs, msa_l, 0, seqFrags);
        
    // clean up
    if (cons_n) {
        for (i = 0; i < cons_n; ++i) {
            free(cons_seq[i]); free(cons_cov[i]);
        }
        free(cons_seq); free(cons_cov); free(cons_l);
    }
    if (msa_l) {
        for (i = 0; i < n_seqs; ++i) {
            free(msa_seq[i]);
        }
        free(msa_seq);
    }

    for (i = 0; i < n_seqs; ++i) {
        free(bseqs[i]);
    }
    
    free(bseqs);
    free(seq_lens);
    abpoa_free(ab, abpt);
    abpoa_free_para(abpt); 
    
    return mA;
}

static inline char toBase(uint8_t n) {
    return "ACGTN-"[n];
}

MultipleAlignment *poaMatrixToMultipleAlignment(uint8_t** msaSeq, int numSeqs, int msaWidth, int score, stList* seqFrags) {

    fprintf(stdout, ">Multiple_sequence_alignment\n");
    for (int i = 0; i < numSeqs; ++i) {
        for (int j = 0; j < msaWidth; ++j) {
            fprintf(stdout, "%c", toBase(msaSeq[i][j]));
        }
        fprintf(stdout, "\n");
    }

    // this is a big waste of memory
    // todo: can we go lower into the poa api to avoid, or perhaps just use an interval index
    int32_t** offsets = st_calloc(numSeqs, sizeof(int32_t*));
    for (int i = 0; i < numSeqs; ++i) {
        offsets[i] = st_calloc(msaWidth, sizeof(int32_t));
    }
    for (int i = 0; i < numSeqs; ++i) {
        int32_t offset = 0;
        for (int j = 0; j < msaWidth; ++j) {
            offsets[i][j] = offset;
            if (toBase(msaSeq[i][j]) != '-') {
                ++offset;
            }
        }
        int64_t seqLength = ((SeqFrag*)stList_get(seqFrags, i))->length;
        assert(offset >= seqLength);
    }
    
    MultipleAlignment *mA = st_calloc(1, sizeof(MultipleAlignment));

    //pairwise alignment pairs, with sequence indices
    mA->alignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct); //pairwise alignment pairs, with sequence indices

    for (int column = 0; column < msaWidth; ++column) {
        int anchor = -1;
        for (int seqIdx = 0; seqIdx < numSeqs; ++seqIdx) {
            if (toBase(msaSeq[seqIdx][column]) != 'N' && toBase(msaSeq[seqIdx][column]) != '-') {
                if (anchor == -1) {
                    anchor = seqIdx;
                } else if (anchor >= 0) {
                    stList_append(mA->alignedPairs, stIntTuple_construct5(
                                      score,
                                      anchor,
                                      offsets[anchor][column],
                                      seqIdx,
                                      offsets[seqIdx][column]));
                }
            }
        }
    }
    return mA;
    
}

