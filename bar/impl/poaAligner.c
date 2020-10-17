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

static inline char toBase(uint8_t n) {
    return "ACGTN-"[n];
}

static inline uint8_t toByte(char c) {
    return nst_nt4_table[(int)c];
}

MultipleAlignment *makePartialOrderAlignment(StateMachine *sM, stList *seqFrags,
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
        seq_lens[i] = frag->length;
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * frag->length);
        for (j = 0; j < frag->length; ++j) {
            // todo: support iupac characters?
            bseqs[i][j] = toByte(frag->seq[j]);
        }
        ++i;
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

    // variables to store resulting MSA;
    uint8_t **msa_seq; int msa_l=0;

    // perform abpoa-msa
    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, NULL, NULL, NULL, NULL, &msa_seq, &msa_l);

    // todo: score
    bool allPairs = 1;
    MultipleAlignment *mA = st_calloc(1, sizeof(MultipleAlignment));
    mA->alignedPairs = poaMatrixToAlignedPairs(msa_seq, n_seqs, msa_l, 1.0, seqFrags, allPairs);
    // todo: refactor the multiple alignment object, as these should be optional
    mA->chosenPairwiseAlignments = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    if(allPairs) {
        for (int64_t i = 0; i < n_seqs; i++) {
            for (int64_t j =i+1; j < n_seqs; j++) {
                stList_append(mA->chosenPairwiseAlignments, stIntTuple_construct2(i, j));
            }
        }
    }
    mA->columns = stSet_construct(); //makeColumns(seqFrags);
        
    // clean up
    for (i = 0; i < n_seqs; ++i) {
      free(msa_seq[i]);
    }
    free(msa_seq);

    for (i = 0; i < n_seqs; ++i) {
        free(bseqs[i]);
    }
    
    free(bseqs);
    free(seq_lens);
    abpoa_free(ab, abpt);
    abpoa_free_para(abpt); 

    // in debug mode, cactus uses the dreaded -Wall -Werror combo.  This line is a hack to allow compilition with these flags
    if (false) SIMDMalloc(0, 0);
    
    return mA;
}

stList *poaMatrixToAlignedPairs(uint8_t** msaSeq, int numSeqs, int msaWidth, int score, stList* seqFrags, bool allPairs) {
/*
    fprintf(stdout, ">Multiple_sequence_alignment\n");
    for (int i = 0; i < numSeqs; ++i) {
        for (int j = 0; j < msaWidth; ++j) {
            fprintf(stdout, "%c", toBase(msaSeq[i][j]));
        }
        fprintf(stdout, "\n");
    }
*/

    stList *alignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    int64_t* offsets = st_calloc(numSeqs, sizeof(int64_t));

    if(allPairs) {
        for (int column = 0; column < msaWidth; ++column) { // For each column
            for (int seqIdx = 0; seqIdx < numSeqs; ++seqIdx) {
                if (toBase(msaSeq[seqIdx][column]) != '-') { // If is not a gap
                    int64_t seqCoordinate = offsets[seqIdx]++;
                    for (int seqIdx2 = seqIdx; seqIdx2 < numSeqs; ++seqIdx2) {
                        if (toBase(msaSeq[seqIdx2][column]) != '-') { // If is not a gap
                            stList_append(alignedPairs, stIntTuple_construct5(
                                    score, seqIdx, seqCoordinate,
                                    seqIdx2, offsets[seqIdx2]));
                        }
                    }
                }
            }
        }
    }
    else {
        for (int column = 0; column < msaWidth; ++column) { // For each column
            int64_t anchor = -1, anchorCoordinate;
            for (int seqIdx = 0; seqIdx < numSeqs; ++seqIdx) {
                if (toBase(msaSeq[seqIdx][column]) != '-') { // If is not a gap
                    int64_t seqCoordinate = offsets[seqIdx]++;
                    if (toBase(msaSeq[seqIdx][column]) != 'N') { // If is not an ambiguous base
                        if (anchor == -1) { // Set as the anchoring base
                            anchor = seqIdx;
                            anchorCoordinate = seqCoordinate;
                        } else { // Otherwise make an aligned pair between the anchor and the base
                            stList_append(alignedPairs, stIntTuple_construct5(
                                    score, anchor, anchorCoordinate,
                                    seqIdx, seqCoordinate));
                        }
                    }
                }
            }
        }
    }


    free(offsets);

    return alignedPairs;
}

