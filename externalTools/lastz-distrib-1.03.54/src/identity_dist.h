//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: identity_dist.h
//
//----------

#ifndef identity_dist_H			// (prevent multiple inclusion)
#define identity_dist_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "sequences.h"			// sequence stuff
#include "segment.h"			// segment table management stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef identity_dist_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef identity_dist_owner
int identity_dist_dbgShowIdentity = false;
#else
global int identity_dist_dbgShowIdentity;
#endif

//----------
//
// data structures and types
//
//----------

#define numIdentityBins       1000
#define identityBinFormat     "%.3f"
#define identityBinLongFormat "%.4f"

// We map identity fraction f to bin b = floor (numBins*f + .5).  So the bins
// run from
//
//    b - 1/2        b + 1/2
//    ------- <= f < -------
//    numBins        numBins
// 
//
// identity_bin(n,d) = numBins * (n/d), rounded to integer
 
#define identity_bin(numer,denom) \
          ((u32)((2*((u64)numer)*(numIdentityBins)+((u64)denom))/(2*(denom))))
#define bin_to_identity(bin)        ( ((float) bin)        / numIdentityBins)
#define bin_bottom_to_identity(bin) ((((float) bin) - 0.5) / numIdentityBins)
#define bin_top_to_identity(bin)    ((((float) bin) + 0.5) / numIdentityBins)

//----------
//
// prototypes for routines in identity_dist.c
//
//----------

alignel* filter_aligns_by_identity     (seq* seq1, seq* seq2, alignel* alignList, 
		                                float minIdentity, float maxIdentity);
void     alignment_identity            (seq* seq1, seq* seq2, alignel* a,
                                        unspos* numer, unspos* denom);
void     filter_segments_by_identity   (seq* seq1, seq* seq2, segtable* st, 
		                                float minIdentity, float maxIdentity);
int      filter_segment_by_identity    (seq* seq1, unspos pos1,
                                        seq* seq2, unspos pos2, unspos length,
		                                float minIdentity, float maxIdentity);
void     segment_identity              (seq* seq1, unspos pos1,
                                        seq* seq2, unspos pos2, unspos length,
                                        unspos* numer, unspos *denom);
unspos   count_substitutions           (seq* seq1, unspos pos1,
                                        seq* seq2, unspos pos2, unspos length,
                                        unspos count[4][4]);

alignel* filter_aligns_by_match_count  (seq* seq1, seq* seq2, alignel* alignList, 
		                                u32 minMatchCount);
void     filter_segments_by_match_count(seq* seq1, seq* seq2, segtable* st, 
		                                u32 minMatchCount);
int      filter_segment_by_match_count (seq* seq1, unspos pos1,
                                        seq* seq2, unspos pos2, unspos length,
		                                u32 minMatchCount);

alignel* filter_aligns_by_mismatch_count  (seq* seq1, seq* seq2, alignel* alignList, 
		                                s32 maxMismatchCount);
void     filter_segments_by_mismatch_count(seq* seq1, seq* seq2, segtable* st, 
		                                s32 maxMismatchCount);
int      filter_segment_by_mismatch_count(seq* seq1, unspos pos1,
                                        seq* seq2, unspos pos2, unspos length,
		                                s32 maxMismatchCount);

void     init_identity_dist_job        (seq* seq1, seq* seq2);
void     print_identity_dist_job       (FILE* f);
void     identity_dist_from_align_list (alignel* alignList,
                                        seq* seq1, seq* seq2);
void     identity_dist_from_match      (seq* seq1, unspos pos1,
                                        seq* seq2, unspos pos2, unspos length);

#undef global
#endif // identity_dist_H
